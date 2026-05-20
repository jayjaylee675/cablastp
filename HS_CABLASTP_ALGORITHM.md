# HS-CaBLASTP algorithm specification

This document describes the **Hierarchical Subsequence CaBLASTP** algorithm as it is actually implemented in `hs_cablastp/`. It is detailed enough to reconstruct the implementation: every data structure, every hyperparameter, every control-flow branch is captured here.

The original conceptual sketch lives in `hs_cablastp.md` and `debug.md`. Where the implementation diverges from those sketches (depth cap, descend-to-best, gapped fallback, pruning score), this document is authoritative.

---

## 1. Motivation

Classical `cablastp` deduplicates whole sequences: a long protein is split into matched segments stored as edit-scripts against a flat coarse database, plus the leftover residues.

HS-CaBLASTP generalises this to a **forest of trees**. A child node represents not just "this sequence is similar to that one" but "this *subsequence* is similar to a *subsequence* of an existing node, whether root or itself a child". Searching the forest needs only the root sequences in BLAST; descendants are pruned in-memory using their edit scripts.

The pipeline is:

```
       Phase 1 (offline)                  Phase 2 + 3 (online)
   ┌───────────────────┐                ┌───────────────────────┐
   │ INPUT FASTA       │                │ QUERY FASTA            │
   │  └─ compress      │                │  └─ blastp vs roots    │
   │     ├─ root #0    │                │     └─ prune subtree   │
   │     │  ├─ child   │                │        └─ reconstruct  │
   │     │  └─ child   │                │           └─ blastp    │
   │     └─ root #1    │                │              └─ hits   │
   └───────────────────┘                └───────────────────────┘
```

---

## 2. Core data structures

All defined in `hs_cablastp/types.py`.

### 2.1 `EditOp`

```python
@dataclass
class EditOp:
    op_type: str    # 'INSERT' | 'DELETE' | 'SUBSTITUTE'
    position: int   # 0-based offset in the parent's *matched segment*
    char: str       # residue (one of 20 amino acids); '' for DELETE
```

`position` is relative to the parent's matched segment `parent[parent_start:parent_end]`, **not** the parent's full sequence and **not** the child's reconstructed sequence.

### 2.2 `TreeNode`

```python
@dataclass
class TreeNode:
    node_id: int
    is_root: bool
    sequence: Optional[str]            # set IFF is_root else None
    parent_id: Optional[int]           # None IFF is_root
    children: List[int] = []
    diff_script: List[EditOp] = []
    depth: int = 0                     # root = 0
    ref_original_seq: str = ""         # "<fasta_id>:<start>-<end>"
    parent_start: int = 0              # range in parent's *full sequence* covered by this child
    parent_end: int = 0
```

Invariants:

- `is_root == True` ↔ `sequence` is the literal residue string and `parent_id is None`.
- `is_root == False` ↔ `sequence is None`, the child is reconstructable as `apply_edit_script(parent.sequence[parent_start:parent_end], diff_script)` (recursing for non-root parents).
- `depth(child) == depth(parent) + 1`.
- `depth(root) == 0` and `depth(node) ≤ MAX_DEPTH`.
- `ref_original_seq` is the slice of the original input FASTA record this node represents, in the form `"<fasta_id>:<start>-<end>"` with half-open `[start, end)` indices.

### 2.3 `SeedTable`

```python
class SeedTable:
    k: int
    table: Dict[str, List[Tuple[int, int]]]   # kmer -> [(root_node_id, offset_in_root), ...]

    def add_root_sequence(node_id, sequence):
        for i in 0 .. len(sequence) - k:
            kmer = sequence[i:i+k]
            if '*' in kmer or '-' in kmer: continue
            table[kmer].append((node_id, i))

    def lookup(kmer): return table.get(kmer, [])
```

**Only root sequences populate the seed table.** Descendants are reached via tree descent (§3.5), not seed lookup.

### 2.4 `CompressedDB`

```python
class CompressedDB:
    forest: Dict[int, TreeNode]
    _next_id: int

    def new_node(**kwargs) -> TreeNode:
        node = TreeNode(node_id=_next_id, **kwargs)
        forest[_next_id] = node
        _next_id += 1
        return node

    def roots(): yield each node where is_root
```

Node IDs are monotonically increasing integers across the whole forest.

---

## 3. Phase 1: Compression

Implemented in `hs_cablastp/compress.py` as `_Compressor.compress_sequence`.

### 3.1 Hyperparameters

Spec defaults (from `compress.py`):

```python
K_MER_SIZE   = 4
MIN_IDENTITY = 0.7
MIN_LENGTH   = 40
MAX_DEPTH    = 5
```

Implementation knobs (also in `compress.py`):

```python
MAX_SEED_HITS = 32   # only the first 32 hits per seed are examined
GAPPED_BAND   = 10   # band width for banded Needleman-Wunsch
EXTEND_DROP   = 12   # X-drop for ungapped extension
```

### 3.2 The driver loop

For each FASTA record, call `compress_sequence(fasta_id, seq)`. Edge cases first:

```
if len(seq) == 0:        return                       # nothing to compress
if len(seq) <  K_MER_SIZE: make_root over [0, len)    # too short to seed
```

Otherwise, slide a window over `seq` of length `K_MER_SIZE`:

```
pos = 0
pending_root_start = 0
while pos + k <= N:
    kmer = seq[pos:pos+k]
    hits = seeds.lookup(kmer)
    match = best_match(seq, pos, hits) if hits else None

    if match is not None:
        # 1. Flush any unmatched residues since the last match as a NEW ROOT.
        if match.q_start > pending_root_start:
            make_root(fasta_id, seq, pending_root_start, match.q_start)

        # 2. Attach the matched region as a child (or new root if depth cap hit).
        if forest[match.target_id].depth < MAX_DEPTH:
            attach_child(match)
        else:
            make_root(fasta_id, seq, match.q_start, match.q_end)

        pos = match.q_end
        pending_root_start = pos
    else:
        pos += 1

# Final flush: residues after the last match also become a root.
if N > pending_root_start:
    make_root(fasta_id, seq, pending_root_start, N)
```

**Critical invariant: no residue of the input is ever discarded.** Every position in `seq` ends up either inside a child node's `diff_script` (matched region) or inside a root node's `sequence` field (unmatched region). The two flushes above are what enforce this — see `debug.md` "Dangling Ends".

### 3.3 `_best_match`: from seed hit to candidate alignment

Given a seed hit `(root_id, root_off)` and the current query position `pos`:

```python
def best_match(seq, pos, hits) -> Optional[_Match]:
    best = None
    for root_id, root_off in hits[:MAX_SEED_HITS]:
        root_seq = forest[root_id].sequence

        # (a) Ungapped X-drop extension around the seed.
        (sa, ea, sb, eb) = ungapped_extend(seq, pos, root_seq, root_off,
                                           K_MER_SIZE, max_drop=EXTEND_DROP)
        q_sub = seq[sa:ea]
        t_sub = root_seq[sb:eb]

        # (b) If the ungapped window is too short, widen symmetrically
        #     until it covers MIN_LENGTH residues and fall through to NW.
        if len(q_sub) < MIN_LENGTH:
            (q_sub, t_sub, sa, ea, sb, eb) = gapped_extend(seq, root_seq, sa, ea, sb, eb)
            if len(q_sub) < MIN_LENGTH:
                continue

        # (c) Banded Needleman-Wunsch on the widened window.
        (parent_aln, child_aln) = needleman_wunsch(t_sub, q_sub, band=GAPPED_BAND)
        identity = alignment_identity(parent_aln, child_aln)
        if identity < MIN_IDENTITY:
            continue

        cand = _Match(target_id=root_id, target_seq=root_seq,
                      target_start=sb, target_end=eb,
                      q_start=sa,  q_end=ea,
                      identity=identity,
                      parent_aln=parent_aln, child_aln=child_aln)

        # (d) Refine: try existing descendants of root that align better.
        cand = descend_to_best(cand)

        # (e) Score: identity × length. Larger wins.
        if best is None or cand.identity * cand.length > best.identity * best.length:
            best = cand
    return best
```

### 3.4 Ungapped X-drop extension

`ungapped_extend(a, ia, b, ib, k, max_drop)` (in `hs_cablastp/alignment.py`):

1. Seed the alignment with the `k`-mer at `(ia, ib)` on the two strings. Score it via BLOSUM62 (`cablastp.blosum.MATRIX62`).
2. Walk right one position at a time: `cur += score(a[pa], b[pb])`. Track the best score and the position where it occurred. Stop when `best - cur > max_drop` or either string ends.
3. Walk left the same way from `(ia-1, ib-1)`.

Returns half-open intervals `(start_a, end_a, start_b, end_b)`. The reported start/end correspond to the best-score frontier, not the final cursor.

### 3.5 `_descend_to_best`: walking into the subtree

After (a)-(c) find a viable match against a root, the compressor tries to find a *descendant* of that root whose reconstructed sequence aligns even better — this is what lets new sequences attach deep in the tree without needing seeds on descendants.

```python
def descend_to_best(match) -> _Match:
    best = match
    current_id = match.target_id
    while True:
        current = forest[current_id]
        if no children or current.depth >= MAX_DEPTH - 1:
            break
        improved = False
        for child_id in current.children:
            child_seq = reconstruct(child_id)
            q_sub = best.child_aln without gaps      # the query slice that produced `best`
            (p_aln, c_aln) = needleman_wunsch(child_seq, q_sub, band=GAPPED_BAND)
            ident = alignment_identity(p_aln, c_aln)
            if ident >= MIN_IDENTITY and ident > best.identity:
                best = _Match(target_id=child_id, target_seq=child_seq,
                              target_start=0, target_end=len(child_seq),
                              q_start=match.q_start, q_end=match.q_end,
                              identity=ident,
                              parent_aln=p_aln, child_aln=c_aln)
                current_id = child_id
                improved = True
                break   # restart descent from new child
        if not improved:
            break
    return best
```

The loop stops at the *first* improving child it finds (depth-first, take-best-or-break-out), so descent is linear in tree depth and is not exhaustive. The depth check `current.depth >= MAX_DEPTH - 1` ensures the new child node (one level deeper) won't exceed `MAX_DEPTH`.

### 3.6 Banded Needleman-Wunsch

`needleman_wunsch(a, b, band)`:

- BLOSUM62 substitution scores; the gap penalty is `MATRIX62[GAP][GAP]` (treating the last row/column as the gap pseudo-residue).
- DP table `(n+1) × (m+1)` initialised to `NEG = -10**9` except `[0][0] = 0`.
- Initialise the first row and column with `gap_penalty * i`, but stop at `band` if banding is on (cells outside the band stay at `NEG`).
- Fill row `i` from `j = max(1, i-band)` to `min(m, i+band)` (if banded), otherwise full row. `table[i][j] = max(diag + sub, up + gap, left + gap)`.
- Traceback from `(n, m)` following the same `max` decisions, breaking ties as **diag > up > left**.

Returns the two aligned strings, each with `'-'` for gaps.

### 3.7 Building and applying edit scripts

`make_edit_script(parent_aln, child_aln)` walks the two aligned strings left-to-right, tracking `parent_pos` = number of non-gap parent residues seen so far:

| parent_aln | child_aln | emit                                  | parent_pos++ |
| ---------- | --------- | ------------------------------------- | ------------ |
| `-`        | `X`       | `INSERT(parent_pos, X)`               | no           |
| `X`        | `-`       | `DELETE(parent_pos, '')`              | yes          |
| `X`        | `Y` (≠ X) | `SUBSTITUTE(parent_pos, Y)`           | yes          |
| `X`        | `X`       | (nothing — match)                     | yes          |

`apply_edit_script(parent_segment, script)` is the inverse:

```python
ops_by_pos = group script by op.position
out = []
for i in 0 .. len(parent_segment):
    for op at position i:
        if op.op_type == INSERT: out.append(op.char)
    if i == len(parent_segment): break
    applied = False
    for op at position i:
        if op.op_type == DELETE:     applied = True; break        # skip the residue
        if op.op_type == SUBSTITUTE: out.append(op.char); applied = True; break
    if not applied: out.append(parent_segment[i])
return ''.join(out)
```

INSERTs at the same position keep their original order from the script. Multiple ops at the same position are unusual but supported — at most one DELETE/SUBSTITUTE plus any number of INSERTs.

### 3.8 Reconstruction cache

`_Compressor.reconstruct(node_id)` is memoised in `_seq_cache: Dict[int, str]`. It is needed by `_descend_to_best` (so descendants get reconstructed once per compression) and by the search phase. Cache eviction is not currently needed because roots are immutable and children are append-only during compression.

### 3.9 `_make_root`

```python
def _make_root(fasta_id, seq, start, end):
    sub = seq[start:end]
    node = forest.new_node(is_root=True, sequence=sub, parent_id=None, depth=0,
                           ref_original_seq=f"{fasta_id}:{start}-{end}")
    seeds.add_root_sequence(node.node_id, sub)
    return node
```

Adding the root to the seed table is what makes that subsequence available as an attachment point for all subsequently compressed sequences. **Children's k-mers are not seeded.** Children are reached only via the descend-to-best refinement step.

---

## 4. Phase 2: Coarse search + tree pruning

Implemented in `hs_cablastp/search.py:search`.

### 4.1 Hyperparameters

```python
COARSE_EVALUE   = 1e-3
FINE_EVALUE     = 1e-10
PRUNE_THRESHOLD = 0.5    # heuristic score in [0, 1]
```

### 4.2 Step 1 — Coarse `blastp` against the root FASTA

The on-disk database (see §6) stores all root sequences as `>node_<id> <ref_original_seq>` records in `coarse.fasta` and a pre-built BLAST DB at `blastdb-coarse.*`.

```python
rows = run_blastp(blastp, query_path, db_dir / "blastdb-coarse",
                  evalue=COARSE_EVALUE,
                  outfmt="6 qseqid sseqid pident length mismatch gapopen "
                         "qstart qend sstart send evalue bitscore")
```

Each row tells us a query HSP against a specific root. We collect:

```python
@dataclass
class _RootHit:
    root_id: int     # parsed from sseqid "node_<id>"
    q_start: int     # 0-based half-open on query
    q_end:   int
    s_start: int     # 0-based half-open on root.sequence
    s_end:   int
    evalue:  float
    bitscore: float
```

BLAST coordinates are 1-based inclusive; the search code converts to 0-based half-open as it parses.

### 4.3 Step 2 — Tree pruning

For each root hit, DFS the subtree of `root_id`, carrying down a "still-covered range" `q_range` (initially `(s_start, s_end)`, expressed in coordinates of the *current* node's sequence). At each child:

1. Compute the intersection of the child's `[parent_start, parent_end)` window with the inherited `q_range`. If empty, prune.
2. Count "disruption" ops in `diff_script` whose `parent_pos = parent_start + op.position` falls inside the intersection. Weights:
   - `DELETE` → 2
   - `SUBSTITUTE`, `INSERT` → 1
3. Score: `score = 1 - disruptions / max(window_len, 1)`. If `score < PRUNE_THRESHOLD`, prune.
4. Otherwise descend with the new range, **translated into the child's coordinate system**: `(isect_a - parent_start, isect_b - parent_start)`.

Pseudocode:

```python
def heuristic_keep(child, q_range_on_parent, threshold):
    isect_a = max(child.parent_start, q_range_on_parent[0])
    isect_b = min(child.parent_end,   q_range_on_parent[1])
    if isect_a >= isect_b: return False, (0, 0)
    window_len = isect_b - isect_a
    disruptions = 0
    for op in child.diff_script:
        abs_pos = child.parent_start + op.position
        if isect_a <= abs_pos < isect_b:
            disruptions += 2 if op.op_type == 'DELETE' else 1
    score = 1.0 - disruptions / max(window_len, 1)
    if score < threshold: return False, (0, 0)
    return True, (isect_a - child.parent_start, isect_b - child.parent_start)

def prune_tree(db, root_id, q_range_on_root, threshold):
    out = []
    stack = [(root_id, q_range_on_root)]
    while stack:
        node_id, q_range = stack.pop()
        out.append(CandidateHit(node_id, root_id, inherited_s_range=q_range))
        for child_id in db.forest[node_id].children:
            keep, new_range = heuristic_keep(db.forest[child_id], q_range, threshold)
            if keep:
                stack.append((child_id, new_range))
    return out
```

**The root itself is always kept** (it's pushed onto `out` before any pruning). Deduplication is by `node_id`, keeping the widest `inherited_s_range`:

```python
candidates: Dict[int, CandidateHit] = {}
for root_hit in root_hits:
    for cand in prune_tree(db, root_hit.root_id, (root_hit.s_start, root_hit.s_end), threshold):
        prev = candidates.get(cand.node_id)
        if prev is None or width(cand) > width(prev):
            candidates[cand.node_id] = cand
```

---

## 5. Phase 3: Reconstruction + fine BLAST

### 5.1 Reconstruction

For each surviving candidate, `reconstruct_sequence(db, node_id)` walks parent pointers up to the root, then applies edit scripts on the way down:

```python
def reconstruct_sequence(db, node_id):
    chain = []
    cur = db.forest[node_id]
    while cur is not None:
        chain.append(cur)
        cur = db.forest[cur.parent_id] if cur.parent_id is not None else None
    chain.reverse()                       # [root, ..., target]
    seq = chain[0].sequence
    for node in chain[1:]:
        segment = seq[node.parent_start:node.parent_end]
        seq = apply_edit_script(segment, node.diff_script)
    return seq
```

The function recomputes from scratch each call (unlike `_Compressor.reconstruct` which is memoised) because `search` runs after `load_db` and there is no `_Compressor` instance around.

### 5.2 Fine BLAST

Write all reconstructed candidate sequences to a temp FASTA:

```
>node_<id> <ref_original_seq>
<60-col wrapped reconstructed sequence>
...
```

Build a BLAST DB with `makeblastdb -dbtype prot`. Run `blastp` with `--fine-evalue` (default 1e-10) and the same tabular `outfmt 6` as Phase 2.

Each fine BLAST row becomes a `SearchHit`:

```python
@dataclass
class SearchHit:
    node_id: int
    fasta_ref: str       # = TreeNode.ref_original_seq
    pident: float
    length: int
    qstart: int          # 1-based BLAST coords (left as-is)
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
```

Note: fine-hit coords stay in BLAST's 1-based inclusive convention; coarse-hit coords were converted to 0-based half-open in §4.2. This asymmetry is by design — fine hits are user-facing.

### 5.3 Result envelope

`search(...)` returns a dict:

```python
{
    "coarse_hits":    int,                # number of coarse blast rows
    "matched_roots":  int,                # distinct root nodes hit
    "candidates":     int,                # candidates after pruning
    "fine_hits":      List[SearchHit],
}
```

The CLI then prints summary counts plus up to 50 fine hits.

---

## 6. Disk format

`hs_cablastp/io.py`. Layout under `<db_dir>/`:

```
forest.pkl          pickle.HIGHEST_PROTOCOL of CompressedDB.forest (Dict[int, TreeNode])
meta.pkl            pickle of {"params": <dict>, "next_id": int}
coarse.fasta        FASTA of root sequences (one record per root, 60-col wrapped)
                    Defline: ">node_<id> <ref_original_seq>"
blastdb-coarse.*    output of `makeblastdb -dbtype prot -in coarse.fasta -out blastdb-coarse`
```

`save_db` always re-writes the FASTA and re-runs `makeblastdb`. `load_db` reads `forest.pkl` and `meta.pkl` only — the BLAST DB is consumed by `blastp` directly.

**`params` content**: `{"k", "min_identity", "min_length", "max_depth"}`. The search phase does not currently re-use these values (the search hyperparameters are independent CLI flags) but storing them is necessary for reproducibility and future variants of the seed/match logic.

---

## 7. Reconstruction guarantee

The forest preserves **every residue of every input sequence**. The proof:

- During Phase 1 compression, every position of `seq` is either part of a flushed root (residues from `pending_root_start..match.q_start` and `pending_root_start..N`) or inside a child node's edit script (residues `match.q_start..match.q_end`, applied to `target.sequence[match.target_start:match.target_end]`).
- Roots store the literal substring; children's `diff_script` against a known parent segment is invertible by `apply_edit_script` (every `INSERT`/`DELETE`/`SUBSTITUTE`/match operation has a definite effect).

So every input sequence can be recovered exactly by enumerating the children whose `ref_original_seq` shares the same `fasta_id`, sorting their `[start, end)` ranges (with the inserted-root ranges) and concatenating their reconstructed segments. (No decompressor for HS-CaBLASTP is implemented yet — the data is sufficient, the function is just absent.)

---

## 8. Hyperparameter reference

| Constant | Module | Default | Used in |
| --- | --- | --- | --- |
| `K_MER_SIZE` | `compress.py` | 4 | seed lookup |
| `MIN_IDENTITY` | `compress.py` | 0.7 | reject attachments below this NW identity |
| `MIN_LENGTH` | `compress.py` | 40 | reject windows shorter than this |
| `MAX_DEPTH` | `compress.py` | 5 | cap on tree depth |
| `MAX_SEED_HITS` | `compress.py` | 32 | examine at most N hits per seed |
| `GAPPED_BAND` | `compress.py` | 10 | NW band radius |
| `EXTEND_DROP` | `compress.py` | 12 | X-drop for ungapped extension |
| `COARSE_EVALUE` | `search.py` | 1e-3 | E-value for `blastp` vs roots |
| `FINE_EVALUE` | `search.py` | 1e-10 | E-value for fine `blastp` |
| `PRUNE_THRESHOLD` | `search.py` | 0.5 | minimum heuristic score to keep a child |
| `BLOSUM62 gap penalty` | `cablastp/blosum.py` | `MATRIX62[GAP][GAP]` | substitution scores |

---

## 9. Implementation notes & known divergences from `debug.md`

1. **Depth cap behaviour.** `debug.md` says "if depth limit reached, create a new root". Implementation: when the *target's* depth equals `MAX_DEPTH`, the matched region becomes a new root (matching the spec). When walking children in `_descend_to_best`, the loop exits one step early (`current.depth >= MAX_DEPTH - 1`) so the new child stays within the cap.
2. **Gapped extension is greedy and symmetric.** `debug.md` mentions "constrained gapped alignment (Needleman-Wunsch)". Implementation: we widen the ungapped window symmetrically to exactly `MIN_LENGTH` residues (`_gapped_extend`) and then run banded NW once. There's no iterative drop-out.
3. **Pruning score.** `debug.md` says "apply a heuristic scoring function using its diff_script against Q". Implementation: it doesn't look at `Q` at all — it counts edit operations in the parent's coordinate window. DELETEs cost double. This is intentional: the query alignment to the parent already covers that window, so the question is "how much does the child diverge from the parent in the matched region".
4. **Seed table is roots-only.** Mentioned in §3.9 above. Descendants are only reachable via `_descend_to_best` during compression and via DFS pruning during search.
5. **No multi-process compression.** Unlike classical `cablastp-compress`, the HS variant is single-threaded today. Parallelism would require partitioning the forest (the seed table is shared mutable state).
6. **No decompressor CLI.** The data is sufficient (§7) but no `hs-cablastp-decompress` script is exposed in `pyproject.toml`.
