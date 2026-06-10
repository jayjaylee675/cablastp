# Why DB size and search speed move the way they do — an algorithmic analysis

This note explains, at the algorithm level, **why** the measured DB size and search
time differ between `cablastp` and `hs-cablastp`, and **why** the differences flip
between the redundant `dense_2k` set and the divergent `etrembl_2k` set. All numbers
are from the final configuration (hs-cablastp with short-flank absorption, no boundary
overlap; both pipelines benchmarked with DB-load excluded and equal coarse e-value).

The single idea that explains everything below:

> **The redundancy structure of the input decides the *shape* of the hs-cablastp
> forest, and the forest shape decides both the DB size and the size of the
> fine-search candidate set.** A deep forest (few roots, many children) compresses
> hard but expands to many fine candidates; a flat forest (many roots, few children)
> barely compresses but expands to few candidates.

---

## 1. What each pipeline actually stores and searches

### cablastp (whole-sequence redundancy removal)
- **Coarse DB** (`coarse.fasta` + BLAST index): a set of *representative whole
  sequences*. Every input sequence is either a representative or is stored as an
  edit-script + link against one (`compressed`, `coarse.links`).
- **Search**: coarse BLAST (query vs representatives) → for every coarse hit,
  *decompress* the linked originals back to full sequences (`expand_blast_hits`,
  pure-Python link walking) → build a fine BLAST DB of those originals → fine BLAST.

### hs-cablastp (hierarchical sub-sequence forest)
- **Coarse DB** (`coarse.fasta` + BLAST index): only the **root** node sequences.
- **Forest blob** (`forest.pkl`): every **child** node, stored as a compact
  `diff_script` (edit ops) against its parent's matched segment — *no literal
  sequence*.
- **Search**: coarse BLAST (query vs roots only) → **tree pruning** (in-memory DFS
  that keeps/drops descendants by a cheap diff-script disruption score) → reconstruct
  the surviving nodes via `apply_edit_script` → build a fine BLAST DB → fine BLAST.

The crucial asymmetry: in hs-cablastp a *child* costs only its diff-script, and a
*root* costs a full literal sequence (plus its BLAST index). **So DB size is dominated
by how many roots there are, and search cost is dominated by how many nodes survive
pruning into the fine DB.** Both of those are set by the forest shape.

---

## 2. Forest shape is data-dependent (measured)

| forest (hs-cablastp) | dense_2k | etrembl_2k |
| --- | ---: | ---: |
| total nodes | 3474 | 2321 |
| roots | **538** | **1816** |
| children | 2936 | 505 |
| child ratio | **0.85** | **0.22** |
| depth histogram (d0..d5) | 538 / 813 / 780 / 746 / 391 / 206 | 1816 / 443 / 55 / 7 / 0 / 0 |

- **dense_2k** is built from near-identical strain proteins. A new sequence almost
  always seed-matches an existing node at ≥ the identity threshold, so it attaches as
  a **child** (a small diff-script) and even nests several levels deep. Result: a
  **deep forest** — 85% of nodes are children, and the tree fills all 5 depth levels.
- **etrembl_2k** is a uniform, divergent sample. New sequences rarely clear the
  identity/length gate against any existing node, so they are flushed as **new roots**.
  Result: a **flat forest** — 78% of nodes are roots, almost nothing nests past depth 1.

Everything in §3 and §4 follows from these two shapes.

---

## 3. DB size — why hs is 81% of cablastp on dense but 97% on etrembl

Measured DB component breakdown:

| component | dense cab | dense hs | etrembl cab | etrembl hs |
| --- | ---: | ---: | ---: | ---: |
| coarse FASTA (literal seqs) | 239 KB | **148 KB** | 499 KB | **569 KB** |
| compressed / forest blob | 328 KB | 230 KB | 266 KB | 140 KB |
| **total DB (incl. BLAST index)** | **1918 KB** | **1547 KB (81%)** | **2467 KB** | **2392 KB (97%)** |
| (raw input FASTA for reference) | ~1102 KB | | ~838 KB | |

**Why hs wins big on dense:** the coarse FASTA stores only roots. With just 538 roots,
dense's coarse FASTA is tiny (148 KB) and its BLAST index (built from that FASTA) is
correspondingly small. The 2936 redundant children cost only diff-scripts (230 KB
blob). The literal-sequence footprint is pushed down hard because the data genuinely
is mostly repeats of a few hundred root sequences.

**Why the win nearly vanishes on etrembl:** with 1816 roots, the coarse FASTA (569 KB)
is essentially *the whole input stored verbatim* — divergent sequences could not be
folded into children, so there is almost nothing to compress. hs's coarse FASTA is in
fact slightly *larger* than cablastp's (569 vs 499 KB); hs only edges out on the total
because cablastp carries more auxiliary index/link files.

**Conclusion:** hs-cablastp's DB-size advantage is *proportional to the child ratio*,
which is *proportional to input redundancy*. It is a real, robust win — but only on
redundant databases.

---

## 4. Search speed — why hs is slower on dense but faster on etrembl

Algorithm-only search time (DB load excluded, equal coarse e-value, mean of 5 runs):

| | dense cab | dense hs | etrembl cab | etrembl hs |
| --- | ---: | ---: | ---: | ---: |
| search time | 5.9 s | **6.4 s (slower)** | 4.0 s | **3.3 s (faster)** |
| fine candidates (seqs in fine DB) | 1417 | **2623** | 119 | 144 |
| hs coarse hits / candidates after prune | — | 363 / 2623 | — | 132 / 144 |

Two opposing algorithmic forces decide the winner; which one dominates depends on the
forest shape.

### Force A — coarse-stage processing (favours hs, always)
cablastp's coarse stage is heavier *per coarse hit*: it BLASTs in XML (`outfmt 5`) and
then runs `expand_blast_hits`, a pure-Python walk that decompresses each linked
original from disk. hs's coarse stage is lean: tabular BLAST against roots, then an
**in-memory** tree-prune (counting diff-script disruptions — no sequence work) before
any reconstruction. When coarse hits are few, this fixed-overhead difference is the
whole story.

### Force B — fine-candidate-set size (favours cablastp on redundant data)
The fine BLAST cost scales with the number of sequences written into the fine DB.
- On **dense**, hs's deep forest is its own enemy here: one coarse hit on a root drags
  in a whole subtree of near-identical descendants, so 363 coarse hits **expand to
  2623 fine candidates** — nearly 2× cablastp's 1417. hs's fine BLAST is therefore the
  heaviest single phase, and it loses the overall race despite the leaner coarse step.
- On **etrembl**, the flat forest barely expands: 132 coarse hits → only 144
  candidates (≈ cablastp's 119). Force B is essentially neutral, so Force A wins and
  hs comes out ahead.

### Net result
- **dense (deep forest): Force B dominates → hs slower (+8%).** The same redundancy
  that shrinks the DB inflates the candidate set, because the forest keeps many
  near-duplicate reconstructions of the same region.
- **etrembl (flat forest): Force A dominates → hs faster (−18%).** Few candidates for
  either pipeline, so hs's lean in-memory coarse/prune beats cablastp's XML +
  link-decompression.

This is also why the earlier (unfair) measurement made hs look uniformly faster: it
included cablastp's one-time multi-file DB open (~0.5 s) and let cablastp use a laxer
coarse e-value. Removing both confounds exposes Force B, which genuinely hurts hs on
redundant data.

---

## 5. The unifying picture

```
            input redundancy
                   │
                   ▼
        ┌──────────────────────┐
        │   hs-cablastp forest │
        │        shape         │
        └──────────────────────┘
          /                    \
   high redundancy        low redundancy
   = DEEP forest          = FLAT forest
   (few roots,            (many roots,
    many children)         few children)
        │                      │
        ├── small coarse FASTA ─┤ large coarse FASTA
        │   → DB much smaller   │ → DB ≈ cablastp
        │                       │
        └── big subtrees →      └── tiny subtrees →
            many fine candidates    few fine candidates
            → search SLOWER         → search FASTER
```

**DB size and search speed are driven by the same lever (forest depth) in opposite
directions.** Redundant data is exactly where hs-cablastp compresses best *and*
searches worst, because both effects are caused by the same thing: many near-identical
sequences collapsing into deep subtrees.

## 6. Practical implications

- hs-cablastp's clear, robust advantage is **DB size on redundant databases** (strain
  collections, clustered families). The more redundant the corpus, the bigger the
  compression win.
- Its search is **not** algorithmically faster in general; it is a data-dependent
  trade-off (slower on redundant, faster on divergent data).
- The obvious lever to fix the dense-search penalty is to **shrink the fine-candidate
  set without losing recall** — e.g. deduplicate identical reconstructions before the
  fine BLAST, or cap/merge near-identical nodes in a hit subtree. That attacks Force B
  directly while leaving the DB-size win intact.
