## HS-CaBLASTP Algorithm Specification

This specification defines the **Hierarchical Subsequence CaBLASTP (HS-CaBLASTP)** algorithm, which extracts subsequences based on local similarity and performs hierarchical compression and search using a multi-tree (Forest) structure.

### 1. Core Data Structures

These are the core classes and structures that must be maintained in memory and on disk to implement the algorithm.

```python
from typing import List, Dict, Tuple, Optional

# Represents an Edit Operation for the difference script
class EditOp:
    op_type: str       # 'INSERT', 'DELETE', or 'SUBSTITUTE'
    position: int      # Position index in the parent sequence
    char: str          # The amino acid character involved

# 1. Tree Node (Element of a subsequence cluster)
class TreeNode:
    node_id: int                    # Unique ID for the node
    is_root: bool                   # True if this is a Root node (stored as an actual sequence in Coarse DB)
    sequence: Optional[str]         # Actual amino acid string (ONLY if is_root == True), else None
    parent_id: Optional[int]        # ID of the parent node (None if Root)
    children: List[int]             # List of child node IDs
    diff_script: List[EditOp]       # Edit operations to transform the parent sequence into this sequence
    depth: int                      # Current depth in the tree (Root = 0)
    ref_original_seq: str           # Original FASTA ID and offset info (for mapping back to the source)

# 2. Seed Table (k-mer Index)
class SeedTable:
    # Key: k-mer string (e.g., 4-amino-acid sequence)
    # Value: List of Tuples containing (node_id, relative offset from the start of the Root)
    table: Dict[str, List[Tuple[int, int]]]

# 3. Compressed Database (Compression Output)
class CompressedDB:
    forest: Dict[int, TreeNode]     # The entire set of nodes, accessed by node_id
    coarse_fasta: List[str]         # A flat list/file of only Root node sequences (for 1st-pass search)

```

---

### 2. Phase 1: Hierarchical Compression & Tree Building

This phase reads the original protein sequence database and generates the tree-based compressed structure.

**Hyperparameters:**

* `K_MER_SIZE` = 4 (Seed length)
* `MIN_IDENTITY` = 0.7 (Minimum sequence identity of 70% between parent and child)
* `MIN_LENGTH` = 40 (Minimum length of the matching subsequence)
* `MAX_DEPTH` = 5 (Maximum tree depth to prevent reconstruction bottlenecks)

**Execution Process:**

1. **Initialization:** Initialize an empty `CompressedDB` and `SeedTable`.
2. **Sequence Traversal:** For an input original sequence $S$, begin scanning using a sliding window.
3. **Seed Matching:** * Check if the k-mer extracted from $S$ exists in the `SeedTable`.
* If it exists, retrieve the target node ($N_{target}$) associated with that k-mer.


4. **Extension & Alignment:**
* Attempt greedy **ungapped extension** to the left and right of the matched position.
* When ungapped extension can no longer proceed, switch to **gapped alignment** using the Needleman-Wunsch algorithm (with a constrained maximum gap limit).


5. **Node Creation Logic:**
* Check if the aligned subsequence meets the `MIN_LENGTH` and `MIN_IDENTITY` criteria.
* **Condition A (Criteria met AND $N_{target}$ depth < `MAX_DEPTH`):**
* Create a new `TreeNode` and register it as a child of $N_{target}$.
* Parse the alignment result to generate the insertion/deletion/substitution record and store it in `diff_script`.
* Do NOT store the actual sequence string (to save memory).


* **Condition B (Criteria NOT met OR depth has reached `MAX_DEPTH` OR no seed match):**
* Create a new Root node (`is_root = True`).
* Store the literal subsequence string in `sequence` and append it to `coarse_fasta`.
* Extract all k-mers from this new sequence and register them into the `SeedTable`.




6. **Dangling Ends Processing:** Repeat from Step 2 for the remaining unaligned segments (dangling ends) of $S$.

---

### 3. Phase 2: Search & Pruning

When a Query sequence ($Q$) is provided, the algorithm utilizes the tree structure to rapidly prune the search space.

**Step 1: Coarse Search (1st-pass Root Node Search)**

* Perform a standard BLASTP search comparing $Q$ against the `coarse_fasta` (the set of Root nodes).
* Use a relaxed E-value threshold (e.g., $10^{-3}$) to ensure potential candidates are not missed.
* Retrieve the list of matched Root nodes, denoted as $H_{root}$.

**Step 2: DFS/BFS Pruning Search (Tree Traversal)**

* For each matched Root node $R \in H_{root}$, traverse its subtree (DFS is recommended).
* Upon visiting a child node $C$, retrieve its `diff_script` (the variance from its parent).
* Perform **Heuristic Scoring** to rapidly estimate the similarity between $Q$ and $C$.
* *Logic:* Verify if the core matched region between $Q$ and the parent node was heavily disrupted (e.g., via massive deletions) by $C$'s `diff_script`.


* If the heuristic score falls below a predefined pruning threshold, **prune** node $C$ and entirely discard its descendant branches (grandchildren).
* Collect all surviving, valid nodes into a final candidate set, $H_{candidate}$.

---

### 4. Phase 3: Reconstruction & Fine Search

This phase perfectly reconstructs the original sequences from the pruned candidate pool and performs the final high-precision search.

**Step 1: Lossless Sequence Reconstruction**

* For each node $N \in H_{candidate}$, trace back its ancestry using `parent_id` until reaching the Root node.
* Collect the `diff_script` chain and apply the operations sequentially from the Root down to node $N$.
* Apply these script patches to the Root node's literal string to **100% losslessly reconstruct** the original amino acid sequence of $N$.

**Step 2: Fine Search (Final Precision Search)**

* Construct a small, temporary FASTA database in memory containing only the reconstructed sequences.
* Perform a standard BLASTP search between $Q$ and this temporary database using a strict, default E-value threshold (e.g., $10^{-10}$).
* Compute the final alignment results and scores. Use the `ref_original_seq` metadata to map the results back to the source proteins in the original database and return them to the user.