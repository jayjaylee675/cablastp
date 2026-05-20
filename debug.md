**Objective:** Implement the "Hierarchical Subsequence CaBLASTP" (HS-CaBLASTP) algorithm. This algorithm compresses protein sequence databases using a forest of trees based on local subsequence similarity, enabling massive pruning during searches and lossless sequence reconstruction.

### 1. Core Data Structures

Please define the following data structures:

* **`EditOp`**: Represents a single edit operation (e.g., `type` (insert/delete/sub), `position`, `char`).
* **`TreeNode`**: Represents a sequence cluster. Fields: `node_id`, `is_root` (boolean), `sequence` (string, ONLY if root, otherwise null), `parent_id`, `children` (list of IDs), `diff_script` (list of `EditOp`s comparing this node to its parent), `depth`, and `ref_original_seq` (metadata mapping).
* **`SeedTable`**: A hash map indexing k-mers to a list of `(node_id, offset_position)` tuples.
* **`CompressedDB`**: Contains `forest` (a dictionary of all `TreeNode` objects by ID) and `coarse_fasta` (a list/file of only the root node sequences).

### 2. Phase 1: Compression (Forest Building)

Implement a sliding-window compression function with the following constraints: `K_MER_SIZE = 4`, `MIN_IDENTITY = 0.7`, `MIN_LENGTH = 40`, `MAX_DEPTH = 5`.

* **Sliding Window & Seeding:** Scan the input sequence $S$. If a k-mer matches the `SeedTable`, fetch the target node ($N_{target}$).
* **Alignment:** Perform greedy ungapped extension. If it stops, fallback to constrained gapped alignment (Needleman-Wunsch).
* **Node Creation Logic:**
* *Condition A (Child):* If alignment length >= `MIN_LENGTH`, identity >= `MIN_IDENTITY`, AND $N_{target}$ depth < `MAX_DEPTH`: Create a child node. Save **only** the `diff_script` (parsed from the alignment). Do NOT store the full string.
* *Condition B (Root):* If conditions fail, depth limit is reached, or no seed matches: Create a new Root node. Store the literal substring in `sequence` and `coarse_fasta`. Update the `SeedTable` with all k-mers from this new root.


* **CRITICAL - Dangling Ends:** Unaligned leftover segments of $S$ MUST NOT be discarded. Process them recursively as new roots to prevent data loss.

### 3. Phase 2: Search & Pruning

Implement the search logic to quickly discard irrelevant branches.

* **Coarse Search:** Run a standard BLASTP search of the Query ($Q$) against `coarse_fasta` (Roots only) using a relaxed E-value (e.g., 1e-3).
* **Tree Pruning (DFS):** For each matched Root, traverse its subtree. At each child node, apply a heuristic scoring function using its `diff_script` against $Q$.
* If the heuristic score falls below a threshold (e.g., major deletions in the matched region), **prune** the node and all its descendants immediately. Collect surviving nodes into a candidate list.

### 4. Phase 3: Lossless Reconstruction & Fine Search

* **Reconstruction:** For each surviving candidate node, trace its `parent_id` back to the Root. Sequentially apply the `diff_script` patches from the Root downwards to 100% losslessly reconstruct the candidate's original amino acid sequence.
* **Fine Search:** Create a temporary in-memory database with only these reconstructed sequences. Perform a strict BLASTP search (E-value 1e-10) against this subset and return the final hits mapped to `ref_original_seq`.