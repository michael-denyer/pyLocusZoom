/-!
# LD-heatmap highlight cells and cell edges

Source: `src/pylocuszoom/backends/composition.py`

* `lower_triangle`           l.172-188 (mask at l.187, `np.triu(..., k=1)`)
* `heatmap_highlight_cells`  l.191-211 (guard l.207, cells l.209-210)
* `cell_edges`               l.214-232 (single branch l.226-227, arithmetic l.228-232)
* `heatmap_highlight_rects`  l.235-258 (indexing l.256)
* `draw_ld_heatmap`          l.261-312 (`y_coords = list(range(len(x_coords)))`, l.287)

## Conventions fixed by reading the backends

Every backend draws `data[row][col]` at `(x_coords[col], y_coords[row])`:
matplotlib `pcolormesh(X, Y, C)` (`matplotlib_backend.py:540-549`), plotly
`go.Heatmap(z, x, y)` (`plotly_backend.py:797-800`), and bokeh's explicit
`xs = x_edges[j]`, `ys = y_edges[i]` for `data[i, j]`
(`bokeh_backend.py:741-746`). `heatmap_highlight_rects` reads `x_edges[x]`,
`y_edges[y]` (l.256), so a highlight cell `(x, y)` is `(column, row)`.
`lower_triangle` masks `column > row`, so a cell is rendered iff `x ≤ y`.

## Number model

* Indices are Python ints. Inputs to `heatmap_highlight_cells` are `Int`
  because the guard tests `snp_idx < 0` and `n_snps < 1`; cells are `Nat`.
* Coordinates are Python floats. They are modelled as `Int`, and every
  edge is stored DOUBLED so the midpoint `(a + b) / 2` is exact:
  `2 * mid = a + b`. The single-coordinate branch's `± 0.5` is `± 1`
  doubled. A centre `c` is compared against doubled edges as `2 * c`.
  This is exact for the real callers, which pass Python ints (genomic
  positions or `range(n)`) far below 2^52; it does not model float
  rounding of non-integer coordinates.
* A Python exception (`ValueError` l.208, `IndexError` l.229/l.256) is `none`.
-/

/-! ## Transcription -/

/-- `lower_triangle` mask, `composition.py:187`: `np.triu(ones, k=1)` masks
entry `(row, col)` iff `col > row`. -/
def masked (row col : Nat) : Bool := decide (col > row)

/-- `heatmap_highlight_cells`, `composition.py:191-211`. Cells are `(x, y)`.
`range(snp_idx + 1, n_snps)` is `List.range' (s + 1) (n - (s + 1))`. -/
def heatmap_highlight_cells (snp_idx n_snps : Int) : Option (List (Nat × Nat)) :=
  if n_snps < 1 ∨ snp_idx < 0 ∨ snp_idx ≥ n_snps then none            -- l.207-208
  else
    let s := snp_idx.toNat
    let n := n_snps.toNat
    some ((List.range (s + 1)).map (fun j => (j, s))                   -- l.209
          ++ (List.range' (s + 1) (n - (s + 1))).map (fun i => (s, i))) -- l.210

/-- Doubled midpoints, `composition.py:228`:
`mids = [(a + b) / 2 for a, b in zip(coords, coords[1:])]`. -/
def mids (coords : List Int) : List Int :=
  (coords.zip coords.tail).map (fun p => p.1 + p.2)

/-- `cell_edges`, `composition.py:214-232`, edges doubled.
The empty list reaches `coords[0]` / `mids[0]` at l.229 and raises
`IndexError`. No plot call reaches that arm: `prepare_ld_matrix`
(`_ld_matrix.py:23-26`) rejects a matrix with no SNPs at intake, and
`HeatmapPanel.from_matrix` (`panels/heatmap.py:63-66`) rejects a regional
heatmap with no SNP in the region. In the last arm `mids` is non-empty (`mids_length`), so the
`getD` defaults are never used. -/
def cell_edges (coords : List Int) : Option (List (Int × Int)) :=
  match coords with
  | [] => none                                                         -- l.229
  | [c] => some [(2 * c - 1, 2 * c + 1)]                               -- l.226-227
  | c0 :: c1 :: t =>
    let coords := c0 :: c1 :: t
    let m := mids coords                                               -- l.228
    let first := 2 * c0 - (m.getD 0 0 - 2 * c0)                        -- l.229
    let cl := coords.getD (coords.length - 1) 0                        -- coords[-1]
    let last := 2 * cl + (2 * cl - m.getD (m.length - 1) 0)            -- l.230
    let bounds := first :: (m ++ [last])                               -- l.231
    some (bounds.zip bounds.tail)                                      -- l.232

/-- Run `f` over a list; `none` as soon as one call raises. -/
def allOrNone {α β : Type} (f : α → Option β) : List α → Option (List β)
  | [] => some []
  | a :: t =>
    match f a, allOrNone f t with
    | some b, some bs => some (b :: bs)
    | _, _ => none

/-- Loop body of `heatmap_highlight_rects`, `composition.py:256-257`.
Result is `(x0, y0, width, height)`, all doubled. A cell index is a `Nat`,
so Python's negative indexing cannot occur; out of range is `IndexError`. -/
def rect_of (x_edges y_edges : List (Int × Int)) (cell : Nat × Nat) :
    Option (Int × Int × Int × Int) :=
  match x_edges[cell.1]?, y_edges[cell.2]? with
  | some (x0, x1), some (y0, y1) => some (x0, y0, x1 - x0, y1 - y0)
  | _, _ => none

/-- `heatmap_highlight_rects`, `composition.py:235-258`. -/
def heatmap_highlight_rects (snp_idx : Int) (x_coords y_coords : List Int) :
    Option (List (Int × Int × Int × Int)) :=
  match heatmap_highlight_cells snp_idx x_coords.length,               -- l.252
        cell_edges x_coords, cell_edges y_coords with                  -- l.253
  | some cells, some x_edges, some y_edges => allOrNone (rect_of x_edges y_edges) cells
  | _, _, _ => none

/-! ## Properties -/

/-- (a) What a correct highlight of SNP `s` among `n` looks like. -/
def CellsGood (s n : Nat) (cells : List (Nat × Nat)) : Prop :=
  cells.length = n ∧ cells.Nodup ∧ (s, s) ∈ cells ∧
  (∀ c ∈ cells, c.1 < n ∧ c.2 < n ∧ masked c.2 c.1 = false) ∧
  (∀ x, x < n → ∀ y, y < n → ((x, y) ∈ cells ↔ masked y x = false ∧ (x = s ∨ y = s)))

instance (s n : Nat) (cells : List (Nat × Nat)) : Decidable (CellsGood s n cells) := by
  unfold CellsGood; infer_instance

/-- Strictly ascending: every coordinate is below its successor. -/
def Asc (c : List Int) : Prop := ∀ i, i < c.length - 1 → c.getD i 0 < c.getD (i + 1) 0

instance (c : List Int) : Decidable (Asc c) := by unfold Asc; infer_instance

/-- (c) One cell per coordinate; positive width; the centre strictly inside;
each cell's right edge is the next cell's left edge. Edges are doubled. -/
def EdgesGood (c : List Int) (e : List (Int × Int)) : Prop :=
  e.length = c.length ∧
  (∀ i, i < c.length →
    (e.getD i (0, 0)).1 < (e.getD i (0, 0)).2 ∧
    (e.getD i (0, 0)).1 < 2 * c.getD i 0 ∧ 2 * c.getD i 0 < (e.getD i (0, 0)).2) ∧
  (∀ i, i < c.length - 1 → (e.getD i (0, 0)).2 = (e.getD (i + 1) (0, 0)).1)

instance (c : List Int) (e : List (Int × Int)) : Decidable (EdgesGood c e) := by
  unfold EdgesGood; infer_instance

/-- (d) A rectangle for cell `(x, y)`: positive width and height, and it
strictly contains the cell centre `(xc[x], yc[y])` (all doubled). -/
def RectFor (xc yc : List Int) (cell : Nat × Nat) (r : Int × Int × Int × Int) : Prop :=
  0 < r.2.2.1 ∧ 0 < r.2.2.2 ∧
  r.1 < 2 * xc.getD cell.1 0 ∧ 2 * xc.getD cell.1 0 < r.1 + r.2.2.1 ∧
  r.2.1 < 2 * yc.getD cell.2 0 ∧ 2 * yc.getD cell.2 0 < r.2.1 + r.2.2.2

instance (xc yc : List Int) (cell : Nat × Nat) (r : Int × Int × Int × Int) :
    Decidable (RectFor xc yc cell r) := by unfold RectFor; infer_instance

/-- (d) One rectangle per cell, each good for some highlighted cell. -/
def RectsGood (xc yc : List Int) (cells : List (Nat × Nat))
    (rects : List (Int × Int × Int × Int)) : Prop :=
  rects.length = cells.length ∧ ∀ r ∈ rects, ∃ cell ∈ cells, RectFor xc yc cell r

instance (xc yc : List Int) (cells : List (Nat × Nat)) (rects : List (Int × Int × Int × Int)) :
    Decidable (RectsGood xc yc cells rects) := by unfold RectsGood; infer_instance

/-! ## Bounded exhaustive checks -/

def intRange (lo : Int) (count : Nat) : List Int := (List.range count).map (fun (k : Nat) => lo + Int.ofNat k)

/-- (a)+(b): `(snp_idx, n_snps)` in `[-3, 12]²` where a valid input gives a
bad or missing result, or an invalid input gives any result. -/
def badCells : List (Int × Int) :=
  (intRange (-3) 16).flatMap fun s =>
    (intRange (-3) 16).flatMap fun n =>
      let valid : Bool := decide (0 ≤ s ∧ s < n)
      match heatmap_highlight_cells s n with
      | some cells =>
        if valid && decide (CellsGood s.toNat n.toNat cells) then [] else [(s, n)]
      | none => if valid then [(s, n)] else []

#eval badCells   -- []
#guard badCells.isEmpty

def listsOfLen (vals : List Int) : Nat → List (List Int)
  | 0 => [[]]
  | k + 1 => (listsOfLen vals k).flatMap fun l => vals.map (· :: l)

/-- Every list of length `≤ maxLen` over `vals`. -/
def allLists (vals : List Int) (maxLen : Nat) : List (List Int) :=
  (List.range (maxLen + 1)).flatMap (listsOfLen vals)

/-- The coordinate lists searched: length 0..5 over `-2..5` (37449 lists). -/
def searchLists : List (List Int) := allLists (intRange (-2) 8) 5

-- The search is not vacuous: 37449 lists, 218 of them non-empty and ascending.
#guard searchLists.length == 37449
#guard (searchLists.filter fun c => decide (1 ≤ c.length ∧ Asc c)).length == 218

/-- (c): non-empty strictly ascending lists whose edges are missing or bad. -/
def badEdges : List (List Int) :=
  searchLists.filter fun c =>
    decide (1 ≤ c.length ∧ Asc c) &&
      (match cell_edges c with
       | some e => !decide (EdgesGood c e)
       | none => true)

#eval badEdges   -- []
#guard badEdges.isEmpty

/-- (c) converse: lists of length ≥ 2 that are NOT strictly ascending yet
get good edges. Empty means strict ascent is exactly the precondition. -/
def nonAscButGood : List (List Int) :=
  searchLists.filter fun c =>
    decide (2 ≤ c.length ∧ ¬ Asc c) &&
      (match cell_edges c with
       | some e => decide (EdgesGood c e)
       | none => false)

#eval nonAscButGood   -- []
#guard nonAscButGood.isEmpty

/-- (d): `(snp_idx, x_coords)` with ascending `x_coords` of length 1..4 over
`-2..5`, `y_coords = range(n)` as `draw_ld_heatmap` builds them (l.287), and
every `snp_idx` in `[-2, n + 1]`: a valid index must give good rectangles,
an invalid one must raise. -/
def badRects : List (Int × List Int) :=
  (allLists (intRange (-2) 8) 4).flatMap fun xc =>
    if decide (1 ≤ xc.length ∧ Asc xc) then
      let yc := intRange 0 xc.length
      (intRange (-2) (xc.length + 4)).flatMap fun (s : Int) =>
        let valid : Bool := decide (0 ≤ s ∧ s < xc.length)
        match heatmap_highlight_rects s xc yc, heatmap_highlight_cells s xc.length with
        | some rects, some cells =>
          if valid && decide (RectsGood xc yc cells rects) then [] else [(s, xc)]
        | none, _ => if valid then [(s, xc)] else []
        | some _, none => [(s, xc)]
    else []

#eval badRects   -- []
#guard badRects.isEmpty

/-! ### Findings outside the precondition (code as written)

These are counter-examples to `EdgesGood` for inputs that are not strictly
ascending. They are pinned so the build fails if the model stops showing them. -/

-- Empty list: `IndexError` at l.229. Unreachable from a plot call since
-- `prepare_ld_matrix` rejects a `(0, 0)` matrix (`_ld_matrix.py:23-26`).
#guard cell_edges [] == none
-- Duplicate at the low end: first cell has zero width. Python: [(5.0, 5.0), (5.0, 6.0), (6.0, 8.0)].
#guard cell_edges [5, 5, 7] == some [(10, 10), (10, 12), (12, 16)]
-- Interior duplicate: widths stay positive but both centres sit on a shared edge.
#guard cell_edges [1, 5, 5, 9] == some [(-2, 6), (6, 10), (10, 14), (14, 22)]
-- Descending: every cell has low > high (negative width).
#guard cell_edges [3, 2, 1] == some [(7, 5), (5, 3), (3, 1)]
-- Unsorted: the last cell is inverted and cell 1 does not contain its centre 10.
#guard cell_edges [0, 10, 5] == some [(-10, 10), (10, 15), (15, 5)]
-- A zero-width edge becomes a zero-width highlight rectangle.
#guard heatmap_highlight_rects 0 [5, 5, 7] [0, 1, 2]
  == some [(10, -1, 0, 2), (10, 1, 0, 2), (10, 3, 0, 2)]
-- `y_coords` shorter than `x_coords`: `IndexError` at l.256.
#guard heatmap_highlight_rects 0 [1, 2, 3] [0, 1] == none

/-! ## Theorems (all sizes) -/

/-- The cells for a valid input, in closed form. -/
theorem cells_eq (s n : Nat) (h : s < n) :
    heatmap_highlight_cells s n =
      some ((List.range (s + 1)).map (fun j => (j, s))
            ++ (List.range' (s + 1) (n - (s + 1))).map (fun i => (s, i))) := by
  unfold heatmap_highlight_cells
  have hg : ¬ ((n : Int) < 1 ∨ (s : Int) < 0 ∨ (s : Int) ≥ n) := by omega
  simp only [hg, ite_false, Int.toNat_natCast]

/-- (a) For `0 ≤ s < n` the function returns, and its cells are in bounds,
rendered, `n` in number, distinct, include `(s, s)`, and are exactly the
rendered cells in SNP `s`'s row or column. -/
theorem cells_good (s n : Nat) (h : s < n) :
    ∃ cells, heatmap_highlight_cells s n = some cells ∧ CellsGood s n cells := by
  refine ⟨_, cells_eq s n h, ?_, ?_, ?_, ?_, ?_⟩
  · simp only [List.length_append, List.length_map, List.length_range, List.length_range']
    omega
  · rw [List.nodup_append]
    refine ⟨?_, ?_, ?_⟩
    · unfold List.Nodup
      rw [List.pairwise_map]
      exact List.Pairwise.imp (fun hab heq => hab (by simpa using heq)) List.nodup_range
    · unfold List.Nodup
      rw [List.pairwise_map]
      exact List.Pairwise.imp (fun hab heq => hab (by simpa using heq)) List.nodup_range'
    · intro a ha b hb
      simp only [List.mem_map, List.mem_range, List.mem_range'_1] at ha hb
      obtain ⟨j, hj, rfl⟩ := ha
      obtain ⟨i, hi, rfl⟩ := hb
      simp only [ne_eq, Prod.mk.injEq, not_and]
      omega
  · simp only [List.mem_append, List.mem_map, List.mem_range, List.mem_range'_1]
    exact Or.inl ⟨s, by omega, rfl⟩
  · intro c hc
    simp only [List.mem_append, List.mem_map, List.mem_range, List.mem_range'_1] at hc
    simp only [masked, decide_eq_false_iff_not]
    rcases hc with ⟨j, hj, rfl⟩ | ⟨i, hi, rfl⟩ <;> dsimp only <;> omega
  · intro x hx y hy
    simp only [masked, decide_eq_false_iff_not, List.mem_append, List.mem_map, List.mem_range,
      List.mem_range'_1, Prod.mk.injEq]
    constructor
    · rintro (⟨j, hj, rfl, rfl⟩ | ⟨i, hi, rfl, rfl⟩)
      · exact ⟨by omega, Or.inr rfl⟩
      · exact ⟨by omega, Or.inl rfl⟩
    · rintro ⟨hm, hs | hs⟩
      · by_cases hxy : y = s
        · exact Or.inl ⟨x, by omega, rfl, hxy.symm⟩
        · exact Or.inr ⟨y, by omega, hs.symm, rfl⟩
      · exact Or.inl ⟨x, by omega, rfl, hs.symm⟩

/-- (b) The guard: the function returns iff `0 ≤ snp_idx < n_snps`. -/
theorem cells_guard (s n : Int) :
    (heatmap_highlight_cells s n).isSome = true ↔ 0 ≤ s ∧ s < n := by
  unfold heatmap_highlight_cells
  by_cases hg : n < 1 ∨ s < 0 ∨ s ≥ n
  · simp only [hg, ite_true, Option.isSome_none, Bool.false_eq_true, false_iff]; omega
  · simp only [hg, ite_false, Option.isSome_some, true_iff]; omega

theorem mids_cons_cons (a b : Int) (t : List Int) :
    mids (a :: b :: t) = (a + b) :: mids (b :: t) := by
  simp [mids]

/-- `mids` has one entry fewer than `coords`, so `mids[0]` and `mids[-1]`
exist whenever `len(coords) ≥ 2`. -/
theorem mids_length : ∀ c : List Int, (mids c).length = c.length - 1
  | [] => by simp [mids]
  | [_] => by simp [mids]
  | a :: b :: t => by
    rw [mids_cons_cons, List.length_cons, mids_length (b :: t)]
    simp

theorem mids_getD : ∀ (c : List Int) (i : Nat), i + 1 < c.length →
    (mids c).getD i 0 = c.getD i 0 + c.getD (i + 1) 0
  | [], _, h => by simp at h
  | [_], _, h => by simp at h
  | a :: b :: t, 0, _ => by simp [mids_cons_cons]
  | a :: b :: t, i + 1, h => by
    have ih := mids_getD (b :: t) i (by simpa using h)
    rw [mids_cons_cons]
    simpa using ih

theorem getD_zip (l₁ l₂ : List Int) (i : Nat) (h₁ : i < l₁.length) (h₂ : i < l₂.length) :
    (l₁.zip l₂).getD i (0, 0) = (l₁.getD i 0, l₂.getD i 0) := by
  have hz : i < (l₁.zip l₂).length := by simp only [List.length_zip]; omega
  rw [List.getD_eq_getElem?_getD, List.getElem?_eq_getElem hz, List.getElem_zip]
  simp [List.getD_eq_getElem?_getD, h₁, h₂]

theorem getD_append_left (l₁ l₂ : List Int) (i : Nat) (h : i < l₁.length) :
    (l₁ ++ l₂).getD i 0 = l₁.getD i 0 := by
  simp [List.getD_eq_getElem?_getD, List.getElem?_append_left h]

theorem getD_append_right (l₁ l₂ : List Int) (i : Nat) (h : l₁.length ≤ i) :
    (l₁ ++ l₂).getD i 0 = l₂.getD (i - l₁.length) 0 := by
  simp [List.getD_eq_getElem?_getD, List.getElem?_append_right h]

/-- Closed form of every cell for `len(coords) ≥ 2`, ascending or not:
cell `i` runs from the midpoint below to the midpoint above, and the two
outer cells mirror their inner gap. Edges doubled. -/
theorem cell_edges_spec (c : List Int) (h : 2 ≤ c.length) :
    ∃ e, cell_edges c = some e ∧ e.length = c.length ∧
      ∀ i, i < c.length → e.getD i (0, 0) =
        (if i = 0 then 2 * c.getD 0 0 - (c.getD 0 0 + c.getD 1 0 - 2 * c.getD 0 0)
         else c.getD (i - 1) 0 + c.getD i 0,
         if i + 1 < c.length then c.getD i 0 + c.getD (i + 1) 0
         else 2 * c.getD i 0 + (2 * c.getD i 0 - (c.getD (i - 1) 0 + c.getD i 0))) := by
  match c, h with
  | c0 :: c1 :: t, _ =>
    have hml := mids_length (c0 :: c1 :: t)
    simp only [List.length_cons] at hml
    refine ⟨_, rfl, ?_, ?_⟩
    · simp only [List.tail_cons, List.length_zip, List.length_cons, List.length_append,
        List.length_nil, hml]
      omega
    · intro i hi
      simp only [List.length_cons] at hi
      simp only [List.tail_cons]
      rw [getD_zip _ _ _ (by simp only [List.length_cons, List.length_append, hml,
            List.length_nil]; omega)
          (by simp only [List.length_append, hml, List.length_cons, List.length_nil]; omega)]
      have hlast : t.length + 1 + 1 - 1 - 1 + 1 < (c0 :: c1 :: t).length := by
        simp only [List.length_cons]; omega
      have hm0 := mids_getD (c0 :: c1 :: t) 0 (by simp only [List.length_cons]; omega)
      have hmlast := mids_getD (c0 :: c1 :: t) (t.length + 1 + 1 - 1 - 1) hlast
      congr 1
      · cases i with
        | zero => simp only [List.getD_cons_zero, ite_true, hm0]
        | succ k =>
          have hk := mids_getD (c0 :: c1 :: t) k (by simp only [List.length_cons]; omega)
          have : k < (mids (c0 :: c1 :: t)).length := by rw [hml]; omega
          simp only [List.getD_cons_succ, getD_append_left _ _ _ this, hk,
            Nat.add_sub_cancel, Nat.succ_ne_zero, ite_false]
      · by_cases hlt : i + 1 < t.length + 1 + 1
        · have hk := mids_getD (c0 :: c1 :: t) i (by simp only [List.length_cons]; exact hlt)
          have : i < (mids (c0 :: c1 :: t)).length := by rw [hml]; omega
          simp only [getD_append_left _ _ _ this, hk, List.length_cons, hlt, ite_true]
        · have hi' : i = t.length + 1 := by omega
          subst hi'
          have hge : (mids (c0 :: c1 :: t)).length ≤ t.length + 1 := by rw [hml]; omega
          rw [getD_append_right _ _ _ hge]
          simp only [hml, List.length_cons, hlt, ite_false] at hmlast ⊢
          have e1 : t.length + 1 + 1 - 1 - 1 = t.length := by omega
          have e2 : t.length + 1 + 1 - 1 = t.length + 1 := by omega
          rw [e1] at hmlast
          simp only [e2, hmlast, List.getD_cons_zero, Nat.add_sub_cancel, Nat.sub_self]

/-- (c) Strictly ascending coordinates of length ≥ 2 give good edges. -/
theorem cell_edges_good (c : List Int) (h : 2 ≤ c.length) (hasc : Asc c) :
    ∃ e, cell_edges c = some e ∧ EdgesGood c e := by
  obtain ⟨e, he, hlen, hspec⟩ := cell_edges_spec c h
  refine ⟨e, he, hlen, ?_, ?_⟩
  · intro i hi
    rw [hspec i hi]
    by_cases h0 : i = 0
    · subst h0
      have a01 := hasc 0 (by omega)
      have : 0 + 1 < c.length := by omega
      simp only [ite_true, this, Nat.zero_add] at a01 ⊢
      omega
    · have alo := hasc (i - 1) (by omega)
      have e1 : i - 1 + 1 = i := by omega
      rw [e1] at alo
      by_cases hlt : i + 1 < c.length
      · have ahi := hasc i (by omega)
        simp only [h0, hlt, ite_false, ite_true]
        omega
      · simp only [h0, hlt, ite_false]
        omega
  · intro i hi
    have hi' : i + 1 < c.length := by omega
    rw [hspec i (by omega), hspec (i + 1) hi']
    simp only [hi', ite_true, Nat.succ_ne_zero, ite_false, Nat.add_sub_cancel]

/-- (c) converse: for length ≥ 2, good edges force strict ascent. So strict
ascent is exactly the precondition, and any duplicate or descent breaks
`EdgesGood`. -/
theorem cell_edges_good_only_if_asc (c : List Int) (h : 2 ≤ c.length)
    (e : List (Int × Int)) (he : cell_edges c = some e) (hgood : EdgesGood c e) : Asc c := by
  obtain ⟨e', he', _, hspec⟩ := cell_edges_spec c h
  rw [he] at he'
  cases he'
  intro i hi
  have hi' : i + 1 < c.length := by omega
  have hc := (hgood.2.1 (i + 1) hi').2.1
  rw [hspec (i + 1) hi'] at hc
  simp only [Nat.succ_ne_zero, ite_false, Nat.add_sub_cancel] at hc
  omega

/-- (c) the single-coordinate branch (l.226-227): one unit-wide cell. -/
theorem cell_edges_single (a : Int) : EdgesGood [a] [(2 * a - 1, 2 * a + 1)] ∧
    cell_edges [a] = some [(2 * a - 1, 2 * a + 1)] := by
  refine ⟨⟨rfl, ?_, ?_⟩, rfl⟩
  · intro i hi
    have : i = 0 := by simpa using hi
    subst this
    simp only [List.getD_cons_zero]
    omega
  · intro i hi
    simp at hi

/-- (c) for every non-empty strictly ascending list. -/
theorem cell_edges_good' (c : List Int) (h : 1 ≤ c.length) (hasc : Asc c) :
    ∃ e, cell_edges c = some e ∧ EdgesGood c e := by
  match c, h with
  | [a], _ => exact ⟨_, (cell_edges_single a).2, (cell_edges_single a).1⟩
  | a :: b :: t, _ => exact cell_edges_good _ (by simp only [List.length_cons]; omega) hasc

/-- The empty list raises (`IndexError`, l.229). `cell_edges` keeps this
precondition; `prepare_ld_matrix` (`_ld_matrix.py:23-26`) enforces it at intake. -/
theorem cell_edges_nil : cell_edges [] = none := rfl

theorem allOrNone_spec {α β : Type} (f : α → Option β) (P : α → β → Prop) :
    ∀ l : List α, (∀ a ∈ l, ∃ b, f a = some b ∧ P a b) →
      ∃ bs, allOrNone f l = some bs ∧ bs.length = l.length ∧ ∀ b ∈ bs, ∃ a ∈ l, P a b
  | [], _ => ⟨[], rfl, rfl, by intro b hb; cases hb⟩
  | a :: t, hall => by
    obtain ⟨b, hb, hP⟩ := hall a (List.mem_cons_self ..)
    obtain ⟨bs, hbs, hlen, hmem⟩ :=
      allOrNone_spec f P t (fun a' ha' => hall a' (List.mem_cons_of_mem _ ha'))
    refine ⟨b :: bs, ?_, ?_, ?_⟩
    · simp only [allOrNone, hb, hbs]
    · simp only [List.length_cons, hlen]
    · intro b' hb'
      rcases List.mem_cons.mp hb' with rfl | hb'
      · exact ⟨a, List.mem_cons_self .., hP⟩
      · obtain ⟨a', ha', hP'⟩ := hmem b' hb'
        exact ⟨a', List.mem_cons_of_mem _ ha', hP'⟩

theorem getD_eq_of_lt {α : Type} (l : List α) (i : Nat) (d : α) (h : i < l.length) :
    l[i]? = some (l.getD i d) := by
  simp [List.getD_eq_getElem?_getD, h]

/-- (d) With strictly ascending `x_coords` and `y_coords` of equal non-zero
length and `0 ≤ s < n`, every index at l.256 is in bounds (the call returns),
there is one rectangle per cell, and each has positive width and height and
strictly contains its cell's centre. -/
theorem rects_good (s : Nat) (xc yc : List Int) (hs : s < xc.length)
    (hlen : yc.length = xc.length) (hx : Asc xc) (hy : Asc yc) :
    ∃ cells rects, heatmap_highlight_cells s xc.length = some cells ∧
      heatmap_highlight_rects s xc yc = some rects ∧ RectsGood xc yc cells rects := by
  obtain ⟨cells, hcells, _, _, _, hin, _⟩ := cells_good s xc.length hs
  obtain ⟨xe, hxe, hxl, hxg, _⟩ := cell_edges_good' xc (by omega) hx
  obtain ⟨ye, hye, hyl, hyg, _⟩ := cell_edges_good' yc (by omega) hy
  have hall : ∀ cell ∈ cells, ∃ r, rect_of xe ye cell = some r ∧ RectFor xc yc cell r := by
    intro cell hc
    obtain ⟨h1, h2, _⟩ := hin cell hc
    have gx := hxg cell.1 h1
    have gy := hyg cell.2 (by omega)
    refine ⟨((xe.getD cell.1 (0, 0)).1, (ye.getD cell.2 (0, 0)).1,
      (xe.getD cell.1 (0, 0)).2 - (xe.getD cell.1 (0, 0)).1,
      (ye.getD cell.2 (0, 0)).2 - (ye.getD cell.2 (0, 0)).1), ?_, ?_⟩
    · simp only [rect_of, getD_eq_of_lt xe cell.1 (0, 0) (by omega),
        getD_eq_of_lt ye cell.2 (0, 0) (by omega)]
    · simp only [RectFor]
      omega
  obtain ⟨rects, hr, hrl, hrm⟩ := allOrNone_spec (rect_of xe ye) (RectFor xc yc) cells hall
  refine ⟨cells, rects, hcells, ?_, hrl, hrm⟩
  simp only [heatmap_highlight_rects, hcells, hxe, hye, hr]

/-- (d) as `draw_ld_heatmap` calls it (l.287, l.300): `y_coords = range(n)`
is strictly ascending by construction, so only `x_coords` needs a hypothesis. -/
theorem range_asc (n : Nat) : Asc (intRange 0 n) := by
  intro i hi
  have hn : (intRange 0 n).length = n := by simp [intRange]
  rw [hn] at hi
  have h0 : i < n := by omega
  have h1 : i + 1 < n := by omega
  simp [intRange, List.getD_eq_getElem?_getD, h0, h1]
  omega

/-- (b) for the rectangles: an out-of-range `snp_idx` raises. -/
theorem rects_guard (s : Int) (xc yc : List Int) (h : ¬ (0 ≤ s ∧ s < xc.length)) :
    heatmap_highlight_rects s xc yc = none := by
  have hnone : heatmap_highlight_cells s xc.length = none := by
    have := (not_congr (cells_guard s xc.length)).mpr h
    simpa using this
  simp only [heatmap_highlight_rects, hnone]
