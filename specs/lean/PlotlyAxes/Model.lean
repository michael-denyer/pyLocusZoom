/-!
# Plotly subplot index, axis names and pixel splitting

Model of the integer arithmetic in

* `src/pylocuszoom/backends/plotly_layout.py:15-78` (`_Panel.subplot_idx`,
  `axis`, `ref`, `secondary_ref`) and `:102-113` (`secondary_axis_key`);
* `src/pylocuszoom/backends/plotly_backend.py:100-106` and `:137-151` (which
  `row`, `col`, `n_cols` values `create_figure` / `create_figure_grid` build)
  and `:507-536` (`create_twin_axis`, where `secondary_ref` becomes a layout
  key unless a subplot already has that key);
* `src/pylocuszoom/backends/_coerce.py:49-63` (`split_pixels`).

Modelling choices:

* Python `int` is unbounded, so `row`, `col`, `n_cols` and `total` are `Int`.
  Nothing in the code bounds them, and `Int` lets the model state what
  `col = 0` or `col > n_cols` does.
* An axis name is its numeric suffix: `none` is the unsuffixed name (`"yaxis"`,
  `"y"`), `some k` is `f"{kind}{k}"`. Two names of the same kind are equal as
  Python strings exactly when their suffixes are equal. That step is assumed,
  not proved: it rests on `str(int)` being injective and never empty.
* `secondary_axis_key` only swaps the prefix `"y"` for `"yaxis"`, so a
  secondary reference and a primary layout key collide exactly when their
  suffixes are equal.
* A figure is its subplot count `N`: its layout holds the primary keys
  `yaxis`, `yaxis2`, ..., `yaxisN` and nothing else that `create_twin_axis`
  could mistake for one. The guard's `overlaying is None` clause, which tells
  a subplot's own axis from a secondary axis written by an earlier call, is
  assumed to do that and is not modelled.
* `split_pixels` has a float arm. The even arm (`total // n`) is modelled
  exactly: `Int` `/` is floor division for a positive divisor, as `//` is.
  The ratio arm is modelled only for non-negative integer ratios and a
  non-negative total, where `int(total * r / denominator)` is
  `total * r / denominator` in `Nat`; float rounding is outside the model.
* `none` stands for the Python exception (`ZeroDivisionError` in
  `split_pixels`, `ValidationError` in `create_twin_axis`).
-/

namespace PlotlyAxes

/-! ## Transcription -/

/-- `_Panel.subplot_idx`, `plotly_layout.py:31`:
`(self.row - 1) * self.n_cols + self.col`. -/
def subplot_idx (row col n_cols : Int) : Int := (row - 1) * n_cols + col

/-- `_Panel.axis`, `plotly_layout.py:42-43`: `f"{kind}{idx}" if idx > 1 else kind`.
`none` is the unsuffixed layout key. -/
def axis (idx : Int) : Option Int := if idx > 1 then some idx else none

/-- `_Panel.ref`, `plotly_layout.py:67-68`: the same rule as `axis`. -/
def ref (idx : Int) : Option Int := if idx > 1 then some idx else none

/-- `_Panel.secondary_ref`, `plotly_layout.py:56`: `f"y{100 + self.subplot_idx - 1}"`.
Always suffixed. `secondary_axis_key` (`plotly_layout.py:111-112`) keeps the
suffix, so this is also the suffix of the layout key `create_twin_axis` writes
(`plotly_backend.py:514-525`). -/
def secondary_ref (idx : Int) : Option Int := some (100 + idx - 1)

/-- `PlotlyBackend.create_twin_axis`, `plotly_backend.py:514-536`, on a figure
of `N` subplots: the suffix of the secondary axis it creates. `none` is the
`ValidationError` of `plotly_backend.py:517-522`, raised when the key is in
the layout as a subplot's own y-axis. The secondary key is always suffixed, so
it is one of `yaxis2 .. yaxisN` exactly when its suffix is in `2 .. N`. -/
def twin_axis (N idx : Int) : Option Int :=
  if 1 < 100 + idx - 1 ∧ 100 + idx - 1 ≤ N then none else some (100 + idx - 1)

/-- The cell is inside an `n_rows` by `n_cols` grid, 1-based. These are the
only cells `create_figure` (`plotly_backend.py:100-106`, `n_cols = 1`) and
`create_figure_grid` (`plotly_backend.py:137-151`) build. -/
def InGrid (n_rows n_cols row col : Int) : Prop :=
  1 ≤ row ∧ row ≤ n_rows ∧ 1 ≤ col ∧ col ≤ n_cols

instance (n_rows n_cols row col : Int) : Decidable (InGrid n_rows n_cols row col) := by
  unfold InGrid; infer_instance

/-- Even arm of `split_pixels`, `_coerce.py:60-61`: `[total // n] * n`.
`none` is the `ZeroDivisionError` raised when `n == 0`. -/
def split_pixels_even (total : Int) (n : Nat) : Option (List Int) :=
  if n = 0 then none else some (List.replicate n (total / (n : Int)))

/-- Ratio arm of `split_pixels`, `_coerce.py:62-63`, for non-negative integer
ratios: `[int(total * r / denominator) for r in ratios]`. An empty list
divides nothing and returns `[]`; otherwise a zero denominator raises. -/
def split_pixels_ratio (total : Nat) (ratios : List Nat) : Option (List Nat) :=
  if ratios = [] then some []
  else if ratios.sum = 0 then none
  else some (ratios.map fun r => total * r / ratios.sum)

/-! ## Bounded exhaustive checks -/

/-- Cells in the order `create_figure_grid` returns them
(`plotly_backend.py:147-151`): row-major, 1-based. -/
def cells (n_rows n_cols : Nat) : List (Int × Int) :=
  (List.range n_rows).flatMap fun (r : Nat) =>
    (List.range n_cols).map fun (c : Nat) => ((r : Int) + 1, (c : Int) + 1)

def idxOf (n_cols : Nat) (p : Int × Int) : Int := subplot_idx p.1 p.2 n_cols

/-- Grids up to `maxR` by `maxC` whose indices are not exactly `1, 2, ..., N`
in row-major order. Empty means `subplot_idx` is a bijection on each grid. -/
def badGrids (maxR maxC : Nat) : List (Nat × Nat) :=
  (List.range (maxR + 1)).flatMap fun r =>
    (List.range (maxC + 1)).filterMap fun c =>
      if (cells r c).map (idxOf c) = (List.range (r * c)).map (fun (k : Nat) => (k : Int) + 1)
      then none else some (r, c)

/-- Ordered pairs of distinct cells of one grid that satisfy `hit`. -/
def cellPairs (n_rows n_cols : Nat) (hit : Int → Int → Bool) :
    List ((Int × Int) × (Int × Int)) :=
  let cs := cells n_rows n_cols
  cs.flatMap fun p => cs.filterMap fun q =>
    if hit (idxOf n_cols p) (idxOf n_cols q) then some (p, q) else none

/-- Distinct panels with the same primary axis name. -/
def primaryDups (n_rows n_cols : Nat) : List ((Int × Int) × (Int × Int)) :=
  (cellPairs n_rows n_cols fun i j => axis i == axis j).filter fun pq => pq.1 != pq.2

/-- Distinct panels with the same secondary axis name. -/
def secondaryDups (n_rows n_cols : Nat) : List ((Int × Int) × (Int × Int)) :=
  (cellPairs n_rows n_cols fun i j => secondary_ref i == secondary_ref j).filter
    fun pq => pq.1 != pq.2

/-- `(p, q)`: the secondary axis of panel `p` has the name of the primary
axis of panel `q`. Worst case: every panel is given a secondary axis. -/
def crossHits (n_rows n_cols : Nat) : List ((Int × Int) × (Int × Int)) :=
  cellPairs n_rows n_cols fun i j => secondary_ref i == axis j

/-- `(p, q)`: `create_twin_axis` accepts panel `p` and the axis it creates
has the name of the primary axis of panel `q`. Empty means the guard holds. -/
def acceptedHits (n_rows n_cols : Nat) : List ((Int × Int) × (Int × Int)) :=
  cellPairs n_rows n_cols fun i j =>
    match twin_axis ((n_rows * n_cols : Nat) : Int) i with
    | none => false
    | some s => some s == axis j

/-- Panels of one grid that `create_twin_axis` refuses. -/
def refusedCells (n_rows n_cols : Nat) : List (Int × Int) :=
  (cells n_rows n_cols).filter fun p =>
    (twin_axis ((n_rows * n_cols : Nat) : Int) (idxOf n_cols p)).isNone

/-- Grids up to `maxR` by `maxC` with any axis-name collision. -/
def collidingGrids (maxR maxC : Nat) : List (Nat × Nat) :=
  (List.range (maxR + 1)).flatMap fun r =>
    (List.range (maxC + 1)).filterMap fun c =>
      if (primaryDups r c).isEmpty && (secondaryDups r c).isEmpty && (crossHits r c).isEmpty
      then none else some (r, c)

-- Bound: 0..26 rows by 0..5 columns, so up to 130 subplots, past the offset 100.
#eval badGrids 26 5   -- []
#guard (badGrids 26 5).isEmpty

-- No grid in the bound has two panels sharing a primary name or a secondary name.
#guard ((List.range 27).all fun r => (List.range 6).all fun c =>
  (primaryDups r c).isEmpty && (secondaryDups r c).isEmpty)

-- Every colliding grid has at least 100 subplots, and every grid with at
-- least 100 subplots collides: the exact threshold in the bound is N = 100.
#eval (collidingGrids 26 5).map fun rc => (rc, rc.1 * rc.2)
#guard (collidingGrids 26 5).all fun rc => rc.1 * rc.2 ≥ 100
#guard ((List.range 27).all fun r => (List.range 6).all fun c =>
  r * c < 100 || (collidingGrids 26 5).contains (r, c))

-- Smallest colliding configurations: 100 subplots, panel 1's secondary axis
-- `y100` against the primary axis of subplot 100.
#eval crossHits 100 1   -- [((1, 1), (100, 1))]
#eval crossHits 50 2    -- [((1, 1), (50, 2))]
#eval crossHits 20 5    -- [((1, 1), (20, 5))]
#guard crossHits 100 1 = [((1, 1), (100, 1))]
#guard crossHits 50 2 = [((1, 1), (50, 2))]
#guard (crossHits 99 1).isEmpty
#guard (crossHits 33 3).isEmpty
#guard secondary_ref 1 = some 100

-- The guard refuses exactly the panels whose secondary name is taken, and
-- leaves every other panel its `secondary_ref` name.
#eval refusedCells 100 1   -- [(1, 1)]
#guard twin_axis 100 1 = none
#guard twin_axis 100 2 = some 101
#guard twin_axis 99 1 = some 100
#guard refusedCells 100 1 = [(1, 1)]
#guard refusedCells 50 2 = [(1, 1)]
#guard (refusedCells 26 5).length = 31
#guard ((List.range 27).all fun r => (List.range 6).all fun c =>
  (acceptedHits r c).isEmpty &&
    refusedCells r c == (crossHits r c).map Prod.fst)

-- Outside the grid two panels share an index (`col = n_cols + 1`, `col = 0`).
#eval (subplot_idx 1 3 2, subplot_idx 2 1 2)   -- (3, 3)
#eval (subplot_idx 2 0 2, subplot_idx 1 2 2)   -- (2, 2)
#guard subplot_idx 1 3 2 = subplot_idx 2 1 2
#guard subplot_idx 2 0 2 = subplot_idx 1 2 2
-- An index of 0 gets the unsuffixed name of subplot 1.
#guard axis (subplot_idx 1 0 2) = axis (subplot_idx 1 1 2)

/-- The property of the even arm, for one input. -/
def EvenSplitOk (total : Int) (n : Nat) : Prop :=
  match split_pixels_even total n with
  | none => n = 0
  | some parts =>
    parts.length = n ∧ (0 ≤ total → ∀ p ∈ parts, 0 ≤ p) ∧
      parts.sum ≤ total ∧ total - parts.sum < n

instance (total : Int) (n : Nat) : Decidable (EvenSplitOk total n) := by
  unfold EvenSplitOk; split <;> infer_instance

/-- Inputs `(total, n)` with `-8 ≤ total ≤ maxT`, `n ≤ maxN` that break the
even-arm property. -/
def badEvenSplits (maxT maxN : Nat) : List (Int × Nat) :=
  (List.range (maxT + 9)).flatMap fun (t : Nat) =>
    (List.range (maxN + 1)).filterMap fun n =>
      let total : Int := (t : Int) - 8
      if EvenSplitOk total n then none else some (total, n)

#eval badEvenSplits 64 12   -- []
#guard (badEvenSplits 64 12).isEmpty
#guard split_pixels_even 800 0 = none
#guard split_pixels_even 800 3 = some [266, 266, 266]

/-- Ratio lists of length at most 3 with entries at most 4, and totals at most
`maxT`, where the ratio arm returns parts that sum above `total` or fall short
by the number of parts or more. Bounded check only for the shortfall. -/
def badRatioSplits (maxT : Nat) : List (Nat × List Nat) :=
  let rs := List.range 5
  let lists : List (List Nat) :=
    [[]] ++ rs.map (fun a => [a]) ++ rs.flatMap (fun a => rs.map fun b => [a, b]) ++
      rs.flatMap (fun a => rs.flatMap fun b => rs.map fun c => [a, b, c])
  (List.range (maxT + 1)).flatMap fun total =>
    lists.filterMap fun ratios =>
      match split_pixels_ratio total ratios with
      | none => if ratios ≠ [] ∧ ratios.sum = 0 then none else some (total, ratios)
      | some parts =>
        if parts.sum ≤ total ∧ (ratios = [] ∨ total - parts.sum < ratios.length)
        then none else some (total, ratios)

#eval badRatioSplits 40   -- []
#guard (badRatioSplits 40).isEmpty
#guard split_pixels_ratio 800 [3, 1] = some [600, 200]
#guard split_pixels_ratio 800 [0, 0] = none
#guard split_pixels_ratio 800 [] = some []

/-! ## (a) `subplot_idx` is a bijection from the grid onto `1 .. n_rows * n_cols` -/

private theorem row_step (n a b : Int) (hn : 0 ≤ n) (h : a < b) : a * n + n ≤ b * n := by
  have h1 : (a + 1) * n ≤ b * n := Int.mul_le_mul_of_nonneg_right (by omega) hn
  rw [Int.add_mul, Int.one_mul] at h1
  exact h1

/-- Injective on the grid. Only the column bounds are used. -/
theorem subplot_idx_injective (n_cols row₁ col₁ row₂ col₂ : Int)
    (h₁ : 1 ≤ col₁ ∧ col₁ ≤ n_cols) (h₂ : 1 ≤ col₂ ∧ col₂ ≤ n_cols)
    (h : subplot_idx row₁ col₁ n_cols = subplot_idx row₂ col₂ n_cols) :
    row₁ = row₂ ∧ col₁ = col₂ := by
  unfold subplot_idx at h
  have hn : 0 ≤ n_cols := by omega
  have hlt := row_step n_cols (row₁ - 1) (row₂ - 1) hn
  have hgt := row_step n_cols (row₂ - 1) (row₁ - 1) hn
  have hrow : row₁ = row₂ := by omega
  subst hrow
  omega

/-- In range on the grid. -/
theorem subplot_idx_in_range (n_rows n_cols row col : Int)
    (h : InGrid n_rows n_cols row col) :
    1 ≤ subplot_idx row col n_cols ∧ subplot_idx row col n_cols ≤ n_rows * n_cols := by
  obtain ⟨hr1, hr2, hc1, hc2⟩ := h
  unfold subplot_idx
  have hn : 0 ≤ n_cols := by omega
  have h0 : 0 ≤ (row - 1) * n_cols := Int.mul_nonneg (by omega) hn
  have h1 : (row - 1 + 1) * n_cols ≤ n_rows * n_cols :=
    Int.mul_le_mul_of_nonneg_right (by omega) hn
  rw [Int.add_mul, Int.one_mul] at h1
  omega

/-- Onto `1 .. n_rows * n_cols`. -/
theorem subplot_idx_surjective (n_rows n_cols k : Int)
    (hk : 1 ≤ k ∧ k ≤ n_rows * n_cols) (hc : 0 < n_cols) :
    ∃ row col, InGrid n_rows n_cols row col ∧ subplot_idx row col n_cols = k := by
  refine ⟨(k - 1) / n_cols + 1, (k - 1) % n_cols + 1, ?_, ?_⟩
  · have hq0 : 0 ≤ (k - 1) / n_cols := Int.ediv_nonneg (by omega) (by omega)
    have hq1 : (k - 1) / n_cols < n_rows := Int.ediv_lt_of_lt_mul hc (by omega)
    have hm0 : 0 ≤ (k - 1) % n_cols := Int.emod_nonneg _ (by omega)
    have hm1 : (k - 1) % n_cols < n_cols := Int.emod_lt_of_pos _ hc
    unfold InGrid
    omega
  · unfold subplot_idx
    have h := Int.mul_ediv_add_emod (k - 1) n_cols
    have hcomm : ((k - 1) / n_cols + 1 - 1) * n_cols = n_cols * ((k - 1) / n_cols) := by
      rw [Int.add_sub_cancel, Int.mul_comm]
    omega

/-- `col = n_cols + 1` lands on the first cell of the next row. -/
theorem subplot_idx_col_overflow (row n_cols : Int) :
    subplot_idx row (n_cols + 1) n_cols = subplot_idx (row + 1) 1 n_cols := by
  unfold subplot_idx
  have h : (row + 1 - 1) * n_cols = (row - 1) * n_cols + n_cols := by
    have : row + 1 - 1 = row - 1 + 1 := by omega
    rw [this, Int.add_mul, Int.one_mul]
  omega

/-- `col = 0` lands on the last cell of the previous row. -/
theorem subplot_idx_col_zero (row n_cols : Int) :
    subplot_idx (row + 1) 0 n_cols = subplot_idx row n_cols n_cols := by
  unfold subplot_idx
  have h : (row + 1 - 1) * n_cols = (row - 1) * n_cols + n_cols := by
    have : row + 1 - 1 = row - 1 + 1 := by omega
    rw [this, Int.add_mul, Int.one_mul]
  omega

/-! ## (b) primary axis names are distinct for distinct panels -/

/-- For indices from 1 up, equal names mean equal indices: the unsuffixed name
of index 1 is no other index's name. -/
theorem axis_injective (i j : Int) (hi : 1 ≤ i) (hj : 1 ≤ j) (h : axis i = axis j) :
    i = j := by
  unfold axis at h
  split at h <;> split at h <;> simp at h <;> omega

theorem ref_eq_axis (i : Int) : ref i = axis i := rfl

/-- Below 1 the rule does collide: every index at most 1 is unsuffixed. -/
theorem axis_unsuffixed (i : Int) (hi : i ≤ 1) : axis i = none := by
  unfold axis
  split
  · omega
  · rfl

/-- Two grid panels with the same primary axis name are the same panel. -/
theorem primary_names_distinct (n_rows n_cols row₁ col₁ row₂ col₂ : Int)
    (h₁ : InGrid n_rows n_cols row₁ col₁) (h₂ : InGrid n_rows n_cols row₂ col₂)
    (h : axis (subplot_idx row₁ col₁ n_cols) = axis (subplot_idx row₂ col₂ n_cols)) :
    row₁ = row₂ ∧ col₁ = col₂ := by
  have r₁ := subplot_idx_in_range _ _ _ _ h₁
  have r₂ := subplot_idx_in_range _ _ _ _ h₂
  have hidx := axis_injective _ _ r₁.1 r₂.1 h
  exact subplot_idx_injective n_cols row₁ col₁ row₂ col₂ ⟨h₁.2.2.1, h₁.2.2.2⟩
    ⟨h₂.2.2.1, h₂.2.2.2⟩ hidx

/-! ## (c) secondary axis names -/

/-- Two panels never share a secondary axis name, for any indices. -/
theorem secondary_ref_injective (i j : Int) (h : secondary_ref i = secondary_ref j) :
    i = j := by
  unfold secondary_ref at h
  simp at h
  omega

/-- The exact collision condition: panel `i`'s secondary axis is panel `j`'s
primary axis exactly when `j = i + 99`. -/
theorem secondary_hits_primary_iff (i j : Int) (hi : 1 ≤ i) :
    secondary_ref i = axis j ↔ j = i + 99 := by
  unfold secondary_ref axis
  constructor
  · intro h
    split at h <;> simp at h
    omega
  · intro h
    subst h
    have : i + 99 > 1 := by omega
    simp [this]
    omega

/-- With at most 99 subplots, no secondary axis has a primary axis's name. -/
theorem no_collision_upto_99 (N i j : Int) (hN : N ≤ 99)
    (hi : 1 ≤ i ∧ i ≤ N) (hj : 1 ≤ j ∧ j ≤ N) : secondary_ref i ≠ axis j := by
  intro h
  have := (secondary_hits_primary_iff i j hi.1).mp h
  omega

/-- The first panel's secondary axis is `y100`, which is subplot 100's
primary axis. -/
theorem collision_at_100 : secondary_ref 1 = axis 100 := by decide

/-- 99 is the largest subplot count that is collision-free when any panel may
carry a secondary axis. -/
theorem collision_free_iff (N : Int) :
    (∀ i j, 1 ≤ i ∧ i ≤ N → 1 ≤ j ∧ j ≤ N → secondary_ref i ≠ axis j) ↔ N ≤ 99 := by
  constructor
  · intro h
    apply Int.not_lt.mp
    intro hN
    exact h 1 100 (by omega) (by omega) collision_at_100
  · intro hN i j hi hj
    exact no_collision_upto_99 N i j hN hi hj

/-- The same bound for a grid: the count that matters is `n_rows * n_cols`,
not the number of rows. -/
theorem grid_no_collision (n_rows n_cols row₁ col₁ row₂ col₂ : Int)
    (hN : n_rows * n_cols ≤ 99)
    (h₁ : InGrid n_rows n_cols row₁ col₁) (h₂ : InGrid n_rows n_cols row₂ col₂) :
    secondary_ref (subplot_idx row₁ col₁ n_cols) ≠ axis (subplot_idx row₂ col₂ n_cols) :=
  no_collision_upto_99 (n_rows * n_cols) _ _ hN
    (subplot_idx_in_range _ _ _ _ h₁) (subplot_idx_in_range _ _ _ _ h₂)

/-- A 50 by 2 grid has only 50 rows and still collides. -/
theorem grid_50x2_collides :
    InGrid 50 2 1 1 ∧ InGrid 50 2 50 2 ∧
      secondary_ref (subplot_idx 1 1 2) = axis (subplot_idx 50 2 2) := by decide

/-! ## (d) the `create_twin_axis` guard -/

/-- The guard refuses exactly when some subplot of the figure has the name. -/
theorem twin_axis_refuses_iff (N i : Int) :
    twin_axis N i = none ↔ ∃ j, 1 ≤ j ∧ j ≤ N ∧ secondary_ref i = axis j := by
  unfold twin_axis secondary_ref axis
  constructor
  · intro h
    split at h
    · refine ⟨100 + i - 1, by omega, by omega, ?_⟩
      have : 100 + i - 1 > 1 := by omega
      simp [this]
    · simp at h
  · intro ⟨j, hj1, hjN, h⟩
    split at h
    · simp at h
      have : 1 < 100 + i - 1 ∧ 100 + i - 1 ≤ N := by omega
      simp [this]
    · simp at h

/-- An accepted secondary axis keeps the `secondary_ref` name. -/
theorem twin_axis_accepted_name (N i s : Int) (h : twin_axis N i = some s) :
    secondary_ref i = some s := by
  unfold twin_axis at h
  unfold secondary_ref
  split at h
  · simp at h
  · exact h

/-- No accepted secondary axis has the name of a subplot's primary axis, for
any figure size and any index. -/
theorem twin_axis_accepted_no_collision (N i s j : Int) (h : twin_axis N i = some s)
    (hj : 1 ≤ j ∧ j ≤ N) : some s ≠ axis j := by
  unfold twin_axis at h
  unfold axis
  split at h
  · simp at h
  · simp at h
    intro heq
    split at heq
    · simp at heq
      omega
    · simp at heq

/-- With at most 99 subplots the guard refuses nothing, so the names are the
ones `secondary_ref` gave before the guard existed. -/
theorem twin_axis_unchanged_upto_99 (N i : Int) (hN : N ≤ 99) (hi : 1 ≤ i) :
    twin_axis N i = secondary_ref i := by
  unfold twin_axis secondary_ref
  have : ¬(1 < 100 + i - 1 ∧ 100 + i - 1 ≤ N) := by omega
  simp [this]

/-- The same for a grid of any size: an accepted secondary axis is no grid
cell's primary axis. -/
theorem grid_twin_axis_no_collision (n_rows n_cols row₁ col₁ row₂ col₂ s : Int)
    (h : twin_axis (n_rows * n_cols) (subplot_idx row₁ col₁ n_cols) = some s)
    (h₂ : InGrid n_rows n_cols row₂ col₂) :
    some s ≠ axis (subplot_idx row₂ col₂ n_cols) :=
  twin_axis_accepted_no_collision _ _ _ _ h (subplot_idx_in_range _ _ _ _ h₂)

/-! ## (e) `split_pixels` -/

private theorem sum_replicate (n : Nat) (x : Int) :
    (List.replicate n x).sum = (n : Int) * x := by
  induction n with
  | zero => simp
  | succ k ih =>
    rw [List.replicate_succ, List.sum_cons, ih]
    have : ((k + 1 : Nat) : Int) * x = (k : Int) * x + x := by
      rw [Int.natCast_succ, Int.add_mul, Int.one_mul]
    omega

/-- `n = 0` raises (`total // 0`). -/
theorem split_pixels_even_zero (total : Int) : split_pixels_even total 0 = none := rfl

/-- For `n > 0` the even arm returns `n` parts that sum to at most `total`,
with a shortfall below `n`, and no part is negative when `total` is not. -/
theorem split_pixels_even_ok (total : Int) (n : Nat) (hn : 0 < n) :
    ∃ parts, split_pixels_even total n = some parts ∧ parts.length = n ∧
      (0 ≤ total → ∀ p ∈ parts, 0 ≤ p) ∧ parts.sum ≤ total ∧ total - parts.sum < n := by
  refine ⟨List.replicate n (total / (n : Int)), ?_, ?_, ?_, ?_, ?_⟩
  · unfold split_pixels_even
    have : n ≠ 0 := by omega
    simp [this]
  · simp
  · intro ht p hp
    rw [List.eq_of_mem_replicate hp]
    exact Int.ediv_nonneg ht (by omega)
  · rw [sum_replicate]
    have h := Int.mul_ediv_add_emod total (n : Int)
    have hm : 0 ≤ total % (n : Int) := Int.emod_nonneg _ (by omega)
    omega
  · rw [sum_replicate]
    have h := Int.mul_ediv_add_emod total (n : Int)
    have hm : total % (n : Int) < (n : Int) := Int.emod_lt_of_pos _ (by omega)
    omega

private theorem ratio_sum_mul_le (total d : Nat) (ratios : List Nat) :
    (ratios.map fun r => total * r / d).sum * d ≤ total * ratios.sum := by
  induction ratios with
  | nil => simp
  | cons r rs ih =>
    rw [List.map_cons, List.sum_cons, List.sum_cons, Nat.add_mul, Nat.mul_add]
    exact Nat.add_le_add (Nat.div_mul_le_self _ _) ih

/-- The ratio arm, for non-negative integer ratios with a positive sum,
returns parts that sum to at most `total`. -/
theorem split_pixels_ratio_sum_le (total : Nat) (ratios parts : List Nat)
    (h : split_pixels_ratio total ratios = some parts) : parts.sum ≤ total := by
  unfold split_pixels_ratio at h
  split at h
  · cases h; simp
  · split at h
    · cases h
    · cases h
      have hd : 0 < ratios.sum := by omega
      exact Nat.le_of_mul_le_mul_right (ratio_sum_mul_le total ratios.sum ratios) hd

end PlotlyAxes
