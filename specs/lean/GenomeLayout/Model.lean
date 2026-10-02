/-!
# Cumulative chromosome offsets on the shared Manhattan x axis

Source: `src/pylocuszoom/manhattan.py`
* `GenomeLayout.from_frames`, l.114-184: per-chromosome maximum l.150, the
  offset loop l.152-158, cumulative x l.166, `total_length` l.183.
* `_apply_genome_layout`, l.350-367: `map(layout.offsets) + pos` at l.363-365.
* `src/pylocuszoom/panels/miami.py` l.99-109: a highlight region is drawn from
  `offsets[chrom] + start` to `offsets[chrom] + end`.
* `CHROMOSOME_GAP = 1_000_000` (`_plotter_utils.py:47`), user-settable as
  `GenomeWideStyle.chrom_gap` with `ge=0` (`config.py:629-631`).

Model. The input is the list of *present* chromosomes in display order, each
with the positions the pooled frames carry for it. The loop at l.155-158 skips
a chromosome with no rows, so dropping absent chromosomes loses nothing.
Positions, the gap and the offsets are `Int`: Python `int` is unbounded and the
plot-time contract (`schemas.gwas_plot_spec`, l.54-73) puts no lower bound on a
position, so 0 and negative values are representable on purpose.

Outside the model (assumptions, see the report):
* the display order has no duplicate names (a duplicate would overwrite
  `offsets[chrom]` at l.157 and add its length twice at l.158). Enforced for a
  user order by `GenomeWideConfig.validate_custom_chrom_order`
  (`config.py:534-548`), which rejects a name repeated after normalisation;
  the built-in species orders (`species.py:31-64`) carry none;
* positions are integers (l.158 truncates a float maximum with `int()`);
* every x fits int64, the dtype pandas gives `map(offsets) + pos`
  (`totalFrom_le` bounds the largest x by `n * (M + gap)`).
-/

/-! ## Transcription -/

/-- `pooled.groupby("_chrom_str")["_pos"].max()`, `manhattan.py:150`. A group
always has a row; the empty list is given 0 only to make the function total. -/
def max_by_chrom : List Int → Int
  | [] => 0
  | q :: qs => qs.foldl max q

/-- The offset loop, `manhattan.py:152-158`: `offsets[chrom] = cumulative`
then `cumulative += int(max_by_chrom[chrom]) + gap`. -/
def offsetsFrom (gap : Int) : Int → List Int → List Int
  | _, [] => []
  | cumulative, m :: rest => cumulative :: offsetsFrom gap (cumulative + (m + gap)) rest

/-- The value `cumulative` holds when the loop at `manhattan.py:155-158` ends,
returned as `total_length` at l.183. -/
def totalFrom (gap : Int) : Int → List Int → Int
  | cumulative, [] => cumulative
  | cumulative, m :: rest => totalFrom gap (cumulative + (m + gap)) rest

/-- The arithmetic fields of `GenomeLayout`, `manhattan.py:97-102`. `offsets`
is indexed by a present chromosome's rank in display order. -/
structure GenomeLayout where
  offsets : List Int
  total_length : Int
deriving Repr

/-- `GenomeLayout.from_frames`, `manhattan.py:114-184`. `chroms` holds the
pooled positions of each present chromosome, in display order. -/
def from_frames (chroms : List (List Int)) (gap : Int) : GenomeLayout :=
  let maxes := chroms.map max_by_chrom
  { offsets := offsetsFrom gap 0 maxes, total_length := totalFrom gap 0 maxes }

/-- `_cumulative_pos`, `manhattan.py:166` and `:363-365`:
`map(layout.offsets) + pos`. `none` is the NaN pandas gives a chromosome the
layout has no offset for. -/
def cumulative_pos (L : GenomeLayout) (c : Nat) (p : Int) : Option Int :=
  match L.offsets[c]? with
  | some o => some (o + p)
  | none => none

/-- A Miami highlight, `panels/miami.py:101-103`: `(offsets[chrom] + start,
offsets[chrom] + end)`, skipped when the chromosome has no offset (l.108). -/
def highlight (L : GenomeLayout) (c : Nat) (start stop : Int) : Option (Int × Int) :=
  match L.offsets[c]? with
  | some o => some (o + start, o + stop)
  | none => none

/-! ## Properties, as decidable `Prop`s -/

/-- Two points `(i, p)` and `(j, q)` of one genome. -/
structure Case where
  gap : Int
  chroms : List (List Int)
  i : Nat
  p : Int
  j : Nat
  q : Int
deriving Repr

def Case.layout (k : Case) : GenomeLayout := from_frames k.chroms k.gap
def Case.xp (k : Case) : Int := (cumulative_pos k.layout k.i k.p).getD 0
def Case.xq (k : Case) : Int := (cumulative_pos k.layout k.j k.q).getD 0

/-- Both points get an x (no NaN). -/
def Defined (k : Case) : Prop :=
  (cumulative_pos k.layout k.i k.p).isSome = true ∧ (cumulative_pos k.layout k.j k.q).isSome = true

/-- (a) Every point of an earlier chromosome is strictly left of every point
of a later one. -/
def OrderOK (k : Case) : Prop := k.i < k.j → k.xp < k.xq

/-- (b) Distinct `(chrom, pos)` get distinct x. -/
def InjectiveOK (k : Case) : Prop := (k.i ≠ k.j ∨ k.p ≠ k.q) → k.xp ≠ k.xq

/-- (c) Points of different chromosomes are more than `gap` apart. -/
def SeparationOK (k : Case) : Prop := k.i < k.j → k.xp + k.gap < k.xq

/-- (d) `total_length` is one gap past the last chromosome's end, one gap or
more past every point, and the first chromosome starts at 0. -/
def TotalOK (k : Case) : Prop :=
  k.xp + k.gap ≤ k.layout.total_length ∧
  k.layout.total_length
    = k.layout.offsets.getLast?.getD 0 + (k.chroms.map max_by_chrom).getLast?.getD 0 + k.gap ∧
  k.layout.offsets.head? = some 0

instance (k : Case) : Decidable (Defined k) := by unfold Defined; infer_instance
instance (k : Case) : Decidable (OrderOK k) := by unfold OrderOK; infer_instance
instance (k : Case) : Decidable (InjectiveOK k) := by unfold InjectiveOK; infer_instance
instance (k : Case) : Decidable (SeparationOK k) := by unfold SeparationOK; infer_instance
instance (k : Case) : Decidable (TotalOK k) := by unfold TotalOK; infer_instance

/-- A highlight region `(i, start, stop)` on one genome. -/
structure Region where
  gap : Int
  chroms : List (List Int)
  i : Nat
  start : Int
  stop : Int
deriving Repr

def Region.layout (r : Region) : GenomeLayout := from_frames r.chroms r.gap
def Region.span (r : Region) : Int × Int := (highlight r.layout r.i r.start r.stop).getD (0, 0)
def Region.offset (r : Region) (j : Nat) : Int := r.layout.offsets.getD j 0
def Region.maxPos (r : Region) (j : Nat) : Int := (r.chroms.map max_by_chrom).getD j 0

/-- (e) The span stays inside its own chromosome's extent
`[offset + 1, offset + max_pos]` and misses every other chromosome's extent. -/
def HighlightOK (r : Region) : Prop :=
  (highlight r.layout r.i r.start r.stop).isSome = true ∧
  r.offset r.i + 1 ≤ r.span.1 ∧ r.span.2 ≤ r.offset r.i + r.maxPos r.i ∧
  ∀ j, j < r.chroms.length → j ≠ r.i →
    r.span.2 < r.offset j + 1 ∨ r.offset j + r.maxPos j < r.span.1

instance (r : Region) : Decidable (HighlightOK r) := by unfold HighlightOK; infer_instance

/-! ## Bounded exhaustive check -/

def intRange (lo hi : Int) : List Int :=
  (List.range (hi - lo + 1).toNat).map fun (k : Nat) => lo + Int.ofNat k

/-- Every list over `alpha` of length at most `n`. -/
def listsUpTo {α : Type} (alpha : List α) : Nat → List (List α)
  | 0 => [[]]
  | n + 1 => [] :: alpha.flatMap fun a => (listsUpTo alpha n).map (a :: ·)

/-- Every genome of 1..`nChrom` chromosomes, each with 1..`nPos` positions in
`[posLo, posHi]`. -/
def genomes (posLo posHi : Int) (nChrom nPos : Nat) : List (List (List Int)) :=
  let chromSets := (listsUpTo (intRange posLo posHi) nPos).filter (!·.isEmpty)
  (listsUpTo chromSets nChrom).filter (!·.isEmpty)

def points (chroms : List (List Int)) : List (Nat × Int) :=
  (List.range chroms.length).flatMap fun i => (chroms.getD i []).map fun p => (i, p)

/-- Every pair of data points of every genome, for every gap in
`[gapLo, gapHi]`, that fails `ok`. -/
def badPairs (ok : Case → Bool) (gapLo gapHi posLo posHi : Int) (nChrom nPos : Nat) : List Case :=
  (intRange gapLo gapHi).flatMap fun gap =>
    (genomes posLo posHi nChrom nPos).flatMap fun chroms =>
      (points chroms).flatMap fun (i, p) =>
        (points chroms).filterMap fun (j, q) =>
          let k : Case := ⟨gap, chroms, i, p, j, q⟩
          if ok k then none else some k

/-- Every highlight `(i, start, stop)` with `1 ≤ start ≤ stop ≤ stopHi`, and
`stop ≤ max_pos[i]` when `clip`, that fails `HighlightOK`. -/
def badRegions (clip : Bool) (gapLo gapHi posLo posHi stopHi : Int) (nChrom nPos : Nat) :
    List Region :=
  (intRange gapLo gapHi).flatMap fun gap =>
    (genomes posLo posHi nChrom nPos).flatMap fun chroms =>
      (List.range chroms.length).flatMap fun i =>
        (intRange 1 stopHi).flatMap fun start =>
          (intRange start stopHi).filterMap fun stop =>
            let r : Region := ⟨gap, chroms, i, start, stop⟩
            if clip && decide (r.maxPos i < stop) then none
            else if HighlightOK r then none else some r

def allOK (k : Case) : Bool :=
  decide (Defined k) && decide (OrderOK k) && decide (InjectiveOK k) &&
  decide (SeparationOK k) && decide (TotalOK k)

-- Under the hypotheses of the theorems (gap ≥ 0, positions ≥ 1): up to 3
-- chromosomes, 2 positions each in 1..3, gap 0..3. No counter-example.
#eval (badPairs (fun _ => false) 0 3 1 3 3 2).length   -- cases examined
#eval (badPairs allOK 0 3 1 3 3 2).length   -- 0
#guard (badPairs allOK 0 3 1 3 3 2).isEmpty

-- 0-based positions with a positive gap are also safe for (a), (b) and (d).
#guard (badPairs (fun k => decide (OrderOK k) && decide (InjectiveOK k) && decide (TotalOK k))
  1 3 0 2 3 2).isEmpty

-- COUNTER-EXAMPLE 1 (positions ≥ 1 dropped): gap = 0 with a position 0.
-- chroms = [[1], [0]]: x(chrom 0, pos 1) = 1 = x(chrom 1, pos 0).
#eval (badPairs (fun k => decide (InjectiveOK k)) 0 0 0 1 2 1).take 2
#guard !(badPairs (fun k => decide (InjectiveOK k)) 0 0 0 1 2 1).isEmpty
#guard !(badPairs (fun k => decide (OrderOK k)) 0 0 0 1 2 1).isEmpty
-- (c) is the only property a position 0 breaks under a positive gap: the
-- distance is exactly `gap`, not more.
#guard !(badPairs (fun k => decide (SeparationOK k)) 1 1 0 1 2 1).isEmpty

-- COUNTER-EXAMPLE 2 (gap ≥ 0 dropped): the model accepts a negative gap, the
-- code does not (`chrom_gap` has `ge=0`); order breaks.
#guard !(badPairs (fun k => decide (OrderOK k)) (-2) (-2) 1 2 2 1).isEmpty

-- The real constant, `CHROMOSOME_GAP = 1_000_000`, on human-sized chromosomes.
#guard allOK ⟨1000000, [[10583, 248946422], [10019, 242183529], [9, 198235559]], 0, 248946422, 1, 10019⟩
#guard allOK ⟨1000000, [[10583, 248946422], [10019, 242183529], [9, 198235559]], 1, 242183529, 2, 9⟩

-- Highlights clipped to the data (`stop ≤ max_pos[i]`): no counter-example.
#eval (badRegions true 0 2 1 3 6 3 1).length   -- 0
#guard (badRegions true 0 2 1 3 6 3 1).isEmpty

-- COUNTER-EXAMPLE 3 (unclipped, as `miami.py:101-103` is written): a region
-- wider than the data leaves its chromosome, and past `max_pos + gap` it
-- covers the next chromosome's points.
#eval ((badRegions false 1 1 1 2 4 2 1).filter (·.chroms.length = 2)).take 2
-- gap = 1, chroms = [[1], [1]], region (0, 1, 3): span 1..3 covers x = 3, the
-- point of chromosome 1.
#guard !decide (HighlightOK ⟨1, [[1], [1]], 0, 1, 3⟩)
#guard !(badRegions false 0 2 1 3 6 3 1).isEmpty

/-! ## Lemmas about the loop -/

theorem init_le_foldl_max (ps : List Int) (init : Int) : init ≤ ps.foldl max init := by
  induction ps generalizing init with
  | nil => exact Int.le_refl _
  | cons q qs ih =>
    simp only [List.foldl_cons]
    have h1 := ih (max init q)
    have h2 := Int.le_max_left init q
    omega

theorem le_foldl_max (ps : List Int) (init p : Int) (hp : p ∈ ps) : p ≤ ps.foldl max init := by
  induction ps generalizing init with
  | nil => cases hp
  | cons q qs ih =>
    simp only [List.foldl_cons]
    rcases List.mem_cons.mp hp with h | h
    · subst h
      have h1 := init_le_foldl_max qs (max init p)
      have h2 := Int.le_max_right init p
      omega
    · exact ih _ h

/-- `max_by_chrom` is an upper bound of the chromosome's positions. -/
theorem le_max_by_chrom (ps : List Int) (p : Int) (hp : p ∈ ps) : p ≤ max_by_chrom ps := by
  cases ps with
  | nil => cases hp
  | cons q qs =>
    simp only [max_by_chrom]
    rcases List.mem_cons.mp hp with h | h
    · subst h; exact init_le_foldl_max qs p
    · exact le_foldl_max qs q p h

theorem max_by_chrom_nonneg (ps : List Int) (h : ∀ p ∈ ps, 0 ≤ p) : 0 ≤ max_by_chrom ps := by
  cases ps with
  | nil => simp [max_by_chrom]
  | cons q qs =>
    have h1 := le_max_by_chrom (q :: qs) q (by simp)
    have h2 := h q (by simp)
    omega

theorem offsetsFrom_length (gap c : Int) (maxes : List Int) :
    (offsetsFrom gap c maxes).length = maxes.length := by
  induction maxes generalizing c with
  | nil => simp [offsetsFrom]
  | cons m rest ih => simp [offsetsFrom, ih]

/-- No offset is left of the running total it started from. -/
theorem offsetsFrom_ge (gap c : Int) (maxes : List Int) (h : ∀ m ∈ maxes, 0 ≤ m + gap)
    (j : Nat) (o : Int) (ho : (offsetsFrom gap c maxes)[j]? = some o) : c ≤ o := by
  induction maxes generalizing c j with
  | nil => simp [offsetsFrom] at ho
  | cons m rest ih =>
    cases j with
    | zero =>
      simp [offsetsFrom] at ho
      omega
    | succ j =>
      simp only [offsetsFrom, List.getElem?_cons_succ] at ho
      have h1 := ih (c + (m + gap)) (fun x hx => h x (by simp [hx])) j ho
      have h2 := h m (by simp)
      omega

/-- The next chromosome starts exactly `max_pos + gap` after this one. -/
theorem offsets_succ (gap c : Int) (maxes : List Int) (i : Nat) (oi oj mi : Int)
    (hi : (offsetsFrom gap c maxes)[i]? = some oi)
    (hj : (offsetsFrom gap c maxes)[i + 1]? = some oj)
    (hm : maxes[i]? = some mi) : oj = oi + mi + gap := by
  induction maxes generalizing c i with
  | nil => simp at hm
  | cons m rest ih =>
    cases i with
    | zero =>
      cases rest with
      | nil => simp [offsetsFrom] at hj
      | cons m2 rest2 =>
        simp [offsetsFrom] at hi hj hm
        omega
    | succ i =>
      simp only [offsetsFrom, List.getElem?_cons_succ] at hi hj hm
      exact ih (c + (m + gap)) i hi hj hm

/-- A later chromosome starts at least `max_pos + gap` after an earlier one. -/
theorem offsets_step (gap c : Int) (maxes : List Int) (h : ∀ m ∈ maxes, 0 ≤ m + gap)
    (i j : Nat) (hij : i < j) (oi oj mi : Int)
    (hi : (offsetsFrom gap c maxes)[i]? = some oi)
    (hj : (offsetsFrom gap c maxes)[j]? = some oj)
    (hm : maxes[i]? = some mi) : oi + mi + gap ≤ oj := by
  induction maxes generalizing c i j with
  | nil => simp at hm
  | cons m rest ih =>
    cases j with
    | zero => omega
    | succ j =>
      have hrest : ∀ x ∈ rest, 0 ≤ x + gap := fun x hx => h x (by simp [hx])
      simp only [offsetsFrom, List.getElem?_cons_succ] at hj
      cases i with
      | zero =>
        simp [offsetsFrom] at hi hm
        have h1 := offsetsFrom_ge gap (c + (m + gap)) rest hrest j oj hj
        omega
      | succ i =>
        simp only [offsetsFrom, List.getElem?_cons_succ] at hi hm
        exact ih (c + (m + gap)) hrest i j (by omega) hi hj hm

theorem totalFrom_ge (gap c : Int) (maxes : List Int) (h : ∀ m ∈ maxes, 0 ≤ m + gap) :
    c ≤ totalFrom gap c maxes := by
  induction maxes generalizing c with
  | nil => simp [totalFrom]
  | cons m rest ih =>
    simp only [totalFrom]
    have h1 := ih (c + (m + gap)) (fun x hx => h x (by simp [hx]))
    have h2 := h m (by simp)
    omega

/-- `total_length` is at least one gap past every chromosome's end. -/
theorem totalFrom_ge_end (gap c : Int) (maxes : List Int) (h : ∀ m ∈ maxes, 0 ≤ m + gap)
    (i : Nat) (oi mi : Int)
    (hi : (offsetsFrom gap c maxes)[i]? = some oi) (hm : maxes[i]? = some mi) :
    oi + mi + gap ≤ totalFrom gap c maxes := by
  induction maxes generalizing c i with
  | nil => simp at hm
  | cons m rest ih =>
    have hrest : ∀ x ∈ rest, 0 ≤ x + gap := fun x hx => h x (by simp [hx])
    cases i with
    | zero =>
      simp [offsetsFrom] at hi hm
      simp only [totalFrom]
      have h1 := totalFrom_ge gap (c + (m + gap)) rest hrest
      omega
    | succ i =>
      simp only [offsetsFrom, List.getElem?_cons_succ] at hi hm
      simp only [totalFrom]
      exact ih (c + (m + gap)) hrest i hi hm

/-- `total_length` is exactly one gap past the last chromosome's end. -/
theorem totalFrom_eq_last (gap c : Int) (maxes : List Int) (i : Nat) (oi mi : Int)
    (hlast : i + 1 = maxes.length)
    (hi : (offsetsFrom gap c maxes)[i]? = some oi) (hm : maxes[i]? = some mi) :
    totalFrom gap c maxes = oi + mi + gap := by
  induction maxes generalizing c i with
  | nil => simp at hm
  | cons m rest ih =>
    cases i with
    | zero =>
      cases rest with
      | nil =>
        simp [offsetsFrom] at hi hm
        simp only [totalFrom]
        omega
      | cons m2 rest2 => simp at hlast
    | succ i =>
      simp only [offsetsFrom, List.getElem?_cons_succ] at hi hm
      simp only [totalFrom]
      exact ih (c + (m + gap)) i (by simp at hlast; omega) hi hm

/-- Width bound: `n` chromosomes of at most `M` bases end before
`n * (M + gap)`. Compare with `2^63 - 1` for the int64 x column. -/
theorem totalFrom_le (gap c M : Int) (maxes : List Int) (h : ∀ m ∈ maxes, m ≤ M) :
    totalFrom gap c maxes ≤ c + (maxes.length : Int) * (M + gap) := by
  induction maxes generalizing c with
  | nil => simp [totalFrom]
  | cons m rest ih =>
    simp only [totalFrom]
    have h1 := ih (c + (m + gap)) (fun x hx => h x (by simp [hx]))
    have h2 := h m (by simp)
    have h3 : (((m :: rest).length : Nat) : Int) * (M + gap)
        = (rest.length : Int) * (M + gap) + (M + gap) := by
      simp only [List.length_cons, Int.natCast_add, Int.add_mul]
      simp
    omega

/-- 1000 chromosomes of 150 Gbp each (the largest known genome is about
150 Gbp in total) with the default gap stay far inside int64. -/
example : (1000 : Int) * (150000000000 + 1000000) < 9223372036854775807 := by decide

/-! ## Glue between the loop lemmas and `from_frames` -/

theorem cumulative_pos_some {L : GenomeLayout} {c : Nat} {p x : Int}
    (h : cumulative_pos L c p = some x) : ∃ o, L.offsets[c]? = some o ∧ x = o + p := by
  unfold cumulative_pos at h
  split at h
  · rename_i o ho
    exact ⟨o, ho, by cases h; rfl⟩
  · cases h

theorem highlight_some {L : GenomeLayout} {c : Nat} {s e xs xe : Int}
    (h : highlight L c s e = some (xs, xe)) :
    ∃ o, L.offsets[c]? = some o ∧ xs = o + s ∧ xe = o + e := by
  unfold highlight at h
  split at h
  · rename_i o ho
    exact ⟨o, ho, by cases h; exact ⟨rfl, rfl⟩⟩
  · cases h

theorem maxes_get {chroms : List (List Int)} {i : Nat} {ps : List Int}
    (h : chroms[i]? = some ps) : (chroms.map max_by_chrom)[i]? = some (max_by_chrom ps) := by
  simp [h]

theorem maxes_nonneg (chroms : List (List Int)) (gap : Int) (hgap : 0 ≤ gap)
    (hpos : ∀ ps ∈ chroms, ∀ p ∈ ps, 0 ≤ p) :
    ∀ m ∈ chroms.map max_by_chrom, 0 ≤ m + gap := by
  intro m hm
  obtain ⟨ps, hps, rfl⟩ := List.mem_map.mp hm
  have := max_by_chrom_nonneg ps (hpos ps hps)
  omega

/-! ## Theorems, for every genome size

Named hypotheses:
* `hgap`  : `0 ≤ gap`. Enforced by `GenomeWideStyle.chrom_gap` (`ge=0`).
* `hpos`  : every plotted position is `≥ 0`. NOT enforced at plot time.
* `hone`  : `1 ≤ gap + q` for every plotted position `q`, that is positions
            `≥ 1` (1-based coordinates), or positions `≥ 0` with `gap ≥ 1`.
* `hci`/`hp` : the point belongs to a frame the layout was built from, which
            `prepare_manhattan_frames` (l.331-347) guarantees.
-/

/-- Core of (a), (b), (c): a point of a later chromosome is at least
`gap + q` right of any point of an earlier one. -/
theorem separation (chroms : List (List Int)) (gap : Int) (hgap : 0 ≤ gap)
    (hpos : ∀ ps ∈ chroms, ∀ p ∈ ps, 0 ≤ p)
    (i j : Nat) (hij : i < j) (pi : List Int) (hci : chroms[i]? = some pi)
    (p q : Int) (hp : p ∈ pi) (xp xq : Int)
    (hxp : cumulative_pos (from_frames chroms gap) i p = some xp)
    (hxq : cumulative_pos (from_frames chroms gap) j q = some xq) :
    xp + gap + q ≤ xq := by
  obtain ⟨oi, hoi, rfl⟩ := cumulative_pos_some hxp
  obtain ⟨oj, hoj, rfl⟩ := cumulative_pos_some hxq
  simp only [from_frames] at hoi hoj
  have h1 := offsets_step gap 0 _ (maxes_nonneg chroms gap hgap hpos) i j hij oi oj _
    hoi hoj (maxes_get hci)
  have h2 := le_max_by_chrom pi p hp
  omega

/-- (a) Order: every point of an earlier chromosome is strictly left of every
point of a later one. -/
theorem order_strict (chroms : List (List Int)) (gap : Int) (hgap : 0 ≤ gap)
    (hpos : ∀ ps ∈ chroms, ∀ p ∈ ps, 0 ≤ p)
    (i j : Nat) (hij : i < j) (pi : List Int) (hci : chroms[i]? = some pi)
    (p q : Int) (hp : p ∈ pi) (hone : 1 ≤ gap + q) (xp xq : Int)
    (hxp : cumulative_pos (from_frames chroms gap) i p = some xp)
    (hxq : cumulative_pos (from_frames chroms gap) j q = some xq) :
    xp < xq := by
  have := separation chroms gap hgap hpos i j hij pi hci p q hp xp xq hxp hxq
  omega

/-- (c) Separation with 1-based positions: points of different chromosomes
are more than `gap` apart (`gap + 1` is attained). -/
theorem separation_gap (chroms : List (List Int)) (gap : Int) (hgap : 0 ≤ gap)
    (hpos : ∀ ps ∈ chroms, ∀ p ∈ ps, 0 ≤ p)
    (i j : Nat) (hij : i < j) (pi : List Int) (hci : chroms[i]? = some pi)
    (p q : Int) (hp : p ∈ pi) (hq : 1 ≤ q) (xp xq : Int)
    (hxp : cumulative_pos (from_frames chroms gap) i p = some xp)
    (hxq : cumulative_pos (from_frames chroms gap) j q = some xq) :
    xp + gap < xq := by
  have := separation chroms gap hgap hpos i j hij pi hci p q hp xp xq hxp hxq
  omega

/-- (b) Injectivity: distinct `(chrom, pos)` get distinct x. -/
theorem injective (chroms : List (List Int)) (gap : Int) (hgap : 0 ≤ gap)
    (hpos : ∀ ps ∈ chroms, ∀ p ∈ ps, 0 ≤ p)
    (hone : ∀ ps ∈ chroms, ∀ p ∈ ps, 1 ≤ gap + p)
    (i j : Nat) (pi pj : List Int) (hci : chroms[i]? = some pi) (hcj : chroms[j]? = some pj)
    (p q : Int) (hp : p ∈ pi) (hq : q ∈ pj) (xp xq : Int)
    (hxp : cumulative_pos (from_frames chroms gap) i p = some xp)
    (hxq : cumulative_pos (from_frames chroms gap) j q = some xq)
    (hne : i ≠ j ∨ p ≠ q) : xp ≠ xq := by
  rcases Nat.lt_trichotomy i j with hlt | heq | hgt
  · have := order_strict chroms gap hgap hpos i j hlt pi hci p q hp
      (hone pj (List.mem_of_getElem? hcj) q hq) xp xq hxp hxq
    omega
  · subst heq
    obtain ⟨oi, hoi, rfl⟩ := cumulative_pos_some hxp
    obtain ⟨oj, hoj, rfl⟩ := cumulative_pos_some hxq
    rw [hoi] at hoj
    cases hoj
    rcases hne with h | h
    · exact absurd rfl h
    · omega
  · have := order_strict chroms gap hgap hpos j i hgt pj hcj q p hq
      (hone pi (List.mem_of_getElem? hci) p hp) xq xp hxq hxp
    omega

/-- Every chromosome the layout was built from has an offset, so no point of
a frame the layout was built from gets NaN. -/
theorem defined (chroms : List (List Int)) (gap : Int) (i : Nat) (p : Int)
    (hi : i < chroms.length) :
    (cumulative_pos (from_frames chroms gap) i p).isSome = true := by
  have hlen : i < (offsetsFrom gap 0 (chroms.map max_by_chrom)).length := by
    rw [offsetsFrom_length]; simpa using hi
  simp [cumulative_pos, from_frames, List.getElem?_eq_getElem hlen]

/-- (d) `total_length` is at least one gap right of every point. -/
theorem total_ge_x (chroms : List (List Int)) (gap : Int) (hgap : 0 ≤ gap)
    (hpos : ∀ ps ∈ chroms, ∀ p ∈ ps, 0 ≤ p)
    (i : Nat) (pi : List Int) (hci : chroms[i]? = some pi) (p : Int) (hp : p ∈ pi) (xp : Int)
    (hxp : cumulative_pos (from_frames chroms gap) i p = some xp) :
    xp + gap ≤ (from_frames chroms gap).total_length := by
  obtain ⟨oi, hoi, rfl⟩ := cumulative_pos_some hxp
  simp only [from_frames] at hoi ⊢
  have h1 := totalFrom_ge_end gap 0 _ (maxes_nonneg chroms gap hgap hpos) i oi _
    hoi (maxes_get hci)
  have h2 := le_max_by_chrom pi p hp
  omega

/-- (d) `total_length` is the last chromosome's end plus one gap. Needs no
sign hypothesis. -/
theorem total_eq_last_end (chroms : List (List Int)) (gap : Int)
    (i : Nat) (hlast : i + 1 = chroms.length) (pi : List Int) (hci : chroms[i]? = some pi)
    (oi : Int) (hoi : (from_frames chroms gap).offsets[i]? = some oi) :
    (from_frames chroms gap).total_length = oi + max_by_chrom pi + gap := by
  simp only [from_frames] at hoi ⊢
  exact totalFrom_eq_last gap 0 _ i oi _ (by simpa using hlast) hoi (maxes_get hci)

/-- The first present chromosome starts at x = 0, so its points are drawn at
their own position. -/
theorem first_offset_zero (ps : List Int) (rest : List (List Int)) (gap : Int) :
    (from_frames (ps :: rest) gap).offsets[0]? = some 0 := by
  simp [from_frames, offsetsFrom]

/-- (e) A highlight with `1 ≤ start` and `stop ≤ max_pos[i]` lies inside
chromosome `i`'s extent. -/
theorem highlight_within (chroms : List (List Int)) (gap : Int)
    (i : Nat) (pi : List Int) (start stop : Int)
    (hs : 1 ≤ start) (he : stop ≤ max_by_chrom pi) (oi xs xe : Int)
    (hoi : (from_frames chroms gap).offsets[i]? = some oi)
    (hh : highlight (from_frames chroms gap) i start stop = some (xs, xe)) :
    oi + 1 ≤ xs ∧ xe ≤ oi + max_by_chrom pi := by
  obtain ⟨o, ho, rfl, rfl⟩ := highlight_some hh
  rw [hoi] at ho
  cases ho
  omega

/-- (e) Such a highlight misses every other chromosome's extent
`[offset[j] + 1, offset[j] + max_pos[j]]`. -/
theorem highlight_clear (chroms : List (List Int)) (gap : Int) (hgap : 0 ≤ gap)
    (hpos : ∀ ps ∈ chroms, ∀ p ∈ ps, 0 ≤ p)
    (i j : Nat) (hne : j ≠ i) (pi pj : List Int)
    (hci : chroms[i]? = some pi) (hcj : chroms[j]? = some pj)
    (start stop : Int) (hs : 1 ≤ start) (he : stop ≤ max_by_chrom pi) (oj xs xe : Int)
    (hoj : (from_frames chroms gap).offsets[j]? = some oj)
    (hh : highlight (from_frames chroms gap) i start stop = some (xs, xe)) :
    xe < oj + 1 ∨ oj + max_by_chrom pj < xs := by
  obtain ⟨oi, hoi, rfl, rfl⟩ := highlight_some hh
  simp only [from_frames] at hoi hoj
  have hm := maxes_nonneg chroms gap hgap hpos
  rcases Nat.lt_trichotomy i j with hlt | heq | hgt
  · have := offsets_step gap 0 _ hm i j hlt oi oj _ hoi hoj (maxes_get hci)
    left; omega
  · exact absurd heq.symm hne
  · have := offsets_step gap 0 _ hm j i hgt oj oi _ hoj hoi (maxes_get hcj)
    right; omega

/-- (e), the failure, as `panels/miami.py:101-103` is written: the highlight
reaches the next chromosome's first base exactly when `stop > max_pos + gap`.
Nothing in the code bounds `stop`. -/
theorem highlight_bleeds_iff (chroms : List (List Int)) (gap : Int)
    (i : Nat) (pi : List Int) (hci : chroms[i]? = some pi)
    (start stop : Int) (onext xs xe : Int)
    (hnext : (from_frames chroms gap).offsets[i + 1]? = some onext)
    (hh : highlight (from_frames chroms gap) i start stop = some (xs, xe)) :
    onext + 1 ≤ xe ↔ max_by_chrom pi + gap < stop := by
  obtain ⟨oi, hoi, rfl, rfl⟩ := highlight_some hh
  simp only [from_frames] at hoi hnext
  have := offsets_succ gap 0 _ i oi onext _ hoi hnext (maxes_get hci)
  omega
