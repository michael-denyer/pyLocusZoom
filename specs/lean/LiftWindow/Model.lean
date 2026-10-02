/-!
# Lifted-window arithmetic and 0/1-based conversion in liftover

Source: `src/pylocuszoom/_liftover.py`
* `_lift_one`        l.170-181  (`pos - 1` in, `+ 1` out)
* `liftover_region`  l.235-302  (`start`/`end` = min/max of lifted positions, l.295-296;
                                 the lead is lifted on its own, l.291-293)
* `lift_window`      l.342-450  (window arithmetic l.415-432, aggregation l.435;
                                 lead checks l.416-429 and l.436-443)
* `filter_by_region` `src/pylocuszoom/utils.py` l.208 (inclusive bounds)

Callers: `LocusZoomPlotter._lift`, `src/pylocuszoom/plotter.py` l.288-331, reached from
`plot()` (l.272) and `plot_stacked()` (l.405). `RegionConfig`
(`src/pylocuszoom/config.py` l.104-122) guarantees `1 ≤ start < end`;
`_require_lead_in_region` (config.py l.405-411) guarantees `start ≤ lead_pos ≤ end`.

Types: Python ints are unbounded, so every position is an `Int`. Positions are
1-based; the lifter speaks 0-based on both sides. A frame is the list of its
position column, already restricted to the region's chromosome; the other
columns do not enter the arithmetic.

What is proved for every size (no bound): `lift_window_ordered` (a),
`lift_window_contains` (b), `window_start_margin` / `window_end_margin` (c),
`liftOne_shift` / `liftOne_ge_one` (d), `lift_window_lead_inside` (e: every
lead `lift_window` returns lies in the returned window, row of its frame or
not), `lead_in_frame_inside_window` (a lead that is a row lifts inside the
window, so the check of l.436-443 does not drop it) and `stray_lead_off_rows`
(a lead that is no kept row never shares a kept row's lifted position). What
is a concrete example: `lead_not_a_row_is_dropped` and
`stray_lead_on_a_row_is_dropped` (the two leads the checks drop) and
`negative_lifter_escapes` (why (b) needs `NonnegLifter`). Everything else in
this file is a bounded search.
-/

/-! ## Transcription -/

/-- `_Outcome`, `_liftover.py` l.163-167. -/
inductive Outcome where
  | lifted
  | unmapped
  | multiMapped
  | crossChrom
deriving Repr, DecidableEq

/-- `CoordinateLifter.convert_coordinate` for one fixed chromosome
(`_liftover.py` l.39-45): a 0-based query gives a list of hits, each
`(hit is on the query chromosome, 0-based target position)`. Python's `None`
(chromosome unknown) and `[]` are both falsy at l.175, so both are `[]`. -/
abbrev Lifter := Int → List (Bool × Int)

/-- `_lift_one`, `_liftover.py` l.170-181. -/
def liftOne (lifter : Lifter) (pos : Int) : Outcome × Option Int :=
  match lifter (pos - 1) with                       -- l.174
  | [] => (.unmapped, none)                         -- l.175-176
  | [(same, p)] =>
    if same then (.lifted, some (p + 1))            -- l.181
    else (.crossChrom, none)                        -- l.179-180
  | _ :: _ :: _ => (.multiMapped, none)             -- l.177-178

/-- `filter_by_region`, `utils.py` l.208: `(pos >= start) & (pos <= end)`. -/
def filterByRegion (start end_ : Int) (frame : List Int) : List Int :=
  frame.filter fun p => decide (start ≤ p) && decide (p ≤ end_)

/-- The kept rows of `liftover_region`, `_liftover.py` l.273-280, as
`(source position, lifted position)` pairs in input order. -/
def liftoverRegion (lifter : Lifter) (selected : List Int) : List (Int × Int) :=
  selected.filterMap fun p => (liftOne lifter p).2.map fun q => (p, q)

/-- `min` of a non-empty column, head first (`.min()`, l.295, l.430, l.435). -/
def minOf (a : Int) (l : List Int) : Int := l.foldl min a

/-- `max` of a non-empty column, head first (`.max()`, l.296, l.432, l.435). -/
def maxOf (a : Int) (l : List Int) : Int := l.foldl max a

/-- `_liftover.py` l.430: `max(1, lift.start - int(source_pos.min() - start))`. -/
def windowStart (start srcMin liftStart : Int) : Int :=
  max 1 (liftStart - (srcMin - start))

/-- `_liftover.py` l.432:
`max(lift.end + int(end - source_pos.max()), window_start + 1)`. -/
def windowEnd (end_ srcMax liftEnd ws : Int) : Int :=
  max (liftEnd + (end_ - srcMax)) (ws + 1)

/-- One iteration of the loop in `lift_window`, `_liftover.py` l.392-432, from
the kept rows. `none` is the `ValidationError` of l.392-396. `source_pos`
(l.415) is the source column of the kept rows; `lift.start` / `lift.end`
(l.295-296) are the min / max of their lifted column. -/
def frameWindow (start end_ : Int) : List (Int × Int) → Option (Int × Int)
  | [] => none
  | (p, q) :: t =>
    let ws := windowStart start (minOf p (t.map Prod.fst)) (minOf q (t.map Prod.snd))
    some (ws, windowEnd end_ (maxOf p (t.map Prod.fst)) (maxOf q (t.map Prod.snd)) ws)

/-- Filter, lift and window one frame: `_liftover.py` l.381-432. -/
def panelWindow (lifter : Lifter) (start end_ : Int) (frame : List Int) : Option (Int × Int) :=
  frameWindow start end_ (liftoverRegion lifter (filterByRegion start end_ frame))

/-- The `starts` / `ends` lists of `lift_window`, l.379-434; `none` as soon as
one frame raises. -/
def frameWindows (lifter : Lifter) (start end_ : Int) : List (List Int) → Option (List (Int × Int))
  | [] => some []
  | f :: fs =>
    match panelWindow lifter start end_ f, frameWindows lifter start end_ fs with
    | some w, some ws => some (w :: ws)
    | _, _ => none

/-- `lift_window`, `_liftover.py` l.342-450, returning `(start, end)`:
`(min(starts), max(ends))`, l.435. An empty frame list is `none`, since
`min([])` raises. -/
def liftWindow (lifter : Lifter) (start end_ : Int) (frames : List (List Int)) : Option (Int × Int) :=
  match frameWindows lifter start end_ frames with
  | some (w :: ws) => some (minOf w.1 (ws.map Prod.fst), maxOf w.2 (ws.map Prod.snd))
  | _ => none

/-- The lifted lead of `liftover_region`, `_liftover.py` l.291-293: lifted by
its own `_lift_one` call, whether or not it is a row of the frame. -/
def liftLead (lifter : Lifter) (lead : Option Int) : Option Int :=
  lead.bind fun p => (liftOne lifter p).2

/-- The lead one frame hands on, `_liftover.py` l.416-429: the lifted lead,
dropped when its source position is no kept row of the frame yet its lifted
position is a kept row's. `kept` is the `(source, lifted)` pairs of l.415. -/
def frameLead (kept : List (Int × Int)) (lifter : Lifter) (p : Int) : Option Int :=
  match (liftOne lifter p).2 with
  | some q =>
    if p ∉ kept.map Prod.fst ∧ q ∈ kept.map Prod.snd then none else some q
  | none => none

/-- `_liftover.py` l.436-443: a lead outside the final window is dropped. -/
def keepInside (w : Int × Int) : Option Int → Option Int
  | some q => if w.1 ≤ q ∧ q ≤ w.2 then some q else none
  | none => none

/-- `lift_window`, `_liftover.py` l.342-450, returning
`((start, end), lead_positions)`. Frames and leads are paired by `zip`
(l.380). -/
def liftWindowLeads (lifter : Lifter) (start end_ : Int) (frames : List (List Int))
    (leads : List (Option Int)) : Option ((Int × Int) × List (Option Int)) :=
  match liftWindow lifter start end_ frames with
  | some w =>
    some (w, (List.zipWith
      (fun f lead =>
        lead.bind (frameLead (liftoverRegion lifter (filterByRegion start end_ f)) lifter))
      frames leads).map (keepInside w))
  | none => none

/-! ## Properties -/

/-- (a) The window is a valid `RegionConfig`: `1 ≤ start < end`. -/
def WindowOrdered (w : Int × Int) : Prop := 1 ≤ w.1 ∧ w.1 < w.2

instance (w : Int × Int) : Decidable (WindowOrdered w) := by
  unfold WindowOrdered; infer_instance

/-- (b), (e) A position lies in the closed window. -/
def Inside (w : Int × Int) (q : Int) : Prop := w.1 ≤ q ∧ q ≤ w.2

instance (w : Int × Int) (q : Int) : Decidable (Inside w q) := by
  unfold Inside; infer_instance

/-- Hypothesis on the real lifter, not on the model: every hit is a 0-based
coordinate, so it is `≥ 0`. pyliftover returns chain coordinates, which are
non-negative; `InMemoryLifter` returns whatever its dict holds. -/
def NonnegLifter (lifter : Lifter) : Prop := ∀ x, ∀ h ∈ lifter x, 0 ≤ h.2

/-! ## Bounded exhaustive checks -/

def sublists : List Int → List (List Int)
  | [] => [[]]
  | x :: xs => let r := sublists xs; r ++ r.map (x :: ·)

/-- Every function from `n` consecutive 0-based positions to `cod`. -/
def tables (cod : List (List (Bool × Int))) : Nat → List (List (List (Bool × Int)))
  | 0 => [[]]
  | n + 1 => cod.flatMap fun c => (tables cod n).map (c :: ·)

def tableLifter (tbl : List (List (Bool × Int))) : Lifter :=
  fun x => if 0 ≤ x then tbl.getD x.toNat [] else []

/-- Hit lists a lifter may return: unmapped, three single hits, a hit on
another chromosome, an ambiguous double hit. -/
def fullCod : List (List (Bool × Int)) :=
  [[], [(true, 0)], [(true, 2)], [(true, 6)], [(false, 1)], [(true, 1), (true, 3)]]

def smallCod : List (List (Bool × Int)) := [[], [(true, 0)], [(true, 2)], [(true, 6)]]

/-- Regions `1 ≤ start < end ≤ n + 1`, as `RegionConfig` allows. -/
def regions (n : Nat) : List (Int × Int) :=
  (List.range n).flatMap fun a =>
    (List.range (n + 1)).filterMap fun b =>
      if a < b then some (Int.ofNat a + 1, Int.ofNat b + 1) else none

def positions (n : Nat) : List Int := (List.range n).map fun i => Int.ofNat i + 1

/-- (a) and (b) on one input: the window is ordered and holds every selected
row that lifted. `true` when `lift_window` raises. -/
def windowHolds (lifter : Lifter) (start end_ : Int) (frames : List (List Int)) : Bool :=
  match liftWindow lifter start end_ frames with
  | none => true
  | some w =>
    decide (WindowOrdered w) &&
      frames.all fun f =>
        f.all fun p =>
          -- The inclusive region test is restated here, not taken from
          -- `filterByRegion`, so a wrong filter bound is a counter-example.
          if start ≤ p ∧ p ≤ end_ then
            match (liftOne lifter p).2 with
            | some q => decide (Inside w q)
            | none => true
          else true

/-- Counter-examples to (a)/(b) with one frame: positions `1..n`, every lifter
table over `cod`, every region. -/
def badWindows1 (cod : List (List (Bool × Int))) (n : Nat) :
    List (List (List (Bool × Int)) × (Int × Int) × List Int) :=
  (tables cod n).flatMap fun tbl =>
    (regions n).flatMap fun r =>
      (sublists (positions n)).filterMap fun f =>
        if windowHolds (tableLifter tbl) r.1 r.2 [f] then none else some (tbl, r, f)

/-- Counter-examples to (a)/(b) with two frames, exercising the `min(starts)`
/ `max(ends)` aggregation. -/
def badWindows2 (cod : List (List (Bool × Int))) (n : Nat) :
    List (List (List (Bool × Int)) × (Int × Int) × List Int × List Int) :=
  (tables cod n).flatMap fun tbl =>
    (regions n).flatMap fun r =>
      (sublists (positions n)).flatMap fun f =>
        (sublists (positions n)).filterMap fun g =>
          if windowHolds (tableLifter tbl) r.1 r.2 [f, g] then none else some (tbl, r, f, g)

#eval (badWindows1 fullCod 4).length   -- 0
#guard (badWindows1 fullCod 4).isEmpty
#eval (badWindows2 smallCod 4).length  -- 0
#guard (badWindows2 smallCod 4).isEmpty

/-- (c) on one frame: where a clamp does not bind, the margin beside the
outermost lifted row equals the source margin. -/
def marginHolds (start end_ : Int) (kept : List (Int × Int)) : Bool :=
  match kept, frameWindow start end_ kept with
  | (p, q) :: t, some w =>
    let srcMin := minOf p (t.map Prod.fst)
    let srcMax := maxOf p (t.map Prod.fst)
    let liftStart := minOf q (t.map Prod.snd)
    let liftEnd := maxOf q (t.map Prod.snd)
    (decide (liftStart - (srcMin - start) < 1) || decide (liftStart - w.1 = srcMin - start)) &&
      (decide (liftEnd + (end_ - srcMax) < w.1 + 1) || decide (w.2 - liftEnd = end_ - srcMax))
  | _, _ => true

def badMargins (n : Nat) : List (List (List (Bool × Int)) × (Int × Int) × List Int) :=
  (tables fullCod n).flatMap fun tbl =>
    (regions n).flatMap fun r =>
      (sublists (positions n)).filterMap fun f =>
        if marginHolds r.1 r.2 (liftoverRegion (tableLifter tbl) (filterByRegion r.1 r.2 f))
        then none else some (tbl, r, f)

#eval (badMargins 4).length   -- 0
#guard (badMargins 4).isEmpty

/-- (d) Shift lifters `x ↦ x + k`, `|k| ≤ n`, on positions `1..2n`: any
1-based position that does not come back as `pos + k`. -/
def badShifts (n : Nat) : List (Int × Int) :=
  (List.range (2 * n + 1)).flatMap fun i =>
    let k : Int := Int.ofNat i - Int.ofNat n
    (positions (2 * n)).filterMap fun pos =>
      if liftOne (fun x => [(true, x + k)]) pos = (.lifted, some (pos + k)) then none
      else some (k, pos)

#eval badShifts 4   -- []
#guard (badShifts 4).isEmpty

/-- (d) A lifter over a non-negative codomain never yields a position `< 1`. -/
def badFloors (n : Nat) : List (List (List (Bool × Int)) × Int) :=
  (tables fullCod n).flatMap fun tbl =>
    (positions n).filterMap fun pos =>
      match (liftOne (tableLifter tbl) pos).2 with
      | some q => if 1 ≤ q then none else some (tbl, pos)
      | none => none

#guard (badFloors 4).isEmpty

/-- (e) Leads inside `[start, end]`, one frame, that `check` rejects given the
window and the lead `lift_window` returns for them. `onlyRows` restricts the
lead to a row of the frame. -/
def badLeads (onlyRows : Bool) (check : Int × Int → Option Int → Option Int → Bool) (n : Nat) :
    List (List (List (Bool × Int)) × (Int × Int) × List Int × Int) :=
  (tables smallCod n).flatMap fun tbl =>
    (regions n).flatMap fun r =>
      (sublists (positions n)).flatMap fun f =>
        (positions n).filterMap fun lead =>
          if decide (r.1 ≤ lead) && decide (lead ≤ r.2) && (!onlyRows || f.contains lead) then
            match liftWindowLeads (tableLifter tbl) r.1 r.2 [f] [some lead] with
            | some (w, [out]) =>
              if check w (liftLead (tableLifter tbl) (some lead)) out then none
              else some (tbl, r, f, lead)
            | some _ => some (tbl, r, f, lead)
            | none => none
          else none

-- No returned lead is outside the window, whether or not it is a row.
#eval (badLeads false (fun w _ out => out.all fun q => decide (Inside w q)) 4).length   -- 0
#guard (badLeads false (fun w _ out => out.all fun q => decide (Inside w q)) 4).isEmpty
-- A lead that is a row of its frame is returned as it lifted: neither check
-- drops it.
#eval (badLeads true (fun _ lifted out => out == lifted) 4).length   -- 0
#guard (badLeads true (fun _ lifted out => out == lifted) 4).isEmpty
-- The checks are not vacuous: some lead that is no row lifts and is dropped.
#guard !(badLeads false (fun _ lifted out => out == lifted) 4).isEmpty

/-! ## Proofs for every size -/

theorem minOf_le_head (a : Int) (l : List Int) : minOf a l ≤ a := by
  induction l generalizing a with
  | nil => exact Int.le_refl a
  | cons b t ih =>
    show minOf (min a b) t ≤ a
    have := ih (min a b)
    omega

theorem minOf_le_mem (a : Int) (l : List Int) (x : Int) (hx : x ∈ l) : minOf a l ≤ x := by
  induction l generalizing a with
  | nil => cases hx
  | cons b t ih =>
    show minOf (min a b) t ≤ x
    rcases List.mem_cons.mp hx with h | h
    · have := minOf_le_head (min a b) t
      omega
    · exact ih (min a b) h

theorem le_minOf (c a : Int) (l : List Int) (ha : c ≤ a) (hl : ∀ x ∈ l, c ≤ x) :
    c ≤ minOf a l := by
  induction l generalizing a with
  | nil => exact ha
  | cons b t ih =>
    show c ≤ minOf (min a b) t
    have hb : c ≤ b := hl b (List.mem_cons.mpr (Or.inl rfl))
    exact ih (min a b) (by omega) fun x hx => hl x (List.mem_cons.mpr (Or.inr hx))

theorem head_le_maxOf (a : Int) (l : List Int) : a ≤ maxOf a l := by
  induction l generalizing a with
  | nil => exact Int.le_refl a
  | cons b t ih =>
    show a ≤ maxOf (max a b) t
    have := ih (max a b)
    omega

theorem mem_le_maxOf (a : Int) (l : List Int) (x : Int) (hx : x ∈ l) : x ≤ maxOf a l := by
  induction l generalizing a with
  | nil => cases hx
  | cons b t ih =>
    show x ≤ maxOf (max a b) t
    rcases List.mem_cons.mp hx with h | h
    · have := head_le_maxOf (max a b) t
      omega
    · exact ih (max a b) h

theorem maxOf_le (c a : Int) (l : List Int) (ha : a ≤ c) (hl : ∀ x ∈ l, x ≤ c) :
    maxOf a l ≤ c := by
  induction l generalizing a with
  | nil => exact ha
  | cons b t ih =>
    show maxOf (max a b) t ≤ c
    have hb : b ≤ c := hl b (List.mem_cons.mpr (Or.inl rfl))
    exact ih (max a b) (by omega) fun x hx => hl x (List.mem_cons.mpr (Or.inr hx))

/-- (a), per frame. No hypothesis: the two clamps alone give it. -/
theorem frame_window_ordered (start end_ : Int) (kept : List (Int × Int)) (w : Int × Int)
    (h : frameWindow start end_ kept = some w) : WindowOrdered w := by
  cases kept with
  | nil => cases h
  | cons hd t =>
    obtain ⟨p, q⟩ := hd
    simp only [frameWindow, Option.some.injEq] at h
    subst h
    simp only [WindowOrdered, windowStart, windowEnd]
    omega

/-- (b), per frame. Hypothesis `hk`: every kept row has its source position in
`[start, end]` (the region filter, `utils.py` l.208) and a lifted position
`≥ 1` (a non-negative lifter plus the `+ 1` of l.181). -/
theorem frame_window_contains (start end_ : Int) (kept : List (Int × Int)) (w : Int × Int)
    (hk : ∀ pq ∈ kept, start ≤ pq.1 ∧ pq.1 ≤ end_ ∧ 1 ≤ pq.2)
    (h : frameWindow start end_ kept = some w) :
    ∀ pq ∈ kept, Inside w pq.2 := by
  cases kept with
  | nil => cases h
  | cons hd t =>
    obtain ⟨p, q⟩ := hd
    simp only [frameWindow, Option.some.injEq] at h
    subst h
    have hhd := hk (p, q) (List.mem_cons.mpr (Or.inl rfl))
    have htl : ∀ y ∈ t, start ≤ y.1 ∧ y.1 ≤ end_ ∧ 1 ≤ y.2 :=
      fun y hy => hk y (List.mem_cons.mpr (Or.inr hy))
    have hsmin : start ≤ minOf p (t.map Prod.fst) :=
      le_minOf start p _ hhd.1 fun x hx => by
        obtain ⟨y, hy, rfl⟩ := List.mem_map.mp hx
        exact (htl y hy).1
    have hsmax : maxOf p (t.map Prod.fst) ≤ end_ :=
      maxOf_le end_ p _ hhd.2.1 fun x hx => by
        obtain ⟨y, hy, rfl⟩ := List.mem_map.mp hx
        exact (htl y hy).2.1
    have hl1 : 1 ≤ minOf q (t.map Prod.snd) :=
      le_minOf 1 q _ hhd.2.2 fun x hx => by
        obtain ⟨y, hy, rfl⟩ := List.mem_map.mp hx
        exact (htl y hy).2.2
    intro pq hpq
    have hlo : minOf q (t.map Prod.snd) ≤ pq.2 := by
      rcases List.mem_cons.mp hpq with h | h
      · subst h; exact minOf_le_head q _
      · exact minOf_le_mem q _ pq.2 (List.mem_map.mpr ⟨pq, h, rfl⟩)
    have hhi : pq.2 ≤ maxOf q (t.map Prod.snd) := by
      rcases List.mem_cons.mp hpq with h | h
      · subst h; exact head_le_maxOf q _
      · exact mem_le_maxOf q _ pq.2 (List.mem_map.mpr ⟨pq, h, rfl⟩)
    simp only [Inside, windowStart, windowEnd]
    omega

theorem frameWindows_of_window (lifter : Lifter) (start end_ : Int) :
    ∀ (fs : List (List Int)) (ws : List (Int × Int)),
      frameWindows lifter start end_ fs = some ws →
      ∀ w ∈ ws, ∃ f ∈ fs, panelWindow lifter start end_ f = some w := by
  intro fs
  induction fs with
  | nil =>
    intro ws h w hw
    simp only [frameWindows, Option.some.injEq] at h
    subst h
    cases hw
  | cons f fs ih =>
    intro ws h w hw
    simp only [frameWindows] at h
    split at h
    · rename_i w0 ws0 h1 h2
      cases h
      rcases List.mem_cons.mp hw with hw' | hw'
      · subst hw'
        exact ⟨f, List.mem_cons.mpr (Or.inl rfl), h1⟩
      · obtain ⟨g, hg, hgw⟩ := ih ws0 h2 w hw'
        exact ⟨g, List.mem_cons.mpr (Or.inr hg), hgw⟩
    · cases h

theorem frameWindows_of_frame (lifter : Lifter) (start end_ : Int) :
    ∀ (fs : List (List Int)) (ws : List (Int × Int)),
      frameWindows lifter start end_ fs = some ws →
      ∀ f ∈ fs, ∃ w ∈ ws, panelWindow lifter start end_ f = some w := by
  intro fs
  induction fs with
  | nil => intro ws _ f hf; cases hf
  | cons f0 fs ih =>
    intro ws h f hf
    simp only [frameWindows] at h
    split at h
    · rename_i w0 ws0 h1 h2
      cases h
      rcases List.mem_cons.mp hf with hf' | hf'
      · subst hf'
        exact ⟨w0, List.mem_cons.mpr (Or.inl rfl), h1⟩
      · obtain ⟨w, hw, hfw⟩ := ih ws0 h2 f hf'
        exact ⟨w, List.mem_cons.mpr (Or.inr hw), hfw⟩
    · cases h

/-- `liftWindow` succeeded: the per-frame windows, and the result as their
min / max. -/
theorem liftWindow_some (lifter : Lifter) (start end_ : Int) (frames : List (List Int))
    (w : Int × Int) (h : liftWindow lifter start end_ frames = some w) :
    ∃ w0 ws, frameWindows lifter start end_ frames = some (w0 :: ws) ∧
      w = (minOf w0.1 (ws.map Prod.fst), maxOf w0.2 (ws.map Prod.snd)) := by
  unfold liftWindow at h
  split at h
  · rename_i w0 ws hfw
    cases h
    exact ⟨w0, ws, hfw, rfl⟩
  · cases h

/-- (a) `lift_window` returns `1 ≤ start < end`, for any lifter, any frames,
any `start` / `end`. No hypothesis beyond "it returned": an empty frame list
or a frame with no lifted row raises instead. -/
theorem lift_window_ordered (lifter : Lifter) (start end_ : Int) (frames : List (List Int))
    (w : Int × Int) (h : liftWindow lifter start end_ frames = some w) : WindowOrdered w := by
  obtain ⟨w0, ws, hfw, rfl⟩ := liftWindow_some lifter start end_ frames w h
  have hall : ∀ v ∈ w0 :: ws, WindowOrdered v := fun v hv => by
    obtain ⟨f, _, hf⟩ := frameWindows_of_window lifter start end_ frames _ hfw v hv
    exact frame_window_ordered start end_ _ v hf
  have h0 := hall w0 (List.mem_cons.mpr (Or.inl rfl))
  have hmin1 : 1 ≤ minOf w0.1 (ws.map Prod.fst) :=
    le_minOf 1 w0.1 _ h0.1 fun x hx => by
      obtain ⟨y, hy, rfl⟩ := List.mem_map.mp hx
      exact (hall y (List.mem_cons.mpr (Or.inr hy))).1
  have hminle := minOf_le_head w0.1 (ws.map Prod.fst)
  have hmaxge := head_le_maxOf w0.2 (ws.map Prod.snd)
  have h02 := h0.2
  simp only [WindowOrdered]
  omega

/-- A lifter with non-negative hits lifts to a position `≥ 1` (d, second half). -/
theorem liftOne_ge_one (lifter : Lifter) (hNonneg : NonnegLifter lifter) (pos q : Int)
    (h : (liftOne lifter pos).2 = some q) : 1 ≤ q := by
  unfold liftOne at h
  split at h
  · cases h
  · rename_i same p heq
    have hp : 0 ≤ p := hNonneg (pos - 1) (same, p) (by rw [heq]; exact List.mem_cons.mpr (Or.inl rfl))
    split at h
    · simp only [Option.some.injEq] at h
      omega
    · cases h
  · cases h

/-- (d) 1-based round trip: under the 0-based shift lifter `x ↦ x + k`,
`_lift_one` returns exactly `pos + k`; the `- 1` of l.174 and the `+ 1` of
l.181 cancel. `k = 0` is the identity lifter. -/
theorem liftOne_shift (lifter : Lifter) (k : Int) (hShift : ∀ x, lifter x = [(true, x + k)])
    (pos : Int) : liftOne lifter pos = (.lifted, some (pos + k)) := by
  unfold liftOne
  rw [hShift]
  simp only [ite_true, Prod.mk.injEq, Option.some.injEq, true_and]
  omega

/-- (b) Containment. Every row of every frame that the region filter selects
and that lifts lies in the returned window.

Hypotheses, all about the real code rather than the model:
* `hNonneg`: the lifter returns 0-based coordinates `≥ 0`, so `1 ≤ lift.start`.
* `hlo`, `hhi`: the row is inside `[start, end]`; `filter_by_region` selects
  exactly those, which gives `start ≤ src_min ≤ src_max ≤ end`.
* `hf`: the frame is one of the frames passed in (so the list is non-empty). -/
theorem lift_window_contains (lifter : Lifter) (hNonneg : NonnegLifter lifter)
    (start end_ : Int) (frames : List (List Int)) (w : Int × Int)
    (h : liftWindow lifter start end_ frames = some w)
    (f : List Int) (hf : f ∈ frames) (p q : Int) (hp : p ∈ f)
    (hlo : start ≤ p) (hhi : p ≤ end_) (hq : (liftOne lifter p).2 = some q) :
    Inside w q := by
  obtain ⟨w0, ws, hfw, rfl⟩ := liftWindow_some lifter start end_ frames w h
  obtain ⟨v, hv, hpanel⟩ := frameWindows_of_frame lifter start end_ frames _ hfw f hf
  unfold panelWindow at hpanel
  have hkept : ∀ pq ∈ liftoverRegion lifter (filterByRegion start end_ f),
      start ≤ pq.1 ∧ pq.1 ≤ end_ ∧ 1 ≤ pq.2 := by
    intro pq hpq
    obtain ⟨a, ha, hmap⟩ := List.mem_filterMap.mp hpq
    obtain ⟨b, hb, hpair⟩ := Option.map_eq_some_iff.mp hmap
    subst hpair
    have hsel := (List.mem_filter.mp ha).2
    simp only [Bool.and_eq_true, decide_eq_true_eq] at hsel
    exact ⟨hsel.1, hsel.2, liftOne_ge_one lifter hNonneg a b hb⟩
  have hmem : (p, q) ∈ liftoverRegion lifter (filterByRegion start end_ f) := by
    refine List.mem_filterMap.mpr ⟨p, List.mem_filter.mpr ⟨hp, ?_⟩, ?_⟩
    · simp only [Bool.and_eq_true, decide_eq_true_eq]
      exact ⟨hlo, hhi⟩
    · rw [hq]; rfl
  have hin := frame_window_contains start end_ _ v hkept hpanel (p, q) hmem
  have hvlo : minOf w0.1 (ws.map Prod.fst) ≤ v.1 := by
    rcases List.mem_cons.mp hv with h | h
    · subst h; exact minOf_le_head _ _
    · exact minOf_le_mem _ _ _ (List.mem_map.mpr ⟨v, h, rfl⟩)
  have hvhi : v.2 ≤ maxOf w0.2 (ws.map Prod.snd) := by
    rcases List.mem_cons.mp hv with h | h
    · subst h; exact head_le_maxOf _ _
    · exact mem_le_maxOf _ _ _ (List.mem_map.mpr ⟨v, h, rfl⟩)
  simp only [Inside] at hin ⊢
  omega

/-- (c) Left margin. `hNoClamp`: the `max(1, …)` of l.430 does not bind. Then
the gap before the first lifted row equals the source gap `src_min - start`. -/
theorem window_start_margin (start srcMin liftStart : Int)
    (hNoClamp : 1 ≤ liftStart - (srcMin - start)) :
    liftStart - windowStart start srcMin liftStart = srcMin - start := by
  simp only [windowStart]
  omega

/-- (c) Right margin. `hNoFloor`: the `window_start + 1` floor of l.432 does
not bind. Then the gap after the last lifted row equals `end - src_max`. -/
theorem window_end_margin (end_ srcMax liftEnd ws : Int)
    (hNoFloor : ws + 1 ≤ liftEnd + (end_ - srcMax)) :
    windowEnd end_ srcMax liftEnd ws - liftEnd = end_ - srcMax := by
  simp only [windowEnd]
  omega

/-- A lead that is a row of its own frame, inside `[start, end]`, and that
lifts, lies in the window, so the check of l.436-443 keeps it. It is (b)
applied to the lead's row. -/
theorem lead_in_frame_inside_window (lifter : Lifter) (hNonneg : NonnegLifter lifter)
    (start end_ : Int) (frames : List (List Int)) (w : Int × Int)
    (h : liftWindow lifter start end_ frames = some w)
    (f : List Int) (hf : f ∈ frames) (lead q : Int) (hRow : lead ∈ f)
    (hlo : start ≤ lead) (hhi : lead ≤ end_) (hq : liftLead lifter (some lead) = some q) :
    Inside w q :=
  lift_window_contains lifter hNonneg start end_ frames w h f hf lead q hRow hlo hhi hq

/-- `keepInside` returns only positions of the window. -/
theorem keepInside_inside (w : Int × Int) (lead : Option Int) (q : Int)
    (h : keepInside w lead = some q) : Inside w q := by
  cases lead with
  | none => cases h
  | some x =>
    simp only [keepInside] at h
    split at h
    · rename_i hin
      cases h
      exact hin
    · cases h

/-- (e) Every lead `lift_window` returns lies in the window it returns. No
hypothesis on the lifter, the frames or the leads: a lead need not be a row of
its frame, nor inside `[start, end]`, and the lifter may return anything. The
check of l.436-443 is what gives it. -/
theorem lift_window_lead_inside (lifter : Lifter) (start end_ : Int)
    (frames : List (List Int)) (leads : List (Option Int))
    (w : Int × Int) (out : List (Option Int))
    (h : liftWindowLeads lifter start end_ frames leads = some (w, out))
    (q : Int) (hq : some q ∈ out) : Inside w q := by
  unfold liftWindowLeads at h
  split at h
  · cases h
    obtain ⟨lead, _, hlead⟩ := List.mem_map.mp hq
    exact keepInside_inside _ lead q hlead
  · cases h

/-- A lead that is no kept row of its frame is never handed on at a kept row's
lifted position (l.416-429), so it cannot make that row's SNP the lead. -/
theorem stray_lead_off_rows (kept : List (Int × Int)) (lifter : Lifter) (p q : Int)
    (hStray : p ∉ kept.map Prod.fst) (h : frameLead kept lifter p = some q) :
    q ∉ kept.map Prod.snd := by
  unfold frameLead at h
  split at h
  · split at h
    · cases h
    · rename_i hc
      cases h
      exact fun hmem => hc ⟨hStray, hmem⟩
  · cases h

/-! ## Concrete examples -/

/-- Rows at 100 and 200 lift to themselves; position 150, which is no row of
the frame, lifts to 5000. -/
def escapeLifter : Lifter := fun x =>
  if x = 99 then [(true, 99)]
  else if x = 149 then [(true, 4999)]
  else if x = 199 then [(true, 199)]
  else []

/-- Why (e) needs the check of l.436-443: `start = 50 ≤ lead = 150 ≤ end = 250`,
the lifter is non-negative, the window is `(50, 250)` and the lead lifts to
`5000`, because it is lifted on its own (`_liftover.py` l.291-293) and never
enters `lift.start` / `lift.end` (l.295-296). `lift_window` returns no lead
for it. -/
theorem lead_not_a_row_is_dropped :
    liftLead escapeLifter (some 150) = some 5000 ∧
      ¬ Inside (50, 250) 5000 ∧
      liftWindowLeads escapeLifter 50 250 [[100, 200]] [some 150] =
        some ((50, 250), [none]) := by
  decide

/-- Row 100 lifts to itself; row 200 and position 150, which is no row of the
frame, both lift to 300. -/
def collideLifter : Lifter := fun x =>
  if x = 99 then [(true, 99)]
  else if x = 149 then [(true, 299)]
  else if x = 199 then [(true, 299)]
  else []

/-- Why the check of l.416-429 exists: lead 150 lifts to `300`, inside the
window `(50, 350)` and on the lifted position of row 200, which a lookup by
position would take for the lead. `lift_window` returns no lead for it, and
still returns `300` when the lead is row 200 itself. -/
theorem stray_lead_on_a_row_is_dropped :
    liftLead collideLifter (some 150) = some 300 ∧
      Inside (50, 350) 300 ∧
      liftWindowLeads collideLifter 50 250 [[100, 200]] [some 150] =
        some ((50, 350), [none]) ∧
      liftWindowLeads collideLifter 50 250 [[100, 200]] [some 200] =
        some ((50, 350), [some 300]) := by
  decide

/-- (b) needs `NonnegLifter`: a hit at 0-based `-5` lifts row 3 to `-4`, and
the `max(1, …)` clamp of l.430 puts the window start at 1, right of it. -/
theorem negative_lifter_escapes :
    liftWindow (fun _ => [(true, -5)]) 1 5 [[3]] = some (1, 2) ∧
      (liftOne (fun _ => [(true, -5)]) 3).2 = some (-4) ∧
      ¬ Inside (1, 2) (-4) := by
  decide
