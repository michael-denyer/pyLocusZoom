/-!
# Greedy gene-row assignment

Model of `src/pylocuszoom/gene_track.py::assign_gene_positions`, as called by
`src/pylocuszoom/panels/genes.py::GenePanel.from_genes`.

Modelling decisions:

* Coordinates are Python unbounded ints, so `Int`. Rows are `Nat`.
* `label_buffer = region_width * 0.08` is a float. The
  model uses the exact rational 2/25 scaled by 25: the loop test
  `row_ends[row] > gene_start - label_buffer` becomes
  `25 * rowEnd > 25 * geneStart - 2 * regionWidth`. The float product can
  differ from the exact value by a rounding error; that gap is measured from
  Python and reported separately, it is not part of this model.
* `row_ends: dict[int, int]` is a `List Int` indexed by
  row. The loop only ever inserts at the first row that is absent, so the key
  set is always `{0, .., k-1}`; `row in row_ends` is `row < rowEnds.length`.
* The frame is a `List Gene` in iteration order. Nothing is assumed about
  order, about `start ≤ end`, or about the region; each theorem names the
  hypotheses it needs. The caller does enforce `start ≤ end`: see `WellFormed`.
-/

namespace GeneRows

/-- One row of `genes_df`: the `start` and `end` columns (`end` is a Lean
keyword, hence `stop`). -/
structure Gene where
  start : Int
  stop : Int
deriving Repr, DecidableEq

/-- `gene_start = max(gene["start"], start)` in `assign_gene_positions`. -/
def clipStart (S : Int) (g : Gene) : Int := max g.start S

/-- `gene_end = min(gene["end"], end)` in `assign_gene_positions`. -/
def clipEnd (E : Int) (g : Gene) : Int := min g.stop E

/-- `row_ends[row] > gene_start - label_buffer`, the `while` test of
`assign_gene_positions`, with `label_buffer = W * 0.08` scaled by 25. -/
def Collides (W rowEnd geneStart : Int) : Prop := 25 * rowEnd > 25 * geneStart - 2 * W

instance (W rowEnd geneStart : Int) : Decidable (Collides W rowEnd geneStart) := by
  unfold Collides; infer_instance

/-- The `while` loop of `assign_gene_positions`: first row that is absent or whose
recorded end does not collide. -/
def findRow (W geneStart : Int) : List Int → Nat
  | [] => 0
  | e :: es => if Collides W e geneStart then findRow W geneStart es + 1 else 0

/-- `row_ends[row] = gene_end` in `assign_gene_positions`: overwrite, or append
when the row is new. -/
def setRow : List Int → Nat → Int → List Int
  | [], _, v => [v]
  | _ :: es, 0, v => v :: es
  | e :: es, r + 1, v => e :: setRow es r v

/-- The `for` loop body of `assign_gene_positions`, threaded over `row_ends`. -/
def go (S E : Int) : List Int → List Gene → List Nat
  | _, [] => []
  | rowEnds, g :: rest =>
    findRow (E - S) (clipStart S g) rowEnds ::
      go S E (setRow rowEnds (findRow (E - S) (clipStart S g) rowEnds) (clipEnd E g)) rest

/-- `assign_gene_positions(genes_df, start, end)` in `gene_track.py`. -/
def assign_gene_positions (genes : List Gene) (S E : Int) : List Nat := go S E [] genes

/-! ## Properties -/

/-- Two placed genes do not collide when they share a row: the later one's
clipped start minus the label buffer is at or past the earlier one's clipped
end. -/
def Clear (S E : Int) (p q : Nat × Gene) : Prop :=
  p.1 = q.1 → 25 * clipEnd E p.2 ≤ 25 * clipStart S q.2 - 2 * (E - S)

instance (S E : Int) (p q : Nat × Gene) : Decidable (Clear S E p q) := by
  unfold Clear; infer_instance

/-- (a) Non-overlap: every earlier gene in a row is clear of every later gene
in that row, not only the most recent one. -/
def NonOverlap (S E : Int) (genes : List Gene) (rows : List Nat) : Prop :=
  (rows.zip genes).Pairwise (Clear S E)

instance (S E : Int) (genes : List Gene) (rows : List Nat) :
    Decidable (NonOverlap S E genes rows) := by
  unfold NonOverlap; infer_instance

/-- (b1) Gene `i` gets a row `≤ i`. -/
def RowLeIndex (rows : List Nat) : Prop := ∀ x ∈ rows.zipIdx, x.1 ≤ x.2

instance (rows : List Nat) : Decidable (RowLeIndex rows) := by
  unfold RowLeIndex; infer_instance

/-- (b2) and (c): for each gene, every lower row is in use by an earlier gene
that collides with it, so the gene's row is the lowest it could take.
`earlier` holds the genes already placed. The collision arithmetic is written
out here, independent of `Collides`, so a wrong comparison in the model's
loop test cannot satisfy the property by changing both sides at once. -/
def Greedy (S E : Int) : List (Nat × Gene) → List (Nat × Gene) → Prop
  | _, [] => True
  | earlier, p :: rest =>
    (∀ r, r < p.1 →
      ∃ q ∈ earlier, q.1 = r ∧ 25 * clipEnd E q.2 > 25 * clipStart S p.2 - 2 * (E - S)) ∧
    Greedy S E (earlier ++ [p]) rest

instance decGreedy (S E : Int) :
    ∀ earlier rest, Decidable (Greedy S E earlier rest)
  | _, [] => isTrue trivial
  | earlier, p :: rest => by
    unfold Greedy
    have := decGreedy S E (earlier ++ [p]) rest
    infer_instance

/-- Hypothesis of (a): the gene's clipped end is not left of its clipped start
by more than the label buffer, so writing it to `row_ends` cannot lower the
row's recorded end. -/
def NoShrink (S E : Int) (g : Gene) : Prop :=
  25 * clipStart S g - 2 * (E - S) ≤ 25 * clipEnd E g

instance (S E : Int) (g : Gene) : Decidable (NoShrink S E g) := by
  unfold NoShrink; infer_instance

/-- What the real caller supplies: `start ≤ end` per gene, enforced by the
`ordering` rule of `schemas.py::GENES_PLOT` that `GenePanel.from_genes` checks
before the row assignment (`check(genes_df, GENES_PLOT)`, rule run by
`validation.py::check`), and the gene intersects the region
(`gene_track.py::filter_genes_by_region`). -/
def WellFormed (S E : Int) (g : Gene) : Prop :=
  g.start ≤ g.stop ∧ S ≤ g.stop ∧ g.start ≤ E

instance (S E : Int) (g : Gene) : Decidable (WellFormed S E g) := by
  unfold WellFormed; infer_instance

/-- Alternative hypothesis of (a): the frame is sorted by `start`, as
`GenePanel.from_genes` does with `.sort_values("start")`. -/
def SortedByStart (genes : List Gene) : Prop := genes.Pairwise fun a b => a.start ≤ b.start

instance (genes : List Gene) : Decidable (SortedByStart genes) := by
  unfold SortedByStart; infer_instance

/-- The x-extent `GenePanel.draw` paints for a gene: a rectangle anchored at
`gene_start` with width `gene_end - gene_start` (`panels/genes.py::_gene_band`),
so it covers `[min, max]` of the two clipped ends. Two glyphs in
one row are disjoint when the earlier one's right edge is at or left of the
later one's left edge. -/
def GlyphClear (S E : Int) (p q : Nat × Gene) : Prop :=
  p.1 = q.1 →
    max (clipStart S p.2) (clipEnd E p.2) ≤ min (clipStart S q.2) (clipEnd E q.2)

instance (S E : Int) (p q : Nat × Gene) : Decidable (GlyphClear S E p q) := by
  unfold GlyphClear; infer_instance

def GlyphsDisjoint (S E : Int) (genes : List Gene) (rows : List Nat) : Prop :=
  (rows.zip genes).Pairwise (GlyphClear S E)

instance (S E : Int) (genes : List Gene) (rows : List Nat) :
    Decidable (GlyphsDisjoint S E genes rows) := by
  unfold GlyphsDisjoint; infer_instance

/-! ## Bounded exhaustive check -/

/-- Coordinates around both region edges and on the buffer boundary
(buffer is 2 for width 25 and 4 for width 50). -/
def coords (S E : Int) : List Int := [S - 1, S, S + 2, S + 4, E - 2, E, E + 1]

/-- Every gene over `coords`, including `start > end`. -/
def allGenes (cs : List Int) : List Gene :=
  cs.flatMap fun a => cs.map fun b => ⟨a, b⟩

/-- Every list over `xs` of length at most `n`. -/
def listsUpTo (xs : List Gene) : Nat → List (List Gene)
  | 0 => [[]]
  | n + 1 => [] :: xs.flatMap fun x => (listsUpTo xs n).map (x :: ·)

/-- Regions: widths 25 and 50 (integer buffer), 3 (fractional buffer), 0, and
a reversed region. `RegionConfig` only admits `1 ≤ start < end`. -/
def regions : List (Int × Int) := [(0, 25), (0, 50), (7, 10), (3, 3), (30, 5)]

/-- Inputs (region, genes, rows) that satisfy `pre` and break `post`. -/
def bad (rs : List (Int × Int)) (cs : Int → Int → List Int) (n : Nat)
    (pre : Int → Int → List Gene → Bool) (post : Int → Int → List Gene → List Nat → Bool) :
    List ((Int × Int) × List Gene × List Nat) :=
  rs.flatMap fun (S, E) =>
    (listsUpTo (allGenes (cs S E)) n).filterMap fun genes =>
      let rows := assign_gene_positions genes S E
      if pre S E genes && !post S E genes rows then some ((S, E), genes, rows) else none

def preNoShrink (S E : Int) (genes : List Gene) : Bool := genes.all (NoShrink S E ·)
def preWellFormed (S E : Int) (genes : List Gene) : Bool :=
  decide (S ≤ E) && genes.all (WellFormed S E ·)
def preSorted (_ _ : Int) (genes : List Gene) : Bool := decide (SortedByStart genes)
def postGlyphs (S E : Int) (genes : List Gene) (rows : List Nat) : Bool :=
  decide (GlyphsDisjoint S E genes rows)
def preReversedRegion (S E : Int) (_ : List Gene) : Bool := decide (E ≤ S)
def postNonOverlap (S E : Int) (genes : List Gene) (rows : List Nat) : Bool :=
  decide (NonOverlap S E genes rows)
def postTight (S E : Int) (genes : List Gene) (rows : List Nat) : Bool :=
  decide (RowLeIndex rows) && decide (Greedy S E [] (rows.zip genes)) &&
    rows.length == genes.length

-- (a) holds under `NoShrink`, under `WellFormed`, and on a reversed or empty region.
#guard (bad regions coords 3 preNoShrink postNonOverlap).isEmpty
#guard (bad regions coords 3 preWellFormed postNonOverlap).isEmpty
#guard (bad regions coords 3 preReversedRegion postNonOverlap).isEmpty
-- (b) and (c) hold with no hypothesis at all.
#guard (bad regions coords 3 (fun _ _ _ => true) postTight).isEmpty
-- Length 4 on a smaller coordinate set.
#guard (bad regions (fun S E => [S, S + 2, S + 4, E]) 4 preNoShrink postNonOverlap).isEmpty
#guard (bad regions (fun S E => [S, S + 2, S + 4, E]) 4 (fun _ _ _ => true) postTight).isEmpty

-- (a) FAILS with no hypothesis: a gene with `end < start` lowers the row's
-- recorded end and hides the earlier, longer gene. First counter-example:
#eval (bad [(0, 25)] coords 3 (fun _ _ _ => true) postNonOverlap).head?
#guard !(bad [(0, 25)] coords 3 (fun _ _ _ => true) postNonOverlap).isEmpty
-- The same shape at the scale of the Python reproduction: region 1..1001
-- (buffer 80), genes (1,500), (600,100), (200,300) all land in row 0.
#guard assign_gene_positions [⟨1, 500⟩, ⟨600, 100⟩, ⟨200, 300⟩] 1 1001 = [0, 0, 0]
#guard !decide (NonOverlap 1 1001 [⟨1, 500⟩, ⟨600, 100⟩, ⟨200, 300⟩] [0, 0, 0])
-- No counter-example has fewer than three genes.
#guard (bad regions coords 2 (fun _ _ _ => true) postNonOverlap).isEmpty
-- (a) also holds for any frame sorted by start, reversed genes included.
#guard (bad regions coords 3 preSorted postNonOverlap).isEmpty

-- (a) is a statement about (clipStart, clipEnd). `GenePanel.draw` paints the
-- rectangle from `gene_start` with width `gene_end - gene_start`
-- (`panels/genes.py::_gene_band`), which for a gene with `end < start` covers
-- `[gene_end, gene_start]`. Sorted input: (1,500) and (600,100) share row 0,
-- (a) holds, and the glyphs 1..500 and 100..600 overlap. `GenePanel.from_genes`
-- now rejects this frame (the `ordering` rule of `GENES_PLOT`).
#guard assign_gene_positions [⟨1, 500⟩, ⟨600, 100⟩] 1 1001 = [0, 0]
#guard decide (NonOverlap 1 1001 [⟨1, 500⟩, ⟨600, 100⟩] [0, 0])
#guard !decide (GlyphsDisjoint 1 1001 [⟨1, 500⟩, ⟨600, 100⟩] [0, 0])
-- With `WellFormed` genes the glyphs in a row are disjoint.
#guard (bad regions coords 3 preWellFormed postGlyphs).isEmpty

/-! ## Proofs for every size -/

theorem findRow_le (W gs : Int) : ∀ re, findRow W gs re ≤ re.length
  | [] => by simp [findRow]
  | e :: es => by
    have := findRow_le W gs es
    simp only [findRow, List.length_cons]
    split <;> omega

/-- The row found is absent or does not collide. -/
theorem findRow_free (W gs : Int) :
    ∀ re e, re[findRow W gs re]? = some e → ¬ Collides W e gs := by
  intro re
  induction re with
  | nil => intro e h; simp [findRow] at h
  | cons a es ih =>
    intro e h
    simp only [findRow] at h
    split at h
    · simp only [List.getElem?_cons_succ] at h
      exact ih e h
    · simp only [List.getElem?_cons_zero, Option.some.injEq] at h
      subst h; assumption

/-- Every row below the one found is present and collides. -/
theorem findRow_lower (W gs : Int) :
    ∀ re r, r < findRow W gs re → ∃ e, re[r]? = some e ∧ Collides W e gs := by
  intro re
  induction re with
  | nil => intro r h; simp [findRow] at h
  | cons a es ih =>
    intro r h
    simp only [findRow] at h
    split at h
    · cases r with
      | zero => exact ⟨a, by simp, by assumption⟩
      | succ r =>
        have ⟨e, he, hc⟩ := ih r (by omega)
        exact ⟨e, by simpa using he, hc⟩
    · omega

theorem setRow_self : ∀ (re : List Int) (r : Nat) (v : Int),
    r ≤ re.length → (setRow re r v)[r]? = some v := by
  intro re
  induction re with
  | nil => intro r v h; simp at h; subst h; simp [setRow]
  | cons a es ih =>
    intro r v h
    cases r with
    | zero => simp [setRow]
    | succ r => simp only [setRow, List.getElem?_cons_succ]; exact ih r v (by simpa using h)

theorem setRow_other : ∀ (re : List Int) (r r' : Nat) (v : Int),
    r ≤ re.length → r' ≠ r → (setRow re r v)[r']? = re[r']? := by
  intro re
  induction re with
  | nil =>
    intro r r' v h hne
    simp at h; subst h
    cases r' with
    | zero => exact absurd rfl hne
    | succ r' => simp [setRow]
  | cons a es ih =>
    intro r r' v h hne
    cases r with
    | zero =>
      cases r' with
      | zero => exact absurd rfl hne
      | succ r' => simp [setRow]
    | succ r =>
      cases r' with
      | zero => simp [setRow]
      | succ r' =>
        simp only [setRow, List.getElem?_cons_succ]
        exact ih r r' v (by simpa using h) (by omega)

theorem setRow_length : ∀ (re : List Int) (r : Nat) (v : Int),
    (setRow re r v).length ≤ re.length + 1 := by
  intro re
  induction re with
  | nil => intro r v; simp [setRow]
  | cons a es ih =>
    intro r v
    cases r with
    | zero => simp [setRow]
    | succ r => have := ih r v; simp only [setRow, List.length_cons]; omega

/-- Under `NoShrink`, every gene placed in a row is clear of the end that row
recorded before the run: the recorded end never decreases. -/
theorem go_clears_state (S E : Int) : ∀ (genes : List Gene) (re : List Int),
    (∀ g ∈ genes, NoShrink S E g) →
    ∀ p ∈ (go S E re genes).zip genes, ∀ e, re[p.1]? = some e →
      25 * e ≤ 25 * clipStart S p.2 - 2 * (E - S) := by
  intro genes
  induction genes with
  | nil => intro re _ p hp; simp [go] at hp
  | cons g rest ih =>
    intro re hv p hp e he
    simp only [go, List.zip_cons_cons, List.mem_cons] at hp
    have hfree := findRow_free (E - S) (clipStart S g) re
    have hle := findRow_le (E - S) (clipStart S g) re
    have hrest : ∀ g' ∈ rest, NoShrink S E g' := fun g' h => hv g' (List.mem_cons_of_mem _ h)
    rcases hp with rfl | hp
    · have := hfree e he
      simp only [Collides] at this ⊢
      omega
    · by_cases hr : p.1 = findRow (E - S) (clipStart S g) re
      · have h1 := ih _ hrest p hp (clipEnd E g) (by rw [hr]; exact setRow_self re _ _ hle)
        have h2 := hfree e (by rw [← hr]; exact he)
        have h3 : NoShrink S E g := hv g (List.mem_cons_self ..)
        simp only [Collides, NoShrink] at h2 h3
        omega
      · exact ih _ hrest p hp e (by rw [setRow_other re _ _ _ hle hr]; exact he)

theorem go_nonOverlap (S E : Int) : ∀ (genes : List Gene) (re : List Int),
    (∀ g ∈ genes, NoShrink S E g) → NonOverlap S E genes (go S E re genes) := by
  intro genes
  induction genes with
  | nil => intro re _; simp [NonOverlap, go]
  | cons g rest ih =>
    intro re hv
    have hle := findRow_le (E - S) (clipStart S g) re
    have hrest : ∀ g' ∈ rest, NoShrink S E g' := fun g' h => hv g' (List.mem_cons_of_mem _ h)
    simp only [NonOverlap, go, List.zip_cons_cons, List.pairwise_cons]
    refine ⟨?_, ih _ hrest⟩
    intro q hq heq
    exact go_clears_state S E rest _ hrest q hq (clipEnd E g)
      (by rw [← heq]; exact setRow_self re _ _ hle)

/-- (a), main form. Hypothesis `hv` (`NoShrink`) is a fact about the input
frame, not about the code. No sortedness and no region hypothesis is needed. -/
theorem assign_nonOverlap (S E : Int) (genes : List Gene)
    (hv : ∀ g ∈ genes, NoShrink S E g) :
    NonOverlap S E genes (assign_gene_positions genes S E) :=
  go_nonOverlap S E genes [] hv

/-- With a frame sorted by start, every gene placed in a row is clear of the
end that row recorded before the run. No `NoShrink` needed: the next occupant
of the row cleared that end, and later starts are no smaller. -/
theorem go_clears_state_sorted (S E : Int) : ∀ (genes : List Gene) (re : List Int),
    SortedByStart genes →
    ∀ p ∈ (go S E re genes).zip genes, ∀ e, re[p.1]? = some e →
      25 * e ≤ 25 * clipStart S p.2 - 2 * (E - S) := by
  intro genes
  induction genes with
  | nil => intro re _ p hp; simp [go] at hp
  | cons g rest ih =>
    intro re hs p hp e he
    simp only [go, List.zip_cons_cons, List.mem_cons] at hp
    have hfree := findRow_free (E - S) (clipStart S g) re
    have hle := findRow_le (E - S) (clipStart S g) re
    simp only [SortedByStart, List.pairwise_cons] at hs
    rcases hp with rfl | hp
    · have := hfree e he
      simp only [Collides] at this ⊢
      omega
    · by_cases hr : p.1 = findRow (E - S) (clipStart S g) re
      · have h2 := hfree e (by rw [← hr]; exact he)
        have h3 := hs.1 p.2 (List.of_mem_zip hp).2
        simp only [Collides, clipStart] at h2 ⊢
        omega
      · exact ih _ hs.2 p hp e (by rw [setRow_other re _ _ _ hle hr]; exact he)

/-- (a) for a frame sorted by start. No hypothesis on `start ≤ end` per gene
or on the region. -/
theorem assign_nonOverlap_sorted (S E : Int) (genes : List Gene) (hs : SortedByStart genes) :
    NonOverlap S E genes (assign_gene_positions genes S E) := by
  unfold assign_gene_positions
  generalize ([] : List Int) = re
  induction genes generalizing re with
  | nil => simp [NonOverlap, go]
  | cons g rest ih =>
    have hle := findRow_le (E - S) (clipStart S g) re
    simp only [SortedByStart, List.pairwise_cons] at hs
    simp only [NonOverlap, go, List.zip_cons_cons, List.pairwise_cons]
    refine ⟨?_, ih hs.2 _⟩
    intro q hq heq
    exact go_clears_state_sorted S E rest _ hs.2 q hq (clipEnd E g)
      (by rw [← heq]; exact setRow_self re _ _ hle)

/-- (a) implies disjoint glyphs when the region has `start ≤ end` and every
gene is `WellFormed`. Without `WellFormed` it does not: see the `#guard`s. -/
theorem glyphsDisjoint_of_nonOverlap (S E : Int) (genes : List Gene) (rows : List Nat)
    (hSE : S ≤ E) (hv : ∀ g ∈ genes, WellFormed S E g) (h : NonOverlap S E genes rows) :
    GlyphsDisjoint S E genes rows := by
  unfold GlyphsDisjoint
  unfold NonOverlap at h
  refine List.Pairwise.imp_of_mem ?_ h
  intro p q hp hq hc heq
  have hp' := hv p.2 (List.of_mem_zip hp).2
  have hq' := hv q.2 (List.of_mem_zip hq).2
  have := hc heq
  simp only [WellFormed, clipStart, clipEnd] at *
  omega

theorem wellFormed_noShrink (S E : Int) (g : Gene) (hSE : S ≤ E) (h : WellFormed S E g) :
    NoShrink S E g := by
  simp only [WellFormed, NoShrink, clipStart, clipEnd] at *
  omega

/-- (a) for the real caller's contract: region `start ≤ end`, every gene has
`start ≤ end` and intersects the region. -/
theorem assign_nonOverlap_wellFormed (S E : Int) (genes : List Gene) (hSE : S ≤ E)
    (hv : ∀ g ∈ genes, WellFormed S E g) :
    NonOverlap S E genes (assign_gene_positions genes S E) :=
  assign_nonOverlap S E genes fun g hg => wellFormed_noShrink S E g hSE (hv g hg)

/-- (a) on an empty or reversed region (`end ≤ start`): clipping alone forces
it, for any rows whatsoever and any genes. -/
theorem nonOverlap_of_reversed_region (S E : Int) (genes : List Gene) (rows : List Nat)
    (hSE : E ≤ S) : NonOverlap S E genes rows := by
  apply List.Pairwise.imp (R := fun _ _ => True)
  · intro p q _ _
    simp only [clipStart, clipEnd]
    omega
  · exact List.pairwise_of_forall (fun _ _ => trivial)

theorem go_length (S E : Int) : ∀ (genes : List Gene) (re : List Int),
    (go S E re genes).length = genes.length := by
  intro genes
  induction genes with
  | nil => intro re; simp [go]
  | cons g rest ih => intro re; simp [go, ih]

/-- One row per gene. -/
theorem assign_length (S E : Int) (genes : List Gene) :
    (assign_gene_positions genes S E).length = genes.length := go_length S E genes []

theorem go_rowLeIndex (S E : Int) : ∀ (genes : List Gene) (re : List Int) (k : Nat),
    ∀ x ∈ (go S E re genes).zipIdx k, x.1 + k ≤ re.length + x.2 := by
  intro genes
  induction genes with
  | nil => intro re k x hx; simp [go] at hx
  | cons g rest ih =>
    intro re k x hx
    simp only [go, List.zipIdx_cons, List.mem_cons] at hx
    rcases hx with rfl | hx
    · have := findRow_le (E - S) (clipStart S g) re
      simp only
      omega
    · have h1 := ih _ (k + 1) x hx
      have h2 := setRow_length re (findRow (E - S) (clipStart S g) re) (clipEnd E g)
      omega

/-- (b1) Gene `i` gets row `≤ i`. No hypothesis. -/
theorem assign_rowLeIndex (S E : Int) (genes : List Gene) :
    RowLeIndex (assign_gene_positions genes S E) := by
  intro x hx
  have := go_rowLeIndex S E genes [] 0 x hx
  simpa using this

theorem go_greedy (S E : Int) : ∀ (genes : List Gene) (re : List Int)
    (earlier : List (Nat × Gene)),
    (∀ r e, re[r]? = some e → ∃ q ∈ earlier, q.1 = r ∧ clipEnd E q.2 = e) →
    Greedy S E earlier ((go S E re genes).zip genes) := by
  intro genes
  induction genes with
  | nil => intro re earlier _; simp [go, Greedy]
  | cons g rest ih =>
    intro re earlier H
    have hle := findRow_le (E - S) (clipStart S g) re
    simp only [go, List.zip_cons_cons, Greedy]
    refine ⟨?_, ih _ _ ?_⟩
    · intro r hr
      have ⟨e, he, hc⟩ := findRow_lower (E - S) (clipStart S g) re r hr
      have ⟨q, hq, hq1, hq2⟩ := H r e he
      exact ⟨q, hq, hq1, by rw [hq2]; exact hc⟩
    · intro r e he
      by_cases hr : r = findRow (E - S) (clipStart S g) re
      · subst hr
        rw [setRow_self re _ _ hle] at he
        exact ⟨_, List.mem_append_right _ (List.mem_singleton.mpr rfl), rfl,
          Option.some.inj he⟩
      · rw [setRow_other re _ _ _ hle hr] at he
        have ⟨q, hq, hq1, hq2⟩ := H r e he
        exact ⟨q, List.mem_append_left _ hq, hq1, hq2⟩

/-- (b2) and (c): every row below a gene's row holds an earlier gene that
collides with it, so rows are gap-free and each gene takes the lowest row whose
recorded end it clears. No hypothesis. -/
theorem assign_greedy (S E : Int) (genes : List Gene) :
    Greedy S E [] ((assign_gene_positions genes S E).zip genes) :=
  go_greedy S E genes [] [] (by intro r e h; simp at h)

end GeneRows
