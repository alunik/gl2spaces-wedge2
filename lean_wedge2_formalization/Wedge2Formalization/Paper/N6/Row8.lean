import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N4
import Wedge2Formalization.Paper.N5
import Wedge2Formalization.N4PaperSummary
import Wedge2Formalization.N6Summary
import Wedge2Formalization.N6WeightedTwoPoint
import Wedge2Formalization.N6DoublePureSingular
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.LinearAlgebra.Matrix.Rank

namespace Wedge2Formalization
namespace Paper
namespace N6

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 6`, row 8.
Representative `⟨e₂∧e₄ + e₅∧e₆, e₃∧e₄⟩`.
Divisor `[a]`.
Claimed stabilizer:
`K_L = U_9 \rtimes (G_m(k) \times G_m(k) \times SL_2(k))`, exact quotient
family `Q_L = B`.
-/
namespace Row8

def rep₁ :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  Wedge2Formalization.N6SimplePoint.radSimpleRep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  Wedge2Formalization.N6SimplePoint.radSimpleRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 8. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 8. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  rep₂ (k := k)

/-- The row-8 paper representative already agrees with the internal working pair. -/
def paperChange :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6SimplePoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6SimplePoint.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6SimplePoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6SimplePoint.ActBivector]

/-- Exact pointwise kernel family for the radical simple-point row. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6SimplePoint.V
        Wedge2Formalization.N6SimplePoint.V k) :=
  { g |
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N6SimplePoint.W
            Wedge2Formalization.N6SimplePoint.I k,
          ∃ D : Matrix Wedge2Formalization.N6SimplePoint.W
              Wedge2Formalization.N6SimplePoint.W k,
            D ∈ Wedge2Formalization.Paper.N5.Row3.K (k := k) ∧
              g =
                Matrix.fromBlocks
                  (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) u)
                  0
                  C
                  D }

/-- The embedded simple-point core, i.e. the `U_4 \rtimes (G_m(k) \times SL_2(k))`
factor from Appendix A, `n = 5`, row 3. -/
def SimplePointCore :
    Set
      (Matrix Wedge2Formalization.N6SimplePoint.W
        Wedge2Formalization.N6SimplePoint.W k) :=
  Wedge2Formalization.Paper.N5.Row3.K (k := k)

/-- The displayed radical unipotent family. -/
def U
    (C : Matrix Wedge2Formalization.N6SimplePoint.W
      Wedge2Formalization.N6SimplePoint.I k)
    (x y p q : k) :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) (1 : k))
    0
    C
    (Wedge2Formalization.Paper.N5.Row3.U (k := k) x y p q)

/-- The displayed Levi/core family. -/
def Levi
    (u a : k)
    (E : Matrix Wedge2Formalization.N5SimplePoint.I
      Wedge2Formalization.N5SimplePoint.I k) :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) u)
    0
    0
    (Wedge2Formalization.Paper.N5.Row3.Levi (k := k) a E)

/-- Exact coefficient-side projective Borel family. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

/-- Chosen Borel lift. -/
def lift
    (a b : k) :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k :=
  let G : Matrix Wedge2Formalization.N5SimplePoint.W
      Wedge2Formalization.N5SimplePoint.W k :=
    Wedge2Formalization.N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
  let E : Matrix Wedge2Formalization.N5SimplePoint.I
      Wedge2Formalization.N5SimplePoint.I k := !![a, 0; 0, 1]
  let h : Matrix Wedge2Formalization.N6SimplePoint.W
      Wedge2Formalization.N6SimplePoint.W k := Matrix.fromBlocks G 0 0 E
  Matrix.fromBlocks
    (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) (1 : k))
    0
    0
    h

private theorem scalarBlock_eq
    (A : Matrix Wedge2Formalization.N6SimplePoint.I
      Wedge2Formalization.N6SimplePoint.I k) :
    A = Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) (A 0 0) := by
  ext i j
  fin_cases i <;> fin_cases j
  simp [Wedge2Formalization.N6SimplePoint.scalarBlock]

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6SimplePoint.FixesRadSimplePairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    have hDfix :
        Wedge2Formalization.N5SimplePoint.FixesPairBivector (k := k) g.toBlocks₂₂ :=
      Wedge2Formalization.N6SimplePoint.fixesRadSimplePair_lowerRight
        (k := k)
        (A := g.toBlocks₁₁)
        (B := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)
        (by simpa [Matrix.fromBlocks_toBlocks] using hg)
    have hDmem :
        g.toBlocks₂₂ ∈ Wedge2Formalization.Paper.N5.Row3.K (k := k) :=
      (Wedge2Formalization.Paper.N5.Row3.pointwise_stabilizer
        (k := k)
        g.toBlocks₂₂).1 hDfix
    rcases
      (Wedge2Formalization.Paper.N5.Row3.mem_K_iff
        (k := k)
        g.toBlocks₂₂).1 hDmem with
      ⟨a, x, y, p, q, E, ha, hE, hshape⟩
    let g' :
        Matrix Wedge2Formalization.N6SimplePoint.V
          Wedge2Formalization.N6SimplePoint.V k :=
      Matrix.fromBlocks
        (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) (g.toBlocks₁₁ 0 0))
        g.toBlocks₁₂
        g.toBlocks₂₁
        (Matrix.fromBlocks
          (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
          (!![(0 : k), 0; 0, 0; p, q])
          (!![a * (p * E 0 1 - q * E 0 0), 0, 0;
             a * (p * E 1 1 - q * E 1 0), 0, 0])
          E)
    have hgEq : g' = g := by
      calc
        g' =
            Matrix.fromBlocks
              (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) (g.toBlocks₁₁ 0 0))
              g.toBlocks₁₂
              g.toBlocks₂₁
              g.toBlocks₂₂ := by
                simp [g', hshape]
        _ =
            Matrix.fromBlocks
              g.toBlocks₁₁
              g.toBlocks₁₂
              g.toBlocks₂₁
              g.toBlocks₂₂ := by
                ext i j
                rcases i with i | i <;> rcases j with j | j
                ·
                  fin_cases i
                  fin_cases j
                  simp [Wedge2Formalization.N6SimplePoint.scalarBlock]
                · simp [Matrix.fromBlocks]
                · simp [Matrix.fromBlocks]
                · simp [Matrix.fromBlocks]
        _ = g := Matrix.fromBlocks_toBlocks g
    have hg' :
        Wedge2Formalization.N6SimplePoint.FixesRadSimplePairBivector
          g' := by
      rw [hgEq]
      exact hg
    have hshapeN6 :=
      (Wedge2Formalization.N6Summary.simplePoint_nested_iff_shape
        (k := k)
        (u := g.toBlocks₁₁ 0 0)
        (B0 := g.toBlocks₁₂)
        (C0 := g.toBlocks₂₁)
        (a := a)
        (x := x)
        (y := y)
        (B := !![(0 : k), 0; 0, 0; p, q])
        (C := !![a * (p * E 0 1 - q * E 0 0), 0, 0;
                 a * (p * E 1 1 - q * E 1 0), 0, 0])
        (E := E)
        ha).1 hg'
    rcases hshapeN6 with ⟨hB0, hdetE, _⟩
    refine ⟨g.toBlocks₁₁ 0 0, g.toBlocks₂₁, g.toBlocks₂₂, hDmem, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ =
          Matrix.fromBlocks
            g.toBlocks₁₁
            0
            g.toBlocks₂₁
            g.toBlocks₂₂ := by simp [hB0]
      _ =
          Matrix.fromBlocks
            (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) (g.toBlocks₁₁ 0 0))
            0
            g.toBlocks₂₁
            g.toBlocks₂₂ := by
              ext i j
              rcases i with i | i <;> rcases j with j | j
              ·
                fin_cases i
                fin_cases j
                simp [Wedge2Formalization.N6SimplePoint.scalarBlock, Matrix.fromBlocks]
              · simp [Matrix.fromBlocks]
              · simp [Matrix.fromBlocks]
              · simp [Matrix.fromBlocks]
  · rintro ⟨u, C, D, hD, rfl⟩
    have hDfix :
        Wedge2Formalization.N5SimplePoint.FixesPairBivector (k := k) D :=
      (Wedge2Formalization.Paper.N5.Row3.pointwise_stabilizer
        (k := k)
        D).2 hD
    exact
      (Wedge2Formalization.N6SimplePoint.fixesRadSimplePair_fromBlocks_zeroUpperRight_iff
        (k := k)
        (A := Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) u)
        (C := C)
        (D := D)).2 hDfix

theorem mem_K_shape_iff
    (g :
      Matrix Wedge2Formalization.N6SimplePoint.V
        Wedge2Formalization.N6SimplePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N6SimplePoint.W
            Wedge2Formalization.N6SimplePoint.I k,
          ∃ a x y p q : k,
            ∃ E : Matrix Wedge2Formalization.N5SimplePoint.I
                Wedge2Formalization.N5SimplePoint.I k,
              a ≠ 0 ∧
              E.det = 1 ∧
              g =
                Matrix.fromBlocks
                  (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) u)
                  0
                  C
                  (Matrix.fromBlocks
                    (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
                    (!![(0 : k), 0; 0, 0; p, q])
                    (!![a * (p * E 0 1 - q * E 0 0), 0, 0;
                       a * (p * E 1 1 - q * E 1 0), 0, 0])
                    E) := by
  constructor
  · rintro ⟨u, C, D, hD, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N5.Row3.mem_K_iff
        (k := k)
        D).1 hD with
      ⟨a, x, y, p, q, E, ha, hE, rfl⟩
    refine ⟨u, C, a, x, y, p, q, E, ha, hE, rfl⟩
  · rintro ⟨u, C, a, x, y, p, q, E, ha, hE, rfl⟩
    refine ⟨u, C, Matrix.fromBlocks
      (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
      (!![(0 : k), 0; 0, 0; p, q])
      (!![a * (p * E 0 1 - q * E 0 0), 0, 0;
         a * (p * E 1 1 - q * E 1 0), 0, 0])
      E, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N5.Row3.mem_K_iff
        (k := k)
        (Matrix.fromBlocks
          (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
          (!![(0 : k), 0; 0, 0; p, q])
          (!![a * (p * E 0 1 - q * E 0 0), 0, 0;
             a * (p * E 1 1 - q * E 1 0), 0, 0])
          E)).2 ⟨a, x, y, p, q, E, ha, hE, rfl⟩

theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6SimplePoint.V
        Wedge2Formalization.N6SimplePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N6SimplePoint.W
            Wedge2Formalization.N6SimplePoint.I k,
          ∃ D : Matrix Wedge2Formalization.N6SimplePoint.W
              Wedge2Formalization.N6SimplePoint.W k,
            D ∈ Wedge2Formalization.Paper.N5.Row3.K (k := k) ∧
              g =
                Matrix.fromBlocks
                  (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) u)
                  0
                  C
                  D := by
  rfl

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_9 \rtimes (G_m(k) \times G_m(k) \times SL_2(k))`, written in the
explicit local coordinates of this row rather than through a nested `n = 5`
core theorem. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6SimplePoint.V
        Wedge2Formalization.N6SimplePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N6SimplePoint.W
            Wedge2Formalization.N6SimplePoint.I k,
          ∃ a x y p q : k,
            ∃ E : Matrix Wedge2Formalization.N5SimplePoint.I
                Wedge2Formalization.N5SimplePoint.I k,
              a ≠ 0 ∧
              E.det = 1 ∧
              g =
                Matrix.fromBlocks
                  (Wedge2Formalization.N6SimplePoint.scalarBlock (k := k) u)
                  0
                  C
                  (Matrix.fromBlocks
                    (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
                    (!![(0 : k), 0; 0, 0; p, q])
                    (!![a * (p * E 0 1 - q * E 0 0), 0, 0;
                       a * (p * E 1 1 - q * E 1 0), 0, 0])
                    E) :=
  mem_K_shape_iff (k := k) (g := g)

theorem U_pointwise
    (C : Matrix Wedge2Formalization.N6SimplePoint.W
      Wedge2Formalization.N6SimplePoint.I k)
    (x y p q : k) :
    Wedge2Formalization.N6SimplePoint.FixesRadSimplePairBivector
      (U (k := k) C x y p q) := by
  simpa [U, Wedge2Formalization.Paper.N5.Row3.U] using
    Wedge2Formalization.N6SimplePoint.radSimple_pointwise_unipotent_family
      (k := k) (u := (1 : k)) (C := C) (x := x) (y := y) (p := p) (q := q)

theorem Levi_pointwise
    (u a : k)
    (ha : a ≠ 0)
    (E : Matrix Wedge2Formalization.N5SimplePoint.I
      Wedge2Formalization.N5SimplePoint.I k)
    (hE : E.det = 1) :
    Wedge2Formalization.N6SimplePoint.FixesRadSimplePairBivector
      (Levi (k := k) u a E) := by
  simpa [Levi] using
    Wedge2Formalization.N6SimplePoint.radSimple_pointwise_family
      (k := k) (u := u) (a := a) ha (C := 0) (E := E) hE

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6SimplePoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.simplePoint_borel_lift_action
        (k := k)
        (u := (1 : k))
        (a := a)
        (b := b)
        ha
        (C := 0)).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.N4PaperSummary.row1_borelCoeff] using
      (Wedge2Formalization.N6Summary.simplePoint_borel_lift_action
        (k := k)
        (u := (1 : k))
        (a := a)
        (b := b)
        ha
        (C := 0)).2

private def bottomProj :
    Matrix Wedge2Formalization.N6SimplePoint.W
      Wedge2Formalization.N6SimplePoint.V k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => if i = j then 1 else 0

private def bottomIncl :
    Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.W k :=
  fun i j =>
    match i with
    | Sum.inl _ => 0
    | Sum.inr i => if i = j then 1 else 0

private theorem bottomProj_mul_combo_mul_bottomIncl
    (a b : k) :
    bottomProj (k := k) * (a • rep₁ (k := k) + b • rep₂ (k := k)) * bottomIncl (k := k) =
      a • Wedge2Formalization.N5SimplePoint.rep₁ (k := k) +
        b • Wedge2Formalization.N5SimplePoint.rep₂ (k := k) := by
  ext i j
  change
    (bottomProj (k := k) * (a • rep₁ (k := k) + b • rep₂ (k := k)) * bottomIncl (k := k)) i j
      =
    (a • Wedge2Formalization.N5SimplePoint.rep₁ (k := k) +
      b • Wedge2Formalization.N5SimplePoint.rep₂ (k := k)) i j
  simp [bottomProj, bottomIncl, rep₁, rep₂,
    Wedge2Formalization.N6SimplePoint.radSimpleRep₁,
    Wedge2Formalization.N6SimplePoint.radSimpleRep₂,
    Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type, Matrix.add_apply,
    Matrix.smul_apply]

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • Wedge2Formalization.N5SimplePoint.rep₁ (k := k) +
          b • Wedge2Formalization.N5SimplePoint.rep₂ (k := k) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6SimplePoint.radSimpleRep₁,
      Wedge2Formalization.N6SimplePoint.radSimpleRep₂] using h'
  have h01 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h''
  have ha : a = 0 := by
    simpa [Wedge2Formalization.N5SimplePoint.rep₁,
      Wedge2Formalization.N5SimplePoint.rep₂, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h01
  have h12 := congrArg (fun M => M (Sum.inl 1) (Sum.inl 2)) h''
  have hb : b = 0 := by
    simpa [ha, Wedge2Formalization.N5SimplePoint.rep₁,
      Wedge2Formalization.N5SimplePoint.rep₂,
      Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N3PureSingular.ω23,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h12
  exact ⟨ha, hb⟩

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h12 := congrArg (fun M => M (Sum.inr (Sum.inl 1)) (Sum.inr (Sum.inl 2))) hzero
  simpa [rep₂, Wedge2Formalization.N6SimplePoint.radSimpleRep₂,
    Wedge2Formalization.N5SimplePoint.rep₂, Wedge2Formalization.N3PureSingular.ω23,
    Matrix.fromBlocks] using h12

private theorem rank_ge_four_of_leftCoeff_ne_zero_aux
    (a b : k)
    (ha : a ≠ 0) :
    4 ≤ (a • rep₁ (k := k) + b • rep₂ (k := k)).rank := by
  calc
    4 ≤ (a • Wedge2Formalization.N5SimplePoint.rep₁ (k := k) +
          b • Wedge2Formalization.N5SimplePoint.rep₂ (k := k)).rank := by
          exact Wedge2Formalization.Paper.N5.Row3.rank_ge_four_of_leftCoeff_ne_zero
            (k := k) a b ha
    _ = (bottomProj (k := k) * (a • rep₁ (k := k) + b • rep₂ (k := k)) *
          bottomIncl (k := k)).rank := by
          rw [bottomProj_mul_combo_mul_bottomIncl (k := k) (a := a) (b := b)]
    _ ≤ (a • rep₁ (k := k) + b • rep₂ (k := k)).rank := by
          exact (Matrix.rank_mul_le_left _ _).trans <| Matrix.rank_mul_le_right _ _

private def rows₂ : Fin 2 → Wedge2Formalization.N6SimplePoint.V
  | 0 => Sum.inr (Sum.inl 1)
  | 1 => Sum.inr (Sum.inl 2)

private def pick₂ :
    Matrix Wedge2Formalization.N6SimplePoint.V (Fin 2) k :=
  fun i j =>
    match i, j with
    | Sum.inr (Sum.inl 2), 0 => 1
    | Sum.inr (Sum.inl 1), 1 => -1
    | _, _ => 0

private theorem rows₂_mul_pick₂ :
    (rep₂ (k := k)).submatrix rows₂ (Equiv.refl _) * pick₂ (k := k) =
      (1 : Matrix (Fin 2) (Fin 2) k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rows₂, pick₂, rep₂, Wedge2Formalization.N6SimplePoint.radSimpleRep₂,
      Wedge2Formalization.N5SimplePoint.rep₂, Wedge2Formalization.N3PureSingular.ω23,
      Matrix.submatrix, Matrix.mul_apply, Matrix.fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two]

private def rep₂_left :
    Matrix Wedge2Formalization.N6SimplePoint.V (Fin 2) k :=
  fun i j =>
    match i, j with
    | Sum.inr (Sum.inl 1), 0 => 1
    | Sum.inr (Sum.inl 2), 1 => 1
    | _, _ => 0

private def rep₂_right :
    Matrix (Fin 2) Wedge2Formalization.N6SimplePoint.V k :=
  fun i j =>
    match i, j with
    | 0, Sum.inr (Sum.inl 2) => 1
    | 1, Sum.inr (Sum.inl 1) => -1
    | _, _ => 0

private theorem rep₂_factor :
    rep₂ (k := k) = rep₂_left (k := k) * rep₂_right (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N6SimplePoint.radSimpleRep₂,
        Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Matrix.fromBlocks,
        Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]
  ·
    fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N6SimplePoint.radSimpleRep₂,
        Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Matrix.fromBlocks,
        Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]
  ·
    rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N6SimplePoint.radSimpleRep₂,
        Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Matrix.fromBlocks,
        Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]
  ·
    rcases i with i | i <;> rcases j with j | j <;> fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N6SimplePoint.radSimpleRep₂,
        Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Matrix.fromBlocks,
        Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]

private theorem rep₂_rank_eq_two :
    (rep₂ (k := k)).rank = 2 := by
  have hle : (rep₂ (k := k)).rank ≤ 2 := by
    rw [rep₂_factor (k := k)]
    exact (Matrix.rank_mul_le_left (rep₂_left (k := k)) (rep₂_right (k := k))).trans <| by
      simpa using (Matrix.rank_le_card_width (rep₂_left (k := k)))
  have hge : 2 ≤ (rep₂ (k := k)).rank := by
    calc
      2 = (1 : Matrix (Fin 2) (Fin 2) k).rank := by
            simpa using (Matrix.rank_one (R := k) (n := Fin 2))
      _ =
          (((rep₂ (k := k)).submatrix rows₂ (Equiv.refl _)) * pick₂ (k := k)).rank := by
            rw [rows₂_mul_pick₂ (k := k)]
      _ ≤ ((rep₂ (k := k)).submatrix rows₂ (Equiv.refl _)).rank := by
            exact Matrix.rank_mul_le_left _ _
      _ ≤ (rep₂ (k := k)).rank := by
            exact Matrix.rank_submatrix_le rows₂ (Equiv.refl _) (rep₂ (k := k))
  exact le_antisymm hle hge

private theorem act_rank_eq
    (Ω : Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k)
    (g : Matrix Wedge2Formalization.N6SimplePoint.V
      Wedge2Formalization.N6SimplePoint.V k)
    (hg : Matrix.det g ≠ 0) :
    (Wedge2Formalization.N6SimplePoint.ActBivector Ω g).rank = Ω.rank := by
  have hg_unit : IsUnit (Matrix.det g) := isUnit_iff_ne_zero.mpr hg
  have hgt_unit : IsUnit (Matrix.det gᵀ) := by
    simpa [Matrix.det_transpose] using hg_unit
  calc
    (Wedge2Formalization.N6SimplePoint.ActBivector Ω g).rank = (g * Ω * gᵀ).rank := by rfl
    _ = (g * Ω).rank := by
          simpa [Wedge2Formalization.N6SimplePoint.ActBivector, Matrix.mul_assoc] using
            (Matrix.rank_mul_eq_left_of_isUnit_det (A := gᵀ) (B := g * Ω) hgt_unit)
    _ = Ω.rank := by
          simpa using
            (Matrix.rank_mul_eq_right_of_isUnit_det (A := g) (B := Ω) hg_unit)

private theorem act_rep₂_rank_eq_two
    (g :
      Matrix Wedge2Formalization.N6SimplePoint.V
        Wedge2Formalization.N6SimplePoint.V k)
    (hg : Matrix.det g ≠ 0) :
    (Wedge2Formalization.N6SimplePoint.ActBivector (rep₂ (k := k)) g).rank = 2 := by
  calc
    (Wedge2Formalization.N6SimplePoint.ActBivector (rep₂ (k := k)) g).rank
        = (rep₂ (k := k)).rank := act_rank_eq (k := k) (Ω := rep₂ (k := k)) (g := g) hg
    _ = 2 := rep₂_rank_eq_two (k := k)

/-- The rank obstruction on the embedded `[a]` quotient line again isolates the repeated
support point. -/
theorem rank_ge_four_of_leftCoeff_ne_zero
    (a b : k)
    (ha : a ≠ 0) :
    4 ≤ (a • rep₁ (k := k) + b • rep₂ (k := k)).rank :=
  rank_ge_four_of_leftCoeff_ne_zero_aux (k := k) a b ha

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6SimplePoint.V
        Wedge2Formalization.N6SimplePoint.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N6SimplePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6SimplePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      Wedge2Formalization.N6SimplePoint.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      Wedge2Formalization.N6SimplePoint.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have hγ : M 1 0 = 0 := by
    by_contra hγnz
    have hbound :
        4 ≤ (Wedge2Formalization.N6SimplePoint.ActBivector (rep₂ (k := k)) g).rank := by
      simpa [h2] using
        rank_ge_four_of_leftCoeff_ne_zero (k := k) (M 1 0) (M 1 1) hγnz
    rw [act_rep₂_rank_eq_two (k := k) (g := g) hg] at hbound
    omega
  have hδ : M 1 1 ≠ 0 := by
    intro hδ0
    have hzero :
        Wedge2Formalization.N6SimplePoint.ActBivector (rep₂ (k := k)) g = 0 := by
      simpa [hγ, hδ0] using h2
    have horig :=
      (actBivector_eq_zero_iff_of_det_ne_zero
        (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
    exact rep₂_ne_zero (k := k) horig
  have hα : M 0 0 ≠ 0 := by
    intro hα0
    have hzero :
        Wedge2Formalization.N6SimplePoint.ActBivector
          (M 1 1 • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k)) g = 0 := by
      calc
        Wedge2Formalization.N6SimplePoint.ActBivector
            (M 1 1 • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k)) g
            =
          M 1 1 • Wedge2Formalization.N6SimplePoint.ActBivector (rep₁ (k := k)) g +
            (-(M 0 1)) • Wedge2Formalization.N6SimplePoint.ActBivector (rep₂ (k := k)) g := by
              simp [Wedge2Formalization.N6SimplePoint.ActBivector, Matrix.mul_add, Matrix.add_mul]
        _ = 0 := by
              rw [h1, h2]
              ext i j
              simp [hα0, hγ, Matrix.add_apply, Matrix.smul_apply]
              ring
    have horig :=
      (actBivector_eq_zero_iff_of_det_ne_zero
        (k := k)
        (Ω := M 1 1 • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k))
        (g := g)
        hg).1 hzero
    have hcomb_ne :
        M 1 1 • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k) ≠ 0 := by
      intro hcomb
      rcases rep_pair_independent (k := k) hcomb with ⟨hδ', hβ'⟩
      exact hδ hδ'
    exact hcomb_ne horig
  refine ⟨M, ?_, hM⟩
  constructor
  · simpa [hγ]
  · have hdet :
        Matrix.det M = M 0 0 * M 1 1 := by
      simp [Matrix.det_fin_two, hγ]
    simpa [hdet] using mul_ne_zero hα hδ

end Row8

end N6
end Paper
end Wedge2Formalization
