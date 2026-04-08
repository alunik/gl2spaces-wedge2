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

/-! Appendix A, `n = 6`, row 4.
Representative `⟨e₁∧e₂ + e₃∧e₄ + e₅∧e₆, e₅∧e₆⟩`.
Divisor `2[a]+[b]`.
Claimed stabilizer:
`K_L = \Sp_4(k) \times \SL_2(k)`, exact projective quotient family
`Q_L = T`.
-/
namespace Row4

def rep₁ :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 4. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 4. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  rep₂ (k := k)

/-- The row-4 paper representative already agrees with the internal working pair. -/
def paperChange :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6WeightedTwoPoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6WeightedTwoPoint.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6WeightedTwoPoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6WeightedTwoPoint.ActBivector]

/-- The displayed block-diagonal pointwise family. -/
def Levi
    (A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
      Wedge2Formalization.N6WeightedTwoPoint.W k)
    (D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
      Wedge2Formalization.N6WeightedTwoPoint.I k) :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  Matrix.fromBlocks A 0 0 D

/-- Exact pointwise kernel for the weighted two-point row. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6WeightedTwoPoint.V
        Wedge2Formalization.N6WeightedTwoPoint.V k) :=
  { g |
      ∃ A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
          Wedge2Formalization.N6WeightedTwoPoint.W k,
        ∃ D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
            Wedge2Formalization.N6WeightedTwoPoint.I k,
          Wedge2Formalization.N4.FixesBivector
              (Wedge2Formalization.N4.splitRep₁ (k := k)) A ∧
            D.det = 1 ∧
            g = Levi (k := k) A D }

/-- The symplectic factor on the `2[a]+[b]` row, written as the stabilizer of the
symplectic form `e₁∧e₂ + e₃∧e₄`. -/
def Core :
    Set
      (Matrix Wedge2Formalization.N6WeightedTwoPoint.W
        Wedge2Formalization.N6WeightedTwoPoint.W k) :=
  { A | Wedge2Formalization.N4.FixesBivector (Wedge2Formalization.N4.splitRep₁ (k := k)) A }

/-- Public name for the `\Sp_4(k)` factor in Appendix A, row 4. -/
def Sp4Core :
    Set
      (Matrix Wedge2Formalization.N6WeightedTwoPoint.W
        Wedge2Formalization.N6WeightedTwoPoint.W k) :=
  Core (k := k)

/-- Exact coefficient-side torus family in the ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ M 0 0 + M 0 1 = M 1 1 ∧ Matrix.det M ≠ 0 }

/-- Chosen torus lift. -/
def lift
    (u v : k) :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  let A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k := !![u, 0; 0, 1]
  let H : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
      Wedge2Formalization.N6WeightedTwoPoint.W k := Matrix.fromBlocks A 0 0 A
  let E : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
      Wedge2Formalization.N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
  Matrix.fromBlocks H 0 0 E

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6WeightedTwoPoint.FixesPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N6Summary.weightedTwoPoint_pointwise_bivector_iff
        (k := k) (A := g.toBlocks₁₁) (B := g.toBlocks₁₂)
        (C := g.toBlocks₂₁) (D := g.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hg) with
      ⟨hB, hC, hA, hD⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₂, hA, hD, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Levi (k := k) g.toBlocks₁₁ g.toBlocks₂₂ := by
            simp [Levi, hB, hC]
  · rintro ⟨A, D, hA, hD, rfl⟩
    exact
      (Wedge2Formalization.N6Summary.weightedTwoPoint_pointwise_bivector_iff
        (k := k) (A := A) (B := 0) (C := 0) (D := D)).2 ⟨rfl, rfl, hA, hD⟩

theorem mem_K_shape_iff
    (g :
      Matrix Wedge2Formalization.N6WeightedTwoPoint.V
        Wedge2Formalization.N6WeightedTwoPoint.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
          Wedge2Formalization.N6WeightedTwoPoint.W k,
        ∃ D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
            Wedge2Formalization.N6WeightedTwoPoint.I k,
          Wedge2Formalization.N4.FixesBivector
              (Wedge2Formalization.N4.splitRep₁ (k := k)) A ∧
            D.det = 1 ∧
            g = Levi (k := k) A D := by
  rfl

theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6WeightedTwoPoint.V
        Wedge2Formalization.N6WeightedTwoPoint.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
          Wedge2Formalization.N6WeightedTwoPoint.W k,
        ∃ D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
            Wedge2Formalization.N6WeightedTwoPoint.I k,
          A ∈ Core (k := k) ∧
            D.det = 1 ∧
            g = Levi (k := k) A D := by
  constructor
  · rintro ⟨A, D, hA, hD, rfl⟩
    exact ⟨A, D, hA, hD, rfl⟩
  · rintro ⟨A, D, hA, hD, rfl⟩
    exact (mem_K_shape_iff (k := k) (g := Levi (k := k) A D)).2 ⟨A, D, hA, hD, rfl⟩

/-- Table-facing kernel statement for Appendix A, row 4:
`K_L = \Sp_4(k) \times \SL_2(k)`. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6WeightedTwoPoint.V
        Wedge2Formalization.N6WeightedTwoPoint.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
          Wedge2Formalization.N6WeightedTwoPoint.W k,
        ∃ D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
            Wedge2Formalization.N6WeightedTwoPoint.I k,
          A ∈ Sp4Core (k := k) ∧
            D.det = 1 ∧
            g = Levi (k := k) A D :=
  mem_K_iff (k := k) (g := g)

theorem Levi_pointwise
    (A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
      Wedge2Formalization.N6WeightedTwoPoint.W k)
    (D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
      Wedge2Formalization.N6WeightedTwoPoint.I k)
    (hA : Wedge2Formalization.N4.FixesBivector (Wedge2Formalization.N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1) :
    Wedge2Formalization.N6WeightedTwoPoint.FixesPairBivector (Levi (k := k) A D) := by
  simpa [Levi] using
    Wedge2Formalization.N6Summary.weightedTwoPoint_pointwise_family
      (k := k) (A := A) (D := D) hA hD

theorem coeff_mem_Qproj
    (u v : k)
    (hu : u ≠ 0)
    (hv : v ≠ 0) :
    Wedge2Formalization.N4PaperSummary.row2_torusCoeff (k := k) u v ∈ Qproj (k := k) := by
  constructor
  · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff]
  · constructor
    · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff]
    · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, Matrix.det_fin_two, hu, hv]

theorem quotient_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6WeightedTwoPoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_torusCoeff (k := k) u v) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.weightedTwoPoint_torus_lift_action
        (k := k) (u := u) (v := v)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.N4PaperSummary.row2_torusCoeff] using
      (Wedge2Formalization.N6Summary.weightedTwoPoint_torus_lift_action
        (k := k) (u := u) (v := v)).2

private def rep₂_left :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.I k :=
  fun i j =>
    match i, j with
    | Sum.inr i, j => if i = j then 1 else 0
    | _, _ => 0

private def rep₂_right :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.I
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => Wedge2Formalization.N4.J (k := k) i j

private theorem rep₂_factor :
    rep₂ (k := k) = rep₂_left (k := k) * rep₂_right (k := k) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₂,
      rep₂_left, rep₂_right, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.mul_apply, Fin.sum_univ_two]

private def bottomProj :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.I
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => if i = j then 1 else 0

private def bottomIncl :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.I k :=
  fun i j =>
    match i with
    | Sum.inl _ => 0
    | Sum.inr i => if i = j then 1 else 0

private theorem bottomProj_mul_rep₂_mul_bottomIncl :
    bottomProj (k := k) * rep₂ (k := k) * bottomIncl (k := k) =
      Wedge2Formalization.N4.J (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simp [bottomProj, bottomIncl, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₂,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply,
      Fintype.sum_sum_type, Fin.sum_univ_four, Fin.sum_univ_two]
  ·
    simp [bottomProj, bottomIncl, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₂,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply,
      Fintype.sum_sum_type, Fin.sum_univ_four, Fin.sum_univ_two]
  ·
    simp [bottomProj, bottomIncl, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₂,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply,
      Fintype.sum_sum_type, Fin.sum_univ_four, Fin.sum_univ_two]
  ·
    simp [bottomProj, bottomIncl, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₂,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply,
      Fintype.sum_sum_type, Fin.sum_univ_four, Fin.sum_univ_two]

private theorem rep₂_rank_le_two :
    (rep₂ (k := k)).rank ≤ 2 := by
  rw [rep₂_factor (k := k)]
  exact (Matrix.rank_mul_le_left (rep₂_left (k := k)) (rep₂_right (k := k))).trans <| by
    simpa using (Matrix.rank_le_card_width (rep₂_left (k := k)))

private theorem rep₂_rank_ge_two :
    2 ≤ (rep₂ (k := k)).rank := by
  calc
    2 = (Wedge2Formalization.N4.J (k := k)).rank := by
      have hUnit : IsUnit (Matrix.det (Wedge2Formalization.N4.J (k := k))) := by
        simpa [Wedge2Formalization.N4.J_det] using (isUnit_one : IsUnit (1 : k))
      simpa [Wedge2Formalization.N6WeightedTwoPoint.I] using
        (Matrix.rank_of_isUnit
          (Wedge2Formalization.N4.J (k := k))
          ((Wedge2Formalization.N4.J (k := k)).isUnit_iff_isUnit_det.mpr hUnit)).symm
    _ = (bottomProj (k := k) * rep₂ (k := k) * bottomIncl (k := k)).rank := by
      rw [bottomProj_mul_rep₂_mul_bottomIncl (k := k)]
    _ ≤ (rep₂ (k := k)).rank := by
      exact (Matrix.rank_mul_le_left _ _).trans <| Matrix.rank_mul_le_right _ _

private theorem rep₂_rank_eq_two :
    (rep₂ (k := k)).rank = 2 := by
  exact le_antisymm rep₂_rank_le_two (rep₂_rank_ge_two (k := k))

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have hcoord := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) hzero
  simpa [rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₂,
    Wedge2Formalization.N4.J, Matrix.fromBlocks] using hcoord

private def topProj :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.W
      Wedge2Formalization.N6WeightedTwoPoint.V k :=
  fun i j =>
    match j with
    | Sum.inl j => if i = j then 1 else 0
    | Sum.inr _ => 0

private def topIncl :
    Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.W k :=
  fun i j =>
    match i with
    | Sum.inl i => if i = j then 1 else 0
    | Sum.inr _ => 0

private theorem topProj_mul_topRep_mul_topIncl :
    topProj (k := k) * (rep₁ (k := k) - rep₂ (k := k)) * topIncl (k := k) =
      Wedge2Formalization.N4.splitRep₁ (k := k) := by
  ext i j
  cases i with
  | inl i =>
      cases j with
      | inl j =>
          fin_cases i <;> fin_cases j <;>
            simp [topProj, topIncl, rep₁, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₁,
              Wedge2Formalization.N6WeightedTwoPoint.rep₂, Wedge2Formalization.N4.splitRep₁,
              Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
              Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type]
      | inr j =>
          fin_cases i <;> fin_cases j <;>
            simp [topProj, topIncl, rep₁, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₁,
              Wedge2Formalization.N6WeightedTwoPoint.rep₂, Wedge2Formalization.N4.splitRep₁,
              Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
              Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type]
  | inr i =>
      cases j with
      | inl j =>
          fin_cases i <;> fin_cases j <;>
            simp [topProj, topIncl, rep₁, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₁,
              Wedge2Formalization.N6WeightedTwoPoint.rep₂, Wedge2Formalization.N4.splitRep₁,
              Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
              Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type]
      | inr j =>
          fin_cases i <;> fin_cases j <;>
            simp [topProj, topIncl, rep₁, rep₂, Wedge2Formalization.N6WeightedTwoPoint.rep₁,
              Wedge2Formalization.N6WeightedTwoPoint.rep₂, Wedge2Formalization.N4.splitRep₁,
              Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
              Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type]

private theorem topCombo_rank_ge_four
    (a : k)
    (ha : a ≠ 0) :
    4 ≤ (a • (rep₁ (k := k) - rep₂ (k := k))).rank := by
  have hUnit : IsUnit (Matrix.det (a • Wedge2Formalization.N4.splitRep₁ (k := k))) := by
    refine isUnit_iff_ne_zero.mpr ?_
    rw [Matrix.det_smul, Wedge2Formalization.N6WeightedTwoPoint.det_splitRep₁]
    simp [ha]
  calc
    4 = (a • Wedge2Formalization.N4.splitRep₁ (k := k)).rank := by
      simpa [Wedge2Formalization.N6WeightedTwoPoint.W] using
        (Matrix.rank_of_isUnit
          (a • Wedge2Formalization.N4.splitRep₁ (k := k))
          ((a • Wedge2Formalization.N4.splitRep₁ (k := k)).isUnit_iff_isUnit_det.mpr hUnit)).symm
    _ = (topProj (k := k) * (a • (rep₁ (k := k) - rep₂ (k := k))) * topIncl (k := k)).rank := by
      have hEq :
          topProj (k := k) * (a • (rep₁ (k := k) - rep₂ (k := k))) * topIncl (k := k) =
            a • Wedge2Formalization.N4.splitRep₁ (k := k) := by
        calc
          topProj (k := k) * (a • (rep₁ (k := k) - rep₂ (k := k))) * topIncl (k := k)
              = a • (topProj (k := k) * (rep₁ (k := k) - rep₂ (k := k)) * topIncl (k := k)) := by
                  simp [Matrix.mul_assoc]
          _ = a • Wedge2Formalization.N4.splitRep₁ (k := k) := by
                rw [topProj_mul_topRep_mul_topIncl (k := k)]
      rw [hEq]
    _ ≤ (a • (rep₁ (k := k) - rep₂ (k := k))).rank := by
      exact (Matrix.rank_mul_le_left _ _).trans <| Matrix.rank_mul_le_right _ _

private theorem act_rank_eq
    (Ω : Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k)
    (g : Matrix Wedge2Formalization.N6WeightedTwoPoint.V
      Wedge2Formalization.N6WeightedTwoPoint.V k)
    (hg : Matrix.det g ≠ 0) :
    (Wedge2Formalization.N6WeightedTwoPoint.ActBivector Ω g).rank = Ω.rank := by
  have hg_unit : IsUnit (Matrix.det g) := isUnit_iff_ne_zero.mpr hg
  have hgt_unit : IsUnit (Matrix.det gᵀ) := by
    simpa [Matrix.det_transpose] using hg_unit
  calc
    (Wedge2Formalization.N6WeightedTwoPoint.ActBivector Ω g).rank
        = (g * Ω * gᵀ).rank := by rfl
    _ = (g * Ω).rank := by
          simpa [Wedge2Formalization.N6WeightedTwoPoint.ActBivector, Matrix.mul_assoc] using
            (Matrix.rank_mul_eq_left_of_isUnit_det (A := gᵀ) (B := g * Ω) hgt_unit)
    _ = Ω.rank := by
          simpa using
            (Matrix.rank_mul_eq_right_of_isUnit_det (A := g) (B := Ω) hg_unit)

private theorem combo_rank_le_two_of_leftCoeff_zero
    (c : k) :
    ((0 : k) • rep₁ (k := k) + c • rep₂ (k := k)).rank ≤ 2 := by
  rw [zero_smul, zero_add]
  have hfactor :
      c • rep₂ (k := k) = (c • rep₂_left (k := k)) * rep₂_right (k := k) := by
    rw [rep₂_factor (k := k)]
    ext i j
    simp [Matrix.mul_apply, Fin.sum_univ_two]
    ring
  rw [hfactor]
  exact (Matrix.rank_mul_le_left (c • rep₂_left (k := k)) (rep₂_right (k := k))).trans <| by
    simpa using (Matrix.rank_le_card_width (c • rep₂_left (k := k)))

theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • rep₁ (k := k) + b • rep₂ (k := k)) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  Wedge2Formalization.N6Summary.weightedTwoPoint_det_zero_iff (k := k) (a := a) (b := b)

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6WeightedTwoPoint.V
        Wedge2Formalization.N6WeightedTwoPoint.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N6WeightedTwoPoint.PreservesSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6WeightedTwoPoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N6WeightedTwoPoint.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N6WeightedTwoPoint.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  have hγ : γ = 0 := by
    by_contra hγnz
    have hrep₂_det :
        Matrix.det (rep₂ (k := k)) = 0 := by
      simpa using
        ((det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl))
    have hrep₂_det_zero :
        Matrix.det (γ • rep₁ (k := k) + δ • rep₂ (k := k)) = 0 := by
      rw [← h2', Wedge2Formalization.N6WeightedTwoPoint.ActBivector, Matrix.det_mul,
        Matrix.det_mul, Matrix.det_transpose, hrep₂_det]
      simp
    have hsum : γ + δ = 0 := by
      rcases (det_zero_iff (k := k) (a := γ) (b := δ)).1 hrep₂_det_zero with hγ0 | hsum
      · exact (hγnz hγ0).elim
      · exact hsum
    have hlarge :
        4 ≤ (Wedge2Formalization.N6WeightedTwoPoint.ActBivector (rep₂ (k := k)) g).rank := by
      have hrewrite :
          γ • rep₁ (k := k) + δ • rep₂ (k := k) =
            γ • (rep₁ (k := k) - rep₂ (k := k)) := by
        have hδeq : δ = -γ := by
          rw [eq_neg_iff_add_eq_zero]
          simpa [add_comm] using hsum
        ext i j
        simp [hδeq, sub_eq_add_neg]
        ring
      rw [h2', hrewrite]
      exact topCombo_rank_ge_four (k := k) γ hγnz
    have hsmall :
        (Wedge2Formalization.N6WeightedTwoPoint.ActBivector (rep₂ (k := k)) g).rank = 2 := by
      rw [act_rank_eq (k := k) (Ω := rep₂ (k := k)) (g := g) hg]
      exact rep₂_rank_eq_two (k := k)
    omega
  have hδ : δ ≠ 0 := by
    intro hδ0
    have hzero :
        Wedge2Formalization.N6WeightedTwoPoint.ActBivector (rep₂ (k := k)) g = 0 := by
      simpa [hγ, hδ0] using h2'
    have horig :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
    exact rep₂_ne_zero (k := k) horig
  have htop_action :
      Wedge2Formalization.N6WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g =
        α • rep₁ (k := k) + (β - δ) • rep₂ (k := k) := by
    calc
      Wedge2Formalization.N6WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g
          =
        Wedge2Formalization.N6WeightedTwoPoint.ActBivector
          (rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) g := by
            simp [sub_eq_add_neg]
      _ =
        Wedge2Formalization.N6WeightedTwoPoint.ActBivector (rep₁ (k := k)) g +
          (-1 : k) • Wedge2Formalization.N6WeightedTwoPoint.ActBivector (rep₂ (k := k)) g := by
            rw [Wedge2Formalization.N6WeightedTwoPoint.actBivector_add,
              Wedge2Formalization.N6WeightedTwoPoint.actBivector_smul]
      _ = (α • rep₁ (k := k) + β • rep₂ (k := k)) +
            (-1 : k) • (γ • rep₁ (k := k) + δ • rep₂ (k := k)) := by
              rw [h1', h2']
      _ = α • rep₁ (k := k) + (β - δ) • rep₂ (k := k) := by
            ext i j
            simp [hγ, sub_eq_add_neg]
            ring
  have hα : α ≠ 0 := by
    intro hα0
    have hlarge :
        4 ≤ (Wedge2Formalization.N6WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g).rank := by
      rw [act_rank_eq (k := k) (Ω := rep₁ (k := k) - rep₂ (k := k)) (g := g) hg]
      simpa using topCombo_rank_ge_four (k := k) (1 : k) one_ne_zero
    have hsmall :
        (Wedge2Formalization.N6WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g).rank ≤ 2 := by
      rw [htop_action, hα0]
      simpa using combo_rank_le_two_of_leftCoeff_zero (k := k) (β - δ)
    omega
  have halpha_plus : α + β = δ := by
    have htop_det_zero :
        Matrix.det (rep₁ (k := k) - rep₂ (k := k)) = 0 := by
      simpa [sub_eq_add_neg] using
        ((det_zero_iff (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (by ring)))
    have hcombo_zero :
        Matrix.det (α • rep₁ (k := k) + (β - δ) • rep₂ (k := k)) = 0 := by
      rw [← htop_action, Wedge2Formalization.N6WeightedTwoPoint.ActBivector, Matrix.det_mul,
        Matrix.det_mul, Matrix.det_transpose, htop_det_zero]
      simp
    rcases (det_zero_iff (k := k) (a := α) (b := β - δ)).1 hcombo_zero with hα0 | hsum0
    · exact False.elim (hα hα0)
    · calc
        α + β = α + (β - δ) + δ := by ring
        _ = 0 + δ := by rw [hsum0]
        _ = δ := by ring
  refine ⟨!![α, β; γ, δ], ?_, ?_⟩
  · constructor
    · simpa [hγ]
    · constructor
      · simpa using halpha_plus
      · simpa [Matrix.det_fin_two, hγ] using mul_ne_zero hα hδ
  · exact ⟨h1', h2'⟩

end Row4

end N6
end Paper
end Wedge2Formalization
