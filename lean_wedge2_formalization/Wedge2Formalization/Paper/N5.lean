import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N4
import Wedge2Formalization.N5Summary
import Wedge2Formalization.N5PaperSummary
import Wedge2Formalization.N5OnePoint
import Wedge2Formalization.N5PureSingular
import Wedge2Formalization.N5SimplePoint
import Wedge2Formalization.N5PureSingularLong
import Mathlib.LinearAlgebra.Matrix.Rank

namespace Wedge2Formalization
namespace Paper
namespace N5

open Matrix

variable {k : Type*} [Field k]

/-!
This module exposes the paper-facing `n = 5` rows whose kernel and quotient
calculations are complete in the internal development. The public surface is
organized row-by-row in the appendix order, with each row namespace exposing
the representative pair, the named kernel family, the exact coefficient-side
quotient family, and the main pointwise / quotient theorems.
-/

/-! Appendix A, `n = 5`, row 1.
Representative `⟨e₁∧e₄ + e₂∧e₅, e₂∧e₄ + e₃∧e₅⟩`.
Divisor `0`.
Claimed stabilizer:
`K_L = U_4 \rtimes G_m(k)`, exact quotient family `Q_L = GL_2(k)`.
-/
namespace Row1

def rep₁ :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  Wedge2Formalization.N5PureSingularLong.rep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  Wedge2Formalization.N5PureSingularLong.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 1. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 1. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  rep₂ (k := k)

/-- The row-1 paper representative already agrees with the internal working pair. -/
def paperChange :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N5PureSingularLong.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N5PureSingularLong.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N5PureSingularLong.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N5PureSingularLong.ActBivector]

/-- The displayed unipotent family `U_4`. -/
def U
    (u v w z : k) :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N5PureSingularLong.W Wedge2Formalization.N5PureSingularLong.W k)
    0
    (Wedge2Formalization.N5PureSingularLong.lowerHankel (k := k) u v w z)
    (1 : Matrix Wedge2Formalization.N5PureSingularLong.I Wedge2Formalization.N5PureSingularLong.I k)

/-- The displayed Levi torus `G_m`. -/
def Levi
    (a : k) :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  Matrix.fromBlocks
    (a • (1 : Matrix Wedge2Formalization.N5PureSingularLong.W Wedge2Formalization.N5PureSingularLong.W k))
    0
    0
    (a⁻¹ • (1 : Matrix Wedge2Formalization.N5PureSingularLong.I Wedge2Formalization.N5PureSingularLong.I k))

/-- Exact pointwise kernel family for the pure singular length-two row, written in the
displayed paper form `U_4 ⋊ G_m`. -/
def K :
    Set (Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k) :=
  { g | ∃ a u v w z : k, a ≠ 0 ∧ g = U (k := k) u v w z * Levi (k := k) a }

/-- Exact coefficient-side projective quotient family `GL_2(k)` in the ordered basis
`(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

/-- The coefficient matrix produced by the chosen determinant-scaled `GL_2` lift. -/
def coeff
    (a b c d : k) :
    Matrix (Fin 2) (Fin 2) k :=
  !![(Wedge2Formalization.N5PureSingularLong.Delta a b c d) * a,
     (Wedge2Formalization.N5PureSingularLong.Delta a b c d) * c;
     (Wedge2Formalization.N5PureSingularLong.Delta a b c d) * b,
     (Wedge2Formalization.N5PureSingularLong.Delta a b c d) * d]

/-- Chosen determinant-scaled `GL_2(k)` lift. -/
def lift
    (a b c d : k) :
    Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5PureSingularLong.sym2Raw (k := k) a b c d)
    0
    0
    (Wedge2Formalization.N5PureSingularLong.adjointT (k := k) a b c d)

/-- The raw pointwise matrix shape on the pure singular length-two row factors as the
displayed `U_4 ⋊ G_m` product. -/
theorem shape_eq_product
    (a u v w z : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N5PureSingularLong.pointwiseShape (k := k) a u v w z =
      U (k := k) (a⁻¹ * u) (a⁻¹ * v) (a⁻¹ * w) (a⁻¹ * z) * Levi (k := k) a := by
  rw [U, Levi, Wedge2Formalization.N5PureSingularLong.pointwiseShape, Matrix.fromBlocks_multiply]
  ext i j
  cases i with
  | inl i =>
      cases j with
      | inl j =>
          fin_cases i <;> fin_cases j <;>
            simp [Wedge2Formalization.N5PureSingularLong.lowerHankel,
              Matrix.fromBlocks, Matrix.smul_apply, ha]
      | inr j =>
          fin_cases i <;> fin_cases j <;>
            simp [Wedge2Formalization.N5PureSingularLong.lowerHankel,
              Matrix.fromBlocks, Matrix.smul_apply, ha]
  | inr i =>
      cases j with
      | inl j =>
          fin_cases i <;> fin_cases j <;>
            simp [Wedge2Formalization.N5PureSingularLong.lowerHankel,
              Matrix.fromBlocks, Matrix.smul_apply, ha]
      | inr j =>
          fin_cases i <;> fin_cases j <;>
            simp [Wedge2Formalization.N5PureSingularLong.lowerHankel,
              Matrix.fromBlocks, Matrix.smul_apply, ha]

theorem U_pointwise
    (u v w z : k) :
    Wedge2Formalization.N5PureSingularLong.FixesPairBivector (U (k := k) u v w z) := by
  simpa [U] using
    Wedge2Formalization.N5Summary.pureSingularLong_lowerHankel_pointwise_family
      (k := k) u v w z

theorem Levi_pointwise
    (a : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N5PureSingularLong.FixesPairBivector (Levi (k := k) a) := by
  simpa [Levi, inv_inv] using
    Wedge2Formalization.N5Summary.pureSingularLong_pointwise_scale_family
      (k := k) (t := a⁻¹) (inv_ne_zero ha)

private def radVecRep₁ : Wedge2Formalization.N5PureSingularLong.V → k
  | Sum.inl 2 => 1
  | _ => 0

private def radVecRep₂ : Wedge2Formalization.N5PureSingularLong.V → k
  | Sum.inl 0 => 1
  | _ => 0

private def radVecRep12 : Wedge2Formalization.N5PureSingularLong.V → k
  | Sum.inl 0 => 1
  | Sum.inl 1 => -1
  | Sum.inl 2 => 1
  | _ => 0

private theorem rep₁_mulVec_radVecRep₁ :
    rep₁ (k := k) *ᵥ radVecRep₁ (k := k) = 0 := by
  funext i
  cases i with
  | inl i =>
      fin_cases i <;>
        simp [rep₁, Wedge2Formalization.N5PureSingularLong.rep₁,
          radVecRep₁, Wedge2Formalization.N5PureSingularLong.mulX, Matrix.mulVec,
          Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]
  | inr i =>
      fin_cases i <;>
        simp [rep₁, Wedge2Formalization.N5PureSingularLong.rep₁,
          radVecRep₁, Wedge2Formalization.N5PureSingularLong.mulX, Matrix.mulVec,
          Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]

private theorem rep₂_mulVec_radVecRep₂ :
    rep₂ (k := k) *ᵥ radVecRep₂ (k := k) = 0 := by
  funext i
  cases i with
  | inl i =>
      fin_cases i <;>
        simp [rep₂, Wedge2Formalization.N5PureSingularLong.rep₂,
          radVecRep₂, Wedge2Formalization.N5PureSingularLong.mulY, Matrix.mulVec,
          Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]
  | inr i =>
      fin_cases i <;>
        simp [rep₂, Wedge2Formalization.N5PureSingularLong.rep₂,
          radVecRep₂, Wedge2Formalization.N5PureSingularLong.mulY, Matrix.mulVec,
          Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]

private theorem rep12_mulVec_radVecRep12 :
    (rep₁ (k := k) + rep₂ (k := k)) *ᵥ radVecRep12 (k := k) = 0 := by
  funext i
  cases i with
  | inl i =>
      fin_cases i <;>
        simp [rep₁, rep₂, Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.rep₂, radVecRep12,
          Wedge2Formalization.N5PureSingularLong.mulX, Wedge2Formalization.N5PureSingularLong.mulY,
          Matrix.mulVec, Matrix.fromBlocks, Matrix.add_apply, dotProduct, Fintype.sum_sum_type,
          Fin.sum_univ_three, Fin.sum_univ_two] <;>
        ring
  | inr i =>
      fin_cases i <;>
        simp [rep₁, rep₂, Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.rep₂, radVecRep12,
          Wedge2Formalization.N5PureSingularLong.mulX, Wedge2Formalization.N5PureSingularLong.mulY,
          Matrix.mulVec, Matrix.fromBlocks, Matrix.add_apply, dotProduct, Fintype.sum_sum_type,
          Fin.sum_univ_three, Fin.sum_univ_two] <;>
        ring

private theorem upperRight_zero_of_pointwise
    (g : Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k)
    (hg : Matrix.det g ≠ 0)
    (hfix : Wedge2Formalization.N5PureSingularLong.FixesPairBivector (k := k) g) :
    g.toBlocks₁₂ = 0 := by
  rcases hfix with ⟨h₁, h₂⟩
  have h₁' :
      Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k)) g = rep₁ (k := k) := by
    simpa [rep₁] using h₁
  have h₂' :
      Wedge2Formalization.N5PureSingularLong.ActBivector (rep₂ (k := k)) g = rep₂ (k := k) := by
    simpa [rep₂] using h₂
  have h12 :
      Wedge2Formalization.N5PureSingularLong.ActBivector
          (rep₁ (k := k) + rep₂ (k := k)) g =
        rep₁ (k := k) + rep₂ (k := k) := by
    calc
      Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k) + rep₂ (k := k)) g
          =
        Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k)) g +
          Wedge2Formalization.N5PureSingularLong.ActBivector (rep₂ (k := k)) g := by
            rw [Wedge2Formalization.N5PureSingularLong.actBivector_add]
      _ = rep₁ (k := k) + rep₂ (k := k) := by
            simpa [rep₁, rep₂] using congrArg₂ HAdd.hAdd h₁ h₂
  have hpre₁ :
      rep₁ (k := k) *ᵥ (gᵀ *ᵥ radVecRep₁ (k := k)) = 0 := by
    have hzero :
        Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k)) g *ᵥ
            radVecRep₁ (k := k) = 0 := by
      rw [h₁']
      exact rep₁_mulVec_radVecRep₁ (k := k)
    have hzero' :
        g *ᵥ (rep₁ (k := k) *ᵥ (gᵀ *ᵥ radVecRep₁ (k := k))) = 0 := by
      simpa [Wedge2Formalization.N5PureSingularLong.ActBivector, Matrix.mulVec_mulVec,
        Matrix.mul_assoc] using hzero
    exact Matrix.eq_zero_of_mulVec_eq_zero hg hzero'
  have hpre₂ :
      rep₂ (k := k) *ᵥ (gᵀ *ᵥ radVecRep₂ (k := k)) = 0 := by
    have hzero :
        Wedge2Formalization.N5PureSingularLong.ActBivector (rep₂ (k := k)) g *ᵥ
            radVecRep₂ (k := k) = 0 := by
      rw [h₂']
      exact rep₂_mulVec_radVecRep₂ (k := k)
    have hzero' :
        g *ᵥ (rep₂ (k := k) *ᵥ (gᵀ *ᵥ radVecRep₂ (k := k))) = 0 := by
      simpa [Wedge2Formalization.N5PureSingularLong.ActBivector, Matrix.mulVec_mulVec,
        Matrix.mul_assoc] using hzero
    exact Matrix.eq_zero_of_mulVec_eq_zero hg hzero'
  have hpre₁₂ :
      (rep₁ (k := k) + rep₂ (k := k)) *ᵥ (gᵀ *ᵥ radVecRep12 (k := k)) = 0 := by
    have hzero :
        Wedge2Formalization.N5PureSingularLong.ActBivector
            (rep₁ (k := k) + rep₂ (k := k)) g *ᵥ radVecRep12 (k := k) = 0 := by
      rw [h12]
      exact rep12_mulVec_radVecRep12 (k := k)
    have hzero' :
        g *ᵥ ((rep₁ (k := k) + rep₂ (k := k)) *ᵥ (gᵀ *ᵥ radVecRep12 (k := k))) = 0 := by
      simpa [Wedge2Formalization.N5PureSingularLong.ActBivector, Matrix.mulVec_mulVec,
        Matrix.mul_assoc] using hzero
    exact Matrix.eq_zero_of_mulVec_eq_zero hg hzero'
  have hB20 : g.toBlocks₁₂ 2 0 = 0 := by
    have hcoord := congrArg (fun x => x (Sum.inl 0)) hpre₁
    have hcoord' :
        (gᵀ *ᵥ radVecRep₁ (k := k)) (Sum.inr 0) = 0 := by
      simpa [rep₁, Wedge2Formalization.N5PureSingularLong.rep₁, Wedge2Formalization.N5PureSingularLong.mulX,
        Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
        Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
    simpa [Matrix.mulVec, Matrix.transpose_apply, radVecRep₁, Matrix.toBlocks₁₂,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] using hcoord'
  have hB21 : g.toBlocks₁₂ 2 1 = 0 := by
    have hcoord := congrArg (fun x => x (Sum.inl 1)) hpre₁
    have hcoord' :
        (gᵀ *ᵥ radVecRep₁ (k := k)) (Sum.inr 1) = 0 := by
      simpa [rep₁, Wedge2Formalization.N5PureSingularLong.rep₁, Wedge2Formalization.N5PureSingularLong.mulX,
        Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
        Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
    simpa [Matrix.mulVec, Matrix.transpose_apply, radVecRep₁, Matrix.toBlocks₁₂,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] using hcoord'
  have hB00 : g.toBlocks₁₂ 0 0 = 0 := by
    have hcoord := congrArg (fun x => x (Sum.inl 1)) hpre₂
    have hcoord' :
        (gᵀ *ᵥ radVecRep₂ (k := k)) (Sum.inr 0) = 0 := by
      simpa [rep₂, Wedge2Formalization.N5PureSingularLong.rep₂, Wedge2Formalization.N5PureSingularLong.mulY,
        Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
        Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
    simpa [Matrix.mulVec, Matrix.transpose_apply, radVecRep₂, Matrix.toBlocks₁₂,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] using hcoord'
  have hB01 : g.toBlocks₁₂ 0 1 = 0 := by
    have hcoord := congrArg (fun x => x (Sum.inl 2)) hpre₂
    have hcoord' :
        (gᵀ *ᵥ radVecRep₂ (k := k)) (Sum.inr 1) = 0 := by
      simpa [rep₂, Wedge2Formalization.N5PureSingularLong.rep₂, Wedge2Formalization.N5PureSingularLong.mulY,
        Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
        Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
    simpa [Matrix.mulVec, Matrix.transpose_apply, radVecRep₂, Matrix.toBlocks₁₂,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] using hcoord'
  have hrow10 :
      g.toBlocks₁₂ 0 0 - g.toBlocks₁₂ 1 0 + g.toBlocks₁₂ 2 0 = 0 := by
    have hcoord := congrArg (fun x => x (Sum.inl 0)) hpre₁₂
    have hcoord' :
        (gᵀ *ᵥ radVecRep12 (k := k)) (Sum.inr 0) = 0 := by
      simpa [rep₁, rep₂, Wedge2Formalization.N5PureSingularLong.rep₁,
        Wedge2Formalization.N5PureSingularLong.rep₂, Wedge2Formalization.N5PureSingularLong.mulX,
        Wedge2Formalization.N5PureSingularLong.mulY, Matrix.mulVec, Matrix.fromBlocks,
        Matrix.add_apply, dotProduct, Fintype.sum_sum_type,
        Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
    simpa [Matrix.mulVec, Matrix.transpose_apply, radVecRep12, Matrix.toBlocks₁₂,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
      sub_eq_add_neg, add_assoc] using hcoord'
  have hrow11 :
      g.toBlocks₁₂ 0 1 - g.toBlocks₁₂ 1 1 + g.toBlocks₁₂ 2 1 = 0 := by
    have hcoord := congrArg (fun x => x (Sum.inl 2)) hpre₁₂
    have hcoord' :
        (gᵀ *ᵥ radVecRep12 (k := k)) (Sum.inr 1) = 0 := by
      simpa [rep₁, rep₂, Wedge2Formalization.N5PureSingularLong.rep₁,
        Wedge2Formalization.N5PureSingularLong.rep₂, Wedge2Formalization.N5PureSingularLong.mulX,
        Wedge2Formalization.N5PureSingularLong.mulY, Matrix.mulVec, Matrix.fromBlocks,
        Matrix.add_apply, dotProduct, Fintype.sum_sum_type,
        Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
    simpa [Matrix.mulVec, Matrix.transpose_apply, radVecRep12, Matrix.toBlocks₁₂,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
      sub_eq_add_neg, add_assoc] using hcoord'
  have hB10 : g.toBlocks₁₂ 1 0 = 0 := by
    have : -g.toBlocks₁₂ 1 0 = 0 := by simpa [hB00, hB20] using hrow10
    simpa using neg_eq_zero.mp this
  have hB11 : g.toBlocks₁₂ 1 1 = 0 := by
    have : -g.toBlocks₁₂ 1 1 = 0 := by simpa [hB01, hB21] using hrow11
    simpa using neg_eq_zero.mp this
  ext i j
  fin_cases i <;> fin_cases j <;> simp [hB00, hB01, hB10, hB11, hB20, hB21]

theorem pointwise_stabilizer
    (g : Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k)
    (hg : Matrix.det g ≠ 0) :
    Wedge2Formalization.N5PureSingularLong.FixesPairBivector (k := k) g ↔
      g ∈ K (k := k) := by
  constructor
  · intro hfix
    have hB : g.toBlocks₁₂ = 0 := upperRight_zero_of_pointwise (k := k) g hg hfix
    have hfix' :
        Wedge2Formalization.N5PureSingularLong.FixesPairBivector
          (Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂) := by
      have hEq : Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ = g := by
        calc
          Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂
              = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
                  simpa [hB]
          _ = g := Matrix.fromBlocks_toBlocks g
      simpa [hEq] using hfix
    rcases
      (Wedge2Formalization.N5PureSingularLong.fixesPair_fromBlocks_zeroUpperRight_iff_shape
        (k := k) (A := g.toBlocks₁₁) (C := g.toBlocks₂₁) (D := g.toBlocks₂₂)).1 hfix' with
      ⟨ha, hshape⟩
    refine ⟨g.toBlocks₁₁ 0 0,
      (g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 0 0,
      (g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 0 1,
      (g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 0 2,
      (g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 1 2,
      ha, ?_⟩
    calc
      g
          = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := by
              calc
                g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
                      exact (Matrix.fromBlocks_toBlocks g).symm
                _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := by
                      simpa [hB]
      _ = Wedge2Formalization.N5PureSingularLong.pointwiseShape
            (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₁ 0 1)
            (g.toBlocks₂₁ 0 2) (g.toBlocks₂₁ 1 2) := hshape
      _ =
          U (k := k) ((g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 0 0)
            ((g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 0 1)
            ((g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 0 2)
            ((g.toBlocks₁₁ 0 0)⁻¹ * g.toBlocks₂₁ 1 2) *
          Levi (k := k) (g.toBlocks₁₁ 0 0) := by
            rw [shape_eq_product (ha := ha)]
  · rintro ⟨a, u, v, w, z, ha, rfl⟩
    have hU := U_pointwise (k := k) u v w z
    have hL := Levi_pointwise (k := k) a ha
    have hL₁ :
        Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k)) (Levi (k := k) a) =
          rep₁ (k := k) := by
      simpa [rep₁] using hL.1
    have hL₂ :
        Wedge2Formalization.N5PureSingularLong.ActBivector (rep₂ (k := k)) (Levi (k := k) a) =
          rep₂ (k := k) := by
      simpa [rep₂] using hL.2
    have hU₁ :
        Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k)) (U (k := k) u v w z) =
          rep₁ (k := k) := by
      simpa [rep₁] using hU.1
    have hU₂ :
        Wedge2Formalization.N5PureSingularLong.ActBivector (rep₂ (k := k)) (U (k := k) u v w z) =
          rep₂ (k := k) := by
      simpa [rep₂] using hU.2
    constructor
    ·
      calc
        Wedge2Formalization.N5PureSingularLong.ActBivector
            (rep₁ (k := k)) (U (k := k) u v w z * Levi (k := k) a)
            =
          Wedge2Formalization.N5PureSingularLong.ActBivector
            (Wedge2Formalization.N5PureSingularLong.ActBivector
              (rep₁ (k := k)) (Levi (k := k) a))
            (U (k := k) u v w z) := by
              rw [Wedge2Formalization.N5PureSingularLong.actBivector_mul]
        _ = Wedge2Formalization.N5PureSingularLong.ActBivector
              (rep₁ (k := k)) (U (k := k) u v w z) := by rw [hL₁]
        _ = rep₁ (k := k) := hU₁
    ·
      calc
        Wedge2Formalization.N5PureSingularLong.ActBivector
            (rep₂ (k := k)) (U (k := k) u v w z * Levi (k := k) a)
            =
          Wedge2Formalization.N5PureSingularLong.ActBivector
            (Wedge2Formalization.N5PureSingularLong.ActBivector
              (rep₂ (k := k)) (Levi (k := k) a))
            (U (k := k) u v w z) := by
              rw [Wedge2Formalization.N5PureSingularLong.actBivector_mul]
        _ = Wedge2Formalization.N5PureSingularLong.ActBivector
              (rep₂ (k := k)) (U (k := k) u v w z) := by rw [hL₂]
        _ = rep₂ (k := k) := hU₂

theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k) :
    g ∈ K (k := k) ↔
      ∃ a u v w z : k, a ≠ 0 ∧ g = U (k := k) u v w z * Levi (k := k) a := by
  rfl

/-- Table-facing kernel statement for Appendix A, row 1:
`K_L = U_4 \rtimes G_m(k)`. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N5PureSingularLong.V
      Wedge2Formalization.N5PureSingularLong.V k) :
    g ∈ K (k := k) ↔
      ∃ a u v w z : k, a ≠ 0 ∧ g = U (k := k) u v w z * Levi (k := k) a :=
  mem_K_iff (k := k) (g := g)

theorem coeff_mem_Qproj
    (a b c d : k)
    (hΔ : Wedge2Formalization.N5PureSingularLong.Delta a b c d ≠ 0) :
    coeff (k := k) a b c d ∈ Qproj (k := k) := by
  have hdet :
      Matrix.det (coeff (k := k) a b c d) =
        Wedge2Formalization.N5PureSingularLong.Delta a b c d ^ 3 := by
    simp [coeff, Wedge2Formalization.N5PureSingularLong.Delta, Matrix.det_fin_two]
    ring
  change Matrix.det (coeff (k := k) a b c d) ≠ 0
  rw [hdet]
  exact pow_ne_zero 3 hΔ

theorem quotient_action
    (a b c d : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N5PureSingularLong.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b c d)
      (coeff (k := k) a b c d) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift, coeff,
      Wedge2Formalization.N5PureSingularLong.Delta, smul_smul, mul_assoc]
      using (Wedge2Formalization.N5Summary.pureSingularLong_GL2_scaled_lift_action
        (k := k) (a := a) (b := b) (c := c) (d := d)).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift, coeff,
      Wedge2Formalization.N5PureSingularLong.Delta, smul_smul, mul_assoc]
      using (Wedge2Formalization.N5Summary.pureSingularLong_GL2_scaled_lift_action
        (k := k) (a := a) (b := b) (c := c) (d := d)).2

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h03 := congrArg (fun M => M (Sum.inl 0) (Sum.inr 0)) h
  have ha : a = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₂, Wedge2Formalization.N5PureSingularLong.mulX,
      Wedge2Formalization.N5PureSingularLong.mulY, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h03
  have h13 := congrArg (fun M => M (Sum.inl 1) (Sum.inr 0)) h
  have hb : b = 0 := by
    simpa [ha, rep₁, rep₂, Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₂, Wedge2Formalization.N5PureSingularLong.mulX,
      Wedge2Formalization.N5PureSingularLong.mulY, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h13
  exact ⟨ha, hb⟩

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h13 := congrArg (fun M => M (Sum.inl 1) (Sum.inr 0)) hzero
  simpa [rep₂, Wedge2Formalization.N5PureSingularLong.rep₂,
    Wedge2Formalization.N5PureSingularLong.mulY, Matrix.fromBlocks] using h13

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have h03 := congrArg (fun M => M (Sum.inl 0) (Sum.inr 0)) hzero
  simpa [rep₁, Wedge2Formalization.N5PureSingularLong.rep₁,
    Wedge2Formalization.N5PureSingularLong.mulX, Matrix.fromBlocks] using h03

theorem quotient_image
    (g : Matrix Wedge2Formalization.N5PureSingularLong.V Wedge2Formalization.N5PureSingularLong.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N5PureSingularLong.PreservesSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N5PureSingularLong.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N5PureSingularLong.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  have hdetcoeff : Matrix.det (!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k) ≠ 0 := by
    intro hdet
    by_cases hleft : α = 0 ∧ γ = 0
    · by_cases hright : β = 0 ∧ δ = 0
      · have hz :
          Wedge2Formalization.N5PureSingularLong.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [hleft.1, hleft.2, hright.1] using h1'
        have hrep₁zero :=
          (Wedge2Formalization.N5PureSingularLong.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := rep₁ (k := k)) (g := g) hg).1 hz
        exact rep₁_ne_zero (k := k) hrep₁zero
      · have hzero :
          Wedge2Formalization.N5PureSingularLong.ActBivector
              (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g = 0 := by
          rw [Wedge2Formalization.N5PureSingularLong.actBivector_add,
            Wedge2Formalization.N5PureSingularLong.actBivector_smul,
            Wedge2Formalization.N5PureSingularLong.actBivector_smul,
            h1', h2']
          ext i j <;> simp [hleft.1, hleft.2, Matrix.add_apply, Matrix.smul_apply] <;> ring
        have hcomb_zero :=
          (Wedge2Formalization.N5PureSingularLong.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) (g := g) hg).1 hzero
        have hcomb_ne :
            δ • rep₁ (k := k) + (-β) • rep₂ (k := k) ≠ 0 := by
          intro hcomb
          rcases rep_pair_independent (k := k) hcomb with ⟨hδ, hβ'⟩
          exact hright ⟨by simpa using neg_eq_zero.mp hβ', hδ⟩
        exact hcomb_ne hcomb_zero
    · have hzero :
        Wedge2Formalization.N5PureSingularLong.ActBivector
            (γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) g = 0 := by
        have hdet' : α * δ - β * γ = 0 := by
          simpa [Matrix.det_fin_two] using hdet
        have hcoeff : γ * β = α * δ := by
          simpa [mul_comm] using (sub_eq_zero.mp hdet').symm
        rw [Wedge2Formalization.N5PureSingularLong.actBivector_add,
          Wedge2Formalization.N5PureSingularLong.actBivector_smul,
          Wedge2Formalization.N5PureSingularLong.actBivector_smul,
          h1', h2']
        ext i j
        simp [Matrix.add_apply, Matrix.smul_apply]
        have hexpr :
            γ * (α * rep₁ (k := k) i j) + γ * (β * rep₂ (k := k) i j) +
                (-(α * (γ * rep₁ (k := k) i j)) + -(α * (δ * rep₂ (k := k) i j))) =
              γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ := by
          ring
        have hcoeff_zero :
            γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ = 0 := by
          calc
            γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ
                = (γ * β - α * δ) * rep₂ (k := k) i j := by ring
            _ = 0 := by simp [hcoeff]
        calc
          γ * (α * rep₁ (k := k) i j) + γ * (β * rep₂ (k := k) i j) +
              (-(α * (γ * rep₁ (k := k) i j)) + -(α * (δ * rep₂ (k := k) i j)))
              = γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ := hexpr
          _ = 0 := hcoeff_zero
      have hcomb_zero :=
        (Wedge2Formalization.N5PureSingularLong.actBivector_eq_zero_iff_of_det_ne_zero
          (k := k) (Ω := γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) (g := g) hg).1 hzero
      have hcomb_ne :
          γ • rep₁ (k := k) + (-α) • rep₂ (k := k) ≠ 0 := by
        intro hcomb
        rcases rep_pair_independent (k := k) hcomb with ⟨hγ, hα'⟩
        exact hleft ⟨by simpa using neg_eq_zero.mp hα', hγ⟩
      exact hcomb_ne hcomb_zero
  refine ⟨!![α, β; γ, δ], ?_, ?_⟩
  · exact hdetcoeff
  · exact ⟨h1', h2'⟩

end Row1

/-! Appendix A, `n = 5`, row 2.
Representative `⟨e₁∧e₃, e₂∧e₃⟩`.
Divisor `0`.
Claimed stabilizer:
`K_L = U_8 \rtimes (GL_2(k) \times G_m(k))`, exact quotient family `Q_L = GL_2(k)`.
-/
namespace Row2

def rep₁ : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  Wedge2Formalization.N5PureSingular.radPureSingularRep₁ (k := k)

def rep₂ : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  Wedge2Formalization.N5PureSingular.radPureSingularRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 2:
`e₁∧e₃`. -/
def paperRep₁ : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  let R : Matrix Wedge2Formalization.N5PureSingular.I Wedge2Formalization.N5PureSingular.W k :=
    !![(0 : k), 1, 0, 0]
  Matrix.fromBlocks 0 R (-Rᵀ) 0

/-- Literal second basis vector from Appendix A, row 2:
`e₂∧e₃`. -/
def paperRep₂ : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.N4PureSingular.ω12 (k := k))

/-- Explicit basis change sending the paper row-2 representative to the internal
embedded pair `ω12, ω13`. -/
def paperChange : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  let TR : Matrix Wedge2Formalization.N5PureSingular.I Wedge2Formalization.N5PureSingular.W k :=
    !![(0 : k), 0, 1, 0]
  let BL : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k :=
    !![(0 : k); 1; 0; 0]
  let BR : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.W k :=
    !![(0 : k), -1, 0, 0;
       0, 0, 0, 0;
       1, 0, 0, 0;
       0, 0, 0, 1]
  Matrix.fromBlocks 0 TR BL BR

/-- Transport of the literal paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N5PureSingular.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [paperRep₁, paperChange, rep₁,
      Wedge2Formalization.N5PureSingular.radPureSingularRep₁,
      Wedge2Formalization.N4PureSingular.ω12,
      Wedge2Formalization.N5PureSingular.ActBivector,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_one, Fin.sum_univ_four]

/-- Transport of the literal paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N5PureSingular.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [paperRep₂, paperChange, rep₂,
      Wedge2Formalization.N5PureSingular.radPureSingularRep₂,
      Wedge2Formalization.N4PureSingular.ω12,
      Wedge2Formalization.N4PureSingular.ω13,
      Wedge2Formalization.N5PureSingular.ActBivector,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_one, Fin.sum_univ_four]

/-- Exact pointwise kernel family for the radical pure singular row. -/
def K : Set (Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k) :=
  { g |
      ∃ a : Matrix Wedge2Formalization.N5PureSingular.I Wedge2Formalization.N5PureSingular.I k,
        ∃ R : Matrix Wedge2Formalization.N5PureSingular.I Wedge2Formalization.N5PureSingular.W k,
          ∃ C : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k,
            ∃ D : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k,
              Wedge2Formalization.N5PaperSummary.row2_pointwiseShape (k := k) R D ∧
              g = Matrix.fromBlocks a R C D }

/-- The displayed radical unipotent family `U_6` on the radical pure singular row. -/
def U
    (x : k)
    (C : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k) :
    Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) (1 : k))
    (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
    C
    (1 : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k)

/-- The displayed Levi/core family on the radical pure singular row. -/
def Levi
    (u : k)
    (D : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) :
    Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
    0
    0
    D

/-- Exact coefficient-side quotient family `GL_2(k)` in the chosen basis. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

/-- Chosen `GL_2(k)` lift. -/
def lift
    (α β γ δ : k) :
    Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) (1 : k))
    0
    0
    (Wedge2Formalization.N4PureSingular.pureSingularSetwiseShape
      (k := k) 1 0 0 0 α γ β δ 0 0 1)

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N5PureSingular.FixesRadPureSingularPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    refine ⟨g.toBlocks₁₁, g.toBlocks₁₂, g.toBlocks₂₁, g.toBlocks₂₂, ?_, ?_⟩
    · have hg' :
        Wedge2Formalization.N5PureSingular.FixesRadPureSingularPairBivector
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) := by
        simpa [Matrix.fromBlocks_toBlocks] using hg
      exact
        (Wedge2Formalization.N5PaperSummary.row2_pointwise_iff
          (k := k)
          (a := g.toBlocks₁₁)
          (R := g.toBlocks₁₂)
          (C := g.toBlocks₂₁)
          (D := g.toBlocks₂₂)).1 hg'
    · exact (Matrix.fromBlocks_toBlocks g).symm
  · rintro ⟨a, R, C, D, hshape, rfl⟩
    exact
      (Wedge2Formalization.N5PaperSummary.row2_pointwise_iff
        (k := k) (a := a) (R := R) (C := C) (D := D)).2 hshape

theorem U_pointwise
    (x : k)
    (C : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k) :
    Wedge2Formalization.N5PureSingular.FixesRadPureSingularPairBivector (U (k := k) x C) := by
  have hI :
      (1 : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) =
        (!![(1 : k), 0, 0, 0;
            0, 1, 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1] :
          Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp
  rw [U, hI]
  simpa [U, Wedge2Formalization.N5PureSingular.scalarBlock,
    Wedge2Formalization.N4PureSingular.pureSingularShape] using
    Wedge2Formalization.N5PureSingular.radPureSingular_pointwise_family_with_upperRight
      (k := k) (a := (1 : k)) (b := 0) (c := 0) (d := 0) (e := 0) (f := 0) (t := (1 : k))
      (u := (1 : k)) (x := x) one_ne_zero C

theorem Levi_pointwise
    (u a b c d e f t : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N5PureSingular.FixesRadPureSingularPairBivector
      (Levi (k := k) u
        (Wedge2Formalization.N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
  simpa [Levi, Wedge2Formalization.N5PureSingular.scalarBlock] using
    Wedge2Formalization.N5PureSingular.radPureSingular_pointwise_family
      (k := k) (a := a) (b := b) (c := c) (d := d) (e := e) (f := f) (t := t) (u := u) ha 0

theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k) :
    g ∈ K (k := k) ↔
      ∃ u x : k,
        ∃ C : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k,
          ∃ a b c d e f t : k,
            a ≠ 0 ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
                (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
                C
                (Wedge2Formalization.N4PureSingular.pureSingularShape
                  (k := k) a b c d e f t) := by
  constructor
  · rintro ⟨aM, R, C, D, hshape, hg⟩
    rcases hshape with ⟨hR, hD00, hDshape⟩
    let u : k := aM 0 0
    let x : k := R 0 3
    let a : k := D 0 0
    let b : k := D 0 1
    let c : k := D 0 2
    let d : k := D 0 3
    let e : k := D 1 3
    let f : k := D 2 3
    let t : k := D 3 3
    refine ⟨u, x, C, a, b, c, d, e, f, t, hD00, ?_⟩
    have haM :
        aM = Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u := by
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [u, Wedge2Formalization.N5PureSingular.scalarBlock]
    calc
      g = Matrix.fromBlocks aM R C D := hg
      _ =
          Matrix.fromBlocks
            (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
            R
            C
            D := by rw [haM]
      _ =
          Matrix.fromBlocks
            (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
            (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
            C
            D := by rw [show R = Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x by
              simpa [x] using hR]
      _ =
          Matrix.fromBlocks
            (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
            (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
            C
            (Wedge2Formalization.N4PureSingular.pureSingularShape
              (k := k) a b c d e f t) := by
            rw [show D = Wedge2Formalization.N4PureSingular.pureSingularShape (k := k) a b c d e f t by
              simpa [a, b, c, d, e, f, t] using hDshape]
  · rintro ⟨u, x, C, a, b, c, d, e, f, t, ha, rfl⟩
    refine ⟨Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u,
      Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x, C,
      Wedge2Formalization.N4PureSingular.pureSingularShape (k := k) a b c d e f t,
      ?_, rfl⟩
    exact ⟨rfl, ha, rfl⟩

/-- The embedded pure singular core, i.e. the `U_5 \rtimes (G_m(k) \times G_m(k))`
factor from Appendix A, `n = 4`, row 3. -/
def PureSingularCore :
    Set
      (Matrix Wedge2Formalization.N4PureSingular.I
        Wedge2Formalization.N4PureSingular.I k) :=
  Wedge2Formalization.Paper.N4.Row3.K (k := k)

/-- Secondary kernel-membership statement through the embedded `n = 4`, row 3 core. -/
theorem mem_K_core_iff
    (g : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k) :
    g ∈ K (k := k) ↔
      ∃ u x : k,
        ∃ C : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k,
          ∃ D : Matrix Wedge2Formalization.N4PureSingular.I
              Wedge2Formalization.N4PureSingular.I k,
            D ∈ PureSingularCore (k := k) ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
                (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
                C
                D := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨u, x, C, a, b, c, d, e, f, t, ha, rfl⟩
    refine ⟨u, x, C,
      Wedge2Formalization.N4PureSingular.pureSingularShape (k := k) a b c d e f t, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N4.Row3.mem_K_iff
        (k := k)
        (g := Wedge2Formalization.N4PureSingular.pureSingularShape (k := k) a b c d e f t)).2
        ⟨a, b, c, d, e, f, t, ha, rfl⟩
  · rintro ⟨u, x, C, D, hD, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row3.mem_K_iff
        (k := k)
        (g := D)).1 hD with
      ⟨a, b, c, d, e, f, t, ha, rfl⟩
    exact
      (mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks
          (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
          (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
          C
          (Wedge2Formalization.N4PureSingular.pureSingularShape (k := k) a b c d e f t))).2
        ⟨u, x, C, a, b, c, d, e, f, t, ha, rfl⟩

/-- Public kernel-membership statement in literal paper coordinates. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k) :
    g ∈ K (k := k) ↔
      ∃ u x : k,
        ∃ C : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k,
          ∃ D : Matrix Wedge2Formalization.N4PureSingular.I
              Wedge2Formalization.N4PureSingular.I k,
            D ∈ PureSingularCore (k := k) ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
                (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
                C
                D := by
  exact mem_K_core_iff (k := k) (g := g)

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_8 \rtimes (GL_2(k) \times G_m(k))`, written in the explicit local
coordinates of this row rather than through a nested `n = 4` core theorem. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k) :
    g ∈ K (k := k) ↔
      ∃ u x : k,
        ∃ C : Matrix Wedge2Formalization.N5PureSingular.W Wedge2Formalization.N5PureSingular.I k,
          ∃ a b c d e f t : k,
            a ≠ 0 ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5PureSingular.scalarBlock (k := k) u)
                (Wedge2Formalization.N5PureSingular.upperRightLast (k := k) x)
                C
                (Wedge2Formalization.N4PureSingular.pureSingularShape
                  (k := k) a b c d e f t) :=
  mem_K_shape_iff (k := k) (g := g)

theorem quotient_action
    (α β γ δ : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N5PureSingular.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) α β γ δ)
      (!![α, β; γ, δ]) := by
  simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
    Wedge2Formalization.N5Summary.radPureSingular_GL2_lift_action
      (k := k) (α := α) (β := β) (γ := γ) (δ := δ) (u := (1 : k)) (C := 0)

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h01 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h
  have ha : a = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N5PureSingular.radPureSingularRep₁,
      Wedge2Formalization.N5PureSingular.radPureSingularRep₂,
      Wedge2Formalization.N4PureSingular.ω12, Wedge2Formalization.N4PureSingular.ω13,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h01
  have h02 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 2)) h
  have hb : b = 0 := by
    simpa [ha, rep₁, rep₂, Wedge2Formalization.N5PureSingular.radPureSingularRep₁,
      Wedge2Formalization.N5PureSingular.radPureSingularRep₂,
      Wedge2Formalization.N4PureSingular.ω12, Wedge2Formalization.N4PureSingular.ω13,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h02
  exact ⟨ha, hb⟩

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h02 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 2)) hzero
  simpa [rep₂, Wedge2Formalization.N5PureSingular.radPureSingularRep₂,
    Wedge2Formalization.N4PureSingular.ω13, Matrix.fromBlocks] using h02

theorem quotient_image
    (g : Matrix Wedge2Formalization.N5PureSingular.V Wedge2Formalization.N5PureSingular.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N5PureSingular.PreservesRadPureSingularSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N5PureSingular.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N5PureSingular.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N5PureSingular.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  have hdetcoeff : Matrix.det (!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k) ≠ 0 := by
    intro hdet
    by_cases hleft : α = 0 ∧ γ = 0
    · by_cases hright : β = 0 ∧ δ = 0
      · have hz :
          Wedge2Formalization.N5PureSingular.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [hleft.1, hleft.2, hright.1] using h1
        have hrep₁zero :=
          (Wedge2Formalization.N5PureSingular.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := rep₁ (k := k)) (g := g) hg).1 hz
        have hrep₁nz : rep₁ (k := k) ≠ 0 := by
          intro hzero
          have h01 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) hzero
          simpa [rep₁, Wedge2Formalization.N5PureSingular.radPureSingularRep₁,
            Wedge2Formalization.N4PureSingular.ω12, Matrix.fromBlocks] using h01
        exact hrep₁nz hrep₁zero
      · have hzero :
          Wedge2Formalization.N5PureSingular.ActBivector
              (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g = 0 := by
          rw [Wedge2Formalization.N5PureSingular.actBivector_add,
            Wedge2Formalization.N5PureSingular.actBivector_smul,
            Wedge2Formalization.N5PureSingular.actBivector_smul,
            h1', h2']
          ext i j <;> simp [hleft.1, hleft.2] <;> ring
        have hcomb_zero :=
          (Wedge2Formalization.N5PureSingular.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) (g := g) hg).1 hzero
        have hcomb_ne :
            δ • rep₁ (k := k) + (-β) • rep₂ (k := k) ≠ 0 := by
          intro hcomb
          rcases rep_pair_independent (k := k) hcomb with ⟨hδ, hβ'⟩
          exact hright ⟨by simpa using neg_eq_zero.mp hβ', hδ⟩
        exact hcomb_ne hcomb_zero
    · have hzero :
        Wedge2Formalization.N5PureSingular.ActBivector
            (γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) g = 0 := by
        have hdet' : α * δ - β * γ = 0 := by
          simpa [Matrix.det_fin_two] using hdet
        have hcoeff : γ * β = α * δ := by
          simpa [mul_comm] using (sub_eq_zero.mp hdet').symm
        rw [Wedge2Formalization.N5PureSingular.actBivector_add,
          Wedge2Formalization.N5PureSingular.actBivector_smul,
          Wedge2Formalization.N5PureSingular.actBivector_smul,
          h1', h2']
        ext i j
        simp [Matrix.add_apply, Matrix.smul_apply]
        have hexpr :
            γ * (α * rep₁ (k := k) i j) + γ * (β * rep₂ (k := k) i j) +
                (-(α * (γ * rep₁ (k := k) i j)) + -(α * (δ * rep₂ (k := k) i j))) =
              γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ := by
          ring
        have hcoeff_zero :
            γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ = 0 := by
          calc
            γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ
                = (γ * β - α * δ) * rep₂ (k := k) i j := by ring
            _ = 0 := by simp [hcoeff]
        calc
          γ * (α * rep₁ (k := k) i j) + γ * (β * rep₂ (k := k) i j) +
              (-(α * (γ * rep₁ (k := k) i j)) + -(α * (δ * rep₂ (k := k) i j)))
              =
            γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ := hexpr
          _ = 0 := hcoeff_zero
      have hcomb_zero :=
        (Wedge2Formalization.N5PureSingular.actBivector_eq_zero_iff_of_det_ne_zero
          (k := k) (Ω := γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) (g := g) hg).1 hzero
      have hcomb_ne :
          γ • rep₁ (k := k) + (-α) • rep₂ (k := k) ≠ 0 := by
        intro hcomb
        rcases rep_pair_independent (k := k) hcomb with ⟨hγ, hα'⟩
        exact hleft ⟨by simpa using neg_eq_zero.mp hα', hγ⟩
      exact hcomb_ne hcomb_zero
  refine ⟨!![α, β; γ, δ], ?_, ?_⟩
  · exact hdetcoeff
  · simpa [ActsOnOrderedPair] using And.intro h1 h2

end Row2

/-! Appendix A, `n = 5`, row 3.
Representative `⟨e₁∧e₃ + e₄∧e₅, e₂∧e₃⟩`.
Divisor `[a]`.
Claimed stabilizer:
`K_L = U_4 \rtimes (G_m(k) \times SL_2(k))`, exact quotient family `Q_L = B`.
-/
namespace Row3

def rep₁ : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  Wedge2Formalization.N5SimplePoint.rep₁ (k := k)

def rep₂ : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  Wedge2Formalization.N5SimplePoint.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 3. -/
def paperRep₁ : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 3. -/
def paperRep₂ : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  rep₂ (k := k)

/-- The row-3 paper representative already agrees with the internal working pair. -/
def paperChange : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N5SimplePoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N5SimplePoint.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N5SimplePoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N5SimplePoint.ActBivector]

/-- Exact pointwise kernel family for the direct-sum `[a]` row, packaged in the
paper's explicit matrix-shape language. -/
def K :
    Set (Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k) :=
  { g |
      ∃ a x y u v : k,
        ∃ E : Matrix Wedge2Formalization.N5SimplePoint.I Wedge2Formalization.N5SimplePoint.I k,
          a ≠ 0 ∧
          E.det = 1 ∧
          g =
            Matrix.fromBlocks
              (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
              (!![(0 : k), 0; 0, 0; u, v])
              (!![a * (u * E 0 1 - v * E 0 0), 0, 0;
                 a * (u * E 1 1 - v * E 1 0), 0, 0])
              E }

/-- The displayed unipotent family `U_4`. -/
def U
    (x y p q : k) :
    Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  Wedge2Formalization.N5SimplePoint.pointwiseUnipotent (k := k) x y p q

/-- The displayed Levi family `G_m x SL_2(k)`. -/
def Levi
    (a : k)
    (E : Matrix Wedge2Formalization.N5SimplePoint.I Wedge2Formalization.N5SimplePoint.I k) :
    Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N3PureSingular.pointwiseScale (k := k) a)
    0
    0
    E

/-- Exact coefficient-side projective Borel family in the ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

/-- A convenient explicit section inside the projective Borel family. -/
def Q : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | ∃ a b : k, a ≠ 0 ∧ M = Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b }

/-- Chosen Borel lift. -/
def lift
    (a b : k) :
    Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k :=
  let G : Matrix Wedge2Formalization.N5SimplePoint.W Wedge2Formalization.N5SimplePoint.W k :=
    Wedge2Formalization.N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
  let E : Matrix Wedge2Formalization.N5SimplePoint.I Wedge2Formalization.N5SimplePoint.I k :=
    !![a, 0; 0, 1]
  Matrix.fromBlocks G 0 0 E

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N5SimplePoint.FixesPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    have hg' :
        Wedge2Formalization.N5SimplePoint.FixesPairBivector
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) := by
      simpa [Matrix.fromBlocks_toBlocks] using hg
    rcases
      (Wedge2Formalization.N5SimplePoint.fixesPair_fromBlocks_iff_shape
        (k := k) (A := g.toBlocks₁₁) (B := g.toBlocks₁₂) (C := g.toBlocks₂₁) (D := g.toBlocks₂₂)).1 hg'
      with ⟨ha, hdetE, hshape⟩
    refine ⟨g.toBlocks₁₁ 0 0, g.toBlocks₁₁ 2 0, g.toBlocks₁₁ 2 1,
      g.toBlocks₁₂ 2 0, g.toBlocks₁₂ 2 1, g.toBlocks₂₂, ha, hdetE, ?_⟩
    rw [show g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ by
      exact (Matrix.fromBlocks_toBlocks g).symm]
    exact hshape
  · rintro ⟨a, x, y, u, v, E, ha, hdetE, rfl⟩
    exact
      (Wedge2Formalization.N5SimplePoint.fixesPair_fromBlocks_iff_shape
        (k := k)
        (A := Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
        (B := !![(0 : k), 0; 0, 0; u, v])
        (C := !![a * (u * E 0 1 - v * E 0 0), 0, 0;
                 a * (u * E 1 1 - v * E 1 0), 0, 0])
        (D := E)).2 ⟨ha, hdetE, rfl⟩

theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ a x y u v : k,
        ∃ E : Matrix Wedge2Formalization.N5SimplePoint.I Wedge2Formalization.N5SimplePoint.I k,
          a ≠ 0 ∧
          E.det = 1 ∧
          g =
            Matrix.fromBlocks
              (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
              (!![(0 : k), 0; 0, 0; u, v])
              (!![a * (u * E 0 1 - v * E 0 0), 0, 0;
                 a * (u * E 1 1 - v * E 1 0), 0, 0])
              E := by
  rfl

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_4 \rtimes (G_m(k) \times SL_2(k))`, written in the explicit public
coordinates of the direct-sum `[a]` family. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ a x y u v : k,
        ∃ E : Matrix Wedge2Formalization.N5SimplePoint.I
            Wedge2Formalization.N5SimplePoint.I k,
          a ≠ 0 ∧
          E.det = 1 ∧
          g =
            Matrix.fromBlocks
              (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
              (!![(0 : k), 0; 0, 0; u, v])
              (!![a * (u * E 0 1 - v * E 0 0), 0, 0;
                 a * (u * E 1 1 - v * E 1 0), 0, 0])
              E :=
  mem_K_iff (k := k) (g := g)

/-- The coefficient matrix attached to the standard lift belongs to the explicit Borel
section. -/
theorem coeff_mem_Q
    (a b : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b ∈ Q (k := k) :=
  ⟨a, b, ha, rfl⟩

/-- The coefficient matrix attached to the standard lift lies in the exact projective Borel
family. -/
theorem coeff_mem_Qproj
    (a b : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b ∈ Qproj (k := k) := by
  constructor
  · simp [Wedge2Formalization.N4PaperSummary.row1_borelCoeff]
  · simp [Wedge2Formalization.N4PaperSummary.row1_borelCoeff, Matrix.det_fin_two, ha]

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N5SimplePoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N5Summary.simplePoint_borel_lift_action
        (k := k) (a := a) (b := b) ha).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift, Wedge2Formalization.N4PaperSummary.row1_borelCoeff] using
      (Wedge2Formalization.N5Summary.simplePoint_borel_lift_action
        (k := k) (a := a) (b := b) ha).2

private theorem act_eq_zero_iff_of_det_ne_zero
    (Ω : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k)
    (g : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k)
    (hg : Matrix.det g ≠ 0) :
    Wedge2Formalization.N5SimplePoint.ActBivector Ω g = 0 ↔ Ω = 0 := by
  letI : Invertible (Matrix.det g) := invertibleOfNonzero hg
  letI : Invertible g := Matrix.invertibleOfDetInvertible g
  have hunit : IsUnit (Matrix.det g) := isUnit_of_invertible (Matrix.det g)
  constructor
  · intro hzero
    have h' := congrArg (fun M => g⁻¹ * M * (g⁻¹)ᵀ) hzero
    simp [Wedge2Formalization.N5SimplePoint.ActBivector, Matrix.mul_assoc,
      Matrix.nonsing_inv_mul _ hunit, Matrix.mul_nonsing_inv _ hunit,
      Matrix.transpose_nonsing_inv] at h'
    exact h'
  · intro hΩ
    simp [hΩ, Wedge2Formalization.N5SimplePoint.ActBivector]

private def rows₂ : Fin 2 → Wedge2Formalization.N5SimplePoint.V
  | 0 => Sum.inl 1
  | 1 => Sum.inl 2

private def pick₂ : Matrix Wedge2Formalization.N5SimplePoint.V (Fin 2) k :=
  fun i j =>
    match i, j with
    | Sum.inl 2, 0 => 1
    | Sum.inl 1, 1 => -1
    | _, _ => 0

private theorem rows₂_mul_pick₂ :
    (rep₂ (k := k)).submatrix rows₂ (Equiv.refl _) * pick₂ (k := k) =
      (1 : Matrix (Fin 2) (Fin 2) k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rows₂, pick₂, rep₂, Wedge2Formalization.N5SimplePoint.rep₂,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.submatrix, Matrix.mul_apply,
      Matrix.fromBlocks, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]

private def rep₂_left : Matrix Wedge2Formalization.N5SimplePoint.V (Fin 2) k :=
  fun i j =>
    match i, j with
    | Sum.inl 1, 0 => 1
    | Sum.inl 2, 1 => 1
    | _, _ => 0

private def rep₂_right : Matrix (Fin 2) Wedge2Formalization.N5SimplePoint.V k :=
  fun i j =>
    match i, j with
    | 0, Sum.inl 2 => 1
    | 1, Sum.inl 1 => -1
    | _, _ => 0

private theorem rep₂_factor :
    rep₂ (k := k) = rep₂_left (k := k) * rep₂_right (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Fin.sum_univ_two,
        Matrix.fromBlocks]
  ·
    fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Fin.sum_univ_two,
        Matrix.fromBlocks]
  ·
    fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Fin.sum_univ_two,
        Matrix.fromBlocks]
  ·
    fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N5SimplePoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Fin.sum_univ_two,
        Matrix.fromBlocks]

private theorem rep₂_rank_eq_two : (rep₂ (k := k)).rank = 2 := by
  have hle : (rep₂ (k := k)).rank ≤ 2 := by
    rw [rep₂_factor (k := k)]
    exact (Matrix.rank_mul_le_left (rep₂_left (k := k)) (rep₂_right (k := k))).trans <| by
      simpa using (Matrix.rank_le_card_width (rep₂_left (k := k)))
  have hge : 2 ≤ (rep₂ (k := k)).rank := by
    calc
      2 = (1 : Matrix (Fin 2) (Fin 2) k).rank := by
        simpa using (Matrix.rank_one (R := k) (n := Fin 2))
      _ = (((rep₂ (k := k)).submatrix rows₂ (Equiv.refl _)) * pick₂ (k := k)).rank := by
        rw [rows₂_mul_pick₂ (k := k)]
      _ ≤ ((rep₂ (k := k)).submatrix rows₂ (Equiv.refl _)).rank := by
        exact Matrix.rank_mul_le_left _ _
      _ ≤ (rep₂ (k := k)).rank := by
        exact Matrix.rank_submatrix_le rows₂ (Equiv.refl _) (rep₂ (k := k))
  exact le_antisymm hle hge

private def rows₄ : Fin 4 → Wedge2Formalization.N5SimplePoint.V
  | 0 => Sum.inl 0
  | 1 => Sum.inl 2
  | 2 => Sum.inr 0
  | 3 => Sum.inr 1

private def pick₄ (γ : k) : Matrix Wedge2Formalization.N5SimplePoint.V (Fin 4) k :=
  fun i j =>
    match i, j with
    | Sum.inl 2, 0 => γ⁻¹
    | Sum.inl 0, 1 => -γ⁻¹
    | Sum.inr 1, 2 => γ⁻¹
    | Sum.inr 0, 3 => -γ⁻¹
    | _, _ => 0

private theorem rows₄_mul_pick₄
    (γ δ : k)
    (hγ : γ ≠ 0) :
    ((γ • rep₁ (k := k) + δ • rep₂ (k := k)).submatrix rows₄ (Equiv.refl _)) * pick₄ (k := k) γ =
      (1 : Matrix (Fin 4) (Fin 4) k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [rows₄, pick₄, rep₁, rep₂, Wedge2Formalization.N5SimplePoint.rep₁,
      Wedge2Formalization.N5SimplePoint.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Wedge2Formalization.N4.J, Matrix.submatrix,
      Matrix.mul_apply, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply, hγ,
      Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two]

private theorem combo_rank_ge_four
    (γ δ : k)
    (hγ : γ ≠ 0) :
    4 ≤ (γ • rep₁ (k := k) + δ • rep₂ (k := k)).rank := by
  calc
    4 = (1 : Matrix (Fin 4) (Fin 4) k).rank := by
      simpa using (Matrix.rank_one (R := k) (n := Fin 4))
    _ = ((((γ • rep₁ (k := k) + δ • rep₂ (k := k)).submatrix rows₄ (Equiv.refl _)) *
          pick₄ (k := k) γ)).rank := by
      rw [rows₄_mul_pick₄ (k := k) γ δ hγ]
    _ ≤ ((γ • rep₁ (k := k) + δ • rep₂ (k := k)).submatrix rows₄ (Equiv.refl _)).rank := by
      exact Matrix.rank_mul_le_left _ _
    _ ≤ (γ • rep₁ (k := k) + δ • rep₂ (k := k)).rank := by
      exact Matrix.rank_submatrix_le rows₄ (Equiv.refl _) _

private theorem act_rep₂_rank_eq_two
    (g : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k)
    (hg : Matrix.det g ≠ 0) :
    (Wedge2Formalization.N5SimplePoint.ActBivector (rep₂ (k := k)) g).rank = 2 := by
  have hg_unit : IsUnit (Matrix.det g) := isUnit_iff_ne_zero.mpr hg
  have hgt_unit : IsUnit (Matrix.det gᵀ) := by
    simpa [Matrix.det_transpose] using hg_unit
  calc
    (Wedge2Formalization.N5SimplePoint.ActBivector (rep₂ (k := k)) g).rank
        = (g * rep₂ (k := k) * gᵀ).rank := by
          rfl
    _ = (g * rep₂ (k := k)).rank := by
          simpa [Wedge2Formalization.N5SimplePoint.ActBivector, Matrix.mul_assoc] using
            (Matrix.rank_mul_eq_left_of_isUnit_det (A := gᵀ) (B := g * rep₂ (k := k)) hgt_unit)
    _ = (rep₂ (k := k)).rank := by
          simpa using
            (Matrix.rank_mul_eq_right_of_isUnit_det (A := g) (B := rep₂ (k := k)) hg_unit)
    _ = 2 := rep₂_rank_eq_two (k := k)

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h01 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h
  have ha : a = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N5SimplePoint.rep₁,
      Wedge2Formalization.N5SimplePoint.rep₂, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h01
  have h12 := congrArg (fun M => M (Sum.inl 1) (Sum.inl 2)) h
  have hb : b = 0 := by
    simpa [ha, rep₁, rep₂, Wedge2Formalization.N5SimplePoint.rep₁,
      Wedge2Formalization.N5SimplePoint.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h12
  exact ⟨ha, hb⟩

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h12 := congrArg (fun M => M (Sum.inl 1) (Sum.inl 2)) hzero
  simpa [rep₂, Wedge2Formalization.N5SimplePoint.rep₂, Wedge2Formalization.N3PureSingular.ω23,
    Matrix.fromBlocks] using h12

/-- The rank obstruction on the quotient pencil occurs exactly at the unique support point
`[rep₂]`. This is the concrete invariant used to force the quotient image into the Borel. -/
theorem rank_ge_four_of_leftCoeff_ne_zero
    (a b : k)
    (ha : a ≠ 0) :
    4 ≤ (a • rep₁ (k := k) + b • rep₂ (k := k)).rank :=
  combo_rank_ge_four (k := k) a b ha

theorem quotient_image
    (g : Matrix Wedge2Formalization.N5SimplePoint.V Wedge2Formalization.N5SimplePoint.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N5SimplePoint.PreservesSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N5SimplePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N5SimplePoint.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N5SimplePoint.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  have hγ : γ = 0 := by
    by_contra hγnz
    have hbound :
        4 ≤ (Wedge2Formalization.N5SimplePoint.ActBivector (rep₂ (k := k)) g).rank := by
      simpa [h2'] using combo_rank_ge_four (k := k) γ δ hγnz
    rw [act_rep₂_rank_eq_two (k := k) (g := g) hg] at hbound
    omega
  have hδ : δ ≠ 0 := by
    intro hδ0
    have hzero :
        Wedge2Formalization.N5SimplePoint.ActBivector (rep₂ (k := k)) g = 0 := by
      simpa [hγ, hδ0] using h2'
    have horig :=
      (act_eq_zero_iff_of_det_ne_zero (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
    exact rep₂_ne_zero (k := k) horig
  have hα : α ≠ 0 := by
    intro hα0
    have hzero :
        Wedge2Formalization.N5SimplePoint.ActBivector
          (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g = 0 := by
      calc
        Wedge2Formalization.N5SimplePoint.ActBivector
            (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g
            =
          δ • Wedge2Formalization.N5SimplePoint.ActBivector (rep₁ (k := k)) g +
            (-β) • Wedge2Formalization.N5SimplePoint.ActBivector (rep₂ (k := k)) g := by
              simp [Wedge2Formalization.N5SimplePoint.ActBivector, Matrix.mul_add, Matrix.add_mul,
                Matrix.mul_smul, smul_mul_assoc]
        _ = 0 := by
              rw [h1', h2', hγ, hα0]
              ext i j
              simp [Matrix.add_apply, Matrix.smul_apply]
              ring
    have horig :=
      (act_eq_zero_iff_of_det_ne_zero
        (k := k) (Ω := δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) (g := g) hg).1 hzero
    have hcomb_ne :
        δ • rep₁ (k := k) + (-β) • rep₂ (k := k) ≠ 0 := by
      intro hcomb
      rcases rep_pair_independent (k := k) hcomb with ⟨hδ', hβ'⟩
      exact hδ hδ'
    exact hcomb_ne horig
  refine ⟨!![α, β; 0, δ], ?_, ?_⟩
  · constructor
    · simp
    · simp [Matrix.det_fin_two, hα, hδ]
  · constructor
    · simpa [ActsOnOrderedPair] using h1
    · simpa [ActsOnOrderedPair, hγ] using h2

end Row3

private theorem act_embedded_toBlocks₂₂
    (Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    Matrix.toBlocks₂₂
        (Wedge2Formalization.N5.ActBivector (Matrix.fromBlocks 0 0 0 Ω) g) =
      Wedge2Formalization.N4.ActBivector Ω g.toBlocks₂₂ := by
  rw [← Matrix.fromBlocks_toBlocks g]
  simp [Wedge2Formalization.N5.ActBivector, Wedge2Formalization.N4.ActBivector,
    Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

/-! Appendix A, `n = 5`, row 4.
Representative `⟨e₂∧e₄ + e₃∧e₅, e₂∧e₅⟩`.
Divisor `2[a]`.
Claimed stabilizer:
`K_L = U_7 \rtimes (G_m(k) \times SL_2(k))`, exact quotient family `Q_L = B`.
-/
namespace Row4

def rep₁ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)

def rep₂ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 4:
`e₂∧e₄ + e₃∧e₅`. -/
def paperRep₁ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))

/-- Literal second basis vector from Appendix A, row 4:
`e₂∧e₅`. -/
def paperRep₂ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))

/-- Explicit basis change sending the literal paper representative to the internal
working representative pair. -/
def paperChange : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.I k)
    0
    0
    (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))

/-- Transport of the literal paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N5.ActBivector (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  calc
    Wedge2Formalization.N5.ActBivector (paperRep₁ (k := k)) (paperChange (k := k)) =
        Matrix.fromBlocks 0 0 0
          (Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))
            (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))) := by
          simpa [paperRep₁, paperChange] using
            (Wedge2Formalization.N5.act_embedded_fromBlocks_zeroUpperRight
              (k := k)
              (Ω := Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))
              (a := (1 : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.I k))
              (C := (0 : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k))
              (D := Wedge2Formalization.Paper.N4.Row1.paperChange (k := k)))
    _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.N4.onePointRep₁ (k := k)) := by
          simpa [Wedge2Formalization.Paper.N4.Row1.rep₁] using
            congrArg
              (fun Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k =>
                Matrix.fromBlocks
                  (0 : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.I k)
                  (0 : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N4.V k)
                  (0 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N5.I k)
                  Ω)
              (Wedge2Formalization.Paper.N4.Row1.paperRep₁_transport (k := k))
    _ = rep₁ (k := k) := by
          simp [rep₁, Wedge2Formalization.N5OnePoint.radOnePointRep₁]

/-- Transport of the literal paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N5.ActBivector (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  calc
    Wedge2Formalization.N5.ActBivector (paperRep₂ (k := k)) (paperChange (k := k)) =
        Matrix.fromBlocks 0 0 0
          (Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))
            (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))) := by
          simpa [paperRep₂, paperChange] using
            (Wedge2Formalization.N5.act_embedded_fromBlocks_zeroUpperRight
              (k := k)
              (Ω := Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))
              (a := (1 : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.I k))
              (C := (0 : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k))
              (D := Wedge2Formalization.Paper.N4.Row1.paperChange (k := k)))
    _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
          simpa [Wedge2Formalization.Paper.N4.Row1.rep₂] using
            congrArg
              (fun Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k =>
                Matrix.fromBlocks
                  (0 : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.I k)
                  (0 : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N4.V k)
                  (0 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N5.I k)
                  Ω)
              (Wedge2Formalization.Paper.N4.Row1.paperRep₂_transport (k := k))
    _ = rep₂ (k := k) := by
          simp [rep₂, Wedge2Formalization.N5OnePoint.radOnePointRep₂]

/-- Exact pointwise kernel family for the radical one-point row. -/
def K : Set (Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :=
  { g |
      ∃ a : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.I k,
        ∃ R : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.W k,
          ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
            ∃ A B₁ C₁ D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
              Wedge2Formalization.N5PaperSummary.row4_pointwiseShape (k := k) R A B₁ C₁ D ∧
              g = Matrix.fromBlocks a R C (Matrix.fromBlocks A B₁ C₁ D) }

/-- The displayed unipotent radical on the radical one-point row. -/
def U
    (x y z : k)
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) (1 : k))
    0
    C
    (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z)

/-- The displayed Levi/core family on the radical one-point row. -/
def Levi
    (u : k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) u)
    0
    0
    (Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)

/-- The scalar factor on the radical line. -/
def Scale
    (u : k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) u)
    0
    0
    (1 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)

/-- The embedded `SL₂(k)` Levi factor from the repeated-support `n = 4` core. -/
def CoreLevi
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) (1 : k))
    0
    0
    (Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)

/-- Exact coefficient-side projective Borel family in the ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

/-- A convenient explicit section inside the projective Borel family. -/
def Q : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | ∃ a b : k, a ≠ 0 ∧ M = Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b }

/-- Chosen Borel lift. -/
def lift
    (a b : k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  let B : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k :=
    !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) (1 : k))
    0
    0
    (Wedge2Formalization.N4.onePointScale (k := k) a *
      Wedge2Formalization.N4.onePointUpperShear (k := k) B)

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N5OnePoint.FixesRadOnePointPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    refine
      ⟨g.toBlocks₁₁, g.toBlocks₁₂, g.toBlocks₂₁,
        g.toBlocks₂₂.toBlocks₁₁, g.toBlocks₂₂.toBlocks₁₂,
        g.toBlocks₂₂.toBlocks₂₁, g.toBlocks₂₂.toBlocks₂₂, ?_, ?_⟩
    · have hg' :
        Wedge2Formalization.N5OnePoint.FixesRadOnePointPairBivector
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁
            (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
              g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂)) := by
        have hEq22 :
            Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
              g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂ = g.toBlocks₂₂ := by
          simpa using Matrix.fromBlocks_toBlocks g.toBlocks₂₂
        have hEq :
            Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁
                (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
                g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂) = g := by
          rw [hEq22]
          exact Matrix.fromBlocks_toBlocks g
        simpa [hEq] using hg
      exact
        (Wedge2Formalization.N5PaperSummary.row4_pointwise_iff
          (k := k)
          (a := g.toBlocks₁₁)
          (R := g.toBlocks₁₂)
          (C := g.toBlocks₂₁)
          (A := g.toBlocks₂₂.toBlocks₁₁)
          (B₁ := g.toBlocks₂₂.toBlocks₁₂)
          (C₁ := g.toBlocks₂₂.toBlocks₂₁)
          (D := g.toBlocks₂₂.toBlocks₂₂)).1 hg'
    · have hEq22 :
          Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
            g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂ = g.toBlocks₂₂ := by
        simpa using Matrix.fromBlocks_toBlocks g.toBlocks₂₂
      rw [hEq22]
      exact (Matrix.fromBlocks_toBlocks g).symm
  · rintro ⟨a, R, C, A, B₁, C₁, D, hshape, rfl⟩
    exact
      (Wedge2Formalization.N5PaperSummary.row4_pointwise_iff
        (k := k) (a := a) (R := R) (C := C) (A := A) (B₁ := B₁) (C₁ := C₁) (D := D)).2 hshape

theorem U_pointwise
    (x y z : k)
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k) :
    Wedge2Formalization.N5OnePoint.FixesRadOnePointPairBivector (U (k := k) x y z C) := by
  have htrace' :
      (!![x, y; z, -x] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) *
          Wedge2Formalization.N4.J +
        Wedge2Formalization.N4.J *
          (!![x, y; z, -x] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)ᵀ =
        (0 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) := by
    ext i j
    fin_cases i <;> fin_cases j
    ·
      have hij :=
        congrArg
          (fun M =>
            M 0 0)
          (Wedge2Formalization.N4.mul_J_add_J_mul_transpose
            (k := k)
            (B := (!![x, y; z, -x] :
              Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)))
      simpa using hij
    ·
      have hij :=
        congrArg
          (fun M =>
            M 0 1)
          (Wedge2Formalization.N4.mul_J_add_J_mul_transpose
            (k := k)
            (B := (!![x, y; z, -x] :
              Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)))
      simpa using hij
    ·
      have hij :=
        congrArg
          (fun M =>
            M 1 0)
          (Wedge2Formalization.N4.mul_J_add_J_mul_transpose
            (k := k)
            (B := (!![x, y; z, -x] :
              Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)))
      simpa using hij
    ·
      have hij :=
        congrArg
          (fun M =>
            M 1 1)
          (Wedge2Formalization.N4.mul_J_add_J_mul_transpose
            (k := k)
            (B := (!![x, y; z, -x] :
              Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)))
      simpa using hij
  have htrace :
      (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) *
          Wedge2Formalization.N4.J * (!![x, y; z, -x] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)ᵀ +
        (!![x, y; z, -x] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) *
          Wedge2Formalization.N4.J *
          (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)ᵀ =
        (0 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) := by
    simpa [Matrix.one_mul, Matrix.mul_one, add_comm] using htrace'
  simpa [U, Wedge2Formalization.N5.scalarBlock] using
    Wedge2Formalization.N5OnePoint.radOnePoint_pointwise_family
      (k := k) (a := (1 : k)) (C := C)
      (A := (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))
      (B := (!![x, y; z, -x] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))
      (hA := by simp) (hB := htrace)

theorem Levi_pointwise
    (u : k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1) :
    Wedge2Formalization.N5OnePoint.FixesRadOnePointPairBivector (Levi (k := k) u A) := by
  simpa [Levi, Wedge2Formalization.N5.scalarBlock] using
    Wedge2Formalization.N5OnePoint.radOnePoint_pointwise_family
      (k := k) (a := u) (C := 0) (A := A) (B := 0) (hA := hA) (hB := by simp)

theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ A B₁ : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            A * Wedge2Formalization.N4.J * B₁ᵀ + B₁ * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5.scalarBlock (k := k) u)
                0
                C
                (Matrix.fromBlocks A B₁ 0 A) := by
  constructor
  · rintro ⟨aM, R, C, A, B₁, C₁, D, hshape, hg⟩
    rcases hshape with ⟨hR, hA, hC, hD, hrel⟩
    refine ⟨aM 0 0, C, A, B₁, hA, hrel, ?_⟩
    rw [hg, hR, hC, hD]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [Wedge2Formalization.N5.scalarBlock]
  · rintro ⟨u, C, A, B₁, hA, hrel, rfl⟩
    refine ⟨Wedge2Formalization.N5.scalarBlock (k := k) u, 0, C, A, B₁, 0, A, ?_, rfl⟩
    exact ⟨rfl, hA, rfl, rfl, hrel⟩

/-- The radical scalar, radical unipotent family, and embedded `n = 4` Levi multiply to
the raw pointwise block shape. -/
private theorem scale_mul_U
    (u x y z : k)
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k) :
    Scale (k := k) u * U (k := k) x y z C =
      Matrix.fromBlocks
        (Wedge2Formalization.N5.scalarBlock (k := k) u)
        0
        C
        (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z) := by
  rw [Scale, U, Matrix.fromBlocks_multiply]
  simp [Wedge2Formalization.N5.scalarBlock]

private theorem embedded_mul_coreLevi
    (u : k)
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k)
    (H : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix.fromBlocks
        (Wedge2Formalization.N5.scalarBlock (k := k) u)
        0
        C
        H *
      CoreLevi (k := k) A =
      Matrix.fromBlocks
        (Wedge2Formalization.N5.scalarBlock (k := k) u)
        0
        C
        (H * Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) := by
  have hC :
      C * Wedge2Formalization.N5.scalarBlock (k := k) (1 : k) = C := by
    ext i j
    fin_cases j
    simp [Wedge2Formalization.N5.scalarBlock, Matrix.mul_apply]
  rw [CoreLevi, Matrix.fromBlocks_multiply, hC]
  simp [Wedge2Formalization.N5.scalarBlock]

private theorem scale_mul_U_mul_coreLevi
    (u x y z : k)
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Scale (k := k) u * U (k := k) x y z C * CoreLevi (k := k) A =
      Matrix.fromBlocks
        (Wedge2Formalization.N5.scalarBlock (k := k) u)
        0
        C
        (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
          Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) := by
  calc
    Scale (k := k) u * U (k := k) x y z C * CoreLevi (k := k) A =
        Matrix.fromBlocks
          (Wedge2Formalization.N5.scalarBlock (k := k) u)
          0
          C
          (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z) *
        CoreLevi (k := k) A := by
          rw [scale_mul_U]
    _ =
        Matrix.fromBlocks
          (Wedge2Formalization.N5.scalarBlock (k := k) u)
          0
          C
          ((Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z) *
            Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) := by
          rw [embedded_mul_coreLevi]

/-- The embedded repeated-support core, i.e. the `U_3 \rtimes SL_2(k)` factor
from Appendix A, `n = 4`, row 1. -/
def OnePointCore :
    Set
      (Matrix Wedge2Formalization.N4.V
        Wedge2Formalization.N4.V k) :=
  Wedge2Formalization.Paper.N4.Row1.K (k := k)

/-- Secondary kernel-membership statement through the embedded `n = 4`, row 1 core. -/
theorem mem_K_core_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ H : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k,
            H ∈ OnePointCore (k := k) ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5.scalarBlock (k := k) u)
                0
                C
                H := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨u, C, A, B₁, hA, hrel, rfl⟩
    refine ⟨u, C, Matrix.fromBlocks A B₁ 0 A, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks A B₁ 0 A)).2
        ⟨A, B₁, hA, hrel, rfl⟩
  · rintro ⟨u, C, H, hH, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := H)).1 hH with
      ⟨A, B₁, hA, hrel, rfl⟩
    exact
      (mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks
          (Wedge2Formalization.N5.scalarBlock (k := k) u)
          0
          C
          (Matrix.fromBlocks A B₁ 0 A))).2
        ⟨u, C, A, B₁, hA, hrel, rfl⟩

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u x y z : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            g =
              Scale (k := k) u * U (k := k) x y z C * CoreLevi (k := k) A := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨u, C, A, B₁, hA, hrel, hEq⟩
    have hcore :
        Matrix.fromBlocks A B₁ 0 A ∈ Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks A B₁ 0 A)).2
        ⟨A, B₁, hA, hrel, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row1.mem_K_iff
        (k := k)
        (g := Matrix.fromBlocks A B₁ 0 A)).1 hcore with
      ⟨A', x, y, z, hA', hcoreEq⟩
    refine ⟨u, x, y, z, C, A', hA', ?_⟩
    calc
      g =
        Matrix.fromBlocks
          (Wedge2Formalization.N5.scalarBlock (k := k) u)
          0
          C
          (Matrix.fromBlocks A B₁ 0 A) := hEq
      _ =
        Matrix.fromBlocks
          (Wedge2Formalization.N5.scalarBlock (k := k) u)
          0
          C
          (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
            Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A') := by
              rw [hcoreEq]
      _ = Scale (k := k) u * U (k := k) x y z C * CoreLevi (k := k) A' := by
              rw [scale_mul_U_mul_coreLevi]
  · rintro ⟨u, x, y, z, C, A, hA, hEq⟩
    have hcore :
        Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
            Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A ∈
          Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
      (Wedge2Formalization.Paper.N4.Row1.mem_K_iff
        (k := k)
        (g := Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
          Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)).2
        ⟨A, x, y, z, hA, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
          Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)).1 hcore with
      ⟨A', B₁, hA', hrel, hcoreEq⟩
    apply (mem_K_shape_iff (k := k) (g := g)).2
    refine ⟨u, C, A', B₁, hA', hrel, ?_⟩
    rw [hEq, scale_mul_U_mul_coreLevi, hcoreEq]

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_7 \rtimes (G_m(k) \times SL_2(k))`, written in the displayed
`Scale * U * CoreLevi` factorization. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u x y z : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            g = Scale (k := k) u * U (k := k) x y z C * CoreLevi (k := k) A :=
  mem_K_iff (k := k) (g := g)

/-- The coefficient matrix attached to the standard lift belongs to the Borel family. -/
theorem coeff_mem_Q
    (a b : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b ∈ Q (k := k) :=
  ⟨a, b, ha, rfl⟩

/-- The coefficient matrix attached to the standard lift lies in the exact projective Borel
family. -/
theorem coeff_mem_Qproj
    (a b : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b ∈ Qproj (k := k) := by
  constructor
  · simp [Wedge2Formalization.N4PaperSummary.row1_borelCoeff]
  · simp [Wedge2Formalization.N4PaperSummary.row1_borelCoeff, Matrix.det_fin_two, ha]

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N5.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N5Summary.radOnePoint_borel_lift_action
        (k := k) (a := a) (b := b) (t := (1 : k)) ha (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift, Wedge2Formalization.N4PaperSummary.row1_borelCoeff] using
      (Wedge2Formalization.N5Summary.radOnePoint_borel_lift_action
        (k := k) (a := a) (b := b) (t := (1 : k)) ha (C := 0)).2

/-- Rank drop on the repeated-support quotient pencil occurs exactly at the repeated point. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det
      (a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
        b • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 ↔
      a = 0 :=
  Wedge2Formalization.N4Summary.onePoint_det_zero_iff (k := k) (a := a) (b := b)

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          b • (Wedge2Formalization.N4.onePointRep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N5OnePoint.radOnePointRep₁,
      Wedge2Formalization.N5OnePoint.radOnePointRep₂] using h'
  have ha :
      a = 0 := by
    have hdet :
        Matrix.det
          (a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
            b • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 := by
      rw [h'']
      simpa using (Matrix.det_zero (n := Wedge2Formalization.N4.V) (R := k))
    exact
      (Wedge2Formalization.N4Summary.onePoint_det_zero_iff (k := k) (a := a) (b := b)).1 hdet
  have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h''
  have hb :
      b = 0 := by
    simpa [ha, Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.J,
      Matrix.add_apply, Matrix.smul_apply] using h01
  exact ⟨ha, hb⟩

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  rcases rep_pair_independent (k := k) hcomb with ⟨_, h1⟩
  exact one_ne_zero h1

private theorem lowerRight_action
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k)
    {α β γ δ : k}
    (h1 :
      Wedge2Formalization.N5.ActBivector
        (Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)) g
        =
      α • (Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)) +
        β • (Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k)))
    (h2 :
      Wedge2Formalization.N5.ActBivector
        (Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k)) g
        =
      γ • (Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)) +
        δ • (Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k))) :
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.onePointRep₁ (k := k)) g.toBlocks₂₂ =
        α • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.onePointRep₂ (k := k)) ∧
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.onePointRep₂ (k := k)) g.toBlocks₂₂ =
        γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
  have h1' := congrArg Matrix.toBlocks₂₂ h1
  have h2' := congrArg Matrix.toBlocks₂₂ h2
  have hrhs1 :
      Matrix.toBlocks₂₂
        (α • (Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)) +
          β • (Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k))) =
        α • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
    ext i j
    change
      (α • (Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)) +
        β • (Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (α • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
        β • (Wedge2Formalization.N4.onePointRep₂ (k := k))) i j
    simp [Wedge2Formalization.N5OnePoint.radOnePointRep₁, Wedge2Formalization.N5OnePoint.radOnePointRep₂,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
  have hrhs2 :
      Matrix.toBlocks₂₂
        (γ • (Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k))) =
        γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
    ext i j
    change
      (γ • (Wedge2Formalization.N5OnePoint.radOnePointRep₁ (k := k)) +
        δ • (Wedge2Formalization.N5OnePoint.radOnePointRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
        δ • (Wedge2Formalization.N4.onePointRep₂ (k := k))) i j
    simp [Wedge2Formalization.N5OnePoint.radOnePointRep₁, Wedge2Formalization.N5OnePoint.radOnePointRep₂,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
  constructor
  · simpa [act_embedded_toBlocks₂₂, Wedge2Formalization.N4.ActBivector,
      Wedge2Formalization.N5OnePoint.radOnePointRep₁, Wedge2Formalization.N5OnePoint.radOnePointRep₂,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Matrix.add_apply, Matrix.smul_apply] using h1'
  · simpa [act_embedded_toBlocks₂₂, Wedge2Formalization.N4.ActBivector,
      Wedge2Formalization.N5OnePoint.radOnePointRep₁, Wedge2Formalization.N5OnePoint.radOnePointRep₂,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Matrix.add_apply, Matrix.smul_apply] using h2'

theorem quotient_image
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N5OnePoint.PreservesRadOnePointSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N5.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N5.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N5.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  rcases lowerRight_action (k := k) (g := g) h1 h2 with ⟨hD1, hD2⟩
  have hrep₂_det :
      Matrix.det (Wedge2Formalization.N4.onePointRep₂ (k := k)) = 0 := by
    rw [Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Matrix.det_fromBlocks_zero₂₁]
    simpa using (Matrix.det_zero (n := Wedge2Formalization.N4.I) (R := k))
  have hrep₂_det_zero :
      Matrix.det
        (γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 := by
    rw [← hD2, Wedge2Formalization.N4Summary.actBivector_det, hrep₂_det]
    simp
  have hγ : γ = 0 :=
    Wedge2Formalization.N4Summary.onePoint_action_case_of_rankDrop
      (k := k) (γ := γ) (δ := δ) hrep₂_det_zero
  have hδ : δ ≠ 0 := by
    intro hδ0
    have hzero :
        Wedge2Formalization.N5.ActBivector (rep₂ (k := k)) g = 0 := by
      simpa [hγ, hδ0] using h2'
    have horig :=
      (Wedge2Formalization.N5.actBivector_eq_zero_iff_of_det_ne_zero
        (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
    exact rep₂_ne_zero (k := k) horig
  have hα : α ≠ 0 := by
    intro hα0
    have hzero :
        Wedge2Formalization.N5.ActBivector
          (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g = 0 := by
      calc
        Wedge2Formalization.N5.ActBivector
            (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g
            =
          δ • Wedge2Formalization.N5.ActBivector (rep₁ (k := k)) g +
            (-β) • Wedge2Formalization.N5.ActBivector (rep₂ (k := k)) g := by
              rw [Wedge2Formalization.N5.actBivector_add,
                Wedge2Formalization.N5.actBivector_smul,
                Wedge2Formalization.N5.actBivector_smul]
        _ = 0 := by
              rw [h1', h2']
              ext i j
              simp [rep₁, rep₂, hα0, hγ]
              ring
    have horig :=
      (Wedge2Formalization.N5.actBivector_eq_zero_iff_of_det_ne_zero
        (k := k) (Ω := δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) (g := g) hg).1 hzero
    have hcomb_ne :
        δ • rep₁ (k := k) + (-β) • rep₂ (k := k) ≠ 0 := by
      intro hcomb
      rcases rep_pair_independent (k := k) hcomb with ⟨hδ', hβ'⟩
      exact hδ hδ'
    exact hcomb_ne horig
  refine ⟨!![α, β; 0, δ], ?_, ?_⟩
  · constructor
    · simp
    · simp [Matrix.det_fin_two, hα, hδ]
  · constructor
    · simpa [ActsOnOrderedPair] using h1
    · simpa [ActsOnOrderedPair, hγ] using h2

end Row4

/-! Appendix A, `n = 5`, row 5.
Representative `⟨e₂∧e₃ + e₄∧e₅, e₄∧e₅⟩`.
Divisor `[a] + [b]`.
Claimed stabilizer:
`K_L = U_4 \rtimes (G_m(k) \times SL_2(k) \times SL_2(k))`, exact quotient family `Q_L = N_T`.
-/
namespace Row5

def rep₁ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Wedge2Formalization.N5.radSplitRep₁ (k := k)

def rep₂ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Wedge2Formalization.N5.radSplitRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 5:
`e₂∧e₃ + e₄∧e₅`. -/
def paperRep₁ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row2.rep₁ (k := k))

/-- Literal second basis vector from Appendix A, row 5:
`e₄∧e₅`. -/
def paperRep₂ : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row2.rep₂ (k := k))

/-- The paper representative already agrees with the internal working representative. -/
def paperChange : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  1

/-- Transport of the literal paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N5.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N5.radSplitRep₁,
    Wedge2Formalization.Paper.N4.Row2.rep₁, Wedge2Formalization.N4.splitRep₁,
    Wedge2Formalization.N5.ActBivector]

/-- Transport of the literal paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N5.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N5.radSplitRep₂,
    Wedge2Formalization.Paper.N4.Row2.rep₂, Wedge2Formalization.N4.splitRep₂,
    Wedge2Formalization.N5.ActBivector]

/-- Exact pointwise kernel family for the radical split row. -/
def K : Set (Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :=
  { g |
      ∃ a : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.I k,
        ∃ R : Matrix Wedge2Formalization.N5.I Wedge2Formalization.N5.W k,
          ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
            ∃ A B₁ C₁ D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
              Wedge2Formalization.N5PaperSummary.row5_pointwiseShape (k := k) R A B₁ C₁ D ∧
              g = Matrix.fromBlocks a R C (Matrix.fromBlocks A B₁ C₁ D) }

/-- The displayed radical unipotent family on the radical split row. -/
def U
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) (1 : k))
    0
    C
    (1 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)

/-- The displayed Levi/core family on the radical split row. -/
def Levi
    (u : k)
    (A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) u)
    0
    0
    (Matrix.fromBlocks A 0 0 E)

/-- Exact coefficient-side split normalizer family in the ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N4.Row2.Qproj (k := k)

/-- Chosen torus lift. -/
def torusLift
    (u v : k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) (1 : k))
    0
    0
    (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
      0 0 (!![v, 0; 0, 1] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))

/-- Chosen swap-coset lift. -/
def swapLift
    (u v : k) :
    Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k :=
  let A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k := !![u, 0; 0, 1]
  let D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k := !![v, 0; 0, 1]
  Matrix.fromBlocks
    (Wedge2Formalization.N5.scalarBlock (k := k) (1 : k))
    0
    0
    (Matrix.fromBlocks A 0 0 D * Wedge2Formalization.N4.splitSwap (k := k))

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N5.FixesRadSplitPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    refine
      ⟨g.toBlocks₁₁, g.toBlocks₁₂, g.toBlocks₂₁,
        g.toBlocks₂₂.toBlocks₁₁, g.toBlocks₂₂.toBlocks₁₂,
        g.toBlocks₂₂.toBlocks₂₁, g.toBlocks₂₂.toBlocks₂₂, ?_, ?_⟩
    · have hg' :
        Wedge2Formalization.N5.FixesRadSplitPairBivector
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁
            (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
              g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂)) := by
        have hEq22 :
            Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
              g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂ = g.toBlocks₂₂ := by
          simpa using Matrix.fromBlocks_toBlocks g.toBlocks₂₂
        have hEq :
            Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁
                (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
                g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂) = g := by
          rw [hEq22]
          exact Matrix.fromBlocks_toBlocks g
        simpa [hEq] using hg
      exact
        (Wedge2Formalization.N5PaperSummary.row5_pointwise_iff
          (k := k)
          (a := g.toBlocks₁₁)
          (R := g.toBlocks₁₂)
          (C := g.toBlocks₂₁)
          (A := g.toBlocks₂₂.toBlocks₁₁)
          (B₁ := g.toBlocks₂₂.toBlocks₁₂)
          (C₁ := g.toBlocks₂₂.toBlocks₂₁)
          (D := g.toBlocks₂₂.toBlocks₂₂)).1 hg'
    · have hEq22 :
          Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
            g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂ = g.toBlocks₂₂ := by
        simpa using Matrix.fromBlocks_toBlocks g.toBlocks₂₂
      rw [hEq22]
      exact (Matrix.fromBlocks_toBlocks g).symm
  · rintro ⟨a, R, C, A, B₁, C₁, D, hshape, rfl⟩
    exact
      (Wedge2Formalization.N5PaperSummary.row5_pointwise_iff
        (k := k) (a := a) (R := R) (C := C) (A := A) (B₁ := B₁) (C₁ := C₁) (D := D)).2 hshape

theorem U_pointwise
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k) :
    Wedge2Formalization.N5.FixesRadSplitPairBivector (U (k := k) C) := by
  simpa [U, Wedge2Formalization.N5.scalarBlock] using
    Wedge2Formalization.N5.radSplit_pointwise_family
      (k := k) (a := (1 : k)) (C := C)
      (A := (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))
      (E := (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))
      (hA := by simp) (hE := by simp)

theorem Levi_pointwise
    (u : k)
    (A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1)
    (hE : E.det = 1) :
    Wedge2Formalization.N5.FixesRadSplitPairBivector (Levi (k := k) u A E) := by
  simpa [Levi, Wedge2Formalization.N5.scalarBlock] using
    Wedge2Formalization.N5.radSplit_pointwise_family
      (k := k) (a := u) (C := 0) (A := A) (E := E) hA hE

theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5.scalarBlock (k := k) u)
                0
                C
                (Matrix.fromBlocks A 0 0 E) := by
  constructor
  · rintro ⟨aM, R, C, A, B₁, C₁, D, hshape, hg⟩
    rcases hshape with ⟨hR, hA, hB, hC, hD⟩
    refine ⟨aM 0 0, C, A, D, hA, hD, ?_⟩
    rw [hg, hR, hB, hC]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [Wedge2Formalization.N5.scalarBlock]
  · rintro ⟨u, C, A, E, hA, hE, rfl⟩
    refine ⟨Wedge2Formalization.N5.scalarBlock (k := k) u, 0, C, A, 0, 0, E, ?_, rfl⟩
    exact ⟨rfl, hA, rfl, rfl, hE⟩

/-- The displayed Levi family and radical family multiply to the raw pointwise block
shape. -/
private theorem Levi_mul_U_eq_shape
    (u : k)
    (C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k)
    (A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Levi (k := k) u A E * U (k := k) C =
      Matrix.fromBlocks
        (Wedge2Formalization.N5.scalarBlock (k := k) u)
        0
        (Matrix.fromBlocks A 0 0 E * C)
        (Matrix.fromBlocks A 0 0 E) := by
  rw [U, Levi, Matrix.fromBlocks_multiply]
  ext i j
  fin_cases i <;> fin_cases j <;> simp [Wedge2Formalization.N5.scalarBlock]

/-- The embedded split core, i.e. the `SL_2(k) \times SL_2(k)` factor from
Appendix A, `n = 4`, row 2. -/
def SplitCore :
    Set
      (Matrix Wedge2Formalization.N4.V
        Wedge2Formalization.N4.V k) :=
  Wedge2Formalization.Paper.N4.Row2.K (k := k)

/-- Secondary kernel-membership statement using the embedded `n = 4`, row 2 core. -/
theorem mem_K_core_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ H : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k,
            H ∈ SplitCore (k := k) ∧
            g =
              Matrix.fromBlocks
                (Wedge2Formalization.N5.scalarBlock (k := k) u)
                0
                C
                H := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨u, C, A, E, hA, hE, rfl⟩
    refine ⟨u, C, Matrix.fromBlocks A 0 0 E, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N4.Row2.mem_K_iff
        (k := k)
        (g := Matrix.fromBlocks A 0 0 E)).2
        ⟨A, E, hA, hE, rfl⟩
  · rintro ⟨u, C, H, hH, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row2.mem_K_iff
        (k := k)
        (g := H)).1 hH with
      ⟨A, E, hA, hE, rfl⟩
    exact
      (mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks
          (Wedge2Formalization.N5.scalarBlock (k := k) u)
          0
          C
          (Matrix.fromBlocks A 0 0 E))).2
        ⟨u, C, A, E, hA, hE, rfl⟩

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g = Levi (k := k) u A E * U (k := k) C := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with ⟨u, C, A, E, hA, hE, hEq⟩
    let H : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
      Matrix.fromBlocks A 0 0 E
    have hdetH : Matrix.det H = 1 := by
      dsimp [H]
      rw [Matrix.det_fromBlocks_zero₂₁]
      simp [hA, hE]
    have hunitH : IsUnit (Matrix.det H) := by rw [hdetH]; exact isUnit_one
    have hHC : H * (H⁻¹ * C) = C := by
      simpa using Matrix.mul_nonsing_inv_cancel_left H C hunitH
    refine ⟨u, H⁻¹ * C, A, E, hA, hE, ?_⟩
    rw [Levi_mul_U_eq_shape (k := k) (u := u) (C := H⁻¹ * C) (A := A) (E := E)]
    simpa [H, hHC] using hEq
  · rintro ⟨u, C, A, E, hA, hE, hEq⟩
    apply (mem_K_shape_iff (k := k) (g := g)).2
    refine ⟨u, Matrix.fromBlocks A 0 0 E * C, A, E, hA, hE, ?_⟩
    rw [hEq, Levi_mul_U_eq_shape (k := k) (u := u) (C := C) (A := A) (E := E)]

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_4 \rtimes (G_m(k) \times SL_2(k) \times SL_2(k))`, written in the
displayed `Levi * U` factorization. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k) :
    g ∈ K (k := k) ↔
      ∃ u : k,
        ∃ C : Matrix Wedge2Formalization.N5.W Wedge2Formalization.N5.I k,
          ∃ A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g = Levi (k := k) u A E * U (k := k) C :=
  mem_K_iff (k := k) (g := g)

theorem torus_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N5.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (torusLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_torusCoeff (k := k) u v) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, torusLift] using
      (Wedge2Formalization.N5Summary.radSplit_torus_lift_action
        (k := k) (a := (1 : k)) (u := u) (v := v) (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, torusLift, Wedge2Formalization.N4PaperSummary.row2_torusCoeff] using
      (Wedge2Formalization.N5Summary.radSplit_torus_lift_action
        (k := k) (a := (1 : k)) (u := u) (v := v) (C := 0)).2

theorem swap_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N5.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (swapLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_swapCoeff (k := k) u v) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, swapLift] using
      (Wedge2Formalization.N5Summary.radSplit_swapCoset_lift_action
        (k := k) (a := (1 : k)) (u := u) (v := v) (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, swapLift, Wedge2Formalization.N4PaperSummary.row2_swapCoeff] using
      (Wedge2Formalization.N5Summary.radSplit_swapCoset_lift_action
        (k := k) (a := (1 : k)) (u := u) (v := v) (C := 0)).2

/-- Rank drop on the split-support quotient pencil occurs exactly at the two support points. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det
      (a • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        b • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  Wedge2Formalization.N4Summary.split_det_zero_iff (k := k) (a := a) (b := b)

private theorem rep_pair_linearCombination_ne_zero
    {a b : k}
    (h : a ≠ 0 ∨ b ≠ 0) :
    a • rep₁ (k := k) + b • rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h' := congrArg Matrix.toBlocks₂₂ hzero
  have h'' :
      a • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        b • (Wedge2Formalization.N4.splitRep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N5.radSplitRep₁, Wedge2Formalization.N5.radSplitRep₂] using h'
  have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h''
  have ha : a = 0 := by
    simpa [Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J, Matrix.add_apply,
      Matrix.smul_apply] using h01
  have h23 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h''
  have hb : b = 0 := by
    simpa [ha, Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J, Matrix.add_apply,
      Matrix.smul_apply] using h23
  exact h.elim (fun ha' => ha' ha) (fun hb' => hb' hb)

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  have h :=
    rep_pair_linearCombination_ne_zero (k := k) (a := (0 : k)) (b := (1 : k)) (Or.inr one_ne_zero)
  simpa using h

private theorem lowerRight_action
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k)
    {α β γ δ : k}
    (h1 :
      Wedge2Formalization.N5.ActBivector
        (Wedge2Formalization.N5.radSplitRep₁ (k := k)) g
        =
      α • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
        β • (Wedge2Formalization.N5.radSplitRep₂ (k := k)))
    (h2 :
      Wedge2Formalization.N5.ActBivector
        (Wedge2Formalization.N5.radSplitRep₂ (k := k)) g
        =
      γ • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
        δ • (Wedge2Formalization.N5.radSplitRep₂ (k := k))) :
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₁ (k := k)) g.toBlocks₂₂ =
        α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.splitRep₂ (k := k)) ∧
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₂ (k := k)) g.toBlocks₂₂ =
        γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
  have h1' := congrArg Matrix.toBlocks₂₂ h1
  have h2' := congrArg Matrix.toBlocks₂₂ h2
  have hrhs1 :
      Matrix.toBlocks₂₂
        (α • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
          β • (Wedge2Formalization.N5.radSplitRep₂ (k := k))) =
        α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    ext i j
    change
      (α • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
        β • (Wedge2Formalization.N5.radSplitRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        β • (Wedge2Formalization.N4.splitRep₂ (k := k))) i j
    simp [Wedge2Formalization.N5.radSplitRep₁, Wedge2Formalization.N5.radSplitRep₂,
      Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
    ring
  have hrhs2 :
      Matrix.toBlocks₂₂
        (γ • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N5.radSplitRep₂ (k := k))) =
        γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    ext i j
    change
      (γ • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
        δ • (Wedge2Formalization.N5.radSplitRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        δ • (Wedge2Formalization.N4.splitRep₂ (k := k))) i j
    simp [Wedge2Formalization.N5.radSplitRep₁, Wedge2Formalization.N5.radSplitRep₂,
      Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
    ring
  constructor
  · calc
      Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₁ (k := k)) g.toBlocks₂₂
          =
        Matrix.toBlocks₂₂
          (Wedge2Formalization.N5.ActBivector
            (Wedge2Formalization.N5.radSplitRep₁ (k := k)) g) := by
              symm
              simpa [Wedge2Formalization.N5.radSplitRep₁] using
                act_embedded_toBlocks₂₂
                  (k := k)
                  (Ω := Wedge2Formalization.N4.splitRep₁ (k := k))
                  (g := g)
      _ = Matrix.toBlocks₂₂
            (α • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
              β • (Wedge2Formalization.N5.radSplitRep₂ (k := k))) := by
                simpa using h1'
      _ =
          α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            β • (Wedge2Formalization.N4.splitRep₂ (k := k)) := hrhs1
  · calc
      Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₂ (k := k)) g.toBlocks₂₂
          =
        Matrix.toBlocks₂₂
          (Wedge2Formalization.N5.ActBivector
            (Wedge2Formalization.N5.radSplitRep₂ (k := k)) g) := by
              symm
              simpa [Wedge2Formalization.N5.radSplitRep₂] using
                act_embedded_toBlocks₂₂
                  (k := k)
                  (Ω := Wedge2Formalization.N4.splitRep₂ (k := k))
                  (g := g)
      _ = Matrix.toBlocks₂₂
            (γ • (Wedge2Formalization.N5.radSplitRep₁ (k := k)) +
              δ • (Wedge2Formalization.N5.radSplitRep₂ (k := k))) := by
                simpa using h2'
      _ =
          γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            δ • (Wedge2Formalization.N4.splitRep₂ (k := k)) := hrhs2

theorem quotient_image
    (g : Matrix Wedge2Formalization.N5.V Wedge2Formalization.N5.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N5.PreservesRadSplitSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N5.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N5.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N5.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  rcases lowerRight_action (k := k) (g := g) h1 h2 with ⟨hD1, hD2⟩
  have hrep2_det :
      Matrix.det (Wedge2Formalization.N4.splitRep₂ (k := k)) = 0 := by
    simpa using
      (Wedge2Formalization.N4Summary.split_det_zero_iff
        (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl)
  have hrep2_det_zero :
      Matrix.det
        (γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    rw [← hD2, Wedge2Formalization.N4Summary.actBivector_det, hrep2_det]
    simp
  have hω12_det :
      Matrix.det
        ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
          (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    have h :=
      (Wedge2Formalization.N4Summary.split_det_zero_iff
        (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (by ring))
    simpa [sub_eq_add_neg] using h
  have hω12_action :
      Wedge2Formalization.N4.ActBivector
        ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
          (Wedge2Formalization.N4.splitRep₂ (k := k))) g.toBlocks₂₂
        =
      (α - γ) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        (β - δ) • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    calc
      Wedge2Formalization.N4.ActBivector
          ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
            (Wedge2Formalization.N4.splitRep₂ (k := k))) g.toBlocks₂₂
          =
        Wedge2Formalization.N4.ActBivector
          ((Wedge2Formalization.N4.splitRep₁ (k := k)) +
            (-1 : k) • (Wedge2Formalization.N4.splitRep₂ (k := k))) g.toBlocks₂₂ := by
              simp [sub_eq_add_neg]
      _ =
          Wedge2Formalization.N4.ActBivector
              (Wedge2Formalization.N4.splitRep₁ (k := k)) g.toBlocks₂₂ +
            (-1 : k) •
              Wedge2Formalization.N4.ActBivector
                (Wedge2Formalization.N4.splitRep₂ (k := k)) g.toBlocks₂₂ := by
                  rw [Wedge2Formalization.N4Summary.actBivector_add,
                    Wedge2Formalization.N4Summary.actBivector_smul]
      _ =
          (α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            β • (Wedge2Formalization.N4.splitRep₂ (k := k))) +
          (-1 : k) •
            (γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
              δ • (Wedge2Formalization.N4.splitRep₂ (k := k))) := by
                rw [hD1, hD2]
      _ =
          (α - γ) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            (β - δ) • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
              ext i j
              simp [sub_eq_add_neg]
              ring
  have hω12_det_zero :
      Matrix.det
        ((α - γ) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          (β - δ) • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    rw [← hω12_action, Wedge2Formalization.N4Summary.actBivector_det, hω12_det]
    simp
  have hcoeff_det_ne : α * δ - β * γ ≠ 0 := by
    intro hdet
    by_cases hbd : β = 0 ∧ δ = 0
    · by_cases hca : γ = 0 ∧ α = 0
      · have hrep2_zero :
          Wedge2Formalization.N5.ActBivector (rep₂ (k := k)) g = 0 := by
          simpa [hca.1, hbd.2] using h2'
        have horig :=
          (Wedge2Formalization.N5.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hrep2_zero
        exact rep₂_ne_zero (k := k) horig
      · have hnonzero : γ ≠ 0 ∨ α ≠ 0 := by
          by_cases hγ : γ = 0
          · right
            intro hα
            exact hca ⟨hγ, hα⟩
          · exact Or.inl hγ
        have hzero2 :
            Wedge2Formalization.N5.ActBivector
              (γ • rep₁ (k := k) - α • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N5.ActBivector
                (γ • rep₁ (k := k) - α • rep₂ (k := k)) g
                =
              γ • Wedge2Formalization.N5.ActBivector (rep₁ (k := k)) g +
                (-α) • Wedge2Formalization.N5.ActBivector (rep₂ (k := k)) g := by
                  rw [sub_eq_add_neg, Wedge2Formalization.N5.actBivector_add,
                    Wedge2Formalization.N5.actBivector_smul]
                  rw [show -(α • rep₂ (k := k)) = (-α) • rep₂ (k := k) by simp]
                  rw [Wedge2Formalization.N5.actBivector_smul]
            _ =
                γ • (α • rep₁ (k := k) + β • rep₂ (k := k)) +
                  (-α) • (γ • rep₁ (k := k) + δ • rep₂ (k := k)) := by
                    rw [h1', h2']
            _ = 0 := by
                  ext i j
                  simp [hbd.1, hbd.2, sub_eq_add_neg]
                  ring
        have horig :=
          (Wedge2Formalization.N5.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := γ • rep₁ (k := k) - α • rep₂ (k := k)) (g := g) hg).1 hzero2
        have hcomb_ne :
            γ • rep₁ (k := k) - α • rep₂ (k := k) ≠ 0 := by
          have hcomb_ne' :
              γ ≠ 0 ∨ -α ≠ 0 := by
            rcases hnonzero with hγ | hα
            · exact Or.inl hγ
            · exact Or.inr (neg_ne_zero.mpr hα)
          have := rep_pair_linearCombination_ne_zero (k := k) (a := γ) (b := -α) hcomb_ne'
          simpa [sub_eq_add_neg] using this
        exact hcomb_ne horig
    · have hnonzero : β ≠ 0 ∨ δ ≠ 0 := by
        by_cases hβ : β = 0
        · right
          intro hδ
          exact hbd ⟨hβ, hδ⟩
        · exact Or.inl hβ
      have hzero1 :
          Wedge2Formalization.N5.ActBivector
            (δ • rep₁ (k := k) - β • rep₂ (k := k)) g = 0 := by
        calc
          Wedge2Formalization.N5.ActBivector
              (δ • rep₁ (k := k) - β • rep₂ (k := k)) g
              =
            δ • Wedge2Formalization.N5.ActBivector (rep₁ (k := k)) g +
              (-β) • Wedge2Formalization.N5.ActBivector (rep₂ (k := k)) g := by
                rw [sub_eq_add_neg, Wedge2Formalization.N5.actBivector_add,
                  Wedge2Formalization.N5.actBivector_smul]
                rw [show -(β • rep₂ (k := k)) = (-β) • rep₂ (k := k) by simp]
                rw [Wedge2Formalization.N5.actBivector_smul]
          _ =
              δ • (α • rep₁ (k := k) + β • rep₂ (k := k)) +
                (-β) • (γ • rep₁ (k := k) + δ • rep₂ (k := k)) := by
                  rw [h1', h2']
          _ = (α * δ - β * γ) • rep₁ (k := k) := by
                ext i j
                simp [sub_eq_add_neg]
                ring
          _ = 0 := by
                simp [hdet]
      have horig :=
        (Wedge2Formalization.N5.actBivector_eq_zero_iff_of_det_ne_zero
          (k := k) (Ω := δ • rep₁ (k := k) - β • rep₂ (k := k)) (g := g) hg).1 hzero1
      have hcomb_ne :
          δ • rep₁ (k := k) - β • rep₂ (k := k) ≠ 0 := by
        have hcomb_ne' :
            δ ≠ 0 ∨ -β ≠ 0 := by
          rcases hnonzero with hβ | hδ
          · exact Or.inr (neg_ne_zero.mpr hβ)
          · exact Or.inl hδ
        have := rep_pair_linearCombination_ne_zero (k := k) (a := δ) (b := -β) hcomb_ne'
        simpa [sub_eq_add_neg] using this
      exact hcomb_ne horig
  have hM_mem :
      ((!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k)) ∈ Qproj (k := k) := by
    rcases
      Wedge2Formalization.N4Summary.split_action_cases_of_rankDrop
        (k := k)
        (α := α)
        (β := β)
        (γ := γ)
        (δ := δ)
        hcoeff_det_ne
        hrep2_det_zero
        hω12_det_zero with
      hdiag | hswap
    · rcases hdiag with ⟨hγ, hsum⟩
      left
      have hbeta : β = δ - α := by
        calc
          β = α + β - α := by ring
          _ = δ - α := by rw [hsum]
      refine ⟨α, δ, ?_⟩
      ext i j
      fin_cases i <;> fin_cases j
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, hbeta]
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, hbeta]
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, hγ]
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff]
    · rcases hswap with ⟨hαγ, hγδ⟩
      right
      refine ⟨α, α + β, ?_⟩
      have hδ : δ = -α := by
        calc
          δ = γ + δ - γ := by ring
          _ = 0 - γ := by rw [hγδ]
          _ = -γ := by ring
          _ = -α := by rw [hαγ]
      ext i j
      fin_cases i <;> fin_cases j
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff]
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff]
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff, hαγ]
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff, hδ]
  refine ⟨!![α, β; γ, δ], hM_mem, ?_⟩
  constructor
  · simpa [ActsOnOrderedPair] using h1
  · simpa [ActsOnOrderedPair] using h2

end Row5

end N5
end Paper
end Wedge2Formalization
