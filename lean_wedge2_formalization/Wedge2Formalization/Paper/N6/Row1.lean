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

/-! Appendix A, `n = 6`, row 1.
Representative `⟨e₁∧e₄ + e₂∧e₅ + e₃∧e₆, e₁∧e₅ + e₂∧e₆⟩`.
Divisor `3[a]`.
Claimed stabilizer:
`K_L = U_6 \rtimes SL_2(k)`, exact quotient family `Q_L = B`.
-/
namespace Row1

abbrev I := Wedge2Formalization.N6OnePointLong.I
abbrev V := Wedge2Formalization.N6OnePointLong.V

def rep₁ :
    Matrix V V k :=
  Wedge2Formalization.N6OnePointLong.rep₁ (k := k)

def rep₂ :
    Matrix V V k :=
  Wedge2Formalization.N6OnePointLong.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 1. -/
def paperRep₁ : Matrix V V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 1. -/
def paperRep₂ : Matrix V V k :=
  rep₂ (k := k)

/-- The row-1 paper representative already agrees with the internal working pair. -/
def paperChange : Matrix V V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6OnePointLong.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6OnePointLong.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6OnePointLong.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6OnePointLong.ActBivector]

private def diagN : Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N6OnePointLong.N (k := k))
    0
    0
    ((Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ)

/-- The anti-diagonal permutation matrix on the local `3`-space. -/
def flip : Matrix I I k :=
  !![(0 : k), 0, 1;
     0, 1, 0;
     1, 0, 0]

/-- Upper Toeplitz block corresponding to a truncated polynomial in `k[t]/(t^3)`. -/
def upperToeplitz
    (a₀ a₁ a₂ : k) : Matrix I I k :=
  !![a₀, a₁, a₂;
     0,  a₀, a₁;
     0,  0,  a₀]

/-- Lower Toeplitz block corresponding to a truncated polynomial in `k[t]/(t^3)`. -/
def lowerToeplitz
    (d₀ d₁ d₂ : k) : Matrix I I k :=
  !![d₀, 0,  0;
     d₁, d₀, 0;
     d₂, d₁, d₀]

private theorem upperToeplitz_transpose
    (a₀ a₁ a₂ : k) :
    (upperToeplitz (k := k) a₀ a₁ a₂)ᵀ = lowerToeplitz (k := k) a₀ a₁ a₂ := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [upperToeplitz, lowerToeplitz, Matrix.transpose_apply]

private theorem lowerToeplitz_transpose
    (d₀ d₁ d₂ : k) :
    (lowerToeplitz (k := k) d₀ d₁ d₂)ᵀ = upperToeplitz (k := k) d₀ d₁ d₂ := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [upperToeplitz, lowerToeplitz, Matrix.transpose_apply]

private theorem flip_transpose :
    (flip (k := k))ᵀ = flip (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [flip, Matrix.transpose_apply]

private theorem flip_mul_flip :
    flip (k := k) * flip (k := k) = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [flip, Matrix.mul_apply, Fin.sum_univ_three]

private theorem N_mul_flip :
    Wedge2Formalization.N6OnePointLong.N (k := k) * flip (k := k) =
      flip (k := k) * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [flip, Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply,
      Matrix.transpose_apply, Fin.sum_univ_three]

private theorem flip_mul_N_transpose_mul_flip :
    flip (k := k) * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * flip (k := k) =
      Wedge2Formalization.N6OnePointLong.N (k := k) := by
  calc
    flip (k := k) * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * flip (k := k) =
      (Wedge2Formalization.N6OnePointLong.N (k := k) * flip (k := k)) * flip (k := k) := by
        rw [show flip (k := k) * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ =
              Wedge2Formalization.N6OnePointLong.N (k := k) * flip (k := k) by
              simpa using (N_mul_flip (k := k)).symm]
    _ = Wedge2Formalization.N6OnePointLong.N (k := k) * (flip (k := k) * flip (k := k)) := by
        rw [Matrix.mul_assoc]
    _ = Wedge2Formalization.N6OnePointLong.N (k := k) * 1 := by
        rw [flip_mul_flip (k := k)]
    _ = Wedge2Formalization.N6OnePointLong.N (k := k) := by
        rw [Matrix.mul_one]

private theorem neg_rep₁_mul_rep₁ :
    (-(rep₁ (k := k))) * rep₁ (k := k) = 1 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [rep₁, Wedge2Formalization.N6OnePointLong.rep₁, Matrix.fromBlocks,
      Matrix.mul_apply, Fin.sum_univ_three]

private theorem rep₁_mul_neg_rep₁ :
    rep₁ (k := k) * (-(rep₁ (k := k))) = 1 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [rep₁, Wedge2Formalization.N6OnePointLong.rep₁, Matrix.fromBlocks,
      Matrix.mul_apply, Fin.sum_univ_three]

set_option maxHeartbeats 1000000 in
private theorem rep₂_mul_neg_rep₁ :
    rep₂ (k := k) * (-(rep₁ (k := k))) = diagN (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [rep₁, rep₂, diagN, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Wedge2Formalization.N6OnePointLong.N,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]

private theorem diagN_mul_rep₁ :
    diagN (k := k) * rep₁ (k := k) = rep₂ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [rep₁, rep₂, diagN, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Wedge2Formalization.N6OnePointLong.N,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]

def U
    (x₁ x₂ y₁ y₂ z₁ z₂ : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (upperToeplitz (k := k) 1 x₁ x₂)
    (Wedge2Formalization.N6OnePointLong.upperHankel (k := k) y₂ y₁ 0)
    (Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) 0 z₁ z₂)
    (lowerToeplitz (k := k) 1 (-x₁) (x₁ * x₁ - x₂ + y₁ * z₁))

def Levi
    (a b c d : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (a • (1 : Matrix I I k))
    (b • flip (k := k))
    (c • flip (k := k))
    (d • (1 : Matrix I I k))

/-- The exact pointwise kernel family on the indecomposable `3[a]` row.  This is the
explicit `SL_2(k[t]/(t^3))`-type block family written in coordinates. -/
def K :
    Set (Matrix V V k) :=
  { g | ∃ a₀ a₁ a₂ b₀ b₁ b₂ c₀ c₁ c₂ d₀ d₁ d₂ : k,
      a₀ * d₀ - b₀ * c₀ = 1 ∧
      a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁ = 0 ∧
      a₀ * d₂ + a₁ * d₁ + a₂ * d₀ - b₂ * c₀ - b₁ * c₁ - b₀ * c₂ = 0 ∧
      g =
        Matrix.fromBlocks
          (upperToeplitz (k := k) a₀ a₁ a₂)
          (Wedge2Formalization.N6OnePointLong.upperHankel (k := k) b₂ b₁ b₀)
          (Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) c₀ c₁ c₂)
          (lowerToeplitz (k := k) d₀ d₁ d₂) }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

def lift
    (a b : k) :
    Matrix V V k :=
  let gS : Matrix V V k :=
    Matrix.fromBlocks
      (Wedge2Formalization.N6OnePointLong.scaleTop (k := k) a)
      0
      0
      (Wedge2Formalization.N6OnePointLong.scaleBottom (k := k) a)
  let gU : Matrix V V k :=
    Matrix.fromBlocks
      (Wedge2Formalization.N6OnePointLong.shearTop (k := k) b)
      0
      0
      (Wedge2Formalization.N6OnePointLong.shearBottom (k := k) b)
  gS * gU

theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • rep₁ (k := k) + b • rep₂ (k := k)) = 0 ↔
      a = 0 :=
  Wedge2Formalization.N6Summary.onePointLong_det_zero_iff
    (k := k) (a := a) (b := b)

theorem mem_K_shape_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ a₀ a₁ a₂ b₀ b₁ b₂ c₀ c₁ c₂ d₀ d₁ d₂ : k,
        a₀ * d₀ - b₀ * c₀ = 1 ∧
        a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁ = 0 ∧
        a₀ * d₂ + a₁ * d₁ + a₂ * d₀ - b₂ * c₀ - b₁ * c₁ - b₀ * c₂ = 0 ∧
        g =
          Matrix.fromBlocks
            (upperToeplitz (k := k) a₀ a₁ a₂)
            (Wedge2Formalization.N6OnePointLong.upperHankel (k := k) b₂ b₁ b₀)
            (Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) c₀ c₁ c₂)
            (lowerToeplitz (k := k) d₀ d₁ d₂) := by
  rfl

private theorem U_mul_Levi_eq_shape
    (x₁ x₂ y₁ y₂ z₁ z₂ a b c d : k) :
    U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ * Levi (k := k) a b c d =
      Matrix.fromBlocks
        (upperToeplitz (k := k) a (a * x₁ + c * y₁) (a * x₂ + c * y₂))
        (Wedge2Formalization.N6OnePointLong.upperHankel
          (k := k) (b * x₂ + d * y₂) (b * x₁ + d * y₁) b)
        (Wedge2Formalization.N6OnePointLong.lowerHankel
          (k := k) c (-c * x₁ + a * z₁) ((x₁ * x₁ - x₂ + y₁ * z₁) * c + a * z₂))
        (lowerToeplitz (k := k) d (-d * x₁ + b * z₁)
          ((x₁ * x₁ - x₂ + y₁ * z₁) * d + b * z₂)) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [U, Levi, upperToeplitz, lowerToeplitz, flip,
      Wedge2Formalization.N6OnePointLong.upperHankel,
      Wedge2Formalization.N6OnePointLong.lowerHankel,
      Matrix.fromBlocks_multiply, Matrix.mul_apply, Fin.sum_univ_three] <;>
    ring

private theorem shape_eq_U_mul_Levi
    (a₀ a₁ a₂ b₀ b₁ b₂ c₀ c₁ c₂ d₀ d₁ d₂ : k)
    (h₀ : a₀ * d₀ - b₀ * c₀ = 1)
    (h₁ : a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁ = 0)
    (h₂ : a₀ * d₂ + a₁ * d₁ + a₂ * d₀ - b₂ * c₀ - b₁ * c₁ - b₀ * c₂ = 0) :
    Matrix.fromBlocks
      (upperToeplitz (k := k) a₀ a₁ a₂)
      (Wedge2Formalization.N6OnePointLong.upperHankel (k := k) b₂ b₁ b₀)
      (Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) c₀ c₁ c₂)
      (lowerToeplitz (k := k) d₀ d₁ d₂) =
      U (k := k)
        (a₁ * d₀ - b₁ * c₀)
        (a₂ * d₀ - b₂ * c₀)
        (a₀ * b₁ - a₁ * b₀)
        (a₀ * b₂ - a₂ * b₀)
        (c₁ * d₀ - c₀ * d₁)
        (c₂ * d₀ - c₀ * d₂) *
      Levi (k := k) a₀ b₀ c₀ d₀ := by
  have ha1 :
      a₀ * (a₁ * d₀ - b₁ * c₀) + c₀ * (a₀ * b₁ - a₁ * b₀) = a₁ := by
    calc
      a₀ * (a₁ * d₀ - b₁ * c₀) + c₀ * (a₀ * b₁ - a₁ * b₀)
          = a₁ * (a₀ * d₀ - b₀ * c₀) := by ring
      _ = a₁ := by rw [h₀]; ring
  have ha2 :
      a₀ * (a₂ * d₀ - b₂ * c₀) + c₀ * (a₀ * b₂ - a₂ * b₀) = a₂ := by
    calc
      a₀ * (a₂ * d₀ - b₂ * c₀) + c₀ * (a₀ * b₂ - a₂ * b₀)
          = a₂ * (a₀ * d₀ - b₀ * c₀) := by ring
      _ = a₂ := by rw [h₀]; ring
  have hb1 :
      b₀ * (a₁ * d₀ - b₁ * c₀) + d₀ * (a₀ * b₁ - a₁ * b₀) = b₁ := by
    calc
      b₀ * (a₁ * d₀ - b₁ * c₀) + d₀ * (a₀ * b₁ - a₁ * b₀)
          = b₁ * (a₀ * d₀ - b₀ * c₀) := by ring
      _ = b₁ := by rw [h₀]; ring
  have hb2 :
      b₀ * (a₂ * d₀ - b₂ * c₀) + d₀ * (a₀ * b₂ - a₂ * b₀) = b₂ := by
    calc
      b₀ * (a₂ * d₀ - b₂ * c₀) + d₀ * (a₀ * b₂ - a₂ * b₀)
          = b₂ * (a₀ * d₀ - b₀ * c₀) := by ring
      _ = b₂ := by rw [h₀]; ring
  have hc1 :
      -c₀ * (a₁ * d₀ - b₁ * c₀) + a₀ * (c₁ * d₀ - c₀ * d₁) = c₁ := by
    have hzero :
        -c₀ * (a₁ * d₀ - b₁ * c₀) + a₀ * (c₁ * d₀ - c₀ * d₁) - c₁ =
          c₁ * (a₀ * d₀ - b₀ * c₀ - 1) - c₀ * (a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁) := by
      ring
    have hzero' :
        -c₀ * (a₁ * d₀ - b₁ * c₀) + a₀ * (c₁ * d₀ - c₀ * d₁) - c₁ = 0 := by
      rw [hzero, h₀, h₁]
      ring
    exact sub_eq_zero.mp hzero'
  have hd1 :
      -d₀ * (a₁ * d₀ - b₁ * c₀) + b₀ * (c₁ * d₀ - c₀ * d₁) = d₁ := by
    have hzero :
        -d₀ * (a₁ * d₀ - b₁ * c₀) + b₀ * (c₁ * d₀ - c₀ * d₁) - d₁ =
          d₁ * (a₀ * d₀ - b₀ * c₀ - 1) - d₀ * (a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁) := by
      ring
    have hzero' :
        -d₀ * (a₁ * d₀ - b₁ * c₀) + b₀ * (c₁ * d₀ - c₀ * d₁) - d₁ = 0 := by
      rw [hzero, h₀, h₁]
      ring
    exact sub_eq_zero.mp hzero'
  have hc2 :
      ((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
          (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * c₀ +
        a₀ * (c₂ * d₀ - c₀ * d₂) = c₂ := by
    have hzero :
        (((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
              (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * c₀ +
            a₀ * (c₂ * d₀ - c₀ * d₂) - c₂) =
          (c₂ - c₀ * (a₁ * d₁ - b₁ * c₁)) * (a₀ * d₀ - b₀ * c₀ - 1) +
            c₀ * (a₁ * d₀ - b₁ * c₀) * (a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁) -
            c₀ * (a₀ * d₂ + a₁ * d₁ + a₂ * d₀ - b₂ * c₀ - b₁ * c₁ - b₀ * c₂) := by
      ring
    have hzero' :
        ((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
            (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * c₀ +
          a₀ * (c₂ * d₀ - c₀ * d₂) - c₂ = 0 := by
      rw [hzero, h₀, h₁, h₂]
      ring
    exact sub_eq_zero.mp hzero'
  have hd2 :
      ((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
          (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * d₀ +
        b₀ * (c₂ * d₀ - c₀ * d₂) = d₂ := by
    have hzero :
        (((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
              (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * d₀ +
            b₀ * (c₂ * d₀ - c₀ * d₂) - d₂) =
          (d₂ - d₀ * (a₁ * d₁ - b₁ * c₁)) * (a₀ * d₀ - b₀ * c₀ - 1) +
            d₀ * (a₁ * d₀ - b₁ * c₀) * (a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁) -
            d₀ * (a₀ * d₂ + a₁ * d₁ + a₂ * d₀ - b₂ * c₀ - b₁ * c₁ - b₀ * c₂) := by
      ring
    have hzero' :
        ((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
            (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * d₀ +
          b₀ * (c₂ * d₀ - c₀ * d₂) - d₂ = 0 := by
      rw [hzero, h₀, h₁, h₂]
      ring
    exact sub_eq_zero.mp hzero'
  apply Eq.symm
  calc
    U (k := k)
        (a₁ * d₀ - b₁ * c₀)
        (a₂ * d₀ - b₂ * c₀)
        (a₀ * b₁ - a₁ * b₀)
        (a₀ * b₂ - a₂ * b₀)
        (c₁ * d₀ - c₀ * d₁)
        (c₂ * d₀ - c₀ * d₂) *
      Levi (k := k) a₀ b₀ c₀ d₀
        =
      Matrix.fromBlocks
        (upperToeplitz (k := k) a₀
          (a₀ * (a₁ * d₀ - b₁ * c₀) + c₀ * (a₀ * b₁ - a₁ * b₀))
          (a₀ * (a₂ * d₀ - b₂ * c₀) + c₀ * (a₀ * b₂ - a₂ * b₀)))
        (Wedge2Formalization.N6OnePointLong.upperHankel
          (k := k)
          (b₀ * (a₂ * d₀ - b₂ * c₀) + d₀ * (a₀ * b₂ - a₂ * b₀))
          (b₀ * (a₁ * d₀ - b₁ * c₀) + d₀ * (a₀ * b₁ - a₁ * b₀))
          b₀)
        (Wedge2Formalization.N6OnePointLong.lowerHankel
          (k := k)
          c₀
          (-c₀ * (a₁ * d₀ - b₁ * c₀) + a₀ * (c₁ * d₀ - c₀ * d₁))
          (((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
              (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * c₀ +
            a₀ * (c₂ * d₀ - c₀ * d₂)))
        (lowerToeplitz (k := k)
          d₀
          (-d₀ * (a₁ * d₀ - b₁ * c₀) + b₀ * (c₁ * d₀ - c₀ * d₁))
          (((a₁ * d₀ - b₁ * c₀) * (a₁ * d₀ - b₁ * c₀) - (a₂ * d₀ - b₂ * c₀) +
              (a₀ * b₁ - a₁ * b₀) * (c₁ * d₀ - c₀ * d₁)) * d₀ +
            b₀ * (c₂ * d₀ - c₀ * d₂))) := by
        simpa using U_mul_Levi_eq_shape (k := k)
          (a := a₀) (b := b₀) (c := c₀) (d := d₀)
          (x₁ := a₁ * d₀ - b₁ * c₀)
          (x₂ := a₂ * d₀ - b₂ * c₀)
          (y₁ := a₀ * b₁ - a₁ * b₀)
          (y₂ := a₀ * b₂ - a₂ * b₀)
          (z₁ := c₁ * d₀ - c₀ * d₁)
          (z₂ := c₂ * d₀ - c₀ * d₂)
    _ =
      Matrix.fromBlocks
        (upperToeplitz (k := k) a₀ a₁ a₂)
        (Wedge2Formalization.N6OnePointLong.upperHankel (k := k) b₂ b₁ b₀)
        (Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) c₀ c₁ c₂)
        (lowerToeplitz (k := k) d₀ d₁ d₂) := by
        rw [ha1, ha2, hb1, hb2, hc1, hc2, hd1, hd2]

set_option maxHeartbeats 20000000 in
private theorem raw_family_pointwise
    (a₀ a₁ a₂ b₀ b₁ b₂ c₀ c₁ c₂ d₀ d₁ d₂ : k)
    (h₀ : a₀ * d₀ - b₀ * c₀ = 1)
    (h₁ : a₀ * d₁ + a₁ * d₀ - b₁ * c₀ - b₀ * c₁ = 0)
    (h₂ : a₀ * d₂ + a₁ * d₁ + a₂ * d₀ - b₂ * c₀ - b₁ * c₁ - b₀ * c₂ = 0) :
    Wedge2Formalization.N6OnePointLong.FixesPairBivector
      (Matrix.fromBlocks
        (upperToeplitz (k := k) a₀ a₁ a₂)
        (Wedge2Formalization.N6OnePointLong.upperHankel (k := k) b₂ b₁ b₀)
        (Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) c₀ c₁ c₂)
        (lowerToeplitz (k := k) d₀ d₁ d₂)) := by
  let A : Matrix I I k := upperToeplitz (k := k) a₀ a₁ a₂
  let B : Matrix I I k := Wedge2Formalization.N6OnePointLong.upperHankel (k := k) b₂ b₁ b₀
  let C : Matrix I I k := Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) c₀ c₁ c₂
  let D : Matrix I I k := lowerToeplitz (k := k) d₀ d₁ d₂
  have hrep1 :
      Wedge2Formalization.N6OnePointLong.FixesBivector
        (rep₁ (k := k))
        (Matrix.fromBlocks A B C D) := by
    rw [Wedge2Formalization.N6OnePointLong.FixesBivector]
    have hTL : -B * (1 : Matrix I I k)ᵀ * Aᵀ + A * (1 : Matrix I I k) * Bᵀ = 0 := by
      rw [Wedge2Formalization.N6OnePointLong.upperHankel_transpose, upperToeplitz_transpose]
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [A, B, upperToeplitz, lowerToeplitz,
          Wedge2Formalization.N6OnePointLong.upperHankel,
          Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] <;>
        ring_nf
    have hTR : -B * (1 : Matrix I I k)ᵀ * Cᵀ + A * (1 : Matrix I I k) * Dᵀ = (1 : Matrix I I k) := by
      rw [Wedge2Formalization.N6OnePointLong.lowerHankel_transpose, lowerToeplitz_transpose]
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [A, B, C, D, upperToeplitz, lowerToeplitz,
          Wedge2Formalization.N6OnePointLong.upperHankel,
          Wedge2Formalization.N6OnePointLong.lowerHankel,
          Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three,
          h₀, h₁, h₂] <;>
        ring_nf at h₀ h₁ h₂ ⊢ <;>
        first | assumption | linarith [h₀, h₁, h₂] | norm_num
    have hBL0 :
        -(C * (1 : Matrix I I k)ᵀ * Bᵀ) + D * (1 : Matrix I I k) * Aᵀ = (1 : Matrix I I k) := by
      simpa [Matrix.transpose_add, Matrix.transpose_mul, Matrix.transpose_neg,
        Wedge2Formalization.N6OnePointLong.upperHankel_transpose,
        Wedge2Formalization.N6OnePointLong.lowerHankel_transpose,
        upperToeplitz_transpose, lowerToeplitz_transpose] using congrArg Matrix.transpose hTR
    have hBL : -D * (1 : Matrix I I k)ᵀ * Aᵀ + C * (1 : Matrix I I k) * Bᵀ = -(1 : Matrix I I k) := by
      have hBL0' : -(C * (1 : Matrix I I k) * Bᵀ) + D * (1 : Matrix I I k) * Aᵀ = 1 := by
        simpa [Matrix.transpose_one] using hBL0
      rw [Matrix.transpose_one]
      rw [eq_neg_iff_add_eq_zero]
      have hsum :
          (C * (1 : Matrix I I k) * Bᵀ + -(D * (1 : Matrix I I k) * Aᵀ)) +
            (-(C * (1 : Matrix I I k) * Bᵀ) + D * (1 : Matrix I I k) * Aᵀ) = 0 := by
        abel_nf
      rw [hBL0'] at hsum
      simpa [add_comm, add_left_comm, add_assoc] using hsum
    have hBR : -D * (1 : Matrix I I k)ᵀ * Cᵀ + C * (1 : Matrix I I k) * Dᵀ = 0 := by
      rw [Wedge2Formalization.N6OnePointLong.lowerHankel_transpose, lowerToeplitz_transpose]
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [A, C, D, upperToeplitz, lowerToeplitz,
          Wedge2Formalization.N6OnePointLong.lowerHankel,
          Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] <;>
        ring_nf
    calc
      Wedge2Formalization.N6OnePointLong.ActBivector
          (rep₁ (k := k))
          (Matrix.fromBlocks A B C D) =
        Matrix.fromBlocks
          (-B * (1 : Matrix I I k)ᵀ * Aᵀ + A * (1 : Matrix I I k) * Bᵀ)
          (-B * (1 : Matrix I I k)ᵀ * Cᵀ + A * (1 : Matrix I I k) * Dᵀ)
          (-D * (1 : Matrix I I k)ᵀ * Aᵀ + C * (1 : Matrix I I k) * Bᵀ)
          (-D * (1 : Matrix I I k)ᵀ * Cᵀ + C * (1 : Matrix I I k) * Dᵀ) := by
            simpa [A, B, C, D, rep₁, Wedge2Formalization.N6OnePointLong.rep₁] using
              (Wedge2Formalization.N6OnePointLong.act_fromBlocks
                (k := k) (A := A) (B := B) (C := C) (D := D)
                (M := (1 : Matrix I I k)))
      _ = rep₁ (k := k) := by
            rw [hTL, hTR, hBL, hBR]
            simp [rep₁, Wedge2Formalization.N6OnePointLong.rep₁]
  have hAcomm :
      A * Wedge2Formalization.N6OnePointLong.N (k := k) =
        Wedge2Formalization.N6OnePointLong.N (k := k) * A := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [A, upperToeplitz, Wedge2Formalization.N6OnePointLong.N,
        Matrix.mul_apply, Fin.sum_univ_three]
  have hBcomm :
      Wedge2Formalization.N6OnePointLong.N (k := k) * B =
        B * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ := by
    simpa [B] using
      Wedge2Formalization.N6OnePointLong.upperHankel_comm
        (k := k) b₂ b₁ b₀
  have hCcomm :
      C * Wedge2Formalization.N6OnePointLong.N (k := k) =
        (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * C := by
    simpa [C] using
      Wedge2Formalization.N6OnePointLong.lowerHankel_comm
        (k := k) c₀ c₁ c₂
  have hDcomm :
      D * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ =
        (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * D := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [D, lowerToeplitz, Wedge2Formalization.N6OnePointLong.N,
        Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]
  have hcomm :
      Matrix.fromBlocks A B C D * diagN (k := k) =
        diagN (k := k) * Matrix.fromBlocks A B C D := by
    ext i j
    rcases i with i | i <;> rcases j with j | j
    · simpa [diagN, Matrix.fromBlocks_multiply] using congrArg (fun M => M i j) hAcomm
    · simpa [diagN, Matrix.fromBlocks_multiply] using congrArg (fun M => M i j) hBcomm.symm
    · simpa [diagN, Matrix.fromBlocks_multiply] using congrArg (fun M => M i j) hCcomm
    · simpa [diagN, Matrix.fromBlocks_multiply] using congrArg (fun M => M i j) hDcomm
  have hrep2 :
      Wedge2Formalization.N6OnePointLong.FixesBivector
        (rep₂ (k := k))
        (Matrix.fromBlocks A B C D) := by
    rw [Wedge2Formalization.N6OnePointLong.FixesBivector]
    calc
      Wedge2Formalization.N6OnePointLong.ActBivector
          (rep₂ (k := k))
          (Matrix.fromBlocks A B C D)
          = Matrix.fromBlocks A B C D * (diagN (k := k) * rep₁ (k := k)) *
              (Matrix.fromBlocks A B C D)ᵀ := by
              simp [Wedge2Formalization.N6OnePointLong.ActBivector, diagN_mul_rep₁]
      _ = diagN (k := k) *
            (Matrix.fromBlocks A B C D * rep₁ (k := k) *
              (Matrix.fromBlocks A B C D)ᵀ) := by
              simpa [Matrix.mul_assoc] using
                congrArg
                  (fun M => M * (rep₁ (k := k) * (Matrix.fromBlocks A B C D)ᵀ))
                  hcomm
      _ = diagN (k := k) * rep₁ (k := k) := by
              simpa [Wedge2Formalization.N6OnePointLong.FixesBivector] using congrArg (fun M => diagN (k := k) * M) hrep1
      _ = rep₂ (k := k) := diagN_mul_rep₁ (k := k)
  exact ⟨hrep1, hrep2⟩

private theorem upperToeplitz_zero_zero
    (a : k) :
    upperToeplitz (k := k) a 0 0 = a • (1 : Matrix I I k) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [upperToeplitz]

private theorem lowerToeplitz_zero_zero
    (d : k) :
    lowerToeplitz (k := k) d 0 0 = d • (1 : Matrix I I k) := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [lowerToeplitz]

theorem Levi_pointwise
    (a b c d : k)
    (hdet : a * d - b * c = 1) :
    Wedge2Formalization.N6OnePointLong.FixesPairBivector
      (Levi (k := k) a b c d) := by
  simpa [Levi, upperToeplitz_zero_zero, lowerToeplitz_zero_zero, flip,
    Wedge2Formalization.N6OnePointLong.upperHankel,
    Wedge2Formalization.N6OnePointLong.lowerHankel]
    using raw_family_pointwise (k := k)
      a 0 0 b 0 0 c 0 0 d 0 0
      hdet
      (by ring)
      (by ring)

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ a b c d x₁ x₂ y₁ y₂ z₁ z₂ : k,
        a * d - b * c = 1 ∧
        g = U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ * Levi (k := k) a b c d := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨a₀, a₁, a₂, b₀, b₁, b₂, c₀, c₁, c₂, d₀, d₁, d₂, h₀, h₁, h₂, rfl⟩
    refine ⟨a₀, b₀, c₀, d₀,
      a₁ * d₀ - b₁ * c₀,
      a₂ * d₀ - b₂ * c₀,
      a₀ * b₁ - a₁ * b₀,
      a₀ * b₂ - a₂ * b₀,
      c₁ * d₀ - c₀ * d₁,
      c₂ * d₀ - c₀ * d₂,
      h₀, ?_⟩
    simpa using shape_eq_U_mul_Levi (k := k)
      a₀ a₁ a₂ b₀ b₁ b₂ c₀ c₁ c₂ d₀ d₁ d₂ h₀ h₁ h₂
  · rintro ⟨a, b, c, d, x₁, x₂, y₁, y₂, z₁, z₂, hdet, rfl⟩
    refine (mem_K_shape_iff (k := k)
      (g := U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ * Levi (k := k) a b c d)).2 ?_
    refine ⟨a, a * x₁ + c * y₁, a * x₂ + c * y₂,
      b, b * x₁ + d * y₁, b * x₂ + d * y₂,
      c, -c * x₁ + a * z₁, (x₁ * x₁ - x₂ + y₁ * z₁) * c + a * z₂,
      d, -d * x₁ + b * z₁, (x₁ * x₁ - x₂ + y₁ * z₁) * d + b * z₂, ?_, ?_, ?_, ?_⟩
    · exact hdet
    · ring
    · ring
    · simpa using (U_mul_Levi_eq_shape (k := k)
        (x₁ := x₁) (x₂ := x₂) (y₁ := y₁) (y₂ := y₂) (z₁ := z₁) (z₂ := z₂)
        (a := a) (b := b) (c := c) (d := d))

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_6 \rtimes SL_2(k)`. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ a b c d x₁ x₂ y₁ y₂ z₁ z₂ : k,
        a * d - b * c = 1 ∧
        g = U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ * Levi (k := k) a b c d :=
  mem_K_iff (k := k) (g := g)

theorem pointwise_of_mem_K
    (g : Matrix V V k)
    (hg : g ∈ K (k := k)) :
    Wedge2Formalization.N6OnePointLong.FixesPairBivector g := by
  rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
    ⟨a₀, a₁, a₂, b₀, b₁, b₂, c₀, c₁, c₂, d₀, d₁, d₂, h₀, h₁, h₂, rfl⟩
  exact raw_family_pointwise (k := k) a₀ a₁ a₂ b₀ b₁ b₂ c₀ c₁ c₂ d₀ d₁ d₂ h₀ h₁ h₂

theorem U_pointwise
    (x₁ x₂ y₁ y₂ z₁ z₂ : k) :
    Wedge2Formalization.N6OnePointLong.FixesPairBivector
      (U (k := k) x₁ x₂ y₁ y₂ z₁ z₂) := by
  simpa [U, upperToeplitz, lowerToeplitz]
    using raw_family_pointwise (k := k)
      1 x₁ x₂ 0 y₁ y₂ 0 z₁ z₂ 1 (-x₁) (x₁ * x₁ - x₂ + y₁ * z₁)
      (by ring)
      (by ring)
      (by ring)

/-- Exact pointwise criterion on the intrinsic unipotent kernel cell before adjoining
the pointwise torus. -/
theorem cell_pointwise_iff
    (B C : Matrix Wedge2Formalization.N6OnePointLong.I
      Wedge2Formalization.N6OnePointLong.I k) :
    Wedge2Formalization.N6OnePointLong.FixesPairBivector
      (Matrix.fromBlocks
        ((1 : Matrix Wedge2Formalization.N6OnePointLong.I
            Wedge2Formalization.N6OnePointLong.I k) + B * C)
        B
        C
        (1 : Matrix Wedge2Formalization.N6OnePointLong.I
          Wedge2Formalization.N6OnePointLong.I k)) ↔
      B =
        Wedge2Formalization.N6OnePointLong.upperHankel
          (k := k) (B 0 0) (B 0 1) (B 0 2) ∧
        C =
          Wedge2Formalization.N6OnePointLong.lowerHankel
            (k := k) (C 0 2) (C 1 2) (C 2 2) :=
  Wedge2Formalization.N6Summary.onePointLong_productCell_iff
    (k := k) (B := B) (C := C)

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6OnePointLong.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (!![a, b; 0, 1]) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.onePointLong_borel_lift_action
        (k := k) (a := a) (b := b) ha).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.onePointLong_borel_lift_action
        (k := k) (a := a) (b := b) ha).2

private theorem rep_pair_independent
    {a b : k}
    (h :
      a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h03 := congrArg
    (fun M => M (Sum.inl 0) (Sum.inr 0)) h
  have ha : a = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Wedge2Formalization.N6OnePointLong.N,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h03
  have h04 := congrArg
    (fun M => M (Sum.inl 0) (Sum.inr 1)) h
  have hb : b = 0 := by
    simpa [ha, rep₁, rep₂, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Wedge2Formalization.N6OnePointLong.N,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h04
  exact ⟨ha, hb⟩

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  exact one_ne_zero (rep_pair_independent (k := k) hcomb).1

private theorem rep₁_det_ne_zero : Matrix.det (rep₁ (k := k)) ≠ 0 := by
  intro hdet
  have h' : (1 : k) = 0 :=
    (det_zero_iff (k := k) (a := (1 : k)) (b := (0 : k))).1 (by simpa using hdet)
  exact one_ne_zero h'

private theorem rep₂_det_zero : Matrix.det (rep₂ (k := k)) = 0 := by
  have h :
      Matrix.det ((0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k)) = 0 :=
    (det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 rfl
  simpa using h

theorem det_ne_zero_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6OnePointLong.V
        Wedge2Formalization.N6OnePointLong.V k)
    (hfix :
      Wedge2Formalization.N6OnePointLong.FixesPairBivector g) :
    Matrix.det g ≠ 0 := by
  intro hg0
  have hdet := congrArg Matrix.det hfix.1
  rw [Wedge2Formalization.N6OnePointLong.ActBivector, Matrix.det_mul, Matrix.det_mul,
    Matrix.det_transpose] at hdet
  have hrep0 : Matrix.det (rep₁ (k := k)) = 0 := by
    simpa [hg0] using hdet.symm
  exact rep₁_det_ne_zero (k := k) hrep0

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6OnePointLong.V
        Wedge2Formalization.N6OnePointLong.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N6OnePointLong.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6OnePointLong.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h2' :
      Wedge2Formalization.N6OnePointLong.ActBivector
          (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    exact hM.2
  have h20 :
      M 1 0 = 0 := by
    have hrep₂det :
        Matrix.det
          (Wedge2Formalization.N6OnePointLong.ActBivector
            (rep₂ (k := k)) g) = 0 := by
      rw [Wedge2Formalization.N6OnePointLong.ActBivector, Matrix.det_mul, Matrix.det_mul,
        Matrix.det_transpose, rep₂_det_zero (k := k)]
      simp
    exact (det_zero_iff (k := k) (a := M 1 0) (b := M 1 1)).1 (h2'.symm ▸ hrep₂det)
  refine ⟨M, ?_, hM⟩
  constructor
  · simpa [h20]
  ·
    change Matrix.det M ≠ 0
    by_contra hdet
    have hdet_formula : M 0 0 * M 1 1 = 0 := by
      simpa [Matrix.det_fin_two, h20] using hdet
    rcases mul_eq_zero.mp hdet_formula with h00 | h11
    · have hsing :
          Matrix.det
            (Wedge2Formalization.N6OnePointLong.ActBivector
              (rep₁ (k := k)) g) = 0 := by
        rcases hM with ⟨h1, _⟩
        rw [h1]
        exact (det_zero_iff (k := k) (a := M 0 0) (b := M 0 1)).2 h00
      have hnonsing :
          Matrix.det
            (Wedge2Formalization.N6OnePointLong.ActBivector
              (rep₁ (k := k)) g) ≠ 0 := by
        rw [Wedge2Formalization.N6OnePointLong.ActBivector, Matrix.det_mul, Matrix.det_mul,
          Matrix.det_transpose]
        simp [rep₁_det_ne_zero (k := k), hg]
      exact hnonsing hsing
    · have hzero :
          Wedge2Formalization.N6OnePointLong.ActBivector
              (rep₂ (k := k)) g = 0 := by
        rcases hM with ⟨_, h2⟩
        simpa [ActsOnOrderedPair, h20, h11] using h2
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
      have hrep₂nonzero : rep₂ (k := k) ≠ 0 := by
        intro hzero'
        have hcomb : (0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k) = 0 := by
          simpa using hzero'
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).2
      exact hrep₂nonzero horig

private theorem commutes_diagN_of_pointwise
    (g : Matrix V V k)
    (hfix : Wedge2Formalization.N6OnePointLong.FixesPairBivector g) :
    g * diagN (k := k) = diagN (k := k) * g := by
  have hg : Matrix.det g ≠ 0 := det_ne_zero_of_pointwise (k := k) g hfix
  letI : Invertible (Matrix.det g) := invertibleOfNonzero hg
  letI : Invertible g := Matrix.invertibleOfDetInvertible g
  letI : Invertible (Matrix.det gᵀ) := by
    simpa [Matrix.det_transpose] using (inferInstance : Invertible (Matrix.det g))
  letI : Invertible gᵀ := Matrix.invertibleOfDetInvertible gᵀ
  have hunit : IsUnit (Matrix.det g) := isUnit_iff_ne_zero.mpr hg
  have hunitT : IsUnit (Matrix.det gᵀ) := by
    simpa [Matrix.det_transpose] using hunit
  have h1act :
      g * rep₁ (k := k) * gᵀ = rep₁ (k := k) := by
    simpa [Wedge2Formalization.N6OnePointLong.ActBivector] using hfix.1
  have h1inv0 := congrArg
    (fun M =>
      (-(rep₁ (k := k))) * M * (gᵀ)⁻¹) h1act
  have h1inv :
      (gᵀ)⁻¹ =
        (-(rep₁ (k := k))) * g * rep₁ (k := k) := by
    simpa [Wedge2Formalization.N6OnePointLong.ActBivector, Matrix.mul_assoc,
      neg_rep₁_mul_rep₁ (k := k), Matrix.mul_nonsing_inv _ hunitT] using h1inv0.symm
  have h2act :
      g * rep₂ (k := k) * gᵀ = rep₂ (k := k) := by
    simpa [Wedge2Formalization.N6OnePointLong.ActBivector] using hfix.2
  have h2eq0 := congrArg (fun M => M * (gᵀ)⁻¹) h2act
  have h2eq :
      g * rep₂ (k := k) = rep₂ (k := k) * (gᵀ)⁻¹ := by
    simpa [Wedge2Formalization.N6OnePointLong.ActBivector, Matrix.mul_assoc,
      Matrix.mul_nonsing_inv _ hunitT] using h2eq0
  calc
    g * diagN (k := k) = g * (rep₂ (k := k) * (-(rep₁ (k := k)))) := by
      rw [rep₂_mul_neg_rep₁ (k := k)]
    _ = (g * rep₂ (k := k)) * (-(rep₁ (k := k))) := by
      rw [Matrix.mul_assoc]
    _ = (rep₂ (k := k) * (gᵀ)⁻¹) * (-(rep₁ (k := k))) := by
      rw [h2eq]
    _ = (rep₂ (k := k) * (((-(rep₁ (k := k))) * g * rep₁ (k := k)))) * (-(rep₁ (k := k))) := by
      rw [h1inv]
    _ = rep₂ (k := k) * (((-(rep₁ (k := k))) * g) * (rep₁ (k := k) * (-(rep₁ (k := k))))) := by
      simp [Matrix.mul_assoc]
    _ = rep₂ (k := k) * ((-(rep₁ (k := k))) * g) := by
      simp [Matrix.mul_assoc, rep₁_mul_neg_rep₁ (k := k)]
    _ = rep₂ (k := k) * (-(rep₁ (k := k))) * g := by
      simp [Matrix.mul_assoc]
    _ = diagN (k := k) * g := by
      rw [rep₂_mul_neg_rep₁ (k := k)]

private theorem eq_upperToeplitz_of_comm
    (A : Matrix I I k)
    (h : A * Wedge2Formalization.N6OnePointLong.N (k := k) =
      Wedge2Formalization.N6OnePointLong.N (k := k) * A) :
    A = upperToeplitz (k := k) (A 0 0) (A 0 1) (A 0 2) := by
  have h10 : A 1 0 = 0 := by
    have h' := congrArg (fun M => M 0 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Fin.sum_univ_three] using h'.symm
  have h20 : A 2 0 = 0 := by
    have h' := congrArg (fun M => M 1 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Fin.sum_univ_three] using h'.symm
  have h21 : A 2 1 = 0 := by
    have h' := congrArg (fun M => M 2 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Fin.sum_univ_three] using h'
  have h11 : A 1 1 = A 0 0 := by
    have h' := congrArg (fun M => M 0 1) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Fin.sum_univ_three] using h'.symm
  have h12 : A 1 2 = A 0 1 := by
    have h' := congrArg (fun M => M 0 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Fin.sum_univ_three] using h'.symm
  have h22 : A 2 2 = A 0 0 := by
    have h' := congrArg (fun M => M 1 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Fin.sum_univ_three, h11] using h'.symm
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [upperToeplitz, h10, h20, h21, h11, h12, h22]

private theorem eq_upperHankel_of_comm
    (B : Matrix I I k)
    (h : Wedge2Formalization.N6OnePointLong.N (k := k) * B =
      B * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ) :
    B = Wedge2Formalization.N6OnePointLong.upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) := by
  have h10 : B 1 0 = B 0 1 := by
    have h' := congrArg (fun M => M 0 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h11 : B 1 1 = B 0 2 := by
    have h' := congrArg (fun M => M 0 1) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h12 : B 1 2 = 0 := by
    have h' := congrArg (fun M => M 0 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h20 : B 2 0 = B 0 2 := by
    have h' := congrArg (fun M => M 1 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three, h11] using h'
  have h21 : B 2 1 = 0 := by
    have h' := congrArg (fun M => M 1 1) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three, h12] using h'
  have h22 : B 2 2 = 0 := by
    have h' := congrArg (fun M => M 1 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Wedge2Formalization.N6OnePointLong.upperHankel, h10, h11, h12, h20, h21, h22]

private theorem eq_lowerHankel_of_comm
    (C : Matrix I I k)
    (h : C * Wedge2Formalization.N6OnePointLong.N (k := k) =
      (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * C) :
    C = Wedge2Formalization.N6OnePointLong.lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) := by
  have h00 : C 0 0 = 0 := by
    have h' := congrArg (fun M => M 0 1) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h01 : C 0 1 = 0 := by
    have h' := congrArg (fun M => M 0 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h10 : C 1 0 = 0 := by
    have h' := congrArg (fun M => M 2 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'.symm
  have h11 : C 1 1 = C 0 2 := by
    have h' := congrArg (fun M => M 1 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h20 : C 2 0 = C 0 2 := by
    have h' := congrArg (fun M => M 2 1) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three, h11] using h'
  have h21 : C 2 1 = C 1 2 := by
    have h' := congrArg (fun M => M 2 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Wedge2Formalization.N6OnePointLong.lowerHankel, h00, h01, h10, h11, h20, h21]

private theorem eq_lowerToeplitz_of_comm
    (D : Matrix I I k)
    (h : D * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ =
      (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * D) :
    D = lowerToeplitz (k := k) (D 0 0) (D 1 0) (D 2 0) := by
  have h01 : D 0 1 = 0 := by
    have h' := congrArg (fun M => M 0 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h02 : D 0 2 = 0 := by
    have h' := congrArg (fun M => M 0 1) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h11 : D 1 1 = D 0 0 := by
    have h' := congrArg (fun M => M 1 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'
  have h12 : D 1 2 = 0 := by
    have h' := congrArg (fun M => M 2 2) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] using h'.symm
  have h21 : D 2 1 = D 1 0 := by
    have h' := congrArg (fun M => M 2 0) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three, h11] using h'
  have h22 : D 2 2 = D 0 0 := by
    have h' := congrArg (fun M => M 2 1) h
    simpa [Wedge2Formalization.N6OnePointLong.N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three, h11] using h'
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [lowerToeplitz, h01, h02, h11, h12, h21, h22]

set_option maxRecDepth 4096 in
theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6OnePointLong.FixesPairBivector)
      (K (k := k)) := by
  intro g
  constructor
  · intro hfix
    have hcomm := commutes_diagN_of_pointwise (k := k) g hfix
    have hcomm' :
        Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ * diagN (k := k) =
          diagN (k := k) * Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
      simpa [Matrix.fromBlocks_toBlocks] using hcomm
    have hAcomm : g.toBlocks₁₁ * Wedge2Formalization.N6OnePointLong.N (k := k) =
        Wedge2Formalization.N6OnePointLong.N (k := k) * g.toBlocks₁₁ := by
      have h' := congrArg Matrix.toBlocks₁₁ hcomm'
      simpa [diagN, Matrix.fromBlocks_multiply] using h'
    have hBcomm : Wedge2Formalization.N6OnePointLong.N (k := k) * g.toBlocks₁₂ =
        g.toBlocks₁₂ * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ := by
      have h' := congrArg Matrix.toBlocks₁₂ hcomm'
      simpa [diagN, Matrix.fromBlocks_multiply] using h'.symm
    have hCcomm : g.toBlocks₂₁ * Wedge2Formalization.N6OnePointLong.N (k := k) =
        (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * g.toBlocks₂₁ := by
      have h' := congrArg Matrix.toBlocks₂₁ hcomm'
      simpa [diagN, Matrix.fromBlocks_multiply] using h'
    have hDcomm : g.toBlocks₂₂ * (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ =
        (Wedge2Formalization.N6OnePointLong.N (k := k))ᵀ * g.toBlocks₂₂ := by
      have h' := congrArg Matrix.toBlocks₂₂ hcomm'
      simpa [diagN, Matrix.fromBlocks_multiply] using h'
    have hAshape := eq_upperToeplitz_of_comm (k := k) g.toBlocks₁₁ hAcomm
    have hBshape := eq_upperHankel_of_comm (k := k) g.toBlocks₁₂ hBcomm
    have hCshape := eq_lowerHankel_of_comm (k := k) g.toBlocks₂₁ hCcomm
    have hDshape := eq_lowerToeplitz_of_comm (k := k) g.toBlocks₂₂ hDcomm
    have hfix1' :
        Matrix.fromBlocks
            (-g.toBlocks₁₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₁₁ᵀ +
              g.toBlocks₁₁ * (1 : Matrix I I k) * g.toBlocks₁₂ᵀ)
            (-g.toBlocks₁₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₂₁ᵀ +
              g.toBlocks₁₁ * (1 : Matrix I I k) * g.toBlocks₂₂ᵀ)
            (-g.toBlocks₂₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₁₁ᵀ +
              g.toBlocks₂₁ * (1 : Matrix I I k) * g.toBlocks₁₂ᵀ)
            (-g.toBlocks₂₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₂₁ᵀ +
              g.toBlocks₂₁ * (1 : Matrix I I k) * g.toBlocks₂₂ᵀ)
          = rep₁ (k := k) := by
      have hfix1'' :
          Wedge2Formalization.N6OnePointLong.ActBivector
              (rep₁ (k := k))
              (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) =
            rep₁ (k := k) := by
        simpa [Matrix.fromBlocks_toBlocks] using hfix.1
      have hact1 :
          Wedge2Formalization.N6OnePointLong.ActBivector
              (rep₁ (k := k))
              (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) =
            Matrix.fromBlocks
              (-g.toBlocks₁₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₁₁ᵀ +
                g.toBlocks₁₁ * (1 : Matrix I I k) * g.toBlocks₁₂ᵀ)
              (-g.toBlocks₁₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₂₁ᵀ +
                g.toBlocks₁₁ * (1 : Matrix I I k) * g.toBlocks₂₂ᵀ)
              (-g.toBlocks₂₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₁₁ᵀ +
                g.toBlocks₂₁ * (1 : Matrix I I k) * g.toBlocks₁₂ᵀ)
              (-g.toBlocks₂₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₂₁ᵀ +
                g.toBlocks₂₁ * (1 : Matrix I I k) * g.toBlocks₂₂ᵀ) := by
        simpa [rep₁, Wedge2Formalization.N6OnePointLong.rep₁] using
          (Wedge2Formalization.N6OnePointLong.act_fromBlocks
            (k := k)
            (A := g.toBlocks₁₁) (B := g.toBlocks₁₂)
            (C := g.toBlocks₂₁) (D := g.toBlocks₂₂)
            (M := (1 : Matrix I I k)))
      rw [hact1] at hfix1''
      exact hfix1''
    have hTR : -g.toBlocks₁₂ * (1 : Matrix I I k)ᵀ * g.toBlocks₂₁ᵀ +
          g.toBlocks₁₁ * (1 : Matrix I I k) * g.toBlocks₂₂ᵀ = (1 : Matrix I I k) := by
      have h' := congrArg Matrix.toBlocks₁₂ hfix1'
      simpa [rep₁, Wedge2Formalization.N6OnePointLong.rep₁] using h'
    have h0eq := congrArg (fun M => M 0 0) hTR
    have h1eq := congrArg (fun M => M 0 1) hTR
    have h2eq := congrArg (fun M => M 0 2) hTR
    rw [hAshape, hBshape, hCshape, hDshape,
      lowerToeplitz_transpose,
      Wedge2Formalization.N6OnePointLong.lowerHankel_transpose] at h0eq h1eq h2eq
    ring_nf at h0eq h1eq h2eq
    have h0eq' :
        -(g.toBlocks₁₂ 0 2 * g.toBlocks₂₁ 0 2) + g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 = 1 := by
      simpa [upperToeplitz, Wedge2Formalization.N6OnePointLong.upperHankel,
        Wedge2Formalization.N6OnePointLong.lowerHankel, Matrix.mul_apply, Fin.sum_univ_three] using h0eq
    have h1eq' :
        -(g.toBlocks₁₂ 0 1 * g.toBlocks₂₁ 0 2) + -(g.toBlocks₁₂ 0 2 * g.toBlocks₂₁ 1 2) +
            (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 1 0 + g.toBlocks₁₁ 0 1 * g.toBlocks₂₂ 0 0) = 0 := by
      simpa [upperToeplitz, Wedge2Formalization.N6OnePointLong.upperHankel,
        Wedge2Formalization.N6OnePointLong.lowerHankel, Matrix.mul_apply, Fin.sum_univ_three] using h1eq
    have h2eq' :
        -(g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 2) +
            (-(g.toBlocks₁₂ 0 1 * g.toBlocks₂₁ 1 2) + -(g.toBlocks₁₂ 0 2 * g.toBlocks₂₁ 2 2)) +
          (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 2 0 +
            (g.toBlocks₁₁ 0 1 * g.toBlocks₂₂ 1 0 + g.toBlocks₁₁ 0 2 * g.toBlocks₂₂ 0 0)) = 0 := by
      simpa [upperToeplitz, Wedge2Formalization.N6OnePointLong.upperHankel,
        Wedge2Formalization.N6OnePointLong.lowerHankel, Matrix.mul_apply, Fin.sum_univ_three] using h2eq
    refine ⟨g.toBlocks₁₁ 0 0, g.toBlocks₁₁ 0 1, g.toBlocks₁₁ 0 2,
      g.toBlocks₁₂ 0 2, g.toBlocks₁₂ 0 1, g.toBlocks₁₂ 0 0,
      g.toBlocks₂₁ 0 2, g.toBlocks₂₁ 1 2, g.toBlocks₂₁ 2 2,
      g.toBlocks₂₂ 0 0, g.toBlocks₂₂ 1 0, g.toBlocks₂₂ 2 0, ?_, ?_, ?_, ?_⟩
    · simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using h0eq'
    · simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using h1eq'
    · simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using h2eq'
    · calc
        g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
              exact (Matrix.fromBlocks_toBlocks g).symm
        _ =
            Matrix.fromBlocks
              (upperToeplitz (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₁ 0 1) (g.toBlocks₁₁ 0 2))
              (Wedge2Formalization.N6OnePointLong.upperHankel
                (k := k) (g.toBlocks₁₂ 0 0) (g.toBlocks₁₂ 0 1) (g.toBlocks₁₂ 0 2))
              (Wedge2Formalization.N6OnePointLong.lowerHankel
                (k := k) (g.toBlocks₂₁ 0 2) (g.toBlocks₂₁ 1 2) (g.toBlocks₂₁ 2 2))
              (lowerToeplitz (k := k) (g.toBlocks₂₂ 0 0) (g.toBlocks₂₂ 1 0) (g.toBlocks₂₂ 2 0)) := by
              rw [hAshape, hBshape, hCshape, hDshape]
              simp [upperToeplitz, lowerToeplitz,
                Wedge2Formalization.N6OnePointLong.upperHankel,
                Wedge2Formalization.N6OnePointLong.lowerHankel]
  · exact pointwise_of_mem_K (k := k) g

end Row1

end N6
end Paper
end Wedge2Formalization
