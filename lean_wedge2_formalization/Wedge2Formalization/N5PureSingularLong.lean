import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Tactic

open Matrix

namespace Wedge2Formalization
namespace N5PureSingularLong

variable {k : Type*} [Field k]

/-- The `3`-dimensional block. -/
abbrev W := Fin 3

/-- The `2`-dimensional block. -/
abbrev I := Fin 2

/-- The ambient `5`-dimensional space. -/
abbrev V := W ⊕ I

/-- The natural action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- Exact stabilizer of a bivector. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  ActBivector Ω g = Ω

/-- The bivector action is linear in the tensor variable. -/
theorem actBivector_add
    (Ω₁ Ω₂ : Matrix V V k) (g : Matrix V V k) :
    ActBivector (Ω₁ + Ω₂) g = ActBivector Ω₁ g + ActBivector Ω₂ g := by
  simp [ActBivector, Matrix.mul_add, Matrix.add_mul]

/-- The bivector action commutes with scalar multiplication in the tensor variable. -/
theorem actBivector_smul
    (c : k) (Ω : Matrix V V k) (g : Matrix V V k) :
    ActBivector (c • Ω) g = c • ActBivector Ω g := by
  simp [ActBivector, Matrix.mul_smul, smul_mul_assoc]

/-- The bivector action is multiplicative in the acting matrix. -/
theorem actBivector_mul
    (Ω : Matrix V V k)
    (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  simp [ActBivector, Matrix.mul_assoc]

/-- An invertible matrix cannot send a bivector to zero under the natural action. -/
theorem actBivector_eq_zero_iff_of_det_ne_zero
    (Ω : Matrix V V k) (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0) :
    ActBivector Ω g = 0 ↔ Ω = 0 := by
  letI : Invertible (Matrix.det g) := invertibleOfNonzero hg
  letI : Invertible g := Matrix.invertibleOfDetInvertible g
  have hunit : IsUnit (Matrix.det g) := isUnit_of_invertible (Matrix.det g)
  constructor
  · intro hzero
    have h' := congrArg (fun M => g⁻¹ * M * (g⁻¹)ᵀ) hzero
    simp [ActBivector, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hunit,
      Matrix.mul_nonsing_inv _ hunit, Matrix.transpose_nonsing_inv] at h'
    exact h'
  · intro hΩ
    simp [hΩ, ActBivector]

/-- The two multiplication matrices on `Sym¹ -> Sym²`. -/
def mulX : Matrix W I k := !![(1 : k), 0; 0, 1; 0, 0]

def mulY : Matrix W I k := !![(0 : k), 0; 1, 0; 0, 1]

/-- The pure singular length-two representative. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 (mulX (k := k)) (-(mulX (k := k))ᵀ) 0

/-- The second basis vector of the pure singular length-two representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 (mulY (k := k)) (-(mulY (k := k))ᵀ) 0

/-- Exact stabilizer of the pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rep₁ g ∧ FixesBivector rep₂ g

/-- Setwise stabilizer of the `2`-space in the chosen basis. -/
def PreservesSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector rep₁ g = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ g = γ • rep₁ + δ • rep₂

/-- The lower-left Hankel block appearing in the pointwise unipotent radical. -/
def lowerHankel (u v w z : k) : Matrix I W k :=
  !![u, v, w;
     v, w, z]

/-- Acting by a block diagonal matrix keeps the representatives in block form. -/
theorem act_blockDiagonal
    (P : Matrix W W k)
    (Q : Matrix I I k)
    (M : Matrix W I k) :
    ActBivector
        (Matrix.fromBlocks 0 M (-Mᵀ) 0)
        (Matrix.fromBlocks P 0 0 Q) =
      Matrix.fromBlocks 0 (P * M * Qᵀ) (-(P * M * Qᵀ)ᵀ) 0 := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
      Matrix.neg_mul, Matrix.mul_neg, Matrix.transpose_neg, Matrix.mul_assoc]

/-- Acting by a lower-left unipotent preserves the off-diagonal block and adds the
expected skew correction in the lower-right block. -/
theorem act_lowerLeftUnipotent
    (C : Matrix I W k)
    (M : Matrix W I k) :
    ActBivector
        (Matrix.fromBlocks 0 M (-Mᵀ) 0)
        (Matrix.fromBlocks (1 : Matrix W W k) 0 C (1 : Matrix I I k)) =
      Matrix.fromBlocks 0 M (-Mᵀ) (C * M - Mᵀ * Cᵀ) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
      Matrix.neg_mul, Matrix.mul_neg, sub_eq_add_neg, Matrix.mul_assoc]
    <;> ring_nf

/-- A scalar torus fixes the pair pointwise. -/
theorem pointwise_scale_family
    (t : k)
    (ht : t ≠ 0) :
    FixesPairBivector
      (Matrix.fromBlocks
        ((t⁻¹) • (1 : Matrix W W k))
        0
        0
        (t • (1 : Matrix I I k))) := by
  constructor
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₁
          (Matrix.fromBlocks ((t⁻¹) • (1 : Matrix W W k)) 0 0 (t • (1 : Matrix I I k))) =
        Matrix.fromBlocks
          0
          (((t⁻¹) • (1 : Matrix W W k)) * mulX * ((t • (1 : Matrix I I k)))ᵀ)
          ( -((((t⁻¹) • (1 : Matrix W W k)) * mulX * ((t • (1 : Matrix I I k)))ᵀ))ᵀ)
          0 := by
            simpa [rep₁, mulX] using
              act_blockDiagonal
                (k := k)
                (((t⁻¹) • (1 : Matrix W W k)))
                ((t • (1 : Matrix I I k)))
                (mulX (k := k))
      _ = rep₁ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, mulX, ht, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₂
          (Matrix.fromBlocks ((t⁻¹) • (1 : Matrix W W k)) 0 0 (t • (1 : Matrix I I k))) =
        Matrix.fromBlocks
          0
          (((t⁻¹) • (1 : Matrix W W k)) * mulY * ((t • (1 : Matrix I I k)))ᵀ)
          ( -((((t⁻¹) • (1 : Matrix W W k)) * mulY * ((t • (1 : Matrix I I k)))ᵀ))ᵀ)
          0 := by
            simpa [rep₂, mulY] using
              act_blockDiagonal
                (k := k)
                (((t⁻¹) • (1 : Matrix W W k)))
                ((t • (1 : Matrix I I k)))
                (mulY (k := k))
      _ = rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₂, mulY, ht, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]

/-- The Hankel lower-left block gives symmetric products with `mulX`. -/
theorem lowerHankel_mulX_symmetric
    (u v w z : k) :
    lowerHankel (k := k) u v w z * mulX (k := k) =
      (lowerHankel (k := k) u v w z * mulX (k := k))ᵀ := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [lowerHankel, mulX, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three]

/-- The Hankel lower-left block gives symmetric products with `mulY`. -/
theorem lowerHankel_mulY_symmetric
    (u v w z : k) :
    lowerHankel (k := k) u v w z * mulY (k := k) =
      (lowerHankel (k := k) u v w z * mulY (k := k))ᵀ := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [lowerHankel, mulY, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three]

/-- The lower-left Hankel family fixes the pure singular length-two pair pointwise. -/
theorem lowerHankel_pointwise_family
    (u v w z : k) :
    FixesPairBivector
      (Matrix.fromBlocks
        (1 : Matrix W W k)
        0
        (lowerHankel (k := k) u v w z)
        (1 : Matrix I I k)) := by
  constructor <;> rw [FixesBivector]
  ·
    calc
      ActBivector rep₁
          (Matrix.fromBlocks
            (1 : Matrix W W k)
            0
            (lowerHankel (k := k) u v w z)
            (1 : Matrix I I k)) =
        Matrix.fromBlocks
          0
          (mulX (k := k))
          (-(mulX (k := k))ᵀ)
          (lowerHankel (k := k) u v w z * mulX (k := k) -
            (mulX (k := k))ᵀ * (lowerHankel (k := k) u v w z)ᵀ) := by
          simpa [rep₁, mulX] using
            act_lowerLeftUnipotent (k := k) (lowerHankel (k := k) u v w z) (mulX (k := k))
      _ = rep₁ := by
          rw [lowerHankel_mulX_symmetric (k := k) u v w z]
          simp [rep₁, mulX]
  ·
    calc
      ActBivector rep₂
          (Matrix.fromBlocks
            (1 : Matrix W W k)
            0
            (lowerHankel (k := k) u v w z)
            (1 : Matrix I I k)) =
        Matrix.fromBlocks
          0
          (mulY (k := k))
          (-(mulY (k := k))ᵀ)
          (lowerHankel (k := k) u v w z * mulY (k := k) -
            (mulY (k := k))ᵀ * (lowerHankel (k := k) u v w z)ᵀ) := by
          simpa [rep₂, mulY] using
            act_lowerLeftUnipotent (k := k) (lowerHankel (k := k) u v w z) (mulY (k := k))
      _ = rep₂ := by
          rw [lowerHankel_mulY_symmetric (k := k) u v w z]
          simp [rep₂, mulY]

/-- Combining the scalar torus and the lower-left Hankel family gives a concrete
`5`-parameter pointwise subgroup. -/
theorem scaleLowerHankel_product_pointwise
    (t u v w z : k)
    (ht : t ≠ 0) :
    let gScale : Matrix V V k :=
      Matrix.fromBlocks
        ((t⁻¹) • (1 : Matrix W W k))
        0
        0
        (t • (1 : Matrix I I k))
    let gH : Matrix V V k :=
      Matrix.fromBlocks
        (1 : Matrix W W k)
        0
        (lowerHankel (k := k) u v w z)
        (1 : Matrix I I k)
    FixesPairBivector (gH * gScale) := by
  intro gScale gH
  have hScale := pointwise_scale_family (k := k) t ht
  have hH := lowerHankel_pointwise_family (k := k) u v w z
  constructor <;> rw [FixesBivector]
  ·
    calc
      ActBivector rep₁ (gH * gScale) =
          ActBivector (ActBivector rep₁ gScale) gH := by
            rw [actBivector_mul]
      _ = ActBivector rep₁ gH := by rw [hScale.1]
      _ = rep₁ := by simpa [FixesBivector] using hH.1
  ·
    calc
      ActBivector rep₂ (gH * gScale) =
          ActBivector (ActBivector rep₂ gScale) gH := by
            rw [actBivector_mul]
      _ = ActBivector rep₂ gH := by rw [hScale.2]
      _ = rep₂ := by simpa [FixesBivector] using hH.2

/-- The explicit pointwise stabilizer shape for the pure singular length-two pair. -/
def pointwiseShape (a u v w z : k) : Matrix V V k :=
  Matrix.fromBlocks
    (a • (1 : Matrix W W k))
    0
    (lowerHankel (k := k) u v w z)
    (a⁻¹ • (1 : Matrix I I k))

/-- Every matrix of the expected pointwise shape fixes the pure singular length-two
pair. -/
theorem pointwiseShape_fixes
    (a u v w z : k)
    (ha : a ≠ 0) :
    FixesPairBivector (pointwiseShape (k := k) a u v w z) := by
  have h :=
    scaleLowerHankel_product_pointwise
      (k := k)
      (t := a⁻¹)
      (u := a⁻¹ * u)
      (v := a⁻¹ * v)
      (w := a⁻¹ * w)
      (z := a⁻¹ * z)
      (ht := inv_ne_zero ha)
  simpa [pointwiseShape, lowerHankel, Matrix.fromBlocks_multiply, ha, inv_inv] using h

-- When the upper-right block vanishes, pointwise fixing forces the full scalar-plus-
-- Hankel shape predicted by the singular local model.
set_option maxHeartbeats 1200000 in
theorem fixesPair_fromBlocks_zeroUpperRight_iff_shape
    (A : Matrix W W k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks A 0 C D) ↔
      A 0 0 ≠ 0 ∧
        Matrix.fromBlocks A 0 C D =
          pointwiseShape (k := k) (A 0 0) (C 0 0) (C 0 1) (C 0 2) (C 1 2) := by
  constructor
  · rintro ⟨h₁, h₂⟩
    have hur1 : A * mulX (k := k) * Dᵀ = mulX (k := k) := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₁
      simpa [FixesBivector, ActBivector, rep₁, mulX, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have hur2 : A * mulY (k := k) * Dᵀ = mulY (k := k) := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₂
      simpa [FixesBivector, ActBivector, rep₂, mulY, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have hlr1 : C * mulX (k := k) * Dᵀ - D * (mulX (k := k))ᵀ * Cᵀ = 0 := by
      ext i j
      fin_cases i <;> fin_cases j
      ·
        have hij := congrArg (fun M => M (Sum.inr 0) (Sum.inr 0)) h₁
        simpa [FixesBivector, ActBivector, rep₁, mulX, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
      ·
        have hij := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h₁
        simpa [FixesBivector, ActBivector, rep₁, mulX, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
      ·
        have hij := congrArg (fun M => M (Sum.inr 1) (Sum.inr 0)) h₁
        simpa [FixesBivector, ActBivector, rep₁, mulX, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
      ·
        have hij := congrArg (fun M => M (Sum.inr 1) (Sum.inr 1)) h₁
        simpa [FixesBivector, ActBivector, rep₁, mulX, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
    have hlr2 : C * mulY (k := k) * Dᵀ - D * (mulY (k := k))ᵀ * Cᵀ = 0 := by
      ext i j
      fin_cases i <;> fin_cases j
      ·
        have hij := congrArg (fun M => M (Sum.inr 0) (Sum.inr 0)) h₂
        simpa [FixesBivector, ActBivector, rep₂, mulY, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
      ·
        have hij := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h₂
        simpa [FixesBivector, ActBivector, rep₂, mulY, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
      ·
        have hij := congrArg (fun M => M (Sum.inr 1) (Sum.inr 0)) h₂
        simpa [FixesBivector, ActBivector, rep₂, mulY, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
      ·
        have hij := congrArg (fun M => M (Sum.inr 1) (Sum.inr 1)) h₂
        simpa [FixesBivector, ActBivector, rep₂, mulY, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm,
          add_assoc] using hij
    have h00 : A 0 0 * D 0 0 + A 0 1 * D 0 1 = 1 := by
      have hij := congrArg (fun M => M 0 0) hur1
      simpa [mulX, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h01 : A 0 0 * D 1 0 + A 0 1 * D 1 1 = 0 := by
      have hij := congrArg (fun M => M 0 1) hur1
      simpa [mulX, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h10 : A 1 0 * D 0 0 + A 1 1 * D 0 1 = 0 := by
      have hij := congrArg (fun M => M 1 0) hur1
      simpa [mulX, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h11 : A 1 0 * D 1 0 + A 1 1 * D 1 1 = 1 := by
      have hij := congrArg (fun M => M 1 1) hur1
      simpa [mulX, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h20 : A 2 0 * D 0 0 + A 2 1 * D 0 1 = 0 := by
      have hij := congrArg (fun M => M 2 0) hur1
      simpa [mulX, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h21 : A 2 0 * D 1 0 + A 2 1 * D 1 1 = 0 := by
      have hij := congrArg (fun M => M 2 1) hur1
      simpa [mulX, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h00' : A 0 1 * D 0 0 + A 0 2 * D 0 1 = 0 := by
      have hij := congrArg (fun M => M 0 0) hur2
      simpa [mulY, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h01' : A 0 1 * D 1 0 + A 0 2 * D 1 1 = 0 := by
      have hij := congrArg (fun M => M 0 1) hur2
      simpa [mulY, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h10' : A 1 1 * D 0 0 + A 1 2 * D 0 1 = 1 := by
      have hij := congrArg (fun M => M 1 0) hur2
      simpa [mulY, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h11' : A 1 1 * D 1 0 + A 1 2 * D 1 1 = 0 := by
      have hij := congrArg (fun M => M 1 1) hur2
      simpa [mulY, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h20' : A 2 1 * D 0 0 + A 2 2 * D 0 1 = 0 := by
      have hij := congrArg (fun M => M 2 0) hur2
      simpa [mulY, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have h21' : A 2 1 * D 1 0 + A 2 2 * D 1 1 = 1 := by
      have hij := congrArg (fun M => M 2 1) hur2
      simpa [mulY, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have hC01 : C 0 0 * D 1 0 + C 0 1 * D 1 1 - D 0 0 * C 1 0 - D 0 1 * C 1 1 = 0 := by
      have hij := congrArg (fun M => M 0 1) hlr1
      simpa [mulX, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
        Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hij
    have hC02 : C 0 1 * D 1 0 + C 0 2 * D 1 1 - D 0 0 * C 1 1 - D 0 1 * C 1 2 = 0 := by
      have hij := congrArg (fun M => M 0 1) hlr2
      simpa [mulY, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
        Fin.sum_univ_three, sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hij
    have hDetNonzero : D 0 0 * D 1 1 - D 0 1 * D 1 0 ≠ 0 := by
      intro hDet
      have hDetEq :
          (A 0 0 * A 1 1 - A 0 1 * A 1 0) * (D 0 0 * D 1 1 - D 0 1 * D 1 0) = 1 := by
        calc
          (A 0 0 * A 1 1 - A 0 1 * A 1 0) * (D 0 0 * D 1 1 - D 0 1 * D 1 0) =
              (A 0 0 * D 0 0 + A 0 1 * D 0 1) * (A 1 0 * D 1 0 + A 1 1 * D 1 1) -
              (A 0 0 * D 1 0 + A 0 1 * D 1 1) * (A 1 0 * D 0 0 + A 1 1 * D 0 1) := by
                ring
          _ = 1 := by rw [h00, h11, h01, h10]; ring
      rw [hDet, mul_zero] at hDetEq
      exact zero_ne_one hDetEq
    have hA01Mul : A 0 1 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) = 0 := by
      calc
        A 0 1 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) =
            (A 0 1 * D 0 0 + A 0 2 * D 0 1) * D 1 1 -
              (A 0 1 * D 1 0 + A 0 2 * D 1 1) * D 0 1 := by
                ring
        _ = 0 := by rw [h00', h01']; ring
    have hA01 : A 0 1 = 0 := by
      exact (mul_eq_zero.mp hA01Mul).resolve_right hDetNonzero
    have hA02Mul : A 0 2 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) = 0 := by
      calc
        A 0 2 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) =
            (A 0 1 * D 1 0 + A 0 2 * D 1 1) * D 0 0 -
              (A 0 1 * D 0 0 + A 0 2 * D 0 1) * D 1 0 := by
                ring
        _ = 0 := by rw [h00', h01']; ring
    have hA02 : A 0 2 = 0 := by
      exact (mul_eq_zero.mp hA02Mul).resolve_right hDetNonzero
    have hA20Mul : A 2 0 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) = 0 := by
      calc
        A 2 0 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) =
            (A 2 0 * D 0 0 + A 2 1 * D 0 1) * D 1 1 -
              (A 2 0 * D 1 0 + A 2 1 * D 1 1) * D 0 1 := by
                ring
        _ = 0 := by rw [h20, h21]; ring
    have hA20 : A 2 0 = 0 := by
      exact (mul_eq_zero.mp hA20Mul).resolve_right hDetNonzero
    have hA21Mul : A 2 1 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) = 0 := by
      calc
        A 2 1 * (D 0 0 * D 1 1 - D 0 1 * D 1 0) =
            (A 2 0 * D 1 0 + A 2 1 * D 1 1) * D 0 0 -
              (A 2 0 * D 0 0 + A 2 1 * D 0 1) * D 1 0 := by
                ring
        _ = 0 := by rw [h20, h21]; ring
    have hA21 : A 2 1 = 0 := by
      exact (mul_eq_zero.mp hA21Mul).resolve_right hDetNonzero
    have hA00D00 : A 0 0 * D 0 0 = 1 := by simpa [hA01] using h00
    have hA00 : A 0 0 ≠ 0 := by
      intro h0
      have : (0 : k) = 1 := by simpa [h0] using hA00D00
      exact zero_ne_one this
    have hD00 : D 0 0 ≠ 0 := by
      intro h0
      have : (0 : k) = 1 := by simpa [h0] using hA00D00
      exact zero_ne_one this
    have hD10 : D 1 0 = 0 := by
      have : A 0 0 * D 1 0 = 0 := by simpa [hA01] using h01
      exact (mul_eq_zero.mp this).resolve_left hA00
    have hDetSimple : D 0 0 * D 1 1 ≠ 0 := by
      simpa [hD10] using hDetNonzero
    have hD11 : D 1 1 ≠ 0 := by
      intro h0
      exact hDetSimple (by simp [h0])
    have hA12 : A 1 2 = 0 := by
      have : A 1 2 * D 1 1 = 0 := by simpa [hD10] using h11'
      exact (mul_eq_zero.mp this).resolve_right hD11
    have hA11D00 : A 1 1 * D 0 0 = 1 := by simpa [hA12] using h10'
    have hA11eq : A 1 1 = A 0 0 := by
      have hsub : (A 1 1 - A 0 0) * D 0 0 = 0 := by
        calc
          (A 1 1 - A 0 0) * D 0 0 = A 1 1 * D 0 0 - A 0 0 * D 0 0 := by ring
          _ = 0 := by rw [hA11D00, hA00D00]; ring
      exact sub_eq_zero.mp ((mul_eq_zero.mp hsub).resolve_right hD00)
    have hA22D11 : A 2 2 * D 1 1 = 1 := by simpa [hA21, hD10] using h21'
    have hD01 : D 0 1 = 0 := by
      have : A 2 2 * D 0 1 = 0 := by simpa [hA21] using h20'
      have hA22 : A 2 2 ≠ 0 := by
        intro h0
        have : (0 : k) = 1 := by simpa [h0] using hA22D11
        exact zero_ne_one this
      exact (mul_eq_zero.mp this).resolve_left hA22
    have hA10 : A 1 0 = 0 := by
      have : A 1 0 * D 0 0 = 0 := by simpa [hD01] using h10
      exact (mul_eq_zero.mp this).resolve_right hD00
    have hA11D11 : A 1 1 * D 1 1 = 1 := by simpa [hA10] using h11
    have hA22eq : A 2 2 = A 0 0 := by
      have hsub : (A 2 2 - A 1 1) * D 1 1 = 0 := by
        calc
          (A 2 2 - A 1 1) * D 1 1 = A 2 2 * D 1 1 - A 1 1 * D 1 1 := by ring
          _ = 0 := by rw [hA22D11, hA11D11]; ring
      exact (sub_eq_zero.mp ((mul_eq_zero.mp hsub).resolve_right hD11)).trans hA11eq
    have hD00eq : D 0 0 = D 1 1 := by
      have hsub : A 0 0 * (D 0 0 - D 1 1) = 0 := by
        calc
          A 0 0 * (D 0 0 - D 1 1) = A 0 0 * D 0 0 - A 0 0 * D 1 1 := by ring
          _ = 0 := by rw [hA00D00, ← hA11eq, hA11D11]; ring
      exact sub_eq_zero.mp ((mul_eq_zero.mp hsub).resolve_left hA00)
    have hD11inv : D 1 1 = (A 0 0)⁻¹ := by
      calc
        D 1 1 = 1 * D 1 1 := by ring
        _ = ((A 0 0)⁻¹ * A 0 0) * D 1 1 := by rw [inv_mul_cancel₀ hA00]
        _ = (A 0 0)⁻¹ * (A 0 0 * D 1 1) := by ring
        _ = (A 0 0)⁻¹ := by rw [← hA11eq, hA11D11]; simp
    have hC10eq : C 0 1 = C 1 0 := by
      have : (C 0 1 - C 1 0) * D 0 0 = 0 := by
        calc
          (C 0 1 - C 1 0) * D 0 0 = C 0 1 * D 1 1 - C 1 0 * D 0 0 := by
            rw [hD00eq]
            ring
          _ = 0 := by
            rw [hD10, hD01] at hC01
            simpa [mul_comm] using hC01
      exact sub_eq_zero.mp ((mul_eq_zero.mp this).resolve_right hD00)
    have hC11eq : C 0 2 = C 1 1 := by
      have : (C 0 2 - C 1 1) * D 0 0 = 0 := by
        calc
          (C 0 2 - C 1 1) * D 0 0 = C 0 2 * D 1 1 - C 1 1 * D 0 0 := by
            rw [hD00eq]
            ring
          _ = 0 := by
            rw [hD10, hD01] at hC02
            simpa [mul_comm] using hC02
      exact sub_eq_zero.mp ((mul_eq_zero.mp this).resolve_right hD00)
    have hAshape : A = (A 0 0) • (1 : Matrix W W k) := by
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [Matrix.smul_apply, hA01, hA02, hA10, hA11eq, hA12, hA20, hA21, hA22eq]
    have hCshape : C = lowerHankel (k := k) (C 0 0) (C 0 1) (C 0 2) (C 1 2) := by
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [lowerHankel, hC10eq, hC11eq]
    have hDshape : D = (A 0 0)⁻¹ • (1 : Matrix I I k) := by
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [Matrix.smul_apply, hD01, hD10, hD00eq, hD11inv]
    refine ⟨hA00, ?_⟩
    rw [pointwiseShape, hAshape, hCshape, hDshape]
    simp [lowerHankel, Matrix.smul_apply]
  · rintro ⟨ha, hshape⟩
    rw [hshape]
    exact pointwiseShape_fixes (k := k) (A 0 0) (C 0 0) (C 0 1) (C 0 2) (C 1 2) ha

/-- The determinant parameter for a `2 x 2` matrix. -/
def Delta (a b c d : k) : k := a * d - b * c

/-- The raw symmetric-square lift on the `3`-dimensional block. -/
def sym2Raw (a b c d : k) : Matrix W W k :=
  !![a * a, a * b, b * b;
     2 * a * c, a * d + b * c, 2 * b * d;
     c * c, c * d, d * d]

/-- The transposed adjugate on the `2`-dimensional block. -/
def adjointT (a b c d : k) : Matrix I I k :=
  !![d, -c;
      -b, a]

/-- The raw symmetric-square lift acts on `mulX` and `mulY` with the expected
determinant-scaled `GL₂` coefficients. -/
theorem sym2Raw_action
    (a b c d : k) :
    sym2Raw (k := k) a b c d * (mulX (k := k)) * (adjointT (k := k) a b c d)ᵀ =
        (Delta a b c d) • (a • (mulX (k := k)) + c • (mulY (k := k))) ∧
      sym2Raw (k := k) a b c d * (mulY (k := k)) * (adjointT (k := k) a b c d)ᵀ =
        (Delta a b c d) • (b • (mulX (k := k)) + d • (mulY (k := k))) := by
  constructor
  ·
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [sym2Raw, adjointT, mulX, mulY, Delta, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_two, Fin.sum_univ_three, Matrix.smul_apply] <;>
      ring
  ·
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [sym2Raw, adjointT, mulX, mulY, Delta, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_two, Fin.sum_univ_three, Matrix.smul_apply] <;>
      ring

/-- Recombining the off-diagonal block expression gives the expected linear combination
of the two chosen basis vectors. -/
theorem linearCombination_fromBlocks
    (a b : k) :
    Matrix.fromBlocks
      0
      (a • (mulX (k := k)) + b • (mulY (k := k)))
      (-((a • (mulX (k := k)) + b • (mulY (k := k)))ᵀ))
      0 =
      a • rep₁ + b • rep₂ := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    simp [rep₁, rep₂, mulX, mulY, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
  ·
    fin_cases i <;> fin_cases j <;>
      simp [rep₁, rep₂, mulX, mulY, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] <;>
      ring
  ·
    fin_cases i <;> fin_cases j <;>
      simp [rep₁, rep₂, mulX, mulY, Matrix.fromBlocks, Matrix.add_apply,
        Matrix.smul_apply, Matrix.transpose_apply] <;>
      ring
  ·
    simp [rep₁, rep₂, mulX, mulY, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]

/-- Pulling out a scalar from the off-diagonal block matches scalar multiplication
of the corresponding bivector. -/
theorem smul_fromBlocks_offDiag
    (t : k)
    (M : Matrix W I k) :
    Matrix.fromBlocks 0 (t • M) (-((t • M)ᵀ)) 0 =
      t • Matrix.fromBlocks 0 M (-Mᵀ) 0 := by
  ext i j
  cases i <;> cases j <;>
    simp [Matrix.fromBlocks, Matrix.smul_apply, Matrix.transpose_smul]

/-- A determinant-scaled `GL₂` lift for the pure singular length-two block. -/
theorem GL2_scaled_lift_action
    (a b c d : k) :
    let P : Matrix W W k := sym2Raw (k := k) a b c d
    let Q : Matrix I I k := adjointT (k := k) a b c d
    let g : Matrix V V k := Matrix.fromBlocks P 0 0 Q
    ActBivector rep₁ g = (Delta a b c d) • (a • rep₁ + c • rep₂) ∧
      ActBivector rep₂ g = (Delta a b c d) • (b • rep₁ + d • rep₂) := by
  intro P Q g
  have hPQ := sym2Raw_action (k := k) a b c d
  constructor
  ·
    calc
      ActBivector rep₁ g =
        Matrix.fromBlocks
          0
          (P * (mulX (k := k)) * Qᵀ)
          (-(P * (mulX (k := k)) * Qᵀ)ᵀ)
          0 := by
            simpa [g, P, Q, rep₁, mulX] using
              act_blockDiagonal (k := k) P Q (mulX (k := k))
      _ = Matrix.fromBlocks
            0
            ((Delta a b c d) • (a • (mulX (k := k)) + c • (mulY (k := k))))
            (-(((Delta a b c d) • (a • (mulX (k := k)) + c • (mulY (k := k))))ᵀ))
            0 := by
            simpa [P, Q] using congrArg (fun M => Matrix.fromBlocks 0 M (-Mᵀ) 0) hPQ.1
      _ = (Delta a b c d) • (a • rep₁ + c • rep₂) := by
            rw [smul_fromBlocks_offDiag (k := k) (Delta a b c d) (a • (mulX (k := k)) + c • (mulY (k := k)))]
            exact congrArg (fun M => (Delta a b c d) • M)
              (linearCombination_fromBlocks (k := k) a c)
  ·
    calc
      ActBivector rep₂ g =
        Matrix.fromBlocks
          0
          (P * (mulY (k := k)) * Qᵀ)
          (-(P * (mulY (k := k)) * Qᵀ)ᵀ)
          0 := by
            simpa [g, P, Q, rep₂, mulY] using
              act_blockDiagonal (k := k) P Q (mulY (k := k))
      _ = Matrix.fromBlocks
            0
            ((Delta a b c d) • (b • (mulX (k := k)) + d • (mulY (k := k))))
            (-(((Delta a b c d) • (b • (mulX (k := k)) + d • (mulY (k := k))))ᵀ))
            0 := by
            simpa [P, Q] using congrArg (fun M => Matrix.fromBlocks 0 M (-Mᵀ) 0) hPQ.2
      _ = (Delta a b c d) • (b • rep₁ + d • rep₂) := by
            rw [smul_fromBlocks_offDiag (k := k) (Delta a b c d) (b • (mulX (k := k)) + d • (mulY (k := k)))]
            exact congrArg (fun M => (Delta a b c d) • M)
              (linearCombination_fromBlocks (k := k) b d)

end N5PureSingularLong
end Wedge2Formalization
