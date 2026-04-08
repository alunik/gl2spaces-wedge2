import Mathlib.Data.Matrix.Reflection
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

open Matrix

namespace Wedge2Formalization
namespace N3PureSingular

variable {k : Type*} [Field k]

/-- The ambient `3`-dimensional space. -/
abbrev I := Fin 3

/-- The simple bivector `e₁ ∧ e₃`. -/
def ω13 : Matrix I I k :=
  !![(0 : k), 0, 1;
      0, 0, 0;
     -1, 0, 0]

/-- The simple bivector `e₂ ∧ e₃`. -/
def ω23 : Matrix I I k :=
  !![(0 : k), 0, 0;
      0, 0, 1;
      0, -1, 0]

/-- The natural action on bivectors. -/
def ActBivector (Ω : Matrix I I k) (g : Matrix I I k) : Matrix I I k :=
  g * Ω * gᵀ

/-- Exact stabilizer of a bivector. -/
def FixesBivector (Ω : Matrix I I k) (g : Matrix I I k) : Prop :=
  ActBivector Ω g = Ω

/-- Exact stabilizer of the pure singular pair. -/
def FixesPairBivector (g : Matrix I I k) : Prop :=
  FixesBivector ω13 g ∧ FixesBivector ω23 g

/-- Setwise stabilizer of the pure singular `2`-space in the chosen basis. -/
def PreservesSubspaceBivector (g : Matrix I I k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector ω13 g = α • ω13 + β • ω23 ∧
      ActBivector ω23 g = γ • ω13 + δ • ω23

/-- A concrete block upper triangular lift for the quotient action. -/
def pureLift (a b c d t : k) : Matrix I I k :=
  !![a, c, 0;
      b, d, 0;
      0, 0, t]

/-- The one-parameter pointwise family. -/
def pointwiseScale (a : k) : Matrix I I k :=
  !![a⁻¹, 0, 0;
      0, a⁻¹, 0;
      0, 0, a]

/-- The full pointwise stabilizer shape of the pure singular `3`-dimensional pair. -/
def pureSingularShape (a x y : k) : Matrix I I k :=
  !![a, 0, 0;
      0, a, 0;
      x, y, a⁻¹]

theorem pureLift_act_ω13
    (a b c d t : k) :
    ActBivector ω13 (pureLift (k := k) a b c d t) =
      (a * t) • ω13 + (b * t) • ω23 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [ActBivector, ω13, ω23, pureLift, Matrix.mul_apply, Fin.sum_univ_three,
      Matrix.add_apply, Matrix.smul_apply] <;> ring_nf

theorem pureLift_act_ω23
    (a b c d t : k) :
    ActBivector ω23 (pureLift (k := k) a b c d t) =
      (c * t) • ω13 + (d * t) • ω23 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [ActBivector, ω13, ω23, pureLift, Matrix.mul_apply, Fin.sum_univ_three,
      Matrix.add_apply, Matrix.smul_apply] <;> ring_nf

/-- Every concrete upper triangular lift preserves the pure singular `2`-space. -/
theorem pureLift_preserves
    (a b c d t : k) :
    PreservesSubspaceBivector (pureLift (k := k) a b c d t) := by
  refine ⟨a * t, b * t, c * t, d * t, ?_, ?_⟩
  · exact pureLift_act_ω13 (k := k) a b c d t
  · exact pureLift_act_ω23 (k := k) a b c d t

/-- The determinant of the concrete pure singular lift is `t(ad-bc)`. -/
theorem pureLift_det
    (a b c d t : k) :
    Matrix.det (pureLift (k := k) a b c d t) = t * (a * d - b * c) := by
  simp [pureLift, Matrix.det_fin_three]
  ring

/-- The determinant of the full pure singular pointwise shape is the scalar `a`. -/
theorem pureSingularShape_det
    (a x y : k) :
    Matrix.det (pureSingularShape (k := k) a x y) = a := by
  by_cases ha : a = 0
  · simp [pureSingularShape, Matrix.det_fin_three, ha]
  · simp [pureSingularShape, Matrix.det_fin_three, ha, mul_assoc]

/-- A convenient `GL₂` lift for the quotient action. -/
theorem GL2_lift_action
    (α β γ δ : k) :
    let g := pureLift (k := k) α β γ δ (1 : k)
    ActBivector ω13 g = α • ω13 + β • ω23 ∧
      ActBivector ω23 g = γ • ω13 + δ • ω23 := by
  intro g
  constructor
  · simpa using pureLift_act_ω13 (k := k) α β γ δ (1 : k)
  · simpa using pureLift_act_ω23 (k := k) α β γ δ (1 : k)

/-- The basic two-parameter unipotent family in the pointwise stabilizer. -/
def pointwiseUnipotent (x y : k) : Matrix I I k :=
  !![(1 : k), 0, 0;
      0, 1, 0;
      x, y, 1]

/-- The obvious `G_m` family fixes the pair pointwise. -/
theorem pointwise_scale_family
    (a : k)
    (ha : a ≠ 0) :
    FixesPairBivector (pointwiseScale (k := k) a) := by
  constructor
  ·
    rw [FixesBivector, pointwiseScale]
    simpa [ha] using pureLift_act_ω13 (k := k) a⁻¹ 0 0 a⁻¹ a
  ·
    rw [FixesBivector, pointwiseScale]
    simpa [ha] using pureLift_act_ω23 (k := k) a⁻¹ 0 0 a⁻¹ a

/-- The basic two-parameter unipotent family fixes the pure singular pair pointwise. -/
theorem pointwise_unipotent_family
    (x y : k) :
    FixesPairBivector (pointwiseUnipotent (k := k) x y) := by
  constructor
  ·
    rw [FixesBivector, pointwiseUnipotent]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [ActBivector, ω13, pointwiseUnipotent, Matrix.mul_apply, Fin.sum_univ_three]
      <;> ring_nf
  ·
    rw [FixesBivector, pointwiseUnipotent]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [ActBivector, ω23, pointwiseUnipotent, Matrix.mul_apply, Fin.sum_univ_three]
      <;> ring_nf

/-- Every matrix of the expected pointwise shape fixes the pure singular pair. -/
theorem pureSingularShape_fixes
    (a x y : k)
    (ha : a ≠ 0) :
    FixesPairBivector (pureSingularShape (k := k) a x y) := by
  constructor
  ·
    rw [FixesBivector]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [ActBivector, ω13, pureSingularShape, Matrix.mul_apply, Fin.sum_univ_three, ha]
      <;> ring_nf
  ·
    rw [FixesBivector]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [ActBivector, ω23, pureSingularShape, Matrix.mul_apply, Fin.sum_univ_three, ha]
      <;> ring_nf

/-- The pointwise stabilizer of the pure singular pair has the expected
`G_a^2 ⋊ G_m` matrix shape. -/
theorem fixesPairBivector_iff_shape
    (g : Matrix I I k) :
    FixesPairBivector g ↔
      g 0 0 ≠ 0 ∧
        g = pureSingularShape (k := k) (g 0 0) (g 2 0) (g 2 1) := by
  constructor
  · rintro ⟨h13, h23⟩
    have h13_01 := congrArg (fun M => M 0 1) h13
    have h13_02 := congrArg (fun M => M 0 2) h13
    have h13_12 := congrArg (fun M => M 1 2) h13
    have h23_01 := congrArg (fun M => M 0 1) h23
    have h23_02 := congrArg (fun M => M 0 2) h23
    have h23_12 := congrArg (fun M => M 1 2) h23
    have h1301 : g 0 0 * g 1 2 - g 0 2 * g 1 0 = 0 := by
      simp [FixesBivector, ActBivector, ω13, Matrix.mul_apply, Fin.sum_univ_three] at h13_01
      ring_nf at h13_01 ⊢
      exact h13_01
    have h1302 : g 0 0 * g 2 2 - g 0 2 * g 2 0 = 1 := by
      simp [FixesBivector, ActBivector, ω13, Matrix.mul_apply, Fin.sum_univ_three] at h13_02
      ring_nf at h13_02 ⊢
      exact h13_02
    have h1312 : g 1 0 * g 2 2 - g 1 2 * g 2 0 = 0 := by
      simp [FixesBivector, ActBivector, ω13, Matrix.mul_apply, Fin.sum_univ_three] at h13_12
      ring_nf at h13_12 ⊢
      exact h13_12
    have h2301 : g 0 1 * g 1 2 - g 0 2 * g 1 1 = 0 := by
      simp [FixesBivector, ActBivector, ω23, Matrix.mul_apply, Fin.sum_univ_three] at h23_01
      ring_nf at h23_01 ⊢
      exact h23_01
    have h2302 : g 0 1 * g 2 2 - g 0 2 * g 2 1 = 0 := by
      simp [FixesBivector, ActBivector, ω23, Matrix.mul_apply, Fin.sum_univ_three] at h23_02
      ring_nf at h23_02 ⊢
      exact h23_02
    have h2312 : g 1 1 * g 2 2 - g 1 2 * g 2 1 = 1 := by
      simp [FixesBivector, ActBivector, ω23, Matrix.mul_apply, Fin.sum_univ_three] at h23_12
      ring_nf at h23_12 ⊢
      exact h23_12
    have h02 : g 0 2 = 0 := by
      calc
        g 0 2 =
            g 0 2 * (g 1 1 * g 2 2 - g 1 2 * g 2 1) -
              g 1 2 * (g 0 1 * g 2 2 - g 0 2 * g 2 1) := by
              rw [h2312, h2302]
              ring
        _ = g 2 2 * (g 0 2 * g 1 1 - g 1 2 * g 0 1) := by ring
        _ = 0 := by
              calc
                g 2 2 * (g 0 2 * g 1 1 - g 1 2 * g 0 1) =
                    g 2 2 * (-(g 0 1 * g 1 2 - g 0 2 * g 1 1)) := by ring
                _ = 0 := by rw [h2301]; ring
    have h12 : g 1 2 = 0 := by
      calc
        g 1 2 =
            g 1 2 * (g 0 0 * g 2 2 - g 0 2 * g 2 0) -
              g 0 2 * (g 1 0 * g 2 2 - g 1 2 * g 2 0) := by
              rw [h1302, h1312]
              ring
        _ = g 2 2 * (g 1 2 * g 0 0 - g 0 2 * g 1 0) := by ring
        _ = 0 := by
              calc
                g 2 2 * (g 1 2 * g 0 0 - g 0 2 * g 1 0) =
                    g 2 2 * (g 0 0 * g 1 2 - g 0 2 * g 1 0) := by ring
                _ = 0 := by rw [h1301]; ring
    have h00h22 : g 0 0 * g 2 2 = 1 := by simpa [h02] using h1302
    have h11h22 : g 1 1 * g 2 2 = 1 := by simpa [h12] using h2312
    have h00 : g 0 0 ≠ 0 := by
      intro h0
      have : (0 : k) = 1 := by simpa [h0] using h00h22
      exact zero_ne_one this
    have h22 : g 2 2 = (g 0 0)⁻¹ := by
      calc
        g 2 2 = 1 * g 2 2 := by ring
        _ = ((g 0 0)⁻¹ * g 0 0) * g 2 2 := by rw [inv_mul_cancel₀ h00]
        _ = (g 0 0)⁻¹ * (g 0 0 * g 2 2) := by ring
        _ = (g 0 0)⁻¹ := by rw [h00h22]; simp
    have h22nz : g 2 2 ≠ 0 := by
      rw [h22]
      exact inv_ne_zero h00
    have h01 : g 0 1 = 0 := by
      have : g 0 1 * g 2 2 = 0 := by simpa [h02] using h2302
      exact (mul_eq_zero.mp this).resolve_right h22nz
    have h10 : g 1 0 = 0 := by
      have : g 1 0 * g 2 2 = 0 := by simpa [h12] using h1312
      exact (mul_eq_zero.mp this).resolve_right h22nz
    have h11inv : g 1 1 = (g 2 2)⁻¹ := by
      have hleft : g 1 1 * g 2 2 * (g 2 2)⁻¹ = g 1 1 := by
        calc
          g 1 1 * g 2 2 * (g 2 2)⁻¹ = g 1 1 * (g 2 2 * (g 2 2)⁻¹) := by ring
          _ = g 1 1 * 1 := by rw [mul_inv_cancel₀ h22nz]
          _ = g 1 1 := by ring
      have hright : g 1 1 * g 2 2 * (g 2 2)⁻¹ = (g 2 2)⁻¹ := by
        rw [h11h22]
        ring
      calc
        g 1 1 = g 1 1 * g 2 2 * (g 2 2)⁻¹ := by
              symm
              exact hleft
        _ = (g 2 2)⁻¹ := hright
    have h11 : g 1 1 = g 0 0 := by
      calc
        g 1 1 = (g 2 2)⁻¹ := h11inv
        _ = g 0 0 := by rw [h22, inv_inv]
    refine ⟨h00, ?_⟩
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [pureSingularShape, h01, h02, h10, h11, h12, h22]
  · rintro ⟨h00, hshape⟩
    rw [hshape]
    exact pureSingularShape_fixes (k := k) (g 0 0) (g 2 0) (g 2 1) h00

end N3PureSingular
end Wedge2Formalization
