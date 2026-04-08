import Mathlib.Data.Matrix.Reflection
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

open Matrix

namespace Wedge2Formalization
namespace N4PureSingular

variable {k : Type*} [Field k]

abbrev I := Fin 4

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix I I k) (g : Matrix I I k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix I I k) (g : Matrix I I k) : Matrix I I k :=
  g * Ω * gᵀ

/-- The simple bivector `e₁ ∧ e₂`. -/
def ω12 : Matrix I I k :=
  !![(0 : k), 1, 0, 0;
    -1, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0]

/-- The simple bivector `e₁ ∧ e₃`. -/
def ω13 : Matrix I I k :=
  !![(0 : k), 0, 1, 0;
     0, 0, 0, 0;
    -1, 0, 0, 0;
     0, 0, 0, 0]

/-- The pure singular `n = 4` representative `⟨e₁∧e₂, e₁∧e₃⟩`. -/
def FixesPureSingularPair (g : Matrix I I k) : Prop :=
  FixesBivector ω12 g ∧ FixesBivector ω13 g

/-- Setwise stabilizer of the pure singular `2`-space in the chosen basis. -/
def PreservesPureSingularSubspaceBivector (g : Matrix I I k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector ω12 g = α • ω12 + β • ω13 ∧
    ActBivector ω13 g = γ • ω12 + δ • ω13

def pureSingularSetwiseShape (a b c d p q r s e f t : k) : Matrix I I k :=
  !![a, b, c, d;
     0, p, q, e;
     0, r, s, f;
     0, 0, 0, t]

def pureSingularShape (a b c d e f t : k) : Matrix I I k :=
  !![a, b, c, d;
     0, a⁻¹, 0, e;
     0, 0, a⁻¹, f;
     0, 0, 0, t]

lemma pureSingularSetwiseShape_act_ω12 (a b c d p q r s e f t : k) :
    ActBivector ω12 (pureSingularSetwiseShape (k := k) a b c d p q r s e f t) =
      (a * p) • ω12 + (a * r) • ω13 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [ActBivector, ω12, ω13, pureSingularSetwiseShape, Matrix.mul_apply, Fin.sum_univ_four,
      Matrix.add_apply, Matrix.smul_apply] <;> ring_nf

lemma pureSingularSetwiseShape_act_ω13 (a b c d p q r s e f t : k) :
    ActBivector ω13 (pureSingularSetwiseShape (k := k) a b c d p q r s e f t) =
      (a * q) • ω12 + (a * s) • ω13 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [ActBivector, ω12, ω13, pureSingularSetwiseShape, Matrix.mul_apply, Fin.sum_univ_four,
      Matrix.add_apply, Matrix.smul_apply] <;> ring_nf

theorem pureSingularSetwiseShape_preserves
    (a b c d p q r s e f t : k) :
    PreservesPureSingularSubspaceBivector
      (pureSingularSetwiseShape (k := k) a b c d p q r s e f t) := by
  refine ⟨a * p, a * r, a * q, a * s, ?_, ?_⟩
  · exact pureSingularSetwiseShape_act_ω12 (k := k) a b c d p q r s e f t
  · exact pureSingularSetwiseShape_act_ω13 (k := k) a b c d p q r s e f t

lemma pureSingularShape_fixes (a b c d e f t : k) (ha : a ≠ 0) :
    FixesPureSingularPair (pureSingularShape (k := k) a b c d e f t) := by
  constructor
  ·
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [FixesBivector, ω12, pureSingularShape, Matrix.mul_apply, Fin.sum_univ_four, ha] <;>
      ring_nf
  ·
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [FixesBivector, ω13, pureSingularShape, Matrix.mul_apply, Fin.sum_univ_four, ha] <;>
      ring_nf

theorem fixesPureSingularPair_iff_shape (g : Matrix I I k) :
    FixesPureSingularPair g ↔
      g 0 0 ≠ 0 ∧
        g = pureSingularShape (k := k) (g 0 0) (g 0 1) (g 0 2) (g 0 3) (g 1 3) (g 2 3) (g 3 3) := by
  constructor
  · rintro ⟨h12, h13⟩
    have h12_01 := congrArg (fun M => M 0 1) h12
    have h12_02 := congrArg (fun M => M 0 2) h12
    have h12_03 := congrArg (fun M => M 0 3) h12
    have h12_12 := congrArg (fun M => M 1 2) h12
    have h12_13 := congrArg (fun M => M 1 3) h12
    have h13_01 := congrArg (fun M => M 0 1) h13
    have h13_02 := congrArg (fun M => M 0 2) h13
    have h13_03 := congrArg (fun M => M 0 3) h13
    have h13_12 := congrArg (fun M => M 1 2) h13
    have h1201 : g 0 0 * g 1 1 - g 0 1 * g 1 0 = 1 := by
      simp [FixesBivector, ω12, Matrix.mul_apply, Fin.sum_univ_four] at h12_01
      ring_nf at h12_01 ⊢
      exact h12_01
    have h1202 : g 0 0 * g 2 1 - g 0 1 * g 2 0 = 0 := by
      simp [FixesBivector, ω12, Matrix.mul_apply, Fin.sum_univ_four] at h12_02
      ring_nf at h12_02 ⊢
      exact h12_02
    have h1203 : g 0 0 * g 3 1 - g 0 1 * g 3 0 = 0 := by
      simp [FixesBivector, ω12, Matrix.mul_apply, Fin.sum_univ_four] at h12_03
      ring_nf at h12_03 ⊢
      exact h12_03
    have h1212 : g 1 0 * g 2 1 - g 1 1 * g 2 0 = 0 := by
      simp [FixesBivector, ω12, Matrix.mul_apply, Fin.sum_univ_four] at h12_12
      ring_nf at h12_12 ⊢
      exact h12_12
    have h1213 : g 1 0 * g 3 1 - g 1 1 * g 3 0 = 0 := by
      simp [FixesBivector, ω12, Matrix.mul_apply, Fin.sum_univ_four] at h12_13
      ring_nf at h12_13 ⊢
      exact h12_13
    have h1301 : g 0 0 * g 1 2 - g 0 2 * g 1 0 = 0 := by
      simp [FixesBivector, ω13, Matrix.mul_apply, Fin.sum_univ_four] at h13_01
      ring_nf at h13_01 ⊢
      exact h13_01
    have h1302 : g 0 0 * g 2 2 - g 0 2 * g 2 0 = 1 := by
      simp [FixesBivector, ω13, Matrix.mul_apply, Fin.sum_univ_four] at h13_02
      ring_nf at h13_02 ⊢
      exact h13_02
    have h1303 : g 0 0 * g 3 2 - g 0 2 * g 3 0 = 0 := by
      simp [FixesBivector, ω13, Matrix.mul_apply, Fin.sum_univ_four] at h13_03
      ring_nf at h13_03 ⊢
      exact h13_03
    have h1312 : g 1 0 * g 2 2 - g 1 2 * g 2 0 = 0 := by
      simp [FixesBivector, ω13, Matrix.mul_apply, Fin.sum_univ_four] at h13_12
      ring_nf at h13_12 ⊢
      exact h13_12
    have h1212' : g 2 0 * g 1 1 - g 1 0 * g 2 1 = 0 := by
      calc
        g 2 0 * g 1 1 - g 1 0 * g 2 1 = -(g 1 0 * g 2 1 - g 1 1 * g 2 0) := by ring
        _ = 0 := by rw [h1212]; ring
    have h1213' : g 3 0 * g 1 1 - g 1 0 * g 3 1 = 0 := by
      calc
        g 3 0 * g 1 1 - g 1 0 * g 3 1 = -(g 1 0 * g 3 1 - g 1 1 * g 3 0) := by ring
        _ = 0 := by rw [h1213]; ring
    have h1312' : g 1 0 * g 2 2 - g 2 0 * g 1 2 = 0 := by
      calc
        g 1 0 * g 2 2 - g 2 0 * g 1 2 = g 1 0 * g 2 2 - g 1 2 * g 2 0 := by ring
        _ = 0 := h1312
    have h20 : g 2 0 = 0 := by
      calc
        g 2 0 = g 2 0 * (g 0 0 * g 1 1 - g 0 1 * g 1 0) - g 1 0 * (g 0 0 * g 2 1 - g 0 1 * g 2 0) := by
          rw [h1201, h1202]
          ring
        _ = g 0 0 * (g 2 0 * g 1 1 - g 1 0 * g 2 1) := by ring
        _ = 0 := by rw [h1212']; ring
    have h30 : g 3 0 = 0 := by
      calc
        g 3 0 = g 3 0 * (g 0 0 * g 1 1 - g 0 1 * g 1 0) - g 1 0 * (g 0 0 * g 3 1 - g 0 1 * g 3 0) := by
          rw [h1201, h1203]
          ring
        _ = g 0 0 * (g 3 0 * g 1 1 - g 1 0 * g 3 1) := by ring
        _ = 0 := by rw [h1213']; ring
    have h10 : g 1 0 = 0 := by
      calc
        g 1 0 = g 1 0 * (g 0 0 * g 2 2 - g 0 2 * g 2 0) - g 2 0 * (g 0 0 * g 1 2 - g 0 2 * g 1 0) := by
          rw [h1302, h1301]
          ring
        _ = g 0 0 * (g 1 0 * g 2 2 - g 2 0 * g 1 2) := by ring
        _ = 0 := by rw [h1312']; ring
    have h00 : g 0 0 ≠ 0 := by
      intro h0
      have : (0 : k) = 1 := by simpa [h0, h10] using h1201
      exact zero_ne_one this
    have h1201' : g 0 0 * g 1 1 = 1 := by simpa [h10] using h1201
    have h1302' : g 0 0 * g 2 2 = 1 := by simpa [h20] using h1302
    have h21 : g 2 1 = 0 := by
      have : g 0 0 * g 2 1 = 0 := by simpa [h20] using h1202
      exact mul_eq_zero.mp (by simpa using this) |>.resolve_left h00
    have h31 : g 3 1 = 0 := by
      have : g 0 0 * g 3 1 = 0 := by simpa [h30] using h1203
      exact mul_eq_zero.mp (by simpa using this) |>.resolve_left h00
    have h12z : g 1 2 = 0 := by
      have : g 0 0 * g 1 2 = 0 := by simpa [h10] using h1301
      exact mul_eq_zero.mp (by simpa using this) |>.resolve_left h00
    have h32 : g 3 2 = 0 := by
      have : g 0 0 * g 3 2 = 0 := by simpa [h30] using h1303
      exact mul_eq_zero.mp (by simpa using this) |>.resolve_left h00
    have h11 : g 1 1 = (g 0 0)⁻¹ := by
      calc
        g 1 1 = 1 * g 1 1 := by simp
        _ = ((g 0 0)⁻¹ * g 0 0) * g 1 1 := by rw [inv_mul_cancel₀ h00]
        _ = (g 0 0)⁻¹ * (g 0 0 * g 1 1) := by ring
        _ = (g 0 0)⁻¹ := by rw [h1201']; simp
    have h22 : g 2 2 = (g 0 0)⁻¹ := by
      calc
        g 2 2 = 1 * g 2 2 := by simp
        _ = ((g 0 0)⁻¹ * g 0 0) * g 2 2 := by rw [inv_mul_cancel₀ h00]
        _ = (g 0 0)⁻¹ * (g 0 0 * g 2 2) := by ring
        _ = (g 0 0)⁻¹ := by rw [h1302']; simp
    refine ⟨h00, ?_⟩
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [pureSingularShape, h10, h20, h30, h21, h31, h12z, h32, h11, h22]
  · rintro ⟨h00, hshape⟩
    rw [hshape]
    exact pureSingularShape_fixes (k := k) _ _ _ _ _ _ _ h00

end N4PureSingular
end Wedge2Formalization
