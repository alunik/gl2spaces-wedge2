import Wedge2Formalization.PaperCore
import Mathlib.LinearAlgebra.Matrix.SchurComplement

namespace Wedge2Formalization
namespace Paper

open Matrix

variable {V : Type*} [Fintype V] [DecidableEq V]
variable {k : Type*} [Field k]

/-!
# Paper Core

This module is the public paper-facing vocabulary for the row-by-row formalization.
It intentionally stays thin: the public layer should expose exact matrix-family
statements, not abstract algebraic-group packaging.
-/

/-- Public alias for the literal pointwise stabilizer attached to a fixing predicate. -/
abbrev PointwiseStabilizer
    (Fixes : Matrix V V k → Prop) : Set (Matrix V V k) :=
  Wedge2Formalization.PaperCore.PointwiseStabilizer Fixes

/-- Public alias for exactness of a named matrix family for a fixing predicate. -/
abbrev ExactFamily
    (Fixes : Matrix V V k → Prop)
    (K : Set (Matrix V V k)) : Prop :=
  Wedge2Formalization.PaperCore.ExactFamily Fixes K

/-- Public alias for the coefficient action on an ordered basis of a `2`-space. -/
abbrev ActsOnOrderedPair
    (Act : Matrix V V k → Matrix V V k → Matrix V V k)
    (Ω₁ Ω₂ : Matrix V V k)
    (g : Matrix V V k)
    (M : Matrix (Fin 2) (Fin 2) k) : Prop :=
  Wedge2Formalization.PaperCore.ActsOnOrderedPair Act Ω₁ Ω₂ g M

/-- Public alias for setwise preservation of the span of an ordered pair. -/
abbrev PreservesSpanPair
    (Act : Matrix V V k → Matrix V V k → Matrix V V k)
    (Ω₁ Ω₂ : Matrix V V k)
    (g : Matrix V V k) : Prop :=
  Wedge2Formalization.PaperCore.PreservesSpanPair Act Ω₁ Ω₂ g

/-- For the natural `g Ω gᵀ` action on bivectors, an invertible matrix cannot send a
nonzero bivector to zero. -/
theorem actBivector_eq_zero_iff_of_det_ne_zero
    (Ω : Matrix V V k)
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0) :
    g * Ω * gᵀ = 0 ↔ Ω = 0 := by
  letI : Invertible (Matrix.det g) := invertibleOfNonzero hg
  letI : Invertible g := Matrix.invertibleOfDetInvertible g
  have hunit : IsUnit (Matrix.det g) := isUnit_of_invertible (Matrix.det g)
  constructor
  · intro hzero
    have h' := congrArg (fun M => g⁻¹ * M * (g⁻¹)ᵀ) hzero
    simp [Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hunit,
      Matrix.mul_nonsing_inv _ hunit, Matrix.transpose_nonsing_inv] at h'
    exact h'
  · intro hΩ
    simp [hΩ]

end Paper
end Wedge2Formalization
