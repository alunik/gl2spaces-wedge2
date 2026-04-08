import Wedge2Formalization.N5Summary
import Wedge2Formalization.N3PureSingular

namespace Wedge2Formalization
namespace N5PaperSummary

open Matrix

variable {k : Type*} [Field k]

/-!
# Paper-Facing `n = 5` Summary

This file repackages the finished geometric `n = 5` results in the row order of
`paper_wedge2_rewrite/sections/geometric_orbit_tables_generated.tex`.

The aim is the same as for `N4PaperSummary`: expose the reader-facing kernel and
quotient pieces with names that match the paper, instead of making the reader
reverse-engineer them from construction-heavy theorem names.
-/

/-! ## Row 1

Representative
`⟨e₁ ∧ e₄ + e₂ ∧ e₅, e₂ ∧ e₄ + e₃ ∧ e₅⟩`,
divisor type `0`,
paper stabilizer `K_L = U₄ ⋊ G_m`, `Q_L = GL₂(k)`.

Lean model: `N5PureSingularLong`.
-/

/-- The `U₄` Hankel family on the pure singular length-two row. -/
abbrev row1_U4 := N5PureSingularLong.lowerHankel (k := k)

/-- The pointwise torus on the pure singular length-two row. -/
def row1_torus (t : k) : Matrix N5PureSingularLong.V N5PureSingularLong.V k :=
  Matrix.fromBlocks
    ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
    0
    0
    (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))

/-- Concrete `U₄ ⋊ G_m` pointwise family on the pure singular length-two row. -/
abbrev row1_pointwise_family :=
  N5Summary.pureSingularLong_scaleLowerHankel_product_pointwise (k := k)

/-- The full `GL₂(k)` quotient action on the pure singular length-two row. -/
abbrev row1_GL2_lift_action := N5Summary.pureSingularLong_GL2_scaled_lift_action (k := k)

/-! ## Row 2

Representative
`⟨e₁ ∧ e₃, e₂ ∧ e₃⟩`,
divisor type `0`,
paper stabilizer `K_L = U₆ ⋊ (GL₂(k) × G_m)`, `Q_L = GL₂(k)`.

Lean model: `N5PureSingular`.
-/

/-- Paper-facing shape predicate for the radical pure singular row. -/
def row2_pointwiseShape
    (R : Matrix N5PureSingular.I N5PureSingular.W k)
    (D : Matrix N4PureSingular.I N4PureSingular.I k) : Prop :=
  R = N5PureSingular.upperRightLast (k := k) (R 0 3) ∧
    D 0 0 ≠ 0 ∧
    D =
      N4PureSingular.pureSingularShape
        (k := k) (D 0 0) (D 0 1) (D 0 2) (D 0 3) (D 1 3) (D 2 3) (D 3 3)

/-- Exact pointwise criterion for the radical pure singular row. -/
theorem row2_pointwise_iff
    (a : Matrix N5PureSingular.I N5PureSingular.I k)
    (R : Matrix N5PureSingular.I N5PureSingular.W k)
    (C : Matrix N5PureSingular.W N5PureSingular.I k)
    (D : Matrix N4PureSingular.I N4PureSingular.I k) :
    N5PureSingular.FixesRadPureSingularPairBivector
      (Matrix.fromBlocks a R C D) ↔
      row2_pointwiseShape (k := k) R D :=
  N5Summary.radPureSingular_pointwise_bivector_iff
    (k := k) (a := a) (R := R) (C := C) (D := D)

/-- The full `GL₂(k)` quotient action on the radical pure singular row. -/
abbrev row2_GL2_lift_action := N5Summary.radPureSingular_GL2_lift_action (k := k)

/-! ## Row 3

Representative
`⟨e₁ ∧ e₃ + e₄ ∧ e₅, e₂ ∧ e₃⟩`,
divisor type `[a]`,
paper stabilizer `K_L = U₄ ⋊ (G_m × SL₂(k))`, `Q_L = B`.

Lean model: `N5SimplePoint`.
-/

/-- The `U₄` unipotent family on the simple-point row. -/
abbrev row3_U4 := N5SimplePoint.pointwiseUnipotent (k := k)

/-- Concrete pointwise product family on the simple-point row. -/
abbrev row3_pointwise_family := N5Summary.simplePoint_pointwise_product_family (k := k)

/-- The exact shape theorem once the pure singular block is normalized. -/
abbrev row3_pointwise_iff := N5Summary.simplePoint_pureShape_iff (k := k)

/-- The Borel quotient action on the simple-point row. -/
abbrev row3_borel_lift_action := N5Summary.simplePoint_borel_lift_action (k := k)

/-! ## Row 4

Representative
`⟨e₂ ∧ e₄ + e₃ ∧ e₅, e₂ ∧ e₅⟩`,
divisor type `2[a]`,
paper stabilizer `K_L = U₇ ⋊ (G_m × SL₂(k))`, `Q_L = B`.

Lean model: `N5OnePoint`.
-/

/-- Paper-facing shape predicate for the radical one-point row. -/
def row4_pointwiseShape
    (R : Matrix N5.I N5.W k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) : Prop :=
  R = 0 ∧ A.det = 1 ∧ C₁ = 0 ∧ D = A ∧ A * N4.J * B₁ᵀ + B₁ * N4.J * Aᵀ = 0

/-- Exact pointwise criterion for the radical one-point row. -/
theorem row4_pointwise_iff
    (a : Matrix N5.I N5.I k)
    (R : Matrix N5.I N5.W k)
    (C : Matrix N5.W N5.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N5OnePoint.FixesRadOnePointPairBivector
      (Matrix.fromBlocks a R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      row4_pointwiseShape (k := k) R A B₁ C₁ D :=
  N5Summary.radOnePoint_pointwise_bivector_iff
    (k := k) (a := a) (R := R) (C := C) (A := A) (B₁ := B₁) (C₁ := C₁) (D := D)

/-- The Borel quotient action on the radical one-point row. -/
abbrev row4_borel_lift_action := N5Summary.radOnePoint_borel_lift_action (k := k)

/-! ## Row 5

Representative
`⟨e₂ ∧ e₃ + e₄ ∧ e₅, e₄ ∧ e₅⟩`,
divisor type `[a]+[b]`,
paper stabilizer `K_L = U₄ ⋊ (GL₁(k) × SL₂(k) × SL₂(k))`, `Q_L = N_T`.

Lean model: `N5`.
-/

/-- Paper-facing shape predicate for the radical split row. -/
def row5_pointwiseShape
    (R : Matrix N5.I N5.W k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) : Prop :=
  R = 0 ∧ A.det = 1 ∧ B₁ = 0 ∧ C₁ = 0 ∧ D.det = 1

/-- Exact pointwise criterion for the radical split row. -/
theorem row5_pointwise_iff
    (a : Matrix N5.I N5.I k)
    (R : Matrix N5.I N5.W k)
    (C : Matrix N5.W N5.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N5.FixesRadSplitPairBivector
      (Matrix.fromBlocks a R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      row5_pointwiseShape (k := k) R A B₁ C₁ D :=
  N5Summary.radSplit_pointwise_bivector_iff
    (k := k) (a := a) (R := R) (C := C) (A := A) (B₁ := B₁) (C₁ := C₁) (D := D)

/-- The split torus action on the quotient `P(L)` for the radical split row. -/
abbrev row5_torus_lift_action := N5Summary.radSplit_torus_lift_action (k := k)

/-- The swap coset generating the split-torus normalizer on the radical split row. -/
abbrev row5_swapCoset_lift_action := N5Summary.radSplit_swapCoset_lift_action (k := k)

end N5PaperSummary
end Wedge2Formalization
