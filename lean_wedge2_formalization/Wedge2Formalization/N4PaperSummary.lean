import Wedge2Formalization.PaperCore
import Wedge2Formalization.N4Summary
import Wedge2Formalization.N3PureSingular

namespace Wedge2Formalization
namespace N4PaperSummary

open Matrix

variable {k : Type*} [Field k]

/-!
# Paper-Facing `n = 4` Summary

This file repackages the finished geometric `n = 4` results in the row order of
`paper_wedge2_rewrite/sections/geometric_orbit_tables_generated.tex`.

The point is readability. Each row now has:

- a paper-facing pointwise-shape predicate when an exact `iff` theorem is available;
- named concrete families for the unipotent / torus / Levi pieces used in the paper;
- a short alias for the quotient-action theorem.
-/

/-! ## Row 1

Representative
`⟨e₁ ∧ e₃ + e₂ ∧ e₄, e₁ ∧ e₄⟩`,
divisor type `2[a]`,
paper stabilizer `K_L = U₃ ⋊ SL₂(k)`, `Q_L = B`.

Lean model: `N4`.
-/

/-- Paper-facing shape predicate for the exact pointwise stabilizer on the `2[a]` row. -/
def row1_pointwiseShape
    (A B C D : Matrix N4.I N4.I k) : Prop :=
  A.det = 1 ∧ C = 0 ∧ D = A ∧ A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0

/-- The exact pointwise kernel for the `2[a]` row, packaged as a named family. -/
def row1_pointwiseKernel : Set (Matrix N4.V N4.V k) :=
  { g | ∃ A B : Matrix N4.I N4.I k,
      A.det = 1 ∧
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0 ∧
      g = Matrix.fromBlocks A B 0 A }

/-- The upper-block shear used for the unipotent part of the repeated-support row. -/
abbrev row1_upperShear := N4.onePointUpperShear (k := k)

/-- The diagonal scaling used for the Levi `SL₂` part of the repeated-support row. -/
abbrev row1_scaling := N4.onePointScale (k := k)

/-- A concrete three-parameter `U₃` family on the repeated-support row. -/
def row1_U3 (x y z : k) : Matrix N4.V N4.V k :=
  N4.onePointUpperShear (k := k) !![x, y; z, -x]

/-- The block-diagonal `SL₂` family on the repeated-support row. -/
def row1_levi (A : Matrix N4.I N4.I k) : Matrix N4.V N4.V k :=
  Matrix.fromBlocks A 0 0 A

/-- A standard coefficient-side section for the projective Borel quotient on the `2[a]` row. -/
def row1_borelCoeff (a b : k) : Matrix (Fin 2) (Fin 2) k :=
  !![a, b; 0, a * a]

/-- A standard lift of the projective Borel section on the repeated-support row. -/
def row1_borelLift (a b : k) : Matrix N4.V N4.V k :=
  let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
  N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B

/-- Exact pointwise criterion for the repeated-support row, in paper-facing form. -/
theorem row1_pointwise_iff
    (A B C D : Matrix N4.I N4.I k) :
    N4.FixesOnePointPairBivector (Matrix.fromBlocks A B C D) ↔
      row1_pointwiseShape (k := k) A B C D :=
  N4Summary.onePoint_pointwise_bivector_iff (k := k) (A := A) (B := B) (C := C) (D := D)

/-- The named kernel family for the repeated-support row is exact. -/
theorem row1_pointwise_stabilizer_iff_mem
    (g : Matrix N4.V N4.V k) :
    N4.FixesOnePointPairBivector g ↔ g ∈ row1_pointwiseKernel (k := k) := by
  let A : Matrix N4.I N4.I k := Matrix.toBlocks₁₁ g
  let B : Matrix N4.I N4.I k := Matrix.toBlocks₁₂ g
  let C : Matrix N4.I N4.I k := Matrix.toBlocks₂₁ g
  let D : Matrix N4.I N4.I k := Matrix.toBlocks₂₂ g
  have hg : g = Matrix.fromBlocks A B C D := by
    simpa [A, B, C, D] using (Matrix.fromBlocks_toBlocks g).symm
  constructor
  · intro hfix
    have hfix' :
        N4.FixesOnePointPairBivector (Matrix.fromBlocks A B C D) := by
      simpa [hg] using hfix
    rcases (row1_pointwise_iff (k := k) (A := A) (B := B) (C := C) (D := D)).1 hfix' with
      ⟨hA, hC, hD, hrel⟩
    refine ⟨A, B, hA, hrel, ?_⟩
    calc
      g = Matrix.fromBlocks A B C D := hg
      _ = Matrix.fromBlocks A B 0 A := by simp [hC, hD]
  · rintro ⟨A, B, hA, hrel, rfl⟩
    exact (row1_pointwise_iff (k := k) (A := A) (B := B) (C := 0) (D := A)).2
      ⟨hA, rfl, rfl, hrel⟩

/-- Paper-core exactness statement for the repeated-support row. -/
theorem row1_pointwise_stabilizer :
    PaperCore.ExactFamily
      (N4.FixesOnePointPairBivector (k := k))
      (row1_pointwiseKernel (k := k)) := by
  intro g
  exact row1_pointwise_stabilizer_iff_mem (k := k) g

/-- The concrete `U₃` family fixes the repeated-support pair pointwise. -/
theorem row1_U3_pointwise
    (x y z : k) :
    N4.FixesOnePointPairBivector (row1_U3 (k := k) x y z) := by
  constructor
  · have h := (N4Summary.onePoint_upperShear_action (k := k) (B := !![x, y; z, -x])).1
    simpa [row1_U3] using h
  · have h := (N4Summary.onePoint_upperShear_action (k := k) (B := !![x, y; z, -x])).2
    simpa [row1_U3] using h

/-- The block-diagonal `SL₂` family fixes the repeated-support pair pointwise. -/
theorem row1_levi_pointwise
    (A : Matrix N4.I N4.I k)
    (hA : A.det = 1) :
    N4.FixesOnePointPairBivector (row1_levi (k := k) A) := by
  exact (row1_pointwise_iff (k := k) (A := A) (B := 0) (C := 0) (D := A)).2
    ⟨hA, rfl, rfl, by simp⟩

/-- The standard projective-Borel section acts on the repeated-support row as expected. -/
abbrev row1_borel_lift_action := N4Summary.onePoint_borel_lift_action (k := k)

/-- Paper-facing quotient-action statement for the repeated-support row. -/
theorem row1_borel_lift
    (a b : k)
    (ha : a ≠ 0) :
    PaperCore.ActsOnOrderedPair
      (N4.ActBivector (k := k))
      (N4.onePointRep₁ (k := k))
      (N4.onePointRep₂ (k := k))
      (row1_borelLift (k := k) a b)
      (row1_borelCoeff (k := k) a b) := by
  simpa [PaperCore.ActsOnOrderedPair, row1_borelLift, row1_borelCoeff] using
    row1_borel_lift_action (k := k) (a := a) (b := b) ha

/-! ## Row 2

Representative
`⟨e₁ ∧ e₂ + e₃ ∧ e₄, e₃ ∧ e₄⟩`,
divisor type `[a]+[b]`,
paper stabilizer `K_L = SL₂(k) × SL₂(k)`, `Q_L = N_T`.

Lean model: `N4`.
-/

/-- Paper-facing shape predicate for the exact pointwise stabilizer on the split row. -/
def row2_pointwiseShape
    (A B C D : Matrix N4.I N4.I k) : Prop :=
  A.det = 1 ∧ B = 0 ∧ C = 0 ∧ D.det = 1

/-- The exact pointwise kernel for the split row, packaged as a named family. -/
def row2_pointwiseKernel : Set (Matrix N4.V N4.V k) :=
  { g | ∃ A D : Matrix N4.I N4.I k,
      A.det = 1 ∧
      D.det = 1 ∧
      g = Matrix.fromBlocks A 0 0 D }

/-- Coefficient-side split torus matrix in the chosen ordered basis. -/
def row2_torusCoeff (u v : k) : Matrix (Fin 2) (Fin 2) k :=
  !![u, v - u; 0, v]

/-- The standard split-torus lift on the split row. -/
def row2_torusLift (u v : k) : Matrix N4.V N4.V k :=
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  Matrix.fromBlocks A 0 0 D

/-- Coefficient-side swap-coset matrix in the chosen ordered basis. -/
def row2_swapCoeff (u v : k) : Matrix (Fin 2) (Fin 2) k :=
  !![u, v - u; u, -u]

/-- The standard swap-coset lift on the split row. -/
def row2_swapLift (u v : k) : Matrix N4.V N4.V k :=
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k)

/-- Exact pointwise criterion for the split row, in paper-facing form. -/
theorem row2_pointwise_iff
    (A B C D : Matrix N4.I N4.I k) :
    N4.FixesSplitPairBivector (Matrix.fromBlocks A B C D) ↔
      row2_pointwiseShape (k := k) A B C D :=
  N4Summary.split_pointwise_bivector_iff (k := k) (A := A) (B := B) (C := C) (D := D)

/-- The named kernel family for the split row is exact. -/
theorem row2_pointwise_stabilizer_iff_mem
    (g : Matrix N4.V N4.V k) :
    N4.FixesSplitPairBivector g ↔ g ∈ row2_pointwiseKernel (k := k) := by
  let A : Matrix N4.I N4.I k := Matrix.toBlocks₁₁ g
  let B : Matrix N4.I N4.I k := Matrix.toBlocks₁₂ g
  let C : Matrix N4.I N4.I k := Matrix.toBlocks₂₁ g
  let D : Matrix N4.I N4.I k := Matrix.toBlocks₂₂ g
  have hg : g = Matrix.fromBlocks A B C D := by
    simpa [A, B, C, D] using (Matrix.fromBlocks_toBlocks g).symm
  constructor
  · intro hfix
    have hfix' :
        N4.FixesSplitPairBivector (Matrix.fromBlocks A B C D) := by
      simpa [hg] using hfix
    rcases (row2_pointwise_iff (k := k) (A := A) (B := B) (C := C) (D := D)).1 hfix' with
      ⟨hA, hB, hC, hD⟩
    refine ⟨A, D, hA, hD, ?_⟩
    calc
      g = Matrix.fromBlocks A B C D := hg
      _ = Matrix.fromBlocks A 0 0 D := by simp [hB, hC]
  · rintro ⟨A, D, hA, hD, rfl⟩
    exact (row2_pointwise_iff (k := k) (A := A) (B := 0) (C := 0) (D := D)).2
      ⟨hA, rfl, rfl, hD⟩

/-- Paper-core exactness statement for the split row. -/
theorem row2_pointwise_stabilizer :
    PaperCore.ExactFamily
      (N4.FixesSplitPairBivector (k := k))
      (row2_pointwiseKernel (k := k)) := by
  intro g
  exact row2_pointwise_stabilizer_iff_mem (k := k) g

/-- The split torus action on the quotient `P(L)`. -/
abbrev row2_torus_lift_action := N4Summary.split_torus_lift_action (k := k)

/-- The swap coset generating the normalizer of the split torus. -/
abbrev row2_swapCoset_lift_action := N4Summary.split_swapCoset_lift_action (k := k)

/-- Paper-facing quotient-action statement for the split torus on the split row. -/
theorem row2_torus_lift
    (u v : k) :
    PaperCore.ActsOnOrderedPair
      (N4.ActBivector (k := k))
      (N4.splitRep₁ (k := k))
      (N4.splitRep₂ (k := k))
      (row2_torusLift (k := k) u v)
      (row2_torusCoeff (k := k) u v) := by
  simpa [PaperCore.ActsOnOrderedPair, row2_torusLift, row2_torusCoeff] using
    row2_torus_lift_action (k := k) (u := u) (v := v)

/-- Paper-facing quotient-action statement for the swap coset on the split row. -/
theorem row2_swap_lift
    (u v : k) :
    PaperCore.ActsOnOrderedPair
      (N4.ActBivector (k := k))
      (N4.splitRep₁ (k := k))
      (N4.splitRep₂ (k := k))
      (row2_swapLift (k := k) u v)
      (row2_swapCoeff (k := k) u v) := by
  simpa [PaperCore.ActsOnOrderedPair, row2_swapLift, row2_swapCoeff] using
    row2_swapCoset_lift_action (k := k) (u := u) (v := v)

/-! ## Row 3

Representative
`⟨e₂ ∧ e₄, e₃ ∧ e₄⟩`,
divisor type `0`,
paper stabilizer `K_L = U₂ ⋊ G_m`, `Q_L = GL₂(k)`.

Lean model: `N4PureSingular`.
-/

/-- Paper-facing shape predicate for the exact pointwise stabilizer on the pure singular row. -/
def row3_pointwiseShape
    (g : Matrix N4PureSingular.I N4PureSingular.I k) : Prop :=
  g 0 0 ≠ 0 ∧
    g =
      N4PureSingular.pureSingularShape (k := k)
        (g 0 0) (g 0 1) (g 0 2) (g 0 3) (g 1 3) (g 2 3) (g 3 3)

/-- The exact pointwise kernel for the pure singular row, packaged as a named family. -/
def row3_pointwiseKernel : Set (Matrix N4PureSingular.I N4PureSingular.I k) :=
  { g | row3_pointwiseShape (k := k) g }

/-- Coefficient-side `GL₂(k)` matrix in the chosen ordered basis on the pure singular row. -/
def row3_GL2Coeff (α β γ δ : k) : Matrix (Fin 2) (Fin 2) k :=
  !![α, β; γ, δ]

/-- The standard `GL₂(k)` lift on the pure singular row. -/
def row3_GL2Lift (α β γ δ : k) : Matrix N4PureSingular.I N4PureSingular.I k :=
  N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1

/-- The unipotent `U₂` family on the pure singular row, viewed on a single `3 × 3` block. -/
abbrev row3_unipotent := N3PureSingular.pointwiseUnipotent (k := k)

/-- The torus `G_m` family on the pure singular row, viewed on a single `3 × 3` block. -/
abbrev row3_torus := N3PureSingular.pointwiseScale (k := k)

/-- Exact pointwise criterion for the pure singular row, in paper-facing form. -/
theorem row3_pointwise_iff
    (g : Matrix N4PureSingular.I N4PureSingular.I k) :
    N4PureSingular.FixesPureSingularPair g ↔ row3_pointwiseShape (k := k) g :=
  N4Summary.pureSingular_pointwise_iff_shape (k := k) (g := g)

/-- The named kernel family for the pure singular row is exact. -/
theorem row3_pointwise_stabilizer_iff_mem
    (g : Matrix N4PureSingular.I N4PureSingular.I k) :
    N4PureSingular.FixesPureSingularPair g ↔ g ∈ row3_pointwiseKernel (k := k) := by
  exact row3_pointwise_iff (k := k) (g := g)

/-- Paper-core exactness statement for the pure singular row. -/
theorem row3_pointwise_stabilizer :
    PaperCore.ExactFamily
      (N4PureSingular.FixesPureSingularPair (k := k))
      (row3_pointwiseKernel (k := k)) := by
  intro g
  exact row3_pointwise_stabilizer_iff_mem (k := k) g

/-- The full `GL₂(k)` quotient action on the pure singular row. -/
abbrev row3_GL2_lift_action := N4Summary.pureSingular_GL2_lift_action (k := k)

/-- Paper-facing quotient-action statement for the pure singular row. -/
theorem row3_GL2_lift
    (α β γ δ : k) :
    PaperCore.ActsOnOrderedPair
      (N4PureSingular.ActBivector (k := k))
      (N4PureSingular.ω12 (k := k))
      (N4PureSingular.ω13 (k := k))
      (row3_GL2Lift (k := k) α β γ δ)
      (row3_GL2Coeff (k := k) α β γ δ) := by
  simpa [PaperCore.ActsOnOrderedPair, row3_GL2Lift, row3_GL2Coeff] using
    row3_GL2_lift_action (k := k) (α := α) (β := β) (γ := γ) (δ := δ)

end N4PaperSummary
end Wedge2Formalization
