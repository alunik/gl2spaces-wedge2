import Wedge2Formalization.N6Summary

namespace Wedge2Formalization
namespace N6PaperSummary

variable {k : Type*} [Field k]

/-!
# Paper-Facing `n = 6` Summary

This file repackages the finished geometric `n = 6` results in the row order of
the paper tables.

The goal is citation convenience. For each paper row we expose a small set of theorem
aliases:

- a pointwise theorem describing the concrete `K_L` family in Lean coordinates;
- one or more quotient-action theorems describing the displayed `Q_L`;
- when useful, a rank-drop theorem matching the divisor picture from the paper.

All row numbers below match the numbering in
`paper_wedge2_rewrite/sections/geometric_orbit_tables_generated.tex` and
`paper_wedge2_rewrite/sections/geometric_stabilizer_tables_generated.tex`.
-/

/-! ## Row 1

Representative
`⟨e₁ ∧ e₄ + e₂ ∧ e₅ + e₃ ∧ e₆, e₁ ∧ e₅ + e₂ ∧ e₆⟩`,
divisor type `3[a]`,
paper stabilizer `K_L = U₆ ⋊ SL₂(k)`, `Q_L = B`.

Lean model: `N6OnePointLong`.
-/

abbrev row1_pointwise_cell_iff := N6Summary.onePointLong_productCell_iff (k := k)
abbrev row1_pointwise_family :=
  N6Summary.onePointLong_scaleUpperLowerHankel_product_pointwise (k := k)
abbrev row1_borel_lift_action := N6Summary.onePointLong_borel_lift_action (k := k)

/-! ## Row 2

Representative
`⟨e₁ ∧ e₃ + e₂ ∧ e₄ + e₅ ∧ e₆, e₁ ∧ e₄⟩`,
divisor type `3[a]`,
paper stabilizer `K_L = U₇ ⋊ (SL₂(k) × SL₂(k))`, `Q_L = B`.

Lean model: `N6MixedOnePoint`.
-/

abbrev row2_unipotent_cell_iff :=
  N6Summary.mixedOnePoint_Block3_unipotent_iff_exists_centralBlock_mul_coupledKernelFrom (k := k)
abbrev row2_coupled_kernel_family := N6Summary.mixedOnePoint_coupledKernelFrom_pointwise (k := k)
abbrev row2_pointwise_levi_family := N6Summary.mixedOnePoint_pointwise_levi_family (k := k)
abbrev row2_borel_lift_action := N6Summary.mixedOnePoint_borel_lift_action (k := k)

/-! ## Row 3

Representative
`⟨e₁ ∧ e₃ + e₂ ∧ e₄ + e₅ ∧ e₆, e₁ ∧ e₄ + e₅ ∧ e₆⟩`,
divisor type `2[a]+[b]`,
paper stabilizer `K_L = U₃ ⋊ (SL₂(k) × SL₂(k))`, `Q_L = T`.

Lean model: `N6OnePointPlusSimple`.
-/

abbrev row3_pointwise_iff := N6Summary.onePointPlusSimple_pointwise_bivector_iff (k := k)
abbrev row3_pointwise_family := N6Summary.onePointPlusSimple_pointwise_levi_family (k := k)
abbrev row3_torus_lift_action := N6Summary.onePointPlusSimple_torus_lift_action (k := k)

/-! ## Row 4

Representative
`⟨e₁ ∧ e₂ + e₃ ∧ e₄ + e₅ ∧ e₆, e₅ ∧ e₆⟩`,
divisor type `2[a]+[b]`,
paper stabilizer `K_L = Sp₄(k) × SL₂(k)`, `Q_L = T`.

Lean model: `N6WeightedTwoPoint`.
-/

abbrev row4_pointwise_iff := N6Summary.weightedTwoPoint_pointwise_bivector_iff (k := k)
abbrev row4_pointwise_family := N6Summary.weightedTwoPoint_pointwise_family (k := k)
abbrev row4_torus_lift_action := N6Summary.weightedTwoPoint_torus_lift_action (k := k)

/-! ## Row 5

Representative
`⟨e₁ ∧ e₂ + e₃ ∧ e₄, e₃ ∧ e₄ + e₅ ∧ e₆⟩`,
divisor type `[a]+[b]+[c]`,
paper stabilizer `K_L = SL₂(k) × SL₂(k) × SL₂(k)`, `Q_L =` preimage of `S₃`.

Lean model: `N6ThreePoint`.
-/

abbrev row5_pointwise_iff := N6Summary.threePoint_pointwise_bivector_iff (k := k)
abbrev row5_pointwise_family := N6Summary.threePoint_pointwise_family (k := k)
abbrev row5_swap_lift_action := N6Summary.threePoint_swap_lift_action (k := k)
abbrev row5_det_zero_iff := N6Summary.threePoint_det_zero_iff (k := k)

/-! ## Row 6

Representative
`⟨e₃ ∧ e₅ + e₄ ∧ e₆, e₃ ∧ e₆⟩`,
divisor type `2[a]`,
paper stabilizer `K_L = U₁₁ ⋊ (GL₂(k) × SL₂(k))`, `Q_L = B`.

Lean model: `N6OnePoint`.
-/

abbrev row6_pointwise_iff := N6Summary.rad2OnePoint_pointwise_bivector_iff (k := k)
abbrev row6_borel_lift_action := N6Summary.rad2OnePoint_borel_lift_action (k := k)

/-! ## Row 7

Representative
`⟨e₃ ∧ e₄ + e₅ ∧ e₆, e₅ ∧ e₆⟩`,
divisor type `[a]+[b]`,
paper stabilizer `K_L = U₈ ⋊ (GL₂(k) × SL₂(k) × SL₂(k))`, `Q_L = N_T`.

Lean model: `N6`.
-/

abbrev row7_pointwise_iff := N6Summary.rad2Split_pointwise_bivector_iff (k := k)
abbrev row7_torus_lift_action := N6Summary.rad2Split_torus_lift_action (k := k)
abbrev row7_swap_lift_action := N6Summary.rad2Split_swapCoset_lift_action (k := k)

/-! ## Row 8

Representative
`⟨e₂ ∧ e₄ + e₅ ∧ e₆, e₃ ∧ e₄⟩`,
divisor type `[a]`,
paper stabilizer `K_L = U₉ ⋊ (G_m × G_m × SL₂(k))`, `Q_L = B`.

Lean model: `N6SimplePoint`.
-/

abbrev row8_pointwise_iff := N6Summary.simplePoint_nested_iff_shape (k := k)
abbrev row8_borel_lift_action := N6Summary.simplePoint_borel_lift_action (k := k)

/-! ## Row 9

Representative
`⟨e₄ ∧ e₆, e₅ ∧ e₆⟩`,
divisor type `0`,
paper stabilizer `K_L = U₁₁ ⋊ (GL₃(k) × G_m)`, `Q_L = GL₂(k)`.

Lean model: `N6PureSingular3`.
-/

abbrev row9_pointwise_iff := N6Summary.pureSingular3_pointwise_bivector_iff (k := k)
abbrev row9_GL2_lift_action := N6Summary.pureSingular3_GL2_lift_action (k := k)

/-! ## Row 10

Representative
`⟨e₂ ∧ e₅ + e₃ ∧ e₆, e₃ ∧ e₅ + e₄ ∧ e₆⟩`,
divisor type `0`,
paper stabilizer `K_L = U₈ ⋊ (G_m × G_m)`, `Q_L = GL₂(k)`.

Lean model: `N6PureSingularLong`.
-/

abbrev row10_pointwise_iff := N6Summary.pureSingularLong_nested_iff_shape (k := k)
abbrev row10_pointwise_family :=
  N6Summary.pureSingularLong_scaleLowerHankel_product_pointwise (k := k)
abbrev row10_GL2_lift_action := N6Summary.pureSingularLong_GL2_lift_action (k := k)

/-! ## Row 11

Representative
`⟨e₁ ∧ e₃ + e₄ ∧ e₆, e₂ ∧ e₃ + e₅ ∧ e₆⟩`,
divisor type `0`,
paper stabilizer `K_L = U₆ ⋊ GL₂(k)`, `Q_L = GL₂(k)`.

Lean model: `N6DoublePureSingular`.
-/

/-- Paper-facing name for the `U₆` family on the direct-sum pure singular row. -/
abbrev row11_U6 := N6DoublePureSingular.coupledUnipotent (k := k)

/-- Paper-facing name for the pointwise `GL₂(k)` family on the direct-sum pure singular row. -/
abbrev row11_pointwiseGL2 := N6DoublePureSingular.magmaPointwiseGL2 (k := k)

/-- Paper-facing name for the full `U₆ ⋊ GL₂(k)` kernel cell on the direct-sum pure singular row. -/
abbrev row11_kernelCell := N6DoublePureSingular.magmaKernelCell (k := k)

/-- A concrete `U₆ ⋊ (G_m × G_m)` pointwise family on the direct-sum pure singular row. -/
abbrev row11_shapeTimesU6_pointwise :=
  N6Summary.doublePureSingular_coupled_shape_product_pointwise (k := k)

/-- Exact pointwise criterion after right-multiplying by the paper-facing `U₆` family. -/
abbrev row11_blockDiagonalTimesU6_iff :=
  N6Summary.doublePureSingular_blockDiagonal_coupled_right_product_iff_shape (k := k)

abbrev row11_pointwise_iff := N6Summary.doublePureSingular_blockDiagonal_iff_shape (k := k)
abbrev row11_magmaPointwiseGL2_family :=
  N6Summary.doublePureSingular_magmaPointwiseGL2_family (k := k)
abbrev row11_magmaKernelCell_pointwise :=
  N6Summary.doublePureSingular_magmaKernelCell_pointwise (k := k)

/-- The quotient `GL₂(k)` action respects the paper-facing `U₆` family on the left. -/
abbrev row11_U6_GL2_product_action :=
  N6Summary.doublePureSingular_coupled_GL2_product_lift_action (k := k)

/-- Right-multiplying the quotient `GL₂(k)` lift by `U₆` does not change the quotient action. -/
abbrev row11_GL2_U6_right_action :=
  N6Summary.doublePureSingular_GL2_coupled_right_product_lift_action (k := k)

abbrev row11_GL2_lift_action := N6Summary.doublePureSingular_GL2_lift_action (k := k)

end N6PaperSummary
end Wedge2Formalization
