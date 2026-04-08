import Wedge2Formalization.N7
import Wedge2Formalization.N7OnePoint
import Wedge2Formalization.N7PureSingular
import Wedge2Formalization.N7ThreePoint
import Wedge2Formalization.N7WeightedTwoPoint
import Wedge2Formalization.N7MixedOnePoint

namespace Wedge2Formalization
namespace N7Summary

open Matrix

variable {k : Type*} [Field k]

/-- Summary form of the first `n = 7` pointwise stabilizer calculation. This is the
three-dimensional radical extension of the split-support `n = 4` orbit. -/
theorem rad3Split_pointwise_bivector_iff
    (A0 : Matrix N7.I N7.I k)
    (R : Matrix N7.I N7.W k)
    (C : Matrix N7.W N7.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N7.FixesRad3SplitPairBivector
      (Matrix.fromBlocks A0 R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      R = 0 ∧ A.det = 1 ∧ B₁ = 0 ∧ C₁ = 0 ∧ D.det = 1 := by
  rw [N7.fixesRad3SplitPair_fromBlocks_iff]
  constructor
  · rintro ⟨hR, hD⟩
    rcases
      (N4Summary.split_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B₁)
        (C := C₁)
        (D := D)).1 (by simpa using hD) with
      ⟨hA, hB, hC, hDdet⟩
    exact ⟨hR, hA, hB, hC, hDdet⟩
  · rintro ⟨hR, hA, hB, hC, hDdet⟩
    refine ⟨hR, ?_⟩
    exact
      (N4Summary.split_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B₁)
        (C := C₁)
        (D := D)).2 ⟨hA, hB, hC, hDdet⟩

/-- The three-dimensional radical extension family preserves the split-support `2`-space
setwise. -/
theorem rad3Split_diag_preserves
    (A0 : Matrix N7.I N7.I k)
    (C : Matrix N7.W N7.I k)
    (A D : Matrix N4.I N4.I k) :
    N7.PreservesRad3SplitSubspaceBivector
      (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D)) :=
  N7.rad3Split_setwise_family (k := k) (A0 := A0) (C := C) (A := A) (D := D)

/-- The standard torus lift on the quotient extends to the three-dimensional radical
extension. -/
theorem rad3Split_torus_lift_action
    (A0 : Matrix N7.I N7.I k)
    (u v : k)
    (C : Matrix N7.W N7.I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix N7.W N7.W k := Matrix.fromBlocks A 0 0 D
    let g : Matrix N7.V N7.V k := Matrix.fromBlocks A0 0 C h
    N7.ActBivector N7.rad3SplitRep₁ g = u • N7.rad3SplitRep₁ + (v - u) • N7.rad3SplitRep₂ ∧
      N7.ActBivector N7.rad3SplitRep₂ g = v • N7.rad3SplitRep₂ :=
  N7.rad3Split_torus_lift_action (k := k) (A0 := A0) (u := u) (v := v) (C := C)

/-- The swap coset in the split normalizer also extends to the three-dimensional radical
extension. -/
theorem rad3Split_swapCoset_lift_action
    (A0 : Matrix N7.I N7.I k)
    (u v : k)
    (C : Matrix N7.W N7.I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix N7.W N7.W k := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k)
    let g : Matrix N7.V N7.V k := Matrix.fromBlocks A0 0 C h
    N7.ActBivector N7.rad3SplitRep₁ g = u • N7.rad3SplitRep₁ + (v - u) • N7.rad3SplitRep₂ ∧
      N7.ActBivector N7.rad3SplitRep₂ g = u • N7.rad3SplitRep₁ + (-u) • N7.rad3SplitRep₂ :=
  N7.rad3Split_swapCoset_lift_action (k := k) (A0 := A0) (u := u) (v := v) (C := C)

/-- Summary form of the second `n = 7` pointwise stabilizer calculation. This is the
three-dimensional radical extension of the repeated-support `n = 4` orbit. -/
theorem rad3OnePoint_pointwise_bivector_iff
    (A0 : Matrix N7.I N7.I k)
    (R : Matrix N7.I N7.W k)
    (C : Matrix N7.W N7.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N7OnePoint.FixesRad3OnePointPairBivector
      (Matrix.fromBlocks A0 R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      R = 0 ∧ A.det = 1 ∧ C₁ = 0 ∧ D = A ∧ A * N4.J * B₁ᵀ + B₁ * N4.J * Aᵀ = 0 := by
  rw [N7OnePoint.fixesRad3OnePointPair_fromBlocks_iff]
  constructor
  · rintro ⟨hR, hD⟩
    rcases
      (N4Summary.onePoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B₁)
        (C := C₁)
        (D := D)).1 (by simpa using hD) with
      ⟨hA, hC, hD', hrel⟩
    exact ⟨hR, hA, hC, hD', hrel⟩
  · rintro ⟨hR, hA, hC, hD, hrel⟩
    refine ⟨hR, ?_⟩
    exact
      (N4Summary.onePoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B₁)
        (C := C₁)
        (D := D)).2 ⟨hA, hC, hD, hrel⟩

/-- The standard upper-triangular quotient family on the repeated-support orbit extends
to the three-dimensional radical extension. -/
theorem rad3OnePoint_borel_lift_action
    (A0 : Matrix N7.I N7.I k)
    (a b : k)
    (ha : a ≠ 0)
    (C : Matrix N7.W N7.I k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let g : Matrix N7.V N7.V k :=
      Matrix.fromBlocks
        A0
        0
        C
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
    N7.ActBivector N7OnePoint.rad3OnePointRep₁ g =
        a • N7OnePoint.rad3OnePointRep₁ + b • N7OnePoint.rad3OnePointRep₂ ∧
      N7.ActBivector N7OnePoint.rad3OnePointRep₂ g =
        (a * a) • N7OnePoint.rad3OnePointRep₂ :=
  N7OnePoint.rad3OnePoint_borel_lift_action
    (k := k) (A0 := A0) (a := a) (b := b) ha C

/-- Summary form of the third `n = 7` pointwise stabilizer calculation. This is the
three-dimensional radical extension of the pure singular `n = 4` orbit. -/
theorem rad3PureSingular_pointwise_bivector_iff
    (A0 : Matrix N7PureSingular.I N7PureSingular.I k)
    (R : Matrix N7PureSingular.I N7PureSingular.W k)
    (C : Matrix N7PureSingular.W N7PureSingular.I k)
    (D : Matrix N4PureSingular.I N4PureSingular.I k) :
    N7PureSingular.FixesRad3PureSingularPairBivector
      (Matrix.fromBlocks A0 R C D) ↔
      R = N7PureSingular.upperRightLast (k := k) (R 0 3) (R 1 3) (R 2 3) ∧
      D 0 0 ≠ 0 ∧
      D =
        N4PureSingular.pureSingularShape
          (k := k) (D 0 0) (D 0 1) (D 0 2) (D 0 3) (D 1 3) (D 2 3) (D 3 3) := by
  exact
    N7PureSingular.fixesRad3PureSingularPair_fromBlocks_iff
      (k := k) (A := A0) (B := R) (C := C) (D := D)

/-- The full `GL₂` quotient lift on the pure singular orbit extends to the three-dimensional
radical extension. -/
theorem rad3PureSingular_GL2_lift_action
    (α β γ δ : k)
    (A0 : Matrix N7PureSingular.I N7PureSingular.I k)
    (C : Matrix N7PureSingular.W N7PureSingular.I k) :
    let h :=
      N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1
    let g : Matrix N7PureSingular.V N7PureSingular.V k := Matrix.fromBlocks A0 0 C h
    N7PureSingular.ActBivector N7PureSingular.rad3PureSingularRep₁ g =
        α • N7PureSingular.rad3PureSingularRep₁ + β • N7PureSingular.rad3PureSingularRep₂ ∧
      N7PureSingular.ActBivector N7PureSingular.rad3PureSingularRep₂ g =
        γ • N7PureSingular.rad3PureSingularRep₁ + δ • N7PureSingular.rad3PureSingularRep₂ :=
  N7PureSingular.rad3PureSingular_GL2_lift_action
    (k := k) (α := α) (β := β) (γ := γ) (δ := δ) (A0 := A0) C

/-- Summary form of the weighted two-point `n = 7` pointwise stabilizer calculation. -/
theorem radWeighted_pointwise_bivector_iff
    (a : Matrix N7WeightedTwoPoint.I N7WeightedTwoPoint.I k)
    (R : Matrix N7WeightedTwoPoint.I N7WeightedTwoPoint.W k)
    (C : Matrix N7WeightedTwoPoint.W N7WeightedTwoPoint.I k)
    (A : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k)
    (B : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.I k)
    (D : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k) :
    N7WeightedTwoPoint.FixesRadWeightedPairBivector
      (Matrix.fromBlocks a R C (Matrix.fromBlocks A B 0 D)) ↔
      R = 0 ∧ B = 0 ∧ N4.FixesBivector (N4.splitRep₁ (k := k)) A ∧ D.det = 1 := by
  rw [N7WeightedTwoPoint.fixesRadWeightedPair_fromBlocks_iff]
  constructor
  · rintro ⟨hR, hD⟩
    rcases
      (N6Summary.weightedTwoPoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B)
        (C := 0)
        (D := D)).1 (by simpa using hD) with
      ⟨hB, hC, hA, hDdet⟩
    exact ⟨hR, hB, hA, hDdet⟩
  · rintro ⟨hR, hB, hA, hDdet⟩
    refine ⟨hR, ?_⟩
    exact
      (N6Summary.weightedTwoPoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B)
        (C := 0)
        (D := D)).2 ⟨hB, rfl, hA, hDdet⟩

/-- The diagonal torus lift on the weighted two-point orbit extends to the
one-dimensional radical extension. -/
theorem radWeighted_torus_lift_action
    (a0 u v : k)
    (C : Matrix N7WeightedTwoPoint.W N7WeightedTwoPoint.I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k := Matrix.fromBlocks A 0 0 A
    let E : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
    let h : Matrix N7WeightedTwoPoint.W N7WeightedTwoPoint.W k := Matrix.fromBlocks H 0 0 E
    let g : Matrix N7WeightedTwoPoint.V N7WeightedTwoPoint.V k :=
      Matrix.fromBlocks (N7WeightedTwoPoint.scalarBlock (k := k) a0) 0 C h
    N7WeightedTwoPoint.ActBivector N7WeightedTwoPoint.radWeightedRep₁ g =
        u • N7WeightedTwoPoint.radWeightedRep₁ + (v - u) • N7WeightedTwoPoint.radWeightedRep₂ ∧
    N7WeightedTwoPoint.ActBivector N7WeightedTwoPoint.radWeightedRep₂ g =
        v • N7WeightedTwoPoint.radWeightedRep₂ :=
  N7WeightedTwoPoint.radWeighted_torus_lift_action (k := k) (a0 := a0) (u := u) (v := v) C

/-- On the embedded weighted two-point line in dimension `7`, the lower-right
determinant vanishes exactly on the two distinguished projective points. -/
theorem radWeighted_lowerRight_det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂
          (a • (N7WeightedTwoPoint.radWeightedRep₁ (k := k)) +
            b • (N7WeightedTwoPoint.radWeightedRep₂ (k := k)))) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  N7WeightedTwoPoint.lowerRight_det_zero_iff (k := k) (a := a) (b := b)

/-- On the embedded mixed one-point line in dimension `7`, the lower-right determinant
vanishes exactly at the repeated-support point. -/
theorem radMixed_lowerRight_det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂
          (a • (N7MixedOnePoint.radMixedRep₁ (k := k)) +
            b • (N7MixedOnePoint.radMixedRep₂ (k := k)))) = (0 : k) ↔
      a = 0 :=
  N7MixedOnePoint.lowerRight_det_zero_iff (k := k) (a := a) (b := b)

/-- A one-dimensional radical extension of the mixed `SL₂ × SL₂` pointwise family. -/
theorem radMixed_pointwise_levi_family
    (a0 : Matrix N7MixedOnePoint.I N7MixedOnePoint.I k)
    (C : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k)
    (A B H : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    N7MixedOnePoint.FixesRadMixedPairBivector
      (Matrix.fromBlocks a0 0 C (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) :=
  N7MixedOnePoint.radMixed_pointwise_levi_family
    (k := k) (a0 := a0) (C := C) (A := A) (B := B) (H := H) hA hB hH

/-- Summary form of the mixed one-point `n = 7` pointwise stabilizer calculation. -/
theorem radMixed_pointwise_bivector_iff
    (a0 : Matrix N7MixedOnePoint.I N7MixedOnePoint.I k)
    (R : Matrix N7MixedOnePoint.I N7MixedOnePoint.W k)
    (C : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k)
    (D : Matrix N7MixedOnePoint.W N7MixedOnePoint.W k) :
    N7MixedOnePoint.FixesRadMixedPairBivector
      (Matrix.fromBlocks a0 R C D) ↔
      R = 0 ∧ N6MixedOnePoint.FixesPairBivector D :=
  N7MixedOnePoint.fixesRadMixedPair_fromBlocks_iff
    (k := k) (a := a0) (B := R) (C := C) (D := D)

/-- A one-dimensional radical extension of the explicit coupled `E_{12}` family on the
mixed one-point orbit. -/
theorem radMixed_coupledE12_pointwise
    (a0 : Matrix N7MixedOnePoint.I N7MixedOnePoint.I k)
    (C : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k)
    (t : k) :
    N7MixedOnePoint.FixesRadMixedPairBivector
      (Matrix.fromBlocks
        a0
        0
        C
        (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ12 (k := k) t)
          (N6MixedOnePoint.coupledR12 (k := k) t)
          1)) :=
  N7MixedOnePoint.radMixed_coupledE12_pointwise (k := k) (a0 := a0) (C := C) t

/-- A one-dimensional radical extension of the explicit coupled `E_{21}` family on the
mixed one-point orbit. -/
theorem radMixed_coupledE21_pointwise
    (a0 : Matrix N7MixedOnePoint.I N7MixedOnePoint.I k)
    (C : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k)
    (t : k) :
    N7MixedOnePoint.FixesRadMixedPairBivector
      (Matrix.fromBlocks
        a0
        0
        C
        (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          1)) :=
  N7MixedOnePoint.radMixed_coupledE21_pointwise (k := k) (a0 := a0) (C := C) t

/-- The product of the two explicit mixed one-point coupled families still fixes the
embedded pair pointwise. -/
theorem radMixed_coupled_product_pointwise
    (s t : k) :
    N7MixedOnePoint.FixesRadMixedPairBivector
      ((Matrix.fromBlocks
          (1 : Matrix N7MixedOnePoint.I N7MixedOnePoint.I k)
          0
          (0 : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k)
          (Matrix.fromBlocks
            1
            (N6MixedOnePoint.coupledQ12 (k := k) s)
            (N6MixedOnePoint.coupledR12 (k := k) s)
            1)) *
       (Matrix.fromBlocks
          (1 : Matrix N7MixedOnePoint.I N7MixedOnePoint.I k)
          0
          (0 : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k)
          (Matrix.fromBlocks
            1
            (N6MixedOnePoint.coupledQ21 (k := k) t)
            (N6MixedOnePoint.coupledR21 (k := k) t)
            1))) :=
  N7MixedOnePoint.radMixed_coupledE12E21_product_pointwise (k := k) s t

/-- The standard upper-triangular quotient family on the mixed one-point orbit extends
to the one-dimensional radical extension. -/
theorem radMixed_borel_lift_action
    (a0 a b : k)
    (ha : a ≠ 0)
    (C : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N7MixedOnePoint.W N7MixedOnePoint.W k := Matrix.fromBlocks
      (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      0 0
      (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g : Matrix N7MixedOnePoint.V N7MixedOnePoint.V k :=
      Matrix.fromBlocks (N7MixedOnePoint.scalarBlock (k := k) a0) 0 C H
    N7MixedOnePoint.ActBivector N7MixedOnePoint.radMixedRep₁ g =
        a • N7MixedOnePoint.radMixedRep₁ + b • N7MixedOnePoint.radMixedRep₂ ∧
      N7MixedOnePoint.ActBivector N7MixedOnePoint.radMixedRep₂ g =
        (a * a) • N7MixedOnePoint.radMixedRep₂ :=
  N7MixedOnePoint.radMixed_borel_lift_action
    (k := k) (a0 := a0) (a := a) (b := b) ha C

/-- The determinant of the standard mixed one-point lift on the `n = 7` radical
extension is the radical scalar times `a^3`. -/
theorem radMixed_borel_lift_det
    (a0 a b : k)
    (C : Matrix N7MixedOnePoint.W N7MixedOnePoint.I k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N7MixedOnePoint.W N7MixedOnePoint.W k := Matrix.fromBlocks
      (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      0 0
      (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g : Matrix N7MixedOnePoint.V N7MixedOnePoint.V k :=
      Matrix.fromBlocks (N7MixedOnePoint.scalarBlock (k := k) a0) 0 C H
    Matrix.det g = a0 * (a * a * a) :=
  N7MixedOnePoint.radMixed_borel_lift_det (k := k) (a0 := a0) (a := a) (b := b) C

/-- Summary form of the three-point `n = 7` pointwise stabilizer calculation. -/
theorem radThreePoint_pointwise_bivector_iff
    (a : Matrix N7ThreePoint.I N7ThreePoint.I k)
    (R : Matrix N7ThreePoint.I N7ThreePoint.W k)
    (C0 : Matrix N7ThreePoint.W N7ThreePoint.I k)
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (B : Matrix N6ThreePoint.W N6ThreePoint.I k)
    (C : Matrix N6ThreePoint.I N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k) :
    N7ThreePoint.FixesRadThreePointPairBivector
      (Matrix.fromBlocks a R C0 (Matrix.fromBlocks A B C D)) ↔
      R = 0 ∧ B = 0 ∧ C = 0 ∧ N4.FixesSplitPairBivector A ∧ D.det = 1 :=
  N7ThreePoint.fixesRadThreePointPair_fromBlocks_iff
    (k := k) (a := a) (R := R) (C0 := C0) (A := A) (B := B) (C := C) (D := D)

/-- On the embedded three-point line in dimension `7`, the lower-right determinant
vanishes exactly on the three distinguished projective points. -/
theorem radThreePoint_lowerRight_det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂
          (a • (N7ThreePoint.radThreePointRep₁ (k := k)) +
            b • (N7ThreePoint.radThreePointRep₂ (k := k)))) = (0 : k) ↔
      a = 0 ∨ a + b = 0 ∨ b = 0 :=
  N7ThreePoint.lowerRight_det_zero_iff (k := k) (a := a) (b := b)

end N7Summary
end Wedge2Formalization
