import Wedge2Formalization.N6
import Wedge2Formalization.N6OnePoint
import Wedge2Formalization.N6OnePointPlusSimple
import Wedge2Formalization.N6DoublePureSingular
import Wedge2Formalization.N6PureSingular
import Wedge2Formalization.N6PureSingular3
import Wedge2Formalization.N6PureSingularLong
import Wedge2Formalization.N6SimplePoint
import Wedge2Formalization.N6ThreePoint
import Wedge2Formalization.N6WeightedTwoPoint
import Wedge2Formalization.N6MixedOnePoint
import Wedge2Formalization.N6OnePointLong

namespace Wedge2Formalization
namespace N6Summary

open Matrix

variable {k : Type*} [Field k]

/-- Summary form of the first `n = 6` pointwise stabilizer calculation. This is the
two-dimensional radical extension of the split-support `n = 4` orbit. -/
theorem rad2Split_pointwise_bivector_iff
    (A0 : Matrix N6.I N6.I k)
    (R : Matrix N6.I N6.W k)
    (C : Matrix N6.W N6.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N6.FixesRad2SplitPairBivector
      (Matrix.fromBlocks A0 R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      R = 0 ∧ A.det = 1 ∧ B₁ = 0 ∧ C₁ = 0 ∧ D.det = 1 := by
  rw [N6.fixesRad2SplitPair_fromBlocks_iff]
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

/-- The two-dimensional radical extension family preserves the split-support `2`-space
setwise. -/
theorem rad2Split_diag_preserves
    (A0 : Matrix N6.I N6.I k)
    (C : Matrix N6.W N6.I k)
    (A D : Matrix N4.I N4.I k) :
    N6.PreservesRad2SplitSubspaceBivector
      (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D)) :=
  N6.rad2Split_setwise_family (k := k) (A0 := A0) (C := C) (A := A) (D := D)

/-- The standard torus lift on the quotient extends to the two-dimensional radical
extension. -/
theorem rad2Split_torus_lift_action
    (A0 : Matrix N6.I N6.I k)
    (u v : k)
    (C : Matrix N6.W N6.I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix N6.W N6.W k := Matrix.fromBlocks A 0 0 D
    let g : Matrix N6.V N6.V k := Matrix.fromBlocks A0 0 C h
    N6.ActBivector N6.rad2SplitRep₁ g = u • N6.rad2SplitRep₁ + (v - u) • N6.rad2SplitRep₂ ∧
      N6.ActBivector N6.rad2SplitRep₂ g = v • N6.rad2SplitRep₂ :=
  N6.rad2Split_torus_lift_action (k := k) (A0 := A0) (u := u) (v := v) (C := C)

/-- The swap coset in the split normalizer also extends to the two-dimensional radical
extension. -/
theorem rad2Split_swapCoset_lift_action
    (A0 : Matrix N6.I N6.I k)
    (u v : k)
    (C : Matrix N6.W N6.I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix N6.W N6.W k := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k)
    let g : Matrix N6.V N6.V k := Matrix.fromBlocks A0 0 C h
    N6.ActBivector N6.rad2SplitRep₁ g = u • N6.rad2SplitRep₁ + (v - u) • N6.rad2SplitRep₂ ∧
      N6.ActBivector N6.rad2SplitRep₂ g = u • N6.rad2SplitRep₁ + (-u) • N6.rad2SplitRep₂ :=
  N6.rad2Split_swapCoset_lift_action (k := k) (A0 := A0) (u := u) (v := v) (C := C)

/-- Summary form of the second `n = 6` pointwise stabilizer calculation. This is the
two-dimensional radical extension of the repeated-support `n = 4` orbit. -/
theorem rad2OnePoint_pointwise_bivector_iff
    (A0 : Matrix N6.I N6.I k)
    (R : Matrix N6.I N6.W k)
    (C : Matrix N6.W N6.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N6OnePoint.FixesRad2OnePointPairBivector
      (Matrix.fromBlocks A0 R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      R = 0 ∧ A.det = 1 ∧ C₁ = 0 ∧ D = A ∧ A * N4.J * B₁ᵀ + B₁ * N4.J * Aᵀ = 0 := by
  rw [N6OnePoint.fixesRad2OnePointPair_fromBlocks_iff]
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
to the two-dimensional radical extension. -/
theorem rad2OnePoint_borel_lift_action
    (A0 : Matrix N6.I N6.I k)
    (a b : k)
    (ha : a ≠ 0)
    (C : Matrix N6.W N6.I k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let g : Matrix N6.V N6.V k :=
      Matrix.fromBlocks
        A0
        0
        C
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
    N6.ActBivector N6OnePoint.rad2OnePointRep₁ g =
        a • N6OnePoint.rad2OnePointRep₁ + b • N6OnePoint.rad2OnePointRep₂ ∧
      N6.ActBivector N6OnePoint.rad2OnePointRep₂ g =
        (a * a) • N6OnePoint.rad2OnePointRep₂ :=
  N6OnePoint.rad2OnePoint_borel_lift_action
    (k := k) (A0 := A0) (a := a) (b := b) ha C

/-- The indecomposable length-three repeated-support orbit admits the explicit upper
Hankel pointwise family. -/
theorem onePointLong_upperHankel_pointwise
    (u v w : k) :
    N6OnePointLong.FixesPairBivector
      (Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        (N6OnePointLong.upperHankel (k := k) u v w)
        0
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) :=
  N6OnePointLong.upperHankel_pointwise (k := k) u v w

/-- On the upper unipotent cell of the indecomposable `3[a]` orbit, pointwise fixing
is equivalent to the upper Hankel shape. -/
theorem onePointLong_upperIdentity_iff
    (B : Matrix N6OnePointLong.I N6OnePointLong.I k) :
    N6OnePointLong.FixesPairBivector
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          B
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) ↔
      B =
        N6OnePointLong.upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) :=
  N6OnePointLong.fixesPair_upperIdentity_iff (k := k) B

/-- The indecomposable length-three repeated-support orbit also admits the explicit lower
Hankel pointwise family. -/
theorem onePointLong_lowerHankel_pointwise
    (u v w : k) :
    N6OnePointLong.FixesPairBivector
      (Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        0
        (N6OnePointLong.lowerHankel (k := k) u v w)
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) :=
  N6OnePointLong.lowerHankel_pointwise (k := k) u v w

/-- On the lower unipotent cell of the indecomposable `3[a]` orbit, pointwise fixing
is equivalent to the lower Hankel shape. -/
theorem onePointLong_lowerIdentity_iff
    (C : Matrix N6OnePointLong.I N6OnePointLong.I k) :
    N6OnePointLong.FixesPairBivector
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          C
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) ↔
      C =
        N6OnePointLong.lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) :=
  N6OnePointLong.fixesPair_lowerIdentity_iff (k := k) C

/-- Combining the upper and lower Hankel families gives an explicit six-parameter
pointwise unipotent subgroup on the indecomposable `3[a]` orbit. -/
theorem onePointLong_upperLowerHankel_product_pointwise
    (u₁ v₁ w₁ u₂ v₂ w₂ : k) :
    let gUpper : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
        0
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
    let gLower : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        0
        (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
    N6OnePointLong.FixesPairBivector (gUpper * gLower) :=
  N6OnePointLong.upperLowerHankel_product_pointwise
    (k := k) (u₁ := u₁) (v₁ := v₁) (w₁ := w₁) (u₂ := u₂) (v₂ := v₂) (w₂ := w₂)

/-- On the Magma local kernel product cell for the indecomposable `3[a]` orbit,
pointwise fixing is exactly the upper/lower Hankel condition on the two kernel
blocks. -/
theorem onePointLong_productCell_iff
    (B C : Matrix N6OnePointLong.I N6OnePointLong.I k) :
    N6OnePointLong.FixesPairBivector
      (Matrix.fromBlocks
        ((1 : Matrix N6OnePointLong.I N6OnePointLong.I k) + B * C)
        B
        C
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) ↔
      B = N6OnePointLong.upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) ∧
        C = N6OnePointLong.lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) :=
  N6OnePointLong.fixesPair_productCell_iff (k := k) B C

/-- Rank drop on the indecomposable `3[a]` line occurs exactly at the repeated-support
projective point. -/
theorem onePointLong_det_zero_iff
    (a b : k) :
    Matrix.det (a • (N6OnePointLong.rep₁ (k := k)) + b • (N6OnePointLong.rep₂ (k := k))) = 0 ↔
      a = 0 :=
  N6OnePointLong.det_zero_iff (k := k) (a := a) (b := b)

/-- The obvious torus on the indecomposable length-three repeated-support orbit fixes
the pair pointwise. -/
theorem onePointLong_scale_pointwise
    (a : k)
    (ha : a ≠ 0) :
    N6OnePointLong.FixesPairBivector
      (Matrix.fromBlocks
        (a • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
      (a⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) :=
  N6OnePointLong.scale_pointwise (k := k) a ha

/-- On the upper Hankel cell followed by the pointwise torus, pointwise fixing is still
equivalent to the upper Hankel shape. -/
theorem onePointLong_upperScale_iff
    (B : Matrix N6OnePointLong.I N6OnePointLong.I k)
    (a : k)
    (ha : a ≠ 0) :
    N6OnePointLong.FixesPairBivector
        ((Matrix.fromBlocks
            (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
            B
            0
            (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
         (Matrix.fromBlocks
            (a • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
            0
            0
            (a⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))) ↔
      B =
        N6OnePointLong.upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) :=
  N6OnePointLong.fixesPair_upperScale_iff (k := k) B a ha

/-- On the lower Hankel cell followed by the pointwise torus, pointwise fixing is still
equivalent to the lower Hankel shape. -/
theorem onePointLong_lowerScale_iff
    (C : Matrix N6OnePointLong.I N6OnePointLong.I k)
    (a : k)
    (ha : a ≠ 0) :
    N6OnePointLong.FixesPairBivector
        ((Matrix.fromBlocks
            (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
            0
            C
            (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
         (Matrix.fromBlocks
            (a • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
            0
            0
            (a⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))) ↔
      C =
        N6OnePointLong.lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) :=
  N6OnePointLong.fixesPair_lowerScale_iff (k := k) C a ha

/-- Combining the torus with the six-parameter unipotent family gives a concrete
pointwise subgroup on the indecomposable `3[a]` orbit. -/
theorem onePointLong_scaleUpperLowerHankel_product_pointwise
    (u₁ v₁ w₁ u₂ v₂ w₂ a : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      (Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
        0
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
      (Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        0
        (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
    let gScale : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (a • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (a⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
    N6OnePointLong.FixesPairBivector (gU * gScale) :=
  N6OnePointLong.scaleUpperLowerHankel_product_pointwise
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (a := a)
    ha

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit has determinant `1`. -/
theorem onePointLong_scaleUpperLowerHankel_product_det
    (u₁ v₁ w₁ u₂ v₂ w₂ a : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      (Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
        0
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
      (Matrix.fromBlocks
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
        0
        (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
        (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
    let gScale : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (a • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (a⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
    Matrix.det (gU * gScale) = 1 :=
  N6OnePointLong.scaleUpperLowerHankel_product_det
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (a := a)
    ha

/-- A concrete quotient shear lift on the indecomposable length-three repeated-support
orbit. -/
theorem onePointLong_shear_lift_action
    (b : k) :
    let g : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ g =
        N6OnePointLong.rep₁ + b • N6OnePointLong.rep₂ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ g =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.shear_lift_action (k := k) (b := b)

/-- A concrete quotient scaling lift on the indecomposable length-three repeated-support
orbit. -/
theorem onePointLong_scale_lift_action
    (a : k)
    (ha : a ≠ 0) :
    let g : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ g =
        a • N6OnePointLong.rep₁ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ g =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.scale_lift_action (k := k) (a := a) ha

/-- The quotient shear lift on the indecomposable `3[a]` orbit has determinant `1`. -/
theorem onePointLong_shear_lift_det
    (b : k) :
    let g : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    Matrix.det g = 1 :=
  N6OnePointLong.shear_lift_det (k := k) (b := b)

/-- The quotient scaling lift on the indecomposable `3[a]` orbit has determinant `a^3`. -/
theorem onePointLong_scale_lift_det
    (a : k)
    (ha : a ≠ 0) :
    let g : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    Matrix.det g = a * a * a :=
  N6OnePointLong.scale_lift_det (k := k) (a := a) ha

/-- A concrete upper-triangular quotient lift on the indecomposable length-three
repeated-support orbit. -/
theorem onePointLong_borel_lift_action
    (a b : k)
    (ha : a ≠ 0) :
    let gS : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ (gS * gU) =
        a • N6OnePointLong.rep₁ + b • N6OnePointLong.rep₂ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ (gS * gU) =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.borel_lift_action (k := k) (a := a) (b := b) ha

/-- The determinant of the standard Borel lift on the indecomposable `3[a]` orbit
is `a^3`. -/
theorem onePointLong_borel_lift_det
    (a b : k)
    (ha : a ≠ 0) :
    let gS : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    Matrix.det (gS * gU) = a * a * a :=
  N6OnePointLong.borel_lift_det (k := k) (a := a) (b := b) ha

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit combines with the
quotient shear lift with the expected action on the chosen basis. -/
theorem onePointLong_pointwise_shear_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ (gU * gL) =
        N6OnePointLong.rep₁ + b • N6OnePointLong.rep₂ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ (gU * gL) =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.pointwise_shear_product_lift_action
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (b := b)
    hc

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit also combines with the
quotient scaling lift with the expected action on the chosen basis. -/
theorem onePointLong_pointwise_scale_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ (gU * gL) =
        a • N6OnePointLong.rep₁ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ (gU * gL) =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.pointwise_scale_product_lift_action
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    hc
    ha

/-- Left-multiplying the quotient shear lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem onePointLong_pointwise_shear_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    Matrix.det (gU * gL) = 1 :=
  N6OnePointLong.pointwise_shear_product_lift_det
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (b := b)
    hc

/-- Left-multiplying the quotient scaling lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem onePointLong_pointwise_scale_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    Matrix.det (gU * gL) = a * a * a :=
  N6OnePointLong.pointwise_scale_product_lift_det
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    hc
    ha

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit combines with the
full upper-triangular quotient lift with the expected action on the chosen basis. -/
theorem onePointLong_pointwise_borel_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      (Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b))
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ (gU * gL) =
        a • N6OnePointLong.rep₁ + b • N6OnePointLong.rep₂ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ (gU * gL) =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.pointwise_borel_product_lift_action
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    (b := b)
    hc
    ha

/-- Left-multiplying the standard Borel lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem onePointLong_pointwise_borel_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      (Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b))
    Matrix.det (gU * gL) = a * a * a :=
  N6OnePointLong.pointwise_borel_product_lift_det
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    (b := b)
    hc
    ha

/-- Right-multiplying the quotient shear lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change the quotient action. -/
theorem onePointLong_shear_pointwise_right_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ (gL * gU) =
        N6OnePointLong.rep₁ + b • N6OnePointLong.rep₂ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ (gL * gU) =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.shear_pointwise_right_product_lift_action
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (b := b)
    hc

/-- Right-multiplying the quotient scaling lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change the quotient action. -/
theorem onePointLong_scale_pointwise_right_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ (gL * gU) =
        a • N6OnePointLong.rep₁ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ (gL * gU) =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.scale_pointwise_right_product_lift_action
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    hc
    ha

/-- Right-multiplying the quotient shear lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem onePointLong_shear_pointwise_right_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b)
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    Matrix.det (gL * gU) = 1 :=
  N6OnePointLong.shear_pointwise_right_product_lift_det
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (b := b)
    hc

/-- Right-multiplying the quotient scaling lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem onePointLong_scale_pointwise_right_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    Matrix.det (gL * gU) = a * a * a :=
  N6OnePointLong.scale_pointwise_right_product_lift_det
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    hc
    ha

/-- Right-multiplying the full upper-triangular quotient lift by the full pointwise
subgroup on the indecomposable `3[a]` orbit does not change the quotient action. -/
theorem onePointLong_borel_pointwise_right_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      (Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b))
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    N6OnePointLong.ActBivector N6OnePointLong.rep₁ (gL * gU) =
        a • N6OnePointLong.rep₁ + b • N6OnePointLong.rep₂ ∧
      N6OnePointLong.ActBivector N6OnePointLong.rep₂ (gL * gU) =
        N6OnePointLong.rep₂ :=
  N6OnePointLong.borel_pointwise_right_product_lift_action
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    (b := b)
    hc
    ha

/-- Right-multiplying the standard Borel lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem onePointLong_borel_pointwise_right_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      (Matrix.fromBlocks
        (N6OnePointLong.scaleTop (k := k) a)
        0
        0
        (N6OnePointLong.scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (N6OnePointLong.shearTop (k := k) b)
        0
        0
        (N6OnePointLong.shearBottom (k := k) b))
    let gU : Matrix N6OnePointLong.V N6OnePointLong.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          (N6OnePointLong.upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)
          0
          (N6OnePointLong.lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k))
        0
        0
        (c⁻¹ • (1 : Matrix N6OnePointLong.I N6OnePointLong.I k)))
    Matrix.det (gL * gU) = a * a * a :=
  N6OnePointLong.borel_pointwise_right_product_lift_det
    (k := k)
    (u₁ := u₁)
    (v₁ := v₁)
    (w₁ := w₁)
    (u₂ := u₂)
    (v₂ := v₂)
    (w₂ := w₂)
    (c := c)
    (a := a)
    (b := b)
    hc
    ha

/-- Summary form of the third `n = 6` pointwise stabilizer calculation. This is the
two-dimensional radical extension of the pure singular `n = 4` orbit. -/
theorem rad2PureSingular_pointwise_bivector_iff
    (A0 : Matrix N6PureSingular.I N6PureSingular.I k)
    (R : Matrix N6PureSingular.I N6PureSingular.W k)
    (C : Matrix N6PureSingular.W N6PureSingular.I k)
    (D : Matrix N4PureSingular.I N4PureSingular.I k) :
    N6PureSingular.FixesRad2PureSingularPairBivector
      (Matrix.fromBlocks A0 R C D) ↔
      R = N6PureSingular.upperRightLast (k := k) (R 0 3) (R 1 3) ∧
      D 0 0 ≠ 0 ∧
      D =
        N4PureSingular.pureSingularShape
          (k := k) (D 0 0) (D 0 1) (D 0 2) (D 0 3) (D 1 3) (D 2 3) (D 3 3) := by
  exact
    N6PureSingular.fixesRad2PureSingularPair_fromBlocks_iff
      (k := k) (A := A0) (B := R) (C := C) (D := D)

/-- The full `GL₂` quotient lift on the pure singular orbit extends to the two-dimensional
radical extension. -/
theorem rad2PureSingular_GL2_lift_action
    (α β γ δ : k)
    (A0 : Matrix N6PureSingular.I N6PureSingular.I k)
    (C : Matrix N6PureSingular.W N6PureSingular.I k) :
    let h :=
      N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1
    let g : Matrix N6PureSingular.V N6PureSingular.V k := Matrix.fromBlocks A0 0 C h
    N6PureSingular.ActBivector N6PureSingular.rad2PureSingularRep₁ g =
        α • N6PureSingular.rad2PureSingularRep₁ + β • N6PureSingular.rad2PureSingularRep₂ ∧
      N6PureSingular.ActBivector N6PureSingular.rad2PureSingularRep₂ g =
        γ • N6PureSingular.rad2PureSingularRep₁ + δ • N6PureSingular.rad2PureSingularRep₂ :=
  N6PureSingular.rad2PureSingular_GL2_lift_action
    (k := k) (α := α) (β := β) (γ := γ) (δ := δ) (A0 := A0) C

/-- The pure singular `3`-dimensional pair with a `3`-dimensional radical admits the
expected pointwise scaling family. -/
theorem pureSingular3_pointwise_scale_family
    (A0 : Matrix N6PureSingular3.I N6PureSingular3.I k)
    (a : k)
    (ha : a ≠ 0)
    (C : Matrix N6PureSingular3.W N6PureSingular3.I k) :
    N6PureSingular3.FixesPairBivector
      (Matrix.fromBlocks
        A0
        0
        C
        (N3PureSingular.pointwiseScale (k := k) a)) :=
  N6PureSingular3.pointwise_scale_family (k := k) (A0 := A0) (a := a) ha C

/-- The basic two-parameter unipotent family on the pure singular `3`-dimensional
quotient also extends through the `3`-dimensional radical. -/
theorem pureSingular3_pointwise_unipotent_family
    (A0 : Matrix N6PureSingular3.I N6PureSingular3.I k)
    (C : Matrix N6PureSingular3.W N6PureSingular3.I k)
    (x y : k) :
    N6PureSingular3.FixesPairBivector
      (Matrix.fromBlocks
        A0
        0
        C
        (N3PureSingular.pointwiseUnipotent (k := k) x y)) :=
  N6PureSingular3.pointwise_unipotent_family
    (k := k) (A0 := A0) (C := C) (x := x) (y := y)

/-- Summary form of the pure singular `3 + 3` pointwise stabilizer calculation. -/
theorem pureSingular3_pointwise_bivector_iff
    (A : Matrix N6PureSingular3.I N6PureSingular3.I k)
    (B : Matrix N6PureSingular3.I N6PureSingular3.W k)
    (C : Matrix N6PureSingular3.W N6PureSingular3.I k)
    (D : Matrix N6PureSingular3.W N6PureSingular3.W k) :
    N6PureSingular3.FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧
      D 0 0 ≠ 0 ∧
      D = N3PureSingular.pureSingularShape (k := k) (D 0 0) (D 2 0) (D 2 1) :=
  N6PureSingular3.fixesPair_fromBlocks_iff
    (k := k) (A := A) (B := B) (C := C) (D := D)

/-- The full `GL₂` quotient lift on the `3`-dimensional pure singular pair extends
through the `3`-dimensional radical. -/
theorem pureSingular3_GL2_lift_action
    (A0 : Matrix N6PureSingular3.I N6PureSingular3.I k)
    (α β γ δ : k)
    (C : Matrix N6PureSingular3.W N6PureSingular3.I k) :
    let h : Matrix N6PureSingular3.W N6PureSingular3.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix N6PureSingular3.V N6PureSingular3.V k := Matrix.fromBlocks A0 0 C h
    N6PureSingular3.ActBivector N6PureSingular3.rep₁ g =
        α • N6PureSingular3.rep₁ + β • N6PureSingular3.rep₂ ∧
      N6PureSingular3.ActBivector N6PureSingular3.rep₂ g =
        γ • N6PureSingular3.rep₁ + δ • N6PureSingular3.rep₂ :=
  N6PureSingular3.GL2_lift_action
    (k := k) (A0 := A0) (α := α) (β := β) (γ := γ) (δ := δ) C

/-- The determinant of the radical-extension `GL₂` lift on the pure singular
`3 + 3` row is the radical determinant times the `GL₂` determinant. -/
theorem pureSingular3_GL2_lift_det
    (A0 : Matrix N6PureSingular3.I N6PureSingular3.I k)
    (α β γ δ : k)
    (C : Matrix N6PureSingular3.W N6PureSingular3.I k) :
    let h : Matrix N6PureSingular3.W N6PureSingular3.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix N6PureSingular3.V N6PureSingular3.V k := Matrix.fromBlocks A0 0 C h
    Matrix.det g = A0.det * (α * δ - β * γ) :=
  N6PureSingular3.GL2_lift_det
    (k := k) (A0 := A0) (α := α) (β := β) (γ := γ) (δ := δ) C

/-- The radical extension of the pure singular length-two `n = 5` orbit admits the
same explicit scalar-plus-Hankel shape on the lower-right block whenever that block
already has the internal zero upper-right form. -/
theorem pureSingularLong_nested_zeroUpperRight_iff_shape
    (A0 : Matrix N6PureSingularLong.I N6PureSingularLong.I k)
    (C0 : Matrix N6PureSingularLong.W N6PureSingularLong.I k)
    (A1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
    (C1 : Matrix N5PureSingularLong.I N5PureSingularLong.W k)
    (D1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k) :
    N6PureSingularLong.FixesPairBivector
      (Matrix.fromBlocks
        A0
        0
        C0
        (Matrix.fromBlocks A1 0 C1 D1)) ↔
      A1 0 0 ≠ 0 ∧
        Matrix.fromBlocks A1 0 C1 D1 =
          N5PureSingularLong.pointwiseShape (k := k) (A1 0 0) (C1 0 0) (C1 0 1) (C1 0 2) (C1 1 2) :=
  N6PureSingularLong.fixesPair_fromBlocks_nested_zeroUpperRight_iff_shape
    (k := k) (A0 := A0) (C0 := C0) (A1 := A1) (C1 := C1) (D1 := D1)

/-- If the lower-right block is already written in the internal zero upper-right form
from the `n = 5` pure singular length-two model, then pointwise fixing forces the
radical upper-right row to vanish as well. -/
theorem pureSingularLong_nested_iff_shape
    (A0 : Matrix N6PureSingularLong.I N6PureSingularLong.I k)
    (B0 : Matrix N6PureSingularLong.I N6PureSingularLong.W k)
    (C0 : Matrix N6PureSingularLong.W N6PureSingularLong.I k)
    (A1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
    (C1 : Matrix N5PureSingularLong.I N5PureSingularLong.W k)
    (D1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k) :
    N6PureSingularLong.FixesPairBivector
      (Matrix.fromBlocks
        A0
        B0
        C0
        (Matrix.fromBlocks A1 0 C1 D1)) ↔
      B0 = 0 ∧
        A1 0 0 ≠ 0 ∧
        Matrix.fromBlocks A1 0 C1 D1 =
          N5PureSingularLong.pointwiseShape (k := k) (A1 0 0) (C1 0 0) (C1 0 1) (C1 0 2) (C1 1 2) :=
  N6PureSingularLong.fixesPair_fromBlocks_nested_iff_shape
    (k := k) (A0 := A0) (B0 := B0) (C0 := C0) (A1 := A1) (C1 := C1) (D1 := D1)

/-- The radical extension of the pure singular length-two `n = 5` orbit admits the
expected pointwise scaling family. -/
theorem pureSingularLong_pointwise_scale_family
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix N6PureSingularLong.W N6PureSingularLong.I k) :
    N6PureSingularLong.FixesPairBivector
      (Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)))) :=
  N6PureSingularLong.pointwise_scale_family (k := k) (u := u) (t := t) ht C

/-- The lower-left Hankel unipotent family on the pure singular length-two `n = 5`
orbit also extends through the radical direction. -/
theorem pureSingularLong_lowerHankel_pointwise_family
    (u : k)
    (C : Matrix N6PureSingularLong.W N6PureSingularLong.I k)
    (a b c d : k) :
    N6PureSingularLong.FixesPairBivector
      (Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) a b c d)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) :=
  N6PureSingularLong.lowerHankel_pointwise_family
    (k := k) (u := u) (C := C) (a := a) (b := b) (c := c) (d := d)

/-- Combining the radical-extension scalar torus and the lifted lower-left Hankel family
gives a concrete `6`-parameter pointwise subgroup on the pure singular length-two
`n = 6` orbit. -/
theorem pureSingularLong_scaleLowerHankel_product_pointwise
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix N6PureSingularLong.W N6PureSingularLong.I k)
    (a b c d : k) :
    let gScale : Matrix N6PureSingularLong.V N6PureSingularLong.V k :=
      Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)))
    let gH : Matrix N6PureSingularLong.V N6PureSingularLong.V k :=
      Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) 1)
        0
        0
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) a b c d)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))
    N6PureSingularLong.FixesPairBivector (gH * gScale) :=
  N6PureSingularLong.scaleLowerHankel_product_pointwise
    (k := k) (u := u) (t := t) ht (C := C) (a := a) (b := b) (c := c) (d := d)

/-- The determinant-scaled `GL₂` lift on the pure singular length-two `n = 5` orbit
extends through the radical direction. -/
theorem pureSingularLong_GL2_lift_action
    (u α β γ δ : k)
    (C : Matrix N6PureSingularLong.W N6PureSingularLong.I k) :
    let P : Matrix N5PureSingularLong.W N5PureSingularLong.W k :=
      N5PureSingularLong.sym2Raw (k := k) α β γ δ
    let Q : Matrix N5PureSingularLong.I N5PureSingularLong.I k :=
      N5PureSingularLong.adjointT (k := k) α β γ δ
    let h : Matrix N6PureSingularLong.W N6PureSingularLong.W k := Matrix.fromBlocks P 0 0 Q
    let g : Matrix N6PureSingularLong.V N6PureSingularLong.V k :=
      Matrix.fromBlocks (N6PureSingularLong.scalarBlock (k := k) u) 0 C h
    N6PureSingularLong.ActBivector N6PureSingularLong.rep₁ g =
        N5PureSingularLong.Delta α β γ δ •
          (α • N6PureSingularLong.rep₁ + γ • N6PureSingularLong.rep₂) ∧
      N6PureSingularLong.ActBivector N6PureSingularLong.rep₂ g =
        N5PureSingularLong.Delta α β γ δ •
          (β • N6PureSingularLong.rep₁ + δ • N6PureSingularLong.rep₂) :=
  N6PureSingularLong.GL2_scaled_lift_action
    (k := k) (u := u) (a := α) (b := β) (c := γ) (d := δ) C

/-- The explicit pointwise subgroup on the radical extension of the pure singular
length-two orbit combines with the determinant-scaled `GL₂` lift exactly as
expected. -/
theorem pureSingularLong_scaleLowerHankel_GL2_product_lift_action
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix N6PureSingularLong.W N6PureSingularLong.I k)
    (x y z w α β γ δ : k) :
    let gU : Matrix N6PureSingularLong.V N6PureSingularLong.V k :=
      (Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) 1)
        0
        0
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) x y z w)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) *
      (Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))))
    let P : Matrix N5PureSingularLong.W N5PureSingularLong.W k :=
      N5PureSingularLong.sym2Raw (k := k) α β γ δ
    let Q : Matrix N5PureSingularLong.I N5PureSingularLong.I k :=
      N5PureSingularLong.adjointT (k := k) α β γ δ
    let h : Matrix N6PureSingularLong.W N6PureSingularLong.W k := Matrix.fromBlocks P 0 0 Q
    let gL : Matrix N6PureSingularLong.V N6PureSingularLong.V k :=
      Matrix.fromBlocks (N6PureSingularLong.scalarBlock (k := k) 1) 0 0 h
    N6PureSingularLong.ActBivector N6PureSingularLong.rep₁ (gU * gL) =
        N5PureSingularLong.Delta α β γ δ •
          (α • N6PureSingularLong.rep₁ + γ • N6PureSingularLong.rep₂) ∧
      N6PureSingularLong.ActBivector N6PureSingularLong.rep₂ (gU * gL) =
        N5PureSingularLong.Delta α β γ δ •
          (β • N6PureSingularLong.rep₁ + δ • N6PureSingularLong.rep₂) :=
  N6PureSingularLong.scaleLowerHankel_GL2_product_lift_action
    (k := k) (u := u) (t := t) ht (C := C) (x := x) (y := y) (z := z) (w := w)
    (a := α) (b := β) (c := γ) (d := δ)

/-- Right-multiplying the determinant-scaled `GL₂` lift by the explicit pointwise
subgroup on the radical extension of the pure singular length-two orbit does not change
the quotient action. -/
theorem pureSingularLong_GL2_scaleLowerHankel_right_product_lift_action
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix N6PureSingularLong.W N6PureSingularLong.I k)
    (x y z w α β γ δ : k) :
    let P : Matrix N5PureSingularLong.W N5PureSingularLong.W k :=
      N5PureSingularLong.sym2Raw (k := k) α β γ δ
    let Q : Matrix N5PureSingularLong.I N5PureSingularLong.I k :=
      N5PureSingularLong.adjointT (k := k) α β γ δ
    let h : Matrix N6PureSingularLong.W N6PureSingularLong.W k := Matrix.fromBlocks P 0 0 Q
    let gL : Matrix N6PureSingularLong.V N6PureSingularLong.V k :=
      Matrix.fromBlocks (N6PureSingularLong.scalarBlock (k := k) 1) 0 0 h
    let gU : Matrix N6PureSingularLong.V N6PureSingularLong.V k :=
      (Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) 1)
        0
        0
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) x y z w)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) *
      (Matrix.fromBlocks
        (N6PureSingularLong.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))))
    N6PureSingularLong.ActBivector N6PureSingularLong.rep₁ (gL * gU) =
        N5PureSingularLong.Delta α β γ δ •
          (α • N6PureSingularLong.rep₁ + γ • N6PureSingularLong.rep₂) ∧
      N6PureSingularLong.ActBivector N6PureSingularLong.rep₂ (gL * gU) =
        N5PureSingularLong.Delta α β γ δ •
          (β • N6PureSingularLong.rep₁ + δ • N6PureSingularLong.rep₂) :=
  N6PureSingularLong.GL2_scaleLowerHankel_right_product_lift_action
    (k := k) (u := u) (t := t) ht (C := C) (x := x) (y := y) (z := z) (w := w)
    (a := α) (b := β) (c := γ) (d := δ)

/-- Independent pointwise scalings on the two pure singular `3`-dimensional summands
fix the direct-sum pair. -/
theorem doublePureSingular_pointwise_scale_family
    (a b : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    N6DoublePureSingular.FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        (N3PureSingular.pointwiseScale (k := k) b)) :=
  N6DoublePureSingular.pointwise_scale_family (k := k) a b ha hb

/-- For block diagonal matrices, fixing the direct-sum pure singular pair is exactly
the same as fixing the pure singular pair on each summand. -/
theorem doublePureSingular_blockDiagonal_iff
    (G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k)
    (H : Matrix N6DoublePureSingular.I N6DoublePureSingular.I k) :
    N6DoublePureSingular.FixesPairBivector (Matrix.fromBlocks G 0 0 H) ↔
      N3PureSingular.FixesPairBivector G ∧ N3PureSingular.FixesPairBivector H :=
  N6DoublePureSingular.fixesPair_blockDiagonal_iff (k := k) G H

/-- Shape form of the block diagonal pointwise criterion for the direct-sum pure
singular row. -/
theorem doublePureSingular_blockDiagonal_iff_shape
    (G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k)
    (H : Matrix N6DoublePureSingular.I N6DoublePureSingular.I k) :
    N6DoublePureSingular.FixesPairBivector (Matrix.fromBlocks G 0 0 H) ↔
      G 0 0 ≠ 0 ∧ H 0 0 ≠ 0 ∧
        G = N3PureSingular.pureSingularShape (k := k) (G 0 0) (G 2 0) (G 2 1) ∧
        H = N3PureSingular.pureSingularShape (k := k) (H 0 0) (H 2 0) (H 2 1) :=
  N6DoublePureSingular.fixesPair_blockDiagonal_iff_shape (k := k) G H

/-- The full pure singular pointwise shape on each summand gives an explicit block
diagonal pointwise family on the direct-sum `S_2^2` row. -/
theorem doublePureSingular_pointwise_shape_family
    (a b x y r s : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    N6DoublePureSingular.FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b r s)) :=
  N6DoublePureSingular.pointwise_shape_family (k := k) a b x y r s ha hb

/-- The direct-sum pure singular `S_2^2` row admits a concrete six-parameter unipotent
pointwise family, matching the `\operatorname{Sym}_2 \oplus \operatorname{Sym}_2`
kernel visible in the Magma local model. -/
theorem doublePureSingular_coupled_unipotent_family
    (x y z w r s : k) :
    N6DoublePureSingular.FixesPairBivector
      (N6DoublePureSingular.coupledUnipotent (k := k) x y z w r s) :=
  N6DoublePureSingular.coupled_unipotent_family (k := k) x y z w r s

/-- The coupled unipotent family combines with the independent pointwise scalings on
the two pure singular summands to give a larger explicit pointwise subgroup on the
direct-sum `S_2^2` row. -/
theorem doublePureSingular_coupled_scale_product_pointwise
    (x y z w r s a b : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gU : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      N6DoublePureSingular.coupledUnipotent (k := k) x y z w r s
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        (N3PureSingular.pointwiseScale (k := k) b)
    N6DoublePureSingular.FixesPairBivector (gU * gS) :=
  N6DoublePureSingular.coupled_scale_product_pointwise
    (k := k) x y z w r s a b ha hb

/-- The coupled unipotent family also combines with the full pointwise pure singular
shape on the two summands. -/
theorem doublePureSingular_coupled_shape_product_pointwise
    (u v z w r s a b x y p q : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gU : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      N6DoublePureSingular.coupledUnipotent (k := k) u v z w r s
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    N6DoublePureSingular.FixesPairBivector (gU * gS) :=
  N6DoublePureSingular.coupled_shape_product_pointwise
    (k := k) u v z w r s a b x y p q ha hb

/-- Right-multiplying a block-diagonal pure singular pointwise shape by the coupled
unipotent family does not change the exact pointwise criterion on the direct-sum
`S_2^2` row. -/
theorem doublePureSingular_blockDiagonal_coupled_right_product_iff_shape
    (x y z w r s : k)
    (G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k)
    (H : Matrix N6DoublePureSingular.I N6DoublePureSingular.I k) :
    N6DoublePureSingular.FixesPairBivector
        ((Matrix.fromBlocks G 0 0 H) *
          N6DoublePureSingular.coupledUnipotent (k := k) x y z w r s) ↔
      G 0 0 ≠ 0 ∧ H 0 0 ≠ 0 ∧
        G = N3PureSingular.pureSingularShape (k := k) (G 0 0) (G 2 0) (G 2 1) ∧
        H = N3PureSingular.pureSingularShape (k := k) (H 0 0) (H 2 0) (H 2 1) :=
  N6DoublePureSingular.blockDiagonal_coupled_right_product_iff_shape
    (k := k) x y z w r s G H

/-- The transported Magma `diag(T, T, T^{-T})` pointwise `GL₂` family fixes the
direct-sum pure singular `S_2^2` pair pointwise. -/
theorem doublePureSingular_magmaPointwiseGL2_family
    (a b c d : k)
    (hΔ : a * d - b * c ≠ 0) :
    N6DoublePureSingular.FixesPairBivector
      (N6DoublePureSingular.magmaPointwiseGL2 (k := k) a b c d) :=
  N6DoublePureSingular.magmaPointwiseGL2_family (k := k) a b c d hΔ

/-- The full transported Magma local kernel cell on the direct-sum pure singular
`S_2^2` row fixes the pair pointwise. -/
theorem doublePureSingular_magmaKernelCell_pointwise
    (a b c d p q r s u v : k)
    (hΔ : a * d - b * c ≠ 0) :
    N6DoublePureSingular.FixesPairBivector
      (N6DoublePureSingular.magmaKernelCell (k := k) a b c d p q r s u v) :=
  N6DoublePureSingular.magmaKernelCell_pointwise
    (k := k) a b c d p q r s u v hΔ

/-- Applying the same `GL₂` lift on both pure singular summands produces the expected
quotient action on the direct-sum pair. -/
theorem doublePureSingular_GL2_lift_action
    (α β γ δ : k) :
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₁ g =
        α • N6DoublePureSingular.rep₁ + β • N6DoublePureSingular.rep₂ ∧
      N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₂ g =
        γ • N6DoublePureSingular.rep₁ + δ • N6DoublePureSingular.rep₂ :=
  N6DoublePureSingular.GL2_lift_action (k := k) (α := α) (β := β) (γ := γ) (δ := δ)

/-- The simultaneous `GL₂` lift on the direct-sum pure singular row has determinant
`(αδ-βγ)^2`. -/
theorem doublePureSingular_GL2_lift_det
    (α β γ δ : k) :
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    Matrix.det g = (α * δ - β * γ) ^ 2 :=
  N6DoublePureSingular.GL2_lift_det (k := k) α β γ δ

/-- The six-parameter coupled unipotent family and the simultaneous `GL₂` lift combine
as expected for the semidirect-product structure on the direct-sum pure singular
`S_2^2` row. -/
theorem doublePureSingular_coupled_GL2_product_lift_action
    (x y z w r s α β γ δ : k) :
    let gU : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      N6DoublePureSingular.coupledUnipotent (k := k) x y z w r s
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₁ (gU * gL) =
        α • N6DoublePureSingular.rep₁ + β • N6DoublePureSingular.rep₂ ∧
      N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₂ (gU * gL) =
        γ • N6DoublePureSingular.rep₁ + δ • N6DoublePureSingular.rep₂ :=
  N6DoublePureSingular.coupled_GL2_product_lift_action
    (k := k) x y z w r s α β γ δ

/-- The full block-diagonal pure singular pointwise shape also combines with the
simultaneous `GL₂` lift with the expected quotient action. -/
theorem doublePureSingular_shape_GL2_product_lift_action
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₁ (gS * gL) =
        α • N6DoublePureSingular.rep₁ + β • N6DoublePureSingular.rep₂ ∧
      N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₂ (gS * gL) =
        γ • N6DoublePureSingular.rep₁ + δ • N6DoublePureSingular.rep₂ :=
  N6DoublePureSingular.shape_GL2_product_lift_action
    (k := k) a b x y p q α β γ δ ha hb

/-- Left-multiplying the simultaneous `GL₂` lift by the full block-diagonal pure
singular pointwise shape changes the determinant only by the pointwise factor `ab`. -/
theorem doublePureSingular_shape_GL2_product_lift_det
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    Matrix.det (gS * gL) = (a * b) * (α * δ - β * γ) ^ 2 :=
  N6DoublePureSingular.shape_GL2_product_lift_det
    (k := k) a b x y p q α β γ δ ha hb

/-- Right-multiplying the simultaneous `GL₂` lift by the coupled unipotent family does
not change the quotient action on the direct-sum pure singular pair. -/
theorem doublePureSingular_GL2_coupled_right_product_lift_action
    (x y z w r s α β γ δ : k) :
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    let gU : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      N6DoublePureSingular.coupledUnipotent (k := k) x y z w r s
    N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₁ (gL * gU) =
        α • N6DoublePureSingular.rep₁ + β • N6DoublePureSingular.rep₂ ∧
      N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₂ (gL * gU) =
        γ • N6DoublePureSingular.rep₁ + δ • N6DoublePureSingular.rep₂ :=
  N6DoublePureSingular.GL2_coupled_right_product_lift_action
    (k := k) x y z w r s α β γ δ

/-- Right-multiplying the simultaneous `GL₂` lift by the full block-diagonal pure
singular pointwise shape does not change the quotient action. -/
theorem doublePureSingular_GL2_shape_right_product_lift_action
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₁ (gL * gS) =
        α • N6DoublePureSingular.rep₁ + β • N6DoublePureSingular.rep₂ ∧
      N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₂ (gL * gS) =
        γ • N6DoublePureSingular.rep₁ + δ • N6DoublePureSingular.rep₂ :=
  N6DoublePureSingular.GL2_shape_right_product_lift_action
    (k := k) a b x y p q α β γ δ ha hb

/-- Right-multiplying the simultaneous `GL₂` lift by the full block-diagonal pure
singular pointwise shape changes the determinant only by the pointwise factor `ab`. -/
theorem doublePureSingular_GL2_shape_right_product_lift_det
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    Matrix.det (gL * gS) = (a * b) * (α * δ - β * γ) ^ 2 :=
  N6DoublePureSingular.GL2_shape_right_product_lift_det
    (k := k) a b x y p q α β γ δ ha hb

/-- The full coupled-and-shape pointwise subgroup combines with the simultaneous
`GL₂` lift with the expected quotient action. -/
theorem doublePureSingular_coupledShape_GL2_product_lift_action
    (u v z w r s a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gU : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      N6DoublePureSingular.coupledUnipotent (k := k) u v z w r s
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let gP : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k := gU * gS
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₁ (gP * gL) =
        α • N6DoublePureSingular.rep₁ + β • N6DoublePureSingular.rep₂ ∧
      N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₂ (gP * gL) =
        γ • N6DoublePureSingular.rep₁ + δ • N6DoublePureSingular.rep₂ :=
  N6DoublePureSingular.coupledShape_GL2_product_lift_action
    (k := k) u v z w r s a b x y p q α β γ δ ha hb

/-- Right-multiplying the simultaneous `GL₂` lift by the combined coupled-and-shape
pointwise subgroup does not change the quotient action. -/
theorem doublePureSingular_GL2_coupledShape_right_product_lift_action
    (u v z w r s a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let G : Matrix N6DoublePureSingular.W N6DoublePureSingular.W k :=
      N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    let gU : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      N6DoublePureSingular.coupledUnipotent (k := k) u v z w r s
    let gS : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let gP : Matrix N6DoublePureSingular.V N6DoublePureSingular.V k := gU * gS
    N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₁ (gL * gP) =
        α • N6DoublePureSingular.rep₁ + β • N6DoublePureSingular.rep₂ ∧
      N6DoublePureSingular.ActBivector N6DoublePureSingular.rep₂ (gL * gP) =
        γ • N6DoublePureSingular.rep₁ + δ • N6DoublePureSingular.rep₂ :=
  N6DoublePureSingular.GL2_coupledShape_right_product_lift_action
    (k := k) u v z w r s a b x y p q α β γ δ ha hb

/-- Summary form of the weighted two-point `n = 6` pointwise stabilizer calculation. -/
theorem weightedTwoPoint_pointwise_bivector_iff
    (A : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k)
    (B : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.I k)
    (C : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.W k)
    (D : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k) :
    N6WeightedTwoPoint.FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ C = 0 ∧ N4.FixesBivector (N4.splitRep₁ (k := k)) A ∧ D.det = 1 :=
  N6WeightedTwoPoint.fixesPair_fromBlocks_iff (k := k) (A := A) (B := B) (C := C) (D := D)

/-- The standard torus lift for the weighted two-point `n = 6` orbit. -/
theorem weightedTwoPoint_torus_lift_action
    (u v : k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k := Matrix.fromBlocks A 0 0 A
    let E : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
    let g : Matrix N6WeightedTwoPoint.V N6WeightedTwoPoint.V k := Matrix.fromBlocks H 0 0 E
    N6WeightedTwoPoint.ActBivector N6WeightedTwoPoint.rep₁ g =
        u • N6WeightedTwoPoint.rep₁ + (v - u) • N6WeightedTwoPoint.rep₂ ∧
    N6WeightedTwoPoint.ActBivector N6WeightedTwoPoint.rep₂ g =
        v • N6WeightedTwoPoint.rep₂ :=
  N6WeightedTwoPoint.torus_lift_action (k := k) (u := u) (v := v)

/-- The determinant of the standard torus lift on the weighted two-point `n = 6`
orbit is `u^2 v`. -/
theorem weightedTwoPoint_torus_lift_det
    (u v : k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k := Matrix.fromBlocks A 0 0 A
    let E : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
    let g : Matrix N6WeightedTwoPoint.V N6WeightedTwoPoint.V k := Matrix.fromBlocks H 0 0 E
    Matrix.det g = (u * u) * v :=
  N6WeightedTwoPoint.torus_lift_det (k := k) (u := u) (v := v)

/-- The obvious block-diagonal pointwise family on the weighted two-point orbit. -/
theorem weightedTwoPoint_pointwise_family
    (A : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k)
    (D : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k)
    (hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1) :
    N6WeightedTwoPoint.FixesPairBivector (Matrix.fromBlocks A 0 0 D) :=
  N6WeightedTwoPoint.pointwise_family (k := k) (A := A) (D := D) hA hD

/-- The block-diagonal pointwise subgroup on the weighted two-point orbit combines
with the standard torus lift with the expected quotient action. -/
theorem weightedTwoPoint_pointwise_torus_product_lift_action
    (A : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k)
    (D : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k)
    (hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1)
    (u v : k) :
    let A0 : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k := Matrix.fromBlocks A0 0 0 A0
    let E : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
    let gU : Matrix N6WeightedTwoPoint.V N6WeightedTwoPoint.V k := Matrix.fromBlocks A 0 0 D
    let gL : Matrix N6WeightedTwoPoint.V N6WeightedTwoPoint.V k := Matrix.fromBlocks H 0 0 E
    N6WeightedTwoPoint.ActBivector N6WeightedTwoPoint.rep₁ (gU * gL) =
        u • N6WeightedTwoPoint.rep₁ + (v - u) • N6WeightedTwoPoint.rep₂ ∧
      N6WeightedTwoPoint.ActBivector N6WeightedTwoPoint.rep₂ (gU * gL) =
        v • N6WeightedTwoPoint.rep₂ :=
  N6WeightedTwoPoint.pointwise_torus_product_lift_action
    (k := k) (A := A) (D := D) hA hD (u := u) (v := v)

/-- Right-multiplying the standard torus lift by the block-diagonal pointwise subgroup
on the weighted two-point orbit does not change the quotient action. -/
theorem weightedTwoPoint_torus_pointwise_right_product_lift_action
    (A : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k)
    (D : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k)
    (hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1)
    (u v : k) :
    let A0 : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k := Matrix.fromBlocks A0 0 0 A0
    let E : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
    let gL : Matrix N6WeightedTwoPoint.V N6WeightedTwoPoint.V k := Matrix.fromBlocks H 0 0 E
    let gU : Matrix N6WeightedTwoPoint.V N6WeightedTwoPoint.V k := Matrix.fromBlocks A 0 0 D
    N6WeightedTwoPoint.ActBivector N6WeightedTwoPoint.rep₁ (gL * gU) =
        u • N6WeightedTwoPoint.rep₁ + (v - u) • N6WeightedTwoPoint.rep₂ ∧
      N6WeightedTwoPoint.ActBivector N6WeightedTwoPoint.rep₂ (gL * gU) =
        v • N6WeightedTwoPoint.rep₂ :=
  N6WeightedTwoPoint.torus_pointwise_right_product_lift_action
    (k := k) (A := A) (D := D) hA hD (u := u) (v := v)

/-- On the weighted two-point line in dimension `6`, rank drop occurs exactly on
the two distinguished projective points. -/
theorem weightedTwoPoint_det_zero_iff
    (a b : k) :
    Matrix.det
        (a • (N6WeightedTwoPoint.rep₁ (k := k)) + b • (N6WeightedTwoPoint.rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  N6WeightedTwoPoint.det_zero_iff (k := k) (a := a) (b := b)

/-- On the mixed one-point line in dimension `6`, rank drop occurs exactly at the
repeated-support point. -/
theorem mixedOnePoint_det_zero_iff
    (a b : k) :
    Matrix.det
        (a • (N6MixedOnePoint.rep₁ (k := k)) + b • (N6MixedOnePoint.rep₂ (k := k))) = 0 ↔
      a = 0 :=
  N6MixedOnePoint.det_zero_iff (k := k) (a := a) (b := b)

/-- The obvious `SL₂ × SL₂` pointwise family on the mixed one-point orbit. -/
theorem mixedOnePoint_pointwise_levi_family
    (A B H : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    N6MixedOnePoint.FixesPairBivector
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H) :=
  N6MixedOnePoint.pointwise_levi_family
    (k := k) (A := A) (B := B) (H := H) hA hB hH

/-- A general coupled pointwise family on the mixed one-point orbit. -/
theorem mixedOnePoint_coupled_pointwise_family
    (Q : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k)
    (R : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k)
    (h11 : Q * N4.J * Qᵀ = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k))
    (h12 : N4.onePointRep₁ * Rᵀ + Q * N4.J = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k))
    (h21 : R * N4.onePointRep₁ + N4.J * Qᵀ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k))
    (h22 : R * N4.onePointRep₁ * Rᵀ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))
    (h₂left : R * N4.onePointRep₂ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k))
    (h₂right : N4.onePointRep₂ * Rᵀ = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k)) :
    N6MixedOnePoint.FixesPairBivector (Matrix.fromBlocks 1 Q R 1) :=
  N6MixedOnePoint.coupled_pointwise_family
    (k := k) (Q := Q) (R := R) h11 h12 h21 h22 h₂left h₂right

/-- A simple explicit one-parameter coupled family on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledE12_pointwise
    (t : k) :
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ12 (k := k) t)
        (N6MixedOnePoint.coupledR12 (k := k) t)
        1) :=
  N6MixedOnePoint.coupledE12_pointwise (k := k) t

/-- A second simple explicit one-parameter coupled family on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledE21_pointwise
    (t : k) :
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ21 (k := k) t)
        (N6MixedOnePoint.coupledR21 (k := k) t)
        1) :=
  N6MixedOnePoint.coupledE21_pointwise (k := k) t

/-- A third simple explicit one-parameter coupled family on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledE11_pointwise
    (t : k) :
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) t)
        (N6MixedOnePoint.coupledR11 (k := k) t)
        1) :=
  N6MixedOnePoint.coupledE11_pointwise (k := k) t

/-- A fourth simple explicit one-parameter coupled family on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledE22_pointwise
    (t : k) :
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t)
        (N6MixedOnePoint.coupledR22 (k := k) t)
        1) :=
  N6MixedOnePoint.coupledE22_pointwise (k := k) t

/-- Each basic explicit mixed `E12` generator has determinant `1`. -/
theorem mixedOnePoint_coupledE12_det
    (t : k) :
    Matrix.det
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ12 (k := k) t)
        (N6MixedOnePoint.coupledR12 (k := k) t)
        1) = 1 :=
  N6MixedOnePoint.coupledE12_det (k := k) t

/-- Each basic explicit mixed `E21` generator has determinant `1`. -/
theorem mixedOnePoint_coupledE21_det
    (t : k) :
    Matrix.det
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ21 (k := k) t)
        (N6MixedOnePoint.coupledR21 (k := k) t)
        1) = 1 :=
  N6MixedOnePoint.coupledE21_det (k := k) t

/-- Each basic explicit mixed `E11` generator has determinant `1`. -/
theorem mixedOnePoint_coupledE11_det
    (t : k) :
    Matrix.det
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) t)
        (N6MixedOnePoint.coupledR11 (k := k) t)
        1) = 1 :=
  N6MixedOnePoint.coupledE11_det (k := k) t

/-- Each basic explicit mixed `E22` generator has determinant `1`. -/
theorem mixedOnePoint_coupledE22_det
    (t : k) :
    Matrix.det
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t)
        (N6MixedOnePoint.coupledR22 (k := k) t)
        1) = 1 :=
  N6MixedOnePoint.coupledE22_det (k := k) t

/-- The older hard-coded `crossUnipotent₂` direction is not pointwise on the mixed
one-point orbit in the bivector model unless the parameter vanishes. -/
theorem mixedOnePoint_cross2_not_pointwise
    {t : k}
    (ht : t ≠ 0) :
    ¬ N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (t • N6MixedOnePoint.crossQ₂ (k := k))
        (t • N6MixedOnePoint.crossR₂ (k := k))
        1) :=
  N6MixedOnePoint.cross₂_not_pointwise (k := k) ht

/-- The older hard-coded `crossUnipotent₁` direction is also not pointwise on the mixed
one-point orbit in the bivector model unless the parameter vanishes. -/
theorem mixedOnePoint_cross1_not_pointwise
    {t : k}
    (ht : t ≠ 0) :
    ¬ N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (t • N6MixedOnePoint.crossQ₁ (k := k))
        (t • N6MixedOnePoint.crossR₁ (k := k))
        1) :=
  N6MixedOnePoint.cross₁_not_pointwise (k := k) ht

/-- The subgroup generated by the two explicit mixed coupled families still fixes the
pair pointwise. -/
theorem mixedOnePoint_coupled_product_pointwise
    (s t : k) :
    N6MixedOnePoint.FixesPairBivector
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          1) *
       (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          1)) :=
  N6MixedOnePoint.coupledE12E21_product_pointwise (k := k) s t

/-- The two diagonal explicit mixed coupled families also generate a concrete pointwise
subgroup on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledDiagonal_product_pointwise
    (s t : k) :
    N6MixedOnePoint.FixesPairBivector
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ11 (k := k) s)
          (N6MixedOnePoint.coupledR11 (k := k) s)
          1) *
       (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ22 (k := k) t)
          (N6MixedOnePoint.coupledR22 (k := k) t)
          1)) :=
  N6MixedOnePoint.coupledE11E22_product_pointwise (k := k) s t

/-- All four explicit mixed coupled one-parameter families combine to give a concrete
four-parameter pointwise subgroup. -/
theorem mixedOnePoint_coupledExplicit_product_pointwise
    (s₁ t₁ s₂ t₂ : k) :
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ12 (k := k) s₁)
        (N6MixedOnePoint.coupledR12 (k := k) s₁)
        1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ21 (k := k) t₁)
        (N6MixedOnePoint.coupledR21 (k := k) t₁)
        1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) s₂)
        (N6MixedOnePoint.coupledR11 (k := k) s₂)
        1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t₂)
        (N6MixedOnePoint.coupledR22 (k := k) t₂)
        1
    N6MixedOnePoint.FixesPairBivector ((g12 * g21) * (g11 * g22)) :=
  N6MixedOnePoint.coupledExplicit_product_pointwise (k := k) s₁ t₁ s₂ t₂

/-- The four basic mixed root families combine into a closed `4`-parameter normal form,
matching the `G_a^4` quotient family suggested by the Magma local model. -/
theorem mixedOnePoint_coupledKernelFrom_eq_product
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    N6MixedOnePoint.coupledKernelFrom (k := k) X =
      ((Matrix.fromBlocks
            1
            (N6MixedOnePoint.coupledQ12 (k := k) (X 0 1))
            (N6MixedOnePoint.coupledR12 (k := k) (X 0 1))
            1) *
        (Matrix.fromBlocks
            1
            (N6MixedOnePoint.coupledQ21 (k := k) (X 1 0))
            (N6MixedOnePoint.coupledR21 (k := k) (X 1 0))
            1)) *
        ((Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ11 (k := k) (X 0 0))
              (N6MixedOnePoint.coupledR11 (k := k) (X 0 0))
              1) *
          (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ22 (k := k) (X 1 1))
              (N6MixedOnePoint.coupledR22 (k := k) (X 1 1))
              1)) :=
  N6MixedOnePoint.coupledKernelFrom_eq_product (k := k) X

/-- The closed `4`-parameter mixed quotient family fixes the mixed one-point pair
pointwise. -/
theorem mixedOnePoint_coupledKernelFrom_pointwise
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    N6MixedOnePoint.FixesPairBivector (N6MixedOnePoint.coupledKernelFrom (k := k) X) :=
  N6MixedOnePoint.coupledKernelFrom_pointwise (k := k) X

/-- The closed `4`-parameter mixed quotient family has determinant `1`. -/
theorem mixedOnePoint_coupledKernelFrom_det
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    Matrix.det (N6MixedOnePoint.coupledKernelFrom (k := k) X) = 1 :=
  N6MixedOnePoint.coupledKernelFrom_det (k := k) X

/-- The central upper-middle layer on the mixed Magma unipotent cell fixes the pair
exactly when its trace vanishes. -/
theorem mixedOnePoint_centralBlock_iff_trace_zero
    (M : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    N6MixedOnePoint.FixesPairBivector (N6MixedOnePoint.centralBlock (k := k) M) ↔
      M 0 0 + M 1 1 = 0 :=
  N6MixedOnePoint.centralBlock_pointwise_iff (k := k) M

/-- On the full Magma `3 × 3` unipotent cell, pointwise fixing is equivalent to the
transported mixed-root upper-right block together with the single trace condition on the
central upper-middle block. -/
theorem mixedOnePoint_Block3_unipotent_iff
    (M U X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    N6MixedOnePoint.FixesPairBivector
      (N6MixedOnePoint.Block3
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) M U
        0 (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) 0
        0 X (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) ↔
      U = N6MixedOnePoint.coupledUpperFrom (k := k) X ∧
        M 0 0 + M 1 1 = -X.det :=
  N6MixedOnePoint.fixesPair_Block3_unipotent_iff (k := k) (M := M) (U := U) (X := X)

/-- Equivalently, the full Magma unipotent cell fixes the mixed pair exactly when it
splits as a trace-zero central factor times the concrete `G_a^4` quotient family. -/
theorem mixedOnePoint_Block3_unipotent_iff_exists_centralBlock_mul_coupledKernelFrom
    (M U X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    N6MixedOnePoint.FixesPairBivector
      (N6MixedOnePoint.Block3
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) M U
        0 (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) 0
        0 X (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) ↔
      ∃ Z : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k,
        Z 0 0 + Z 1 1 = 0 ∧
        N6MixedOnePoint.Block3
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) M U
          0 (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) 0
          0 X (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) =
          N6MixedOnePoint.centralBlock (k := k) Z *
            N6MixedOnePoint.coupledKernelFrom (k := k) X :=
  N6MixedOnePoint.fixesPair_Block3_unipotent_iff_exists_centralBlock_mul_coupledKernelFrom
    (k := k) (M := M) (U := U) (X := X)

/-- Combining the explicit coupled unipotent families with the Levi family gives a
concrete pointwise subgroup on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledLevi_product_pointwise
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (Matrix.fromBlocks
        (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
        (N6MixedOnePoint.coupledQ12 (k := k) s)
        (N6MixedOnePoint.coupledR12 (k := k) s)
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) *
      (Matrix.fromBlocks
        (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
        (N6MixedOnePoint.coupledQ21 (k := k) t)
        (N6MixedOnePoint.coupledR21 (k := k) t)
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H
    N6MixedOnePoint.FixesPairBivector (gU * gL) :=
  N6MixedOnePoint.coupledLevi_product_pointwise
    (k := k)
    (s := s)
    (t := t)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH

/-- Combining the diagonal explicit mixed coupled subgroup with the Levi family gives
a concrete pointwise subgroup on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledDiagonalLevi_product_pointwise
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) s)
        (N6MixedOnePoint.coupledR11 (k := k) s)
        1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t)
        (N6MixedOnePoint.coupledR22 (k := k) t)
        1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    N6MixedOnePoint.FixesPairBivector gU :=
  N6MixedOnePoint.coupledDiagonalLevi_product_pointwise
    (k := k)
    (s := s)
    (t := t)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH

/-- Combining the full explicit four-parameter mixed coupled subgroup with the Levi
family gives a concrete pointwise subgroup on the mixed one-point orbit. -/
theorem mixedOnePoint_coupledExplicitLevi_product_pointwise
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ12 (k := k) s₁)
        (N6MixedOnePoint.coupledR12 (k := k) s₁)
        1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ21 (k := k) t₁)
        (N6MixedOnePoint.coupledR21 (k := k) t₁)
        1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) s₂)
        (N6MixedOnePoint.coupledR11 (k := k) s₂)
        1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t₂)
        (N6MixedOnePoint.coupledR22 (k := k) t₂)
        1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    N6MixedOnePoint.FixesPairBivector gU :=
  N6MixedOnePoint.coupledExplicitLevi_product_pointwise
    (k := k)
    (s₁ := s₁)
    (t₁ := t₁)
    (s₂ := s₂)
    (t₂ := t₂)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH

/-- The standard Borel lift for the mixed one-point `n = 6` orbit. -/
theorem mixedOnePoint_borel_lift_action
    (a b : k)
    (ha : a ≠ 0) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let g : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H 0 0 E
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ g =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ g =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.borel_lift_action (k := k) (a := a) (b := b) ha

/-- The determinant of the standard Borel lift on the mixed one-point
`n = 6` orbit is `a^3`. -/
theorem mixedOnePoint_borel_lift_det
    (a b : k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let g : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H 0 0 E
    Matrix.det g = a * a * a :=
  N6MixedOnePoint.borel_lift_det (k := k) (a := a) (b := b)

/-- The explicit block diagonal basis change carrying the Lean mixed representative to the
native Magma local model. -/
theorem mixedOnePoint_magmaChange_act_rep₁ :
    N6MixedOnePoint.ActBivector
      (N6MixedOnePoint.rep₁ (k := k))
      (N6MixedOnePoint.magmaChange (k := k)) =
      N6MixedOnePoint.magmaRep₁ (k := k) :=
  N6MixedOnePoint.magmaChange_act_rep₁ (k := k)

/-- The same explicit basis change also carries the second mixed basis vector to the
native Magma local model. -/
theorem mixedOnePoint_magmaChange_act_rep₂ :
    N6MixedOnePoint.ActBivector
      (N6MixedOnePoint.rep₂ (k := k))
      (N6MixedOnePoint.magmaChange (k := k)) =
      N6MixedOnePoint.magmaRep₂ (k := k) :=
  N6MixedOnePoint.magmaChange_act_rep₂ (k := k)

/-- The explicit mixed pointwise subgroup combines with the standard Borel lift with the
expected action on the chosen basis. -/
theorem mixedOnePoint_coupled_pointwise_iff
    (Q : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k)
    (R : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k) :
    N6MixedOnePoint.FixesPairBivector (Matrix.fromBlocks 1 Q R 1) ↔
      Q * N4.J * Qᵀ = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k) ∧
        N4.onePointRep₁ * Rᵀ + Q * N4.J = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k) ∧
        R * N4.onePointRep₁ + N4.J * Qᵀ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k) ∧
        R * N4.onePointRep₁ * Rᵀ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) ∧
        R * N4.onePointRep₂ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k) ∧
        N4.onePointRep₂ * Rᵀ = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k) :=
  N6MixedOnePoint.fixesPair_identityOffDiagonal_iff (k := k) (Q := Q) (R := R)

/-- Any identity-plus-off-diagonal pointwise fixer on the mixed row already lies on the
intrinsic surviving coupled cell, with rank-drop on the attached `2 × 2` block. -/
theorem mixedOnePoint_identityOffDiagonal_iff_exists_coupledFrom_det_zero
    (Q : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k)
    (R : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k) :
    N6MixedOnePoint.FixesPairBivector (Matrix.fromBlocks 1 Q R 1) ↔
      ∃ X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k,
        Q = N6MixedOnePoint.coupledQFrom (k := k) X ∧
        R = N6MixedOnePoint.coupledRFrom (k := k) X ∧
        X.det = 0 :=
  N6MixedOnePoint.fixesPair_identityOffDiagonal_iff_exists_coupledFrom_det_zero
    (k := k) (Q := Q) (R := R)

/-- In Lean coordinates, one of the native Magma mixed root generators is exactly the
`E22` coupled family in `3 × 3` block form. -/
theorem mixedOnePoint_coupledE22_Block3
    (t : k) :
    Matrix.fromBlocks
      1
      (N6MixedOnePoint.coupledQ22 (k := k) t)
      (N6MixedOnePoint.coupledR22 (k := k) t)
      1 =
      N6MixedOnePoint.Block3
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
        0
        (-t • N6MixedOnePoint.E11 (k := k))
        0
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
        0
        0
        (t • N6MixedOnePoint.E22 (k := k))
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :=
  N6MixedOnePoint.coupledE22_Block3 (k := k) t

/-- In Lean coordinates, the second transported Magma mixed root generator is exactly the
`E12` coupled family in `3 × 3` block form. -/
theorem mixedOnePoint_magmaMixed2_Block3
    (t : k) :
    Matrix.fromBlocks
      1
      (N6MixedOnePoint.magmaMixedQ2 (k := k) t)
      (N6MixedOnePoint.magmaMixedR2 (k := k) t)
      1 =
      N6MixedOnePoint.Block3
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
        0
        (-t • N6MixedOnePoint.E12 (k := k))
        0
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
        0
        0
        (-t • N6MixedOnePoint.E12 (k := k))
        (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :=
  N6MixedOnePoint.magmaMixed2_Block3 (k := k) t

/-- In the current Lean basis, the second transported Magma mixed root subgroup is
literally the existing `E12` coupled family with parameter `-t`. -/
theorem mixedOnePoint_magmaMixed2_eq_coupledE12
    (t : k) :
    Matrix.fromBlocks
      1
      (N6MixedOnePoint.magmaMixedQ2 (k := k) t)
      (N6MixedOnePoint.magmaMixedR2 (k := k) t)
      1 =
    Matrix.fromBlocks
      1
      (N6MixedOnePoint.coupledQ12 (k := k) (-t))
      (N6MixedOnePoint.coupledR12 (k := k) (-t))
      1 :=
  N6MixedOnePoint.magmaMixed2_eq_coupledE12 (k := k) t

/-- The second transported Magma mixed root subgroup fixes the mixed one-point pair
pointwise; in the current Lean coordinates it is the existing `E12` coupled family. -/
theorem mixedOnePoint_magmaMixed2_pointwise
    (t : k) :
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.magmaMixedQ2 (k := k) t)
        (N6MixedOnePoint.magmaMixedR2 (k := k) t)
        1) :=
  N6MixedOnePoint.magmaMixed2_pointwise (k := k) t

/-- On the surviving mixed coupled cell, pointwise fixing is exactly the
rank-drop condition on the attached `2 × 2` block. -/
theorem mixedOnePoint_coupledFrom_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQFrom (k := k) X)
        (N6MixedOnePoint.coupledRFrom (k := k) X)
        1) ↔
      X.det = 0 :=
  N6MixedOnePoint.coupledFrom_pointwise_iff_det_zero (k := k) X

/-- On the surviving mixed coupled cell, pointwise fixing is equivalent to the
explicit rank-one factorization `X = u vᵀ`. -/
theorem mixedOnePoint_coupledFrom_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQFrom (k := k) X)
        (N6MixedOnePoint.coupledRFrom (k := k) X)
        1) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.coupledFrom_pointwise_iff_exists_outer (k := k) X

/-- Every element of the intrinsic mixed coupled cell has determinant `1`. -/
theorem mixedOnePoint_coupledFrom_det
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) :
    Matrix.det
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQFrom (k := k) X)
        (N6MixedOnePoint.coupledRFrom (k := k) X)
        1) = 1 :=
  N6MixedOnePoint.coupledFrom_det (k := k) X

/-- Rank-one `2 × 2` matrices give a clean four-parameter pointwise family on the
intrinsic mixed coupled cell. -/
theorem mixedOnePoint_coupledOuter_pointwise
    (u v : N6MixedOnePoint.I → k) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    N6MixedOnePoint.FixesPairBivector
      (Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQFrom (k := k) X)
        (N6MixedOnePoint.coupledRFrom (k := k) X)
        1) :=
  N6MixedOnePoint.coupledOuter_pointwise (k := k) u v

/-- On the natural coupled-times-Levi cell for the mixed one-point orbit, pointwise
fixing is still controlled exactly by the coupled block equations. -/
theorem mixedOnePoint_coupledLevi_iff
    (Q : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k)
    (R : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    N6MixedOnePoint.FixesPairBivector
      ((Matrix.fromBlocks 1 Q R 1) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) ↔
      Q * N4.J * Qᵀ = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k) ∧
        N4.onePointRep₁ * Rᵀ + Q * N4.J = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k) ∧
        R * N4.onePointRep₁ + N4.J * Qᵀ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k) ∧
        R * N4.onePointRep₁ * Rᵀ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k) ∧
        R * N4.onePointRep₂ = (0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.W k) ∧
        N4.onePointRep₂ * Rᵀ = (0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.I k) :=
  N6MixedOnePoint.fixesPair_coupledLevi_iff
    (k := k) (Q := Q) (R := R) (A := A) (B := B) (H := H) hA hB hH

/-- On the coupled-from-times-Levi cell, pointwise fixing is exactly the
determinant-zero condition on the attached mixed coupled block. -/
theorem mixedOnePoint_coupledFromLevi_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    N6MixedOnePoint.FixesPairBivector
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) ↔
      X.det = 0 :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_det_zero
    (k := k) (X := X) (A := A) (B := B) (H := H) hA hB hH

/-- On the coupled-from-times-Levi cell, pointwise fixing is equivalent to the
explicit rank-one factorization `X = u vᵀ`. -/
theorem mixedOnePoint_coupledFromLevi_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    N6MixedOnePoint.FixesPairBivector
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_exists_outer
    (k := k) (X := X) (A := A) (B := B) (H := H) hA hB hH

/-- The rank-one mixed coupled family remains pointwise after multiplication by a
pointwise Levi element. -/
theorem mixedOnePoint_coupledOuterLevi_product_pointwise
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    N6MixedOnePoint.FixesPairBivector
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) :=
  N6MixedOnePoint.coupledOuterLevi_product_pointwise (k := k) u v A B H hA hB hH

/-- The rank-one intrinsic mixed coupled family has determinant `1` after
multiplication by a pointwise Levi element. -/
theorem mixedOnePoint_coupledOuterLevi_product_det
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hH : H.det = 1) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    Matrix.det
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) = 1 :=
  N6MixedOnePoint.coupledOuterLevi_product_det (k := k) u v A B H hA hH

/-- Multiplying an intrinsic mixed coupled element by a pointwise Levi element does not
change its determinant. -/
theorem mixedOnePoint_coupledFromLevi_product_det
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hH : H.det = 1) :
    Matrix.det
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) = 1 :=
  N6MixedOnePoint.coupledFromLevi_product_det (k := k) X A B H hA hH

/-- On the intrinsic coupled-from-times-Levi cell, left-multiplying the standard Borel
lift gives the expected quotient action whenever `det X = 0`. -/
theorem mixedOnePoint_coupledFromLevi_borel_product_lift_action
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (hX : X.det = 0)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.coupledFromLevi_borel_product_lift_action
    (k := k) X A B H hA hB hH hX a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the standard Borel product has the
expected quotient action exactly when `det X = 0`. -/
theorem mixedOnePoint_coupledFromLevi_borel_product_lift_action_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      X.det = 0 :=
  N6MixedOnePoint.coupledFromLevi_borel_product_lift_action_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the standard Borel product has the
expected quotient action exactly when the coupled block factors as `X = u vᵀ`. -/
theorem mixedOnePoint_coupledFromLevi_borel_product_lift_action_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.coupledFromLevi_borel_product_lift_action_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to the
usual quotient action after left multiplication by the standard mixed Borel lift. -/
theorem mixedOnePoint_coupledFromLevi_iff_left_borel_action
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    N6MixedOnePoint.FixesPairBivector gU ↔
      (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_left_borel_action
    (k := k) X A B H hA hB hH a b ha

/-- The concrete rank-one intrinsic mixed coupled family combines with the standard
Borel lift with the expected quotient action. -/
theorem mixedOnePoint_coupledOuterLevi_borel_product_lift_action
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.coupledOuterLevi_borel_product_lift_action
    (k := k) u v A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, right-multiplying the standard Borel
lift gives the expected quotient action whenever `det X = 0`. -/
theorem mixedOnePoint_borel_coupledFromLevi_right_product_lift_action
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (hX : X.det = 0)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.borel_coupledFromLevi_right_product_lift_action
    (k := k) X A B H hA hB hH hX a b ha

/-- On the intrinsic coupled-from-times-Levi cell, right-multiplying the standard Borel
lift has the expected quotient action exactly when `det X = 0`. -/
theorem mixedOnePoint_borel_coupledFromLevi_right_product_lift_action_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      X.det = 0 :=
  N6MixedOnePoint.borel_coupledFromLevi_right_product_lift_action_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, right-multiplying the standard Borel
lift has the expected quotient action exactly when the coupled block factors as
`X = u vᵀ`. -/
theorem mixedOnePoint_borel_coupledFromLevi_right_product_lift_action_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.borel_coupledFromLevi_right_product_lift_action_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to the
usual quotient action after right multiplication by the standard mixed Borel lift. -/
theorem mixedOnePoint_coupledFromLevi_iff_right_borel_action
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    N6MixedOnePoint.FixesPairBivector gU ↔
      (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_right_borel_action
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left quotient-action
condition is equivalent to `det X = 0`. -/
theorem mixedOnePoint_left_borel_action_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      X.det = 0 :=
  N6MixedOnePoint.coupledFromLevi_borel_product_lift_action_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left quotient-action
condition is equivalent to explicit rank-one factorization. -/
theorem mixedOnePoint_left_borel_action_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.coupledFromLevi_borel_product_lift_action_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right quotient-action
condition is equivalent to `det X = 0`. -/
theorem mixedOnePoint_right_borel_action_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      X.det = 0 :=
  N6MixedOnePoint.borel_coupledFromLevi_right_product_lift_action_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right quotient-action
condition is equivalent to explicit rank-one factorization. -/
theorem mixedOnePoint_right_borel_action_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.borel_coupledFromLevi_right_product_lift_action_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left and right Borel
quotient-action conditions are equivalent. -/
theorem mixedOnePoint_left_borel_action_iff_right_borel_action
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) :=
  N6MixedOnePoint.coupledFromLevi_left_borel_action_iff_right_borel_action
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to
simultaneously getting the usual left and right quotient actions. -/
theorem mixedOnePoint_coupledFromLevi_iff_left_and_right_borel_action
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    N6MixedOnePoint.FixesPairBivector gU ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂)) :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_left_and_right_borel_action
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, simultaneously getting the usual
left and right quotient actions is equivalent to the coupled block satisfying
`det X = 0`. -/
theorem mixedOnePoint_left_and_right_borel_action_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂)) ↔
      X.det = 0 :=
  N6MixedOnePoint.coupledFromLevi_left_and_right_borel_action_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, simultaneously getting the usual
left and right quotient actions is equivalent to explicit rank-one factorization. -/
theorem mixedOnePoint_left_and_right_borel_action_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂)) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.coupledFromLevi_left_and_right_borel_action_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left one-sided
Borel action-and-determinant package is already equivalent to the full
two-sided package. -/
theorem mixedOnePoint_left_borel_package_iff_left_and_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledFromLevi_left_borel_package_iff_left_and_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right one-sided
Borel action-and-determinant package is already equivalent to the full
two-sided package. -/
theorem mixedOnePoint_right_borel_package_iff_left_and_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) ↔
      (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledFromLevi_right_borel_package_iff_left_and_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to the
full two-sided Borel action-and-determinant package. -/
theorem mixedOnePoint_coupledFromLevi_iff_left_and_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    N6MixedOnePoint.FixesPairBivector gU ↔
      (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_left_and_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the full two-sided Borel
action-and-determinant package is equivalent to the single condition `det X = 0`. -/
theorem mixedOnePoint_left_and_right_borel_package_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) ↔
      X.det = 0 :=
  N6MixedOnePoint.coupledFromLevi_left_and_right_borel_package_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the full two-sided Borel
action-and-determinant package is equivalent to explicit rank-one factorization. -/
theorem mixedOnePoint_left_and_right_borel_package_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.coupledFromLevi_left_and_right_borel_package_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient
action is already equivalent to the full two-sided Borel action-and-determinant
package. -/
theorem mixedOnePoint_left_borel_action_iff_left_and_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledFromLevi_left_borel_action_iff_left_and_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient
action is already equivalent to the full two-sided Borel action-and-determinant
package. -/
theorem mixedOnePoint_right_borel_action_iff_left_and_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
              a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
            N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
              (a * a) • N6MixedOnePoint.rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledFromLevi_right_borel_action_iff_left_and_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient
action is already equivalent to the matching one-sided Borel
action-and-determinant package. -/
theorem mixedOnePoint_left_borel_action_iff_left_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledFromLevi_left_borel_action_iff_left_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient
action is already equivalent to the matching one-sided Borel
action-and-determinant package. -/
theorem mixedOnePoint_right_borel_action_iff_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledFromLevi_right_borel_action_iff_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient
action is already equivalent to the opposite one-sided Borel
action-and-determinant package. -/
theorem mixedOnePoint_left_borel_action_iff_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledFromLevi_left_borel_action_iff_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient
action is already equivalent to the opposite one-sided Borel
action-and-determinant package. -/
theorem mixedOnePoint_right_borel_action_iff_left_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂) ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledFromLevi_right_borel_action_iff_left_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient
action together with its determinant condition is equivalent to `det X = 0`. -/
theorem mixedOnePoint_left_borel_package_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      X.det = 0 :=
  N6MixedOnePoint.coupledFromLevi_left_borel_package_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient
action together with its determinant condition is equivalent to rank-one
factorization. -/
theorem mixedOnePoint_left_borel_package_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.coupledFromLevi_left_borel_package_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient
action together with its determinant condition is equivalent to `det X = 0`. -/
theorem mixedOnePoint_right_borel_package_iff_det_zero
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) ↔
      X.det = 0 :=
  N6MixedOnePoint.coupledFromLevi_right_borel_package_iff_det_zero
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient
action together with its determinant condition is equivalent to rank-one
factorization. -/
theorem mixedOnePoint_right_borel_package_iff_exists_outer
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) ↔
      ∃ u v : N6MixedOnePoint.I → k, X = fun i j => u i * v j :=
  N6MixedOnePoint.coupledFromLevi_right_borel_package_iff_exists_outer
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent
to the usual left Borel quotient action together with its determinant condition. -/
theorem mixedOnePoint_fixesPair_iff_left_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    N6MixedOnePoint.FixesPairBivector gU ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_left_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent
to the usual right Borel quotient action together with its determinant condition. -/
theorem mixedOnePoint_fixesPair_iff_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    N6MixedOnePoint.FixesPairBivector gU ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.fixesPair_coupledFromLevi_iff_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- On the intrinsic coupled-from-times-Levi cell, the one-sided left and right
Borel action-and-determinant packages are equivalent. -/
theorem mixedOnePoint_left_borel_package_iff_right_borel_package
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledFromLevi_left_borel_package_iff_right_borel_package
    (k := k) X A B H hA hB hH a b ha

/-- The concrete rank-one intrinsic mixed coupled family also combines with the
standard Borel lift on the right with the expected quotient action. -/
theorem mixedOnePoint_borel_coupledOuterLevi_right_product_lift_action
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.borel_coupledOuterLevi_right_product_lift_action
    (k := k) u v A B H hA hB hH a b ha

/-- Left-multiplying the standard mixed Borel lift by an intrinsic
coupled-from-times-Levi element keeps determinant `a^3`. -/
theorem mixedOnePoint_coupledFromLevi_borel_product_lift_det
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a :=
  N6MixedOnePoint.coupledFromLevi_borel_product_lift_det
    (k := k) X A B H hA hH a b ha

/-- The concrete rank-one intrinsic mixed coupled family keeps determinant `a^3` after
left multiplication by the standard Borel lift. -/
theorem mixedOnePoint_coupledOuterLevi_borel_product_lift_det
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a :=
  N6MixedOnePoint.coupledOuterLevi_borel_product_lift_det
    (k := k) u v A B H hA hH a b ha

/-- Right-multiplying the standard mixed Borel lift by an intrinsic
coupled-from-times-Levi element keeps determinant `a^3`. -/
theorem mixedOnePoint_borel_coupledFromLevi_right_product_lift_det
    (X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    Matrix.det (gL * gU) = a * a * a :=
  N6MixedOnePoint.borel_coupledFromLevi_right_product_lift_det
    (k := k) X A B H hA hH a b ha

/-- The concrete rank-one intrinsic mixed coupled family keeps determinant `a^3` after
right multiplication by the standard Borel lift. -/
theorem mixedOnePoint_borel_coupledOuterLevi_right_product_lift_det
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    Matrix.det (gL * gU) = a * a * a :=
  N6MixedOnePoint.borel_coupledOuterLevi_right_product_lift_det
    (k := k) u v A B H hA hH a b ha

/-- The concrete rank-one intrinsic mixed coupled-from-times-Levi family gives the
expected left Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledOuterLevi_left_borel_package
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledOuterLevi_left_borel_package
    (k := k) u v A B H hA hB hH a b ha

/-- The concrete rank-one intrinsic mixed coupled-from-times-Levi family gives the
expected right Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledOuterLevi_right_borel_package
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledOuterLevi_right_borel_package
    (k := k) u v A B H hA hB hH a b ha

/-- The concrete rank-one intrinsic coupled-from-times-Levi family gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem mixedOnePoint_coupledOuterLevi_left_and_right_borel_package
    (u v : N6MixedOnePoint.I → k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := fun i j => u i * v j
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQFrom (k := k) X)
          (N6MixedOnePoint.coupledRFrom (k := k) X)
          1) *
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H))
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledOuterLevi_left_and_right_borel_package
    (k := k) u v A B H hA hB hH a b ha

/-- The explicit mixed pointwise subgroup combines with the standard Borel lift with the
expected action on the chosen basis. -/
theorem mixedOnePoint_coupledLevi_borel_product_lift_action
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) *
        (Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.coupledLevi_borel_product_lift_action
    (k := k)
    (s := s)
    (t := t)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH
    (a := a)
    (b := b)
    ha

/-- The diagonal explicit mixed coupled subgroup also combines with the standard Borel
lift with the expected action on the chosen basis. -/
theorem mixedOnePoint_coupledDiagonalLevi_borel_product_lift_action
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) s)
        (N6MixedOnePoint.coupledR11 (k := k) s)
        1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t)
        (N6MixedOnePoint.coupledR22 (k := k) t)
        1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.coupledDiagonalLevi_borel_product_lift_action
    (k := k)
    (s := s)
    (t := t)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH
    (a := a)
    (b := b)
    ha

/-- The explicit off-diagonal mixed coupled-plus-Levi subgroup gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledLevi_left_borel_package
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          1) *
       (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          1)) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledLevi_left_borel_package
    (k := k) s t A B H hA hB hH a b ha

/-- The diagonal mixed coupled-plus-Levi subgroup gives the expected left Borel
quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledDiagonalLevi_left_borel_package
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledDiagonalLevi_left_borel_package
    (k := k) s t A B H hA hB hH a b ha

/-- The full explicit four-parameter mixed pointwise subgroup also combines with the
standard Borel lift with the expected action on the chosen basis. -/
theorem mixedOnePoint_coupledExplicitLevi_borel_product_lift_action
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ12 (k := k) s₁)
        (N6MixedOnePoint.coupledR12 (k := k) s₁)
        1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ21 (k := k) t₁)
        (N6MixedOnePoint.coupledR21 (k := k) t₁)
        1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) s₂)
        (N6MixedOnePoint.coupledR11 (k := k) s₂)
        1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t₂)
        (N6MixedOnePoint.coupledR22 (k := k) t₂)
        1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.coupledExplicitLevi_borel_product_lift_action
    (k := k)
    (s₁ := s₁)
    (t₁ := t₁)
    (s₂ := s₂)
    (t₂ := t₂)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH
    (a := a)
    (b := b)
    ha

/-- Right-multiplying the standard Borel lift by the explicit mixed pointwise subgroup
does not change the quotient action. -/
theorem mixedOnePoint_borel_coupledLevi_right_product_lift_action
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) *
       (Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.borel_coupledLevi_right_product_lift_action
    (k := k)
    (s := s)
    (t := t)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH
    (a := a)
    (b := b)
    ha

/-- Right-multiplying the standard Borel lift by the diagonal explicit mixed coupled
subgroup does not change the quotient action. -/
theorem mixedOnePoint_borel_coupledDiagonalLevi_right_product_lift_action
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) s)
        (N6MixedOnePoint.coupledR11 (k := k) s)
        1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t)
        (N6MixedOnePoint.coupledR22 (k := k) t)
        1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.borel_coupledDiagonalLevi_right_product_lift_action
    (k := k)
    (s := s)
    (t := t)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH
    (a := a)
    (b := b)
    ha

/-- The explicit off-diagonal mixed coupled-plus-Levi subgroup gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledLevi_right_borel_package
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          1) *
       (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          1)) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledLevi_right_borel_package
    (k := k) s t A B H hA hB hH a b ha

/-- The diagonal mixed coupled-plus-Levi subgroup gives the expected right Borel
quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledDiagonalLevi_right_borel_package
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledDiagonalLevi_right_borel_package
    (k := k) s t A B H hA hB hH a b ha

/-- The explicit off-diagonal mixed coupled-plus-Levi subgroup gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem mixedOnePoint_coupledLevi_left_and_right_borel_package
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          1) *
       (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          1)) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledLevi_left_and_right_borel_package
    (k := k) s t A B H hA hB hH a b ha

/-- The diagonal mixed coupled-plus-Levi subgroup gives the expected quotient action
on both sides of the standard Borel lift, and both products have determinant `a^3`. -/
theorem mixedOnePoint_coupledDiagonalLevi_left_and_right_borel_package
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledDiagonalLevi_left_and_right_borel_package
    (k := k) s t A B H hA hB hH a b ha

/-- The single mixed `E12` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem mixedOnePoint_coupledE12Levi_left_and_right_borel_package
    (s : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s) (N6MixedOnePoint.coupledR12 (k := k) s) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g12 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledE12Levi_left_and_right_borel_package
    (k := k) s A B H hA hB hH a b ha

/-- The single mixed `E12` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE12Levi_left_borel_package
    (s : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s) (N6MixedOnePoint.coupledR12 (k := k) s) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g12 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledE12Levi_left_borel_package
    (k := k) s A B H hA hB hH a b ha

/-- The single mixed `E12` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE12Levi_right_borel_package
    (s : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s) (N6MixedOnePoint.coupledR12 (k := k) s) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g12 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledE12Levi_right_borel_package
    (k := k) s A B H hA hB hH a b ha

/-- The single mixed `E21` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem mixedOnePoint_coupledE21Levi_left_and_right_borel_package
    (t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t) (N6MixedOnePoint.coupledR21 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g21 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledE21Levi_left_and_right_borel_package
    (k := k) t A B H hA hB hH a b ha

/-- The single mixed `E21` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE21Levi_left_borel_package
    (t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t) (N6MixedOnePoint.coupledR21 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g21 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledE21Levi_left_borel_package
    (k := k) t A B H hA hB hH a b ha

/-- The single mixed `E21` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE21Levi_right_borel_package
    (t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t) (N6MixedOnePoint.coupledR21 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g21 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledE21Levi_right_borel_package
    (k := k) t A B H hA hB hH a b ha

/-- The single mixed `E11` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem mixedOnePoint_coupledE11Levi_left_and_right_borel_package
    (s : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g11 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledE11Levi_left_and_right_borel_package
    (k := k) s A B H hA hB hH a b ha

/-- The single mixed `E11` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE11Levi_left_borel_package
    (s : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g11 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledE11Levi_left_borel_package
    (k := k) s A B H hA hB hH a b ha

/-- The single mixed `E11` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE11Levi_right_borel_package
    (s : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g11 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledE11Levi_right_borel_package
    (k := k) s A B H hA hB hH a b ha

/-- The single mixed `E22` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem mixedOnePoint_coupledE22Levi_left_and_right_borel_package
    (t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g22 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    (((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) :=
  N6MixedOnePoint.coupledE22Levi_left_and_right_borel_package
    (k := k) t A B H hA hB hH a b ha

/-- The single mixed `E22` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE22Levi_left_borel_package
    (t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g22 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledE22Levi_left_borel_package
    (k := k) t A B H hA hB hH a b ha

/-- The single mixed `E22` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledE22Levi_right_borel_package
    (t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      g22 * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
            a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
          N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
            (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledE22Levi_right_borel_package
    (k := k) t A B H hA hB hH a b ha

/-- Right-multiplying the standard Borel lift by the full explicit four-parameter mixed
pointwise subgroup does not change the quotient action. -/
theorem mixedOnePoint_borel_coupledExplicitLevi_right_product_lift_action
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ12 (k := k) s₁)
        (N6MixedOnePoint.coupledR12 (k := k) s₁)
        1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ21 (k := k) t₁)
        (N6MixedOnePoint.coupledR21 (k := k) t₁)
        1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ11 (k := k) s₂)
        (N6MixedOnePoint.coupledR11 (k := k) s₂)
        1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        1
        (N6MixedOnePoint.coupledQ22 (k := k) t₂)
        (N6MixedOnePoint.coupledR22 (k := k) t₂)
        1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
        a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
      N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
        (a * a) • N6MixedOnePoint.rep₂ :=
  N6MixedOnePoint.borel_coupledExplicitLevi_right_product_lift_action
    (k := k)
    (s₁ := s₁)
    (t₁ := t₁)
    (s₂ := s₂)
    (t₂ := t₂)
    (A := A)
    (B := B)
    (H := H)
    hA
    hB
    hH
    (a := a)
    (b := b)
    ha

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected
quotient action on both sides of the standard Borel lift. -/
theorem mixedOnePoint_coupledExplicitLevi_left_and_right_borel_action
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s₁) (N6MixedOnePoint.coupledR12 (k := k) s₁) 1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t₁) (N6MixedOnePoint.coupledR21 (k := k) t₁) 1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s₂) (N6MixedOnePoint.coupledR11 (k := k) s₂) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t₂) (N6MixedOnePoint.coupledR22 (k := k) t₂) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      (N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂)) :=
  N6MixedOnePoint.coupledExplicitLevi_left_and_right_borel_action
    (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH a b ha

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledExplicitLevi_left_borel_package
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s₁) (N6MixedOnePoint.coupledR12 (k := k) s₁) 1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t₁) (N6MixedOnePoint.coupledR21 (k := k) t₁) 1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s₂) (N6MixedOnePoint.coupledR11 (k := k) s₂) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t₂) (N6MixedOnePoint.coupledR22 (k := k) t₂) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) :=
  N6MixedOnePoint.coupledExplicitLevi_left_borel_package
    (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH a b ha

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem mixedOnePoint_coupledExplicitLevi_right_borel_package
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s₁) (N6MixedOnePoint.coupledR12 (k := k) s₁) 1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t₁) (N6MixedOnePoint.coupledR21 (k := k) t₁) 1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s₂) (N6MixedOnePoint.coupledR11 (k := k) s₂) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t₂) (N6MixedOnePoint.coupledR22 (k := k) t₂) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledExplicitLevi_right_borel_package
    (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH a b ha

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem mixedOnePoint_coupledExplicitLevi_left_and_right_borel_package
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s₁) (N6MixedOnePoint.coupledR12 (k := k) s₁) 1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t₁) (N6MixedOnePoint.coupledR21 (k := k) t₁) 1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s₂) (N6MixedOnePoint.coupledR11 (k := k) s₂) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t₂) (N6MixedOnePoint.coupledR22 (k := k) t₂) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gU * gL) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gU * gL) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₁ (gL * gU) =
          a • N6MixedOnePoint.rep₁ + b • N6MixedOnePoint.rep₂ ∧
        N6MixedOnePoint.ActBivector N6MixedOnePoint.rep₂ (gL * gU) =
          (a * a) • N6MixedOnePoint.rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) :=
  N6MixedOnePoint.coupledExplicitLevi_left_and_right_borel_package
    (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH a b ha

/-- The mixed Levi pointwise family has determinant `1`. -/
theorem mixedOnePoint_pointwise_levi_det
    (A B H : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hH : H.det = 1) :
    Matrix.det (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H) = 1 :=
  N6MixedOnePoint.pointwise_levi_det (k := k) A B H hA hH

/-- The mixed coupled-plus-Levi pointwise subgroup has determinant `1`. -/
theorem mixedOnePoint_coupledLevi_product_det
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) *
       (Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    Matrix.det gU = 1 :=
  N6MixedOnePoint.coupledLevi_product_det (k := k) s t A B H hA hB hH

/-- The diagonal mixed coupled-plus-Levi pointwise subgroup has determinant `1`. -/
theorem mixedOnePoint_coupledDiagonalLevi_product_det
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    Matrix.det gU = 1 :=
  N6MixedOnePoint.coupledDiagonalLevi_product_det (k := k) s t A B H hA hB hH

/-- The full explicit four-parameter mixed coupled-plus-Levi pointwise subgroup has
determinant `1`. -/
theorem mixedOnePoint_coupledExplicitLevi_product_det
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s₁) (N6MixedOnePoint.coupledR12 (k := k) s₁) 1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t₁) (N6MixedOnePoint.coupledR21 (k := k) t₁) 1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s₂) (N6MixedOnePoint.coupledR11 (k := k) s₂) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t₂) (N6MixedOnePoint.coupledR22 (k := k) t₂) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    Matrix.det gU = 1 :=
  N6MixedOnePoint.coupledExplicitLevi_product_det (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH

/-- Left-multiplying the mixed Borel lift by the mixed coupled-plus-Levi pointwise
subgroup does not change the determinant. -/
theorem mixedOnePoint_coupledLevi_borel_product_lift_det
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) *
       (Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a :=
  N6MixedOnePoint.coupledLevi_borel_product_lift_det
    (k := k) s t A B H hA hB hH a b ha

/-- Left-multiplying the mixed Borel lift by the diagonal mixed coupled-plus-Levi
pointwise subgroup does not change the determinant. -/
theorem mixedOnePoint_coupledDiagonalLevi_borel_product_lift_det
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a :=
  N6MixedOnePoint.coupledDiagonalLevi_borel_product_lift_det
    (k := k) s t A B H hA hB hH a b ha

/-- Left-multiplying the mixed Borel lift by the full explicit four-parameter mixed
coupled-plus-Levi pointwise subgroup does not change the determinant. -/
theorem mixedOnePoint_coupledExplicitLevi_borel_product_lift_det
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s₁) (N6MixedOnePoint.coupledR12 (k := k) s₁) 1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t₁) (N6MixedOnePoint.coupledR21 (k := k) t₁) 1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s₂) (N6MixedOnePoint.coupledR11 (k := k) s₂) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t₂) (N6MixedOnePoint.coupledR22 (k := k) t₂) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k := !![a, 0; 0, 1]
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a :=
  N6MixedOnePoint.coupledExplicitLevi_borel_product_lift_det
    (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH a b ha

/-- Right-multiplying the mixed Borel lift by the mixed coupled-plus-Levi pointwise
subgroup does not change the determinant. -/
theorem mixedOnePoint_borel_coupledLevi_right_product_lift_det
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ12 (k := k) s)
          (N6MixedOnePoint.coupledR12 (k := k) s)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)) *
       (Matrix.fromBlocks
          (1 : Matrix N6MixedOnePoint.W N6MixedOnePoint.W k)
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          (1 : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) *
      (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    Matrix.det (gL * gU) = a * a * a :=
  N6MixedOnePoint.borel_coupledLevi_right_product_lift_det
    (k := k) s t A B H hA hB hH a b ha

/-- Right-multiplying the mixed Borel lift by the diagonal mixed coupled-plus-Levi
pointwise subgroup does not change the determinant. -/
theorem mixedOnePoint_borel_coupledDiagonalLevi_right_product_lift_det
    (s t : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s) (N6MixedOnePoint.coupledR11 (k := k) s) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t) (N6MixedOnePoint.coupledR22 (k := k) t) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      (g11 * g22) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    Matrix.det (gL * gU) = a * a * a :=
  N6MixedOnePoint.borel_coupledDiagonalLevi_right_product_lift_det
    (k := k) s t A B H hA hB hH a b ha

/-- Right-multiplying the mixed Borel lift by the full explicit four-parameter mixed
coupled-plus-Levi pointwise subgroup does not change the determinant. -/
theorem mixedOnePoint_borel_coupledExplicitLevi_right_product_lift_det
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g12 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ12 (k := k) s₁) (N6MixedOnePoint.coupledR12 (k := k) s₁) 1
    let g21 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ21 (k := k) t₁) (N6MixedOnePoint.coupledR21 (k := k) t₁) 1
    let g11 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ11 (k := k) s₂) (N6MixedOnePoint.coupledR11 (k := k) s₂) 1
    let g22 : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      Matrix.fromBlocks 1 (N6MixedOnePoint.coupledQ22 (k := k) t₂) (N6MixedOnePoint.coupledR22 (k := k) t₂) 1
    let gU : Matrix N6MixedOnePoint.V N6MixedOnePoint.V k :=
      ((g12 * g21) * (g11 * g22)) * (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
    Matrix.det (gL * gU) = a * a * a :=
  N6MixedOnePoint.borel_coupledExplicitLevi_right_product_lift_det
    (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH a b ha

/-- The direct-sum `[a]` family in dimension `6` admits the expected radical-extension
pointwise family. -/
theorem simplePoint_pointwise_family
    (u a : k)
    (ha : a ≠ 0)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    N6SimplePoint.FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a)
          0
          0
          E)) :=
  N6SimplePoint.radSimple_pointwise_family (k := k) (u := u) (a := a) ha C E hE

/-- The explicit `4`-parameter unipotent family on the `5`-dimensional `[a]` quotient
also extends through the radical direction. -/
theorem simplePoint_pointwise_unipotent_family
    (u : k)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (x y p q : k) :
    N6SimplePoint.FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)) :=
  N6SimplePoint.radSimple_pointwise_unipotent_family
    (k := k) (u := u) (C := C) (x := x) (y := y) (p := p) (q := q)

/-- Combining the lifted quotient-unipotent family with the Levi family gives a
concrete pointwise subgroup of the direct-sum `[a]` stabilizer in dimension `6`. -/
theorem simplePoint_pointwise_product_family
    (u : k)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (x y p q a : k)
    (ha : a ≠ 0)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    let gU : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)
    let gL : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) (1 : k))
        0
        0
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a)
          0
          0
          E)
    N6SimplePoint.FixesRadSimplePairBivector (gU * gL) :=
  N6SimplePoint.radSimple_pointwise_product_family
    (k := k)
    (u := u)
    (C := C)
    (x := x)
    (y := y)
    (p := p)
    (q := q)
    (a := a)
    ha
    (E := E)
    hE

/-- The direct-sum `[a]` family in dimension `6` also admits the expected explicit
general pointwise shape on the `5`-dimensional quotient, together with an arbitrary
radical column. -/
theorem simplePoint_pointwise_shape_family
    (u : k)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (a x y p q : k)
    (ha : a ≠ 0)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    let B : Matrix N5SimplePoint.W N5SimplePoint.I k := !![(0 : k), 0; 0, 0; p, q]
    let C0 : Matrix N5SimplePoint.I N5SimplePoint.W k :=
      !![a * (p * E 0 1 - q * E 0 0), 0, 0;
         a * (p * E 1 1 - q * E 1 0), 0, 0]
    N6SimplePoint.FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E)) :=
  N6SimplePoint.radSimple_pointwise_shape_family
    (k := k)
    (u := u)
    (C := C)
    (a := a)
    (x := x)
    (y := y)
    (p := p)
    (q := q)
    ha
    (E := E)
    hE

/-- If the quotient block is already in the pure singular shape, then pointwise fixing
forces the remaining quotient blocks to have the expected direct-sum `[a]` form in
dimension `6`, while the radical scalar and radical column remain free. -/
theorem simplePoint_nested_pureShape_iff
    (u : k)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (a x y : k)
    (B : Matrix N5SimplePoint.W N5SimplePoint.I k)
    (C0 : Matrix N5SimplePoint.I N5SimplePoint.W k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (ha : a ≠ 0) :
    N6SimplePoint.FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E)) ↔
      E.det = 1 ∧
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E =
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          (!![(0 : k), 0; 0, 0; B 2 0, B 2 1])
          (!![a * (B 2 0 * E 0 1 - B 2 1 * E 0 0), 0, 0;
             a * (B 2 0 * E 1 1 - B 2 1 * E 1 0), 0, 0])
          E :=
  N6SimplePoint.fixesRadSimplePair_fromBlocks_nested_pureShape_iff
    (k := k)
    (u := u)
    (C := C)
    (a := a)
    (x := x)
    (y := y)
    (B := B)
    (C0 := C0)
    (E := E)
    ha

/-- For the direct-sum `[a]` row in dimension `6`, pointwise fixing also forces the
upper-right radical row to vanish. After that, the lower-right block has the expected
nested shape coming from the `5`-dimensional quotient. -/
theorem simplePoint_nested_iff_shape
    (u : k)
    (B0 : Matrix N6SimplePoint.I N6SimplePoint.W k)
    (C0 : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (a x y : k)
    (B : Matrix N5SimplePoint.W N5SimplePoint.I k)
    (C : Matrix N5SimplePoint.I N5SimplePoint.W k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (ha : a ≠ 0) :
    N6SimplePoint.FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        B0
        C0
        (Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C
          E)) ↔
      B0 = 0 ∧
        E.det = 1 ∧
          Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            B
            C
            E =
          Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            (!![(0 : k), 0; 0, 0; B 2 0, B 2 1])
            (!![a * (B 2 0 * E 0 1 - B 2 1 * E 0 0), 0, 0;
               a * (B 2 0 * E 1 1 - B 2 1 * E 1 0), 0, 0])
            E :=
  N6SimplePoint.fixesRadSimplePair_fromBlocks_nested_iff_shape
    (k := k)
    (u := u)
    (B0 := B0)
    (C0 := C0)
    (a := a)
    (x := x)
    (y := y)
    (B := B)
    (C := C)
    (E := E)
    ha

/-- A concrete Borel lift on the direct-sum `[a]` family in dimension `6`. -/
theorem simplePoint_borel_lift_action
    (u a b : k)
    (ha : a ≠ 0)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix N6SimplePoint.W N6SimplePoint.W k := Matrix.fromBlocks G 0 0 E
    let g : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      Matrix.fromBlocks (N6SimplePoint.scalarBlock (k := k) u) 0 C h
    N6SimplePoint.ActBivector N6SimplePoint.radSimpleRep₁ g =
        a • N6SimplePoint.radSimpleRep₁ + b • N6SimplePoint.radSimpleRep₂ ∧
      N6SimplePoint.ActBivector N6SimplePoint.radSimpleRep₂ g =
        (a * a) • N6SimplePoint.radSimpleRep₂ :=
  N6SimplePoint.radSimple_borel_lift_action (k := k) (u := u) (a := a) (b := b) ha C

/-- The determinant of the standard Borel lift on the direct-sum `[a]` family in
dimension `6` is the radical scalar times `a^4`. -/
theorem simplePoint_borel_lift_det
    (u a b : k)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix N6SimplePoint.W N6SimplePoint.W k := Matrix.fromBlocks G 0 0 E
    let g : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      Matrix.fromBlocks (N6SimplePoint.scalarBlock (k := k) u) 0 C h
    Matrix.det g = u * ((a * a) * (a * a)) :=
  N6SimplePoint.radSimple_borel_lift_det (k := k) (u := u) (a := a) (b := b) C

/-- The explicit pointwise subgroup on the direct-sum `[a]` family combines with
the standard Borel lift exactly as expected. -/
theorem simplePoint_pointwise_borel_product_lift_action
    (u : k)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (x y p q a₀ : k)
    (ha₀ : a₀ ≠ 0)
    (E0 : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE0 : E0.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)) *
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) (1 : k))
        0
        0
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a₀)
          0
          0
          E0))
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix N6SimplePoint.W N6SimplePoint.W k := Matrix.fromBlocks G 0 0 E
    let gL : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      Matrix.fromBlocks (N6SimplePoint.scalarBlock (k := k) 1) 0 0 h
    N6SimplePoint.ActBivector N6SimplePoint.radSimpleRep₁ (gU * gL) =
        a • N6SimplePoint.radSimpleRep₁ + b • N6SimplePoint.radSimpleRep₂ ∧
      N6SimplePoint.ActBivector N6SimplePoint.radSimpleRep₂ (gU * gL) =
        (a * a) • N6SimplePoint.radSimpleRep₂ :=
  N6SimplePoint.radSimple_pointwise_borel_product_lift_action
    (k := k)
    (u := u)
    (C := C)
    (x := x)
    (y := y)
    (p := p)
    (q := q)
    (a₀ := a₀)
    ha₀
    (E0 := E0)
    hE0
    (a := a)
    (b := b)
    ha

/-- Right-multiplying the standard Borel lift by the explicit pointwise subgroup on the
direct-sum `[a]` family does not change the quotient action. -/
theorem simplePoint_borel_pointwise_right_product_lift_action
    (u : k)
    (C : Matrix N6SimplePoint.W N6SimplePoint.I k)
    (x y p q a₀ : k)
    (ha₀ : a₀ ≠ 0)
    (E0 : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE0 : E0.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix N6SimplePoint.W N6SimplePoint.W k := Matrix.fromBlocks G 0 0 E
    let gL : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      Matrix.fromBlocks (N6SimplePoint.scalarBlock (k := k) 1) 0 0 h
    let gU : Matrix N6SimplePoint.V N6SimplePoint.V k :=
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)) *
      (Matrix.fromBlocks
        (N6SimplePoint.scalarBlock (k := k) (1 : k))
        0
        0
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a₀)
          0
          0
          E0))
    N6SimplePoint.ActBivector N6SimplePoint.radSimpleRep₁ (gL * gU) =
        a • N6SimplePoint.radSimpleRep₁ + b • N6SimplePoint.radSimpleRep₂ ∧
      N6SimplePoint.ActBivector N6SimplePoint.radSimpleRep₂ (gL * gU) =
        (a * a) • N6SimplePoint.radSimpleRep₂ :=
  N6SimplePoint.radSimple_borel_pointwise_right_product_lift_action
    (k := k)
    (u := u)
    (C := C)
    (x := x)
    (y := y)
    (p := p)
    (q := q)
    (a₀ := a₀)
    ha₀
    (E0 := E0)
    hE0
    (a := a)
    (b := b)
    ha

/-- Summary form of the three-point `n = 6` pointwise stabilizer calculation. -/
theorem threePoint_pointwise_bivector_iff
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (B : Matrix N6ThreePoint.W N6ThreePoint.I k)
    (C : Matrix N6ThreePoint.I N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k) :
    N6ThreePoint.FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ C = 0 ∧ N4.FixesSplitPairBivector A ∧ D.det = 1 :=
  N6ThreePoint.fixesPair_fromBlocks_iff (k := k) (A := A) (B := B) (C := C) (D := D)

/-- The obvious block-diagonal pointwise family on the three-point orbit. -/
theorem threePoint_pointwise_family
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    N6ThreePoint.FixesPairBivector (Matrix.fromBlocks A 0 0 D) :=
  N6ThreePoint.pointwise_family (k := k) (A := A) (D := D) hA hD

/-- The full block-diagonal pointwise family on the three-point orbit has determinant `1`. -/
theorem threePoint_pointwise_family_det
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    Matrix.det (Matrix.fromBlocks A 0 0 D) = 1 :=
  N6ThreePoint.pointwise_family_det (k := k) (A := A) (D := D) hA hD

/-- A split swap on the `4`-block together with a determinant `-1` sign flip on the
simple block preserves the three-point orbit with the expected action. -/
theorem threePoint_swap_lift_action :
    let E : Matrix N6ThreePoint.I N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
    let g : Matrix N6ThreePoint.V N6ThreePoint.V k :=
      Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    N6ThreePoint.ActBivector N6ThreePoint.rep₁ g = N6ThreePoint.rep₁ ∧
    N6ThreePoint.ActBivector N6ThreePoint.rep₂ g =
        N6ThreePoint.rep₁ - N6ThreePoint.rep₂ :=
  N6ThreePoint.swap_lift_action (k := k)

/-- The concrete three-point swap lift is an involution. -/
theorem threePoint_swap_lift_sq :
    let E : Matrix N6ThreePoint.I N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
    let g : Matrix N6ThreePoint.V N6ThreePoint.V k :=
      Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    g * g = 1 :=
  N6ThreePoint.swap_lift_sq (k := k)

/-- Left-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change the induced quotient action. -/
theorem threePoint_pointwise_swap_product_lift_action
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix N6ThreePoint.I N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
    let gU : Matrix N6ThreePoint.V N6ThreePoint.V k := Matrix.fromBlocks A 0 0 D
    let gL : Matrix N6ThreePoint.V N6ThreePoint.V k :=
      Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    N6ThreePoint.ActBivector N6ThreePoint.rep₁ (gU * gL) = N6ThreePoint.rep₁ ∧
      N6ThreePoint.ActBivector N6ThreePoint.rep₂ (gU * gL) =
        N6ThreePoint.rep₁ - N6ThreePoint.rep₂ :=
  N6ThreePoint.pointwise_swap_product_lift_action (k := k) (A := A) (D := D) hA hD

/-- Right-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change the induced quotient action. -/
theorem threePoint_swap_pointwise_right_product_lift_action
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix N6ThreePoint.I N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
    let gL : Matrix N6ThreePoint.V N6ThreePoint.V k :=
      Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    let gU : Matrix N6ThreePoint.V N6ThreePoint.V k := Matrix.fromBlocks A 0 0 D
    N6ThreePoint.ActBivector N6ThreePoint.rep₁ (gL * gU) = N6ThreePoint.rep₁ ∧
      N6ThreePoint.ActBivector N6ThreePoint.rep₂ (gL * gU) =
        N6ThreePoint.rep₁ - N6ThreePoint.rep₂ :=
  N6ThreePoint.swap_pointwise_right_product_lift_action
    (k := k) (A := A) (D := D) hA hD

/-- Left-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change its determinant. -/
theorem threePoint_pointwise_swap_product_lift_det
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix N6ThreePoint.I N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
    let gU : Matrix N6ThreePoint.V N6ThreePoint.V k := Matrix.fromBlocks A 0 0 D
    let gL : Matrix N6ThreePoint.V N6ThreePoint.V k :=
      Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    Matrix.det (gU * gL) = Matrix.det gL :=
  N6ThreePoint.pointwise_swap_product_lift_det
    (k := k) (A := A) (D := D) hA hD

/-- Right-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change its determinant. -/
theorem threePoint_swap_pointwise_right_product_lift_det
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix N6ThreePoint.I N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
    let gL : Matrix N6ThreePoint.V N6ThreePoint.V k :=
      Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    let gU : Matrix N6ThreePoint.V N6ThreePoint.V k := Matrix.fromBlocks A 0 0 D
    Matrix.det (gL * gU) = Matrix.det gL :=
  N6ThreePoint.swap_pointwise_right_product_lift_det
    (k := k) (A := A) (D := D) hA hD

/-- On the three-point line in dimension `6`, rank drop occurs exactly on the
three distinguished projective points. -/
theorem threePoint_det_zero_iff
    (a b : k) :
    Matrix.det
        (a • (N6ThreePoint.rep₁ (k := k)) + b • (N6ThreePoint.rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 ∨ b = 0 :=
  N6ThreePoint.det_zero_iff (k := k) (a := a) (b := b)

/-- On the direct-sum `2[a]+[b]` line in dimension `6`, rank drop occurs exactly on the
two distinguished projective points. -/
theorem onePointPlusSimple_det_zero_iff
    (a b : k) :
    Matrix.det
        (a • (N6OnePointPlusSimple.rep₁ (k := k)) + b • (N6OnePointPlusSimple.rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  N6OnePointPlusSimple.det_zero_iff (k := k) (a := a) (b := b)

/-- Sharp pointwise description for the direct-sum `2[a]+[b]` orbit. -/
theorem onePointPlusSimple_pointwise_bivector_iff
    (A : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.W k)
    (B : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.I k)
    (C : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.W k)
    (D : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k) :
    N6OnePointPlusSimple.FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ C = 0 ∧ N4.FixesOnePointPairBivector A ∧ D.det = 1 :=
  N6OnePointPlusSimple.fixesPair_fromBlocks_iff (k := k) (A := A) (B := B) (C := C) (D := D)

/-- The obvious block-diagonal pointwise family on the direct-sum `2[a]+[b]` orbit. -/
theorem onePointPlusSimple_pointwise_levi_family
    (A B C D E : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE : E.det = 1) :
    N6OnePointPlusSimple.FixesPairBivector
      (Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E) :=
  N6OnePointPlusSimple.pointwise_levi_family
    (k := k) (A := A) (B := B) (C := C) (D := D) (E := E) hA hC hD hrel hE

/-- The full block-diagonal pointwise Levi family on the direct-sum `2[a]+[b]` orbit
has determinant `1`. -/
theorem onePointPlusSimple_pointwise_levi_det
    (A B C D E : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE : E.det = 1) :
    Matrix.det (Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E) = 1 :=
  N6OnePointPlusSimple.pointwise_levi_det
    (k := k) (A := A) (B := B) (C := C) (D := D) (E := E) hA hC hD hrel hE

/-- A concrete torus-type lift on the direct-sum `2[a]+[b]` orbit. -/
theorem onePointPlusSimple_torus_lift_action
    (a : k)
    (ha : a ≠ 0) :
    let B : Matrix N4.I N4.I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k := !![a * a, 0; 0, 1]
    let g : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k := Matrix.fromBlocks H 0 0 E
    N6OnePointPlusSimple.ActBivector N6OnePointPlusSimple.rep₁ g =
        a • N6OnePointPlusSimple.rep₁ + (a * a - a) • N6OnePointPlusSimple.rep₂ ∧
      N6OnePointPlusSimple.ActBivector N6OnePointPlusSimple.rep₂ g =
        (a * a) • N6OnePointPlusSimple.rep₂ :=
  N6OnePointPlusSimple.torus_lift_action (k := k) (a := a) ha

/-- The determinant of the standard torus-type lift on the direct-sum
`2[a]+[b]` orbit is `a^4`. -/
theorem onePointPlusSimple_torus_lift_det
    (a : k) :
    let B : Matrix N4.I N4.I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k := !![a * a, 0; 0, 1]
    let g : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k := Matrix.fromBlocks H 0 0 E
    Matrix.det g = (a * a) * (a * a) :=
  N6OnePointPlusSimple.torus_lift_det (k := k) (a := a)

/-- The full block-diagonal pointwise subgroup on the direct-sum `2[a]+[b]` orbit
combines with the torus-type quotient lift with the expected action. -/
theorem onePointPlusSimple_pointwise_torus_product_lift_action
    (A B C D E0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k :=
      Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    let B0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k :=
      !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k := !![a * a, 0; 0, 1]
    let gL : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k := Matrix.fromBlocks H 0 0 E
    N6OnePointPlusSimple.ActBivector N6OnePointPlusSimple.rep₁ (gU * gL) =
        a • N6OnePointPlusSimple.rep₁ + (a * a - a) • N6OnePointPlusSimple.rep₂ ∧
      N6OnePointPlusSimple.ActBivector N6OnePointPlusSimple.rep₂ (gU * gL) =
        (a * a) • N6OnePointPlusSimple.rep₂ :=
  N6OnePointPlusSimple.pointwise_torus_product_lift_action
    (k := k) (A := A) (B := B) (C := C) (D := D) (E0 := E0) hA hC hD hrel hE0 (a := a) ha

/-- Right-multiplying the torus-type quotient lift by the full block-diagonal pointwise
subgroup on the direct-sum `2[a]+[b]` orbit does not change the quotient action. -/
theorem onePointPlusSimple_torus_pointwise_right_product_lift_action
    (A B C D E0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let B0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k :=
      !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k := !![a * a, 0; 0, 1]
    let gL : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k := Matrix.fromBlocks H 0 0 E
    let gU : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k :=
      Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    N6OnePointPlusSimple.ActBivector N6OnePointPlusSimple.rep₁ (gL * gU) =
        a • N6OnePointPlusSimple.rep₁ + (a * a - a) • N6OnePointPlusSimple.rep₂ ∧
      N6OnePointPlusSimple.ActBivector N6OnePointPlusSimple.rep₂ (gL * gU) =
        (a * a) • N6OnePointPlusSimple.rep₂ :=
  N6OnePointPlusSimple.torus_pointwise_right_product_lift_action
    (k := k) (A := A) (B := B) (C := C) (D := D) (E0 := E0) hA hC hD hrel hE0 (a := a) ha

/-- Left-multiplying the torus-type quotient lift by the full block-diagonal pointwise
subgroup on the direct-sum `2[a]+[b]` orbit does not change the determinant. -/
theorem onePointPlusSimple_pointwise_torus_product_lift_det
    (A B C D E0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let gU : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k :=
      Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    let B0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k :=
      !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k := !![a * a, 0; 0, 1]
    let gL : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k := Matrix.fromBlocks H 0 0 E
    Matrix.det (gU * gL) = (a * a) * (a * a) :=
  N6OnePointPlusSimple.pointwise_torus_product_lift_det
    (k := k) (A := A) (B := B) (C := C) (D := D) (E0 := E0) hA hC hD hrel hE0 (a := a) ha

/-- Right-multiplying the torus-type quotient lift by the full block-diagonal pointwise
subgroup on the direct-sum `2[a]+[b]` orbit does not change the determinant. -/
theorem onePointPlusSimple_torus_pointwise_right_product_lift_det
    (A B C D E0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let B0 : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k :=
      !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix N6OnePointPlusSimple.W N6OnePointPlusSimple.W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix N6OnePointPlusSimple.I N6OnePointPlusSimple.I k := !![a * a, 0; 0, 1]
    let gL : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k := Matrix.fromBlocks H 0 0 E
    let gU : Matrix N6OnePointPlusSimple.V N6OnePointPlusSimple.V k :=
      Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    Matrix.det (gL * gU) = (a * a) * (a * a) :=
  N6OnePointPlusSimple.torus_pointwise_right_product_lift_det
    (k := k) (A := A) (B := B) (C := C) (D := D) (E0 := E0) hA hC hD hrel hE0 (a := a) ha

end N6Summary
end Wedge2Formalization
