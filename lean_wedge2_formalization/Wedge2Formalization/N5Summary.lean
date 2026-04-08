import Wedge2Formalization.N5
import Wedge2Formalization.N3PureSingular
import Wedge2Formalization.N5PureSingularLong
import Wedge2Formalization.N5SimplePoint
import Wedge2Formalization.N5OnePoint
import Wedge2Formalization.N5PureSingular

namespace Wedge2Formalization
namespace N5Summary

open Matrix

variable {k : Type*} [Field k]

/-- Summary form of the first `n = 5` pointwise stabilizer calculation. This is the
radical extension of the split-support `n = 4` orbit. -/
theorem radSplit_pointwise_bivector_iff
    (a : Matrix N5.I N5.I k)
    (R : Matrix N5.I N5.W k)
    (C : Matrix N5.W N5.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N5.FixesRadSplitPairBivector
      (Matrix.fromBlocks a R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      R = 0 ∧ A.det = 1 ∧ B₁ = 0 ∧ C₁ = 0 ∧ D.det = 1 := by
  rw [N5.fixesRadSplitPair_fromBlocks_iff]
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

/-- The radical extension family preserves the split-support `2`-space setwise. -/
theorem radSplit_diag_preserves
    (a : k)
    (C : Matrix N5.W N5.I k)
    (A D : Matrix N4.I N4.I k) :
    N5.PreservesRadSplitSubspaceBivector
      (Matrix.fromBlocks (N5.scalarBlock (k := k) a) 0 C (Matrix.fromBlocks A 0 0 D)) :=
  N5.radSplit_setwise_family (k := k) (a := a) (C := C) (A := A) (D := D)

/-- The standard torus lift on the quotient extends to the radical extension. -/
theorem radSplit_torus_lift_action
    (a u v : k)
    (C : Matrix N5.W N5.I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix N5.W N5.W k := Matrix.fromBlocks A 0 0 D
    let g : Matrix N5.V N5.V k := Matrix.fromBlocks (N5.scalarBlock (k := k) a) 0 C h
    N5.ActBivector N5.radSplitRep₁ g = u • N5.radSplitRep₁ + (v - u) • N5.radSplitRep₂ ∧
      N5.ActBivector N5.radSplitRep₂ g = v • N5.radSplitRep₂ :=
  N5.radSplit_torus_lift_action (k := k) (a := a) (u := u) (v := v) (C := C)

/-- The swap coset in the split normalizer also extends to the radical extension. -/
theorem radSplit_swapCoset_lift_action
    (a u v : k)
    (C : Matrix N5.W N5.I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix N5.W N5.W k := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k)
    let g : Matrix N5.V N5.V k := Matrix.fromBlocks (N5.scalarBlock (k := k) a) 0 C h
    N5.ActBivector N5.radSplitRep₁ g = u • N5.radSplitRep₁ + (v - u) • N5.radSplitRep₂ ∧
      N5.ActBivector N5.radSplitRep₂ g = u • N5.radSplitRep₁ + (-u) • N5.radSplitRep₂ :=
  N5.radSplit_swapCoset_lift_action (k := k) (a := a) (u := u) (v := v) (C := C)

/-- Summary form of the second `n = 5` pointwise stabilizer calculation. This is the
radical extension of the repeated-support `n = 4` orbit. -/
theorem radOnePoint_pointwise_bivector_iff
    (a : Matrix N5.I N5.I k)
    (R : Matrix N5.I N5.W k)
    (C : Matrix N5.W N5.I k)
    (A B₁ C₁ D : Matrix N4.I N4.I k) :
    N5OnePoint.FixesRadOnePointPairBivector
      (Matrix.fromBlocks a R C (Matrix.fromBlocks A B₁ C₁ D)) ↔
      R = 0 ∧ A.det = 1 ∧ C₁ = 0 ∧ D = A ∧ A * N4.J * B₁ᵀ + B₁ * N4.J * Aᵀ = 0 := by
  rw [N5OnePoint.fixesRadOnePointPair_fromBlocks_iff]
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
to the radical extension. -/
theorem radOnePoint_borel_lift_action
    (a b t : k)
    (ha : a ≠ 0)
    (C : Matrix N5.W N5.I k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let g : Matrix N5.V N5.V k :=
      Matrix.fromBlocks
        (N5.scalarBlock (k := k) t)
        0
        C
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
    N5.ActBivector N5OnePoint.radOnePointRep₁ g =
        a • N5OnePoint.radOnePointRep₁ + b • N5OnePoint.radOnePointRep₂ ∧
      N5.ActBivector N5OnePoint.radOnePointRep₂ g =
        (a * a) • N5OnePoint.radOnePointRep₂ :=
  N5OnePoint.radOnePoint_borel_lift_action (k := k) (a := a) (b := b) (t := t) ha C

/-- Summary form of the third `n = 5` pointwise stabilizer calculation. This is the
radical extension of the pure singular `n = 4` orbit. -/
theorem radPureSingular_pointwise_bivector_iff
    (a : Matrix N5PureSingular.I N5PureSingular.I k)
    (R : Matrix N5PureSingular.I N5PureSingular.W k)
    (C : Matrix N5PureSingular.W N5PureSingular.I k)
    (D : Matrix N4PureSingular.I N4PureSingular.I k) :
    N5PureSingular.FixesRadPureSingularPairBivector
      (Matrix.fromBlocks a R C D) ↔
      R = N5PureSingular.upperRightLast (k := k) (R 0 3) ∧
      D 0 0 ≠ 0 ∧
      D =
        N4PureSingular.pureSingularShape
          (k := k) (D 0 0) (D 0 1) (D 0 2) (D 0 3) (D 1 3) (D 2 3) (D 3 3) := by
  exact
    N5PureSingular.fixesRadPureSingularPair_fromBlocks_iff
      (k := k) (a := a) (B := R) (C := C) (D := D)

/-- The full `GL₂` quotient lift on the pure singular orbit extends to the radical
extension. -/
theorem radPureSingular_GL2_lift_action
    (α β γ δ u : k)
    (C : Matrix N5PureSingular.W N5PureSingular.I k) :
    let h :=
      N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1
    let g : Matrix N5PureSingular.V N5PureSingular.V k :=
      Matrix.fromBlocks (N5PureSingular.scalarBlock (k := k) u) 0 C h
    N5PureSingular.ActBivector N5PureSingular.radPureSingularRep₁ g =
        α • N5PureSingular.radPureSingularRep₁ + β • N5PureSingular.radPureSingularRep₂ ∧
      N5PureSingular.ActBivector N5PureSingular.radPureSingularRep₂ g =
        γ • N5PureSingular.radPureSingularRep₁ + δ • N5PureSingular.radPureSingularRep₂ :=
  N5PureSingular.radPureSingular_GL2_lift_action
    (k := k) (α := α) (β := β) (γ := γ) (δ := δ) (u := u) C

/-- The pure singular length-two block in dimension `5` has the expected pointwise
scaling family. -/
theorem pureSingularLong_pointwise_scale_family
    (t : k)
    (ht : t ≠ 0) :
    N5PureSingularLong.FixesPairBivector
      (Matrix.fromBlocks
        ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
        0
        0
        (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) :=
  N5PureSingularLong.pointwise_scale_family (k := k) t ht

/-- The pure singular length-two block in dimension `5` also admits the expected
`4`-parameter pointwise unipotent family. -/
theorem pureSingularLong_lowerHankel_pointwise_family
    (u v w z : k) :
    N5PureSingularLong.FixesPairBivector
      (Matrix.fromBlocks
        (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
        0
        (N5PureSingularLong.lowerHankel (k := k) u v w z)
        (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)) :=
  N5PureSingularLong.lowerHankel_pointwise_family (k := k) u v w z

/-- Combining the scalar torus and the lower-left Hankel family gives a concrete
`5`-parameter pointwise subgroup on the pure singular length-two orbit. -/
theorem pureSingularLong_scaleLowerHankel_product_pointwise
    (t u v w z : k)
    (ht : t ≠ 0) :
    let gScale : Matrix N5PureSingularLong.V N5PureSingularLong.V k :=
      Matrix.fromBlocks
        ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
        0
        0
        (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))
    let gH : Matrix N5PureSingularLong.V N5PureSingularLong.V k :=
      Matrix.fromBlocks
        (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
        0
        (N5PureSingularLong.lowerHankel (k := k) u v w z)
        (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)
    N5PureSingularLong.FixesPairBivector (gH * gScale) :=
  N5PureSingularLong.scaleLowerHankel_product_pointwise
    (k := k) (t := t) (u := u) (v := v) (w := w) (z := z) ht

/-- The pure singular length-two block in dimension `5` admits the full `GL₂`
quotient lift. -/
theorem pureSingularLong_GL2_scaled_lift_action
    (a b c d : k) :
    let P : Matrix N5PureSingularLong.W N5PureSingularLong.W k :=
      N5PureSingularLong.sym2Raw (k := k) a b c d
    let Q : Matrix N5PureSingularLong.I N5PureSingularLong.I k :=
      N5PureSingularLong.adjointT (k := k) a b c d
    let g : Matrix N5PureSingularLong.V N5PureSingularLong.V k := Matrix.fromBlocks P 0 0 Q
    N5PureSingularLong.ActBivector N5PureSingularLong.rep₁ g =
        (N5PureSingularLong.Delta a b c d) •
          (a • N5PureSingularLong.rep₁ + c • N5PureSingularLong.rep₂) ∧
      N5PureSingularLong.ActBivector N5PureSingularLong.rep₂ g =
        (N5PureSingularLong.Delta a b c d) •
          (b • N5PureSingularLong.rep₁ + d • N5PureSingularLong.rep₂) :=
  N5PureSingularLong.GL2_scaled_lift_action (k := k) (a := a) (b := b) (c := c) (d := d)

/-- The direct-sum `[a]` family in dimension `5` admits the expected pointwise Levi
subgroup. -/
theorem simplePoint_pointwise_levi_family
    (a : k)
    (ha : a ≠ 0)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    N5SimplePoint.FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        E) :=
  N5SimplePoint.pointwise_levi_family (k := k) a ha E hE

/-- The direct-sum `[a]` family in dimension `5` also has the expected `U_4`
pointwise radical. -/
theorem simplePoint_pointwise_unipotent_family
    (x y p q : k) :
    N5SimplePoint.FixesPairBivector
      (N5SimplePoint.pointwiseUnipotent (k := k) x y p q) :=
  N5SimplePoint.pointwise_unipotent_family (k := k) x y p q

/-- Combining the pointwise unipotent family with the Levi family gives a concrete
pointwise subgroup of the direct-sum `[a]` stabilizer in dimension `5`. -/
theorem simplePoint_pointwise_product_family
    (x y p q a : k)
    (ha : a ≠ 0)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    let gU : Matrix N5SimplePoint.V N5SimplePoint.V k :=
      N5SimplePoint.pointwiseUnipotent (k := k) x y p q
    let gL : Matrix N5SimplePoint.V N5SimplePoint.V k :=
      Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        E
    N5SimplePoint.FixesPairBivector (gU * gL) :=
  N5SimplePoint.pointwise_product_family
    (k := k)
    (x := x)
    (y := y)
    (p := p)
    (q := q)
    (a := a)
    ha
    (E := E)
    hE

/-- The direct-sum `[a]` family in dimension `5` also admits an explicit general
pointwise shape combining the `U_4` radical with the `G_m x SL_2` Levi factor. -/
theorem simplePoint_pointwise_shape_family
    (a x y u v : k)
    (ha : a ≠ 0)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    let B : Matrix N5SimplePoint.W N5SimplePoint.I k := !![(0 : k), 0; 0, 0; u, v]
    let C : Matrix N5SimplePoint.I N5SimplePoint.W k :=
      !![a * (u * E 0 1 - v * E 0 0), 0, 0;
         a * (u * E 1 1 - v * E 1 0), 0, 0]
    N5SimplePoint.FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        B
        C
        E) :=
  N5SimplePoint.pointwise_shape_family
    (k := k) (a := a) (x := x) (y := y) (u := u) (v := v) ha (E := E) hE

/-- If the upper-left block is already in the pure singular shape, then pointwise
fixing forces the remaining blocks on the direct-sum `[a]` orbit to have the expected
form. -/
theorem simplePoint_pureShape_iff
    (a x y : k)
    (ha : a ≠ 0)
    (B : Matrix N5SimplePoint.W N5SimplePoint.I k)
    (C : Matrix N5SimplePoint.I N5SimplePoint.W k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k) :
    N5SimplePoint.FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        B
        C
        E) ↔
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
  N5SimplePoint.fixesPair_fromBlocks_pureShape_iff
    (k := k) (a := a) (x := x) (y := y) ha (B := B) (C := C) (D := E)

/-- A concrete Borel lift on the direct-sum `[a]` family in dimension `5`. -/
theorem simplePoint_borel_lift_action
    (a b : k)
    (ha : a ≠ 0) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let g : Matrix N5SimplePoint.V N5SimplePoint.V k := Matrix.fromBlocks G 0 0 E
    N5SimplePoint.ActBivector N5SimplePoint.rep₁ g =
        a • N5SimplePoint.rep₁ + b • N5SimplePoint.rep₂ ∧
      N5SimplePoint.ActBivector N5SimplePoint.rep₂ g =
        (a * a) • N5SimplePoint.rep₂ :=
  N5SimplePoint.borel_lift_action (k := k) (a := a) (b := b) ha

/-- The determinant of the standard Borel lift on the direct-sum `[a]` family in
dimension `5` is `a^4`. -/
theorem simplePoint_borel_lift_det
    (a b : k) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let g : Matrix N5SimplePoint.V N5SimplePoint.V k := Matrix.fromBlocks G 0 0 E
    Matrix.det g = (a * a) * (a * a) :=
  N5SimplePoint.borel_lift_det (k := k) (a := a) (b := b)

end N5Summary
end Wedge2Formalization
