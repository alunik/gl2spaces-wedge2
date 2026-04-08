import Wedge2Formalization.N6Summary
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

open Matrix

namespace Wedge2Formalization
namespace N7ThreePoint

variable {k : Type*} [Field k]

/-- The one-dimensional radical factor used for the three-point `n = 7` orbit. -/
abbrev I := Fin 1

/-- The `6`-dimensional quotient space carrying the three-point `n = 6` orbit. -/
abbrev W := N6ThreePoint.V

/-- The ambient `7`-dimensional space. -/
abbrev V := I ⊕ W

/-- The embedded three-point representative in dimension `7`. -/
def radThreePointRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N6ThreePoint.rep₁ (k := k))

/-- The second basis vector of the embedded three-point representative. -/
def radThreePointRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N6ThreePoint.rep₂ (k := k))

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- Exact stabilizer of the embedded three-point pair. -/
def FixesRadThreePointPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector radThreePointRep₁ g ∧ FixesBivector radThreePointRep₂ g

/-- Embedding commutes with linear combinations of the three-point basis vectors. -/
theorem embed_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0
        (α • (N6ThreePoint.rep₁ (k := k)) + β • (N6ThreePoint.rep₂ (k := k))) =
      α • radThreePointRep₁ + β • radThreePointRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [radThreePointRep₁, radThreePointRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The embedded three-point line has the same lower-right determinant as the
underlying `n = 6` orbit. -/
theorem lowerRight_det_in_basis
    (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • radThreePointRep₁ + b • radThreePointRep₂)) =
      (a * (a + b) * b) * (a * (a + b) * b) := by
  rw [← embed_linearCombination (k := k) (α := a) (β := b)]
  simp only [Matrix.toBlocks_fromBlocks₂₂]
  exact N6ThreePoint.det_in_basis (k := k) (a := a) (b := b)

/-- On the embedded three-point line, the lower-right determinant vanishes exactly on
the three distinguished projective points. -/
theorem lowerRight_det_zero_iff
    (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • radThreePointRep₁ + b • radThreePointRep₂)) = (0 : k) ↔
      a = 0 ∨ a + b = 0 ∨ b = 0 := by
  rw [lowerRight_det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    have hprod : a * (a + b) * b = 0 := by
      rcases mul_eq_zero.mp h with h0 | h0 <;> exact h0
    rcases mul_eq_zero.mp hprod with hab | hb
    · rcases mul_eq_zero.mp hab with ha | hab'
      · exact Or.inl ha
      · exact Or.inr (Or.inl hab')
    · exact Or.inr (Or.inr hb)
  · rintro (ha | hab | hb)
    · simp [ha]
    · simp [hab]
    · simp [hb]

/-- If a bivector is supported on the lower-right block, then acting by a block upper
triangular matrix with zero upper-right block only changes the lower-right block. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
      Matrix.fromBlocks 0 0 0 (N6ThreePoint.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N6ThreePoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- Pointwise fixing of an embedded bivector reduces to pointwise fixing on the
`6`-dimensional quotient when the upper-right block is zero. -/
theorem fixes_embedded_fromBlocks_zeroUpperRight_iff
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) ↔
      N6ThreePoint.FixesBivector Ω D := by
  constructor
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω at h
    rw [act_embedded_fromBlocks_zeroUpperRight] at h
    have h' := congrArg Matrix.toBlocks₂₂ h
    simpa [N6ThreePoint.FixesBivector, N6ThreePoint.ActBivector] using h'
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω
    rw [act_embedded_fromBlocks_zeroUpperRight]
    simpa [N6ThreePoint.FixesBivector, N6ThreePoint.ActBivector] using h

/-- The pairwise version of the lower-right reduction for the three-point radical
extension. -/
theorem fixesRadThreePointPair_fromBlocks_zeroUpperRight_iff
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadThreePointPairBivector (Matrix.fromBlocks a 0 C D) ↔
      N6ThreePoint.FixesPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6ThreePoint.rep₁ (k := k)) (a := a) (C := C) (D := D)).1 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6ThreePoint.rep₂ (k := k)) (a := a) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6ThreePoint.rep₁ (k := k)) (a := a) (C := C) (D := D)).2 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6ThreePoint.rep₂ (k := k)) (a := a) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded pair forces
the lower-right block to fix the underlying `n = 6` three-point pair. -/
theorem fixesRadThreePointPair_lowerRight
    (a : Matrix I I k)
    (R : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadThreePointPairBivector (Matrix.fromBlocks a R C D) →
      N6ThreePoint.FixesPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, radThreePointRep₁, N6ThreePoint.FixesBivector,
      N6ThreePoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, radThreePointRep₂, N6ThreePoint.FixesBivector,
      N6ThreePoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij

/-- For the three-point radical extension, pointwise fixing forces the upper-right block
to vanish. -/
theorem fixesRadThreePointPair_fromBlocks_iff
    (a : Matrix I I k)
    (R : Matrix I W k)
    (C0 : Matrix W I k)
    (A : Matrix N6ThreePoint.W N6ThreePoint.W k)
    (B : Matrix N6ThreePoint.W N6ThreePoint.I k)
    (C : Matrix N6ThreePoint.I N6ThreePoint.W k)
    (D : Matrix N6ThreePoint.I N6ThreePoint.I k) :
    FixesRadThreePointPairBivector (Matrix.fromBlocks a R C0 (Matrix.fromBlocks A B C D)) ↔
      R = 0 ∧ B = 0 ∧ C = 0 ∧ N4.FixesSplitPairBivector A ∧ D.det = 1 := by
  constructor
  · intro h
    have hD :
        N6ThreePoint.FixesPairBivector (Matrix.fromBlocks A B C D) :=
      fixesRadThreePointPair_lowerRight
        (a := a) (R := R) (C := C0) (D := Matrix.fromBlocks A B C D) h
    rcases
      (N6Summary.threePoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B)
        (C := C)
        (D := D)).1 hD with
      ⟨hB, hC, hA, hDdet⟩
    have hAunit : IsUnit A.det := by
      rcases
        (N4Summary.split_pointwise_bivector_iff
          (k := k)
          (A := A.toBlocks₁₁)
          (B := A.toBlocks₁₂)
          (C := A.toBlocks₂₁)
          (D := A.toBlocks₂₂)).1 (by simpa [Matrix.fromBlocks_toBlocks] using hA) with
        ⟨hA11, hA12, hA21, hA22⟩
      rw [← Matrix.fromBlocks_toBlocks A, hA21, Matrix.det_fromBlocks_zero₂₁]
      simp [hA11, hA12, hA22]
    have hHunit : IsUnit ((Matrix.fromBlocks A B C D).det) := by
      rw [hC, Matrix.det_fromBlocks_zero₂₁]
      simpa [hB, hDdet] using hAunit.mul (isUnit_one : IsUnit (1 : k))
    have hunitT : IsUnit (((Matrix.fromBlocks A B C D)ᵀ).det) := by
      simpa [Matrix.det_transpose] using hHunit
    rcases h with ⟨h₁, h₂⟩
    have h₂tr : R * (N6ThreePoint.rep₂ (k := k)) * (Matrix.fromBlocks A B C D)ᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
      simpa [FixesBivector, radThreePointRep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₁tr : R * (N6ThreePoint.rep₁ (k := k)) * (Matrix.fromBlocks A B C D)ᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
      simpa [FixesBivector, radThreePointRep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₂zero : R * (N6ThreePoint.rep₂ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * ((Matrix.fromBlocks A B C D)ᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have h₁zero : R * (N6ThreePoint.rep₁ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * ((Matrix.fromBlocks A B C D)ᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hRright1 : R 0 (Sum.inr 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 0)) h₂zero
      have hneg : -R 0 (Sum.inr 1) = 0 := by
        simpa [N6ThreePoint.rep₂, N4.J, Matrix.mul_apply, Fintype.sum_sum_type] using hij
      exact neg_eq_zero.mp hneg
    have hRright0 : R 0 (Sum.inr 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 1)) h₂zero
      simpa [N6ThreePoint.rep₂, N4.J, Matrix.mul_apply, Fintype.sum_sum_type] using hij
    have hRmid1 : R 0 (Sum.inl (Sum.inr 1)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inr 0))) h₂zero
      have hneg : -R 0 (Sum.inl (Sum.inr 1)) = 0 := by
        simpa [N6ThreePoint.rep₂, N4.splitRep₂, N4.ω34, N4.J,
          Matrix.mul_apply, Fintype.sum_sum_type, hRright0, hRright1] using hij
      exact neg_eq_zero.mp hneg
    have hRmid0 : R 0 (Sum.inl (Sum.inr 0)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inr 1))) h₂zero
      simpa [N6ThreePoint.rep₂, N4.splitRep₂, N4.ω34, N4.J,
        Matrix.mul_apply, Fintype.sum_sum_type, hRright0, hRright1] using hij
    have hRleft1 : R 0 (Sum.inl (Sum.inl 1)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 0))) h₁zero
      have hneg : -R 0 (Sum.inl (Sum.inl 1)) = 0 := by
        simpa [N6ThreePoint.rep₁, N4.splitRep₁, N4.ω12, N4.ω34, N4.J,
          Matrix.mul_apply, Fintype.sum_sum_type, hRmid0, hRmid1, hRright0, hRright1] using hij
      exact neg_eq_zero.mp hneg
    have hRleft0 : R 0 (Sum.inl (Sum.inl 0)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 1))) h₁zero
      simpa [N6ThreePoint.rep₁, N4.splitRep₁, N4.ω12, N4.ω34, N4.J,
        Matrix.mul_apply, Fintype.sum_sum_type, hRmid0, hRmid1, hRright0, hRright1] using hij
    refine ⟨?_, hB, hC, hA, hDdet⟩
    ext i j
    fin_cases i
    rcases j with j | j
    · rcases j with j | j
      · fin_cases j <;> simp [hRleft0, hRleft1]
      · fin_cases j <;> simp [hRmid0, hRmid1]
    · fin_cases j <;> simp [hRright0, hRright1]
  · rintro ⟨hR, hB, hC, hA, hDdet⟩
    have hpair :
        N6ThreePoint.FixesPairBivector (Matrix.fromBlocks A B C D) := by
      exact
        (N6Summary.threePoint_pointwise_bivector_iff
          (k := k)
          (A := A)
          (B := B)
          (C := C)
          (D := D)).2 ⟨hB, hC, hA, hDdet⟩
    simpa [hR] using
      (fixesRadThreePointPair_fromBlocks_zeroUpperRight_iff
        (k := k) (a := a) (C := C0) (D := Matrix.fromBlocks A B C D)).2 hpair

end N7ThreePoint
end Wedge2Formalization
