import Wedge2Formalization.N6Summary
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

open Matrix

namespace Wedge2Formalization
namespace N7WeightedTwoPoint

variable {k : Type*} [Field k]

/-- The one-dimensional radical factor used for the weighted two-point `n = 7` orbit. -/
abbrev I := Fin 1

/-- The `6`-dimensional quotient space carrying the weighted two-point `n = 6` orbit. -/
abbrev W := N6WeightedTwoPoint.V

/-- The ambient `7`-dimensional space. -/
abbrev V := I ⊕ W

/-- The embedded weighted two-point representative in dimension `7`. -/
def radWeightedRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N6WeightedTwoPoint.rep₁ (k := k))

/-- The second basis vector of the embedded weighted two-point representative. -/
def radWeightedRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N6WeightedTwoPoint.rep₂ (k := k))

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- Exact stabilizer of the embedded weighted two-point pair. -/
def FixesRadWeightedPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector radWeightedRep₁ g ∧ FixesBivector radWeightedRep₂ g

/-- Setwise stabilizer of the embedded weighted two-point `2`-space in the chosen basis. -/
def PreservesRadWeightedSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector radWeightedRep₁ g = α • radWeightedRep₁ + β • radWeightedRep₂ ∧
    ActBivector radWeightedRep₂ g = γ • radWeightedRep₁ + δ • radWeightedRep₂

/-- The `1 × 1` scalar block in the radical direction. -/
def scalarBlock (a : k) : Matrix I I k := !![a]

/-- Embedding commutes with linear combinations of the weighted two-point basis vectors. -/
theorem embed_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0
        (α • (N6WeightedTwoPoint.rep₁ (k := k)) + β • (N6WeightedTwoPoint.rep₂ (k := k))) =
      α • radWeightedRep₁ + β • radWeightedRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [radWeightedRep₁, radWeightedRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The embedded weighted two-point line has the same lower-right determinant as the
underlying `n = 6` orbit. -/
theorem lowerRight_det_in_basis
    (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • radWeightedRep₁ + b • radWeightedRep₂)) =
      (a * a * (a + b)) * (a * a * (a + b)) := by
  rw [← embed_linearCombination (k := k) (α := a) (β := b)]
  simp only [Matrix.toBlocks_fromBlocks₂₂]
  exact N6WeightedTwoPoint.det_in_basis (k := k) (a := a) (b := b)

/-- On the embedded weighted two-point line, the lower-right determinant vanishes
exactly on the two distinguished projective points. -/
theorem lowerRight_det_zero_iff
    (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • radWeightedRep₁ + b • radWeightedRep₂)) = (0 : k) ↔
      a = 0 ∨ a + b = 0 := by
  change Matrix.det (Matrix.toBlocks₂₂ (a • radWeightedRep₁ + b • radWeightedRep₂)) = (0 : k) ↔
      a = 0 ∨ a + b = 0
  rw [lowerRight_det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    have hsq : a * a * (a + b) = 0 := by
      rcases mul_eq_zero.mp h with h0 | h0 <;> exact h0
    rcases mul_eq_zero.mp hsq with haa | hab
    · rcases mul_eq_zero.mp haa with ha | ha
      · exact Or.inl ha
      · exact Or.inl ha
    · exact Or.inr hab
  · rintro (ha | hab)
    · simp [ha]
    · simp [hab]

/-- If a bivector is supported on the lower-right block, then acting by a block upper
triangular matrix with zero upper-right block only changes the lower-right block. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
      Matrix.fromBlocks 0 0 0 (N6WeightedTwoPoint.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N6WeightedTwoPoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- Pointwise fixing of an embedded bivector reduces to pointwise fixing on the
`6`-dimensional quotient when the upper-right block is zero. -/
theorem fixes_embedded_fromBlocks_zeroUpperRight_iff
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) ↔
      N6WeightedTwoPoint.FixesBivector Ω D := by
  constructor
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω at h
    rw [act_embedded_fromBlocks_zeroUpperRight] at h
    have h' := congrArg Matrix.toBlocks₂₂ h
    simpa [N6WeightedTwoPoint.FixesBivector, N6WeightedTwoPoint.ActBivector] using h'
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω
    rw [act_embedded_fromBlocks_zeroUpperRight]
    simpa [N6WeightedTwoPoint.FixesBivector, N6WeightedTwoPoint.ActBivector] using h

/-- The pairwise version of the lower-right reduction for the weighted two-point radical
extension. -/
theorem fixesRadWeightedPair_fromBlocks_zeroUpperRight_iff
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadWeightedPairBivector (Matrix.fromBlocks a 0 C D) ↔
      N6WeightedTwoPoint.FixesPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6WeightedTwoPoint.rep₁ (k := k)) (a := a) (C := C) (D := D)).1 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6WeightedTwoPoint.rep₂ (k := k)) (a := a) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6WeightedTwoPoint.rep₁ (k := k)) (a := a) (C := C) (D := D)).2 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N6WeightedTwoPoint.rep₂ (k := k)) (a := a) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded pair forces
the lower-right block to fix the underlying `n = 6` weighted two-point pair. -/
theorem fixesRadWeightedPair_lowerRight
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadWeightedPairBivector (Matrix.fromBlocks a B C D) →
      N6WeightedTwoPoint.FixesPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, radWeightedRep₁, N6WeightedTwoPoint.FixesBivector,
      N6WeightedTwoPoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, radWeightedRep₂, N6WeightedTwoPoint.FixesBivector,
      N6WeightedTwoPoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij

/-- For the weighted two-point radical extension, pointwise fixing forces the
upper-right block to vanish. -/
theorem fixesRadWeightedPair_fromBlocks_iff
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadWeightedPairBivector (Matrix.fromBlocks a B C D) ↔
      B = 0 ∧ N6WeightedTwoPoint.FixesPairBivector D := by
  constructor
  · intro h
    have hD : N6WeightedTwoPoint.FixesPairBivector D :=
      fixesRadWeightedPair_lowerRight (a := a) (B := B) (C := C) (D := D) h
    have hDblocks :
        N6WeightedTwoPoint.FixesPairBivector
          (Matrix.fromBlocks D.toBlocks₁₁ D.toBlocks₁₂ D.toBlocks₂₁ D.toBlocks₂₂) := by
      simpa [Matrix.fromBlocks_toBlocks] using hD
    rcases
      (N6Summary.weightedTwoPoint_pointwise_bivector_iff
        (k := k)
        (A := D.toBlocks₁₁)
        (B := D.toBlocks₁₂)
        (C := D.toBlocks₂₁)
        (D := D.toBlocks₂₂)).1 hDblocks with
      ⟨hD₁₂, hD₂₁, hD₁₁, hD₂₂⟩
    have hAunit : IsUnit (D.toBlocks₁₁.det) := by
      have hdetEq := congrArg Matrix.det hD₁₁
      simp [N4.FixesBivector, N4.ActBivector, Matrix.det_mul, Matrix.det_transpose,
        N6WeightedTwoPoint.det_splitRep₁] at hdetEq
      exact ⟨⟨D.toBlocks₁₁.det, D.toBlocks₁₁.det, hdetEq, hdetEq⟩, rfl⟩
    have hDunit : IsUnit D.det := by
      rw [← Matrix.fromBlocks_toBlocks D, hD₂₁, Matrix.det_fromBlocks_zero₂₁]
      simpa [hD₂₂] using hAunit.mul (isUnit_one : IsUnit (1 : k))
    have hunitT : IsUnit ((Dᵀ).det) := by
      simpa [Matrix.det_transpose] using hDunit
    rcases h with ⟨h₁, h₂⟩
    have h₂tr : B * (N6WeightedTwoPoint.rep₂ (k := k)) * Dᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
      simpa [FixesBivector, radWeightedRep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₁tr : B * (N6WeightedTwoPoint.rep₁ (k := k)) * Dᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
      simpa [FixesBivector, radWeightedRep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₂zero : B * (N6WeightedTwoPoint.rep₂ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have h₁zero : B * (N6WeightedTwoPoint.rep₁ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hBright1 : B 0 (Sum.inr 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 0)) h₂zero
      have hneg : -B 0 (Sum.inr 1) = 0 := by
        simpa [N6WeightedTwoPoint.rep₂, N4.J, Matrix.mul_apply, Fintype.sum_sum_type] using hij
      exact neg_eq_zero.mp hneg
    have hBright0 : B 0 (Sum.inr 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 1)) h₂zero
      simpa [N6WeightedTwoPoint.rep₂, N4.J, Matrix.mul_apply, Fintype.sum_sum_type] using hij
    have hBleft11 : B 0 (Sum.inl (Sum.inl 1)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 0))) h₁zero
      have hneg : -B 0 (Sum.inl (Sum.inl 1)) = 0 := by
        simpa [N6WeightedTwoPoint.rep₁, N4.splitRep₁, N4.ω12, N4.ω34, N4.J,
          Matrix.mul_apply, Fintype.sum_sum_type, hBright0, hBright1] using hij
      exact neg_eq_zero.mp hneg
    have hBleft10 : B 0 (Sum.inl (Sum.inl 0)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 1))) h₁zero
      simpa [N6WeightedTwoPoint.rep₁, N4.splitRep₁, N4.ω12, N4.ω34, N4.J,
        Matrix.mul_apply, Fintype.sum_sum_type, hBright0, hBright1] using hij
    have hBleft01 : B 0 (Sum.inl (Sum.inr 1)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inr 0))) h₁zero
      have hneg : -B 0 (Sum.inl (Sum.inr 1)) = 0 := by
        simpa [N6WeightedTwoPoint.rep₁, N4.splitRep₁, N4.ω12, N4.ω34, N4.J,
          Matrix.mul_apply, Fintype.sum_sum_type, hBright0, hBright1] using hij
      exact neg_eq_zero.mp hneg
    have hBleft00 : B 0 (Sum.inl (Sum.inr 0)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inr 1))) h₁zero
      simpa [N6WeightedTwoPoint.rep₁, N4.splitRep₁, N4.ω12, N4.ω34, N4.J,
        Matrix.mul_apply, Fintype.sum_sum_type, hBright0, hBright1] using hij
    refine ⟨?_, hD⟩
    ext i j
    fin_cases i
    rcases j with j | j
    · rcases j with j | j
      · fin_cases j <;> simp [hBleft10, hBleft11]
      · fin_cases j <;> simp [hBleft00, hBleft01]
    · fin_cases j <;> simp [hBright0, hBright1]
  · rintro ⟨hB, hD⟩
    simpa [hB] using
      (fixesRadWeightedPair_fromBlocks_zeroUpperRight_iff
        (k := k) (a := a) (C := C) (D := D)).2 hD

/-- A one-dimensional radical-extension family inside the pointwise stabilizer of the
embedded weighted two-point orbit. -/
theorem radWeighted_pointwise_family
    (a : k)
    (C : Matrix W I k)
    (A : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k)
    (D : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k)
    (hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1) :
    FixesRadWeightedPairBivector
      (Matrix.fromBlocks (scalarBlock (k := k) a) 0 C (Matrix.fromBlocks A 0 0 D)) := by
  have hpair : N6WeightedTwoPoint.FixesPairBivector (Matrix.fromBlocks A 0 0 D) := by
    exact (N6WeightedTwoPoint.fixesPair_fromBlocks_iff
      (k := k) (A := A) (B := 0) (C := 0) (D := D)).2 ⟨rfl, rfl, hA, hD⟩
  exact
    (fixesRadWeightedPair_fromBlocks_zeroUpperRight_iff
      (k := k) (a := scalarBlock (k := k) a) (C := C) (D := Matrix.fromBlocks A 0 0 D)).2 hpair

/-- The diagonal torus lift on the weighted two-point orbit extends to the one-dimensional
radical extension. -/
theorem radWeighted_torus_lift_action
    (a0 u v : k)
    (C : Matrix W I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix N6WeightedTwoPoint.W N6WeightedTwoPoint.W k := Matrix.fromBlocks A 0 0 A
    let E : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks H 0 0 E
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C h
    ActBivector radWeightedRep₁ g =
        u • radWeightedRep₁ + (v - u) • radWeightedRep₂ ∧
      ActBivector radWeightedRep₂ g =
        v • radWeightedRep₂ := by
  dsimp
  have htorus := N6WeightedTwoPoint.torus_lift_action (k := k) (u := u) (v := v)
  let h : Matrix W W k :=
    Matrix.fromBlocks
      (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
        (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
      0 0
      (!![v, 0; 0, 1] : Matrix N6WeightedTwoPoint.I N6WeightedTwoPoint.I k)
  have hblock1 :
      ActBivector radWeightedRep₁ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C h) =
        Matrix.fromBlocks 0 0 0
          (N6WeightedTwoPoint.ActBivector (N6WeightedTwoPoint.rep₁ (k := k)) h) := by
    simpa [radWeightedRep₁, h] using
      act_embedded_fromBlocks_zeroUpperRight
        (Ω := N6WeightedTwoPoint.rep₁ (k := k))
        (a := scalarBlock (k := k) a0) (C := C) (D := h)
  have hblock2 :
      ActBivector radWeightedRep₂ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C h) =
        Matrix.fromBlocks 0 0 0
          (N6WeightedTwoPoint.ActBivector (N6WeightedTwoPoint.rep₂ (k := k)) h) := by
    simpa [radWeightedRep₂, h] using
      act_embedded_fromBlocks_zeroUpperRight
        (Ω := N6WeightedTwoPoint.rep₂ (k := k))
        (a := scalarBlock (k := k) a0) (C := C) (D := h)
  constructor
  ·
    calc
      ActBivector radWeightedRep₁ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C h) =
        Matrix.fromBlocks 0 0 0
          (N6WeightedTwoPoint.ActBivector (N6WeightedTwoPoint.rep₁ (k := k)) h) := by
            simpa using hblock1
      _ = Matrix.fromBlocks 0 0 0
          (u • N6WeightedTwoPoint.rep₁ (k := k) + (v - u) • N6WeightedTwoPoint.rep₂ (k := k)) := by
            simp [h, htorus.1]
      _ = u • radWeightedRep₁ + (v - u) • radWeightedRep₂ := by
            simpa using (embed_linearCombination (k := k) (α := u) (β := v - u))
  ·
    calc
      ActBivector radWeightedRep₂ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C h) =
        Matrix.fromBlocks 0 0 0
          (N6WeightedTwoPoint.ActBivector (N6WeightedTwoPoint.rep₂ (k := k)) h) := by
            simpa using hblock2
      _ = Matrix.fromBlocks 0 0 0 (v • N6WeightedTwoPoint.rep₂ (k := k)) := by
            simp [h, htorus.2]
      _ = v • radWeightedRep₂ := by
            simpa using (embed_linearCombination (k := k) (α := 0) (β := v))

end N7WeightedTwoPoint
end Wedge2Formalization
