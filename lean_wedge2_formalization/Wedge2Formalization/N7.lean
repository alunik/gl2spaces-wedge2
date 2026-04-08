import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N7

variable {k : Type*} [Field k]

/-- The three-dimensional radical factor used for the first `n = 7` orbit. -/
abbrev I := Fin 3

/-- The `4`-dimensional quotient space carrying the split-support `n = 4` orbit. -/
abbrev W := N4.V

/-- The ambient `7`-dimensional space for the radical extension. -/
abbrev V := I ⊕ W

/-- The embedded split-support representative
`⟨e₄∧e₅ + e₆∧e₇, e₆∧e₇⟩`. -/
def rad3SplitRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4.splitRep₁ (k := k))

/-- The second basis vector of the embedded split-support representative. -/
def rad3SplitRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4.splitRep₂ (k := k))

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- Exact stabilizer of the embedded split-support pair. -/
def FixesRad3SplitPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rad3SplitRep₁ g ∧ FixesBivector rad3SplitRep₂ g

/-- Setwise stabilizer of the embedded split-support `2`-space in the chosen basis. -/
def PreservesRad3SplitSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector rad3SplitRep₁ g = α • rad3SplitRep₁ + β • rad3SplitRep₂ ∧
    ActBivector rad3SplitRep₂ g = γ • rad3SplitRep₁ + δ • rad3SplitRep₂

/-- Embedding commutes with linear combinations of the split-support basis vectors. -/
theorem embed_split_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0 (α • (N4.splitRep₁ (k := k)) + β • (N4.splitRep₂ (k := k))) =
      α • rad3SplitRep₁ + β • rad3SplitRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [rad3SplitRep₁, rad3SplitRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- If a bivector is supported on the lower-right block, then acting by a block upper
triangular matrix with zero upper-right block only changes the lower-right block. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
      Matrix.fromBlocks 0 0 0 (N4.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N4.ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

/-- Pointwise fixing of an embedded bivector reduces to pointwise fixing on the
`4`-dimensional quotient when the upper-right block is zero. -/
theorem fixes_embedded_fromBlocks_zeroUpperRight_iff
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) ↔
      N4.FixesBivector Ω D := by
  constructor
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω at h
    rw [act_embedded_fromBlocks_zeroUpperRight] at h
    have h' := congrArg Matrix.toBlocks₂₂ h
    simpa using h'
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω
    rw [act_embedded_fromBlocks_zeroUpperRight]
    simpa [N4.FixesBivector] using h

/-- The pairwise version of the lower-right reduction for the three-dimensional radical
extension. -/
theorem fixesRad3SplitPair_fromBlocks_zeroUpperRight_iff
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRad3SplitPairBivector (Matrix.fromBlocks A 0 C D) ↔
      N4.FixesSplitPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₁ (k := k)) (A := A) (C := C) (D := D)).1 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₂ (k := k)) (A := A) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₁ (k := k)) (A := A) (C := C) (D := D)).2 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₂ (k := k)) (A := A) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded pair forces
the lower-right block to fix the underlying `n = 4` split-support pair. -/
theorem fixesRad3SplitPair_lowerRight
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRad3SplitPairBivector (Matrix.fromBlocks A B C D) →
      N4.FixesSplitPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, rad3SplitRep₁, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, rad3SplitRep₂, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- Setwise preservation also reduces to the lower-right block when the upper-right block
vanishes. -/
theorem preservesRad3SplitSubspace_fromBlocks_zeroUpperRight
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hD : N4.PreservesSplitSubspaceBivector D) :
    PreservesRad3SplitSubspaceBivector (Matrix.fromBlocks A 0 C D) := by
  rcases hD with ⟨α, β, γ, δ, h₁, h₂⟩
  refine ⟨α, β, γ, δ, ?_, ?_⟩
  ·
    calc
      ActBivector rad3SplitRep₁ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4.ActBivector (N4.splitRep₁ (k := k)) D) := by
              simpa [rad3SplitRep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4.splitRep₁ (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (α • N4.splitRep₁ + β • N4.splitRep₂) := by simp [h₁]
      _ = α • rad3SplitRep₁ + β • rad3SplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := α) (β := β)
  ·
    calc
      ActBivector rad3SplitRep₂ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4.ActBivector (N4.splitRep₂ (k := k)) D) := by
              simpa [rad3SplitRep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4.splitRep₂ (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (γ • N4.splitRep₁ + δ • N4.splitRep₂) := by simp [h₂]
      _ = γ • rad3SplitRep₁ + δ • rad3SplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := γ) (β := δ)

/-- The three-dimensional radical extension of the split-support pointwise stabilizer.
The radical block is free, the radical columns are free, and the lower-right block is the
pointwise stabilizer from the `n = 4` split-support orbit. -/
theorem rad3Split_pointwise_family
    (A0 : Matrix I I k)
    (C : Matrix W I k)
    (A E : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hE : E.det = 1) :
    FixesRad3SplitPairBivector
      (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 E)) := by
  have hD : N4.FixesSplitPairBivector (Matrix.fromBlocks A 0 0 E) := by
    exact
      (N4Summary.split_pointwise_bivector_iff (k := k) (A := A) (B := 0) (C := 0) (D := E)).2
        ⟨hA, rfl, rfl, hE⟩
  exact
    (fixesRad3SplitPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := A0)
      (C := C)
      (D := Matrix.fromBlocks A 0 0 E)).2 hD

/-- For the three-dimensional radical extension of the split-support orbit, pointwise
fixing forces the upper-right block to vanish. -/
theorem fixesRad3SplitPair_fromBlocks_iff
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRad3SplitPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ N4.FixesSplitPairBivector D := by
  constructor
  · intro h
    have hD : N4.FixesSplitPairBivector D :=
      fixesRad3SplitPair_lowerRight (A := A) (B := B) (C := C) (D := D) h
    have hDblocks :
        N4.FixesSplitPairBivector
          (Matrix.fromBlocks D.toBlocks₁₁ D.toBlocks₁₂ D.toBlocks₂₁ D.toBlocks₂₂) := by
      simpa [Matrix.fromBlocks_toBlocks] using hD
    rcases
        (N4Summary.split_pointwise_bivector_iff
          (k := k)
          (A := D.toBlocks₁₁)
          (B := D.toBlocks₁₂)
          (C := D.toBlocks₂₁)
          (D := D.toBlocks₂₂)).1 hDblocks with
      ⟨hD₁₁, hD₁₂, hD₂₁, hD₂₂⟩
    have hDdet : D.det = 1 := by
      rw [← Matrix.fromBlocks_toBlocks D, hD₂₁, Matrix.det_fromBlocks_zero₂₁]
      simp [hD₁₁, hD₂₂]
    have hunitT : IsUnit ((Dᵀ).det) := by
      simpa [Matrix.det_transpose, hDdet] using (isUnit_one : IsUnit (1 : k))
    rcases h with ⟨h₁, h₂⟩
    have h₂tr : B * (N4.splitRep₂ (k := k)) * Dᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₂
      simpa [FixesBivector, rad3SplitRep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₁tr : B * (N4.splitRep₁ (k := k)) * Dᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₁
      simpa [FixesBivector, rad3SplitRep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₂zero : B * (N4.splitRep₂ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have h₁zero : B * (N4.splitRep₁ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hrow_zero (i : I) : B i = 0 := by
      have hBinl1 : B i (Sum.inl 1) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inl 0)) h₁zero
        have hneg : -B i (Sum.inl 1) = 0 := by
          simpa [N4.splitRep₁, N4.ω12, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
        exact neg_eq_zero.mp hneg
      have hBinl0 : B i (Sum.inl 0) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inl 1)) h₁zero
        simpa [N4.splitRep₁, N4.ω12, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      have hBinr1 : B i (Sum.inr 1) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inr 0)) h₂zero
        have hneg : -B i (Sum.inr 1) = 0 := by
          simpa [N4.splitRep₂, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
        exact neg_eq_zero.mp hneg
      have hBinr0 : B i (Sum.inr 0) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inr 1)) h₂zero
        simpa [N4.splitRep₂, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      ext j
      rcases j with j | j
      · fin_cases j <;> simp [hBinl0, hBinl1]
      · fin_cases j <;> simp [hBinr0, hBinr1]
    refine ⟨?_, hD⟩
    ext i j
    simpa using congrArg (fun row => row j) (hrow_zero i)
  · rintro ⟨hB, hD⟩
    simpa [hB] using
      (fixesRad3SplitPair_fromBlocks_zeroUpperRight_iff
        (k := k)
        (A := A)
        (C := C)
        (D := D)).2 hD

/-- The three-dimensional radical extension also preserves the split-support `2`-space
setwise whenever the lower-right block preserves the underlying split-support `n = 4`
subspace. -/
theorem rad3Split_setwise_family
    (A0 : Matrix I I k)
    (C : Matrix W I k)
    (A D : Matrix N4.I N4.I k) :
    PreservesRad3SplitSubspaceBivector
      (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D)) := by
  apply preservesRad3SplitSubspace_fromBlocks_zeroUpperRight
  exact N4Summary.split_diag_preserves (k := k) (A := A) (D := D)

/-- The standard torus lift on the split-support orbit extends directly to the
three-dimensional radical extension. -/
theorem rad3Split_torus_lift_action
    (A0 : Matrix I I k)
    (u v : k)
    (C : Matrix W I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks A 0 0 D
    let g : Matrix V V k := Matrix.fromBlocks A0 0 C h
    ActBivector rad3SplitRep₁ g = u • rad3SplitRep₁ + (v - u) • rad3SplitRep₂ ∧
      ActBivector rad3SplitRep₂ g = v • rad3SplitRep₂ := by
  dsimp
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  have htorus := N4Summary.split_torus_lift_action (k := k) (u := u) (v := v)
  constructor
  ·
    calc
      ActBivector rad3SplitRep₁ (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₁ (k := k)) (Matrix.fromBlocks A 0 0 D)) := by
            simpa [rad3SplitRep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₁ (k := k))
                (A := A0)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D)
      _ = Matrix.fromBlocks 0 0 0
          (u • (N4.splitRep₁ (k := k)) + (v - u) • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using htorus.1
      _ = u • rad3SplitRep₁ + (v - u) • rad3SplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := u) (β := v - u)
  ·
    calc
      ActBivector rad3SplitRep₂ (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₂ (k := k)) (Matrix.fromBlocks A 0 0 D)) := by
            simpa [rad3SplitRep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₂ (k := k))
                (A := A0)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D)
      _ = Matrix.fromBlocks 0 0 0 (v • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using htorus.2
      _ = v • rad3SplitRep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rad3SplitRep₂, Matrix.fromBlocks]

/-- The swap coset in the split normalizer also extends directly to the three-dimensional
radical extension. -/
theorem rad3Split_swapCoset_lift_action
    (A0 : Matrix I I k)
    (u v : k)
    (C : Matrix W I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k)
    let g : Matrix V V k := Matrix.fromBlocks A0 0 C h
    ActBivector rad3SplitRep₁ g = u • rad3SplitRep₁ + (v - u) • rad3SplitRep₂ ∧
      ActBivector rad3SplitRep₂ g = u • rad3SplitRep₁ + (-u) • rad3SplitRep₂ := by
  dsimp
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  have hswap := N4Summary.split_swapCoset_lift_action (k := k) (u := u) (v := v)
  constructor
  ·
    calc
      ActBivector rad3SplitRep₁
          (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₁ (k := k))
            (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) := by
            simpa [rad3SplitRep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₁ (k := k))
                (A := A0)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))
      _ = Matrix.fromBlocks 0 0 0
          (u • (N4.splitRep₁ (k := k)) + (v - u) • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using hswap.1
      _ = u • rad3SplitRep₁ + (v - u) • rad3SplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := u) (β := v - u)
  ·
    calc
      ActBivector rad3SplitRep₂
          (Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₂ (k := k))
            (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) := by
            simpa [rad3SplitRep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₂ (k := k))
                (A := A0)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))
      _ = Matrix.fromBlocks 0 0 0
          (u • (N4.splitRep₁ (k := k)) + (-u) • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using hswap.2
      _ = u • rad3SplitRep₁ + (-u) • rad3SplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := u) (β := -u)

end N7
end Wedge2Formalization
