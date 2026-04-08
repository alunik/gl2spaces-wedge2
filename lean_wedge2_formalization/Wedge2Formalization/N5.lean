import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N5

variable {k : Type*} [Field k]

/-- The one-dimensional radical factor used for the first `n = 5` orbit. -/
abbrev I := Fin 1

/-- The `4`-dimensional quotient space carrying the split-support `n = 4` orbit. -/
abbrev W := N4.V

/-- The ambient `5`-dimensional space for the radical extension. -/
abbrev V := I ⊕ W

/-- The embedded split-support representative
`⟨e₂∧e₃ + e₄∧e₅, e₄∧e₅⟩`. -/
def radSplitRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4.splitRep₁ (k := k))

/-- The second basis vector of the embedded split-support representative. -/
def radSplitRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4.splitRep₂ (k := k))

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The bivector action is linear in the tensor variable. -/
theorem actBivector_add
    (Ω₁ Ω₂ : Matrix V V k) (g : Matrix V V k) :
    ActBivector (Ω₁ + Ω₂) g = ActBivector Ω₁ g + ActBivector Ω₂ g := by
  simp [ActBivector, Matrix.mul_add, Matrix.add_mul]

/-- The bivector action commutes with scalar multiplication in the tensor variable. -/
theorem actBivector_smul
    (c : k) (Ω : Matrix V V k) (g : Matrix V V k) :
    ActBivector (c • Ω) g = c • ActBivector Ω g := by
  simp [ActBivector, Matrix.mul_smul, smul_mul_assoc]

/-- An invertible matrix cannot send a bivector to zero under the natural action. -/
theorem actBivector_eq_zero_iff_of_det_ne_zero
    (Ω : Matrix V V k) (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0) :
    ActBivector Ω g = 0 ↔ Ω = 0 := by
  letI : Invertible (Matrix.det g) := invertibleOfNonzero hg
  letI : Invertible g := Matrix.invertibleOfDetInvertible g
  have hunit : IsUnit (Matrix.det g) := isUnit_of_invertible (Matrix.det g)
  constructor
  · intro hzero
    have h' := congrArg (fun M => g⁻¹ * M * (g⁻¹)ᵀ) hzero
    simp [ActBivector, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hunit,
      Matrix.mul_nonsing_inv _ hunit, Matrix.transpose_nonsing_inv] at h'
    exact h'
  · intro hΩ
    simp [hΩ, ActBivector]

/-- Exact stabilizer of the embedded split-support pair. -/
def FixesRadSplitPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector radSplitRep₁ g ∧ FixesBivector radSplitRep₂ g

/-- Setwise stabilizer of the embedded split-support `2`-space in the chosen basis. -/
def PreservesRadSplitSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector radSplitRep₁ g = α • radSplitRep₁ + β • radSplitRep₂ ∧
    ActBivector radSplitRep₂ g = γ • radSplitRep₁ + δ • radSplitRep₂

/-- The `1 × 1` scalar block in the radical direction. -/
def scalarBlock (a : k) : Matrix I I k := !![a]

@[simp]
lemma scalarBlock_apply (a : k) : scalarBlock (k := k) a 0 0 = a := by
  simp [scalarBlock]

/-- Embedding commutes with linear combinations of the split-support basis vectors. -/
theorem embed_split_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0 (α • (N4.splitRep₁ (k := k)) + β • (N4.splitRep₂ (k := k))) =
      α • radSplitRep₁ + β • radSplitRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [radSplitRep₁, radSplitRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- If a bivector is supported on the lower-right block, then acting by a block upper
triangular matrix with zero upper-right block only changes the lower-right block. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
      Matrix.fromBlocks 0 0 0 (N4.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N4.ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

/-- Pointwise fixing of an embedded bivector reduces to pointwise fixing on the
`4`-dimensional quotient when the upper-right block is zero. -/
theorem fixes_embedded_fromBlocks_zeroUpperRight_iff
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) ↔
      N4.FixesBivector Ω D := by
  constructor
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω at h
    rw [act_embedded_fromBlocks_zeroUpperRight] at h
    have h' := congrArg Matrix.toBlocks₂₂ h
    simpa using h'
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω
    rw [act_embedded_fromBlocks_zeroUpperRight]
    simpa [N4.FixesBivector] using h

/-- The pairwise version of the lower-right reduction for the radical extension. -/
theorem fixesRadSplitPair_fromBlocks_zeroUpperRight_iff
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadSplitPairBivector (Matrix.fromBlocks a 0 C D) ↔
      N4.FixesSplitPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₁ (k := k)) (a := a) (C := C) (D := D)).1 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₂ (k := k)) (a := a) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₁ (k := k)) (a := a) (C := C) (D := D)).2 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.splitRep₂ (k := k)) (a := a) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded pair forces
the lower-right block to fix the underlying `n = 4` split-support pair. -/
theorem fixesRadSplitPair_lowerRight
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadSplitPairBivector (Matrix.fromBlocks a B C D) →
      N4.FixesSplitPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, radSplitRep₁, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, radSplitRep₂, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- Setwise preservation also reduces to the lower-right block when the upper-right block
vanishes. -/
theorem preservesRadSplitSubspace_fromBlocks_zeroUpperRight
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hD : N4.PreservesSplitSubspaceBivector D) :
    PreservesRadSplitSubspaceBivector (Matrix.fromBlocks a 0 C D) := by
  rcases hD with ⟨α, β, γ, δ, h₁, h₂⟩
  refine ⟨α, β, γ, δ, ?_, ?_⟩
  ·
    calc
      ActBivector radSplitRep₁ (Matrix.fromBlocks a 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4.ActBivector (N4.splitRep₁ (k := k)) D) := by
              simpa [radSplitRep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4.splitRep₁ (k := k)) (a := a) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (α • N4.splitRep₁ + β • N4.splitRep₂) := by simp [h₁]
      _ = α • radSplitRep₁ + β • radSplitRep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [radSplitRep₁, radSplitRep₂, Matrix.fromBlocks, Matrix.add_apply]
  ·
    calc
      ActBivector radSplitRep₂ (Matrix.fromBlocks a 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4.ActBivector (N4.splitRep₂ (k := k)) D) := by
              simpa [radSplitRep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4.splitRep₂ (k := k)) (a := a) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (γ • N4.splitRep₁ + δ • N4.splitRep₂) := by simp [h₂]
      _ = γ • radSplitRep₁ + δ • radSplitRep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [radSplitRep₁, radSplitRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The radical extension of the split-support pointwise stabilizer. The radical column is
free, the scalar on the radical line is free, and the lower-right block is the pointwise
stabilizer from the `n = 4` split-support orbit. -/
theorem radSplit_pointwise_family
    (a : k)
    (C : Matrix W I k)
    (A E : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hE : E.det = 1) :
    FixesRadSplitPairBivector
      (Matrix.fromBlocks (scalarBlock (k := k) a) 0 C (Matrix.fromBlocks A 0 0 E)) := by
  have hD : N4.FixesSplitPairBivector (Matrix.fromBlocks A 0 0 E) := by
    exact
      (N4Summary.split_pointwise_bivector_iff (k := k) (A := A) (B := 0) (C := 0) (D := E)).2
        ⟨hA, rfl, rfl, hE⟩
  exact
    (fixesRadSplitPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (a := scalarBlock (k := k) a)
      (C := C)
      (D := Matrix.fromBlocks A 0 0 E)).2 hD

/-- For the radical extension of the split-support orbit, pointwise fixing forces the
upper-right block to vanish. Equivalently, the only freedom outside the `n = 4` split
pointwise stabilizer is the free radical column and the free scalar on the radical line. -/
theorem fixesRadSplitPair_fromBlocks_iff
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadSplitPairBivector (Matrix.fromBlocks a B C D) ↔
      B = 0 ∧ N4.FixesSplitPairBivector D := by
  constructor
  · intro h
    have hD : N4.FixesSplitPairBivector D :=
      fixesRadSplitPair_lowerRight (a := a) (B := B) (C := C) (D := D) h
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
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
      simpa [FixesBivector, radSplitRep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₁tr : B * (N4.splitRep₁ (k := k)) * Dᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
      simpa [FixesBivector, radSplitRep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₂zero : B * (N4.splitRep₂ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have h₁zero : B * (N4.splitRep₁ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hBinl1 : B 0 (Sum.inl 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 0)) h₁zero
      have hneg : -B 0 (Sum.inl 1) = 0 := by
        simpa [N4.splitRep₁, N4.ω12, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      exact neg_eq_zero.mp hneg
    have hBinl0 : B 0 (Sum.inl 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 1)) h₁zero
      simpa [N4.splitRep₁, N4.ω12, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
    have hBinr1 : B 0 (Sum.inr 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 0)) h₂zero
      have hneg : -B 0 (Sum.inr 1) = 0 := by
        simpa [N4.splitRep₂, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      exact neg_eq_zero.mp hneg
    have hBinr0 : B 0 (Sum.inr 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 1)) h₂zero
      simpa [N4.splitRep₂, N4.ω34, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
    refine ⟨?_, hD⟩
    ext i j
    fin_cases i
    rcases j with j | j
    · fin_cases j <;> simp [hBinl0, hBinl1]
    · fin_cases j <;> simp [hBinr0, hBinr1]
  · rintro ⟨hB, hD⟩
    simpa [hB] using
      (fixesRadSplitPair_fromBlocks_zeroUpperRight_iff
        (k := k)
        (a := a)
        (C := C)
        (D := D)).2 hD

/-- The radical extension also preserves the split-support `2`-space setwise whenever the
lower-right block preserves the underlying split-support `n = 4` subspace. -/
theorem radSplit_setwise_family
    (a : k)
    (C : Matrix W I k)
    (A D : Matrix N4.I N4.I k) :
    PreservesRadSplitSubspaceBivector
      (Matrix.fromBlocks (scalarBlock (k := k) a) 0 C (Matrix.fromBlocks A 0 0 D)) := by
  apply preservesRadSplitSubspace_fromBlocks_zeroUpperRight
  exact N4Summary.split_diag_preserves (k := k) (A := A) (D := D)

/-- The standard torus lift on the split-support orbit extends directly to the radical
extension. The radical column does not affect the action on the chosen `2`-space. -/
theorem radSplit_torus_lift_action
    (a u v : k)
    (C : Matrix W I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks A 0 0 D
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) a) 0 C h
    ActBivector radSplitRep₁ g = u • radSplitRep₁ + (v - u) • radSplitRep₂ ∧
      ActBivector radSplitRep₂ g = v • radSplitRep₂ := by
  dsimp
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  have htorus := N4Summary.split_torus_lift_action (k := k) (u := u) (v := v)
  constructor
  ·
    calc
      ActBivector radSplitRep₁
          (Matrix.fromBlocks (scalarBlock (k := k) a) 0 C (Matrix.fromBlocks A 0 0 D)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₁ (k := k)) (Matrix.fromBlocks A 0 0 D)) := by
            simpa [radSplitRep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₁ (k := k))
                (a := scalarBlock (k := k) a)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D)
      _ = Matrix.fromBlocks 0 0 0
          (u • (N4.splitRep₁ (k := k)) + (v - u) • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using htorus.1
      _ = u • radSplitRep₁ + (v - u) • radSplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := u) (β := v - u)
  ·
    calc
      ActBivector radSplitRep₂
          (Matrix.fromBlocks (scalarBlock (k := k) a) 0 C (Matrix.fromBlocks A 0 0 D)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₂ (k := k)) (Matrix.fromBlocks A 0 0 D)) := by
            simpa [radSplitRep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₂ (k := k))
                (a := scalarBlock (k := k) a)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D)
      _ = Matrix.fromBlocks 0 0 0 (v • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using htorus.2
      _ = v • radSplitRep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [radSplitRep₂, Matrix.fromBlocks]

/-- The standard torus lift in the quotient also gives the swap coset after composing with
the split-support involution, and this extends directly to the radical extension. -/
theorem radSplit_swapCoset_lift_action
    (a u v : k)
    (C : Matrix W I k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k)
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) a) 0 C h
    ActBivector radSplitRep₁ g = u • radSplitRep₁ + (v - u) • radSplitRep₂ ∧
      ActBivector radSplitRep₂ g = u • radSplitRep₁ + (-u) • radSplitRep₂ := by
  dsimp
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  have hswap := N4Summary.split_swapCoset_lift_action (k := k) (u := u) (v := v)
  constructor
  ·
    calc
      ActBivector radSplitRep₁
          (Matrix.fromBlocks (scalarBlock (k := k) a) 0 C
            (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₁ (k := k))
            (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) := by
            simpa [radSplitRep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₁ (k := k))
                (a := scalarBlock (k := k) a)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))
      _ = Matrix.fromBlocks 0 0 0
          (u • (N4.splitRep₁ (k := k)) + (v - u) • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using hswap.1
      _ = u • radSplitRep₁ + (v - u) • radSplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := u) (β := v - u)
  ·
    calc
      ActBivector radSplitRep₂
          (Matrix.fromBlocks (scalarBlock (k := k) a) 0 C
            (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.splitRep₂ (k := k))
            (Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))) := by
            simpa [radSplitRep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.splitRep₂ (k := k))
                (a := scalarBlock (k := k) a)
                (C := C)
                (D := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k))
      _ = Matrix.fromBlocks 0 0 0
          (u • (N4.splitRep₁ (k := k)) + (-u) • (N4.splitRep₂ (k := k))) := by
            simpa [A, D] using hswap.2
      _ = u • radSplitRep₁ + (-u) • radSplitRep₂ := by
            exact embed_split_linearCombination (k := k) (α := u) (β := -u)

end N5
end Wedge2Formalization
