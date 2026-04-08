import Wedge2Formalization.N7

open Matrix

namespace Wedge2Formalization
namespace N7OnePoint

variable {k : Type*} [Field k]

/-- The embedded repeated-support representative in dimension `7`. -/
def rad3OnePointRep₁ : Matrix N7.V N7.V k :=
  Matrix.fromBlocks 0 0 0 (N4.onePointRep₁ (k := k))

/-- The second basis vector of the embedded repeated-support representative. -/
def rad3OnePointRep₂ : Matrix N7.V N7.V k :=
  Matrix.fromBlocks 0 0 0 (N4.onePointRep₂ (k := k))

/-- Exact stabilizer of the embedded repeated-support pair. -/
def FixesRad3OnePointPairBivector (g : Matrix N7.V N7.V k) : Prop :=
  N7.FixesBivector rad3OnePointRep₁ g ∧ N7.FixesBivector rad3OnePointRep₂ g

/-- Setwise stabilizer of the embedded repeated-support `2`-space in the chosen basis. -/
def PreservesRad3OnePointSubspaceBivector (g : Matrix N7.V N7.V k) : Prop :=
  ∃ α β γ δ : k,
    N7.ActBivector rad3OnePointRep₁ g = α • rad3OnePointRep₁ + β • rad3OnePointRep₂ ∧
    N7.ActBivector rad3OnePointRep₂ g = γ • rad3OnePointRep₁ + δ • rad3OnePointRep₂

/-- Embedding commutes with linear combinations of the repeated-support basis vectors. -/
theorem embed_onePoint_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0 (α • (N4.onePointRep₁ (k := k)) + β • (N4.onePointRep₂ (k := k))) =
      α • rad3OnePointRep₁ + β • rad3OnePointRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [rad3OnePointRep₁, rad3OnePointRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- Pointwise fixing of the embedded repeated-support pair reduces to pointwise fixing on
the lower-right `4`-dimensional quotient when the upper-right block is zero. -/
theorem fixesRad3OnePointPair_fromBlocks_zeroUpperRight_iff
    (A : Matrix N7.I N7.I k)
    (C : Matrix N7.W N7.I k)
    (D : Matrix N7.W N7.W k) :
    FixesRad3OnePointPairBivector (Matrix.fromBlocks A 0 C D) ↔
      N4.FixesOnePointPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(N7.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₁ (k := k)) (A := A) (C := C) (D := D)).1 h₁,
        (N7.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₂ (k := k)) (A := A) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(N7.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₁ (k := k)) (A := A) (C := C) (D := D)).2 h₁,
        (N7.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₂ (k := k)) (A := A) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded repeated-support
pair forces the lower-right block to fix the underlying `n = 4` repeated-support pair. -/
theorem fixesRad3OnePointPair_lowerRight
    (A : Matrix N7.I N7.I k)
    (B : Matrix N7.I N7.W k)
    (C : Matrix N7.W N7.I k)
    (D : Matrix N7.W N7.W k) :
    FixesRad3OnePointPairBivector (Matrix.fromBlocks A B C D) →
      N4.FixesOnePointPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [N7.FixesBivector, rad3OnePointRep₁, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [N7.FixesBivector, rad3OnePointRep₂, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- For the three-dimensional radical extension of the repeated-support orbit, pointwise
fixing forces the upper-right block to vanish. -/
theorem fixesRad3OnePointPair_fromBlocks_iff
    (A : Matrix N7.I N7.I k)
    (B : Matrix N7.I N7.W k)
    (C : Matrix N7.W N7.I k)
    (D : Matrix N7.W N7.W k) :
    FixesRad3OnePointPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ N4.FixesOnePointPairBivector D := by
  constructor
  · intro h
    have hD : N4.FixesOnePointPairBivector D :=
      fixesRad3OnePointPair_lowerRight (A := A) (B := B) (C := C) (D := D) h
    have hDblocks :
        N4.FixesOnePointPairBivector
          (Matrix.fromBlocks D.toBlocks₁₁ D.toBlocks₁₂ D.toBlocks₂₁ D.toBlocks₂₂) := by
      simpa [Matrix.fromBlocks_toBlocks] using hD
    rcases
        (N4Summary.onePoint_pointwise_bivector_iff
          (k := k)
          (A := D.toBlocks₁₁)
          (B := D.toBlocks₁₂)
          (C := D.toBlocks₂₁)
          (D := D.toBlocks₂₂)).1 hDblocks with
      ⟨hD₁₁, hD₂₁, hD₂₂, hrel⟩
    have hDdet : D.det = 1 := by
      rw [← Matrix.fromBlocks_toBlocks D, hD₂₁, hD₂₂, Matrix.det_fromBlocks_zero₂₁]
      simp [hD₁₁]
    have hunitT : IsUnit ((Dᵀ).det) := by
      simpa [Matrix.det_transpose, hDdet] using (isUnit_one : IsUnit (1 : k))
    rcases h with ⟨h₁, h₂⟩
    have h₂tr : B * (N4.onePointRep₂ (k := k)) * Dᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₂
      simpa [N7.FixesBivector, rad3OnePointRep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₁tr : B * (N4.onePointRep₁ (k := k)) * Dᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₁
      simpa [N7.FixesBivector, rad3OnePointRep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₂zero : B * (N4.onePointRep₂ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have h₁zero : B * (N4.onePointRep₁ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hrow_zero (i : N7.I) : B i = 0 := by
      have hBinl1 : B i (Sum.inl 1) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inl 0)) h₂zero
        have hneg : -B i (Sum.inl 1) = 0 := by
          simpa [N4.onePointRep₂, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
        exact neg_eq_zero.mp hneg
      have hBinl0 : B i (Sum.inl 0) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inl 1)) h₂zero
        simpa [N4.onePointRep₂, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      have hBinr1 : B i (Sum.inr 1) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inl 0)) h₁zero
        have hneg : -B i (Sum.inr 1) = 0 := by
          simpa [N4.onePointRep₁, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
        exact neg_eq_zero.mp hneg
      have hBinr0 : B i (Sum.inr 0) = 0 := by
        have hij := congrArg (fun M => M i (Sum.inl 1)) h₁zero
        simpa [N4.onePointRep₁, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      ext j
      rcases j with j | j
      · fin_cases j <;> simp [hBinl0, hBinl1]
      · fin_cases j <;> simp [hBinr0, hBinr1]
    refine ⟨?_, hD⟩
    ext i j
    simpa using congrArg (fun row => row j) (hrow_zero i)
  · rintro ⟨hB, hD⟩
    simpa [hB] using
      (fixesRad3OnePointPair_fromBlocks_zeroUpperRight_iff
        (k := k)
        (A := A)
        (C := C)
        (D := D)).2 hD

/-- A three-dimensional radical-extension family inside the pointwise stabilizer of the
embedded repeated-support orbit. -/
theorem rad3OnePoint_pointwise_family
    (A0 : Matrix N7.I N7.I k)
    (C : Matrix N7.W N7.I k)
    (A B : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0) :
    FixesRad3OnePointPairBivector
      (Matrix.fromBlocks
        A0
        0
        C
        (Matrix.fromBlocks A B 0 A)) := by
  have hD : N4.FixesOnePointPairBivector (Matrix.fromBlocks A B 0 A) := by
    exact
      (N4Summary.onePoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B)
        (C := 0)
        (D := A)).2 ⟨hA, rfl, rfl, hB⟩
  exact
    (fixesRad3OnePointPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := A0)
      (C := C)
      (D := Matrix.fromBlocks A B 0 A)).2 hD

/-- The standard upper-triangular quotient family on the repeated-support orbit extends
directly to the three-dimensional radical extension. -/
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
    N7.ActBivector rad3OnePointRep₁ g = a • rad3OnePointRep₁ + b • rad3OnePointRep₂ ∧
      N7.ActBivector rad3OnePointRep₂ g = (a * a) • rad3OnePointRep₂ := by
  dsimp
  let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
  have hborel := N4Summary.onePoint_borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  ·
    calc
      N7.ActBivector rad3OnePointRep₁
          (Matrix.fromBlocks
            A0
            0
            C
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.onePointRep₁ (k := k))
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) := by
            simpa [rad3OnePointRep₁] using
              N7.act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.onePointRep₁ (k := k))
                (A := A0)
                (C := C)
                (D := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      _ = Matrix.fromBlocks 0 0 0
          (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k))) := by
            simpa [B] using hborel.1
      _ = a • rad3OnePointRep₁ + b • rad3OnePointRep₂ := by
            exact embed_onePoint_linearCombination (k := k) (α := a) (β := b)
  ·
    calc
      N7.ActBivector rad3OnePointRep₂
          (Matrix.fromBlocks
            A0
            0
            C
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.onePointRep₂ (k := k))
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) := by
            simpa [rad3OnePointRep₂] using
              N7.act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.onePointRep₂ (k := k))
                (A := A0)
                (C := C)
                (D := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      _ = Matrix.fromBlocks 0 0 0 ((a * a) • (N4.onePointRep₂ (k := k))) := by
            simpa [B] using hborel.2
      _ = (a * a) • rad3OnePointRep₂ := by
            simpa [rad3OnePointRep₁, rad3OnePointRep₂] using
              (embed_onePoint_linearCombination (k := k) (α := (0 : k)) (β := a * a))

end N7OnePoint
end Wedge2Formalization
