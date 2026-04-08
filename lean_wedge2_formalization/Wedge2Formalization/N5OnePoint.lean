import Wedge2Formalization.N5

open Matrix

namespace Wedge2Formalization
namespace N5OnePoint

variable {k : Type*} [Field k]

/-- The embedded repeated-support representative in dimension `5`. -/
def radOnePointRep₁ : Matrix N5.V N5.V k :=
  Matrix.fromBlocks 0 0 0 (N4.onePointRep₁ (k := k))

/-- The second basis vector of the embedded repeated-support representative. -/
def radOnePointRep₂ : Matrix N5.V N5.V k :=
  Matrix.fromBlocks 0 0 0 (N4.onePointRep₂ (k := k))

/-- Exact stabilizer of the embedded repeated-support pair. -/
def FixesRadOnePointPairBivector (g : Matrix N5.V N5.V k) : Prop :=
  N5.FixesBivector radOnePointRep₁ g ∧ N5.FixesBivector radOnePointRep₂ g

/-- Setwise stabilizer of the embedded repeated-support `2`-space in the chosen basis. -/
def PreservesRadOnePointSubspaceBivector (g : Matrix N5.V N5.V k) : Prop :=
  ∃ α β γ δ : k,
    N5.ActBivector radOnePointRep₁ g = α • radOnePointRep₁ + β • radOnePointRep₂ ∧
    N5.ActBivector radOnePointRep₂ g = γ • radOnePointRep₁ + δ • radOnePointRep₂

/-- Embedding commutes with linear combinations of the repeated-support basis vectors. -/
theorem embed_onePoint_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0 (α • (N4.onePointRep₁ (k := k)) + β • (N4.onePointRep₂ (k := k))) =
      α • radOnePointRep₁ + β • radOnePointRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [radOnePointRep₁, radOnePointRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- Pointwise fixing of the embedded repeated-support pair reduces to pointwise fixing on
the lower-right `4`-dimensional quotient when the upper-right block is zero. -/
theorem fixesRadOnePointPair_fromBlocks_zeroUpperRight_iff
    (a : Matrix N5.I N5.I k)
    (C : Matrix N5.W N5.I k)
    (D : Matrix N5.W N5.W k) :
    FixesRadOnePointPairBivector (Matrix.fromBlocks a 0 C D) ↔
      N4.FixesOnePointPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(N5.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₁ (k := k)) (a := a) (C := C) (D := D)).1 h₁,
        (N5.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₂ (k := k)) (a := a) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(N5.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₁ (k := k)) (a := a) (C := C) (D := D)).2 h₁,
        (N5.fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N4.onePointRep₂ (k := k)) (a := a) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded repeated-support
pair forces the lower-right block to fix the underlying `n = 4` repeated-support pair. -/
theorem fixesRadOnePointPair_lowerRight
    (a : Matrix N5.I N5.I k)
    (B : Matrix N5.I N5.W k)
    (C : Matrix N5.W N5.I k)
    (D : Matrix N5.W N5.W k) :
    FixesRadOnePointPairBivector (Matrix.fromBlocks a B C D) →
      N4.FixesOnePointPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [N5.FixesBivector, radOnePointRep₁, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [N5.FixesBivector, radOnePointRep₂, N4.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- For the radical extension of the repeated-support orbit, pointwise fixing forces the
upper-right block to vanish. Equivalently, the only freedom outside the `n = 4`
repeated-support pointwise stabilizer is the free radical column and the free scalar on the
radical line. -/
theorem fixesRadOnePointPair_fromBlocks_iff
    (a : Matrix N5.I N5.I k)
    (B : Matrix N5.I N5.W k)
    (C : Matrix N5.W N5.I k)
    (D : Matrix N5.W N5.W k) :
    FixesRadOnePointPairBivector (Matrix.fromBlocks a B C D) ↔
      B = 0 ∧ N4.FixesOnePointPairBivector D := by
  constructor
  · intro h
    have hD : N4.FixesOnePointPairBivector D :=
      fixesRadOnePointPair_lowerRight (a := a) (B := B) (C := C) (D := D) h
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
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
      simpa [N5.FixesBivector, radOnePointRep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₁tr : B * (N4.onePointRep₁ (k := k)) * Dᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
      simpa [N5.FixesBivector, radOnePointRep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₂zero : B * (N4.onePointRep₂ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have h₁zero : B * (N4.onePointRep₁ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hBinl1 : B 0 (Sum.inl 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 0)) h₂zero
      have hneg : -B 0 (Sum.inl 1) = 0 := by
        simpa [N4.onePointRep₂, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      exact neg_eq_zero.mp hneg
    have hBinl0 : B 0 (Sum.inl 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 1)) h₂zero
      simpa [N4.onePointRep₂, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
    have hBinr1 : B 0 (Sum.inr 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 0)) h₁zero
      have hneg : -B 0 (Sum.inr 1) = 0 := by
        simpa [N4.onePointRep₁, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
      exact neg_eq_zero.mp hneg
    have hBinr0 : B 0 (Sum.inr 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 1)) h₁zero
      simpa [N4.onePointRep₁, N4.ω12, Matrix.mul_apply, Fintype.sum_sum_type, N4.J] using hij
    refine ⟨?_, hD⟩
    ext i j
    fin_cases i
    rcases j with j | j
    · fin_cases j <;> simp [hBinl0, hBinl1]
    · fin_cases j <;> simp [hBinr0, hBinr1]
  · rintro ⟨hB, hD⟩
    simpa [hB] using
      (fixesRadOnePointPair_fromBlocks_zeroUpperRight_iff
        (k := k)
        (a := a)
        (C := C)
        (D := D)).2 hD

/-- A radical-extension family inside the pointwise stabilizer of the embedded
repeated-support orbit. -/
theorem radOnePoint_pointwise_family
    (a : k)
    (C : Matrix N5.W N5.I k)
    (A B : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0) :
    FixesRadOnePointPairBivector
      (Matrix.fromBlocks
        (N5.scalarBlock (k := k) a)
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
    (fixesRadOnePointPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (a := N5.scalarBlock (k := k) a)
      (C := C)
      (D := Matrix.fromBlocks A B 0 A)).2 hD

/-- The standard upper-triangular quotient family on the repeated-support orbit extends
directly to the radical extension. -/
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
    N5.ActBivector radOnePointRep₁ g = a • radOnePointRep₁ + b • radOnePointRep₂ ∧
      N5.ActBivector radOnePointRep₂ g = (a * a) • radOnePointRep₂ := by
  dsimp
  let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
  have hborel := N4Summary.onePoint_borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  ·
    calc
      N5.ActBivector radOnePointRep₁
          (Matrix.fromBlocks
            (N5.scalarBlock (k := k) t)
            0
            C
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.onePointRep₁ (k := k))
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) := by
            simpa [radOnePointRep₁] using
              N5.act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.onePointRep₁ (k := k))
                (a := N5.scalarBlock (k := k) t)
                (C := C)
                (D := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      _ = Matrix.fromBlocks 0 0 0
          (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k))) := by
            simpa [B] using hborel.1
      _ = a • radOnePointRep₁ + b • radOnePointRep₂ := by
            exact embed_onePoint_linearCombination (k := k) (α := a) (β := b)
  ·
    calc
      N5.ActBivector radOnePointRep₂
          (Matrix.fromBlocks
            (N5.scalarBlock (k := k) t)
            0
            C
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) =
        Matrix.fromBlocks 0 0 0
          (N4.ActBivector (N4.onePointRep₂ (k := k))
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)) := by
            simpa [radOnePointRep₂] using
              N5.act_embedded_fromBlocks_zeroUpperRight
                (Ω := N4.onePointRep₂ (k := k))
                (a := N5.scalarBlock (k := k) t)
                (C := C)
                (D := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      _ = Matrix.fromBlocks 0 0 0 ((a * a) • (N4.onePointRep₂ (k := k))) := by
            simpa [B] using hborel.2
      _ = (a * a) • radOnePointRep₂ := by
            simpa [radOnePointRep₁, radOnePointRep₂] using
              (embed_onePoint_linearCombination (k := k) (α := (0 : k)) (β := a * a))

end N5OnePoint
end Wedge2Formalization
