import Wedge2Formalization.N3PureSingular
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N6PureSingular3

variable {k : Type*} [Field k]

/-- The `3`-dimensional radical factor. -/
abbrev I := Fin 3

/-- The `3`-dimensional quotient carrying the pure singular pair. -/
abbrev W := N3PureSingular.I

/-- The ambient `6`-dimensional space. -/
abbrev V := I ⊕ W

/-- Exact stabilizer of a bivector. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The embedded pure singular representative in dimension `6`. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N3PureSingular.ω13 (k := k))

/-- The second basis vector of the embedded pure singular representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N3PureSingular.ω23 (k := k))

/-- Exact stabilizer of the pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rep₁ g ∧ FixesBivector rep₂ g

/-- Setwise stabilizer of the `2`-space in the chosen basis. -/
def PreservesSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector rep₁ g = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ g = γ • rep₁ + δ • rep₂

/-- Acting on an embedded bivector with zero upper-right block only changes the
lower-right block. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
      Matrix.fromBlocks 0 0 0 (N3PureSingular.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N3PureSingular.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- Pointwise fixing of an embedded bivector reduces to the lower-right block when the
upper-right block is zero. -/
theorem fixes_embedded_fromBlocks_zeroUpperRight_iff
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) ↔
      N3PureSingular.FixesBivector Ω D := by
  constructor
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω at h
    rw [act_embedded_fromBlocks_zeroUpperRight] at h
    have h' := congrArg Matrix.toBlocks₂₂ h
    simpa [N3PureSingular.FixesBivector, N3PureSingular.ActBivector] using h'
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω
    rw [act_embedded_fromBlocks_zeroUpperRight]
    simpa [N3PureSingular.FixesBivector, N3PureSingular.ActBivector] using h

/-- Pairwise lower-right reduction for the radical extension of the pure singular
`3`-dimensional pair. -/
theorem fixesPair_fromBlocks_zeroUpperRight_iff
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesPairBivector (Matrix.fromBlocks A 0 C D) ↔
      N3PureSingular.FixesPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N3PureSingular.ω13 (k := k)) (A := A) (C := C) (D := D)).1 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N3PureSingular.ω23 (k := k)) (A := A) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N3PureSingular.ω13 (k := k)) (A := A) (C := C) (D := D)).2 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N3PureSingular.ω23 (k := k)) (A := A) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block vanishes, pointwise fixing forces the
lower-right block to fix the underlying pure singular `3`-dimensional pair. -/
theorem fixesPair_lowerRight
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesPairBivector (Matrix.fromBlocks A B C D) →
      N3PureSingular.FixesPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  ·
    ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, rep₁, N3PureSingular.FixesBivector, N3PureSingular.ActBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  ·
    ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, rep₂, N3PureSingular.FixesBivector, N3PureSingular.ActBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- Setwise preservation also reduces to the lower-right block when the upper-right block
vanishes. -/
theorem preservesSubspace_fromBlocks_zeroUpperRight
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hD : N3PureSingular.PreservesSubspaceBivector D) :
    PreservesSubspaceBivector (Matrix.fromBlocks A 0 C D) := by
  rcases hD with ⟨α, β, γ, δ, h₁, h₂⟩
  refine ⟨α, β, γ, δ, ?_, ?_⟩
  ·
    calc
      ActBivector rep₁ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) D) := by
              simpa [rep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N3PureSingular.ω13 (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (α • N3PureSingular.ω13 + β • N3PureSingular.ω23) := by simp [h₁]
      _ = α • rep₁ + β • rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]
  ·
    calc
      ActBivector rep₂ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) D) := by
              simpa [rep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N3PureSingular.ω23 (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (γ • N3PureSingular.ω13 + δ • N3PureSingular.ω23) := by simp [h₂]
      _ = γ • rep₁ + δ • rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The obvious pointwise family on the radical extension of the pure singular
`3`-dimensional pair. -/
theorem pointwise_scale_family
    (A0 : Matrix I I k)
    (a : k)
    (ha : a ≠ 0)
    (C : Matrix W I k) :
    FixesPairBivector
      (Matrix.fromBlocks
        A0
        0
        C
        (N3PureSingular.pointwiseScale (k := k) a)) := by
  have hD :
      N3PureSingular.FixesPairBivector
        (N3PureSingular.pointwiseScale (k := k) a) := by
    exact N3PureSingular.pointwise_scale_family (k := k) a ha
  exact
    (fixesPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := A0)
      (C := C)
      (D := N3PureSingular.pointwiseScale (k := k) a)).2 hD

/-- The two-parameter pointwise unipotent family on the `3`-dimensional quotient also
extends through the radical direction. -/
theorem pointwise_unipotent_family
    (A0 : Matrix I I k)
    (C : Matrix W I k)
    (x y : k) :
    FixesPairBivector
      (Matrix.fromBlocks
        A0
        0
        C
        (N3PureSingular.pointwiseUnipotent (k := k) x y)) := by
  have hD :
      N3PureSingular.FixesPairBivector
        (N3PureSingular.pointwiseUnipotent (k := k) x y) := by
    exact N3PureSingular.pointwise_unipotent_family (k := k) x y
  exact
    (fixesPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := A0)
      (C := C)
      (D := N3PureSingular.pointwiseUnipotent (k := k) x y)).2 hD

/-- In the `3 + 3` pure singular row, pointwise fixing forces the upper-right block to
vanish and the lower-right block to have the expected pure singular shape. -/
theorem fixesPair_fromBlocks_iff
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧
      D 0 0 ≠ 0 ∧
      D = N3PureSingular.pureSingularShape (k := k) (D 0 0) (D 2 0) (D 2 1) := by
  constructor
  · intro h
    have hD : N3PureSingular.FixesPairBivector D :=
      fixesPair_lowerRight (k := k) (A := A) (B := B) (C := C) (D := D) h
    rcases
      (N3PureSingular.fixesPairBivector_iff_shape (k := k) (g := D)).1 hD with
      ⟨hD00, hshape⟩
    rcases h with ⟨h₁, h₂⟩
    have h₁tr : B * (N3PureSingular.ω13 (k := k)) * Dᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₁
      simpa [FixesBivector, rep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have h₂tr : B * (N3PureSingular.ω23 (k := k)) * Dᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₂
      simpa [FixesBivector, rep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have hrow (i : I) : B i 0 = 0 ∧ B i 1 = 0 ∧ B i 2 = 0 := by
      have hB2cases : B i 2 = 0 ∨ D 0 0 = 0 := by
        have hij := congrArg (fun M => M i 0) h₁tr
        rw [hshape] at hij
        simpa [N3PureSingular.ω13, N3PureSingular.pureSingularShape, Matrix.mul_apply,
          Fin.sum_univ_three] using hij
      have hB2 : B i 2 = 0 := by
        exact hB2cases.resolve_right hD00
      have hB0cases : B i 0 = 0 ∨ D 0 0 = 0 := by
        have hij := congrArg (fun M => M i 2) h₁tr
        rw [hshape] at hij
        simpa [N3PureSingular.ω13, N3PureSingular.pureSingularShape, Matrix.mul_apply,
          Fin.sum_univ_three, hB2] using hij
      have hB0 : B i 0 = 0 := by
        exact hB0cases.resolve_right hD00
      have hB1cases : B i 1 = 0 ∨ D 0 0 = 0 := by
        have hij := congrArg (fun M => M i 2) h₂tr
        rw [hshape] at hij
        simpa [N3PureSingular.ω23, N3PureSingular.pureSingularShape, Matrix.mul_apply,
          Fin.sum_univ_three, hB2] using hij
      have hB1 : B i 1 = 0 := by
        exact hB1cases.resolve_right hD00
      exact ⟨hB0, hB1, hB2⟩
    have hB00 : B 0 0 = 0 := (hrow 0).1
    have hB01 : B 0 1 = 0 := (hrow 0).2.1
    have hB02 : B 0 2 = 0 := (hrow 0).2.2
    have hB10 : B 1 0 = 0 := (hrow 1).1
    have hB11 : B 1 1 = 0 := (hrow 1).2.1
    have hB12 : B 1 2 = 0 := (hrow 1).2.2
    have hB20 : B 2 0 = 0 := (hrow 2).1
    have hB21 : B 2 1 = 0 := (hrow 2).2.1
    have hB22 : B 2 2 = 0 := (hrow 2).2.2
    refine ⟨?_, hD00, hshape⟩
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [hB00, hB01, hB02, hB10, hB11, hB12, hB20, hB21, hB22]
  · rintro ⟨hB, hD00, hshape⟩
    rw [hB]
    exact
      (fixesPair_fromBlocks_zeroUpperRight_iff (k := k) (A := A) (C := C) (D := D)).2
        ((N3PureSingular.fixesPairBivector_iff_shape (k := k) (g := D)).2 ⟨hD00, hshape⟩)

/-- The full `GL₂` quotient lift extends through the radical direction. -/
theorem GL2_lift_action
    (A0 : Matrix I I k)
    (α β γ δ : k)
    (C : Matrix W I k) :
    let h : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix V V k := Matrix.fromBlocks A0 0 C h
    ActBivector rep₁ g = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ g = γ • rep₁ + δ • rep₂ := by
  intro h g
  have hD := N3PureSingular.GL2_lift_action (k := k) α β γ δ
  constructor
  ·
    calc
      ActBivector rep₁ g = Matrix.fromBlocks 0 0 0 (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) h) := by
        simpa [g, h, rep₁] using
          act_embedded_fromBlocks_zeroUpperRight
            (k := k)
            (Ω := N3PureSingular.ω13 (k := k))
            (A := A0)
            (C := C)
            (D := h)
      _ = Matrix.fromBlocks 0 0 0 (α • N3PureSingular.ω13 + β • N3PureSingular.ω23) := by
        simpa [h] using hD.1
      _ = α • rep₁ + β • rep₂ := by
        ext i j
        cases i <;> cases j <;>
          simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]
  ·
    calc
      ActBivector rep₂ g = Matrix.fromBlocks 0 0 0 (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) h) := by
        simpa [g, h, rep₂] using
          act_embedded_fromBlocks_zeroUpperRight
            (k := k)
            (Ω := N3PureSingular.ω23 (k := k))
            (A := A0)
            (C := C)
            (D := h)
      _ = Matrix.fromBlocks 0 0 0 (γ • N3PureSingular.ω13 + δ • N3PureSingular.ω23) := by
        simpa [h] using hD.2
      _ = γ • rep₁ + δ • rep₂ := by
        ext i j
        cases i <;> cases j <;>
          simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The determinant of the radical-extension `GL₂` lift is the radical determinant
times the `GL₂` determinant. -/
theorem GL2_lift_det
    (A0 : Matrix I I k)
    (α β γ δ : k)
    (C : Matrix W I k) :
    let h : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix V V k := Matrix.fromBlocks A0 0 C h
    Matrix.det g = A0.det * (α * δ - β * γ) := by
  intro h g
  dsimp [g, h]
  rw [Matrix.det_fromBlocks_zero₁₂, N3PureSingular.pureLift_det]
  simp

end N6PureSingular3
end Wedge2Formalization
