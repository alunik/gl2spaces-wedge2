import Wedge2Formalization.N5PureSingularLong
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N6PureSingularLong

variable {k : Type*} [Field k]

/-- The one-dimensional radical factor. -/
abbrev I := Fin 1

/-- The `5`-dimensional quotient carrying the pure singular length-two pair. -/
abbrev W := N5PureSingularLong.V

/-- The ambient `6`-dimensional space. -/
abbrev V := I ⊕ W

/-- Exact stabilizer of a bivector. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The `1 × 1` scalar block in the radical direction. -/
def scalarBlock (u : k) : Matrix I I k := !![u]

/-- The bivector action is multiplicative in the acting matrix. -/
theorem actBivector_mul
    (Ω : Matrix V V k)
    (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  simp [ActBivector, Matrix.mul_assoc]

/-- The embedded pure singular length-two representative in dimension `6`. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N5PureSingularLong.rep₁ (k := k))

/-- The second basis vector of the embedded representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N5PureSingularLong.rep₂ (k := k))

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
      Matrix.fromBlocks 0 0 0 (N5PureSingularLong.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N5PureSingularLong.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- Pointwise fixing of an embedded bivector reduces to the lower-right block when the
upper-right block is zero. -/
theorem fixes_embedded_fromBlocks_zeroUpperRight_iff
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) ↔
      N5PureSingularLong.FixesBivector Ω D := by
  constructor
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω at h
    rw [act_embedded_fromBlocks_zeroUpperRight] at h
    have h' := congrArg Matrix.toBlocks₂₂ h
    simpa [N5PureSingularLong.FixesBivector, N5PureSingularLong.ActBivector] using h'
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω
    rw [act_embedded_fromBlocks_zeroUpperRight]
    simpa [N5PureSingularLong.FixesBivector, N5PureSingularLong.ActBivector] using h

/-- Pairwise lower-right reduction for the radical extension of the pure singular
length-two `5`-dimensional pair. -/
theorem fixesPair_fromBlocks_zeroUpperRight_iff
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesPairBivector (Matrix.fromBlocks A 0 C D) ↔
      N5PureSingularLong.FixesPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5PureSingularLong.rep₁ (k := k)) (A := A) (C := C) (D := D)).1 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5PureSingularLong.rep₂ (k := k)) (A := A) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5PureSingularLong.rep₁ (k := k)) (A := A) (C := C) (D := D)).2 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5PureSingularLong.rep₂ (k := k)) (A := A) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block vanishes, pointwise fixing forces the
lower-right block to fix the underlying pure singular length-two `5`-dimensional pair. -/
theorem fixesPair_lowerRight
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesPairBivector (Matrix.fromBlocks A B C D) →
      N5PureSingularLong.FixesPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  ·
    ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, rep₁, N5PureSingularLong.FixesBivector, N5PureSingularLong.ActBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  ·
    ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, rep₂, N5PureSingularLong.FixesBivector, N5PureSingularLong.ActBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- If the lower-right block already has the internal zero upper-right shape from the
`n = 5` pure singular length-two model, then pointwise fixing is equivalent to the
explicit scalar-plus-Hankel shape on that lower-right block. -/
theorem fixesPair_fromBlocks_nested_zeroUpperRight_iff_shape
    (A0 : Matrix I I k)
    (C0 : Matrix W I k)
    (A1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
    (C1 : Matrix N5PureSingularLong.I N5PureSingularLong.W k)
    (D1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k) :
    FixesPairBivector (Matrix.fromBlocks A0 0 C0 (Matrix.fromBlocks A1 0 C1 D1)) ↔
      A1 0 0 ≠ 0 ∧
        Matrix.fromBlocks A1 0 C1 D1 =
          N5PureSingularLong.pointwiseShape (k := k) (A1 0 0) (C1 0 0) (C1 0 1) (C1 0 2) (C1 1 2) := by
  rw [fixesPair_fromBlocks_zeroUpperRight_iff]
  exact
    N5PureSingularLong.fixesPair_fromBlocks_zeroUpperRight_iff_shape
      (k := k) (A := A1) (C := C1) (D := D1)

/-- If the lower-right block is already written in the internal zero upper-right form
from the `n = 5` pure singular length-two model, then pointwise fixing forces the
radical upper-right row to vanish as well. -/
theorem fixesPair_fromBlocks_nested_iff_shape
    (A0 : Matrix I I k)
    (B0 : Matrix I W k)
    (C0 : Matrix W I k)
    (A1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
    (C1 : Matrix N5PureSingularLong.I N5PureSingularLong.W k)
    (D1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k) :
    FixesPairBivector (Matrix.fromBlocks A0 B0 C0 (Matrix.fromBlocks A1 0 C1 D1)) ↔
      B0 = 0 ∧
        A1 0 0 ≠ 0 ∧
        Matrix.fromBlocks A1 0 C1 D1 =
          N5PureSingularLong.pointwiseShape (k := k) (A1 0 0) (C1 0 0) (C1 0 1) (C1 0 2) (C1 1 2) := by
  constructor
  · intro h
    have hD :
        N5PureSingularLong.FixesPairBivector (Matrix.fromBlocks A1 0 C1 D1) :=
      fixesPair_lowerRight
        (k := k)
        (A := A0)
        (B := B0)
        (C := C0)
        (D := Matrix.fromBlocks A1 0 C1 D1)
        h
    rcases
      (N5PureSingularLong.fixesPair_fromBlocks_zeroUpperRight_iff_shape
        (k := k) (A := A1) (C := C1) (D := D1)).1 hD with
      ⟨hA10, hshape⟩
    rcases h with ⟨h₁, h₂⟩
    have h₁tr : B0 * (N5PureSingularLong.rep₁ (k := k)) * (Matrix.fromBlocks A1 0 C1 D1)ᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
      simpa [FixesBivector, rep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have h₂tr : B0 * (N5PureSingularLong.rep₂ (k := k)) * (Matrix.fromBlocks A1 0 C1 D1)ᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
      simpa [FixesBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have hB30eq : -(B0 0 (Sum.inr 0) * A1 0 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 0)) h₁tr
      rw [hshape] at hij
      simpa [N5PureSingularLong.rep₁, N5PureSingularLong.mulX, N5PureSingularLong.pointwiseShape,
        N5PureSingularLong.lowerHankel, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have hB30mul : B0 0 (Sum.inr 0) * A1 0 0 = 0 := by
      simpa using hB30eq
    have hB30 : B0 0 (Sum.inr 0) = 0 := by
      exact (mul_eq_zero.mp hB30mul).resolve_right hA10
    have hB31eq : -(B0 0 (Sum.inr 1) * A1 0 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 1)) h₁tr
      rw [hshape] at hij
      simpa [N5PureSingularLong.rep₁, N5PureSingularLong.mulX, N5PureSingularLong.pointwiseShape,
        N5PureSingularLong.lowerHankel, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two, Fin.sum_univ_three] using hij
    have hB31mul : B0 0 (Sum.inr 1) * A1 0 0 = 0 := by
      simpa using hB31eq
    have hB31 : B0 0 (Sum.inr 1) = 0 := by
      exact (mul_eq_zero.mp hB31mul).resolve_right hA10
    have hB00eq : B0 0 (Sum.inl 0) * (A1 0 0)⁻¹ = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 0)) h₁tr
      rw [hshape] at hij
      simpa [N5PureSingularLong.rep₁, N5PureSingularLong.mulX, N5PureSingularLong.pointwiseShape,
        N5PureSingularLong.lowerHankel, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two, Fin.sum_univ_three, hB30, hB31] using hij
    have hB00 : B0 0 (Sum.inl 0) = 0 := by
      exact (mul_eq_zero.mp hB00eq).resolve_right (inv_ne_zero hA10)
    have hB01eq : B0 0 (Sum.inl 1) * (A1 0 0)⁻¹ = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 1)) h₁tr
      rw [hshape] at hij
      simpa [N5PureSingularLong.rep₁, N5PureSingularLong.mulX, N5PureSingularLong.pointwiseShape,
        N5PureSingularLong.lowerHankel, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two, Fin.sum_univ_three, hB30, hB31] using hij
    have hB01 : B0 0 (Sum.inl 1) = 0 := by
      exact (mul_eq_zero.mp hB01eq).resolve_right (inv_ne_zero hA10)
    have hB02eq : B0 0 (Sum.inl 2) * (A1 0 0)⁻¹ = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 1)) h₂tr
      rw [hshape] at hij
      simpa [N5PureSingularLong.rep₂, N5PureSingularLong.mulY, N5PureSingularLong.pointwiseShape,
        N5PureSingularLong.lowerHankel, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two, Fin.sum_univ_three, hB30, hB31] using hij
    have hB02 : B0 0 (Sum.inl 2) = 0 := by
      exact (mul_eq_zero.mp hB02eq).resolve_right (inv_ne_zero hA10)
    refine ⟨?_, hA10, hshape⟩
    ext i j
    fin_cases i
    cases j with
    | inl j =>
        fin_cases j <;> simp [hB00, hB01, hB02, hB30, hB31]
    | inr j =>
        fin_cases j <;> simp [hB00, hB01, hB02, hB30, hB31]
  · rintro ⟨hB, hA10, hshape⟩
    rw [hB]
    exact
      (fixesPair_fromBlocks_nested_zeroUpperRight_iff_shape
        (k := k)
        (A0 := A0)
        (C0 := C0)
        (A1 := A1)
        (C1 := C1)
        (D1 := D1)).2 ⟨hA10, hshape⟩

/-- Setwise preservation also reduces to the lower-right block when the upper-right block
vanishes. -/
theorem preservesSubspace_fromBlocks_zeroUpperRight
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hD : N5PureSingularLong.PreservesSubspaceBivector D) :
    PreservesSubspaceBivector (Matrix.fromBlocks A 0 C D) := by
  rcases hD with ⟨α, β, γ, δ, h₁, h₂⟩
  refine ⟨α, β, γ, δ, ?_, ?_⟩
  ·
    calc
      ActBivector rep₁ (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0
          (N5PureSingularLong.ActBivector (N5PureSingularLong.rep₁ (k := k)) D) := by
            simpa [rep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N5PureSingularLong.rep₁ (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0
            (α • N5PureSingularLong.rep₁ + β • N5PureSingularLong.rep₂) := by
              simp [h₁]
      _ = α • rep₁ + β • rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]
  ·
    calc
      ActBivector rep₂ (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0
          (N5PureSingularLong.ActBivector (N5PureSingularLong.rep₂ (k := k)) D) := by
            simpa [rep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (Ω := N5PureSingularLong.rep₂ (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0
            (γ • N5PureSingularLong.rep₁ + δ • N5PureSingularLong.rep₂) := by
              simp [h₂]
      _ = γ • rep₁ + δ • rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The obvious pointwise scalar family on the radical extension of the pure singular
length-two orbit. -/
theorem pointwise_scale_family
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix W I k) :
    FixesPairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)))) := by
  have hD :
      N5PureSingularLong.FixesPairBivector
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) := by
    exact N5PureSingularLong.pointwise_scale_family (k := k) t ht
  exact
    (fixesPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := scalarBlock (k := k) u)
      (C := C)
      (D := Matrix.fromBlocks
        ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
        0
        0
        (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)))).2 hD

/-- The lower-left Hankel unipotent family on the `5`-dimensional quotient also extends
through the radical direction. -/
theorem lowerHankel_pointwise_family
    (u : k)
    (C : Matrix W I k)
    (a b c d : k) :
    FixesPairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) a b c d)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) := by
  have hD :
      N5PureSingularLong.FixesPairBivector
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) a b c d)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)) := by
    exact N5PureSingularLong.lowerHankel_pointwise_family (k := k) a b c d
  exact
    (fixesPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := scalarBlock (k := k) u)
      (C := C)
      (D := Matrix.fromBlocks
        (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
        0
        (N5PureSingularLong.lowerHankel (k := k) a b c d)
        (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))).2 hD

/-- Combining the radical-extension scalar torus and the lifted lower-left Hankel family
gives a concrete `6`-parameter pointwise subgroup. -/
theorem scaleLowerHankel_product_pointwise
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix W I k)
    (a b c d : k) :
    let gScale : Matrix V V k :=
      Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k)))
    let gH : Matrix V V k :=
      Matrix.fromBlocks
        (scalarBlock (k := k) 1)
        0
        0
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) a b c d)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))
    FixesPairBivector (gH * gScale) := by
  intro gScale gH
  have hScale := pointwise_scale_family (k := k) (u := u) (t := t) ht C
  have hH := lowerHankel_pointwise_family (k := k) (u := 1) (C := 0) (a := a) (b := b) (c := c) (d := d)
  constructor <;> rw [FixesBivector]
  ·
    calc
      ActBivector rep₁ (gH * gScale) =
          ActBivector (ActBivector rep₁ gScale) gH := by
            rw [actBivector_mul]
      _ = ActBivector rep₁ gH := by
            simpa [FixesBivector] using congrArg (fun M => ActBivector M gH) hScale.1
      _ = rep₁ := by simpa [FixesBivector] using hH.1
  ·
    calc
      ActBivector rep₂ (gH * gScale) =
          ActBivector (ActBivector rep₂ gScale) gH := by
            rw [actBivector_mul]
      _ = ActBivector rep₂ gH := by
            simpa [FixesBivector] using congrArg (fun M => ActBivector M gH) hScale.2
      _ = rep₂ := by simpa [FixesBivector] using hH.2

/-- The full determinant-scaled `GL₂` quotient lift extends through the radical
direction. -/
theorem GL2_scaled_lift_action
    (u a b c d : k)
    (C : Matrix W I k) :
    let P : Matrix N5PureSingularLong.W N5PureSingularLong.W k :=
      N5PureSingularLong.sym2Raw (k := k) a b c d
    let Q : Matrix N5PureSingularLong.I N5PureSingularLong.I k :=
      N5PureSingularLong.adjointT (k := k) a b c d
    let h : Matrix W W k := Matrix.fromBlocks P 0 0 Q
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) u) 0 C h
    ActBivector rep₁ g =
        N5PureSingularLong.Delta a b c d • (a • rep₁ + c • rep₂) ∧
      ActBivector rep₂ g =
        N5PureSingularLong.Delta a b c d • (b • rep₁ + d • rep₂) := by
  intro P Q h g
  have hD := N5PureSingularLong.GL2_scaled_lift_action (k := k) a b c d
  constructor
  ·
    calc
      ActBivector rep₁ g =
        Matrix.fromBlocks 0 0 0
          (N5PureSingularLong.ActBivector (N5PureSingularLong.rep₁ (k := k)) h) := by
            simpa [g, h, rep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N5PureSingularLong.rep₁ (k := k))
                (A := scalarBlock (k := k) u)
                (C := C)
                (D := h)
      _ = Matrix.fromBlocks 0 0 0
            (N5PureSingularLong.Delta a b c d •
              (a • N5PureSingularLong.rep₁ + c • N5PureSingularLong.rep₂)) := by
              simpa [P, Q, h] using hD.1
      _ = N5PureSingularLong.Delta a b c d • (a • rep₁ + c • rep₂) := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply, smul_add, mul_add]
  ·
    calc
      ActBivector rep₂ g =
        Matrix.fromBlocks 0 0 0
          (N5PureSingularLong.ActBivector (N5PureSingularLong.rep₂ (k := k)) h) := by
            simpa [g, h, rep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N5PureSingularLong.rep₂ (k := k))
                (A := scalarBlock (k := k) u)
                (C := C)
                (D := h)
      _ = Matrix.fromBlocks 0 0 0
            (N5PureSingularLong.Delta a b c d •
              (b • N5PureSingularLong.rep₁ + d • N5PureSingularLong.rep₂)) := by
              simpa [P, Q, h] using hD.2
      _ = N5PureSingularLong.Delta a b c d • (b • rep₁ + d • rep₂) := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply, smul_add, mul_add]

/-- The explicit pointwise subgroup on the radical extension of the pure singular
length-two orbit combines with the determinant-scaled `GL₂` lift exactly as expected. -/
theorem scaleLowerHankel_GL2_product_lift_action
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix W I k)
    (x y z w : k)
    (a b c d : k) :
    let gU : Matrix V V k :=
      (Matrix.fromBlocks
        (scalarBlock (k := k) 1)
        0
        0
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) x y z w)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) *
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))))
    let P : Matrix N5PureSingularLong.W N5PureSingularLong.W k :=
      N5PureSingularLong.sym2Raw (k := k) a b c d
    let Q : Matrix N5PureSingularLong.I N5PureSingularLong.I k :=
      N5PureSingularLong.adjointT (k := k) a b c d
    let h : Matrix W W k := Matrix.fromBlocks P 0 0 Q
    let gL : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) 1) 0 0 h
    ActBivector rep₁ (gU * gL) =
        N5PureSingularLong.Delta a b c d • (a • rep₁ + c • rep₂) ∧
      ActBivector rep₂ (gU * gL) =
        N5PureSingularLong.Delta a b c d • (b • rep₁ + d • rep₂) := by
  intro gU P Q h gL
  have hU :=
    scaleLowerHankel_product_pointwise
      (k := k)
      (u := u)
      (t := t)
      ht
      (C := C)
      (a := x)
      (b := y)
      (c := z)
      (d := w)
  have hL := GL2_scaled_lift_action (k := k) (u := 1) (a := a) (b := b) (c := c) (d := d) (C := 0)
  have hUg1 : ActBivector rep₁ gU = rep₁ := by
    simpa [gU] using hU.1
  have hUg2 : ActBivector rep₂ gU = rep₂ := by
    simpa [gU] using hU.2
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ =
        ActBivector
          (N5PureSingularLong.Delta a b c d • (a • rep₁ + c • rep₂)) gU := by
            simpa [gL, h, P, Q] using congrArg (fun M => ActBivector M gU) hL.1
      _ =
        N5PureSingularLong.Delta a b c d •
          (a • ActBivector rep₁ gU + c • ActBivector rep₂ gU) := by
            simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul,
              smul_add, mul_add]
      _ = N5PureSingularLong.Delta a b c d • (a • rep₁ + c • rep₂) := by
        rw [hUg1, hUg2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ =
        ActBivector
          (N5PureSingularLong.Delta a b c d • (b • rep₁ + d • rep₂)) gU := by
            simpa [gL, h, P, Q] using congrArg (fun M => ActBivector M gU) hL.2
      _ =
        N5PureSingularLong.Delta a b c d •
          (b • ActBivector rep₁ gU + d • ActBivector rep₂ gU) := by
            simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul,
              smul_add, mul_add]
      _ = N5PureSingularLong.Delta a b c d • (b • rep₁ + d • rep₂) := by
        rw [hUg1, hUg2]

/-- Right-multiplying the determinant-scaled `GL₂` lift by the explicit pointwise
subgroup on the radical extension of the pure singular length-two orbit does not change
the quotient action. -/
theorem GL2_scaleLowerHankel_right_product_lift_action
    (u t : k)
    (ht : t ≠ 0)
    (C : Matrix W I k)
    (x y z w a b c d : k) :
    let P : Matrix N5PureSingularLong.W N5PureSingularLong.W k :=
      N5PureSingularLong.sym2Raw (k := k) a b c d
    let Q : Matrix N5PureSingularLong.I N5PureSingularLong.I k :=
      N5PureSingularLong.adjointT (k := k) a b c d
    let h : Matrix W W k := Matrix.fromBlocks P 0 0 Q
    let gL : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) 1) 0 0 h
    let gU : Matrix V V k :=
      (Matrix.fromBlocks
        (scalarBlock (k := k) 1)
        0
        0
        (Matrix.fromBlocks
          (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k)
          0
          (N5PureSingularLong.lowerHankel (k := k) x y z w)
          (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))) *
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          ((t⁻¹) • (1 : Matrix N5PureSingularLong.W N5PureSingularLong.W k))
          0
          0
          (t • (1 : Matrix N5PureSingularLong.I N5PureSingularLong.I k))))
    ActBivector rep₁ (gL * gU) =
        N5PureSingularLong.Delta a b c d • (a • rep₁ + c • rep₂) ∧
      ActBivector rep₂ (gL * gU) =
        N5PureSingularLong.Delta a b c d • (b • rep₁ + d • rep₂) := by
  intro P Q h gL gU
  have hL := GL2_scaled_lift_action (k := k) (u := 1) (a := a) (b := b) (c := c) (d := d) (C := 0)
  have hU :=
    scaleLowerHankel_product_pointwise
      (k := k)
      (u := u)
      (t := t)
      ht
      (C := C)
      (a := x)
      (b := y)
      (c := z)
      (d := w)
  have hUg1 : ActBivector rep₁ gU = rep₁ := by
    simpa [gU] using hU.1
  have hUg2 : ActBivector rep₂ gU = rep₂ := by
    simpa [gU] using hU.2
  constructor
  · calc
      ActBivector rep₁ (gL * gU) = ActBivector (ActBivector rep₁ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hUg1
      _ = N5PureSingularLong.Delta a b c d • (a • rep₁ + c • rep₂) := by
        simpa [gL, h, P, Q] using hL.1
  · calc
      ActBivector rep₂ (gL * gU) = ActBivector (ActBivector rep₂ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hUg2
      _ = N5PureSingularLong.Delta a b c d • (b • rep₁ + d • rep₂) := by
        simpa [gL, h, P, Q] using hL.2

end N6PureSingularLong
end Wedge2Formalization
