import Wedge2Formalization.N6MixedOnePoint
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N7MixedOnePoint

variable {k : Type*} [Field k]

/-- The one-dimensional radical factor used for the mixed one-point `n = 7` orbit. -/
abbrev I := Fin 1

/-- The `6`-dimensional quotient space carrying the mixed one-point `n = 6` orbit. -/
abbrev W := N6MixedOnePoint.V

/-- The ambient `7`-dimensional space. -/
abbrev V := I ⊕ W

/-- The embedded mixed one-point representative in dimension `7`. -/
def radMixedRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k))

/-- The second basis vector of the embedded mixed one-point representative. -/
def radMixedRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k))

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The bivector action of a product factors through successive actions. -/
theorem actBivector_mul
    (Ω : Matrix V V k) (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  rw [ActBivector, ActBivector, ActBivector, Matrix.transpose_mul]
  repeat rw [Matrix.mul_assoc]

/-- Exact stabilizer of the embedded mixed one-point pair. -/
def FixesRadMixedPairBivector (g : Matrix V V k) : Prop :=
  ActBivector (radMixedRep₁ (k := k)) g = radMixedRep₁ ∧
    ActBivector (radMixedRep₂ (k := k)) g = radMixedRep₂

/-- The `1 × 1` scalar block in the radical direction. -/
def scalarBlock (a : k) : Matrix I I k := !![a]

/-- Embedding commutes with linear combinations of the mixed one-point basis vectors. -/
theorem embed_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0
        (α • (N6MixedOnePoint.rep₁ (k := k)) + β • (N6MixedOnePoint.rep₂ (k := k))) =
      α • radMixedRep₁ + β • radMixedRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [radMixedRep₁, radMixedRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The embedded mixed one-point line has the same lower-right determinant as the
underlying `n = 6` orbit. -/
theorem lowerRight_det_in_basis
    (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • radMixedRep₁ + b • radMixedRep₂)) =
      (a * a * a) * (a * a * a) := by
  rw [← embed_linearCombination (k := k) (α := a) (β := b)]
  simp only [Matrix.toBlocks_fromBlocks₂₂]
  exact N6MixedOnePoint.det_in_basis (k := k) (a := a) (b := b)

/-- On the embedded mixed one-point line, the lower-right determinant vanishes exactly
at the repeated-support point. -/
theorem lowerRight_det_zero_iff
    (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • radMixedRep₁ + b • radMixedRep₂)) = (0 : k) ↔
      a = 0 := by
  rw [lowerRight_det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    have hcube : a * a * a = 0 := by
      rcases mul_eq_zero.mp h with h0 | h0 <;> exact h0
    rcases mul_eq_zero.mp hcube with haa | ha
    · rcases mul_eq_zero.mp haa with ha' | ha'
      · exact ha'
      · exact ha'
    · exact ha
  · intro ha
    simp [ha]

/-- If a bivector is supported on the lower-right block, then acting by a block upper
triangular matrix with zero upper-right block only changes the lower-right block. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
      Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N6MixedOnePoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- Pointwise fixing of the embedded mixed one-point pair reduces to pointwise fixing on
the lower-right `6`-dimensional quotient when the upper-right block is zero. -/
theorem fixesRadMixedPair_fromBlocks_zeroUpperRight_iff
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadMixedPairBivector (Matrix.fromBlocks a 0 C D) ↔
      N6MixedOnePoint.FixesPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    refine ⟨?_, ?_⟩
    · change
        ActBivector (Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k)))
            (Matrix.fromBlocks a 0 C D) =
          Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k)) at h₁
      rw [act_embedded_fromBlocks_zeroUpperRight] at h₁
      have h₁' := congrArg Matrix.toBlocks₂₂ h₁
      simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.ActBivector] using h₁'
    · change
        ActBivector (Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k)))
            (Matrix.fromBlocks a 0 C D) =
          Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k)) at h₂
      rw [act_embedded_fromBlocks_zeroUpperRight] at h₂
      have h₂' := congrArg Matrix.toBlocks₂₂ h₂
      simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.ActBivector] using h₂'
  · rintro ⟨h₁, h₂⟩
    refine ⟨?_, ?_⟩
    · change
        ActBivector (Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k)))
            (Matrix.fromBlocks a 0 C D) =
          Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k))
      rw [act_embedded_fromBlocks_zeroUpperRight]
      simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.ActBivector] using h₁
    · change
        ActBivector (Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k)))
            (Matrix.fromBlocks a 0 C D) =
          Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k))
      rw [act_embedded_fromBlocks_zeroUpperRight]
      simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.ActBivector] using h₂

/-- Even without assuming the upper-right block is zero, fixing the embedded pair forces
the lower-right block to fix the underlying `n = 6` mixed one-point pair. -/
theorem fixesRadMixedPair_lowerRight
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadMixedPairBivector (Matrix.fromBlocks a B C D) →
      N6MixedOnePoint.FixesPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [ActBivector, radMixedRep₁, N6MixedOnePoint.ActBivector,
      Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [ActBivector, radMixedRep₂, N6MixedOnePoint.ActBivector,
      Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij

private theorem rep₁_det_ne_zero : Matrix.det (N6MixedOnePoint.rep₁ (k := k)) ≠ 0 := by
  intro hdet
  have h' : (1 : k) = 0 :=
    (N6MixedOnePoint.det_zero_iff (k := k) (a := (1 : k)) (b := (0 : k))).1
      (by simpa using hdet)
  exact one_ne_zero h'

private theorem det_ne_zero_of_pointwise
    (D : Matrix W W k)
    (hfix : N6MixedOnePoint.FixesPairBivector D) :
    Matrix.det D ≠ 0 := by
  intro hD0
  have hdet := congrArg Matrix.det hfix.1
  rw [N6MixedOnePoint.ActBivector, Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose] at hdet
  have hrep0 : Matrix.det (N6MixedOnePoint.rep₁ (k := k)) = 0 := by
    simpa [hD0] using hdet.symm
  exact rep₁_det_ne_zero (k := k) hrep0

/-- For the mixed one-point radical extension, pointwise fixing forces the upper-right
block to vanish. -/
theorem fixesRadMixedPair_fromBlocks_iff
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadMixedPairBivector (Matrix.fromBlocks a B C D) ↔
      B = 0 ∧ N6MixedOnePoint.FixesPairBivector D := by
  constructor
  · intro h
    have hD : N6MixedOnePoint.FixesPairBivector D :=
      fixesRadMixedPair_lowerRight (k := k) (a := a) (B := B) (C := C) (D := D) h
    have hunitD : IsUnit (Matrix.det D) :=
      isUnit_iff_ne_zero.mpr (det_ne_zero_of_pointwise (k := k) (D := D) hD)
    have hunitT : IsUnit ((Dᵀ).det) := by
      simpa [Matrix.det_transpose] using hunitD
    rcases h with ⟨h₁, h₂⟩
    have h₂tr : B * (N6MixedOnePoint.rep₂ (k := k)) * Dᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
      simpa [FixesRadMixedPairBivector, radMixedRep₂, ActBivector,
        Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have h₁tr : B * (N6MixedOnePoint.rep₁ (k := k)) * Dᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
      simpa [FixesRadMixedPairBivector, radMixedRep₁, ActBivector,
        Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have h₂zero : B * (N6MixedOnePoint.rep₂ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have h₁zero : B * (N6MixedOnePoint.rep₁ (k := k)) = 0 := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hBleft1 : B 0 (Sum.inl (Sum.inl 1)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 0))) h₂zero
      have hneg : -B 0 (Sum.inl (Sum.inl 1)) = 0 := by
        simpa [N6MixedOnePoint.rep₂, N4.onePointRep₂, N4.ω12, N4.J, Matrix.mul_apply,
          Fintype.sum_sum_type, Fin.sum_univ_two, Fin.sum_univ_four] using hij
      exact neg_eq_zero.mp hneg
    have hBleft0 : B 0 (Sum.inl (Sum.inl 0)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 1))) h₂zero
      simpa [N6MixedOnePoint.rep₂, N4.onePointRep₂, N4.ω12, N4.J, Matrix.mul_apply,
        Fintype.sum_sum_type, Fin.sum_univ_two, Fin.sum_univ_four] using hij
    have hBmid1 : B 0 (Sum.inl (Sum.inr 1)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 0))) h₁zero
      have hneg : -B 0 (Sum.inl (Sum.inr 1)) = 0 := by
        simpa [N6MixedOnePoint.rep₁, N4.onePointRep₁, N4.J, Matrix.mul_apply,
          Fintype.sum_sum_type, Fin.sum_univ_two, Fin.sum_univ_four] using hij
      exact neg_eq_zero.mp hneg
    have hBmid0 : B 0 (Sum.inl (Sum.inr 0)) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl (Sum.inl 1))) h₁zero
      simpa [N6MixedOnePoint.rep₁, N4.onePointRep₁, N4.J, Matrix.mul_apply,
        Fintype.sum_sum_type, Fin.sum_univ_two, Fin.sum_univ_four] using hij
    have hBright1 : B 0 (Sum.inr 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 0)) h₁zero
      have hneg : -B 0 (Sum.inr 1) = 0 := by
        simpa [N6MixedOnePoint.rep₁, N4.onePointRep₁, N4.J, Matrix.mul_apply,
          Fintype.sum_sum_type, Fin.sum_univ_two, Fin.sum_univ_four] using hij
      exact neg_eq_zero.mp hneg
    have hBright0 : B 0 (Sum.inr 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 1)) h₁zero
      simpa [N6MixedOnePoint.rep₁, N4.onePointRep₁, N4.J, Matrix.mul_apply,
        Fintype.sum_sum_type, Fin.sum_univ_two, Fin.sum_univ_four] using hij
    refine ⟨?_, hD⟩
    ext i j
    fin_cases i
    rcases j with j | j
    · rcases j with j | j
      · fin_cases j <;> simp [hBleft0, hBleft1]
      · fin_cases j <;> simp [hBmid0, hBmid1]
    · fin_cases j <;> simp [hBright0, hBright1]
  · rintro ⟨hB, hD⟩
    simpa [hB] using
      (fixesRadMixedPair_fromBlocks_zeroUpperRight_iff
        (k := k) (a := a) (C := C) (D := D)).2 hD

/-- A one-dimensional radical extension of the mixed `SL₂ × SL₂` pointwise family. -/
theorem radMixed_pointwise_levi_family
    (a0 : Matrix I I k)
    (C : Matrix W I k)
    (A B H : Matrix N4.I N4.I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    FixesRadMixedPairBivector
      (Matrix.fromBlocks a0 0 C (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) := by
  have h6 :
      N6MixedOnePoint.FixesPairBivector
        (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H) :=
    N6MixedOnePoint.pointwise_levi_family
      (k := k) (A := A) (B := B) (H := H) hA hB hH
  constructor
  · calc
      ActBivector (radMixedRep₁ (k := k))
          (Matrix.fromBlocks a0 0 C (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector (N6MixedOnePoint.rep₁ (k := k))
            (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) := by
            simpa [radMixedRep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N6MixedOnePoint.rep₁ (k := k))
                (a := a0)
                (C := C)
                (D := N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
      _ = Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k)) := by
            simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.rep₁] using h6.1
      _ = radMixedRep₁ := by simp [radMixedRep₁]
  · calc
      ActBivector (radMixedRep₂ (k := k))
          (Matrix.fromBlocks a0 0 C (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector (N6MixedOnePoint.rep₂ (k := k))
            (N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)) := by
            simpa [radMixedRep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N6MixedOnePoint.rep₂ (k := k))
                (a := a0)
                (C := C)
                (D := N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H)
      _ = Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k)) := by
            simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.rep₂] using h6.2
      _ = radMixedRep₂ := by simp [radMixedRep₂]

/-- A one-dimensional radical extension of the explicit coupled `E_{12}` family on the
mixed one-point orbit. -/
theorem radMixed_coupledE12_pointwise
    (a0 : Matrix I I k)
    (C : Matrix W I k)
    (t : k) :
    FixesRadMixedPairBivector
      (Matrix.fromBlocks
        a0
        0
        C
        (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ12 (k := k) t)
          (N6MixedOnePoint.coupledR12 (k := k) t)
          1)) := by
  have h6 :
      N6MixedOnePoint.FixesPairBivector
        (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ12 (k := k) t)
          (N6MixedOnePoint.coupledR12 (k := k) t)
          1) :=
    N6MixedOnePoint.coupledE12_pointwise (k := k) t
  constructor
  · calc
      ActBivector (radMixedRep₁ (k := k))
          (Matrix.fromBlocks
            a0
            0
            C
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ12 (k := k) t)
              (N6MixedOnePoint.coupledR12 (k := k) t)
              1)) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector
            (N6MixedOnePoint.rep₁ (k := k))
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ12 (k := k) t)
              (N6MixedOnePoint.coupledR12 (k := k) t)
              1)) := by
            simpa [radMixedRep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N6MixedOnePoint.rep₁ (k := k))
                (a := a0)
                (C := C)
                (D := Matrix.fromBlocks
                  1
                  (N6MixedOnePoint.coupledQ12 (k := k) t)
                  (N6MixedOnePoint.coupledR12 (k := k) t)
                  1)
      _ = Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k)) := by
            simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.rep₁] using h6.1
      _ = radMixedRep₁ := by simp [radMixedRep₁]
  · calc
      ActBivector (radMixedRep₂ (k := k))
          (Matrix.fromBlocks
            a0
            0
            C
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ12 (k := k) t)
              (N6MixedOnePoint.coupledR12 (k := k) t)
              1)) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector
            (N6MixedOnePoint.rep₂ (k := k))
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ12 (k := k) t)
              (N6MixedOnePoint.coupledR12 (k := k) t)
              1)) := by
            simpa [radMixedRep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N6MixedOnePoint.rep₂ (k := k))
                (a := a0)
                (C := C)
                (D := Matrix.fromBlocks
                  1
                  (N6MixedOnePoint.coupledQ12 (k := k) t)
                  (N6MixedOnePoint.coupledR12 (k := k) t)
                  1)
      _ = Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k)) := by
            simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.rep₂] using h6.2
      _ = radMixedRep₂ := by simp [radMixedRep₂]

/-- A one-dimensional radical extension of the explicit coupled `E_{21}` family on the
mixed one-point orbit. -/
theorem radMixed_coupledE21_pointwise
    (a0 : Matrix I I k)
    (C : Matrix W I k)
    (t : k) :
    FixesRadMixedPairBivector
      (Matrix.fromBlocks
        a0
        0
        C
        (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          1)) := by
  have h6 :
      N6MixedOnePoint.FixesPairBivector
        (Matrix.fromBlocks
          1
          (N6MixedOnePoint.coupledQ21 (k := k) t)
          (N6MixedOnePoint.coupledR21 (k := k) t)
          1) :=
    N6MixedOnePoint.coupledE21_pointwise (k := k) t
  constructor
  · calc
      ActBivector (radMixedRep₁ (k := k))
          (Matrix.fromBlocks
            a0
            0
            C
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ21 (k := k) t)
              (N6MixedOnePoint.coupledR21 (k := k) t)
              1)) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector
            (N6MixedOnePoint.rep₁ (k := k))
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ21 (k := k) t)
              (N6MixedOnePoint.coupledR21 (k := k) t)
              1)) := by
            simpa [radMixedRep₁] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N6MixedOnePoint.rep₁ (k := k))
                (a := a0)
                (C := C)
                (D := Matrix.fromBlocks
                  1
                  (N6MixedOnePoint.coupledQ21 (k := k) t)
                  (N6MixedOnePoint.coupledR21 (k := k) t)
                  1)
      _ = Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₁ (k := k)) := by
            simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.rep₁] using h6.1
      _ = radMixedRep₁ := by simp [radMixedRep₁]
  · calc
      ActBivector (radMixedRep₂ (k := k))
          (Matrix.fromBlocks
            a0
            0
            C
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ21 (k := k) t)
              (N6MixedOnePoint.coupledR21 (k := k) t)
              1)) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector
            (N6MixedOnePoint.rep₂ (k := k))
            (Matrix.fromBlocks
              1
              (N6MixedOnePoint.coupledQ21 (k := k) t)
              (N6MixedOnePoint.coupledR21 (k := k) t)
              1)) := by
            simpa [radMixedRep₂] using
              act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := N6MixedOnePoint.rep₂ (k := k))
                (a := a0)
                (C := C)
                (D := Matrix.fromBlocks
                  1
                  (N6MixedOnePoint.coupledQ21 (k := k) t)
                  (N6MixedOnePoint.coupledR21 (k := k) t)
                  1)
      _ = Matrix.fromBlocks 0 0 0 (N6MixedOnePoint.rep₂ (k := k)) := by
            simpa [N6MixedOnePoint.FixesPairBivector, N6MixedOnePoint.rep₂] using h6.2
      _ = radMixedRep₂ := by simp [radMixedRep₂]

/-- The product of the two explicit mixed one-point coupled families still fixes the
embedded pair pointwise. -/
theorem radMixed_coupledE12E21_product_pointwise
    (s t : k) :
    FixesRadMixedPairBivector
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (0 : Matrix W I k)
          (Matrix.fromBlocks
            1
            (N6MixedOnePoint.coupledQ12 (k := k) s)
            (N6MixedOnePoint.coupledR12 (k := k) s)
            1)) *
       (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (0 : Matrix W I k)
          (Matrix.fromBlocks
            1
            (N6MixedOnePoint.coupledQ21 (k := k) t)
            (N6MixedOnePoint.coupledR21 (k := k) t)
            1))) := by
  have h12 :=
    radMixed_coupledE12_pointwise (k := k)
      (a0 := (1 : Matrix I I k))
      (C := (0 : Matrix W I k))
      s
  have h21 :=
    radMixed_coupledE21_pointwise (k := k)
      (a0 := (1 : Matrix I I k))
      (C := (0 : Matrix W I k))
      t
  constructor
  · rw [actBivector_mul]
    rw [h21.1]
    exact h12.1
  · rw [actBivector_mul]
    rw [h21.2]
    exact h12.2

/-- The standard upper-triangular quotient family on the mixed one-point orbit extends
to the one-dimensional radical extension. -/
theorem radMixed_borel_lift_action
    (a0 a b : k)
    (ha : a ≠ 0)
    (C : Matrix W I k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k := Matrix.fromBlocks
      (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      0 0
      (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C H
    ActBivector radMixedRep₁ g = a • radMixedRep₁ + b • radMixedRep₂ ∧
      ActBivector radMixedRep₂ g = (a * a) • radMixedRep₂ := by
  dsimp
  have h6 := N6MixedOnePoint.borel_lift_action (k := k) (a := a) (b := b) ha
  have hblock1 :
      ActBivector radMixedRep₁ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C
        (Matrix.fromBlocks
          (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
            !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
          0 0
          (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector (N6MixedOnePoint.rep₁ (k := k))
            (Matrix.fromBlocks
              (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
                !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
              0 0
              (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) := by
    simpa [radMixedRep₁] using
      act_embedded_fromBlocks_zeroUpperRight
        (k := k)
        (Ω := N6MixedOnePoint.rep₁ (k := k))
        (a := scalarBlock (k := k) a0)
        (C := C)
        (D := Matrix.fromBlocks
          (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
            !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
          0 0
          (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))
  have hblock2 :
      ActBivector radMixedRep₂ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C
        (Matrix.fromBlocks
          (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
            !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
          0 0
          (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector (N6MixedOnePoint.rep₂ (k := k))
            (Matrix.fromBlocks
              (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
                !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
              0 0
              (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) := by
    simpa [radMixedRep₂] using
      act_embedded_fromBlocks_zeroUpperRight
        (k := k)
        (Ω := N6MixedOnePoint.rep₂ (k := k))
        (a := scalarBlock (k := k) a0)
        (C := C)
        (D := Matrix.fromBlocks
          (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
            !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
          0 0
          (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))
  constructor
  · calc
      ActBivector radMixedRep₁ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C
        (Matrix.fromBlocks
          (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
            !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
          0 0
          (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector (N6MixedOnePoint.rep₁ (k := k))
            (Matrix.fromBlocks
              (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
                !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
              0 0
              (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) := by
            simpa using hblock1
      _ = Matrix.fromBlocks 0 0 0
            (a • (N6MixedOnePoint.rep₁ (k := k)) + b • (N6MixedOnePoint.rep₂ (k := k))) := by
            simp [h6.1]
      _ = a • radMixedRep₁ + b • radMixedRep₂ := by
            simpa using (embed_linearCombination (k := k) (α := a) (β := b))
  · calc
      ActBivector radMixedRep₂ (Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C
        (Matrix.fromBlocks
          (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
            !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
          0 0
          (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) =
        Matrix.fromBlocks 0 0 0
          (N6MixedOnePoint.ActBivector (N6MixedOnePoint.rep₂ (k := k))
            (Matrix.fromBlocks
              (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
                !![b * (a⁻¹ * a⁻¹), 0; 0, 0])
              0 0
              (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k))) := by
            simpa using hblock2
      _ = Matrix.fromBlocks 0 0 0 ((a * a) • (N6MixedOnePoint.rep₂ (k := k))) := by
            simp [h6.2]
      _ = (a * a) • radMixedRep₂ := by
            simpa using (embed_linearCombination (k := k) (α := 0) (β := a * a))

/-- The determinant of the standard mixed one-point lift on the `n = 7` radical
extension is the radical scalar times `a^3`. -/
theorem radMixed_borel_lift_det
    (a0 a b : k)
    (C : Matrix W I k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k := Matrix.fromBlocks
      (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B)
      0 0
      (!![a, 0; 0, 1] : Matrix N6MixedOnePoint.I N6MixedOnePoint.I k)
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C H
    Matrix.det g = a0 * (a * a * a) := by
  dsimp
  rw [Matrix.det_fromBlocks_zero₁₂, N6MixedOnePoint.borel_lift_det]
  simp [scalarBlock]

end N7MixedOnePoint
end Wedge2Formalization
