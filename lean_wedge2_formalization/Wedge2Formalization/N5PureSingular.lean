import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N5PureSingular

variable {k : Type*} [Field k]

/-- The one-dimensional radical factor for the pure singular `n = 5` orbit. -/
abbrev I := Fin 1

/-- The `4`-dimensional quotient carrying the pure singular `n = 4` orbit. -/
abbrev W := N4PureSingular.I

/-- The ambient `5`-dimensional space for the radical extension. -/
abbrev V := I ⊕ W

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

/-- The `1 × 1` scalar block in the radical direction. -/
def scalarBlock (a : k) : Matrix I I k := !![a]

/-- A row supported only on the last quotient basis vector. -/
def upperRightLast (x : k) : Matrix I W k := !![(0 : k), 0, 0, x]

/-- The embedded pure singular representative in dimension `5`. -/
def radPureSingularRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k))

/-- The second basis vector of the embedded pure singular representative. -/
def radPureSingularRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k))

/-- Exact stabilizer of the embedded pure singular pair. -/
def FixesRadPureSingularPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector radPureSingularRep₁ g ∧ FixesBivector radPureSingularRep₂ g

/-- Setwise stabilizer of the embedded pure singular `2`-space in the chosen basis. -/
def PreservesRadPureSingularSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector radPureSingularRep₁ g = α • radPureSingularRep₁ + β • radPureSingularRep₂ ∧
    ActBivector radPureSingularRep₂ g = γ • radPureSingularRep₁ + δ • radPureSingularRep₂

/-- Embedding commutes with linear combinations of the pure singular basis vectors. -/
theorem embed_pureSingular_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0 (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) =
      α • radPureSingularRep₁ + β • radPureSingularRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [radPureSingularRep₁, radPureSingularRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- Acting by a block upper triangular matrix with zero upper-right block only changes
the lower-right block on embedded bivectors. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a 0 C D) =
      Matrix.fromBlocks 0 0 0 (N4PureSingular.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N4PureSingular.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- If the upper-right block annihilates an embedded bivector on both sides, then the
action again reduces to the lower-right block. -/
theorem act_embedded_fromBlocks_annihilating_upperRight
    (Ω : Matrix W W k)
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hleft : B * Ω = 0)
    (hright : Ω * Bᵀ = 0) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks a B C D) =
      Matrix.fromBlocks 0 0 0 (N4PureSingular.ActBivector Ω D) := by
  have htop : B * Ω * Bᵀ = 0 := by
    have htmp := congrArg (fun M => M * Bᵀ) hleft
    simpa [Matrix.mul_assoc] using htmp
  have hupper : B * Ω * Dᵀ = 0 := by
    have htmp := congrArg (fun M => M * Dᵀ) hleft
    simpa [Matrix.mul_assoc] using htmp
  have hlower : D * Ω * Bᵀ = 0 := by
    have htmp := congrArg (fun M => D * M) hright
    simpa [Matrix.mul_assoc] using htmp
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N4PureSingular.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply, htop, hupper, hlower]

/-- Pointwise fixing of the embedded pure singular pair reduces to pointwise fixing on
the lower-right `4`-dimensional quotient when the upper-right block is zero. -/
theorem fixesRadPureSingularPair_fromBlocks_zeroUpperRight_iff
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadPureSingularPairBivector (Matrix.fromBlocks a 0 C D) ↔
      N4PureSingular.FixesPureSingularPair D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(by
          have h₁' : ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)))
              (Matrix.fromBlocks a 0 C D) =
            Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)) := by
            simpa [radPureSingularRep₁] using h₁
          rw [act_embedded_fromBlocks_zeroUpperRight] at h₁'
          have h'' := congrArg Matrix.toBlocks₂₂ h₁'
          simpa [N4PureSingular.FixesBivector] using h''),
        (by
          have h₂' : ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)))
              (Matrix.fromBlocks a 0 C D) =
            Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)) := by
            simpa [radPureSingularRep₂] using h₂
          rw [act_embedded_fromBlocks_zeroUpperRight] at h₂'
          have h'' := congrArg Matrix.toBlocks₂₂ h₂'
          simpa [N4PureSingular.FixesBivector] using h'')⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(by
          have h₁' :
              ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)))
                (Matrix.fromBlocks a 0 C D) =
              Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)) := by
            rw [act_embedded_fromBlocks_zeroUpperRight]
            simpa [N4PureSingular.FixesBivector] using h₁
          simpa [radPureSingularRep₁] using h₁'),
        (by
          have h₂' :
              ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)))
                (Matrix.fromBlocks a 0 C D) =
              Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)) := by
            rw [act_embedded_fromBlocks_zeroUpperRight]
            simpa [N4PureSingular.FixesBivector] using h₂
          simpa [radPureSingularRep₂] using h₂')⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded pure singular
pair forces the lower-right block to fix the underlying `n = 4` pure singular pair. -/
theorem fixesRadPureSingularPair_lowerRight
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadPureSingularPairBivector (Matrix.fromBlocks a B C D) →
      N4PureSingular.FixesPureSingularPair D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, radPureSingularRep₁, N4PureSingular.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, radPureSingularRep₂, N4PureSingular.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- For the radical extension of the pure singular orbit, pointwise fixing forces the
upper-right block to be supported only on the last quotient basis vector. -/
theorem fixesRadPureSingularPair_upperRight_shape
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (h : FixesRadPureSingularPairBivector (Matrix.fromBlocks a B C D)) :
    B 0 0 = 0 ∧ B 0 1 = 0 ∧ B 0 2 = 0 := by
  have hD : N4PureSingular.FixesPureSingularPair D :=
    fixesRadPureSingularPair_lowerRight (a := a) (B := B) (C := C) (D := D) h
  rcases
    (N4PureSingular.fixesPureSingularPair_iff_shape (k := k) (g := D)).1 hD with
    ⟨hD00, hshape⟩
  rcases h with ⟨h₁, h₂⟩
  have h₁tr : B * (N4PureSingular.ω12 (k := k)) * Dᵀ = 0 := by
    ext i j
    fin_cases i
    have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
    simpa [FixesBivector, radPureSingularRep₁, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij
  have h₂tr : B * (N4PureSingular.ω13 (k := k)) * Dᵀ = 0 := by
    ext i j
    fin_cases i
    have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
    simpa [FixesBivector, radPureSingularRep₂, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij
  have hB00eq : B 0 0 * (D 0 0)⁻¹ = 0 := by
    have hij := congrArg (fun M => M 0 1) h₁tr
    rw [hshape] at hij
    simpa [N4PureSingular.ω12, N4PureSingular.pureSingularShape,
      Matrix.mul_apply, Fin.sum_univ_four] using hij
  have hB00 : B 0 0 = 0 := by
    exact (mul_eq_zero.mp hB00eq).resolve_right (inv_ne_zero hD00)
  have hB01eq : -(B 0 1) * D 0 0 = 0 := by
    have hij := congrArg (fun M => M 0 0) h₁tr
    rw [hshape] at hij
    simpa [N4PureSingular.ω12, N4PureSingular.pureSingularShape,
      Matrix.mul_apply, Fin.sum_univ_four, hB00] using hij
  have hB01 : B 0 1 = 0 := by
    apply neg_eq_zero.mp
    exact (mul_eq_zero.mp hB01eq).resolve_right hD00
  have hB02eq : -(B 0 2) * D 0 0 = 0 := by
    have hij := congrArg (fun M => M 0 0) h₂tr
    rw [hshape] at hij
    simpa [N4PureSingular.ω13, N4PureSingular.pureSingularShape,
      Matrix.mul_apply, Fin.sum_univ_four, hB00] using hij
  have hB02 : B 0 2 = 0 := by
    apply neg_eq_zero.mp
    exact (mul_eq_zero.mp hB02eq).resolve_right hD00
  exact ⟨hB00, hB01, hB02⟩

lemma upperRightLast_mul_ω12_zero (x : k) :
    upperRightLast (k := k) x * N4PureSingular.ω12 (k := k) = 0 := by
  ext i j
  fin_cases i
  fin_cases j <;>
    simp [upperRightLast, N4PureSingular.ω12, Matrix.mul_apply, Fin.sum_univ_four]

lemma upperRightLast_mul_ω13_zero (x : k) :
    upperRightLast (k := k) x * N4PureSingular.ω13 (k := k) = 0 := by
  ext i j
  fin_cases i
  fin_cases j <;>
    simp [upperRightLast, N4PureSingular.ω13, Matrix.mul_apply, Fin.sum_univ_four]

lemma ω12_mul_upperRightLast_transpose_zero (x : k) :
    N4PureSingular.ω12 (k := k) * (upperRightLast (k := k) x)ᵀ = 0 := by
  have h := congrArg Matrix.transpose (upperRightLast_mul_ω12_zero (k := k) x)
  have htr : (N4PureSingular.ω12 (k := k))ᵀ * (upperRightLast (k := k) x)ᵀ = 0 := by
    simpa [Matrix.transpose_mul] using h
  have hω :
      (N4PureSingular.ω12 (k := k))ᵀ = -(N4PureSingular.ω12 (k := k)) := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp [N4PureSingular.ω12]
  have hneg : -(N4PureSingular.ω12 (k := k) * (upperRightLast (k := k) x)ᵀ) = 0 := by
    simpa [hω, neg_mul] using htr
  exact neg_eq_zero.mp hneg

lemma ω13_mul_upperRightLast_transpose_zero (x : k) :
    N4PureSingular.ω13 (k := k) * (upperRightLast (k := k) x)ᵀ = 0 := by
  have h := congrArg Matrix.transpose (upperRightLast_mul_ω13_zero (k := k) x)
  have htr : (N4PureSingular.ω13 (k := k))ᵀ * (upperRightLast (k := k) x)ᵀ = 0 := by
    simpa [Matrix.transpose_mul] using h
  have hω :
      (N4PureSingular.ω13 (k := k))ᵀ = -(N4PureSingular.ω13 (k := k)) := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp [N4PureSingular.ω13]
  have hneg : -(N4PureSingular.ω13 (k := k) * (upperRightLast (k := k) x)ᵀ) = 0 := by
    simpa [hω, neg_mul] using htr
  exact neg_eq_zero.mp hneg

/-- The free upper-right row in the pure singular radical extension may be supported on
the last quotient basis vector without affecting pointwise fixing. -/
theorem radPureSingular_pointwise_family_with_upperRight
    (a b c d e f t u x : k)
    (ha : a ≠ 0)
    (C : Matrix W I k) :
    FixesRadPureSingularPairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        (upperRightLast (k := k) x)
        C
        (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
  have hfix :=
    N4PureSingular.pureSingularShape_fixes (k := k) a b c d e f t ha
  constructor
  ·
    calc
      ActBivector radPureSingularRep₁
          (Matrix.fromBlocks
            (scalarBlock (k := k) u)
            (upperRightLast (k := k) x)
            C
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector
            (N4PureSingular.ω12 (k := k))
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
              simpa [radPureSingularRep₁] using
                act_embedded_fromBlocks_annihilating_upperRight
                  (Ω := N4PureSingular.ω12 (k := k))
                  (a := scalarBlock (k := k) u)
                  (B := upperRightLast (k := k) x)
                  (C := C)
                  (D := N4PureSingular.pureSingularShape (k := k) a b c d e f t)
                  (hleft := upperRightLast_mul_ω12_zero (k := k) x)
                  (hright := ω12_mul_upperRightLast_transpose_zero (k := k) x)
      _ = Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)) := by
            simpa [N4PureSingular.FixesBivector] using hfix.1
      _ = radPureSingularRep₁ := by simp [radPureSingularRep₁]
  ·
    calc
      ActBivector radPureSingularRep₂
          (Matrix.fromBlocks
            (scalarBlock (k := k) u)
            (upperRightLast (k := k) x)
            C
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector
            (N4PureSingular.ω13 (k := k))
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
              simpa [radPureSingularRep₂] using
                act_embedded_fromBlocks_annihilating_upperRight
                  (Ω := N4PureSingular.ω13 (k := k))
                  (a := scalarBlock (k := k) u)
                  (B := upperRightLast (k := k) x)
                  (C := C)
                  (D := N4PureSingular.pureSingularShape (k := k) a b c d e f t)
                  (hleft := upperRightLast_mul_ω13_zero (k := k) x)
                  (hright := ω13_mul_upperRightLast_transpose_zero (k := k) x)
      _ = Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)) := by
            simpa [N4PureSingular.FixesBivector] using hfix.2
      _ = radPureSingularRep₂ := by simp [radPureSingularRep₂]

/-- Sharp pointwise description for the pure singular radical extension in dimension `5`. -/
theorem fixesRadPureSingularPair_fromBlocks_iff
    (a : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadPureSingularPairBivector (Matrix.fromBlocks a B C D) ↔
      B = upperRightLast (k := k) (B 0 3) ∧
      D 0 0 ≠ 0 ∧
      D =
        N4PureSingular.pureSingularShape
          (k := k) (D 0 0) (D 0 1) (D 0 2) (D 0 3) (D 1 3) (D 2 3) (D 3 3) := by
  constructor
  · intro h
    have hD : N4PureSingular.FixesPureSingularPair D :=
      fixesRadPureSingularPair_lowerRight (a := a) (B := B) (C := C) (D := D) h
    rcases
      (N4PureSingular.fixesPureSingularPair_iff_shape (k := k) (g := D)).1 hD with
      ⟨hD00, hshape⟩
    rcases fixesRadPureSingularPair_upperRight_shape
        (a := a) (B := B) (C := C) (D := D) h with
      ⟨hB00, hB01, hB02⟩
    refine ⟨?_, hD00, hshape⟩
    ext i j
    fin_cases i
    fin_cases j <;> simp [upperRightLast, hB00, hB01, hB02]
  · rintro ⟨hB, hD00, hshape⟩
    rw [hB, hshape]
    have ha : a = scalarBlock (k := k) (a 0 0) := by
      ext i j
      fin_cases i
      fin_cases j
      simp [scalarBlock]
    rw [ha]
    exact
      radPureSingular_pointwise_family_with_upperRight
        (k := k)
        (a := D 0 0)
        (b := D 0 1)
        (c := D 0 2)
        (d := D 0 3)
        (e := D 1 3)
        (f := D 2 3)
        (t := D 3 3)
        (u := a 0 0)
        (x := B 0 3)
        hD00
        C

/-- A radical-extension family inside the pointwise stabilizer of the embedded
pure singular orbit. -/
theorem radPureSingular_pointwise_family
    (a b c d e f t u : k)
    (ha : a ≠ 0)
    (C : Matrix W I k) :
    FixesRadPureSingularPairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
  have hD :
      N4PureSingular.FixesPureSingularPair
        (N4PureSingular.pureSingularShape (k := k) a b c d e f t) :=
    N4PureSingular.pureSingularShape_fixes (k := k) a b c d e f t ha
  exact
    (fixesRadPureSingularPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (a := scalarBlock (k := k) u)
      (C := C)
      (D := N4PureSingular.pureSingularShape (k := k) a b c d e f t)).2 hD

/-- Setwise preservation also reduces to the lower-right block when the upper-right block
vanishes. -/
theorem preservesRadPureSingularSubspace_fromBlocks_zeroUpperRight
    (a : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hD : N4PureSingular.PreservesPureSingularSubspaceBivector D) :
    PreservesRadPureSingularSubspaceBivector (Matrix.fromBlocks a 0 C D) := by
  rcases hD with ⟨α, β, γ, δ, h₁, h₂⟩
  refine ⟨α, β, γ, δ, ?_, ?_⟩
  ·
    calc
      ActBivector radPureSingularRep₁ (Matrix.fromBlocks a 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4PureSingular.ActBivector (N4PureSingular.ω12 (k := k)) D) := by
              simpa [radPureSingularRep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω12 (k := k)) (a := a) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0
          (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) := by simp [h₁]
      _ = α • radPureSingularRep₁ + β • radPureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := α) (β := β)
  ·
    calc
      ActBivector radPureSingularRep₂ (Matrix.fromBlocks a 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4PureSingular.ActBivector (N4PureSingular.ω13 (k := k)) D) := by
              simpa [radPureSingularRep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω13 (k := k)) (a := a) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0
          (γ • (N4PureSingular.ω12 (k := k)) + δ • (N4PureSingular.ω13 (k := k))) := by simp [h₂]
      _ = γ • radPureSingularRep₁ + δ • radPureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := γ) (β := δ)

/-- The pure singular setwise family extends directly to the radical extension. -/
theorem radPureSingular_setwise_family
    (a b c d p q r s e f t u : k)
    (C : Matrix W I k) :
    PreservesRadPureSingularSubspaceBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (N4PureSingular.pureSingularSetwiseShape (k := k) a b c d p q r s e f t)) := by
  apply preservesRadPureSingularSubspace_fromBlocks_zeroUpperRight
  exact N4PureSingular.pureSingularSetwiseShape_preserves (k := k) a b c d p q r s e f t

/-- Every `2 × 2` coefficient matrix on the pure singular orbit lifts through the radical
extension. -/
theorem radPureSingular_GL2_lift_action
    (α β γ δ u : k)
    (C : Matrix W I k) :
    let h :=
      N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) u) 0 C h
    ActBivector radPureSingularRep₁ g =
        α • radPureSingularRep₁ + β • radPureSingularRep₂ ∧
      ActBivector radPureSingularRep₂ g =
        γ • radPureSingularRep₁ + δ • radPureSingularRep₂ := by
  dsimp
  have hlift := N4Summary.pureSingular_GL2_lift_action (k := k) (α := α) (β := β) (γ := γ) (δ := δ)
  constructor
  ·
    calc
      ActBivector radPureSingularRep₁
          (Matrix.fromBlocks
            (scalarBlock (k := k) u)
            0
            C
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector (N4PureSingular.ω12 (k := k))
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) := by
              simpa [radPureSingularRep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω12 (k := k))
                  (a := scalarBlock (k := k) u)
                  (C := C)
                  (D := N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)
      _ = Matrix.fromBlocks 0 0 0
          (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) := by
            simpa using hlift.1
      _ = α • radPureSingularRep₁ + β • radPureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := α) (β := β)
  ·
    calc
      ActBivector radPureSingularRep₂
          (Matrix.fromBlocks
            (scalarBlock (k := k) u)
            0
            C
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector (N4PureSingular.ω13 (k := k))
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) := by
              simpa [radPureSingularRep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω13 (k := k))
                  (a := scalarBlock (k := k) u)
                  (C := C)
                  (D := N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)
      _ = Matrix.fromBlocks 0 0 0
          (γ • (N4PureSingular.ω12 (k := k)) + δ • (N4PureSingular.ω13 (k := k))) := by
            simpa using hlift.2
      _ = γ • radPureSingularRep₁ + δ • radPureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := γ) (β := δ)

end N5PureSingular
end Wedge2Formalization
