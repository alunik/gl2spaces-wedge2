import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N6PureSingular

variable {k : Type*} [Field k]

/-- The two-dimensional radical factor for the pure singular `n = 6` orbit. -/
abbrev I := Fin 2

/-- The `4`-dimensional quotient carrying the pure singular `n = 4` orbit. -/
abbrev W := N4PureSingular.I

/-- The ambient `6`-dimensional space for the radical extension. -/
abbrev V := I ⊕ W

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- A `2 × 4` matrix supported only on the last quotient basis vector. -/
def upperRightLast (x y : k) : Matrix I W k := !![(0 : k), 0, 0, x; 0, 0, 0, y]

/-- The embedded pure singular representative in dimension `6`. -/
def rad2PureSingularRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k))

/-- The second basis vector of the embedded pure singular representative. -/
def rad2PureSingularRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k))

/-- Exact stabilizer of the embedded pure singular pair. -/
def FixesRad2PureSingularPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rad2PureSingularRep₁ g ∧ FixesBivector rad2PureSingularRep₂ g

/-- Setwise stabilizer of the embedded pure singular `2`-space in the chosen basis. -/
def PreservesRad2PureSingularSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector rad2PureSingularRep₁ g = α • rad2PureSingularRep₁ + β • rad2PureSingularRep₂ ∧
    ActBivector rad2PureSingularRep₂ g = γ • rad2PureSingularRep₁ + δ • rad2PureSingularRep₂

/-- Embedding commutes with linear combinations of the pure singular basis vectors. -/
theorem embed_pureSingular_linearCombination
    (α β : k) :
    Matrix.fromBlocks 0 0 0 (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) =
      α • rad2PureSingularRep₁ + β • rad2PureSingularRep₂ := by
  ext i j
  cases i <;> cases j <;>
    simp [rad2PureSingularRep₁, rad2PureSingularRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- Acting by a block upper triangular matrix with zero upper-right block only changes
the lower-right block on embedded bivectors. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
      Matrix.fromBlocks 0 0 0 (N4PureSingular.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N4PureSingular.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- If the upper-right block annihilates an embedded bivector on both sides, then the
action again reduces to the lower-right block. -/
theorem act_embedded_fromBlocks_annihilating_upperRight
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hleft : B * Ω = 0)
    (hright : Ω * Bᵀ = 0) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A B C D) =
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
theorem fixesRad2PureSingularPair_fromBlocks_zeroUpperRight_iff
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRad2PureSingularPairBivector (Matrix.fromBlocks A 0 C D) ↔
      N4PureSingular.FixesPureSingularPair D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(by
          have h₁' :
              ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)))
                (Matrix.fromBlocks A 0 C D) =
              Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)) := by
            simpa [rad2PureSingularRep₁] using h₁
          rw [act_embedded_fromBlocks_zeroUpperRight] at h₁'
          have h'' := congrArg Matrix.toBlocks₂₂ h₁'
          simpa [N4PureSingular.FixesBivector] using h''),
        (by
          have h₂' :
              ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)))
                (Matrix.fromBlocks A 0 C D) =
              Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)) := by
            simpa [rad2PureSingularRep₂] using h₂
          rw [act_embedded_fromBlocks_zeroUpperRight] at h₂'
          have h'' := congrArg Matrix.toBlocks₂₂ h₂'
          simpa [N4PureSingular.FixesBivector] using h'')⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(by
          have h₁' :
              ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)))
                (Matrix.fromBlocks A 0 C D) =
              Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)) := by
            rw [act_embedded_fromBlocks_zeroUpperRight]
            simpa [N4PureSingular.FixesBivector] using h₁
          simpa [rad2PureSingularRep₁] using h₁'),
        (by
          have h₂' :
              ActBivector (Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)))
                (Matrix.fromBlocks A 0 C D) =
              Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)) := by
            rw [act_embedded_fromBlocks_zeroUpperRight]
            simpa [N4PureSingular.FixesBivector] using h₂
          simpa [rad2PureSingularRep₂] using h₂')⟩

/-- Even without assuming the upper-right block is zero, fixing the embedded pure singular
pair forces the lower-right block to fix the underlying `n = 4` pure singular pair. -/
theorem fixesRad2PureSingularPair_lowerRight
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRad2PureSingularPairBivector (Matrix.fromBlocks A B C D) →
      N4PureSingular.FixesPureSingularPair D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, rad2PureSingularRep₁, N4PureSingular.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, rad2PureSingularRep₂, N4PureSingular.FixesBivector,
      Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- For the radical extension of the pure singular orbit, pointwise fixing forces the
upper-right block to be supported only on the last quotient basis vector. -/
theorem fixesRad2PureSingularPair_upperRight_shape
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (h : FixesRad2PureSingularPairBivector (Matrix.fromBlocks A B C D)) :
    B 0 0 = 0 ∧ B 0 1 = 0 ∧ B 0 2 = 0 ∧
    B 1 0 = 0 ∧ B 1 1 = 0 ∧ B 1 2 = 0 := by
  have hD : N4PureSingular.FixesPureSingularPair D :=
    fixesRad2PureSingularPair_lowerRight (A := A) (B := B) (C := C) (D := D) h
  rcases
    (N4PureSingular.fixesPureSingularPair_iff_shape (k := k) (g := D)).1 hD with
    ⟨hD00, hshape⟩
  rcases h with ⟨h₁, h₂⟩
  have h₁tr : B * (N4PureSingular.ω12 (k := k)) * Dᵀ = 0 := by
    ext i j
    have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₁
    simpa [FixesBivector, rad2PureSingularRep₁, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij
  have h₂tr : B * (N4PureSingular.ω13 (k := k)) * Dᵀ = 0 := by
    ext i j
    have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₂
    simpa [FixesBivector, rad2PureSingularRep₂, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply] using hij
  have hrow_shape (i : I) : B i 0 = 0 ∧ B i 1 = 0 ∧ B i 2 = 0 := by
    have hB00eq : B i 0 * (D 0 0)⁻¹ = 0 := by
      have hij := congrArg (fun M => M i 1) h₁tr
      rw [hshape] at hij
      simpa [N4PureSingular.ω12, N4PureSingular.pureSingularShape,
        Matrix.mul_apply, Fin.sum_univ_four] using hij
    have hB00 : B i 0 = 0 := by
      exact (mul_eq_zero.mp hB00eq).resolve_right (inv_ne_zero hD00)
    have hB01eq : -(B i 1) * D 0 0 = 0 := by
      have hij := congrArg (fun M => M i 0) h₁tr
      rw [hshape] at hij
      simpa [N4PureSingular.ω12, N4PureSingular.pureSingularShape,
        Matrix.mul_apply, Fin.sum_univ_four, hB00] using hij
    have hB01 : B i 1 = 0 := by
      apply neg_eq_zero.mp
      exact (mul_eq_zero.mp hB01eq).resolve_right hD00
    have hB02eq : -(B i 2) * D 0 0 = 0 := by
      have hij := congrArg (fun M => M i 0) h₂tr
      rw [hshape] at hij
      simpa [N4PureSingular.ω13, N4PureSingular.pureSingularShape,
        Matrix.mul_apply, Fin.sum_univ_four, hB00] using hij
    have hB02 : B i 2 = 0 := by
      apply neg_eq_zero.mp
      exact (mul_eq_zero.mp hB02eq).resolve_right hD00
    exact ⟨hB00, hB01, hB02⟩
  rcases hrow_shape 0 with ⟨h000, h001, h002⟩
  rcases hrow_shape 1 with ⟨h100, h101, h102⟩
  exact ⟨h000, h001, h002, h100, h101, h102⟩

lemma upperRightLast_mul_ω12_zero (x y : k) :
    upperRightLast (k := k) x y * N4PureSingular.ω12 (k := k) = 0 := by
  ext i j <;>
    fin_cases i <;>
    fin_cases j <;>
    simp [upperRightLast, N4PureSingular.ω12, Matrix.mul_apply, Fin.sum_univ_four]

lemma upperRightLast_mul_ω13_zero (x y : k) :
    upperRightLast (k := k) x y * N4PureSingular.ω13 (k := k) = 0 := by
  ext i j <;>
    fin_cases i <;>
    fin_cases j <;>
    simp [upperRightLast, N4PureSingular.ω13, Matrix.mul_apply, Fin.sum_univ_four]

lemma ω12_mul_upperRightLast_transpose_zero (x y : k) :
    N4PureSingular.ω12 (k := k) * (upperRightLast (k := k) x y)ᵀ = 0 := by
  have h := congrArg Matrix.transpose (upperRightLast_mul_ω12_zero (k := k) x y)
  have htr : (N4PureSingular.ω12 (k := k))ᵀ * (upperRightLast (k := k) x y)ᵀ = 0 := by
    simpa [Matrix.transpose_mul] using h
  have hω :
      (N4PureSingular.ω12 (k := k))ᵀ = -(N4PureSingular.ω12 (k := k)) := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp [N4PureSingular.ω12]
  have hneg : -(N4PureSingular.ω12 (k := k) * (upperRightLast (k := k) x y)ᵀ) = 0 := by
    simpa [hω, neg_mul] using htr
  exact neg_eq_zero.mp hneg

lemma ω13_mul_upperRightLast_transpose_zero (x y : k) :
    N4PureSingular.ω13 (k := k) * (upperRightLast (k := k) x y)ᵀ = 0 := by
  have h := congrArg Matrix.transpose (upperRightLast_mul_ω13_zero (k := k) x y)
  have htr : (N4PureSingular.ω13 (k := k))ᵀ * (upperRightLast (k := k) x y)ᵀ = 0 := by
    simpa [Matrix.transpose_mul] using h
  have hω :
      (N4PureSingular.ω13 (k := k))ᵀ = -(N4PureSingular.ω13 (k := k)) := by
    ext i j
    fin_cases i <;> fin_cases j <;> simp [N4PureSingular.ω13]
  have hneg : -(N4PureSingular.ω13 (k := k) * (upperRightLast (k := k) x y)ᵀ) = 0 := by
    simpa [hω, neg_mul] using htr
  exact neg_eq_zero.mp hneg

/-- The free upper-right block in the pure singular radical extension may be supported on
the last quotient basis vector without affecting pointwise fixing. -/
theorem rad2PureSingular_pointwise_family_with_upperRight
    (a b c d e f t x y : k)
    (ha : a ≠ 0)
    (A0 : Matrix I I k)
    (C : Matrix W I k) :
    FixesRad2PureSingularPairBivector
      (Matrix.fromBlocks
        A0
        (upperRightLast (k := k) x y)
        C
        (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
  have hfix :=
    N4PureSingular.pureSingularShape_fixes (k := k) a b c d e f t ha
  constructor
  ·
    calc
      ActBivector rad2PureSingularRep₁
          (Matrix.fromBlocks
            A0
            (upperRightLast (k := k) x y)
            C
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector
            (N4PureSingular.ω12 (k := k))
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
              simpa [rad2PureSingularRep₁] using
                act_embedded_fromBlocks_annihilating_upperRight
                  (Ω := N4PureSingular.ω12 (k := k))
                  (A := A0)
                  (B := upperRightLast (k := k) x y)
                  (C := C)
                  (D := N4PureSingular.pureSingularShape (k := k) a b c d e f t)
                  (hleft := upperRightLast_mul_ω12_zero (k := k) x y)
                  (hright := ω12_mul_upperRightLast_transpose_zero (k := k) x y)
      _ = Matrix.fromBlocks 0 0 0 (N4PureSingular.ω12 (k := k)) := by
            simpa [N4PureSingular.FixesBivector] using hfix.1
      _ = rad2PureSingularRep₁ := by simp [rad2PureSingularRep₁]
  ·
    calc
      ActBivector rad2PureSingularRep₂
          (Matrix.fromBlocks
            A0
            (upperRightLast (k := k) x y)
            C
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector
            (N4PureSingular.ω13 (k := k))
            (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
              simpa [rad2PureSingularRep₂] using
                act_embedded_fromBlocks_annihilating_upperRight
                  (Ω := N4PureSingular.ω13 (k := k))
                  (A := A0)
                  (B := upperRightLast (k := k) x y)
                  (C := C)
                  (D := N4PureSingular.pureSingularShape (k := k) a b c d e f t)
                  (hleft := upperRightLast_mul_ω13_zero (k := k) x y)
                  (hright := ω13_mul_upperRightLast_transpose_zero (k := k) x y)
      _ = Matrix.fromBlocks 0 0 0 (N4PureSingular.ω13 (k := k)) := by
            simpa [N4PureSingular.FixesBivector] using hfix.2
      _ = rad2PureSingularRep₂ := by simp [rad2PureSingularRep₂]

/-- Sharp pointwise description for the pure singular radical extension in dimension `6`. -/
theorem fixesRad2PureSingularPair_fromBlocks_iff
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRad2PureSingularPairBivector (Matrix.fromBlocks A B C D) ↔
      B = upperRightLast (k := k) (B 0 3) (B 1 3) ∧
      D 0 0 ≠ 0 ∧
      D =
        N4PureSingular.pureSingularShape
          (k := k) (D 0 0) (D 0 1) (D 0 2) (D 0 3) (D 1 3) (D 2 3) (D 3 3) := by
  constructor
  · intro h
    have hD : N4PureSingular.FixesPureSingularPair D :=
      fixesRad2PureSingularPair_lowerRight (A := A) (B := B) (C := C) (D := D) h
    rcases
      (N4PureSingular.fixesPureSingularPair_iff_shape (k := k) (g := D)).1 hD with
      ⟨hD00, hshape⟩
    rcases fixesRad2PureSingularPair_upperRight_shape
        (A := A) (B := B) (C := C) (D := D) h with
      ⟨h000, h001, h002, h100, h101, h102⟩
    refine ⟨?_, hD00, hshape⟩
    ext i j
    fin_cases i <;> fin_cases j <;> simp [upperRightLast, h000, h001, h002, h100, h101, h102]
  · rintro ⟨hB, hD00, hshape⟩
    rw [hB, hshape]
    exact
      rad2PureSingular_pointwise_family_with_upperRight
        (k := k)
        (a := D 0 0)
        (b := D 0 1)
        (c := D 0 2)
        (d := D 0 3)
        (e := D 1 3)
        (f := D 2 3)
        (t := D 3 3)
        (x := B 0 3)
        (y := B 1 3)
        hD00
        A
        C

/-- A radical-extension family inside the pointwise stabilizer of the embedded
pure singular orbit. -/
theorem rad2PureSingular_pointwise_family
    (a b c d e f t : k)
    (ha : a ≠ 0)
    (A0 : Matrix I I k)
    (C : Matrix W I k) :
    FixesRad2PureSingularPairBivector
      (Matrix.fromBlocks
        A0
        0
        C
        (N4PureSingular.pureSingularShape (k := k) a b c d e f t)) := by
  have hD :
      N4PureSingular.FixesPureSingularPair
        (N4PureSingular.pureSingularShape (k := k) a b c d e f t) :=
    N4PureSingular.pureSingularShape_fixes (k := k) a b c d e f t ha
  exact
    (fixesRad2PureSingularPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := A0)
      (C := C)
      (D := N4PureSingular.pureSingularShape (k := k) a b c d e f t)).2 hD

/-- Setwise preservation also reduces to the lower-right block when the upper-right block
vanishes. -/
theorem preservesRad2PureSingularSubspace_fromBlocks_zeroUpperRight
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hD : N4PureSingular.PreservesPureSingularSubspaceBivector D) :
    PreservesRad2PureSingularSubspaceBivector (Matrix.fromBlocks A 0 C D) := by
  rcases hD with ⟨α, β, γ, δ, h₁, h₂⟩
  refine ⟨α, β, γ, δ, ?_, ?_⟩
  ·
    calc
      ActBivector rad2PureSingularRep₁ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4PureSingular.ActBivector (N4PureSingular.ω12 (k := k)) D) := by
              simpa [rad2PureSingularRep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω12 (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0
          (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) := by simp [h₁]
      _ = α • rad2PureSingularRep₁ + β • rad2PureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := α) (β := β)
  ·
    calc
      ActBivector rad2PureSingularRep₂ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N4PureSingular.ActBivector (N4PureSingular.ω13 (k := k)) D) := by
              simpa [rad2PureSingularRep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω13 (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0
          (γ • (N4PureSingular.ω12 (k := k)) + δ • (N4PureSingular.ω13 (k := k))) := by simp [h₂]
      _ = γ • rad2PureSingularRep₁ + δ • rad2PureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := γ) (β := δ)

/-- The pure singular setwise family extends directly to the two-dimensional radical
extension. -/
theorem rad2PureSingular_setwise_family
    (a b c d p q r s e f t : k)
    (A0 : Matrix I I k)
    (C : Matrix W I k) :
    PreservesRad2PureSingularSubspaceBivector
      (Matrix.fromBlocks
        A0
        0
        C
        (N4PureSingular.pureSingularSetwiseShape (k := k) a b c d p q r s e f t)) := by
  apply preservesRad2PureSingularSubspace_fromBlocks_zeroUpperRight
  exact N4PureSingular.pureSingularSetwiseShape_preserves (k := k) a b c d p q r s e f t

/-- Every `2 × 2` coefficient matrix on the pure singular orbit lifts through the
two-dimensional radical extension. -/
theorem rad2PureSingular_GL2_lift_action
    (α β γ δ : k)
    (A0 : Matrix I I k)
    (C : Matrix W I k) :
    let h :=
      N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1
    let g : Matrix V V k := Matrix.fromBlocks A0 0 C h
    ActBivector rad2PureSingularRep₁ g =
        α • rad2PureSingularRep₁ + β • rad2PureSingularRep₂ ∧
      ActBivector rad2PureSingularRep₂ g =
        γ • rad2PureSingularRep₁ + δ • rad2PureSingularRep₂ := by
  dsimp
  have hlift := N4Summary.pureSingular_GL2_lift_action (k := k) (α := α) (β := β) (γ := γ) (δ := δ)
  constructor
  ·
    calc
      ActBivector rad2PureSingularRep₁
          (Matrix.fromBlocks
            A0
            0
            C
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector (N4PureSingular.ω12 (k := k))
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) := by
              simpa [rad2PureSingularRep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω12 (k := k))
                  (A := A0)
                  (C := C)
                  (D := N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)
      _ = Matrix.fromBlocks 0 0 0
          (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) := by
            simpa using hlift.1
      _ = α • rad2PureSingularRep₁ + β • rad2PureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := α) (β := β)
  ·
    calc
      ActBivector rad2PureSingularRep₂
          (Matrix.fromBlocks
            A0
            0
            C
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) =
        Matrix.fromBlocks 0 0 0
          (N4PureSingular.ActBivector (N4PureSingular.ω13 (k := k))
            (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)) := by
              simpa [rad2PureSingularRep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N4PureSingular.ω13 (k := k))
                  (A := A0)
                  (C := C)
                  (D := N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1)
      _ = Matrix.fromBlocks 0 0 0
          (γ • (N4PureSingular.ω12 (k := k)) + δ • (N4PureSingular.ω13 (k := k))) := by
            simpa using hlift.2
      _ = γ • rad2PureSingularRep₁ + δ • rad2PureSingularRep₂ := by
            exact embed_pureSingular_linearCombination (k := k) (α := γ) (β := δ)

end N6PureSingular
end Wedge2Formalization
