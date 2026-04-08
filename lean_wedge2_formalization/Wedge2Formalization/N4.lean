import Mathlib.Data.Matrix.Reflection
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

open Matrix

namespace Wedge2Formalization
namespace N4

variable {k : Type*} [Field k]

abbrev I := Fin 2
abbrev V := I ⊕ I

/-- The standard alternating matrix on a 2-dimensional space. -/
def J : Matrix I I k := !![(0 : k), 1; -1, 0]

/-- The rank-2 tensor `e₁ ∧ e₂`, written in block form. -/
def ω12 : Matrix V V k := Matrix.fromBlocks J 0 0 0

/-- The rank-2 tensor `e₃ ∧ e₄`, written in block form. -/
def ω34 : Matrix V V k := Matrix.fromBlocks 0 0 0 J

/-- The split-support representative `⟨e₁∧e₂ + e₃∧e₄, e₃∧e₄⟩`. -/
def splitRep₁ : Matrix V V k := ω12 + ω34

/-- The second basis vector of the split-support representative. -/
def splitRep₂ : Matrix V V k := ω34

/-- Exact stabilizer of a tensor under congruence. -/
def Fixes (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  gᵀ * Ω * g = Ω

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- Exact stabilizer of the split-support pair. -/
def FixesSplitPair (g : Matrix V V k) : Prop :=
  Fixes splitRep₁ g ∧ Fixes splitRep₂ g

/-- Exact stabilizer of the split-support pair for the action on `∧² V`. -/
def FixesSplitPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector splitRep₁ g ∧ FixesBivector splitRep₂ g

/-- Setwise stabilizer of the split-support `2`-space in the chosen basis. -/
def PreservesSplitSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector splitRep₁ g = α • splitRep₁ + β • splitRep₂ ∧
    ActBivector splitRep₂ g = γ • splitRep₁ + δ • splitRep₂

/-- A repeated-support representative of type `2[a]`. -/
def onePointRep₁ : Matrix V V k := Matrix.fromBlocks 0 J J 0

/-- A second basis vector for the repeated-support representative. -/
def onePointRep₂ : Matrix V V k := ω12

/-- Exact stabilizer of the repeated-support pair. -/
def FixesOnePointPair (g : Matrix V V k) : Prop :=
  Fixes onePointRep₁ g ∧ Fixes onePointRep₂ g

/-- Exact stabilizer of the repeated-support pair for the action on `∧² V`. -/
def FixesOnePointPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector onePointRep₁ g ∧ FixesBivector onePointRep₂ g

/-- Setwise stabilizer of the repeated-support `2`-space in the chosen basis. -/
def PreservesOnePointSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector onePointRep₁ g = α • onePointRep₁ + β • onePointRep₂ ∧
    ActBivector onePointRep₂ g = γ • onePointRep₁ + δ • onePointRep₂

/-- The block-swap matrix exchanging the two `2`-dimensional summands. -/
def splitSwap : Matrix V V k :=
  Matrix.fromBlocks 0 (1 : Matrix I I k) (1 : Matrix I I k) 0

/-- A general anti-diagonal block matrix for the split-support family. -/
def splitAntiDiag (B C : Matrix I I k) : Matrix V V k :=
  Matrix.fromBlocks 0 B C 0

/-- An upper block shear used in the repeated-support setwise stabilizer. -/
def onePointUpperShear (B : Matrix I I k) : Matrix V V k :=
  Matrix.fromBlocks 1 B 0 1

/-- A diagonal scaling used in the repeated-support setwise stabilizer. -/
def onePointScale (a : k) : Matrix V V k :=
  Matrix.fromBlocks (a • (1 : Matrix I I k)) 0 0 1

@[simp]
lemma J_det : (J (k := k)).det = 1 := by
  simp [J, Matrix.det_fin_two]

lemma transpose_mul_J_mul (A : Matrix I I k) :
    Aᵀ * J * A = !![(0 : k), A.det; -A.det, 0] := by
  ext i j
  fin_cases i <;> fin_cases j
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf

lemma mul_J_transpose_mul (A : Matrix I I k) :
    A * J * Aᵀ = !![(0 : k), A.det; -A.det, 0] := by
  ext i j
  fin_cases i <;> fin_cases j
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf
  · simp [J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.det_fin_two]
    ring_nf

lemma mul_J_add_J_mul_transpose (B : Matrix I I k) :
    B * J + J * Bᵀ = !![(0 : k), B 0 0 + B 1 1; -(B 0 0 + B 1 1), 0] := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simp [J, Matrix.mul_apply, Fin.sum_univ_two]
    change -B 0 1 + ((![0, 1] : I → k) ⬝ᵥ fun i => (Bᵀ) i 0) = 0
    simp [Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  ·
    simp [J, Matrix.mul_apply, Fin.sum_univ_two]
    change ((![0, 1] : I → k) ⬝ᵥ fun i => (Bᵀ) i 1) = B 1 1
    simp [Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  ·
    simp [J, Matrix.mul_apply, Fin.sum_univ_two]
    change ((![-1, 0] : I → k) ⬝ᵥ fun i => (Bᵀ) i 0) = -B 0 0
    simp [Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  ·
    simp [J, Matrix.mul_apply, Fin.sum_univ_two]
    change B 1 0 + ((![-1, 0] : I → k) ⬝ᵥ fun i => (Bᵀ) i 1) = 0
    simp [Matrix.vecMul, dotProduct, Fin.sum_univ_two]

lemma det_eq_one_of_preserves_J {A : Matrix I I k} (hA : Aᵀ * J * A = J) :
    A.det = 1 := by
  rw [transpose_mul_J_mul] at hA
  have h01 := congrArg (fun M => M 0 1) hA
  simpa [J] using h01

lemma fixes_ω12_fromBlocks_iff_raw
    (A B C D : Matrix I I k) :
    Fixes ω12 (Matrix.fromBlocks A B C D) ↔
      Aᵀ * J * A = J ∧
      Aᵀ * J * B = 0 ∧
      Bᵀ * J * A = 0 ∧
      Bᵀ * J * B = 0 := by
  constructor
  · intro h
    refine ⟨?_, ?_, ?_, ?_⟩
    · ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h
      simpa [Fixes, ω12, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h
      simpa [Fixes, ω12, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inl j)) h
      simpa [Fixes, ω12, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h
      simpa [Fixes, ω12, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · rintro ⟨hAA, hAB, hBA, hBB⟩
    ext i j
    cases i <;> cases j <;>
      simp [ω12, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply, hAA, hAB, hBA, hBB]

lemma fixes_ω34_fromBlocks_iff_raw
    (A B C D : Matrix I I k) :
    Fixes ω34 (Matrix.fromBlocks A B C D) ↔
      Cᵀ * J * C = 0 ∧
      Cᵀ * J * D = 0 ∧
      Dᵀ * J * C = 0 ∧
      Dᵀ * J * D = J := by
  constructor
  · intro h
    refine ⟨?_, ?_, ?_, ?_⟩
    · ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h
      simpa [Fixes, ω34, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h
      simpa [Fixes, ω34, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inl j)) h
      simpa [Fixes, ω34, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h
      simpa [Fixes, ω34, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  · rintro ⟨hCC, hCD, hDC, hDD⟩
    ext i j
    cases i <;> cases j <;>
      simp [ω34, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply, hCC, hCD, hDC, hDD]

lemma preserves_J_left_isUnit_det {A : Matrix I I k} (hA : Aᵀ * J * A = J) :
    IsUnit ((Aᵀ * J).det) := by
  have hdet : A.det = 1 := det_eq_one_of_preserves_J hA
  rw [Matrix.det_mul, Matrix.det_transpose, hdet, J_det, one_mul]
  exact isUnit_one

lemma preserves_J_right_isUnit_det {D : Matrix I I k} (hD : Dᵀ * J * D = J) :
    IsUnit ((Dᵀ * J).det) := by
  have hdet : D.det = 1 := det_eq_one_of_preserves_J hD
  rw [Matrix.det_mul, Matrix.det_transpose, hdet, J_det, one_mul]
  exact isUnit_one

lemma fixes_ω12_fromBlocks_iff
    (A B C D : Matrix I I k) :
    Fixes ω12 (Matrix.fromBlocks A B C D) ↔
      A.det = 1 ∧ B = 0 := by
  constructor
  · intro h
    rcases (fixes_ω12_fromBlocks_iff_raw A B C D).1 h with ⟨hAA, hAB, _, _⟩
    have hdet : A.det = 1 := det_eq_one_of_preserves_J hAA
    have hunit : IsUnit ((Aᵀ * J).det) := preserves_J_left_isUnit_det hAA
    have hB : B = 0 := by
      calc
        B = 1 * B := by simp
        _ = ((Aᵀ * J)⁻¹ * (Aᵀ * J)) * B := by
          rw [Matrix.nonsing_inv_mul _ hunit]
        _ = (Aᵀ * J)⁻¹ * ((Aᵀ * J) * B) := by simp [Matrix.mul_assoc]
        _ = 0 := by
          have htmp := congrArg (fun M => (Aᵀ * J)⁻¹ * M) hAB
          have hleft : (Aᵀ * J)⁻¹ * ((Aᵀ * J) * B) = 0 := by
            simpa [Matrix.mul_assoc, Matrix.mul_zero] using htmp
          exact hleft
    exact ⟨hdet, hB⟩
  · rintro ⟨hdet, hB⟩
    have hAA : Aᵀ * J * A = J := by
      rw [transpose_mul_J_mul, hdet]
      simp [J]
    have hAB : Aᵀ * J * B = 0 := by simp [hB]
    have hBA : Bᵀ * J * A = 0 := by simp [hB]
    have hBB : Bᵀ * J * B = 0 := by simp [hB]
    exact (fixes_ω12_fromBlocks_iff_raw A B C D).2 ⟨hAA, hAB, hBA, hBB⟩

lemma fixes_ω34_fromBlocks_iff
    (A B C D : Matrix I I k) :
    Fixes ω34 (Matrix.fromBlocks A B C D) ↔
      C = 0 ∧ D.det = 1 := by
  constructor
  · intro h
    rcases (fixes_ω34_fromBlocks_iff_raw A B C D).1 h with ⟨_, _, hDC, hDD⟩
    have hdet : D.det = 1 := det_eq_one_of_preserves_J hDD
    have hunit : IsUnit ((Dᵀ * J).det) := preserves_J_right_isUnit_det hDD
    have hC : C = 0 := by
      calc
        C = 1 * C := by simp
        _ = ((Dᵀ * J)⁻¹ * (Dᵀ * J)) * C := by
          rw [Matrix.nonsing_inv_mul _ hunit]
        _ = (Dᵀ * J)⁻¹ * ((Dᵀ * J) * C) := by simp [Matrix.mul_assoc]
        _ = 0 := by
          have htmp := congrArg (fun M => (Dᵀ * J)⁻¹ * M) hDC
          have hleft : (Dᵀ * J)⁻¹ * ((Dᵀ * J) * C) = 0 := by
            simpa [Matrix.mul_assoc, Matrix.mul_zero] using htmp
          exact hleft
    exact ⟨hC, hdet⟩
  · rintro ⟨hC, hdet⟩
    have hCC : Cᵀ * J * C = 0 := by simp [hC]
    have hCD : Cᵀ * J * D = 0 := by simp [hC]
    have hDC : Dᵀ * J * C = 0 := by simp [hC]
    have hDD : Dᵀ * J * D = J := by
      rw [transpose_mul_J_mul, hdet]
      simp [J]
    exact (fixes_ω34_fromBlocks_iff_raw A B C D).2 ⟨hCC, hCD, hDC, hDD⟩

lemma fixes_splitPair_iff_fixes_components (g : Matrix V V k) :
    FixesSplitPair g ↔ Fixes ω12 g ∧ Fixes ω34 g := by
  constructor
  · rintro ⟨hsum, h34⟩
    have hsum' : gᵀ * ω12 * g + gᵀ * ω34 * g = ω12 + ω34 := by
      simpa [FixesSplitPair, Fixes, splitRep₁, Matrix.mul_add, Matrix.add_mul,
        Matrix.mul_assoc] using hsum
    have h34eq : gᵀ * ω34 * g = ω34 := by
      simpa [Fixes] using h34
    have h12 : gᵀ * ω12 * g = ω12 := by
      have hsum'' : gᵀ * ω12 * g + ω34 = ω12 + ω34 := by
        simpa [h34eq] using hsum'
      exact add_right_cancel hsum''
    exact ⟨h12, h34⟩
  · rintro ⟨h12, h34⟩
    have h12eq : gᵀ * ω12 * g = ω12 := by
      simpa [Fixes] using h12
    have h34eq : gᵀ * ω34 * g = ω34 := by
      simpa [Fixes] using h34
    have hsum : Fixes splitRep₁ g := by
      rw [Fixes, splitRep₁]
      calc
        gᵀ * (ω12 + ω34) * g = gᵀ * ω12 * g + gᵀ * ω34 * g := by
          simp [Matrix.mul_add, Matrix.add_mul, Matrix.mul_assoc]
        _ = ω12 + ω34 := by rw [h12eq, h34eq]
    exact ⟨hsum, h34⟩

theorem fixes_splitPair_fromBlocks_iff
    (A B C D : Matrix I I k) :
    FixesSplitPair (Matrix.fromBlocks A B C D) ↔
      A.det = 1 ∧ B = 0 ∧ C = 0 ∧ D.det = 1 := by
  rw [fixes_splitPair_iff_fixes_components, fixes_ω12_fromBlocks_iff, fixes_ω34_fromBlocks_iff]
  constructor
  · rintro ⟨hA, hC⟩
    exact ⟨hA.1, hA.2, hC.1, hC.2⟩
  · rintro ⟨hA, hB, hC, hD⟩
    exact ⟨⟨hA, hB⟩, ⟨hC, hD⟩⟩

theorem fixes_splitPairBivector_fromBlocks_iff
    (A B C D : Matrix I I k) :
    FixesSplitPairBivector (Matrix.fromBlocks A B C D) ↔
      A.det = 1 ∧ B = 0 ∧ C = 0 ∧ D.det = 1 := by
  constructor
  · intro h
    have hT : FixesSplitPair (Matrix.fromBlocks Aᵀ Cᵀ Bᵀ Dᵀ) := by
      simpa [FixesSplitPairBivector, FixesSplitPair, FixesBivector, Fixes,
        Matrix.fromBlocks_transpose] using h
    have hT' :=
      (fixes_splitPair_fromBlocks_iff (A := Aᵀ) (B := Cᵀ) (C := Bᵀ) (D := Dᵀ)).1 hT
    simpa [Matrix.det_transpose, transpose_eq_zero, and_assoc, and_left_comm, and_comm] using hT'
  · intro h
    have hT :
        (Aᵀ).det = 1 ∧ Cᵀ = 0 ∧ Bᵀ = 0 ∧ (Dᵀ).det = 1 := by
      simpa [Matrix.det_transpose, transpose_eq_zero, and_assoc, and_left_comm, and_comm] using h
    have hFix :
        FixesSplitPair (Matrix.fromBlocks Aᵀ Cᵀ Bᵀ Dᵀ) :=
      (fixes_splitPair_fromBlocks_iff (A := Aᵀ) (B := Cᵀ) (C := Bᵀ) (D := Dᵀ)).2 hT
    simpa [FixesSplitPairBivector, FixesSplitPair, FixesBivector, Fixes,
      Matrix.fromBlocks_transpose] using hFix

lemma split_diag_act_ω12 (A D : Matrix I I k) :
    ActBivector ω12 (Matrix.fromBlocks A 0 0 D) =
      A.det • ω12 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    have h := congrArg (fun M => M i j) (mul_J_transpose_mul (A := A))
    fin_cases i <;> fin_cases j <;>
      simpa [ActBivector, ω12, Matrix.fromBlocks_multiply,
        Matrix.fromBlocks_transpose, J, Matrix.smul_apply] using h
  · simp [ActBivector, ω12, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  · simp [ActBivector, ω12, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  · simp [ActBivector, ω12, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]

lemma split_diag_act_ω34 (A D : Matrix I I k) :
    ActBivector ω34 (Matrix.fromBlocks A 0 0 D) =
      D.det • ω34 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · simp [ActBivector, ω34, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  · simp [ActBivector, ω34, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  · simp [ActBivector, ω34, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  ·
    have h := congrArg (fun M => M i j) (mul_J_transpose_mul (A := D))
    fin_cases i <;> fin_cases j <;>
      simpa [ActBivector, ω34, Matrix.fromBlocks_multiply,
        Matrix.fromBlocks_transpose, J, Matrix.smul_apply] using h

lemma split_diag_act_rep1 (A D : Matrix I I k) :
    ActBivector splitRep₁ (Matrix.fromBlocks A 0 0 D) =
      A.det • ω12 + D.det • ω34 := by
  calc
    ActBivector splitRep₁ (Matrix.fromBlocks A 0 0 D)
        = ActBivector ω12 (Matrix.fromBlocks A 0 0 D)
            + ActBivector ω34 (Matrix.fromBlocks A 0 0 D) := by
            simp [ActBivector, splitRep₁, Matrix.mul_add, Matrix.add_mul]
    _ = A.det • ω12 + D.det • ω34 := by
          rw [split_diag_act_ω12, split_diag_act_ω34]

lemma split_diag_act_rep2 (A D : Matrix I I k) :
    ActBivector splitRep₂ (Matrix.fromBlocks A 0 0 D) =
      D.det • ω34 := by
  simpa [splitRep₂] using split_diag_act_ω34 (A := A) (D := D)

theorem split_diag_preserves_subspace (A D : Matrix I I k) :
    PreservesSplitSubspaceBivector (Matrix.fromBlocks A 0 0 D) := by
  refine ⟨A.det, D.det - A.det, 0, D.det, ?_, ?_⟩
  ·
    rw [split_diag_act_rep1]
    ext i j
    rcases i with i | i <;> rcases j with j | j
    · fin_cases i <;> fin_cases j <;>
        simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg, J]
    · simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg]
    · simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg]
    · fin_cases i <;> fin_cases j <;>
        simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg, J] <;> ring
  ·
    rw [split_diag_act_rep2]
    change D.det • ω34 = 0 • splitRep₁ + D.det • splitRep₂
    simp [splitRep₁, splitRep₂]

lemma split_antidiag_act_ω12 (B C : Matrix I I k) :
    ActBivector ω12 (splitAntiDiag (k := k) B C) =
      C.det • ω34 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · simp [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose, Matrix.smul_apply]
  · simp [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose, Matrix.smul_apply]
  · simp [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose, Matrix.smul_apply]
  ·
    have h := congrArg (fun M => M i j) (mul_J_transpose_mul (A := C))
    fin_cases i <;> fin_cases j <;>
      simpa [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
        Matrix.fromBlocks_transpose, J, Matrix.smul_apply] using h

lemma split_antidiag_act_ω34 (B C : Matrix I I k) :
    ActBivector ω34 (splitAntiDiag (k := k) B C) =
      B.det • ω12 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    have h := congrArg (fun M => M i j) (mul_J_transpose_mul (A := B))
    fin_cases i <;> fin_cases j <;>
      simpa [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
        Matrix.fromBlocks_transpose, J, Matrix.smul_apply] using h
  · simp [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose, Matrix.smul_apply]
  · simp [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose, Matrix.smul_apply]
  · simp [ActBivector, splitAntiDiag, ω12, ω34, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose, Matrix.smul_apply]

theorem split_antidiag_preserves_subspace (B C : Matrix I I k) :
    PreservesSplitSubspaceBivector (splitAntiDiag (k := k) B C) := by
  refine ⟨B.det, C.det - B.det, B.det, -B.det, ?_, ?_⟩
  ·
    calc
      ActBivector splitRep₁ (splitAntiDiag (k := k) B C)
          = ActBivector ω12 (splitAntiDiag (k := k) B C)
              + ActBivector ω34 (splitAntiDiag (k := k) B C) := by
                simp [ActBivector, splitRep₁, Matrix.mul_add, Matrix.add_mul]
      _ = C.det • ω34 + B.det • ω12 := by
            rw [split_antidiag_act_ω12, split_antidiag_act_ω34]
      _ = B.det • splitRep₁ + (C.det - B.det) • splitRep₂ := by
            ext i j
            rcases i with i | i <;> rcases j with j | j
            · fin_cases i <;> fin_cases j <;>
                simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg, J]
            · simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg]
            · simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg]
            · fin_cases i <;> fin_cases j <;>
                simp [splitRep₁, splitRep₂, ω12, ω34, smul_add, sub_eq_add_neg, J] <;> ring
  ·
    calc
      ActBivector splitRep₂ (splitAntiDiag (k := k) B C)
          = ActBivector ω34 (splitAntiDiag (k := k) B C) := by
                simp [ActBivector, splitRep₂]
      _ = B.det • ω12 := by rw [split_antidiag_act_ω34]
      _ = B.det • splitRep₁ + (-B.det) • splitRep₂ := by
            simp [splitRep₁, splitRep₂, smul_add, sub_eq_add_neg]

lemma det_split_linearCombination (a b : k) :
    Matrix.det (a • ω12 + b • ω34) = (a * b) * (a * b) := by
  have hblock :
      a • ω12 + b • ω34 =
        Matrix.fromBlocks (a • J) (0 : Matrix I I k) (0 : Matrix I I k) (b • J) := by
    ext i j
    rcases i with i | i <;> rcases j with j | j <;>
      simp [ω12, ω34, Matrix.fromBlocks, Matrix.smul_apply, add_comm, add_left_comm, add_assoc]
  rw [hblock]
  rw [Matrix.det_fromBlocks_zero₂₁]
  rw [Matrix.det_smul, Matrix.det_smul]
  have hcard : Fintype.card I = 2 := by simp [I]
  rw [hcard, J_det]
  ring

lemma splitSwap_transpose :
    (splitSwap (k := k))ᵀ = splitSwap (k := k) := by
  simpa [splitSwap, Matrix.fromBlocks_transpose]

lemma splitSwap_act_ω12 :
    ActBivector ω12 (splitSwap (k := k)) = ω34 := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, splitSwap, ω12, ω34, splitSwap_transpose,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]

lemma splitSwap_act_ω34 :
    ActBivector ω34 (splitSwap (k := k)) = ω12 := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, splitSwap, ω12, ω34, splitSwap_transpose,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]

theorem split_swap_preserves_subspace :
    PreservesSplitSubspaceBivector (splitSwap (k := k)) := by
  have hcoeff1 :
      (1 : k) • (splitRep₁ (k := k)) + (0 : k) • (splitRep₂ (k := k)) = splitRep₁ (k := k) := by
    simp
  have hcoeff2 :
      (1 : k) • (splitRep₁ (k := k)) + (-1 : k) • (splitRep₂ (k := k)) =
        splitRep₁ (k := k) - splitRep₂ (k := k) := by
    simp [sub_eq_add_neg]
  refine ⟨1, 0, 1, -1, ?_, ?_⟩
  ·
    rw [hcoeff1]
    calc
      ActBivector splitRep₁ (splitSwap (k := k))
          = ActBivector ω12 (splitSwap (k := k))
              + ActBivector ω34 (splitSwap (k := k)) := by
                simp [ActBivector, splitRep₁, Matrix.mul_add, Matrix.add_mul]
      _ = ω34 + ω12 := by rw [splitSwap_act_ω12, splitSwap_act_ω34]
      _ = splitRep₁ := by simp [splitRep₁, add_comm]
  ·
    rw [hcoeff2]
    calc
      ActBivector splitRep₂ (splitSwap (k := k))
          = ActBivector ω34 (splitSwap (k := k)) := by
                simp [ActBivector, splitRep₂]
      _ = ω12 := by rw [splitSwap_act_ω34]
      _ = splitRep₁ - splitRep₂ := by simp [splitRep₁, splitRep₂]

lemma onePoint_upperShear_act_rep1 (B : Matrix I I k) :
    ActBivector onePointRep₁ (onePointUpperShear (k := k) B) =
      onePointRep₁ + (B 0 0 + B 1 1) • onePointRep₂ := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    have h := congrArg (fun M => M i j) (mul_J_add_J_mul_transpose (B := B))
    fin_cases i <;> fin_cases j <;>
      simpa [ActBivector, onePointRep₁, onePointRep₂, onePointUpperShear, ω12,
        Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose, J, Matrix.smul_apply,
        Matrix.add_apply] using h
  ·
    simp [ActBivector, onePointRep₁, onePointRep₂, onePointUpperShear, ω12,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose, Matrix.add_apply]
  ·
    simp [ActBivector, onePointRep₁, onePointRep₂, onePointUpperShear, ω12,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose, Matrix.add_apply]
  ·
    simp [ActBivector, onePointRep₁, onePointRep₂, onePointUpperShear, ω12,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose, Matrix.add_apply]

lemma onePoint_upperShear_act_rep2 (B : Matrix I I k) :
    ActBivector onePointRep₂ (onePointUpperShear (k := k) B) = onePointRep₂ := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · simp [ActBivector, onePointRep₂, onePointUpperShear, ω12,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  · simp [ActBivector, onePointRep₂, onePointUpperShear, ω12,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  · simp [ActBivector, onePointRep₂, onePointUpperShear, ω12,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]
  · simp [ActBivector, onePointRep₂, onePointUpperShear, ω12,
      Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose]

theorem onePoint_upperShear_preserves_subspace (B : Matrix I I k) :
    PreservesOnePointSubspaceBivector (onePointUpperShear (k := k) B) := by
  refine ⟨1, B 0 0 + B 1 1, 0, 1, ?_, ?_⟩
  · rw [onePoint_upperShear_act_rep1]
    simp
  · rw [onePoint_upperShear_act_rep2]
    simp

lemma onePoint_scale_act_rep1 (a : k) :
    ActBivector onePointRep₁ (onePointScale (k := k) a) = a • onePointRep₁ := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · simp [ActBivector, onePointRep₁, onePointScale, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose]
  · fin_cases i <;> fin_cases j <;>
      simp [ActBivector, onePointRep₁, onePointScale, Matrix.fromBlocks_multiply,
        Matrix.fromBlocks_transpose, J, Matrix.smul_apply]
  · fin_cases i <;> fin_cases j <;>
      simp [ActBivector, onePointRep₁, onePointScale, Matrix.fromBlocks_multiply,
        Matrix.fromBlocks_transpose, J, Matrix.smul_apply]
  · simp [ActBivector, onePointRep₁, onePointScale, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose]

lemma onePoint_scale_act_rep2 (a : k) :
    ActBivector onePointRep₂ (onePointScale (k := k) a) = a ^ 2 • onePointRep₂ := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · fin_cases i <;> fin_cases j <;>
      simp [ActBivector, onePointRep₂, onePointScale, ω12, Matrix.fromBlocks_multiply,
        Matrix.fromBlocks_transpose, J, Matrix.smul_apply, pow_two] <;> ring
  · simp [ActBivector, onePointRep₂, onePointScale, ω12, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose]
  · simp [ActBivector, onePointRep₂, onePointScale, ω12, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose]
  · simp [ActBivector, onePointRep₂, onePointScale, ω12, Matrix.fromBlocks_multiply,
      Matrix.fromBlocks_transpose]

theorem onePoint_scale_preserves_subspace (a : k) :
    PreservesOnePointSubspaceBivector (onePointScale (k := k) a) := by
  refine ⟨a, 0, 0, a ^ 2, ?_, ?_⟩
  · rw [onePoint_scale_act_rep1]
    simp
  · rw [onePoint_scale_act_rep2]
    simp


lemma fixes_onePointRep₁_fromBlocks_iff_raw
    (A B C D : Matrix I I k) :
    Fixes onePointRep₁ (Matrix.fromBlocks A B C D) ↔
      Aᵀ * J * C + Cᵀ * J * A = 0 ∧
      Aᵀ * J * D + Cᵀ * J * B = J ∧
      Bᵀ * J * C + Dᵀ * J * A = J ∧
      Bᵀ * J * D + Dᵀ * J * B = 0 := by
  constructor
  · intro h
    refine ⟨?_, ?_, ?_, ?_⟩
    · ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h
      simpa [Fixes, onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h
      simpa [Fixes, onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inl j)) h
      simpa [Fixes, onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij
    · ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h
      simpa [Fixes, onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij
  · rintro ⟨hAC, hAD, hDA, hDB⟩
    ext i j
    rcases i with i | i <;> rcases j with j | j
    ·
      have hij := congrFun (congrFun hAC i) j
      simpa [onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij
    ·
      have hij := congrFun (congrFun hAD i) j
      simpa [onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij
    ·
      have hij := congrFun (congrFun hDA i) j
      simpa [onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij
    ·
      have hij := congrFun (congrFun hDB i) j
      simpa [onePointRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.add_apply, add_comm] using hij

theorem fixes_onePointPair_fromBlocks_iff
    (A B C D : Matrix I I k) :
    FixesOnePointPair (Matrix.fromBlocks A B C D) ↔
      A.det = 1 ∧ B = 0 ∧ D = A ∧ Aᵀ * J * C + Cᵀ * J * A = 0 := by
  constructor
  · intro h
    rcases h with ⟨h1, h2⟩
    rcases (fixes_ω12_fromBlocks_iff A B C D).1 h2 with ⟨hdet, hB⟩
    rcases (fixes_onePointRep₁_fromBlocks_iff_raw A B C D).1 h1 with
      ⟨hAC, hAD, hDA, hDB⟩
    have hAA : Aᵀ * J * A = J := by
      rw [transpose_mul_J_mul, hdet]
      simp [J]
    have hunit : IsUnit ((Aᵀ * J).det) := preserves_J_left_isUnit_det (by
      simpa [hAA])
    have hD : D = A := by
      have hEq : Aᵀ * J * D = Aᵀ * J * A := by
        simpa [hB, hAA, Matrix.zero_mul, zero_add] using hAD
      calc
        D = 1 * D := by simp
        _ = ((Aᵀ * J)⁻¹ * (Aᵀ * J)) * D := by
          rw [Matrix.nonsing_inv_mul _ hunit]
        _ = (Aᵀ * J)⁻¹ * ((Aᵀ * J) * D) := by simp [Matrix.mul_assoc]
        _ = (Aᵀ * J)⁻¹ * ((Aᵀ * J) * A) := by simpa [Matrix.mul_assoc] using congrArg (fun M => (Aᵀ * J)⁻¹ * M) hEq
        _ = A := by
          rw [← Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hunit]
          simp
    exact ⟨hdet, hB, hD, hAC⟩
  · rintro ⟨hdet, hB, hD, hAC⟩
    have hAA : Aᵀ * J * A = J := by
      rw [transpose_mul_J_mul, hdet]
      simp [J]
    have hAD : Aᵀ * J * D + Cᵀ * J * B = J := by
      simp [hD, hAA, hB]
    have hDA : Bᵀ * J * C + Dᵀ * J * A = J := by
      simp [hD, hAA, hB]
    have hDB : Bᵀ * J * D + Dᵀ * J * B = 0 := by
      simp [hD, hB]
    exact ⟨(fixes_onePointRep₁_fromBlocks_iff_raw A B C D).2 ⟨hAC, hAD, hDA, hDB⟩,
      (fixes_ω12_fromBlocks_iff A B C D).2 ⟨hdet, hB⟩⟩

theorem fixes_onePointPairBivector_fromBlocks_iff
    (A B C D : Matrix I I k) :
    FixesOnePointPairBivector (Matrix.fromBlocks A B C D) ↔
      A.det = 1 ∧ C = 0 ∧ D = A ∧ A * J * Bᵀ + B * J * Aᵀ = 0 := by
  constructor
  · intro h
    have hT : FixesOnePointPair (Matrix.fromBlocks Aᵀ Cᵀ Bᵀ Dᵀ) := by
      simpa [FixesOnePointPairBivector, FixesOnePointPair, FixesBivector, Fixes,
        Matrix.fromBlocks_transpose] using h
    have hT' :=
      (fixes_onePointPair_fromBlocks_iff (A := Aᵀ) (B := Cᵀ) (C := Bᵀ) (D := Dᵀ)).1 hT
    simpa [Matrix.det_transpose, transpose_eq_zero, Matrix.transpose_mul] using hT'
  · intro h
    have hT :
        (Aᵀ).det = 1 ∧ Cᵀ = 0 ∧ Dᵀ = Aᵀ ∧ A * J * Bᵀ + B * J * Aᵀ = 0 := by
      simpa [Matrix.det_transpose, transpose_eq_zero, Matrix.transpose_mul] using h
    have hFix :
        FixesOnePointPair (Matrix.fromBlocks Aᵀ Cᵀ Bᵀ Dᵀ) :=
      (fixes_onePointPair_fromBlocks_iff (A := Aᵀ) (B := Cᵀ) (C := Bᵀ) (D := Dᵀ)).2 hT
    simpa [FixesOnePointPairBivector, FixesOnePointPair, FixesBivector, Fixes,
      Matrix.fromBlocks_transpose] using hFix


end N4
end Wedge2Formalization
