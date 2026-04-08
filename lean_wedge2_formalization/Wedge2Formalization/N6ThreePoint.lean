import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

open Matrix

namespace Wedge2Formalization
namespace N6ThreePoint

variable {k : Type*} [Field k]

/-- The `4`-dimensional block carrying two simple support points. -/
abbrev W := N4.V

/-- The `2`-dimensional block carrying the third simple support point. -/
abbrev I := N4.I

/-- The ambient `6`-dimensional space. -/
abbrev V := W ⊕ I

/-- Exact stabilizer of a bivector under the natural `GL(V)` action on `∧²V`. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The bivector action of a product factors through successive actions. -/
theorem actBivector_mul
    (Ω : Matrix V V k) (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  rw [ActBivector, ActBivector, ActBivector, Matrix.transpose_mul]
  repeat rw [Matrix.mul_assoc]

/-- The bivector action is additive in the bivector variable. -/
theorem actBivector_add
    (Ω₁ Ω₂ : Matrix V V k) (g : Matrix V V k) :
    ActBivector (Ω₁ + Ω₂) g = ActBivector Ω₁ g + ActBivector Ω₂ g := by
  simp [ActBivector, Matrix.add_mul, Matrix.mul_add]

/-- The bivector action is `k`-linear in the bivector variable. -/
theorem actBivector_smul
    (a : k) (Ω : Matrix V V k) (g : Matrix V V k) :
    ActBivector (a • Ω) g = a • ActBivector Ω g := by
  simp [ActBivector, Matrix.mul_assoc]

/-- Acting by a block diagonal matrix keeps the representative block diagonal. -/
theorem act_blockDiagonal
    (ΩW : Matrix W W k)
    (ΩI : Matrix I I k)
    (H : Matrix W W k)
    (E : Matrix I I k) :
    ActBivector (Matrix.fromBlocks ΩW 0 0 ΩI) (Matrix.fromBlocks H 0 0 E) =
      Matrix.fromBlocks (N4.ActBivector ΩW H) 0 0 (E * ΩI * Eᵀ) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N4.ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

/-- A representative for the three-point orbit. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks (N4.splitRep₁ (k := k)) 0 0 0

/-- The second basis vector of the three-point representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks (N4.splitRep₂ (k := k)) 0 0 (N4.J (k := k))

/-- Exact stabilizer of the three-point pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rep₁ g ∧ FixesBivector rep₂ g

/-- The embedded split-support block is nonsingular. -/
lemma det_splitRep₁ : Matrix.det (N4.splitRep₁ (k := k)) = 1 := by
  simpa [N4.splitRep₁] using
    N4.det_split_linearCombination (k := k) (a := (1 : k)) (b := (1 : k))

/-- In the chosen basis, every vector on the three-point line is block diagonal. -/
theorem basis_linearCombination
    (a b : k) :
    a • (rep₁ (k := k)) + b • (rep₂ (k := k)) =
      Matrix.fromBlocks
        (a • (N4.splitRep₁ (k := k)) + b • (N4.splitRep₂ (k := k)))
        0
        0
        (b • (N4.J (k := k))) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₁, rep₂, Matrix.add_apply]

/-- The determinant on the three-point line factors as the square of
`a(a+b)b`. -/
theorem det_in_basis
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) =
      (a * (a + b) * b) * (a * (a + b) * b) := by
  rw [basis_linearCombination (k := k) (a := a) (b := b)]
  rw [Matrix.det_fromBlocks_zero₂₁]
  have htop :=
    N4Summary.split_det_in_basis (k := k) (a := a) (b := b)
  have hcard : Fintype.card N4.I = 2 := by
    simp [N4.I]
  have hbot : Matrix.det (b • (N4.J (k := k))) = b * b := by
    simp [Matrix.det_smul, hcard, N4.J_det, pow_two]
  rw [htop, hbot]
  ring

/-- Rank drop on the three-point line occurs exactly on the three distinguished
projective points. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 ∨ b = 0 := by
  rw [det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    have hsq : (a * (a + b) * b) = 0 := by
      rcases mul_eq_zero.mp h with h0 | h0 <;> exact h0
    have hprod : (a * (a + b)) * b = 0 := by
      simpa [mul_assoc] using hsq
    rcases mul_eq_zero.mp hprod with hab | hb
    · rcases mul_eq_zero.mp hab with ha | hab0
      · exact Or.inl ha
      · exact Or.inr (Or.inl hab0)
    · exact Or.inr (Or.inr hb)
  · rintro (ha | hab | hb)
    · simp [ha]
    · simp [hab]
    · simp [hb]

/-- Pointwise stabilizer theorem for the three simple support orbit. -/
theorem fixesPair_fromBlocks_iff
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ C = 0 ∧ N4.FixesSplitPairBivector A ∧ D.det = 1 := by
  constructor
  · rintro ⟨h1, h2⟩
    have hA1 : N4.FixesBivector (N4.splitRep₁ (k := k)) A := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h1
      simpa [FixesBivector, rep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        N4.FixesBivector, N4.ActBivector] using hij
    have hAunit : IsUnit (A.det) := by
      have hdetEq := congrArg Matrix.det hA1
      have hsq : A.det * A.det = 1 := by
        simpa [N4.FixesBivector, N4.ActBivector, Matrix.det_mul, Matrix.det_transpose,
          det_splitRep₁] using hdetEq
      have hne : A.det ≠ 0 := by
        intro hzero
        rw [hzero, zero_mul] at hsq
        exact zero_ne_one hsq
      exact isUnit_iff_ne_zero.mpr hne
    have hCA : C * (N4.splitRep₁ (k := k) * Aᵀ) = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inl j)) h1
      simpa [FixesBivector, rep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_assoc] using hij
    have hunitTop : IsUnit ((N4.splitRep₁ (k := k) * Aᵀ).det) := by
      rw [Matrix.det_mul, Matrix.det_transpose, det_splitRep₁]
      simpa using hAunit
    have hC : C = 0 := by
      calc
        C = C * (1 : Matrix W W k) := by simp
        _ = C * ((N4.splitRep₁ (k := k) * Aᵀ) * (N4.splitRep₁ (k := k) * Aᵀ)⁻¹) := by
              rw [Matrix.mul_nonsing_inv _ hunitTop]
        _ = (C * (N4.splitRep₁ (k := k) * Aᵀ)) * (N4.splitRep₁ (k := k) * Aᵀ)⁻¹ := by
              simp [Matrix.mul_assoc]
        _ = 0 := by simpa [hCA]
    have hDD : D * N4.J * Dᵀ = N4.J := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h2
      simpa [FixesBivector, rep₂, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have hDdet : D.det = 1 := by
      rw [N4.mul_J_transpose_mul] at hDD
      have h01 := congrArg (fun M => M 0 1) hDD
      simpa [N4.J] using h01
    have hunitBotT : IsUnit ((Dᵀ).det) := by
      simpa [Matrix.det_transpose, hDdet] using (isUnit_one : IsUnit (1 : k))
    have hBJ : B * N4.J * Dᵀ = (0 : Matrix W I k) := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h2
      simpa [FixesBivector, rep₂, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have hBJ0 : B * N4.J = (0 : Matrix W I k) := by
      have htmp := congrArg (fun M => M * (Dᵀ)⁻¹) hBJ
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitBotT] using htmp
    have hJunit : IsUnit ((N4.J (k := k)).det) := by
      simpa [N4.J_det] using (isUnit_one : IsUnit (1 : k))
    have hB : B = 0 := by
      calc
        B = B * (1 : Matrix I I k) := by simp
        _ = B * (N4.J * N4.J⁻¹) := by rw [Matrix.mul_nonsing_inv _ hJunit]
        _ = (B * N4.J) * N4.J⁻¹ := by simp [Matrix.mul_assoc]
        _ = 0 := by simpa [hBJ0]
    have hA2 : N4.FixesBivector (N4.splitRep₂ (k := k)) A := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h2
      simpa [FixesBivector, rep₂, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        N4.FixesBivector, N4.ActBivector] using hij
    exact ⟨hB, hC, ⟨hA1, hA2⟩, hDdet⟩
  · rintro ⟨hB, hC, hA, hDdet⟩
    refine ⟨?_, ?_⟩
    · ext i j
      cases i with
      | inl i =>
          cases j with
          | inl j =>
              have hij := congrArg (fun M => M i j) hA.1
              simpa [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
                N4.FixesBivector, N4.ActBivector] using hij
          | inr j =>
              simp [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
      | inr i =>
          cases j with
          | inl j =>
              simp [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
          | inr j =>
              simp [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
    ·
      have hDD : D * N4.J * Dᵀ = N4.J := by
        rw [N4.mul_J_transpose_mul, hDdet]
        simp [N4.J]
      ext i j
      cases i with
      | inl i =>
          cases j with
          | inl j =>
              have hij := congrArg (fun M => M i j) hA.2
              simpa [FixesBivector, rep₂, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
                N4.FixesBivector, N4.ActBivector] using hij
          | inr j =>
              simp [FixesBivector, rep₂, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
      | inr i =>
          cases j with
          | inl j =>
              simp [FixesBivector, rep₂, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
          | inr j =>
              have hij := congrArg (fun M => M i j) hDD
              simpa [FixesBivector, rep₂, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- The obvious block-diagonal pointwise family on the three-point orbit. -/
theorem pointwise_family
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    FixesPairBivector (Matrix.fromBlocks A 0 0 D) := by
  exact (fixesPair_fromBlocks_iff (k := k) (A := A) (B := 0) (C := 0) (D := D)).2
    ⟨rfl, rfl, hA, hD⟩

/-- The full block-diagonal pointwise family on the three-point orbit has determinant `1`. -/
theorem pointwise_family_det
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    Matrix.det (Matrix.fromBlocks A 0 0 D) = 1 := by
  have hA' :
      N4.FixesSplitPairBivector
        (Matrix.fromBlocks A.toBlocks₁₁ A.toBlocks₁₂ A.toBlocks₂₁ A.toBlocks₂₂) := by
    simpa [Matrix.fromBlocks_toBlocks] using hA
  rcases
      (N4Summary.split_pointwise_bivector_iff
        (k := k)
        (A := A.toBlocks₁₁)
        (B := A.toBlocks₁₂)
        (C := A.toBlocks₂₁)
        (D := A.toBlocks₂₂)).1 hA' with
    ⟨h11, h12, h21, h22⟩
  rw [← Matrix.fromBlocks_toBlocks A, h12, h21, Matrix.det_fromBlocks_zero₂₁,
    Matrix.det_fromBlocks_zero₂₁, h11, h22, hD]
  ring

/-- A sign flip on the simple block has determinant `-1`. -/
lemma signFlip_det :
    Matrix.det (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k) = (-1 : k) := by
  simp [Matrix.det_fin_two]

/-- The sign flip on the simple block is an involution. -/
lemma signFlip_sq :
    (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k) *
        (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k) = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Matrix.mul_apply, Fin.sum_univ_two]

/-- A split swap on the `4`-block together with a sign flip on the simple block
preserves the three-point `2`-space with the expected action. -/
theorem swap_lift_action :
    let E : Matrix I I k := !![(1 : k), 0; 0, (-1 : k)]
    let g : Matrix V V k := Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    ActBivector rep₁ g = rep₁ ∧
      ActBivector rep₂ g = rep₁ - rep₂ := by
  dsimp
  have hswap := N4Summary.split_swap_action (k := k)
  have hbot :
      (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k) * N4.J *
          (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)ᵀ =
        (-1 : k) • N4.J := by
    rw [N4.mul_J_transpose_mul, signFlip_det (k := k)]
    simp [N4.J]
  constructor
  · calc
      ActBivector rep₁
          (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
            (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) =
        Matrix.fromBlocks (N4.ActBivector (N4.splitRep₁ (k := k)) (N4.splitSwap (k := k))) 0 0 0 := by
          simpa [rep₁] using
            act_blockDiagonal
              (k := k)
              (ΩW := N4.splitRep₁ (k := k))
              (ΩI := (0 : Matrix I I k))
              (H := N4.splitSwap (k := k))
              (E := (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k))
      _ = Matrix.fromBlocks (N4.splitRep₁ (k := k)) 0 0 0 := by
            rw [hswap.1]
      _ = rep₁ := by simp [rep₁]
  · calc
      ActBivector rep₂
          (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
            (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) =
        Matrix.fromBlocks
          (N4.ActBivector (N4.splitRep₂ (k := k)) (N4.splitSwap (k := k))) 0 0
          (((!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k) * N4.J *
            (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)ᵀ)) := by
          simpa [rep₂] using
            act_blockDiagonal
              (k := k)
              (ΩW := N4.splitRep₂ (k := k))
              (ΩI := N4.J (k := k))
              (H := N4.splitSwap (k := k))
              (E := (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k))
      _ = Matrix.fromBlocks (N4.splitRep₁ (k := k) - N4.splitRep₂ (k := k)) 0 0 ((-1 : k) • N4.J) := by
            rw [hswap.2, hbot]
      _ = rep₁ - rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.sub_apply, sub_eq_add_neg]

/-- The concrete three-point swap lift is an involution. -/
theorem swap_lift_sq :
    let E : Matrix I I k := !![(1 : k), 0; 0, (-1 : k)]
    let g : Matrix V V k := Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    g * g = 1 := by
  dsimp
  rw [Matrix.fromBlocks_multiply, N4Summary.split_swap_sq, signFlip_sq (k := k)]
  simp

/-- Left-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change the induced quotient action. -/
theorem pointwise_swap_product_lift_action
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix I I k := !![(1 : k), 0; 0, (-1 : k)]
    let gU : Matrix V V k := Matrix.fromBlocks A 0 0 D
    let gL : Matrix V V k := Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    ActBivector rep₁ (gU * gL) = rep₁ ∧
      ActBivector rep₂ (gU * gL) = rep₁ - rep₂ := by
  dsimp
  have hU := pointwise_family (k := k) (A := A) (D := D) hA hD
  have hU1 : ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) = rep₁ := by
    simpa [FixesBivector, ActBivector] using hU.1
  have hU2 : ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) = rep₂ := by
    simpa [FixesBivector, ActBivector] using hU.2
  have hL := swap_lift_action (k := k)
  constructor
  · calc
      ActBivector rep₁
          ((Matrix.fromBlocks A 0 0 D) *
            Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) =
        ActBivector
          (ActBivector rep₁
            (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)))
          (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_mul]
      _ = ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) := by rw [hL.1]
      _ = rep₁ := hU1
  · calc
      ActBivector rep₂
          ((Matrix.fromBlocks A 0 0 D) *
            Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) =
        ActBivector
          (ActBivector rep₂
            (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)))
          (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_mul]
      _ = ActBivector (rep₁ - rep₂) (Matrix.fromBlocks A 0 0 D) := by rw [hL.2]
      _ = ActBivector (rep₁ + (-1 : k) • rep₂) (Matrix.fromBlocks A 0 0 D) := by
            simp [sub_eq_add_neg]
      _ = ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) +
            ActBivector ((-1 : k) • rep₂) (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_add]
      _ = ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) +
            (-1 : k) • ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_smul]
      _ = rep₁ + (-1 : k) • rep₂ := by simp [hU1, hU2]
      _ = rep₁ - rep₂ := by simp [sub_eq_add_neg]

/-- Right-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change the induced quotient action. -/
theorem swap_pointwise_right_product_lift_action
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix I I k := !![(1 : k), 0; 0, (-1 : k)]
    let gL : Matrix V V k := Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    let gU : Matrix V V k := Matrix.fromBlocks A 0 0 D
    ActBivector rep₁ (gL * gU) = rep₁ ∧
      ActBivector rep₂ (gL * gU) = rep₁ - rep₂ := by
  dsimp
  have hU := pointwise_family (k := k) (A := A) (D := D) hA hD
  have hU1 : ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) = rep₁ := by
    simpa [FixesBivector, ActBivector] using hU.1
  have hU2 : ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) = rep₂ := by
    simpa [FixesBivector, ActBivector] using hU.2
  have hL := swap_lift_action (k := k)
  constructor
  · calc
      ActBivector rep₁
          ((Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) *
            Matrix.fromBlocks A 0 0 D) =
        ActBivector
          (ActBivector rep₁ (Matrix.fromBlocks A 0 0 D))
          (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
            (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) := by
            rw [actBivector_mul]
      _ = ActBivector rep₁
            (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) := by rw [hU1]
      _ = rep₁ := hL.1
  · calc
      ActBivector rep₂
          ((Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) *
            Matrix.fromBlocks A 0 0 D) =
        ActBivector
          (ActBivector rep₂ (Matrix.fromBlocks A 0 0 D))
          (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
            (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) := by
            rw [actBivector_mul]
      _ = ActBivector rep₂
            (Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0
              (!![(1 : k), 0; 0, (-1 : k)] : Matrix I I k)) := by rw [hU2]
      _ = rep₁ - rep₂ := hL.2

/-- Left-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change its determinant. -/
theorem pointwise_swap_product_lift_det
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix I I k := !![(1 : k), 0; 0, (-1 : k)]
    let gU : Matrix V V k := Matrix.fromBlocks A 0 0 D
    let gL : Matrix V V k := Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    Matrix.det (gU * gL) = Matrix.det gL := by
  dsimp
  have hUdet := pointwise_family_det (k := k) (A := A) (D := D) hA hD
  rw [Matrix.det_mul, hUdet, one_mul]

/-- Right-multiplying the three-point swap lift by the full block-diagonal pointwise
subgroup does not change its determinant. -/
theorem swap_pointwise_right_product_lift_det
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    let E : Matrix I I k := !![(1 : k), 0; 0, (-1 : k)]
    let gL : Matrix V V k := Matrix.fromBlocks (N4.splitSwap (k := k)) 0 0 E
    let gU : Matrix V V k := Matrix.fromBlocks A 0 0 D
    Matrix.det (gL * gU) = Matrix.det gL := by
  dsimp
  have hUdet := pointwise_family_det (k := k) (A := A) (D := D) hA hD
  rw [Matrix.det_mul, hUdet, mul_one]

end N6ThreePoint
end Wedge2Formalization
