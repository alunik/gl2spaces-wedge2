import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse

open Matrix

namespace Wedge2Formalization
namespace N6WeightedTwoPoint

variable {k : Type*} [Field k]

/-- The `4`-dimensional block carrying the rank-`4` line. -/
abbrev W := N4.V

/-- The `2`-dimensional block carrying the simple line. -/
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

/-- The weighted two-point representative with divisor `2[a]+[b]`. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks (N4.splitRep₁ (k := k)) 0 0 (N4.J (k := k))

/-- The second basis vector of the weighted two-point representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N4.J (k := k))

/-- Exact stabilizer of the weighted two-point pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rep₁ g ∧ FixesBivector rep₂ g

/-- Setwise stabilizer of the weighted two-point `2`-space in the chosen basis. -/
def PreservesSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector rep₁ g = α • rep₁ + β • rep₂ ∧
    ActBivector rep₂ g = γ • rep₁ + δ • rep₂

/-- The first basis vector minus the second is the embedded rank-`4` line. -/
def topRep : Matrix V V k :=
  Matrix.fromBlocks (N4.splitRep₁ (k := k)) 0 0 0

lemma rep₁_sub_rep₂ : rep₁ (k := k) - rep₂ (k := k) = topRep (k := k) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₁, rep₂, topRep, Matrix.fromBlocks, Matrix.sub_apply]

/-- The rank-`4` block is nonsingular. -/
lemma det_splitRep₁ : Matrix.det (N4.splitRep₁ (k := k)) = 1 := by
  simpa [N4.splitRep₁] using N4.det_split_linearCombination (k := k) (a := (1 : k)) (b := (1 : k))

/-- In the chosen basis, every vector on the weighted two-point line is block diagonal. -/
theorem basis_linearCombination
    (a b : k) :
    a • (rep₁ (k := k)) + b • (rep₂ (k := k)) =
      Matrix.fromBlocks
        (a • (N4.splitRep₁ (k := k)))
        0
        0
        ((a + b) • (N4.J (k := k))) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₁, rep₂, Matrix.add_apply, smul_add] <;> ring

/-- The determinant on the weighted two-point line factors as the square of
`a^2 (a+b)`. -/
theorem det_in_basis
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) =
      (a * a * (a + b)) * (a * a * (a + b)) := by
  rw [basis_linearCombination (k := k) (a := a) (b := b)]
  rw [Matrix.det_fromBlocks_zero₂₁]
  have hcardW : Fintype.card W = 4 := by
    simp [W, N4.V, N4.I]
  have hcardI : Fintype.card I = 2 := by
    simp [I, N4.I]
  have htop : Matrix.det (a • (N4.splitRep₁ (k := k))) = (a * a) * (a * a) := by
    simp [Matrix.det_smul, hcardW, det_splitRep₁, pow_two]
    ring
  have hbot : Matrix.det ((a + b) • (N4.J (k := k))) = (a + b) * (a + b) := by
    simp [Matrix.det_smul, hcardI, N4.J_det, pow_two]
  rw [htop, hbot]
  ring

/-- Rank drop on the weighted two-point line occurs exactly on the two
distinguished projective points. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 := by
  rw [det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    have hsq : a * a * (a + b) = 0 := by
      rcases mul_eq_zero.mp h with h0 | h0 <;> exact h0
    rcases mul_eq_zero.mp hsq with haa | hab
    · rcases mul_eq_zero.mp haa with ha | ha
      · exact Or.inl ha
      · exact Or.inl ha
    · exact Or.inr hab
  · rintro (ha | hab)
    · simp [ha]
    · simp [hab]

/-- Fixing the lower-right simple line forces the upper-right block to vanish and the
lower-right block to have determinant `1`. -/
theorem fixes_rep₂_fromBlocks_iff
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    FixesBivector rep₂ (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ D.det = 1 := by
  constructor
  · intro h
    have h22 : D * N4.J * Dᵀ = N4.J := by
      ext i j
      have h' := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h
      simpa [FixesBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using h'
    have hdet : D.det = 1 := by
      rw [N4.mul_J_transpose_mul] at h22
      have h01 := congrArg (fun M => M 0 1) h22
      simpa [N4.J] using h01
    have hunit : IsUnit ((D * N4.J).det) := by
      rw [Matrix.det_mul, hdet, N4.J_det, one_mul]
      exact isUnit_one
    have h21 : D * N4.J * Bᵀ = 0 := by
      ext i j
      have h' := congrArg (fun M => M (Sum.inr i) (Sum.inl j)) h
      simpa [FixesBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using h'
    have hBt : Bᵀ = 0 := by
      calc
        Bᵀ = (1 : Matrix I I k) * Bᵀ := by simp
        _ = ((D * N4.J)⁻¹ * (D * N4.J)) * Bᵀ := by rw [Matrix.nonsing_inv_mul _ hunit]
        _ = (D * N4.J)⁻¹ * ((D * N4.J) * Bᵀ) := by simp [Matrix.mul_assoc]
        _ = 0 := by
          have htmp := congrArg (fun M => (D * N4.J)⁻¹ * M) h21
          simpa [Matrix.mul_assoc] using htmp
    exact ⟨by simpa [transpose_eq_zero] using hBt, hdet⟩
  · rintro ⟨hB, hdet⟩
    have hDD : D * N4.J * Dᵀ = N4.J := by
      rw [N4.mul_J_transpose_mul, hdet]
      simp [N4.J]
    ext i j
    cases i with
    | inl i =>
        cases j with
        | inl j =>
            simp [FixesBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply, hB]
        | inr j =>
            simp [FixesBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply, hB]
    | inr i =>
        cases j with
        | inl j =>
            simp [FixesBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply, hB]
        | inr j =>
            have hij := congrArg (fun M => M i j) hDD
            simpa [FixesBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- Sharp pointwise description for the weighted two-point orbit. -/
theorem fixesPair_fromBlocks_iff
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ C = 0 ∧ N4.FixesBivector (N4.splitRep₁ (k := k)) A ∧ D.det = 1 := by
  constructor
  · rintro ⟨h1, h2⟩
    rcases (fixes_rep₂_fromBlocks_iff (A := A) (B := B) (C := C) (D := D)).1 h2 with
      ⟨hB, hD⟩
    have hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h1
      simpa [FixesBivector, rep₁, hB, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        N4.FixesBivector, N4.ActBivector] using hij
    have hAunit : IsUnit (A.det) := by
      have hdetEq := congrArg Matrix.det hA
      simp [N4.FixesBivector, N4.ActBivector, Matrix.det_mul, Matrix.det_transpose, det_splitRep₁] at hdetEq
      exact ⟨⟨A.det, A.det, hdetEq, hdetEq⟩, rfl⟩
    have hCA : C * (N4.splitRep₁ (k := k) * Aᵀ) = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inl j)) h1
      simpa [Matrix.mul_assoc, FixesBivector, rep₁, hB, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
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
    exact ⟨hB, hC, hA, hD⟩
  · rintro ⟨hB, hC, hA, hD⟩
    refine ⟨?_, ?_⟩
    ·
      have hDD : D * N4.J * Dᵀ = N4.J := by
        rw [N4.mul_J_transpose_mul, hD]
        simp [N4.J]
      ext i j
      cases i with
      | inl i =>
          cases j with
          | inl j =>
              have hij := congrArg (fun M => M i j) hA
              simpa [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
                N4.FixesBivector, N4.ActBivector] using hij
          | inr j =>
              simp [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
      | inr i =>
          cases j with
          | inl j =>
              simp [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
          | inr j =>
              have hij := congrArg (fun M => M i j) hDD
              simpa [FixesBivector, rep₁, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    · exact (fixes_rep₂_fromBlocks_iff (A := A) (B := B) (C := C) (D := D)).2 ⟨hB, hD⟩

/-- A pointwise family inside the weighted two-point stabilizer. -/
theorem pointwise_family
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1) :
    FixesPairBivector (Matrix.fromBlocks A 0 0 D) := by
  exact (fixesPair_fromBlocks_iff (A := A) (B := 0) (C := 0) (D := D)).2 ⟨rfl, rfl, hA, hD⟩

/-- A diagonal torus lift for the weighted two-point orbit. -/
theorem torus_lift_action
    (u v : k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix W W k := Matrix.fromBlocks A 0 0 A
    let E : Matrix I I k := !![v, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks H 0 0 E
    ActBivector rep₁ g = u • rep₁ + (v - u) • rep₂ ∧
      ActBivector rep₂ g = v • rep₂ := by
  dsimp
  have htop :=
    N4Summary.split_torus_lift_action (k := k) (u := u) (v := u)
  have htop1 :
      N4.ActBivector (N4.splitRep₁ (k := k))
        (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
          (!![u, 0; 0, 1] : Matrix N4.I N4.I k)) =
        u • N4.splitRep₁ := by
    simpa using htop.1
  have hbot : (!![v, 0; 0, 1] : Matrix I I k) * N4.J * (!![v, 0; 0, 1] : Matrix I I k)ᵀ = v • N4.J := by
    have hdet : Matrix.det (!![v, 0; 0, 1] : Matrix I I k) = v := by
      simp [Matrix.det_fin_two]
    rw [N4.mul_J_transpose_mul, hdet]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [N4.J]
  constructor
  ·
    calc
      ActBivector rep₁ (Matrix.fromBlocks (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0 (!![u, 0; 0, 1] : Matrix N4.I N4.I k)) 0 0 (!![v, 0; 0, 1] : Matrix I I k))
          = Matrix.fromBlocks
              (N4.ActBivector (N4.splitRep₁ (k := k))
                (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0 (!![u, 0; 0, 1] : Matrix N4.I N4.I k)))
              0 0
              ((!![v, 0; 0, 1] : Matrix I I k) * N4.J * (!![v, 0; 0, 1] : Matrix I I k)ᵀ) := by
                ext i j
                cases i <;> cases j <;>
                  simp [ActBivector, rep₁, N4.ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
      _ = Matrix.fromBlocks (u • N4.splitRep₁) 0 0 (v • N4.J) := by
            ext i j
            cases i with
            | inl i =>
                cases j with
                | inl j =>
                    have hij := congrArg (fun M => M i j) htop1
                    simpa [Matrix.fromBlocks] using hij
                | inr j =>
                    simp [Matrix.fromBlocks]
            | inr i =>
                cases j with
                | inl j =>
                    simp [Matrix.fromBlocks]
                | inr j =>
                    have hij := congrArg (fun M => M i j) hbot
                    simpa [Matrix.fromBlocks] using hij
      _ = u • rep₁ + (v - u) • rep₂ := by
            ext i j
            cases i with
            | inl i =>
                cases j <;> simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply, sub_eq_add_neg]
            | inr i =>
                cases j with
                | inl j =>
                    simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply, sub_eq_add_neg]
                | inr j =>
                    fin_cases i <;> fin_cases j <;>
                      simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply, sub_eq_add_neg] <;> ring
  ·
    calc
      ActBivector rep₂
          (Matrix.fromBlocks
            (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
              (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
            0 0
            (!![v, 0; 0, 1] : Matrix I I k)) =
        Matrix.fromBlocks 0 0 0
          ((!![v, 0; 0, 1] : Matrix I I k) * N4.J * (!![v, 0; 0, 1] : Matrix I I k)ᵀ) := by
            ext i j
            cases i <;> cases j <;>
              simp [ActBivector, rep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
      _ = Matrix.fromBlocks 0 0 0 (v • N4.J) := by
            ext i j
            cases i with
            | inl i =>
                cases j <;> simp [Matrix.fromBlocks]
            | inr i =>
                cases j with
                | inl j =>
                    simp [Matrix.fromBlocks]
                | inr j =>
                    have hij := congrArg (fun M => M i j) hbot
                    simpa [Matrix.fromBlocks] using hij
      _ = v • rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₂, Matrix.fromBlocks]

/-- The determinant of the standard torus lift on the weighted two-point orbit is
`u^2 v`. -/
theorem torus_lift_det
    (u v : k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix W W k := Matrix.fromBlocks A 0 0 A
    let E : Matrix I I k := !![v, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks H 0 0 E
    Matrix.det g = (u * u) * v := by
  dsimp
  rw [Matrix.det_fromBlocks_zero₂₁, Matrix.det_fromBlocks_zero₂₁]
  simp [Matrix.det_fin_two]

/-- The block-diagonal pointwise subgroup combines with the standard torus lift with
the expected quotient action. -/
theorem pointwise_torus_product_lift_action
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1)
    (u v : k) :
    let A0 : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix W W k := Matrix.fromBlocks A0 0 0 A0
    let E : Matrix I I k := !![v, 0; 0, 1]
    let gU : Matrix V V k := Matrix.fromBlocks A 0 0 D
    let gL : Matrix V V k := Matrix.fromBlocks H 0 0 E
    ActBivector rep₁ (gU * gL) = u • rep₁ + (v - u) • rep₂ ∧
      ActBivector rep₂ (gU * gL) = v • rep₂ := by
  dsimp
  have hU := pointwise_family (k := k) (A := A) (D := D) hA hD
  have hU1 : ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) = rep₁ := by
    simpa [FixesBivector, ActBivector] using hU.1
  have hU2 : ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) = rep₂ := by
    simpa [FixesBivector, ActBivector] using hU.2
  have hL := torus_lift_action (k := k) (u := u) (v := v)
  constructor
  · calc
      ActBivector rep₁
          ((Matrix.fromBlocks A 0 0 D) *
            Matrix.fromBlocks
              (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
                (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
              0 0 (!![v, 0; 0, 1] : Matrix I I k)) =
        ActBivector
          (ActBivector rep₁
            (Matrix.fromBlocks
              (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
                (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
              0 0 (!![v, 0; 0, 1] : Matrix I I k)))
          (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_mul]
      _ = ActBivector (u • rep₁ + (v - u) • rep₂) (Matrix.fromBlocks A 0 0 D) := by
            rw [hL.1]
      _ = ActBivector (u • rep₁) (Matrix.fromBlocks A 0 0 D) +
            ActBivector ((v - u) • rep₂) (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_add]
      _ = u • ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) +
            (v - u) • ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_smul, actBivector_smul]
      _ = u • rep₁ + (v - u) • rep₂ := by
            have hUr1 :
                u • ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) = u • rep₁ := by
              simpa using congrArg (fun M => u • M) hU1
            have hUr2 :
                (v - u) • ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) =
                  (v - u) • rep₂ := by
              simpa using congrArg (fun M => (v - u) • M) hU2
            rw [hUr1, hUr2]
  · calc
      ActBivector rep₂
          ((Matrix.fromBlocks A 0 0 D) *
            Matrix.fromBlocks
              (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
                (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
              0 0 (!![v, 0; 0, 1] : Matrix I I k)) =
        ActBivector
          (ActBivector rep₂
            (Matrix.fromBlocks
              (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
                (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
              0 0 (!![v, 0; 0, 1] : Matrix I I k)))
          (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_mul]
      _ = ActBivector (v • rep₂) (Matrix.fromBlocks A 0 0 D) := by
            rw [hL.2]
      _ = v • ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) := by
            rw [actBivector_smul]
      _ = v • rep₂ := by
            simpa using congrArg (fun M => v • M) hU2

/-- Right-multiplying the standard torus lift by the block-diagonal pointwise subgroup
does not change the quotient action on the weighted two-point orbit. -/
theorem torus_pointwise_right_product_lift_action
    (A : Matrix W W k)
    (D : Matrix I I k)
    (hA : N4.FixesBivector (N4.splitRep₁ (k := k)) A)
    (hD : D.det = 1)
    (u v : k) :
    let A0 : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let H : Matrix W W k := Matrix.fromBlocks A0 0 0 A0
    let E : Matrix I I k := !![v, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H 0 0 E
    let gU : Matrix V V k := Matrix.fromBlocks A 0 0 D
    ActBivector rep₁ (gL * gU) = u • rep₁ + (v - u) • rep₂ ∧
      ActBivector rep₂ (gL * gU) = v • rep₂ := by
  dsimp
  have hU := pointwise_family (k := k) (A := A) (D := D) hA hD
  have hU1 : ActBivector rep₁ (Matrix.fromBlocks A 0 0 D) = rep₁ := by
    simpa [FixesBivector, ActBivector] using hU.1
  have hU2 : ActBivector rep₂ (Matrix.fromBlocks A 0 0 D) = rep₂ := by
    simpa [FixesBivector, ActBivector] using hU.2
  have hL := torus_lift_action (k := k) (u := u) (v := v)
  constructor
  · calc
      ActBivector rep₁
          (Matrix.fromBlocks
              (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
                (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
              0 0 (!![v, 0; 0, 1] : Matrix I I k) *
            Matrix.fromBlocks A 0 0 D) =
        ActBivector (ActBivector rep₁ (Matrix.fromBlocks A 0 0 D))
          (Matrix.fromBlocks
            (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
              (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
            0 0 (!![v, 0; 0, 1] : Matrix I I k)) := by
            rw [actBivector_mul]
      _ = ActBivector rep₁
          (Matrix.fromBlocks
            (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
              (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
            0 0 (!![v, 0; 0, 1] : Matrix I I k)) := by
            rw [hU1]
      _ = u • rep₁ + (v - u) • rep₂ := by
            rw [hL.1]
  · calc
      ActBivector rep₂
          (Matrix.fromBlocks
              (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
                (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
              0 0 (!![v, 0; 0, 1] : Matrix I I k) *
            Matrix.fromBlocks A 0 0 D) =
        ActBivector (ActBivector rep₂ (Matrix.fromBlocks A 0 0 D))
          (Matrix.fromBlocks
            (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
              (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
            0 0 (!![v, 0; 0, 1] : Matrix I I k)) := by
            rw [actBivector_mul]
      _ = ActBivector rep₂
          (Matrix.fromBlocks
            (Matrix.fromBlocks (!![u, 0; 0, 1] : Matrix N4.I N4.I k) 0 0
              (!![u, 0; 0, 1] : Matrix N4.I N4.I k))
            0 0 (!![v, 0; 0, 1] : Matrix I I k)) := by
            rw [hU2]
      _ = v • rep₂ := by
            rw [hL.2]

end N6WeightedTwoPoint
end Wedge2Formalization
