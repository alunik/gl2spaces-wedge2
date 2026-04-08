import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N6OnePointPlusSimple

variable {k : Type*} [Field k]

/-- The `4`-dimensional repeated-support block. -/
abbrev W := N4.V

/-- The `2`-dimensional simple block. -/
abbrev I := N4.I

/-- The ambient `6`-dimensional space. -/
abbrev V := W ⊕ I

/-- The natural `GL(V)` action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The bivector action of a product factors through successive actions. -/
theorem actBivector_mul
    (Ω : Matrix V V k) (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  rw [ActBivector, ActBivector, ActBivector, Matrix.transpose_mul]
  repeat rw [Matrix.mul_assoc]

/-- The direct-sum representative with divisor `2[a]+[b]`. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks (N4.onePointRep₁ (k := k)) 0 0 (N4.J (k := k))

/-- The second basis vector of the direct-sum representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks (N4.onePointRep₂ (k := k)) 0 0 (N4.J (k := k))

/-- Exact stabilizer of the direct-sum pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  ActBivector (rep₁ (k := k)) g = rep₁ ∧
    ActBivector (rep₂ (k := k)) g = rep₂

/-- The difference of the two basis vectors is the embedded nondegenerate top block. -/
def topRep : Matrix V V k :=
  Matrix.fromBlocks (N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) 0 0 0

lemma rep₁_sub_rep₂ : rep₁ (k := k) - rep₂ (k := k) = topRep (k := k) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₁, rep₂, topRep, Matrix.fromBlocks, Matrix.sub_apply]

lemma det_topRep : Matrix.det (N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) = 1 := by
  simpa [sub_eq_add_neg] using
    (N4Summary.onePoint_det_in_basis (k := k) (a := (1 : k)) (b := (-1 : k)))

/-- In the chosen basis, every vector on this line is block diagonal. -/
theorem basis_linearCombination
    (a b : k) :
    a • (rep₁ (k := k)) + b • (rep₂ (k := k)) =
      Matrix.fromBlocks
        (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k)))
        0
        0
        ((a + b) • (N4.J (k := k))) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply, add_smul, smul_add, add_assoc]

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

/-- Fixing the embedded nondegenerate top block forces both off-diagonal blocks to
vanish. -/
theorem fixes_topRep_fromBlocks_iff
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    ActBivector (topRep (k := k)) (Matrix.fromBlocks A B C D) = topRep ↔
      C = 0 ∧
      N4.FixesBivector (N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) A := by
  constructor
  · intro h
    have hA :
        N4.FixesBivector (N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) A := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h
      simpa [ActBivector, topRep, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        N4.FixesBivector, N4.ActBivector] using hij
    have hAunit : IsUnit (A.det) := by
      have hdetEq := congrArg Matrix.det hA
      simp [N4.FixesBivector, N4.ActBivector, Matrix.det_mul, Matrix.det_transpose, det_topRep] at hdetEq
      exact ⟨⟨A.det, A.det, hdetEq, hdetEq⟩, rfl⟩
    have hCA : C * ((N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) * Aᵀ) = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inr i) (Sum.inl j)) h
      simpa [ActBivector, topRep, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_assoc] using hij
    have hunitRight :
        IsUnit (((N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) * Aᵀ).det) := by
      rw [Matrix.det_mul, Matrix.det_transpose, det_topRep]
      simpa using hAunit
    have hC : C = 0 := by
      calc
        C = C * (1 : Matrix W W k) := by simp
        _ =
            C *
              (((N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) * Aᵀ) *
                ((N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) * Aᵀ)⁻¹) := by
              rw [Matrix.mul_nonsing_inv _ hunitRight]
        _ =
            (C * ((N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) * Aᵀ)) *
              ((N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) * Aᵀ)⁻¹ := by
              simp [Matrix.mul_assoc]
        _ = 0 := by simpa [hCA]
    exact ⟨hC, hA⟩
  · rintro ⟨hC, hA⟩
    ext i j
    cases i with
    | inl i =>
        cases j with
        | inl j =>
            have hij := congrArg (fun M => M i j) hA
            simpa [ActBivector, topRep, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
              N4.FixesBivector, N4.ActBivector] using hij
        | inr j =>
            simp [ActBivector, topRep, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
    | inr i =>
        cases j with
        | inl j =>
            simp [ActBivector, topRep, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]
        | inr j =>
            simp [ActBivector, topRep, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

/-- Sharp pointwise description for the direct-sum `2[a]+[b]` orbit. -/
theorem fixesPair_fromBlocks_iff
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      B = 0 ∧ C = 0 ∧ N4.FixesOnePointPairBivector A ∧ D.det = 1 := by
  constructor
  · rintro ⟨h1, h2⟩
    have htop :
        ActBivector (topRep (k := k)) (Matrix.fromBlocks A B C D) = topRep := by
      calc
        ActBivector (topRep (k := k)) (Matrix.fromBlocks A B C D) =
            ActBivector (rep₁ (k := k) - rep₂ (k := k)) (Matrix.fromBlocks A B C D) := by
              rw [rep₁_sub_rep₂]
        _ =
            ActBivector (rep₁ (k := k)) (Matrix.fromBlocks A B C D) -
              ActBivector (rep₂ (k := k)) (Matrix.fromBlocks A B C D) := by
              simp [ActBivector, Matrix.mul_sub, Matrix.sub_mul]
        _ = topRep (k := k) := by rw [h1, h2, rep₁_sub_rep₂]
    rcases (fixes_topRep_fromBlocks_iff (k := k) (A := A) (B := B) (C := C) (D := D)).1 htop with
      ⟨hC, hAtop⟩
    have hD : D.det = 1 := by
      have hDD : D * N4.J * Dᵀ = N4.J := by
        ext i j
        have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h2
        simpa [ActBivector, rep₂, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
      rw [N4.mul_J_transpose_mul] at hDD
      have h01 := congrArg (fun M => M 0 1) hDD
      simpa [N4.J] using h01
    have hBJD : B * N4.J * Dᵀ = (0 : Matrix W I k) := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h2
      simpa [ActBivector, rep₂, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        Matrix.mul_assoc] using hij
    have hunitBot : IsUnit ((N4.J (k := k) * Dᵀ).det) := by
      rw [Matrix.det_mul, Matrix.det_transpose, N4.J_det, hD]
      simp
    have hB : B = 0 := by
      have hBt : B = 0 := by
        calc
          B = B * (1 : Matrix I I k) := by simp
          _ = B * ((N4.J (k := k) * Dᵀ) * (N4.J (k := k) * Dᵀ)⁻¹) := by
                rw [Matrix.mul_nonsing_inv _ hunitBot]
          _ = (B * (N4.J (k := k) * Dᵀ)) * (N4.J (k := k) * Dᵀ)⁻¹ := by
                simp [Matrix.mul_assoc]
          _ = 0 := by
                rw [← Matrix.mul_assoc, hBJD]
                simp
      exact hBt
    have hA2 : N4.FixesBivector (N4.onePointRep₂ (k := k)) A := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h2
      simpa [ActBivector, rep₂, hB, hC, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
        N4.FixesBivector, N4.ActBivector] using hij
    have hA1 : N4.FixesBivector (N4.onePointRep₁ (k := k)) A := by
      have hAtop' :
          N4.ActBivector (N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) A =
            N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k) := by
        simpa [N4.FixesBivector, N4.ActBivector] using hAtop
      have hA2' :
          N4.ActBivector (N4.onePointRep₂ (k := k)) A = N4.onePointRep₂ (k := k) := by
        simpa [N4.FixesBivector, N4.ActBivector] using hA2
      calc
        N4.ActBivector (N4.onePointRep₁ (k := k)) A =
            N4.ActBivector
              ((N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) + N4.onePointRep₂ (k := k)) A := by
              simp [sub_eq_add_neg, add_assoc]
        _ =
            N4.ActBivector (N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) A +
              N4.ActBivector (N4.onePointRep₂ (k := k)) A := by
              ext i j
              fin_cases i <;> fin_cases j <;>
                simp [N4.ActBivector, Matrix.mul_apply, Fin.sum_univ_four, sub_eq_add_neg] <;>
                ring
        _ =
            (N4.onePointRep₁ (k := k) - N4.onePointRep₂ (k := k)) + N4.onePointRep₂ (k := k) := by
              rw [hAtop', hA2']
        _ = N4.onePointRep₁ (k := k) := by simp [sub_eq_add_neg, add_assoc]
    exact ⟨hB, hC, ⟨hA1, hA2⟩, hD⟩
  · rintro ⟨hB, hC, hA, hD⟩
    have hEq : Matrix.fromBlocks A B C D = Matrix.fromBlocks A 0 0 D := by
      simp [hB, hC]
    have hA1 :
        N4.ActBivector (N4.onePointRep₁ (k := k)) A = N4.onePointRep₁ (k := k) := by
      simpa [N4.FixesBivector, N4.ActBivector] using hA.1
    have hA2 :
        N4.ActBivector (N4.onePointRep₂ (k := k)) A = N4.onePointRep₂ (k := k) := by
      simpa [N4.FixesBivector, N4.ActBivector] using hA.2
    constructor
    · rw [hEq]
      calc
        ActBivector (rep₁ (k := k)) (Matrix.fromBlocks A 0 0 D) =
          Matrix.fromBlocks (N4.ActBivector (N4.onePointRep₁ (k := k)) A) 0 0 (D * N4.J * Dᵀ) := by
            simpa [rep₁] using
              act_blockDiagonal
                (k := k)
                (ΩW := N4.onePointRep₁ (k := k))
                (ΩI := N4.J (k := k))
                (H := A)
                (E := D)
        _ = Matrix.fromBlocks (N4.onePointRep₁ (k := k)) 0 0 (N4.J (k := k)) := by
            rw [hA1]
            rw [N4.mul_J_transpose_mul, hD]
            simp [N4.J]
        _ = rep₁ := by simp [rep₁]
    · rw [hEq]
      calc
        ActBivector (rep₂ (k := k)) (Matrix.fromBlocks A 0 0 D) =
          Matrix.fromBlocks (N4.ActBivector (N4.onePointRep₂ (k := k)) A) 0 0 (D * N4.J * Dᵀ) := by
            simpa [rep₂] using
              act_blockDiagonal
                (k := k)
                (ΩW := N4.onePointRep₂ (k := k))
                (ΩI := N4.J (k := k))
                (H := A)
                (E := D)
        _ = Matrix.fromBlocks (N4.onePointRep₂ (k := k)) 0 0 (N4.J (k := k)) := by
            rw [hA2]
            rw [N4.mul_J_transpose_mul, hD]
            simp [N4.J]
        _ = rep₂ := by simp [rep₂]

/-- The obvious block diagonal pointwise family on the direct-sum orbit. -/
theorem pointwise_blockDiagonal_family
    (H : Matrix W W k)
    (E : Matrix I I k)
    (hH : N4.FixesOnePointPairBivector H)
    (hE : E.det = 1) :
    FixesPairBivector (Matrix.fromBlocks H 0 0 E) := by
  have htop1 : N4.ActBivector (N4.onePointRep₁ (k := k)) H = N4.onePointRep₁ (k := k) := by
    simpa [N4.FixesBivector, N4.ActBivector] using hH.1
  have htop2 : N4.ActBivector (N4.onePointRep₂ (k := k)) H = N4.onePointRep₂ (k := k) := by
    simpa [N4.FixesBivector, N4.ActBivector] using hH.2
  have hEJ : E * N4.J * Eᵀ = N4.J := by
    rw [N4.mul_J_transpose_mul, hE]
    simp [N4.J]
  constructor
  · calc
      ActBivector (rep₁ (k := k)) (Matrix.fromBlocks H 0 0 E) =
        Matrix.fromBlocks (N4.ActBivector (N4.onePointRep₁ (k := k)) H) 0 0 (E * N4.J * Eᵀ) := by
          simpa [rep₁] using
            act_blockDiagonal
              (k := k)
              (ΩW := N4.onePointRep₁ (k := k))
                (ΩI := N4.J (k := k))
                (H := H)
                (E := E)
      _ = Matrix.fromBlocks (N4.onePointRep₁ (k := k)) 0 0 (N4.J (k := k)) := by
          rw [htop1, hEJ]
      _ = rep₁ := by simp [rep₁]
  · calc
      ActBivector (rep₂ (k := k)) (Matrix.fromBlocks H 0 0 E) =
        Matrix.fromBlocks (N4.ActBivector (N4.onePointRep₂ (k := k)) H) 0 0 (E * N4.J * Eᵀ) := by
          simpa [rep₂] using
            act_blockDiagonal
              (k := k)
              (ΩW := N4.onePointRep₂ (k := k))
                (ΩI := N4.J (k := k))
                (H := H)
                (E := E)
      _ = Matrix.fromBlocks (N4.onePointRep₂ (k := k)) 0 0 (N4.J (k := k)) := by
          rw [htop2, hEJ]
      _ = rep₂ := by simp [rep₂]

/-- A concrete block-diagonal Levi family for the direct-sum orbit. -/
theorem pointwise_levi_family
    (A B C D E : Matrix I I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE : E.det = 1) :
    FixesPairBivector
      (Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E) := by
  have hH :
      N4.FixesOnePointPairBivector (Matrix.fromBlocks A B C D) := by
    exact
      (N4Summary.onePoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B)
        (C := C)
        (D := D)).2 ⟨hA, hC, hD, hrel⟩
  exact pointwise_blockDiagonal_family (k := k) _ _ hH hE

/-- The concrete block-diagonal Levi family on the direct-sum orbit has determinant `1`. -/
theorem pointwise_levi_det
    (A B C D E : Matrix I I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE : E.det = 1) :
    Matrix.det (Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E) = 1 := by
  rw [Matrix.det_fromBlocks_zero₂₁, hC, Matrix.det_fromBlocks_zero₂₁, hD, hA, hE]
  ring

/-- The determinant on the direct-sum line factors as the square of `a^2 (a+b)`. -/
theorem det_in_basis
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) =
      ((a * a) * (a + b)) * ((a * a) * (a + b)) := by
  rw [basis_linearCombination (k := k) (a := a) (b := b)]
  rw [Matrix.det_fromBlocks_zero₂₁]
  have htop := N4Summary.onePoint_det_in_basis (k := k) (a := a) (b := b)
  have hcard : Fintype.card I = 2 := by
    simp [I, N4.I]
  have hbot : Matrix.det ((a + b) • (N4.J (k := k))) = (a + b) * (a + b) := by
    simp [Matrix.det_smul, hcard, N4.J_det, pow_two]
  rw [htop, hbot]
  ring

/-- Rank drop on the direct-sum line occurs exactly at the two distinguished projective
points. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 := by
  rw [det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h0 | h0
    · rcases mul_eq_zero.mp h0 with haa | hab
      · rcases mul_eq_zero.mp haa with ha | ha
        · exact Or.inl ha
        · exact Or.inl ha
      · exact Or.inr hab
    · rcases mul_eq_zero.mp h0 with haa | hab
      · rcases mul_eq_zero.mp haa with ha | ha
        · exact Or.inl ha
        · exact Or.inl ha
      · exact Or.inr hab
  · rintro (ha | hab)
    · simp [ha]
    · simp [hab]

/-- A one-parameter torus-type lift for the direct-sum `2[a]+[b]` orbit. -/
theorem torus_lift_action
    (a : k)
    (ha : a ≠ 0) :
    let B : Matrix I I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix I I k := !![a * a, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks H 0 0 E
    ActBivector rep₁ g = a • rep₁ + (a * a - a) • rep₂ ∧
      ActBivector rep₂ g = (a * a) • rep₂ := by
  dsimp
  have htop := N4Summary.onePoint_borel_lift_action (k := k) (a := a) (b := a * a - a) ha
  have hE : Matrix.det (!![a * a, 0; 0, 1] : Matrix I I k) = a * a := by
    simp [Matrix.det_fin_two]
  constructor
  · calc
      ActBivector rep₁
          (Matrix.fromBlocks
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
              !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0])
            0 0
            (!![a * a, 0; 0, 1] : Matrix I I k)) =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₁ (k := k))
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
              !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]))
          0 0
          ((!![a * a, 0; 0, 1] : Matrix I I k) * N4.J * (!![a * a, 0; 0, 1] : Matrix I I k)ᵀ) := by
            simpa [rep₁] using
              act_blockDiagonal
                (k := k)
                (ΩW := N4.onePointRep₁ (k := k))
                (ΩI := N4.J (k := k))
                (H := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
                  !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0])
                (E := (!![a * a, 0; 0, 1] : Matrix I I k))
      _ = Matrix.fromBlocks
          (a • (N4.onePointRep₁ (k := k)) + (a * a - a) • (N4.onePointRep₂ (k := k)))
          0 0
          ((a * a) • (N4.J (k := k))) := by
            rcases htop with ⟨htop1, _⟩
            rw [htop1, N4.mul_J_transpose_mul, hE]
            simp [N4.J]
      _ = a • rep₁ + (a * a - a) • rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply, add_smul]
            ring_nf
  · calc
      ActBivector rep₂
          (Matrix.fromBlocks
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
              !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0])
            0 0
            (!![a * a, 0; 0, 1] : Matrix I I k)) =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₂ (k := k))
            (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
              !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]))
          0 0
          ((!![a * a, 0; 0, 1] : Matrix I I k) * N4.J * (!![a * a, 0; 0, 1] : Matrix I I k)ᵀ) := by
            simpa [rep₂] using
              act_blockDiagonal
                (k := k)
                (ΩW := N4.onePointRep₂ (k := k))
                (ΩI := N4.J (k := k))
                (H := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
                  !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0])
                (E := (!![a * a, 0; 0, 1] : Matrix I I k))
      _ = Matrix.fromBlocks
          ((a * a) • (N4.onePointRep₂ (k := k)))
          0 0
          ((a * a) • (N4.J (k := k))) := by
            rcases htop with ⟨_, htop2⟩
            rw [htop2, N4.mul_J_transpose_mul, hE]
            simp [N4.J]
      _ = (a * a) • rep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₂, Matrix.fromBlocks, Matrix.smul_apply]

/-- The determinant of the standard torus-type lift on the direct-sum
`2[a]+[b]` orbit is `a^4`. -/
theorem torus_lift_det
    (a : k) :
    let B : Matrix I I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix I I k := !![a * a, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks H 0 0 E
    Matrix.det g = (a * a) * (a * a) := by
  dsimp
  rw [Matrix.det_fromBlocks_zero₂₁, N4Summary.onePoint_borel_lift_det]
  simp [Matrix.det_fin_two]

/-- The full block-diagonal pointwise subgroup on the direct-sum `2[a]+[b]` orbit
combines with the torus-type quotient lift with the expected action. -/
theorem pointwise_torus_product_lift_action
    (A B C D E0 : Matrix I I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k := Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    let B0 : Matrix I I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix I I k := !![a * a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H 0 0 E
    ActBivector rep₁ (gU * gL) = a • rep₁ + (a * a - a) • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
  intro gU B0 H E gL
  have hU := pointwise_levi_family (k := k) A B C D E0 hA hC hD hrel hE0
  have hL := torus_lift_action (k := k) (a := a) ha
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • rep₁ + (a * a - a) • rep₂) gU := by
        simpa [gL, E, H, B0] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector rep₁ gU + (a * a - a) • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ + (a * a - a) • rep₂ := by
        rw [hU.1, hU.2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector ((a * a) • rep₂) gU := by
        simpa [gL, E, H, B0] using congrArg (fun M => ActBivector M gU) hL.2
      _ = (a * a) • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
      _ = (a * a) • rep₂ := by
        rw [hU.2]

/-- Right-multiplying the torus-type quotient lift by the full block-diagonal pointwise
subgroup on the direct-sum `2[a]+[b]` orbit does not change the quotient action. -/
theorem torus_pointwise_right_product_lift_action
    (A B C D E0 : Matrix I I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let B0 : Matrix I I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix I I k := !![a * a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H 0 0 E
    let gU : Matrix V V k := Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    ActBivector rep₁ (gL * gU) = a • rep₁ + (a * a - a) • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
  intro B0 H E gL gU
  have hL := torus_lift_action (k := k) (a := a) ha
  have hU := pointwise_levi_family (k := k) A B C D E0 hA hC hD hrel hE0
  constructor
  · calc
      ActBivector rep₁ (gL * gU) = ActBivector (ActBivector rep₁ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.1
      _ = a • rep₁ + (a * a - a) • rep₂ := by
        simpa [gL, E, H, B0] using hL.1
  · calc
      ActBivector rep₂ (gL * gU) = ActBivector (ActBivector rep₂ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.2
      _ = (a * a) • rep₂ := by
        simpa [gL, E, H, B0] using hL.2

/-- Left-multiplying the torus-type quotient lift by the full block-diagonal pointwise
subgroup on the direct-sum `2[a]+[b]` orbit does not change the determinant. -/
theorem pointwise_torus_product_lift_det
    (A B C D E0 : Matrix I I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k := Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    let B0 : Matrix I I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix I I k := !![a * a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H 0 0 E
    Matrix.det (gU * gL) = (a * a) * (a * a) := by
  intro gU B0 H E gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using pointwise_levi_det (k := k) A B C D E0 hA hC hD hrel hE0
  have hLdet : Matrix.det gL = (a * a) * (a * a) := by
    simpa [gL, E, H, B0] using torus_lift_det (k := k) (a := a)
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- Right-multiplying the torus-type quotient lift by the full block-diagonal pointwise
subgroup on the direct-sum `2[a]+[b]` orbit does not change the determinant. -/
theorem torus_pointwise_right_product_lift_det
    (A B C D E0 : Matrix I I k)
    (hA : A.det = 1)
    (hC : C = 0)
    (hD : D = A)
    (hrel : A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hE0 : E0.det = 1)
    (a : k)
    (ha : a ≠ 0) :
    let B0 : Matrix I I k := !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k :=
      N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E : Matrix I I k := !![a * a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H 0 0 E
    let gU : Matrix V V k := Matrix.fromBlocks (Matrix.fromBlocks A B C D) 0 0 E0
    Matrix.det (gL * gU) = (a * a) * (a * a) := by
  intro B0 H E gL gU
  have hLdet : Matrix.det gL = (a * a) * (a * a) := by
    simpa [gL, E, H, B0] using torus_lift_det (k := k) (a := a)
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using pointwise_levi_det (k := k) A B C D E0 hA hC hD hrel hE0
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

end N6OnePointPlusSimple
end Wedge2Formalization
