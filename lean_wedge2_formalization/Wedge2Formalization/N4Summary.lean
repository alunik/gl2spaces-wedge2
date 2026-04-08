import Wedge2Formalization.N4
import Wedge2Formalization.N4PureSingular
import Mathlib.LinearAlgebra.Matrix.SchurComplement

namespace Wedge2Formalization
namespace N4Summary

open Matrix

variable {k : Type*} [Field k]

/-- The bivector action is linear in the tensor variable. -/
theorem actBivector_add
    (Ω₁ Ω₂ : Matrix N4.V N4.V k) (g : Matrix N4.V N4.V k) :
    N4.ActBivector (Ω₁ + Ω₂) g = N4.ActBivector Ω₁ g + N4.ActBivector Ω₂ g := by
  simp [N4.ActBivector, Matrix.mul_add, Matrix.add_mul]

/-- The bivector action commutes with scalar multiplication in the tensor variable. -/
theorem actBivector_smul
    (c : k) (Ω : Matrix N4.V N4.V k) (g : Matrix N4.V N4.V k) :
    N4.ActBivector (c • Ω) g = c • N4.ActBivector Ω g := by
  simp [N4.ActBivector, Matrix.mul_smul, smul_mul_assoc]

/-- The bivector action of a product factors through successive actions. -/
theorem actBivector_mul
    (Ω : Matrix N4.V N4.V k) (g h : Matrix N4.V N4.V k) :
    N4.ActBivector Ω (g * h) = N4.ActBivector (N4.ActBivector Ω h) g := by
  rw [N4.ActBivector, N4.ActBivector, N4.ActBivector, Matrix.transpose_mul]
  repeat rw [Matrix.mul_assoc]

/-- The determinant of the bivector action scales by `det(g)^2`. -/
theorem actBivector_det
    (Ω : Matrix N4.V N4.V k) (g : Matrix N4.V N4.V k) :
    Matrix.det (N4.ActBivector Ω g) =
      Matrix.det g * Matrix.det Ω * Matrix.det g := by
  simp [N4.ActBivector, Matrix.det_mul, Matrix.det_transpose, mul_assoc, mul_left_comm, mul_comm]

/-- An invertible matrix cannot send a bivector to zero under the natural action. -/
theorem actBivector_eq_zero_iff_of_det_ne_zero
    (Ω : Matrix N4.V N4.V k) (g : Matrix N4.V N4.V k)
    (hg : Matrix.det g ≠ 0) :
    N4.ActBivector Ω g = 0 ↔ Ω = 0 := by
  letI : Invertible (Matrix.det g) := invertibleOfNonzero hg
  letI : Invertible g := Matrix.invertibleOfDetInvertible g
  have hunit : IsUnit (Matrix.det g) := isUnit_of_invertible (Matrix.det g)
  constructor
  · intro hzero
    have h' := congrArg (fun M => g⁻¹ * M * (g⁻¹)ᵀ) hzero
    simp [N4.ActBivector, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hunit,
      Matrix.mul_nonsing_inv _ hunit, Matrix.transpose_nonsing_inv] at h'
    exact h'
  · intro hΩ
    simp [hΩ, N4.ActBivector]

/-- Summary form of the split-support pointwise stabilizer calculation. -/
theorem split_pointwise_bivector_iff
    (A B C D : Matrix N4.I N4.I k) :
    N4.FixesSplitPairBivector (Matrix.fromBlocks A B C D) ↔
      A.det = 1 ∧ B = 0 ∧ C = 0 ∧ D.det = 1 :=
  N4.fixes_splitPairBivector_fromBlocks_iff (A := A) (B := B) (C := C) (D := D)

/-- Block diagonal matrices preserve the split-support `2`-space. -/
theorem split_diag_preserves
    (A D : Matrix N4.I N4.I k) :
    N4.PreservesSplitSubspaceBivector (Matrix.fromBlocks A 0 0 D) :=
  N4.split_diag_preserves_subspace (A := A) (D := D)

/-- The explicit action of a block diagonal matrix on the split-support basis. -/
theorem split_diag_action
    (A D : Matrix N4.I N4.I k) :
    N4.ActBivector N4.splitRep₁ (Matrix.fromBlocks A 0 0 D) =
        A.det • N4.splitRep₁ + (D.det - A.det) • N4.splitRep₂ ∧
      N4.ActBivector N4.splitRep₂ (Matrix.fromBlocks A 0 0 D) =
        D.det • N4.splitRep₂ := by
  constructor
  ·
    calc
      N4.ActBivector N4.splitRep₁ (Matrix.fromBlocks A 0 0 D)
          = A.det • N4.ω12 + D.det • N4.ω34 := by
              exact N4.split_diag_act_rep1 (A := A) (D := D)
      _ = A.det • N4.splitRep₁ + (D.det - A.det) • N4.splitRep₂ := by
            ext i j
            rcases i with i | i <;> rcases j with j | j
            · fin_cases i <;> fin_cases j <;>
                simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg, N4.J]
            · simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg]
            · simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg]
            · fin_cases i <;> fin_cases j <;>
                simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg, N4.J] <;> ring
  ·
    calc
      N4.ActBivector N4.splitRep₂ (Matrix.fromBlocks A 0 0 D)
          = D.det • N4.ω34 := by
              exact N4.split_diag_act_rep2 (A := A) (D := D)
      _ = D.det • N4.splitRep₂ := by simp [N4.splitRep₂]

/-- The block swap preserves the split-support `2`-space. -/
theorem split_swap_preserves :
    N4.PreservesSplitSubspaceBivector (N4.splitSwap (k := k)) :=
  N4.split_swap_preserves_subspace (k := k)

/-- The explicit action of the block swap on the split-support basis. -/
theorem split_swap_action :
    N4.ActBivector N4.splitRep₁ (N4.splitSwap (k := k)) = N4.splitRep₁ ∧
      N4.ActBivector N4.splitRep₂ (N4.splitSwap (k := k)) =
        N4.splitRep₁ - N4.splitRep₂ := by
  constructor
  ·
    calc
      N4.ActBivector N4.splitRep₁ (N4.splitSwap (k := k))
          = N4.ActBivector N4.ω12 (N4.splitSwap (k := k))
              + N4.ActBivector N4.ω34 (N4.splitSwap (k := k)) := by
                simp [N4.ActBivector, N4.splitRep₁, Matrix.mul_add, Matrix.add_mul]
      _ = N4.ω34 + N4.ω12 := by rw [N4.splitSwap_act_ω12, N4.splitSwap_act_ω34]
      _ = N4.splitRep₁ := by simp [N4.splitRep₁, add_comm]
  ·
    calc
      N4.ActBivector N4.splitRep₂ (N4.splitSwap (k := k))
          = N4.ActBivector N4.ω34 (N4.splitSwap (k := k)) := by
                simp [N4.ActBivector, N4.splitRep₂]
      _ = N4.ω12 := by rw [N4.splitSwap_act_ω34]
      _ = N4.splitRep₁ - N4.splitRep₂ := by simp [N4.splitRep₁, N4.splitRep₂]

/-- Anti-diagonal block matrices preserve the split-support `2`-space. -/
theorem split_antidiag_preserves
    (B C : Matrix N4.I N4.I k) :
    N4.PreservesSplitSubspaceBivector (N4.splitAntiDiag (k := k) B C) :=
  N4.split_antidiag_preserves_subspace (B := B) (C := C)

/-- The explicit action of an anti-diagonal block matrix on the split-support basis. -/
theorem split_antidiag_action
    (B C : Matrix N4.I N4.I k) :
    N4.ActBivector N4.splitRep₁ (N4.splitAntiDiag (k := k) B C) =
        B.det • N4.splitRep₁ + (C.det - B.det) • N4.splitRep₂ ∧
      N4.ActBivector N4.splitRep₂ (N4.splitAntiDiag (k := k) B C) =
        B.det • N4.splitRep₁ + (-B.det) • N4.splitRep₂ := by
  constructor
  ·
    calc
      N4.ActBivector N4.splitRep₁ (N4.splitAntiDiag (k := k) B C)
          = N4.ActBivector N4.ω12 (N4.splitAntiDiag (k := k) B C)
              + N4.ActBivector N4.ω34 (N4.splitAntiDiag (k := k) B C) := by
                simp [N4.ActBivector, N4.splitRep₁, Matrix.mul_add, Matrix.add_mul]
      _ = C.det • N4.ω34 + B.det • N4.ω12 := by
            rw [N4.split_antidiag_act_ω12, N4.split_antidiag_act_ω34]
      _ = B.det • N4.splitRep₁ + (C.det - B.det) • N4.splitRep₂ := by
            ext i j
            rcases i with i | i <;> rcases j with j | j
            · fin_cases i <;> fin_cases j <;>
                simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg, N4.J]
            · simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg]
            · simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg]
            · fin_cases i <;> fin_cases j <;>
                simp [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, smul_add, sub_eq_add_neg, N4.J] <;> ring
  ·
    calc
      N4.ActBivector N4.splitRep₂ (N4.splitAntiDiag (k := k) B C)
          = N4.ActBivector N4.ω34 (N4.splitAntiDiag (k := k) B C) := by
                simp [N4.ActBivector, N4.splitRep₂]
      _ = B.det • N4.ω12 := by rw [N4.split_antidiag_act_ω34]
      _ = B.det • N4.splitRep₁ + (-B.det) • N4.splitRep₂ := by
            simp [N4.splitRep₁, N4.splitRep₂, smul_add, sub_eq_add_neg]

/-- The determinant on the split-support `2`-space factors as the square of `ab`. -/
theorem split_det_linearCombination
    (a b : k) :
    Matrix.det (a • N4.ω12 + b • N4.ω34) = (a * b) * (a * b) :=
  N4.det_split_linearCombination (a := a) (b := b)

/-- Rewriting a linear combination from the chosen split-support basis into the
coordinate basis `ω12, ω34`. -/
theorem split_basis_linearCombination
    (a b : k) :
    a • (N4.splitRep₁ (k := k)) + b • (N4.splitRep₂ (k := k)) =
      a • (N4.ω12 (k := k)) + (a + b) • (N4.ω34 (k := k)) := by
  calc
    a • (N4.splitRep₁ (k := k)) + b • (N4.splitRep₂ (k := k))
        = a • ((N4.ω12 (k := k)) + (N4.ω34 (k := k))) + b • (N4.ω34 (k := k)) := by
            simp [N4.splitRep₁, N4.splitRep₂]
    _ = a • (N4.ω12 (k := k)) + a • (N4.ω34 (k := k)) + b • (N4.ω34 (k := k)) := by
          simp [smul_add, add_assoc]
    _ = a • (N4.ω12 (k := k)) + (a • (N4.ω34 (k := k)) + b • (N4.ω34 (k := k))) := by
          abel_nf
    _ = a • (N4.ω12 (k := k)) + (a + b) • (N4.ω34 (k := k)) := by
          rw [← add_smul]

/-- The determinant on the split-support line in the chosen basis. -/
theorem split_det_in_basis
    (a b : k) :
    Matrix.det (a • (N4.splitRep₁ (k := k)) + b • (N4.splitRep₂ (k := k))) =
      (a * (a + b)) * (a * (a + b)) := by
  rw [split_basis_linearCombination (k := k)]
  simpa [mul_assoc, mul_left_comm, mul_comm] using
    split_det_linearCombination (a := a) (b := a + b)

/-- On the split-support line, rank drop occurs exactly on the two distinguished lines
spanned by `splitRep₂` and `splitRep₁ - splitRep₂`. -/
theorem split_det_zero_iff
    (a b : k) :
    Matrix.det (a • (N4.splitRep₁ (k := k)) + b • (N4.splitRep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 := by
  rw [split_det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    have hprod2 : (a * (a + b)) * (a * (a + b)) = 0 := h
    have hprod : a * (a + b) = 0 := by
      rcases mul_eq_zero.mp hprod2 with hzero | hzero
      · exact hzero
      · exact hzero
    exact mul_eq_zero.mp hprod
  · rintro (ha | hab)
    · simp [ha]
    · simp [hab]

/-- If an invertible coefficient matrix on the split-support `2`-space preserves the two
rank-drop lines, then it is either of diagonal type or of anti-diagonal type in the
chosen basis. -/
theorem split_action_cases_of_rankDrop
    (α β γ δ : k)
    (hinv : α * δ - β * γ ≠ 0)
    (hrep2 : Matrix.det (γ • (N4.splitRep₁ (k := k)) + δ • (N4.splitRep₂ (k := k))) = 0)
    (hω12 :
      Matrix.det ((α - γ) • (N4.splitRep₁ (k := k)) + (β - δ) • (N4.splitRep₂ (k := k))) = 0) :
    (γ = 0 ∧ α + β = δ) ∨ (α = γ ∧ γ + δ = 0) := by
  have hcol2 : γ = 0 ∨ γ + δ = 0 := (split_det_zero_iff (k := k) (a := γ) (b := δ)).1 hrep2
  have hcol12 :
      α - γ = 0 ∨ (α - γ) + (β - δ) = 0 := by
    exact (split_det_zero_iff (k := k) (a := α - γ) (b := β - δ)).1 hω12
  rcases hcol2 with hγ0 | hγδ
  · rcases hcol12 with hαγ0 | hsum
    · exfalso
      have hαeqγ : α = γ := sub_eq_zero.mp hαγ0
      have hα0 : α = 0 := by simpa [hγ0] using hαeqγ
      apply hinv
      simp [hγ0, hα0]
    · left
      constructor
      · exact hγ0
      ·
        have hsum' : (α + β) - δ = 0 := by
          simpa [hγ0, sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hsum
        exact sub_eq_zero.mp hsum'
  · rcases hcol12 with hαγ0 | hsum
    · right
      constructor
      · exact sub_eq_zero.mp hαγ0
      · exact hγδ
    · exfalso
      have hδ : δ = -γ := by
        simpa using eq_neg_of_add_eq_zero_right hγδ
      have hab0 : α + β = 0 := by
        have hsum' : α + β - (γ + δ) = 0 := by
          simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hsum
        simpa [hγδ] using hsum'
      have hβ : β = -α := by
        simpa using eq_neg_of_add_eq_zero_right hab0
      apply hinv
      calc
        α * δ - β * γ = α * (-γ) - (-α) * γ := by simp [hδ, hβ]
        _ = 0 := by ring

/-- The standard torus lift for the split-support orbit. -/
theorem split_torus_lift_action
    (u v : k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let g : Matrix N4.V N4.V k := Matrix.fromBlocks A 0 0 D
    N4.ActBivector N4.splitRep₁ g = u • N4.splitRep₁ + (v - u) • N4.splitRep₂ ∧
      N4.ActBivector N4.splitRep₂ g = v • N4.splitRep₂ := by
  dsimp
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  simpa [A, D, Matrix.det_fin_two] using
    (split_diag_action
      (k := k)
      (A := A)
      (D := D))

/-- The determinant of the standard torus lift is `uv`. -/
theorem split_torus_lift_det
    (u v : k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let g : Matrix N4.V N4.V k := Matrix.fromBlocks A 0 0 D
    Matrix.det g = u * v := by
  dsimp
  rw [Matrix.det_fromBlocks_zero₂₁]
  simp [Matrix.det_fin_two]

/-- The block swap is an involution. -/
theorem split_swap_sq :
    N4.splitSwap (k := k) * N4.splitSwap (k := k) = 1 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [N4.splitSwap, Matrix.fromBlocks, Matrix.mul_apply, Fin.sum_univ_two]

/-- Composing the standard torus lift with the block swap gives the other coset in the split
normalizer. -/
theorem split_swapCoset_lift_action
    (u v : k) :
    let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
    let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
    let g : Matrix N4.V N4.V k := Matrix.fromBlocks A 0 0 D * N4.splitSwap (k := k)
    N4.ActBivector N4.splitRep₁ g = u • N4.splitRep₁ + (v - u) • N4.splitRep₂ ∧
      N4.ActBivector N4.splitRep₂ g = u • N4.splitRep₁ + (-u) • N4.splitRep₂ := by
  dsimp
  let A : Matrix N4.I N4.I k := !![u, 0; 0, 1]
  let D : Matrix N4.I N4.I k := !![v, 0; 0, 1]
  have hdiag :=
    split_diag_action
      (k := k)
      (A := A)
      (D := D)
  have hswap := split_swap_action (k := k)
  constructor
  · calc
      N4.ActBivector N4.splitRep₁
          (Matrix.fromBlocks A 0 0 D *
            N4.splitSwap (k := k)) =
        N4.ActBivector
            (N4.ActBivector N4.splitRep₁ (N4.splitSwap (k := k)))
            (Matrix.fromBlocks A 0 0 D) := by
              rw [actBivector_mul]
      _ = N4.ActBivector N4.splitRep₁
            (Matrix.fromBlocks A 0 0 D) := by
            rw [hswap.1]
      _ = u • N4.splitRep₁ + (v - u) • N4.splitRep₂ := by
            simpa [A, D, Matrix.det_fin_two] using hdiag.1
  · calc
      N4.ActBivector N4.splitRep₂
          (Matrix.fromBlocks A 0 0 D *
            N4.splitSwap (k := k)) =
        N4.ActBivector
            (N4.ActBivector N4.splitRep₂ (N4.splitSwap (k := k)))
            (Matrix.fromBlocks A 0 0 D) := by
              rw [actBivector_mul]
      _ = N4.ActBivector
            (N4.splitRep₁ + (-1 : k) • N4.splitRep₂)
            (Matrix.fromBlocks A 0 0 D) := by
            rw [hswap.2]
            simp [sub_eq_add_neg]
      _ = N4.ActBivector N4.splitRep₁
            (Matrix.fromBlocks A 0 0 D) +
          (-1 : k) •
            N4.ActBivector N4.splitRep₂
              (Matrix.fromBlocks A 0 0 D) := by
              rw [actBivector_add, actBivector_smul]
      _ = (u • N4.splitRep₁ + (v - u) • N4.splitRep₂) + (-1 : k) • (v • N4.splitRep₂) := by
            rw [hdiag.1, hdiag.2]
            simp [A, D, Matrix.det_fin_two]
      _ = u • N4.splitRep₁ + (-u) • N4.splitRep₂ := by
            have hcoeff : (v - u) + (-v) = -u := by ring
            simp [sub_eq_add_neg, add_smul, hcoeff, add_assoc]

/-- A basis permutation used to put the repeated-support family into a Schur-complement form. -/
def onePointDetPerm : N4.V ≃ N4.V where
  toFun
    | Sum.inl 0 => Sum.inl 0
    | Sum.inl 1 => Sum.inr 1
    | Sum.inr 0 => Sum.inl 1
    | Sum.inr 1 => Sum.inr 0
  invFun
    | Sum.inl 0 => Sum.inl 0
    | Sum.inl 1 => Sum.inr 0
    | Sum.inr 0 => Sum.inr 1
    | Sum.inr 1 => Sum.inl 1
  left_inv := by
    intro x
    cases x with
    | inl i =>
        fin_cases i <;> rfl
    | inr i =>
        fin_cases i <;> rfl
  right_inv := by
    intro x
    cases x with
    | inl i =>
        fin_cases i <;> rfl
    | inr i =>
        fin_cases i <;> rfl

/-- Rewriting a linear combination on the repeated-support line in block form. -/
theorem onePoint_basis_linearCombination
    (a b : k) :
    a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k)) =
      Matrix.fromBlocks (b • (N4.J (k := k))) (a • (N4.J (k := k)))
        (a • (N4.J (k := k))) 0 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    simp [N4.onePointRep₁, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.add_apply]

/-- The inverse of the standard alternating `2 × 2` matrix is `-J`. -/
theorem J_inv :
    (N4.J (k := k))⁻¹ = -N4.J (k := k) := by
  apply Matrix.inv_eq_right_inv
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [N4.J, Matrix.mul_apply, Fin.sum_univ_two]

/-- After reindexing the repeated-support line, the determinant is computed by a Schur complement. -/
theorem onePoint_basis_linearCombination_submatrix
    (a b : k) :
    (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k))).submatrix
        onePointDetPerm onePointDetPerm =
      Matrix.fromBlocks
        (a • (N4.J (k := k)))
        (!![b, 0; 0, 0] : Matrix N4.I N4.I k)
        (!![-b, 0; 0, 0] : Matrix N4.I N4.I k)
        (-a • (N4.J (k := k))) := by
  rw [onePoint_basis_linearCombination (k := k) (a := a) (b := b)]
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [onePointDetPerm, N4.J, Matrix.fromBlocks]

/-- The determinant on the repeated-support line depends only on the coefficient of
`onePointRep₁`. -/
theorem onePoint_det_in_basis
    (a b : k) :
    Matrix.det (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k))) =
      (a * a) * (a * a) := by
  by_cases ha : a = 0
  · rw [ha]
    rw [onePoint_basis_linearCombination (k := k) (a := 0) (b := b)]
    simp [zero_smul]
    right
    simpa using (Matrix.det_zero (n := N4.I) (R := k))
  · rw [← Matrix.det_submatrix_equiv_self onePointDetPerm]
    rw [onePoint_basis_linearCombination_submatrix (k := k) (a := a) (b := b)]
    let A : Matrix N4.I N4.I k := a • (N4.J (k := k))
    let B : Matrix N4.I N4.I k := !![b, 0; 0, 0]
    let C : Matrix N4.I N4.I k := !![-b, 0; 0, 0]
    have hcard : Fintype.card N4.I = 2 := by simp [N4.I]
    have hdetA : Matrix.det A ≠ 0 := by
      simp [A, Matrix.det_smul, hcard, N4.J_det, pow_two, ha]
    letI : Invertible (Matrix.det A) := invertibleOfNonzero hdetA
    letI : Invertible A := Matrix.invertibleOfDetInvertible A
    letI : Invertible a := invertibleOfNonzero ha
    have hJunit : IsUnit (Matrix.det (N4.J (k := k))) := by
      simp [N4.J_det]
    have hAinv : A⁻¹ = a⁻¹ • (-N4.J (k := k)) := by
      apply Matrix.inv_eq_right_inv
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [A, N4.J, Matrix.mul_apply, Fin.sum_univ_two, ha]
    have hCB : C * A⁻¹ * B = 0 := by
      rw [hAinv]
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [B, C, N4.J, Matrix.mul_apply, Fin.sum_univ_two]
    have hdetBlock :
        Matrix.det (Matrix.fromBlocks A B C (-A)) =
          Matrix.det A * Matrix.det (-A - C * A⁻¹ * B) := by
      simpa using (Matrix.det_fromBlocks₁₁ A B C (-A))
    have hdetBlock' :
        Matrix.det
            (Matrix.fromBlocks (a • (N4.J (k := k))) (!![b, 0; 0, 0]) (!![-b, 0; 0, 0])
              (-a • (N4.J (k := k)))) =
          Matrix.det A * Matrix.det (-A - C * A⁻¹ * B) := by
      simpa [A, B, C] using hdetBlock
    have hdetAeq : Matrix.det A = a * a := by
      simp [A, Matrix.det_smul, hcard, N4.J_det, pow_two]
    have hdetNegA : Matrix.det (-A) = a * a := by
      rw [Matrix.det_neg, hcard, hdetAeq]
      norm_num
    rw [hdetBlock']
    rw [hCB]
    rw [sub_zero, hdetAeq, hdetNegA]

/-- The repeated-support line has a unique rank-drop point, namely the line spanned by
`onePointRep₂`. -/
theorem onePoint_det_zero_iff
    (a b : k) :
    Matrix.det (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k))) = 0 ↔
      a = 0 := by
  rw [onePoint_det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h0 | h0
    · exact mul_eq_zero.mp h0 |>.elim id id
    · exact mul_eq_zero.mp h0 |>.elim id id
  · intro ha
    simp [ha]

/-- Preserving the unique rank-drop line on the repeated-support orbit forces an upper-triangular
quotient action in the chosen basis. -/
theorem onePoint_action_case_of_rankDrop
    (γ δ : k)
    (hrep2 : Matrix.det (γ • (N4.onePointRep₁ (k := k)) + δ • (N4.onePointRep₂ (k := k))) = 0) :
    γ = 0 :=
  (onePoint_det_zero_iff (k := k) (a := γ) (b := δ)).1 hrep2

/-- Summary form of the repeated-support pointwise stabilizer calculation. -/
theorem onePoint_pointwise_bivector_iff
    (A B C D : Matrix N4.I N4.I k) :
    N4.FixesOnePointPairBivector (Matrix.fromBlocks A B C D) ↔
      A.det = 1 ∧ C = 0 ∧ D = A ∧ A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0 :=
  N4.fixes_onePointPairBivector_fromBlocks_iff (A := A) (B := B) (C := C) (D := D)

/-- Upper block shears preserve the repeated-support `2`-space. -/
theorem onePoint_upperShear_preserves
    (B : Matrix N4.I N4.I k) :
    N4.PreservesOnePointSubspaceBivector (N4.onePointUpperShear (k := k) B) :=
  N4.onePoint_upperShear_preserves_subspace (B := B)

/-- The explicit action of an upper block shear on the repeated-support basis. -/
theorem onePoint_upperShear_action
    (B : Matrix N4.I N4.I k) :
    N4.ActBivector N4.onePointRep₁ (N4.onePointUpperShear (k := k) B) =
        N4.onePointRep₁ + (B 0 0 + B 1 1) • N4.onePointRep₂ ∧
      N4.ActBivector N4.onePointRep₂ (N4.onePointUpperShear (k := k) B) =
        N4.onePointRep₂ := by
  constructor
  · exact N4.onePoint_upperShear_act_rep1 (B := B)
  · exact N4.onePoint_upperShear_act_rep2 (B := B)

/-- Diagonal scalings preserve the repeated-support `2`-space. -/
theorem onePoint_scale_preserves
    (a : k) :
    N4.PreservesOnePointSubspaceBivector (N4.onePointScale (k := k) a) :=
  N4.onePoint_scale_preserves_subspace (a := a)

/-- The explicit action of a diagonal scaling on the repeated-support basis. -/
theorem onePoint_scale_action
    (a : k) :
    N4.ActBivector N4.onePointRep₁ (N4.onePointScale (k := k) a) =
        a • N4.onePointRep₁ ∧
      N4.ActBivector N4.onePointRep₂ (N4.onePointScale (k := k) a) =
        (a * a) • N4.onePointRep₂ := by
  constructor
  · exact N4.onePoint_scale_act_rep1 (a := a)
  ·
    calc
      N4.ActBivector N4.onePointRep₂ (N4.onePointScale (k := k) a)
          = a ^ 2 • N4.onePointRep₂ := by
              exact N4.onePoint_scale_act_rep2 (a := a)
      _ = (a * a) • N4.onePointRep₂ := by simp [pow_two]

/-- The upper shear used in the repeated-support family always has determinant `1`. -/
theorem onePoint_upperShear_det
    (B : Matrix N4.I N4.I k) :
    Matrix.det (N4.onePointUpperShear (k := k) B) = 1 := by
  rw [N4.onePointUpperShear, Matrix.det_fromBlocks_zero₂₁]
  simp

/-- The scaling matrix used in the repeated-support family has determinant `a^2`. -/
theorem onePoint_scale_det
    (a : k) :
    Matrix.det (N4.onePointScale (k := k) a) = a * a := by
  rw [N4.onePointScale, Matrix.det_fromBlocks_zero₂₁]
  simp [Matrix.det_fin_two, pow_two]

/-- The standard upper-triangular quotient family for the repeated-support orbit lifts to an
explicit `GL₄` action. -/
theorem onePoint_borel_lift_action
    (a b : k)
    (ha : a ≠ 0) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let g := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    N4.ActBivector N4.onePointRep₁ g = a • N4.onePointRep₁ + b • N4.onePointRep₂ ∧
      N4.ActBivector N4.onePointRep₂ g = (a * a) • N4.onePointRep₂ := by
  dsimp
  let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
  have hshear :=
    onePoint_upperShear_action
      (k := k) (B := B)
  have hscale := onePoint_scale_action (k := k) (a := a)
  constructor
  · calc
      N4.ActBivector N4.onePointRep₁
          (N4.onePointScale (k := k) a *
            N4.onePointUpperShear (k := k) B) =
        N4.ActBivector
            (N4.ActBivector N4.onePointRep₁
              (N4.onePointUpperShear (k := k) B))
            (N4.onePointScale (k := k) a) := by
              rw [actBivector_mul]
      _ = N4.ActBivector
            (N4.onePointRep₁ + (B 0 0 + B 1 1) • N4.onePointRep₂)
            (N4.onePointScale (k := k) a) := by
              rw [hshear.1]
      _ = N4.ActBivector
            (N4.onePointRep₁ + (b * (a⁻¹ * a⁻¹)) • N4.onePointRep₂)
            (N4.onePointScale (k := k) a) := by
              simp [B]
      _ =
        N4.ActBivector N4.onePointRep₁ (N4.onePointScale (k := k) a) +
          (b * (a⁻¹ * a⁻¹)) •
            N4.ActBivector N4.onePointRep₂ (N4.onePointScale (k := k) a) := by
              rw [actBivector_add, actBivector_smul]
      _ = a • N4.onePointRep₁ +
            (b * (a⁻¹ * a⁻¹)) • ((a * a) • N4.onePointRep₂) := by
              rw [hscale.1, hscale.2]
      _ = a • N4.onePointRep₁ + ((b * (a⁻¹ * a⁻¹)) * (a * a)) • N4.onePointRep₂ := by
            simp [smul_smul, mul_assoc]
      _ = a • N4.onePointRep₁ + b • N4.onePointRep₂ := by
            have hab : (b * (a⁻¹ * a⁻¹)) * (a * a) = b := by
              field_simp [ha]
            simp [hab]
  · calc
      N4.ActBivector N4.onePointRep₂
          (N4.onePointScale (k := k) a *
            N4.onePointUpperShear (k := k) B) =
        N4.ActBivector
            (N4.ActBivector N4.onePointRep₂
              (N4.onePointUpperShear (k := k) B))
            (N4.onePointScale (k := k) a) := by
              rw [actBivector_mul]
      _ = N4.ActBivector N4.onePointRep₂ (N4.onePointScale (k := k) a) := by
            rw [hshear.2]
      _ = (a * a) • N4.onePointRep₂ := by
            rw [hscale.2]

/-- The determinant of the standard repeated-support lift is `a^2`. -/
theorem onePoint_borel_lift_det
    (a b : k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let g := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    Matrix.det g = a * a := by
  dsimp
  rw [Matrix.det_mul, onePoint_scale_det, onePoint_upperShear_det]
  ring

/-- Summary form of the pure singular pointwise stabilizer calculation. -/
theorem pureSingular_pointwise_iff_shape
    (g : Matrix N4PureSingular.I N4PureSingular.I k) :
    N4PureSingular.FixesPureSingularPair g ↔
      g 0 0 ≠ 0 ∧
        g =
          N4PureSingular.pureSingularShape (k := k)
            (g 0 0) (g 0 1) (g 0 2) (g 0 3) (g 1 3) (g 2 3) (g 3 3) :=
  N4PureSingular.fixesPureSingularPair_iff_shape (g := g)

/-- A large explicit family preserving the pure singular `2`-space. -/
theorem pureSingular_setwiseShape_preserves
    (a b c d p q r s e f t : k) :
    N4PureSingular.PreservesPureSingularSubspaceBivector
      (N4PureSingular.pureSingularSetwiseShape (k := k) a b c d p q r s e f t) :=
  N4PureSingular.pureSingularSetwiseShape_preserves
    (k := k) a b c d p q r s e f t

/-- The explicit action of the pure singular preserving family on the chosen basis. -/
theorem pureSingular_setwiseShape_action
    (a b c d p q r s e f t : k) :
    N4PureSingular.ActBivector N4PureSingular.ω12
        (N4PureSingular.pureSingularSetwiseShape (k := k) a b c d p q r s e f t) =
          (a * p) • N4PureSingular.ω12 + (a * r) • N4PureSingular.ω13 ∧
      N4PureSingular.ActBivector N4PureSingular.ω13
        (N4PureSingular.pureSingularSetwiseShape (k := k) a b c d p q r s e f t) =
          (a * q) • N4PureSingular.ω12 + (a * s) • N4PureSingular.ω13 := by
  constructor
  · exact N4PureSingular.pureSingularSetwiseShape_act_ω12 (k := k) a b c d p q r s e f t
  · exact N4PureSingular.pureSingularSetwiseShape_act_ω13 (k := k) a b c d p q r s e f t

/-- The pure singular bivector action is linear in the tensor variable. -/
theorem pureSingular_actBivector_add
    (Ω₁ Ω₂ : Matrix N4PureSingular.I N4PureSingular.I k)
    (g : Matrix N4PureSingular.I N4PureSingular.I k) :
    N4PureSingular.ActBivector (Ω₁ + Ω₂) g =
      N4PureSingular.ActBivector Ω₁ g + N4PureSingular.ActBivector Ω₂ g := by
  simp [N4PureSingular.ActBivector, Matrix.mul_add, Matrix.add_mul]

/-- The pure singular bivector action commutes with scalar multiplication in the tensor
variable. -/
theorem pureSingular_actBivector_smul
    (c : k)
    (Ω : Matrix N4PureSingular.I N4PureSingular.I k)
    (g : Matrix N4PureSingular.I N4PureSingular.I k) :
    N4PureSingular.ActBivector (c • Ω) g = c • N4PureSingular.ActBivector Ω g := by
  simp [N4PureSingular.ActBivector, Matrix.mul_smul, smul_mul_assoc]

/-- The determinant of the pure singular bivector action scales by `det(g)^2`. -/
theorem pureSingular_actBivector_det
    (Ω : Matrix N4PureSingular.I N4PureSingular.I k)
    (g : Matrix N4PureSingular.I N4PureSingular.I k) :
    Matrix.det (N4PureSingular.ActBivector Ω g) =
      Matrix.det g * Matrix.det Ω * Matrix.det g := by
  simp [N4PureSingular.ActBivector, Matrix.det_mul, Matrix.det_transpose,
    mul_assoc, mul_left_comm, mul_comm]

/-- An invertible matrix cannot send a pure singular bivector to zero under the natural
action. -/
theorem pureSingular_actBivector_eq_zero_iff_of_det_ne_zero
    (Ω : Matrix N4PureSingular.I N4PureSingular.I k)
    (g : Matrix N4PureSingular.I N4PureSingular.I k)
    (hg : Matrix.det g ≠ 0) :
    N4PureSingular.ActBivector Ω g = 0 ↔ Ω = 0 := by
  letI : Invertible (Matrix.det g) := invertibleOfNonzero hg
  letI : Invertible g := Matrix.invertibleOfDetInvertible g
  have hunit : IsUnit (Matrix.det g) := isUnit_of_invertible (Matrix.det g)
  constructor
  · intro hzero
    have h' := congrArg (fun M => g⁻¹ * M * (g⁻¹)ᵀ) hzero
    simp [N4PureSingular.ActBivector, Matrix.mul_assoc, Matrix.nonsing_inv_mul _ hunit,
      Matrix.mul_nonsing_inv _ hunit, Matrix.transpose_nonsing_inv] at h'
    exact h'
  · intro hΩ
    simp [hΩ, N4PureSingular.ActBivector]

/-- Every `2 × 2` coefficient matrix lifts to an explicit pure singular setwise stabilizer. -/
theorem pureSingular_GL2_lift_action
    (α β γ δ : k) :
    let g :=
      N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1
    N4PureSingular.ActBivector N4PureSingular.ω12 g =
        α • N4PureSingular.ω12 + β • N4PureSingular.ω13 ∧
      N4PureSingular.ActBivector N4PureSingular.ω13 g =
        γ • N4PureSingular.ω12 + δ • N4PureSingular.ω13 := by
  dsimp
  constructor
  · simpa using
      (N4PureSingular.pureSingularSetwiseShape_act_ω12
        (k := k) (a := 1) (b := 0) (c := 0) (d := 0)
        (p := α) (q := γ) (r := β) (s := δ) (e := 0) (f := 0) (t := 1))
  · simpa using
      (N4PureSingular.pureSingularSetwiseShape_act_ω13
        (k := k) (a := 1) (b := 0) (c := 0) (d := 0)
        (p := α) (q := γ) (r := β) (s := δ) (e := 0) (f := 0) (t := 1))

/-- The determinant of the explicit pure singular lift is the determinant of the coefficient
matrix on the `2`-space. -/
theorem pureSingular_GL2_lift_det
    (α β γ δ : k) :
    Matrix.det
        (N4PureSingular.pureSingularSetwiseShape (k := k) 1 0 0 0 α γ β δ 0 0 1) =
      α * δ - β * γ := by
  change
    Matrix.det
        (Matrix.fromBlocks (!![(1 : k)]) 0 0
          (!![α, γ, 0;
             β, δ, 0;
             0, 0, 1] : Matrix (Fin 3) (Fin 3) k)) =
      α * δ - β * γ
  rw [Matrix.det_fromBlocks_zero₂₁]
  simp [Matrix.det_fin_three]
  ring

end N4Summary
end Wedge2Formalization
