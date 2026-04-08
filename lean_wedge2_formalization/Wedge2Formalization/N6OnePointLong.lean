import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Tactic

open Matrix

namespace Wedge2Formalization
namespace N6OnePointLong

variable {k : Type*} [Field k]

/-- The `3`-dimensional support block. -/
abbrev I := Fin 3

/-- The ambient `6`-dimensional space. -/
abbrev V := I ⊕ I

/-- The natural action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- Exact stabilizer of a bivector. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  ActBivector Ω g = Ω

/-- The bivector action is multiplicative in the acting matrix. -/
theorem actBivector_mul
    (Ω : Matrix V V k)
    (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  simp [ActBivector, Matrix.mul_assoc]

/-- The nilpotent Jordan block of size `3`. -/
def N : Matrix I I k := !![(0 : k), 1, 0; 0, 0, 1; 0, 0, 0]

/-- The indecomposable length-three repeated-support representative. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 (1 : Matrix I I k) (-(1 : Matrix I I k)) 0

/-- The second basis vector of the length-three repeated-support representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 (N (k := k)) (-((N (k := k))ᵀ)) 0

/-- Exact stabilizer of the pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector (rep₁ (k := k)) g ∧ FixesBivector (rep₂ (k := k)) g

/-- In the chosen basis, every vector on the indecomposable `3[a]` line is determined
by the upper-right block `a I + b N`. -/
theorem basis_linearCombination
    (a b : k) :
    a • (rep₁ (k := k)) + b • (rep₂ (k := k)) =
      Matrix.fromBlocks
        0
        (a • (1 : Matrix I I k) + b • (N (k := k)))
        (-((a • (1 : Matrix I I k) + b • (N (k := k)))ᵀ))
        0 := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₁, rep₂, Matrix.add_apply, add_assoc, add_comm, add_left_comm]

/-- Swapping the two column blocks turns the indecomposable `3[a]` line into a block
diagonal determinant calculation. -/
theorem basis_linearCombination_submatrix
    (a b : k) :
    (a • (rep₁ (k := k)) + b • (rep₂ (k := k))).submatrix id (Equiv.sumComm I I) =
      Matrix.fromBlocks
        (a • (1 : Matrix I I k) + b • (N (k := k)))
        0
        0
        (-((a • (1 : Matrix I I k) + b • (N (k := k)))ᵀ)) := by
  rw [basis_linearCombination (k := k) (a := a) (b := b)]
  ext i j
  cases i <;> cases j <;>
    simp [Matrix.fromBlocks]

/-- The determinant of the upper-right block on the indecomposable `3[a]` line is
`a^3`. -/
theorem det_upperBlock
    (a b : k) :
    Matrix.det (a • (1 : Matrix I I k) + b • (N (k := k))) = a * a * a := by
  simp [N, Matrix.det_fin_three, Matrix.add_apply, Matrix.smul_apply]

/-- The determinant on the indecomposable `3[a]` line factors as the square of `a^3`. -/
theorem det_in_basis
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) =
      (a * a * a) * (a * a * a) := by
  have hsign : Equiv.Perm.sign (Equiv.sumComm I I) = (-1 : Units ℤ) := by
    decide
  have hperm :
      Matrix.det ((a • (rep₁ (k := k)) + b • (rep₂ (k := k))).submatrix id (Equiv.sumComm I I)) =
        -Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) := by
    simpa [hsign] using
      (Matrix.det_permute' (Equiv.sumComm I I) (a • (rep₁ (k := k)) + b • (rep₂ (k := k))))
  have hdet :
      Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) =
        -Matrix.det ((a • (rep₁ (k := k)) + b • (rep₂ (k := k))).submatrix id (Equiv.sumComm I I)) := by
    rw [hperm]
    ring
  rw [hdet, basis_linearCombination_submatrix (k := k) (a := a) (b := b)]
  rw [Matrix.det_fromBlocks_zero₂₁]
  rw [det_upperBlock (k := k) (a := a) (b := b)]
  have hbot :
      Matrix.det (-((a • (1 : Matrix I I k) + b • (N (k := k)))ᵀ)) =
        -(a * a * a) := by
    rw [Matrix.det_neg]
    rw [Matrix.det_transpose, det_upperBlock (k := k) (a := a) (b := b)]
    norm_num
  rw [hbot]
  ring

/-- Rank drop on the indecomposable `3[a]` line occurs exactly at the repeated-support
projective point. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) = 0 ↔
      a = 0 := by
  rw [det_in_basis (k := k) (a := a) (b := b)]
  constructor
  · intro h
    rw [mul_eq_zero] at h
    rcases h with h | h
    · rw [mul_eq_zero] at h
      rcases h with h | h
      · rw [mul_eq_zero] at h
        exact h.elim id id
      · exact h
    · rw [mul_eq_zero] at h
      rcases h with h | h
      · rw [mul_eq_zero] at h
        exact h.elim id id
      · exact h
  · rintro rfl
    ring

/-- A symmetric upper Hankel matrix. -/
def upperHankel (u v w : k) : Matrix I I k :=
  !![u, v, w;
     v, w, 0;
     w, 0, 0]

/-- A symmetric lower Hankel matrix. -/
def lowerHankel (u v w : k) : Matrix I I k :=
  !![0, 0, u;
     0, u, v;
     u, v, w]

/-- The upper block used in the quotient shear lift. -/
def shearTop (b : k) : Matrix I I k :=
  !![1, 0, 0;
     0, 1, b;
     0, 0, 1]

/-- The lower block used in the quotient shear lift. -/
def shearBottom (b : k) : Matrix I I k :=
  !![1, 0, 0;
     b, 1, 0;
     0, 0, 1]

/-- The upper block used in the quotient scaling lift. -/
def scaleTop (a : k) : Matrix I I k :=
  !![1, 0, 0;
     0, a, 0;
     0, 0, a * a]

/-- The lower block used in the quotient scaling lift. -/
def scaleBottom (a : k) : Matrix I I k :=
  !![a, 0, 0;
     0, 1, 0;
     0, 0, a⁻¹]

/-- The upper shear block is unipotent. -/
theorem shearTop_det
    (b : k) :
    Matrix.det (shearTop (k := k) b) = 1 := by
  simp [shearTop, Matrix.det_fin_three]

/-- The lower shear block is unipotent. -/
theorem shearBottom_det
    (b : k) :
    Matrix.det (shearBottom (k := k) b) = 1 := by
  simp [shearBottom, Matrix.det_fin_three]

/-- The upper scaling block has determinant `a^3`. -/
theorem scaleTop_det
    (a : k) :
    Matrix.det (scaleTop (k := k) a) = a * a * a := by
  simp [scaleTop, Matrix.det_fin_three]
  ring

/-- The lower scaling block has determinant `1` when `a ≠ 0`. -/
theorem scaleBottom_det
    (a : k)
    (ha : a ≠ 0) :
    Matrix.det (scaleBottom (k := k) a) = 1 := by
  simp [scaleBottom, Matrix.det_fin_three, ha]

/-- Acting on an off-diagonal bivector by a block matrix. -/
theorem act_fromBlocks
    (A B C D M : Matrix I I k) :
    ActBivector
        (Matrix.fromBlocks 0 M (-Mᵀ) 0)
        (Matrix.fromBlocks A B C D) =
      Matrix.fromBlocks
        (-B * Mᵀ * Aᵀ + A * M * Bᵀ)
        (-B * Mᵀ * Cᵀ + A * M * Dᵀ)
        (-D * Mᵀ * Aᵀ + C * M * Bᵀ)
        (-D * Mᵀ * Cᵀ + C * M * Dᵀ) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply, Matrix.mul_assoc]

/-- Recombining the off-diagonal blocks gives the expected linear combination of the
two chosen basis vectors. -/
theorem linearCombination_fromBlocks
    (a b : k) :
    Matrix.fromBlocks
      0
      (a • (1 : Matrix I I k) + b • (N (k := k)))
      (-((a • (1 : Matrix I I k) + b • (N (k := k)))ᵀ))
      0 =
      a • rep₁ + b • rep₂ := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · simp [rep₁, rep₂, N, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
  · fin_cases i <;> fin_cases j <;>
      simp [rep₁, rep₂, N, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] <;>
      ring
  · fin_cases i <;> fin_cases j <;>
      simp [rep₁, rep₂, N, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply,
        Matrix.transpose_add, Matrix.transpose_smul, Matrix.transpose_apply] <;>
      ring
  · simp [rep₁, rep₂, N, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]

/-- The upper Hankel matrix is symmetric. -/
theorem upperHankel_transpose
    (u v w : k) :
    (upperHankel (k := k) u v w)ᵀ = upperHankel (k := k) u v w := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [upperHankel, Matrix.transpose_apply]

/-- The upper Hankel matrix satisfies the commutation relation needed for `rep₂`. -/
theorem upperHankel_comm
    (u v w : k) :
    N (k := k) * upperHankel (k := k) u v w =
      upperHankel (k := k) u v w * (N (k := k))ᵀ := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [N, upperHankel, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]

/-- A symmetric matrix satisfying the upper commutation relation has the expected upper
Hankel shape. -/
theorem eq_upperHankel_of_transpose_eq_of_comm
    (B : Matrix I I k)
    (hSym : Bᵀ = B)
    (hComm : N (k := k) * B = B * (N (k := k))ᵀ) :
    B = upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) := by
  have h10 : B 1 0 = B 0 1 := by
    simpa [Matrix.transpose_apply] using congrArg (fun M => M 0 1) hSym
  have h20 : B 2 0 = B 0 2 := by
    simpa [Matrix.transpose_apply] using congrArg (fun M => M 0 2) hSym
  have h21 : B 2 1 = B 1 2 := by
    simpa [Matrix.transpose_apply] using congrArg (fun M => M 1 2) hSym
  have h11 : B 1 1 = B 0 2 := by
    simpa [N, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using
      congrArg (fun M => M 0 1) hComm
  have h12 : B 1 2 = 0 := by
    simpa [N, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using
      congrArg (fun M => M 0 2) hComm
  have h22 : B 2 2 = 0 := by
    simpa [N, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using
      congrArg (fun M => M 1 2) hComm
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [upperHankel, h10, h20, h21, h11, h12, h22]

/-- The lower Hankel matrix is symmetric. -/
theorem lowerHankel_transpose
    (u v w : k) :
    (lowerHankel (k := k) u v w)ᵀ = lowerHankel (k := k) u v w := by
  ext i j
  fin_cases i <;> fin_cases j <;> simp [lowerHankel, Matrix.transpose_apply]

/-- The lower Hankel matrix satisfies the commutation relation needed for `rep₂`. -/
theorem lowerHankel_comm
    (u v w : k) :
    lowerHankel (k := k) u v w * N (k := k) =
      (N (k := k))ᵀ * lowerHankel (k := k) u v w := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [N, lowerHankel, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three]

/-- A symmetric matrix satisfying the lower commutation relation has the expected lower
Hankel shape. -/
theorem eq_lowerHankel_of_transpose_eq_of_comm
    (C : Matrix I I k)
    (hSym : Cᵀ = C)
    (hComm : C * N (k := k) = (N (k := k))ᵀ * C) :
    C = lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) := by
  have h10 : C 1 0 = C 0 1 := by
    simpa [Matrix.transpose_apply] using congrArg (fun M => M 0 1) hSym
  have h20 : C 2 0 = C 0 2 := by
    simpa [Matrix.transpose_apply] using congrArg (fun M => M 0 2) hSym
  have h21 : C 2 1 = C 1 2 := by
    simpa [Matrix.transpose_apply] using congrArg (fun M => M 1 2) hSym
  have h00 : C 0 0 = 0 := by
    simpa [N, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using
      congrArg (fun M => M 0 1) hComm
  have h01 : C 0 1 = 0 := by
    simpa [N, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using
      congrArg (fun M => M 0 2) hComm
  have h11 : C 1 1 = C 0 2 := by
    simpa [N, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using
      congrArg (fun M => M 1 2) hComm
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [lowerHankel, h10, h20, h21, h00, h01, h11]

/-- The upper Hankel unipotent family fixes the pair pointwise. -/
theorem upperHankel_pointwise
    (u v w : k) :
    FixesPairBivector
      (Matrix.fromBlocks
        (1 : Matrix I I k)
        (upperHankel (k := k) u v w)
        0
        (1 : Matrix I I k)) := by
  constructor
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₁
          (Matrix.fromBlocks
            (1 : Matrix I I k)
            (upperHankel (k := k) u v w)
            0
            (1 : Matrix I I k)) =
        Matrix.fromBlocks
          (-(upperHankel (k := k) u v w) + (upperHankel (k := k) u v w)ᵀ)
          (1 : Matrix I I k)
          (-(1 : Matrix I I k))
          0 := by
            simpa [rep₁] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := upperHankel (k := k) u v w)
                (C := 0)
                (D := (1 : Matrix I I k))
                (M := (1 : Matrix I I k))
      _ = rep₁ := by
            rw [upperHankel_transpose]
            simp [rep₁]
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₂
          (Matrix.fromBlocks
            (1 : Matrix I I k)
            (upperHankel (k := k) u v w)
            0
            (1 : Matrix I I k)) =
        Matrix.fromBlocks
          (-(upperHankel (k := k) u v w) * (N (k := k))ᵀ +
            N (k := k) * (upperHankel (k := k) u v w)ᵀ)
          (N (k := k))
          (-((N (k := k))ᵀ))
          0 := by
            simpa [rep₂, N] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := upperHankel (k := k) u v w)
                (C := 0)
                (D := (1 : Matrix I I k))
                (M := N (k := k))
      _ = rep₂ := by
            rw [upperHankel_transpose, upperHankel_comm]
            simp [rep₂]

/-- On the upper unipotent cell, pointwise fixing is equivalent to the upper Hankel
shape. -/
theorem fixesPair_upperIdentity_iff
    (B : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) ↔
      B = upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) := by
  constructor
  · rintro ⟨h₁, h₂⟩
    have h₁' :
        Matrix.fromBlocks
          (-B + Bᵀ)
          (1 : Matrix I I k)
          (-(1 : Matrix I I k))
          0 =
        Matrix.fromBlocks
          0
          (1 : Matrix I I k)
          (-(1 : Matrix I I k))
          0 := by
      calc
        Matrix.fromBlocks
            (-B + Bᵀ)
            (1 : Matrix I I k)
            (-(1 : Matrix I I k))
            0 =
          ActBivector rep₁ (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) := by
            symm
            simpa [rep₁] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := B)
                (C := (0 : Matrix I I k))
                (D := (1 : Matrix I I k))
                (M := (1 : Matrix I I k))
        _ =
          Matrix.fromBlocks
            0
            (1 : Matrix I I k)
            (-(1 : Matrix I I k))
            0 := by
              simpa [rep₁] using h₁
    have h₂' :
        Matrix.fromBlocks
          (-B * (N (k := k))ᵀ + N (k := k) * Bᵀ)
          (N (k := k))
          (-((N (k := k))ᵀ))
          0 =
        Matrix.fromBlocks
          0
          (N (k := k))
          (-((N (k := k))ᵀ))
          0 := by
      calc
        Matrix.fromBlocks
            (-B * (N (k := k))ᵀ + N (k := k) * Bᵀ)
            (N (k := k))
            (-((N (k := k))ᵀ))
            0 =
          ActBivector rep₂ (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) := by
            symm
            simpa [rep₂, N] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := B)
                (C := (0 : Matrix I I k))
                (D := (1 : Matrix I I k))
                (M := N (k := k))
        _ =
          Matrix.fromBlocks
            0
            (N (k := k))
            (-((N (k := k))ᵀ))
            0 := by
              simpa [rep₂] using h₂
    have hSym0 := congrArg Matrix.toBlocks₁₁ h₁'
    have hComm0 := congrArg Matrix.toBlocks₁₁ h₂'
    have hSym : Bᵀ = B := by
      have hSym0' : -B + Bᵀ = 0 := by
        simpa [sub_eq_add_neg] using hSym0
      have h := congrArg (fun M => B + M) hSym0'
      simpa [add_assoc, add_comm, add_left_comm] using h
    have hComm : N (k := k) * B = B * (N (k := k))ᵀ := by
      rw [hSym] at hComm0
      have hComm0' : 0 = N (k := k) * B + -(B * (N (k := k))ᵀ) := by
        symm
        simpa [sub_eq_add_neg, add_assoc, add_comm, add_left_comm] using hComm0
      have h := congrArg (fun M => M + (B * (N (k := k))ᵀ)) hComm0'
      have h' : B * (N (k := k))ᵀ = N (k := k) * B := by
        simpa [add_assoc, add_comm, add_left_comm] using h
      exact h'.symm
    exact eq_upperHankel_of_transpose_eq_of_comm (k := k) B hSym hComm
  · intro hB
    rw [hB]
    exact upperHankel_pointwise (k := k) (B 0 0) (B 0 1) (B 0 2)

/-- The lower Hankel unipotent family fixes the pair pointwise. -/
theorem lowerHankel_pointwise
    (u v w : k) :
    FixesPairBivector
      (Matrix.fromBlocks
        (1 : Matrix I I k)
        0
        (lowerHankel (k := k) u v w)
        (1 : Matrix I I k)) := by
  constructor
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₁
          (Matrix.fromBlocks
            (1 : Matrix I I k)
            0
            (lowerHankel (k := k) u v w)
            (1 : Matrix I I k)) =
        Matrix.fromBlocks
          0
          (1 : Matrix I I k)
          (-(1 : Matrix I I k))
          (-(lowerHankel (k := k) u v w)ᵀ + lowerHankel (k := k) u v w) := by
            simpa [rep₁] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := (0 : Matrix I I k))
                (C := lowerHankel (k := k) u v w)
                (D := (1 : Matrix I I k))
                (M := (1 : Matrix I I k))
      _ = rep₁ := by
            rw [lowerHankel_transpose]
            simp [rep₁]
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₂
          (Matrix.fromBlocks
            (1 : Matrix I I k)
            0
            (lowerHankel (k := k) u v w)
            (1 : Matrix I I k)) =
        Matrix.fromBlocks
          0
          (N (k := k))
          (-((N (k := k))ᵀ))
          (-(N (k := k))ᵀ * (lowerHankel (k := k) u v w)ᵀ +
            lowerHankel (k := k) u v w * N (k := k)) := by
            simpa [rep₂, N] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := (0 : Matrix I I k))
                (C := lowerHankel (k := k) u v w)
                (D := (1 : Matrix I I k))
                (M := N (k := k))
      _ = rep₂ := by
            rw [lowerHankel_transpose, lowerHankel_comm]
            simp [rep₂]

/-- On the lower unipotent cell, pointwise fixing is equivalent to the lower Hankel
shape. -/
theorem fixesPair_lowerIdentity_iff
    (C : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) ↔
      C = lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) := by
  constructor
  · rintro ⟨h₁, h₂⟩
    have h₁' :
        Matrix.fromBlocks
          0
          (1 : Matrix I I k)
          (-(1 : Matrix I I k))
          (-Cᵀ + C) =
        Matrix.fromBlocks
          0
          (1 : Matrix I I k)
          (-(1 : Matrix I I k))
          0 := by
      calc
        Matrix.fromBlocks
            0
            (1 : Matrix I I k)
            (-(1 : Matrix I I k))
            (-Cᵀ + C) =
          ActBivector rep₁ (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) := by
            symm
            simpa [rep₁] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := (0 : Matrix I I k))
                (C := C)
                (D := (1 : Matrix I I k))
                (M := (1 : Matrix I I k))
        _ =
          Matrix.fromBlocks
            0
            (1 : Matrix I I k)
            (-(1 : Matrix I I k))
            0 := by
              simpa [rep₁] using h₁
    have h₂' :
        Matrix.fromBlocks
          0
          (N (k := k))
          (-((N (k := k))ᵀ))
          (-((N (k := k))ᵀ) * Cᵀ + C * N (k := k)) =
        Matrix.fromBlocks
          0
          (N (k := k))
          (-((N (k := k))ᵀ))
          0 := by
      calc
        Matrix.fromBlocks
            0
            (N (k := k))
            (-((N (k := k))ᵀ))
            (-((N (k := k))ᵀ) * Cᵀ + C * N (k := k)) =
          ActBivector rep₂ (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) := by
            symm
            simpa [rep₂, N] using
              act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k))
                (B := (0 : Matrix I I k))
                (C := C)
                (D := (1 : Matrix I I k))
                (M := N (k := k))
        _ =
          Matrix.fromBlocks
            0
            (N (k := k))
            (-((N (k := k))ᵀ))
            0 := by
              simpa [rep₂] using h₂
    have hSym0 := congrArg Matrix.toBlocks₂₂ h₁'
    have hComm0 := congrArg Matrix.toBlocks₂₂ h₂'
    have hSym : Cᵀ = C := by
      have hSym0' : -Cᵀ + C = 0 := by
        simpa [sub_eq_add_neg] using hSym0
      have h := congrArg (fun M => M + Cᵀ) hSym0'
      have h' : C = Cᵀ := by
        simpa [add_assoc, add_comm, add_left_comm] using h
      exact h'.symm
    have hComm : C * N (k := k) = (N (k := k))ᵀ * C := by
      rw [hSym] at hComm0
      have hComm0' : 0 = C * N (k := k) + -((N (k := k))ᵀ * C) := by
        symm
        simpa [sub_eq_add_neg, add_assoc, add_comm, add_left_comm] using hComm0
      have h := congrArg (fun M => M + ((N (k := k))ᵀ * C)) hComm0'
      have h' : (N (k := k))ᵀ * C = C * N (k := k) := by
        simpa [add_assoc, add_comm, add_left_comm] using h
      exact h'.symm
    exact eq_lowerHankel_of_transpose_eq_of_comm (k := k) C hSym hComm
  · intro hC
    rw [hC]
    exact lowerHankel_pointwise (k := k) (C 0 2) (C 1 2) (C 2 2)

/-- Combining the upper and lower Hankel families gives an explicit six-parameter
pointwise unipotent subgroup. -/
theorem upperLowerHankel_product_pointwise
    (u₁ v₁ w₁ u₂ v₂ w₂ : k) :
    let gUpper : Matrix V V k :=
      Matrix.fromBlocks
        (1 : Matrix I I k)
        (upperHankel (k := k) u₁ v₁ w₁)
        0
        (1 : Matrix I I k)
    let gLower : Matrix V V k :=
      Matrix.fromBlocks
        (1 : Matrix I I k)
        0
        (lowerHankel (k := k) u₂ v₂ w₂)
        (1 : Matrix I I k)
    FixesPairBivector (gUpper * gLower) := by
  intro gUpper gLower
  have hUpper := upperHankel_pointwise (k := k) u₁ v₁ w₁
  have hLower := lowerHankel_pointwise (k := k) u₂ v₂ w₂
  constructor <;> rw [FixesBivector]
  ·
    calc
      ActBivector rep₁ (gUpper * gLower) =
          ActBivector (ActBivector rep₁ gLower) gUpper := by
            rw [actBivector_mul]
      _ = ActBivector rep₁ gUpper := by rw [hLower.1]
      _ = rep₁ := by simpa [FixesBivector] using hUpper.1
  ·
    calc
      ActBivector rep₂ (gUpper * gLower) =
          ActBivector (ActBivector rep₂ gLower) gUpper := by
            rw [actBivector_mul]
      _ = ActBivector rep₂ gUpper := by rw [hLower.2]
      _ = rep₂ := by simpa [FixesBivector] using hUpper.2

/-- The Magma local kernel product cell for the indecomposable `3[a]` row is the
expected block matrix `[[1 + BC, B], [C, 1]]`. -/
theorem productCell_eq
    (B C : Matrix I I k) :
    Matrix.fromBlocks ((1 : Matrix I I k) + B * C) B C (1 : Matrix I I k) =
      (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) *
        (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    simp [Matrix.fromBlocks_multiply, Matrix.add_apply]

/-- On the Magma local kernel product cell, pointwise fixing is exactly the upper/lower
Hankel condition on the two kernel blocks. -/
theorem fixesPair_productCell_iff
    (B C : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks ((1 : Matrix I I k) + B * C) B C (1 : Matrix I I k)) ↔
      B = upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) ∧
        C = lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) := by
  let g : Matrix V V k := Matrix.fromBlocks ((1 : Matrix I I k) + B * C) B C (1 : Matrix I I k)
  let gUpper : Matrix V V k := Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)
  let gLower : Matrix V V k := Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)
  have hg : g = gUpper * gLower := by
    simpa [g, gUpper, gLower] using productCell_eq (k := k) B C
  constructor
  · intro hfix
    have hfixg : FixesPairBivector (gUpper * gLower) := by
      rwa [← hg]
    have hSym0 : -Cᵀ + C = 0 := by
      calc
        -Cᵀ + C = Matrix.toBlocks₂₂ (ActBivector rep₁ g) := by
          symm
          simpa [g, rep₁, Matrix.mul_add, Matrix.add_mul, Matrix.mul_assoc] using
            congrArg Matrix.toBlocks₂₂
              (act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k) + B * C)
                (B := B)
                (C := C)
                (D := (1 : Matrix I I k))
                (M := (1 : Matrix I I k)))
        _ = Matrix.toBlocks₂₂ rep₁ := by
          simpa [FixesBivector] using congrArg Matrix.toBlocks₂₂ hfix.1
        _ = 0 := by simp [rep₁]
    have hComm0 : -((N (k := k))ᵀ) * Cᵀ + C * N (k := k) = 0 := by
      calc
        -((N (k := k))ᵀ) * Cᵀ + C * N (k := k) =
            Matrix.toBlocks₂₂ (ActBivector rep₂ g) := by
          symm
          simpa [g, rep₂, N, Matrix.mul_add, Matrix.add_mul, Matrix.mul_assoc] using
            congrArg Matrix.toBlocks₂₂
              (act_fromBlocks
                (k := k)
                (A := (1 : Matrix I I k) + B * C)
                (B := B)
                (C := C)
                (D := (1 : Matrix I I k))
                (M := N (k := k)))
        _ = Matrix.toBlocks₂₂ rep₂ := by
          simpa [FixesBivector] using congrArg Matrix.toBlocks₂₂ hfix.2
        _ = 0 := by simp [rep₂]
    have hSym : Cᵀ = C := by
      have h := congrArg (fun M => M + Cᵀ) hSym0
      have h' : C = Cᵀ := by
        simpa [add_assoc, add_comm, add_left_comm] using h
      exact h'.symm
    have hComm : C * N (k := k) = (N (k := k))ᵀ * C := by
      rw [hSym] at hComm0
      have h := congrArg (fun M => M + ((N (k := k))ᵀ * C)) hComm0
      have h' : C * N (k := k) = (N (k := k))ᵀ * C := by
        simpa [add_assoc, add_comm, add_left_comm] using h
      exact h'
    have hC :
        C = lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) := by
      exact eq_lowerHankel_of_transpose_eq_of_comm (k := k) C hSym hComm
    have hLower : FixesPairBivector gLower := by
      dsimp [gLower]
      rw [hC]
      exact lowerHankel_pointwise (k := k) (C 0 2) (C 1 2) (C 2 2)
    have hUpper : FixesPairBivector gUpper := by
      constructor
      · calc
          ActBivector rep₁ gUpper =
              ActBivector (ActBivector rep₁ gLower) gUpper := by
                rw [hLower.1]
          _ = ActBivector rep₁ (gUpper * gLower) := by
                rw [actBivector_mul]
          _ = rep₁ := hfixg.1
      · calc
          ActBivector rep₂ gUpper =
              ActBivector (ActBivector rep₂ gLower) gUpper := by
                rw [hLower.2]
          _ = ActBivector rep₂ (gUpper * gLower) := by
                rw [actBivector_mul]
          _ = rep₂ := hfixg.2
    have hB :
        B = upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) := by
      exact (fixesPair_upperIdentity_iff (k := k) B).1 (by simpa [gUpper] using hUpper)
    exact ⟨hB, hC⟩
  · rintro ⟨hB, hC⟩
    have hfixg : FixesPairBivector (gUpper * gLower) := by
      dsimp [gUpper, gLower]
      rw [hB, hC]
      exact
        upperLowerHankel_product_pointwise
          (k := k)
          (u₁ := B 0 0)
          (v₁ := B 0 1)
          (w₁ := B 0 2)
          (u₂ := C 0 2)
          (v₂ := C 1 2)
          (w₂ := C 2 2)
    have hfixg' : FixesPairBivector g := by
      exact hg.symm ▸ hfixg
    simpa [g] using hfixg'

/-- The pointwise torus on the indecomposable length-three repeated-support orbit. -/
theorem torus_mul
    (a : k)
    (ha : a ≠ 0) :
    (a • (1 : Matrix I I k)) * ((a⁻¹ • (1 : Matrix I I k)))ᵀ = 1 := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [ha]

/-- The pointwise torus also fixes the second basis vector. -/
theorem torus_N_mul
    (a : k)
    (ha : a ≠ 0) :
    (a • (1 : Matrix I I k)) * N (k := k) * ((a⁻¹ • (1 : Matrix I I k)))ᵀ = N (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [N, ha, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] <;>
    field_simp [ha] <;>
    ring

/-- The obvious torus in the length-three repeated-support kernel fixes the pair
pointwise. -/
theorem scale_pointwise
    (a : k)
    (ha : a ≠ 0) :
    FixesPairBivector
      (Matrix.fromBlocks
        (a • (1 : Matrix I I k))
        0
        0
        (a⁻¹ • (1 : Matrix I I k))) := by
  constructor
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₁
          (Matrix.fromBlocks
            (a • (1 : Matrix I I k))
            0
            0
            (a⁻¹ • (1 : Matrix I I k))) =
        Matrix.fromBlocks
          0
          (((a • (1 : Matrix I I k)) * (1 : Matrix I I k) *
            ((a⁻¹ • (1 : Matrix I I k)))ᵀ))
          (-(((a⁻¹ • (1 : Matrix I I k)) * (1 : Matrix I I k) *
            ((a • (1 : Matrix I I k)))ᵀ)))
          0 := by
            simpa [rep₁] using
              act_fromBlocks
                (k := k)
                (A := a • (1 : Matrix I I k))
                (B := (0 : Matrix I I k))
                (C := (0 : Matrix I I k))
                (D := a⁻¹ • (1 : Matrix I I k))
                (M := (1 : Matrix I I k))
      _ = rep₁ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, ha]
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₂
          (Matrix.fromBlocks
            (a • (1 : Matrix I I k))
            0
            0
            (a⁻¹ • (1 : Matrix I I k))) =
        Matrix.fromBlocks
          0
          (((a • (1 : Matrix I I k)) * N (k := k) *
            ((a⁻¹ • (1 : Matrix I I k)))ᵀ))
          (-(((a⁻¹ • (1 : Matrix I I k)) * (N (k := k))ᵀ *
            ((a • (1 : Matrix I I k)))ᵀ)))
          0 := by
            simpa [rep₂, N] using
              act_fromBlocks
                (k := k)
                (A := a • (1 : Matrix I I k))
                (B := (0 : Matrix I I k))
                (C := (0 : Matrix I I k))
                (D := a⁻¹ • (1 : Matrix I I k))
                (M := N (k := k))
      _ = rep₂ := by
            simp [Matrix.mul_assoc, torus_N_mul (k := k) a ha, rep₂, ha]

/-- On the upper Hankel cell followed by the pointwise torus, pointwise fixing is still
equivalent to the upper Hankel shape. -/
theorem fixesPair_upperScale_iff
    (B : Matrix I I k)
    (a : k)
    (ha : a ≠ 0) :
    FixesPairBivector
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          B
          0
          (1 : Matrix I I k)) *
       (Matrix.fromBlocks
          (a • (1 : Matrix I I k))
          0
          0
          (a⁻¹ • (1 : Matrix I I k)))) ↔
      B = upperHankel (k := k) (B 0 0) (B 0 1) (B 0 2) := by
  have hScale := scale_pointwise (k := k) a ha
  constructor
  · intro h
    have hUpper : FixesPairBivector
        (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) := by
      constructor
      · have htmp : ActBivector rep₁
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                B
                0
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) = rep₁ := h.1
        rw [actBivector_mul] at htmp
        rw [hScale.1] at htmp
        exact htmp
      · have htmp : ActBivector rep₂
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                B
                0
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) = rep₂ := h.2
        rw [actBivector_mul] at htmp
        rw [hScale.2] at htmp
        exact htmp
    exact (fixesPair_upperIdentity_iff (k := k) B).1 hUpper
  · intro hB
    have hUpper : FixesPairBivector
        (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) :=
      (fixesPair_upperIdentity_iff (k := k) B).2 hB
    constructor
    · calc
        ActBivector rep₁
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                B
                0
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) =
          ActBivector
            (ActBivector rep₁
              (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k))))
            (Matrix.fromBlocks
              (1 : Matrix I I k)
              B
              0
              (1 : Matrix I I k)) := by
                rw [actBivector_mul]
      _ = ActBivector rep₁ (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) := by
            simpa using congrArg
              (fun M => ActBivector M (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)))
              hScale.1
      _ = rep₁ := hUpper.1
    · calc
        ActBivector rep₂
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                B
                0
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) =
          ActBivector
            (ActBivector rep₂
              (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k))))
            (Matrix.fromBlocks
              (1 : Matrix I I k)
              B
              0
              (1 : Matrix I I k)) := by
                rw [actBivector_mul]
      _ = ActBivector rep₂ (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)) := by
            simpa using congrArg
              (fun M => ActBivector M (Matrix.fromBlocks (1 : Matrix I I k) B 0 (1 : Matrix I I k)))
              hScale.2
      _ = rep₂ := hUpper.2

/-- On the lower Hankel cell followed by the pointwise torus, pointwise fixing is still
equivalent to the lower Hankel shape. -/
theorem fixesPair_lowerScale_iff
    (C : Matrix I I k)
    (a : k)
    (ha : a ≠ 0) :
    FixesPairBivector
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          C
          (1 : Matrix I I k)) *
       (Matrix.fromBlocks
          (a • (1 : Matrix I I k))
          0
          0
          (a⁻¹ • (1 : Matrix I I k)))) ↔
      C = lowerHankel (k := k) (C 0 2) (C 1 2) (C 2 2) := by
  have hScale := scale_pointwise (k := k) a ha
  constructor
  · intro h
    have hLower : FixesPairBivector
        (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) := by
      constructor
      · have htmp : ActBivector rep₁
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                0
                C
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) = rep₁ := h.1
        rw [actBivector_mul] at htmp
        rw [hScale.1] at htmp
        exact htmp
      · have htmp : ActBivector rep₂
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                0
                C
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) = rep₂ := h.2
        rw [actBivector_mul] at htmp
        rw [hScale.2] at htmp
        exact htmp
    exact (fixesPair_lowerIdentity_iff (k := k) C).1 hLower
  · intro hC
    have hLower : FixesPairBivector
        (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) :=
      (fixesPair_lowerIdentity_iff (k := k) C).2 hC
    constructor
    · calc
        ActBivector rep₁
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                0
                C
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) =
          ActBivector
            (ActBivector rep₁
              (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k))))
            (Matrix.fromBlocks
              (1 : Matrix I I k)
              0
              C
              (1 : Matrix I I k)) := by
                rw [actBivector_mul]
      _ = ActBivector rep₁ (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) := by
            simpa using congrArg
              (fun M => ActBivector M (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)))
              hScale.1
      _ = rep₁ := hLower.1
    · calc
        ActBivector rep₂
            ((Matrix.fromBlocks
                (1 : Matrix I I k)
                0
                C
                (1 : Matrix I I k)) *
             (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k)))) =
          ActBivector
            (ActBivector rep₂
              (Matrix.fromBlocks
                (a • (1 : Matrix I I k))
                0
                0
                (a⁻¹ • (1 : Matrix I I k))))
            (Matrix.fromBlocks
              (1 : Matrix I I k)
              0
              C
              (1 : Matrix I I k)) := by
                rw [actBivector_mul]
      _ = ActBivector rep₂ (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)) := by
            simpa using congrArg
              (fun M => ActBivector M (Matrix.fromBlocks (1 : Matrix I I k) 0 C (1 : Matrix I I k)))
              hScale.2
      _ = rep₂ := hLower.2

/-- Combining the six-parameter unipotent family with the torus gives a concrete
pointwise subgroup on the indecomposable length-three repeated-support orbit. -/
theorem scaleUpperLowerHankel_product_pointwise
    (u₁ v₁ w₁ u₂ v₂ w₂ a : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      (Matrix.fromBlocks
        (1 : Matrix I I k)
        (upperHankel (k := k) u₁ v₁ w₁)
        0
        (1 : Matrix I I k)) *
      (Matrix.fromBlocks
        (1 : Matrix I I k)
        0
        (lowerHankel (k := k) u₂ v₂ w₂)
        (1 : Matrix I I k))
    let gScale : Matrix V V k :=
      Matrix.fromBlocks
        (a • (1 : Matrix I I k))
        0
        0
        (a⁻¹ • (1 : Matrix I I k))
    FixesPairBivector (gU * gScale) := by
  intro gU gScale
  have hU := upperLowerHankel_product_pointwise
    (k := k) (u₁ := u₁) (v₁ := v₁) (w₁ := w₁) (u₂ := u₂) (v₂ := v₂) (w₂ := w₂)
  have hScale := scale_pointwise (k := k) a ha
  constructor <;> rw [FixesBivector]
  ·
    calc
      ActBivector rep₁ (gU * gScale) = ActBivector (ActBivector rep₁ gScale) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gU := by
        simpa [FixesBivector] using congrArg (fun M => ActBivector M gU) hScale.1
      _ = rep₁ := by simpa [FixesBivector] using hU.1
  ·
    calc
      ActBivector rep₂ (gU * gScale) = ActBivector (ActBivector rep₂ gScale) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by
        simpa [FixesBivector] using congrArg (fun M => ActBivector M gU) hScale.2
      _ = rep₂ := by simpa [FixesBivector] using hU.2

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit has determinant `1`. -/
theorem scaleUpperLowerHankel_product_det
    (u₁ v₁ w₁ u₂ v₂ w₂ a : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      (Matrix.fromBlocks
        (1 : Matrix I I k)
        (upperHankel (k := k) u₁ v₁ w₁)
        0
        (1 : Matrix I I k)) *
      (Matrix.fromBlocks
        (1 : Matrix I I k)
        0
        (lowerHankel (k := k) u₂ v₂ w₂)
        (1 : Matrix I I k))
    let gScale : Matrix V V k :=
      Matrix.fromBlocks
        (a • (1 : Matrix I I k))
        0
        0
        (a⁻¹ • (1 : Matrix I I k))
    Matrix.det (gU * gScale) = 1 := by
  intro gU gScale
  have hUpper :
      Matrix.det
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) = 1 := by
    rw [Matrix.det_fromBlocks_zero₂₁]
    simp
  have hLower :
      Matrix.det
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k)) = 1 := by
    rw [Matrix.det_fromBlocks_zero₁₂]
    simp
  have hcard : Fintype.card I = 3 := by
    simp [I]
  have hScale : Matrix.det gScale = 1 := by
    simp [gScale, Matrix.det_fromBlocks_zero₂₁, Matrix.det_smul, hcard, ha]
  rw [Matrix.det_mul, Matrix.det_mul, hUpper, hLower, hScale]
  ring

/-- The quotient shear blocks produce the expected upper-triangular action. -/
theorem shear_mul
    (b : k) :
    shearTop (k := k) b * (shearBottom (k := k) b)ᵀ =
      (1 : Matrix I I k) + b • (N (k := k)) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [shearTop, shearBottom, N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] <;>
    ring

/-- The quotient shear keeps the second basis vector fixed. -/
theorem shear_N_mul
    (b : k) :
    shearTop (k := k) b * N (k := k) * (shearBottom (k := k) b)ᵀ = N (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [shearTop, shearBottom, N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] <;>
    ring

/-- The lower block in the quotient shear lift matches the transposed upper block. -/
theorem shear_bottom_mul
    (b : k) :
    shearBottom (k := k) b * (shearTop (k := k) b)ᵀ =
      ((1 : Matrix I I k) + b • (N (k := k)))ᵀ := by
  simpa [Matrix.transpose_mul, Matrix.transpose_add, Matrix.transpose_smul] using
    congrArg Matrix.transpose (shear_mul (k := k) b)

/-- The lower block in the quotient shear lift also fixes the transposed second basis
block. -/
theorem shear_bottom_N_mul
    (b : k) :
    shearBottom (k := k) b * (N (k := k))ᵀ * (shearTop (k := k) b)ᵀ = (N (k := k))ᵀ := by
  simpa [Matrix.transpose_mul, Matrix.mul_assoc] using
    congrArg Matrix.transpose (shear_N_mul (k := k) b)

/-- The quotient scaling blocks act on the first basis vector by scalar multiplication. -/
theorem scale_mul
    (a : k)
    (ha : a ≠ 0) :
    scaleTop (k := k) a * (scaleBottom (k := k) a)ᵀ = a • (1 : Matrix I I k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [scaleTop, scaleBottom, Matrix.smul_apply, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] <;>
    field_simp [ha] <;>
    ring

/-- The quotient scaling blocks keep the second basis vector fixed. -/
theorem scale_N_mul
    (a : k)
    (ha : a ≠ 0) :
    scaleTop (k := k) a * N (k := k) * (scaleBottom (k := k) a)ᵀ = N (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [scaleTop, scaleBottom, N, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three] <;>
    field_simp [ha] <;>
    ring

/-- The lower block in the quotient scaling lift matches the same scalar on the
transposed first basis block. -/
theorem scale_bottom_mul
    (a : k)
    (ha : a ≠ 0) :
    scaleBottom (k := k) a * (scaleTop (k := k) a)ᵀ = (a • (1 : Matrix I I k))ᵀ := by
  simpa [Matrix.transpose_mul, Matrix.transpose_smul] using
    congrArg Matrix.transpose (scale_mul (k := k) a ha)

/-- The lower block in the quotient scaling lift also fixes the transposed second basis
block. -/
theorem scale_bottom_N_mul
    (a : k)
    (ha : a ≠ 0) :
    scaleBottom (k := k) a * (N (k := k))ᵀ * (scaleTop (k := k) a)ᵀ = (N (k := k))ᵀ := by
  simpa [Matrix.transpose_mul, Matrix.mul_assoc] using
    congrArg Matrix.transpose (scale_N_mul (k := k) a ha)

/-- The quotient shear lift on the length-three repeated-support orbit. -/
theorem shear_lift_action
    (b : k) :
    let g : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    ActBivector rep₁ g = rep₁ + b • rep₂ ∧
      ActBivector rep₂ g = rep₂ := by
  intro g
  constructor
  ·
    calc
      ActBivector rep₁ g =
        Matrix.fromBlocks
          0
          (shearTop (k := k) b * (shearBottom (k := k) b)ᵀ)
          (-(shearBottom (k := k) b * (shearTop (k := k) b)ᵀ))
          0 := by
            simpa [g, rep₁] using
              act_fromBlocks
                (k := k)
                (A := shearTop (k := k) b)
                (B := (0 : Matrix I I k))
                (C := (0 : Matrix I I k))
                (D := shearBottom (k := k) b)
                (M := (1 : Matrix I I k))
      _ = Matrix.fromBlocks
            0
            ((1 : Matrix I I k) + b • (N (k := k)))
            (-(((1 : Matrix I I k) + b • (N (k := k)))ᵀ))
            0 := by
            rw [shear_mul (k := k) b, shear_bottom_mul (k := k) b]
      _ = rep₁ + b • rep₂ := by
            simpa [Matrix.transpose_add, Matrix.transpose_smul, one_smul, add_assoc, add_left_comm,
              add_comm] using
              (linearCombination_fromBlocks (k := k) (a := (1 : k)) (b := b))
  ·
    calc
      ActBivector rep₂ g =
        Matrix.fromBlocks
          0
          (shearTop (k := k) b * N (k := k) * (shearBottom (k := k) b)ᵀ)
          (-(shearBottom (k := k) b * (N (k := k))ᵀ * (shearTop (k := k) b)ᵀ))
          0 := by
            simpa [g, rep₂, N] using
              act_fromBlocks
                (k := k)
                (A := shearTop (k := k) b)
                (B := (0 : Matrix I I k))
                (C := (0 : Matrix I I k))
                (D := shearBottom (k := k) b)
                (M := N (k := k))
      _ = rep₂ := by
            rw [shear_N_mul (k := k) b, shear_bottom_N_mul (k := k) b]
            simp [rep₂]

/-- The quotient shear lift on the indecomposable `3[a]` orbit has determinant `1`. -/
theorem shear_lift_det
    (b : k) :
    let g : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    Matrix.det g = 1 := by
  intro g
  rw [Matrix.det_fromBlocks_zero₂₁]
  simp [g, shearTop_det, shearBottom_det]

/-- The quotient scaling lift on the length-three repeated-support orbit. -/
theorem scale_lift_action
    (a : k)
    (ha : a ≠ 0) :
    let g : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    ActBivector rep₁ g = a • rep₁ ∧
      ActBivector rep₂ g = rep₂ := by
  intro g
  constructor
  ·
    calc
      ActBivector rep₁ g =
        Matrix.fromBlocks
          0
          (scaleTop (k := k) a * (scaleBottom (k := k) a)ᵀ)
          (-(scaleBottom (k := k) a * (scaleTop (k := k) a)ᵀ))
          0 := by
            simpa [g, rep₁] using
              act_fromBlocks
                (k := k)
                (A := scaleTop (k := k) a)
                (B := (0 : Matrix I I k))
                (C := (0 : Matrix I I k))
                (D := scaleBottom (k := k) a)
                (M := (1 : Matrix I I k))
      _ = Matrix.fromBlocks
            0
            (a • (1 : Matrix I I k))
            (-((a • (1 : Matrix I I k))ᵀ))
            0 := by
            rw [scale_mul (k := k) a ha, scale_bottom_mul (k := k) a ha]
      _ = a • rep₁ := by
            ext i j
            cases i <;> cases j <;>
              simp [rep₁, Matrix.fromBlocks, Matrix.smul_apply, Matrix.transpose_smul]
  ·
    calc
      ActBivector rep₂ g =
        Matrix.fromBlocks
          0
          (scaleTop (k := k) a * N (k := k) * (scaleBottom (k := k) a)ᵀ)
          (-(scaleBottom (k := k) a * (N (k := k))ᵀ * (scaleTop (k := k) a)ᵀ))
          0 := by
            simpa [g, rep₂, N] using
              act_fromBlocks
                (k := k)
                (A := scaleTop (k := k) a)
                (B := (0 : Matrix I I k))
                (C := (0 : Matrix I I k))
                (D := scaleBottom (k := k) a)
                (M := N (k := k))
      _ = rep₂ := by
            rw [scale_N_mul (k := k) a ha, scale_bottom_N_mul (k := k) a ha]
            simp [rep₂]

/-- The quotient scaling lift on the indecomposable `3[a]` orbit has determinant `a^3`. -/
theorem scale_lift_det
    (a : k)
    (ha : a ≠ 0) :
    let g : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    Matrix.det g = a * a * a := by
  intro g
  rw [Matrix.det_fromBlocks_zero₂₁]
  simp [g, scaleTop_det, scaleBottom_det, ha]

/-- The combined upper-triangular quotient lift on the indecomposable `3[a]` orbit. -/
theorem borel_lift_action
    (a b : k)
    (ha : a ≠ 0) :
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    let gU : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    ActBivector rep₁ (gS * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gS * gU) = rep₂ := by
  intro gS gU
  have hS := scale_lift_action (k := k) (a := a) ha
  have hU := shear_lift_action (k := k) (b := b)
  constructor
  · calc
      ActBivector rep₁ (gS * gU) = ActBivector (ActBivector rep₁ gU) gS := by
        rw [actBivector_mul]
      _ = ActBivector (rep₁ + b • rep₂) gS := by
        simpa [gU] using congrArg (fun M => ActBivector M gS) hU.1
      _ = ActBivector rep₁ gS + b • ActBivector rep₂ gS := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ + b • rep₂ := by
        rw [hS.1, hS.2]
  · calc
      ActBivector rep₂ (gS * gU) = ActBivector (ActBivector rep₂ gU) gS := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gS := by
        simpa [gU] using congrArg (fun M => ActBivector M gS) hU.2
      _ = rep₂ := by
        simpa [gS] using hS.2

/-- The standard Borel lift on the indecomposable `3[a]` orbit has determinant `a^3`. -/
theorem borel_lift_det
    (a b : k)
    (ha : a ≠ 0) :
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    let gU : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    Matrix.det (gS * gU) = a * a * a := by
  intro gS gU
  rw [Matrix.det_mul]
  simp [gS, gU, Matrix.det_mul, Matrix.det_fromBlocks_zero₂₁,
    scaleTop_det, scaleBottom_det, shearTop_det, shearBottom_det, ha]

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit combines with the
quotient shear lift exactly as expected. -/
theorem pointwise_shear_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    ActBivector rep₁ (gU * gL) = rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = rep₂ := by
  intro gU gL
  have hU :=
    scaleUpperLowerHankel_product_pointwise
      (k := k)
      (u₁ := u₁)
      (v₁ := v₁)
      (w₁ := w₁)
      (u₂ := u₂)
      (v₂ := v₂)
      (w₂ := w₂)
      (a := c)
      hc
  have hL := shear_lift_action (k := k) (b := b)
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (rep₁ + b • rep₂) gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.1
      _ = ActBivector rep₁ gU + b • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = rep₁ + b • rep₂ := by
        rw [hU.1, hU.2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.2
      _ = rep₂ := hU.2

/-- Left-multiplying the quotient shear lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem pointwise_shear_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    Matrix.det (gU * gL) = 1 := by
  intro gU gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using
      scaleUpperLowerHankel_product_det
        (k := k)
        (u₁ := u₁)
        (v₁ := v₁)
        (w₁ := w₁)
        (u₂ := u₂)
        (v₂ := v₂)
        (w₂ := w₂)
        (a := c)
        hc
  have hLdet : Matrix.det gL = 1 := by
    simpa [gL] using shear_lift_det (k := k) (b := b)
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit also combines with the
quotient scaling lift exactly as expected. -/
theorem pointwise_scale_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    ActBivector rep₁ (gU * gL) = a • rep₁ ∧
      ActBivector rep₂ (gU * gL) = rep₂ := by
  intro gU gL
  have hU :=
    scaleUpperLowerHankel_product_pointwise
      (k := k)
      (u₁ := u₁)
      (v₁ := v₁)
      (w₁ := w₁)
      (u₂ := u₂)
      (v₂ := v₂)
      (w₂ := w₂)
      (a := c)
      hc
  have hL := scale_lift_action (k := k) (a := a) ha
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • rep₁) gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector rep₁ gU := by
        simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ := by
        rw [hU.1]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.2
      _ = rep₂ := hU.2

/-- Left-multiplying the quotient scaling lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem pointwise_scale_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    Matrix.det (gU * gL) = a * a * a := by
  intro gU gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using
      scaleUpperLowerHankel_product_det
        (k := k)
        (u₁ := u₁)
        (v₁ := v₁)
        (w₁ := w₁)
        (u₂ := u₂)
        (v₂ := v₂)
        (w₂ := w₂)
        (a := c)
        hc
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using scale_lift_det (k := k) (a := a) ha
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- Right-multiplying the quotient shear lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change the quotient action. -/
theorem shear_pointwise_right_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    ActBivector rep₁ (gL * gU) = rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = rep₂ := by
  intro gL gU
  have hL := shear_lift_action (k := k) (b := b)
  have hU :=
    scaleUpperLowerHankel_product_pointwise
      (k := k)
      (u₁ := u₁)
      (v₁ := v₁)
      (w₁ := w₁)
      (u₂ := u₂)
      (v₂ := v₂)
      (w₂ := w₂)
      (a := c)
      hc
  constructor
  · calc
      ActBivector rep₁ (gL * gU) = ActBivector (ActBivector rep₁ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.1
      _ = rep₁ + b • rep₂ := by simpa [gL] using hL.1
  · calc
      ActBivector rep₂ (gL * gU) = ActBivector (ActBivector rep₂ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.2
      _ = rep₂ := by simpa [gL] using hL.2

/-- Right-multiplying the quotient shear lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem shear_pointwise_right_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c b : k)
    (hc : c ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    Matrix.det (gL * gU) = 1 := by
  intro gL gU
  have hLdet : Matrix.det gL = 1 := by
    simpa [gL] using shear_lift_det (k := k) (b := b)
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using
      scaleUpperLowerHankel_product_det
        (k := k)
        (u₁ := u₁)
        (v₁ := v₁)
        (w₁ := w₁)
        (u₂ := u₂)
        (v₂ := v₂)
        (w₂ := w₂)
        (a := c)
        hc
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

/-- The full pointwise subgroup on the indecomposable `3[a]` orbit also combines with
the full upper-triangular quotient lift exactly as expected. -/
theorem pointwise_borel_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    let gL : Matrix V V k :=
      (Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b))
    ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = rep₂ := by
  intro gU gL
  have hU :=
    scaleUpperLowerHankel_product_pointwise
      (k := k)
      (u₁ := u₁)
      (v₁ := v₁)
      (w₁ := w₁)
      (u₂ := u₂)
      (v₂ := v₂)
      (w₂ := w₂)
      (a := c)
      hc
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • rep₁ + b • rep₂) gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector rep₁ gU + b • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ + b • rep₂ := by
        rw [hU.1, hU.2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.2
      _ = rep₂ := hU.2

/-- Left-multiplying the standard Borel lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem pointwise_borel_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    let gL : Matrix V V k :=
      (Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b))
    Matrix.det (gU * gL) = a * a * a := by
  intro gU gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using
      scaleUpperLowerHankel_product_det
        (k := k)
        (u₁ := u₁)
        (v₁ := v₁)
        (w₁ := w₁)
        (u₂ := u₂)
        (v₂ := v₂)
        (w₂ := w₂)
        (a := c)
        hc
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using borel_lift_det (k := k) (a := a) (b := b) ha
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- Right-multiplying the quotient scaling lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change the quotient action. -/
theorem scale_pointwise_right_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    ActBivector rep₁ (gL * gU) = a • rep₁ ∧
      ActBivector rep₂ (gL * gU) = rep₂ := by
  intro gL gU
  have hL := scale_lift_action (k := k) (a := a) ha
  have hU :=
    scaleUpperLowerHankel_product_pointwise
      (k := k)
      (u₁ := u₁)
      (v₁ := v₁)
      (w₁ := w₁)
      (u₂ := u₂)
      (v₂ := v₂)
      (w₂ := w₂)
      (a := c)
      hc
  constructor
  · calc
      ActBivector rep₁ (gL * gU) = ActBivector (ActBivector rep₁ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.1
      _ = a • rep₁ := by simpa [gL] using hL.1
  · calc
      ActBivector rep₂ (gL * gU) = ActBivector (ActBivector rep₂ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.2
      _ = rep₂ := by simpa [gL] using hL.2

/-- Right-multiplying the quotient scaling lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem scale_pointwise_right_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    Matrix.det (gL * gU) = a * a * a := by
  intro gL gU
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using scale_lift_det (k := k) (a := a) ha
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using
      scaleUpperLowerHankel_product_det
        (k := k)
        (u₁ := u₁)
        (v₁ := v₁)
        (w₁ := w₁)
        (u₂ := u₂)
        (v₂ := v₂)
        (w₂ := w₂)
        (a := c)
        hc
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

/-- Right-multiplying the full upper-triangular quotient lift by the full pointwise
subgroup on the indecomposable `3[a]` orbit does not change the quotient action. -/
theorem borel_pointwise_right_product_lift_action
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      (Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b))
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = rep₂ := by
  intro gL gU
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  have hU :=
    scaleUpperLowerHankel_product_pointwise
      (k := k)
      (u₁ := u₁)
      (v₁ := v₁)
      (w₁ := w₁)
      (u₂ := u₂)
      (v₂ := v₂)
      (w₂ := w₂)
      (a := c)
      hc
  constructor
  · calc
      ActBivector rep₁ (gL * gU) = ActBivector (ActBivector rep₁ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.1
      _ = a • rep₁ + b • rep₂ := by
        simpa [gL] using hL.1
  · calc
      ActBivector rep₂ (gL * gU) = ActBivector (ActBivector rep₂ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.2
      _ = rep₂ := by
        simpa [gL] using hL.2

/-- Right-multiplying the standard Borel lift by the full pointwise subgroup on the
indecomposable `3[a]` orbit does not change its determinant. -/
theorem borel_pointwise_right_product_lift_det
    (u₁ v₁ w₁ u₂ v₂ w₂ c a b : k)
    (hc : c ≠ 0)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      (Matrix.fromBlocks
        (scaleTop (k := k) a)
        0
        0
        (scaleBottom (k := k) a)) *
      (Matrix.fromBlocks
        (shearTop (k := k) b)
        0
        0
        (shearBottom (k := k) b))
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks
          (1 : Matrix I I k)
          (upperHankel (k := k) u₁ v₁ w₁)
          0
          (1 : Matrix I I k)) *
        (Matrix.fromBlocks
          (1 : Matrix I I k)
          0
          (lowerHankel (k := k) u₂ v₂ w₂)
          (1 : Matrix I I k))) *
      (Matrix.fromBlocks
        (c • (1 : Matrix I I k))
        0
        0
        (c⁻¹ • (1 : Matrix I I k)))
    Matrix.det (gL * gU) = a * a * a := by
  intro gL gU
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using borel_lift_det (k := k) (a := a) (b := b) ha
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using
      scaleUpperLowerHankel_product_det
        (k := k)
        (u₁ := u₁)
        (v₁ := v₁)
        (w₁ := w₁)
        (u₂ := u₂)
        (v₂ := v₂)
        (w₂ := w₂)
        (a := c)
        hc
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

end N6OnePointLong
end Wedge2Formalization
