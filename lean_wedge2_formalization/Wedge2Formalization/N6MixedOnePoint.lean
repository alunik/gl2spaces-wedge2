import Wedge2Formalization.N4Summary
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.SchurComplement

open Matrix

namespace Wedge2Formalization
namespace N6MixedOnePoint

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

/-- The bivector action by a matrix with invertible determinant is injective. -/
theorem actBivector_injective_of_isUnit_det
    (g : Matrix V V k)
    (hg : IsUnit (Matrix.det g)) :
    Function.Injective (fun Ω : Matrix V V k => ActBivector Ω g) := by
  letI := Matrix.invertibleOfIsUnitDet g hg
  have hgT : IsUnit (Matrix.det gᵀ) := by
    simpa [Matrix.det_transpose] using hg
  letI := Matrix.invertibleOfIsUnitDet gᵀ hgT
  intro Ω₁ Ω₂ hΩ
  have h1 : Ω₁ * gᵀ = Ω₂ * gᵀ := by
    apply (Matrix.mul_right_inj_of_invertible (A := g)).mp
    simpa [ActBivector, Matrix.mul_assoc] using hΩ
  exact (Matrix.mul_left_inj_of_invertible (A := gᵀ)).mp h1

/-- A convenient `3 × 3` block matrix with `2 × 2` blocks. -/
def Block3
    (A11 A12 A13 A21 A22 A23 A31 A32 A33 : Matrix I I k) : Matrix V V k :=
  fun i j =>
    match i, j with
    | Sum.inl (Sum.inl i), Sum.inl (Sum.inl j) => A11 i j
    | Sum.inl (Sum.inl i), Sum.inl (Sum.inr j) => A12 i j
    | Sum.inl (Sum.inl i), Sum.inr j => A13 i j
    | Sum.inl (Sum.inr i), Sum.inl (Sum.inl j) => A21 i j
    | Sum.inl (Sum.inr i), Sum.inl (Sum.inr j) => A22 i j
    | Sum.inl (Sum.inr i), Sum.inr j => A23 i j
    | Sum.inr i, Sum.inl (Sum.inl j) => A31 i j
    | Sum.inr i, Sum.inl (Sum.inr j) => A32 i j
    | Sum.inr i, Sum.inr j => A33 i j

/-- The elementary matrix `E_{12}`. -/
def E12 : Matrix I I k := !![(0 : k), 1; 0, 0]

/-- The elementary matrix `E_{11}`. -/
def E11 : Matrix I I k := !![(1 : k), 0; 0, 0]

/-- The elementary matrix `E_{21}`. -/
def E21 : Matrix I I k := !![(0 : k), 0; 1, 0]

/-- The elementary matrix `E_{22}`. -/
def E22 : Matrix I I k := !![(0 : k), 0; 0, 1]

/-- A first coupled unipotent generator for the mixed one-point kernel. -/
def crossUnipotent₁ : Matrix V V k :=
  Block3
    (1 : Matrix I I k) 0 0
    0 (1 : Matrix I I k) (E21 (k := k))
    0 (E22 (k := k)) (1 : Matrix I I k)

/-- A second coupled unipotent generator for the mixed one-point kernel. -/
def crossUnipotent₂ : Matrix V V k :=
  Block3
    (1 : Matrix I I k) 0 0
    0 (1 : Matrix I I k) (E22 (k := k))
    (-E12 (k := k)) 0 (1 : Matrix I I k)

/-- The upper-right block of the first coupled unipotent generator. -/
def crossQ₁ : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl _ => 0
    | Sum.inr i => E21 (k := k) i j

/-- The lower-left block of the first coupled unipotent generator. -/
def crossR₁ : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => E22 (k := k) i j

/-- The upper-right block of the second coupled unipotent generator. -/
def crossQ₂ : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl _ => 0
    | Sum.inr i => E22 (k := k) i j

/-- The lower-left block of the second coupled unipotent generator. -/
def crossR₂ : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl j => (-E12 (k := k)) i j
    | Sum.inr _ => 0

/-- A representative for the mixed repeated-support orbit `J_{a,2} + J_{a,1}`. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks (N4.onePointRep₁ (k := k)) 0 0 (N4.J (k := k))

/-- The second basis vector of the mixed repeated-support representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks (N4.onePointRep₂ (k := k)) 0 0 0

/-- The `4 × 4` change of basis carrying the Lean `n = 4` repeated-support model to the
native Magma local model. -/
def magmaTopChange : Matrix W W k :=
  Matrix.fromBlocks
    (E12 (k := k))
    (E22 (k := k))
    (-E21 (k := k))
    (-E11 (k := k))

/-- The top repeated-support representative used in the native Magma local model. -/
def magmaTopRep₁ : Matrix W W k :=
  Matrix.fromBlocks 0 (1 : Matrix I I k) (- (1 : Matrix I I k)) 0

/-- The second top basis vector used in the native Magma local model. -/
def magmaTopRep₂ : Matrix W W k :=
  Matrix.fromBlocks 0 (E12 (k := k)) (-E21 (k := k)) 0

/-- The ambient `6 × 6` basis change carrying the Lean mixed representative to the native
Magma local model. -/
def magmaChange : Matrix V V k :=
  Matrix.fromBlocks (magmaTopChange (k := k)) 0 0 (1 : Matrix I I k)

/-- The mixed repeated-support representative in the native Magma local model. -/
def magmaRep₁ : Matrix V V k :=
  Matrix.fromBlocks (magmaTopRep₁ (k := k)) 0 0 (N4.J (k := k))

/-- The second mixed basis vector in the native Magma local model. -/
def magmaRep₂ : Matrix V V k :=
  Matrix.fromBlocks (magmaTopRep₂ (k := k)) 0 0 0

/-- Exact stabilizer of the mixed one-point pair for the action on `∧²V`. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  ActBivector (rep₁ (k := k)) g = rep₁ ∧
    ActBivector (rep₂ (k := k)) g = rep₂

/-- In the chosen basis, every vector on the mixed one-point line is block diagonal. -/
theorem basis_linearCombination
    (a b : k) :
    a • (rep₁ (k := k)) + b • (rep₂ (k := k)) =
      Matrix.fromBlocks
        (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k)))
        0
        0
        (a • (N4.J (k := k))) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₁, rep₂, Matrix.add_apply]

/-- Re-embedding the `n = 4` one-point line together with the compatible lower-right
simple block recovers the chosen basis on the mixed one-point line. -/
theorem embed_linearCombination
    (a b : k) :
    Matrix.fromBlocks
      (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k)))
      0
      0
      (a • (N4.J (k := k))) =
    a • (rep₁ (k := k)) + b • (rep₂ (k := k) ) := by
  simpa using (basis_linearCombination (k := k) (a := a) (b := b)).symm

/-- Re-embedding the second basis vector on the top block recovers the second basis
vector on the mixed one-point line. -/
theorem embed_second
    (c : k) :
    Matrix.fromBlocks (c • (N4.onePointRep₂ (k := k))) 0 0 (0 : Matrix I I k) =
      c • (rep₂ (k := k)) := by
  ext i j
  cases i <;> cases j <;>
    simp [rep₂, Matrix.fromBlocks, Matrix.smul_apply]

/-- Acting by a block diagonal matrix keeps the mixed one-point representative block
diagonal. -/
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

/-- Acting by an identity-plus-off-diagonal block matrix on a block diagonal bivector
produces the expected quadratic correction terms. -/
theorem act_identityOffDiagonal
    (ΩW : Matrix W W k)
    (ΩI : Matrix I I k)
    (Q : Matrix W I k)
    (R : Matrix I W k) :
    ActBivector (Matrix.fromBlocks ΩW 0 0 ΩI) (Matrix.fromBlocks 1 Q R 1) =
      Matrix.fromBlocks
        (ΩW + Q * ΩI * Qᵀ)
        (ΩW * Rᵀ + Q * ΩI)
        (R * ΩW + ΩI * Qᵀ)
        (R * ΩW * Rᵀ + ΩI) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
      Matrix.mul_add, Matrix.add_mul, Matrix.mul_assoc, add_assoc, add_left_comm, add_comm]

/-- If the third row and column blocks vanish off the diagonal, the `3 × 3` block matrix
is the same as the nested `2 × 2` block diagonal matrix. -/
theorem Block3_diag_eq
    (A11 A12 A21 A22 A33 : Matrix I I k) :
    Block3 A11 A12 0 A21 A22 0 0 0 A33 =
      Matrix.fromBlocks (Matrix.fromBlocks A11 A12 A21 A22) 0 0 A33 := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;> rfl
  · rcases i with i | i <;> fin_cases i <;> fin_cases j <;> rfl
  · fin_cases i <;> rcases j with j | j <;> fin_cases j <;> rfl
  · fin_cases i <;> fin_cases j <;> rfl

/-- The explicit top change of basis carries the Lean repeated-support `n = 4`
representative to the native Magma local model. -/
theorem magmaTopChange_act_onePointRep₁ :
    N4.ActBivector (N4.onePointRep₁ (k := k)) (magmaTopChange (k := k)) =
      magmaTopRep₁ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [N4.ActBivector, magmaTopChange, magmaTopRep₁, N4.onePointRep₁, N4.J,
      E11, E12, E21, E22, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
      Matrix.mul_apply, Fin.sum_univ_two] <;>
    ring_nf

/-- The same top change of basis also carries the second repeated-support basis vector to
the native Magma local model. -/
theorem magmaTopChange_act_onePointRep₂ :
    N4.ActBivector (N4.onePointRep₂ (k := k)) (magmaTopChange (k := k)) =
      magmaTopRep₂ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [N4.ActBivector, magmaTopChange, magmaTopRep₂, N4.onePointRep₂, N4.ω12, N4.J,
      E11, E12, E21, E22, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
      Matrix.mul_apply, Fin.sum_univ_two] <;>
    ring_nf

/-- The ambient mixed representative in Lean coordinates is carried to the native Magma
local model by the explicit block diagonal basis change. -/
theorem magmaChange_act_rep₁ :
    ActBivector (rep₁ (k := k)) (magmaChange (k := k)) = magmaRep₁ (k := k) := by
  calc
    ActBivector (rep₁ (k := k)) (magmaChange (k := k)) =
      Matrix.fromBlocks
        (N4.ActBivector (N4.onePointRep₁ (k := k)) (magmaTopChange (k := k)))
        0
        0
        ((1 : Matrix I I k) * N4.J (k := k) * (1 : Matrix I I k)ᵀ) := by
          simpa [rep₁, magmaChange] using
            act_blockDiagonal
              (k := k)
              (ΩW := N4.onePointRep₁ (k := k))
              (ΩI := N4.J (k := k))
              (H := magmaTopChange (k := k))
              (E := (1 : Matrix I I k))
    _ = magmaRep₁ (k := k) := by
      simp [magmaRep₁, magmaTopChange_act_onePointRep₁, N4.J]

/-- The same ambient basis change also transports the second mixed basis vector to the
native Magma local model. -/
theorem magmaChange_act_rep₂ :
    ActBivector (rep₂ (k := k)) (magmaChange (k := k)) = magmaRep₂ (k := k) := by
  calc
    ActBivector (rep₂ (k := k)) (magmaChange (k := k)) =
      Matrix.fromBlocks
        (N4.ActBivector (N4.onePointRep₂ (k := k)) (magmaTopChange (k := k)))
        0
        0
        ((1 : Matrix I I k) * (0 : Matrix I I k) * (1 : Matrix I I k)ᵀ) := by
          simpa [rep₂, magmaChange] using
            act_blockDiagonal
              (k := k)
              (ΩW := N4.onePointRep₂ (k := k))
              (ΩI := (0 : Matrix I I k))
              (H := magmaTopChange (k := k))
              (E := (1 : Matrix I I k))
    _ = magmaRep₂ (k := k) := by
      simp [magmaRep₂, magmaTopChange_act_onePointRep₂]

/-- The top repeated-support factor and the bottom simple-support factor act pointwise on
the mixed one-point representative. -/
theorem pointwise_levi_family
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    FixesPairBivector
      (Block3 A B 0 0 A 0 0 0 H) := by
  have htop :
      N4.FixesOnePointPairBivector (Matrix.fromBlocks A B 0 A) := by
    exact
      (N4Summary.onePoint_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := B)
        (C := 0)
        (D := A)).2 ⟨hA, rfl, rfl, hB⟩
  have hbot : H * N4.J * Hᵀ = N4.J := by
    rw [N4.mul_J_transpose_mul, hH]
    simp [N4.J]
  have htop₁ :
      N4.ActBivector (N4.onePointRep₁ (k := k)) (Matrix.fromBlocks A B 0 A) =
        N4.onePointRep₁ (k := k) := by
    simpa [N4.FixesBivector, N4.ActBivector] using htop.1
  have htop₂ :
      N4.ActBivector (N4.onePointRep₂ (k := k)) (Matrix.fromBlocks A B 0 A) =
        N4.onePointRep₂ (k := k) := by
    simpa [N4.FixesBivector, N4.ActBivector] using htop.2
  constructor
  · calc
      ActBivector rep₁ (Block3 A B 0 0 A 0 0 0 H) =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₁ (k := k)) (Matrix.fromBlocks A B 0 A))
          0
          0
          (H * N4.J * Hᵀ) := by
            rw [Block3_diag_eq]
            simpa [rep₁] using
              act_blockDiagonal
                (k := k)
                (ΩW := N4.onePointRep₁ (k := k))
                (ΩI := N4.J (k := k))
                (H := Matrix.fromBlocks A B 0 A)
                (E := H)
      _ = Matrix.fromBlocks (N4.onePointRep₁ (k := k)) 0 0 (N4.J (k := k)) := by
            rw [htop₁, hbot]
      _ = rep₁ := by simp [rep₁]
  · calc
      ActBivector rep₂ (Block3 A B 0 0 A 0 0 0 H) =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₂ (k := k)) (Matrix.fromBlocks A B 0 A))
          0
          0
          (H * (0 : Matrix I I k) * Hᵀ) := by
            rw [Block3_diag_eq]
            simpa [rep₂] using
              act_blockDiagonal
                (k := k)
                (ΩW := N4.onePointRep₂ (k := k))
                (ΩI := (0 : Matrix I I k))
                (H := Matrix.fromBlocks A B 0 A)
                (E := H)
      _ = Matrix.fromBlocks (N4.onePointRep₂ (k := k)) 0 0 (0 : Matrix I I k) := by
            rw [htop₂]
            simp
      _ = rep₂ := by simp [rep₂]

/-- The mixed Levi pointwise family has determinant `1`. -/
theorem pointwise_levi_det
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hH : H.det = 1) :
    Matrix.det (Block3 A B 0 0 A 0 0 0 H) = 1 := by
  rw [Block3_diag_eq, Matrix.det_fromBlocks_zero₂₁, Matrix.det_fromBlocks_zero₂₁, hA, hH]
  ring

/-- A general coupled pointwise family for the mixed one-point orbit. The equations here
are exactly the block conditions forced by the bivector action. -/
theorem coupled_pointwise_family
    (Q : Matrix W I k)
    (R : Matrix I W k)
    (h11 : Q * N4.J * Qᵀ = (0 : Matrix W W k))
    (h12 : N4.onePointRep₁ * Rᵀ + Q * N4.J = (0 : Matrix W I k))
    (h21 : R * N4.onePointRep₁ + N4.J * Qᵀ = (0 : Matrix I W k))
    (h22 : R * N4.onePointRep₁ * Rᵀ = (0 : Matrix I I k))
    (h₂left : R * N4.onePointRep₂ = (0 : Matrix I W k))
    (h₂right : N4.onePointRep₂ * Rᵀ = (0 : Matrix W I k)) :
    FixesPairBivector (Matrix.fromBlocks 1 Q R 1) := by
  constructor
  · calc
      ActBivector rep₁ (Matrix.fromBlocks 1 Q R 1) =
        Matrix.fromBlocks
          (N4.onePointRep₁ + Q * N4.J * Qᵀ)
          (N4.onePointRep₁ * Rᵀ + Q * N4.J)
          (R * N4.onePointRep₁ + N4.J * Qᵀ)
          (R * N4.onePointRep₁ * Rᵀ + N4.J) := by
            simpa [rep₁] using
              act_identityOffDiagonal
                (k := k)
                (ΩW := N4.onePointRep₁ (k := k))
                (ΩI := N4.J (k := k))
                (Q := Q)
                (R := R)
      _ = Matrix.fromBlocks (N4.onePointRep₁ (k := k)) 0 0 (N4.J (k := k)) := by
            simp [h11, h12, h21, h22]
      _ = rep₁ := by simp [rep₁]
  · calc
      ActBivector rep₂ (Matrix.fromBlocks 1 Q R 1) =
        Matrix.fromBlocks
          (N4.onePointRep₂ + Q * (0 : Matrix I I k) * Qᵀ)
          (N4.onePointRep₂ * Rᵀ + Q * (0 : Matrix I I k))
          (R * N4.onePointRep₂ + (0 : Matrix I I k) * Qᵀ)
          (R * N4.onePointRep₂ * Rᵀ + (0 : Matrix I I k)) := by
            simpa [rep₂] using
              act_identityOffDiagonal
                (k := k)
                (ΩW := N4.onePointRep₂ (k := k))
                (ΩI := (0 : Matrix I I k))
                (Q := Q)
                (R := R)
      _ = Matrix.fromBlocks (N4.onePointRep₂ (k := k)) 0 0 (0 : Matrix I I k) := by
            simp [h₂left, h₂right]
      _ = rep₂ := by simp [rep₂]

/-- For the mixed one-point orbit, an identity-plus-off-diagonal block matrix fixes the
pair pointwise if and only if its coupled blocks satisfy the exact equations forced by
the bivector action. -/
theorem fixesPair_identityOffDiagonal_iff
    (Q : Matrix W I k)
    (R : Matrix I W k) :
    FixesPairBivector (Matrix.fromBlocks 1 Q R 1) ↔
      Q * N4.J * Qᵀ = (0 : Matrix W W k) ∧
        N4.onePointRep₁ * Rᵀ + Q * N4.J = (0 : Matrix W I k) ∧
        R * N4.onePointRep₁ + N4.J * Qᵀ = (0 : Matrix I W k) ∧
        R * N4.onePointRep₁ * Rᵀ = (0 : Matrix I I k) ∧
        R * N4.onePointRep₂ = (0 : Matrix I W k) ∧
        N4.onePointRep₂ * Rᵀ = (0 : Matrix W I k) := by
  constructor
  · rintro ⟨h₁, h₂⟩
    have h₁' :
        Matrix.fromBlocks
          (N4.onePointRep₁ + Q * N4.J * Qᵀ)
          (N4.onePointRep₁ * Rᵀ + Q * N4.J)
          (R * N4.onePointRep₁ + N4.J * Qᵀ)
          (R * N4.onePointRep₁ * Rᵀ + N4.J) =
        Matrix.fromBlocks
          (N4.onePointRep₁ (k := k))
          0
          0
          (N4.J (k := k)) := by
      have htmp :
          ActBivector
            (Matrix.fromBlocks (N4.onePointRep₁ (k := k)) 0 0 (N4.J (k := k)))
            (Matrix.fromBlocks 1 Q R 1) =
          Matrix.fromBlocks
            (N4.onePointRep₁ (k := k))
            0
            0
            (N4.J (k := k)) := by
        simpa [rep₁] using h₁
      rw [act_identityOffDiagonal
        (k := k)
        (ΩW := N4.onePointRep₁ (k := k))
        (ΩI := N4.J (k := k))
        (Q := Q)
        (R := R)] at htmp
      simpa [rep₁] using htmp
    have h₂' :
        Matrix.fromBlocks
          (N4.onePointRep₂ + Q * (0 : Matrix I I k) * Qᵀ)
          (N4.onePointRep₂ * Rᵀ + Q * (0 : Matrix I I k))
          (R * N4.onePointRep₂ + (0 : Matrix I I k) * Qᵀ)
          (R * N4.onePointRep₂ * Rᵀ + (0 : Matrix I I k)) =
        Matrix.fromBlocks
          (N4.onePointRep₂ (k := k))
          0
          0
          (0 : Matrix I I k) := by
      have htmp :
          ActBivector
            (Matrix.fromBlocks (N4.onePointRep₂ (k := k)) 0 0 (0 : Matrix I I k))
            (Matrix.fromBlocks 1 Q R 1) =
          Matrix.fromBlocks
            (N4.onePointRep₂ (k := k))
            0
            0
            (0 : Matrix I I k) := by
        simpa [rep₂] using h₂
      rw [act_identityOffDiagonal
        (k := k)
        (ΩW := N4.onePointRep₂ (k := k))
        (ΩI := (0 : Matrix I I k))
        (Q := Q)
        (R := R)] at htmp
      simpa [rep₂] using htmp
    have h11 := congrArg Matrix.toBlocks₁₁ h₁'
    have h12 := congrArg Matrix.toBlocks₁₂ h₁'
    have h21 := congrArg Matrix.toBlocks₂₁ h₁'
    have h22 := congrArg Matrix.toBlocks₂₂ h₁'
    have h₂left := congrArg Matrix.toBlocks₂₁ h₂'
    have h₂right := congrArg Matrix.toBlocks₁₂ h₂'
    refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
    · simpa using h11
    · simpa using h12
    · simpa using h21
    · simpa using h22
    · simpa using h₂left
    · simpa using h₂right
  · rintro ⟨h11, h12, h21, h22, h₂left, h₂right⟩
    exact
      coupled_pointwise_family
        (k := k)
        (Q := Q)
        (R := R)
        h11
        h12
        h21
        h22
        h₂left
        h₂right

/-- The upper-right block for a simple coupled `E_{12}` family. -/
def coupledQ12 (t : k) : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl i => (t • E12 (k := k)) i j
    | Sum.inr _ => 0

/-- The lower-left block for a simple coupled `E_{12}` family. -/
def coupledR12 (t : k) : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => (t • E12 (k := k)) i j

/-- The upper-right block for a simple coupled `E_{21}` family. -/
def coupledQ21 (t : k) : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl i => (t • E21 (k := k)) i j
    | Sum.inr _ => 0

/-- The lower-left block for a simple coupled `E_{21}` family. -/
def coupledR21 (t : k) : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => (t • E21 (k := k)) i j

/-- The upper-right block for a simple coupled `E_{11}` family. -/
def coupledQ11 (t : k) : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl i => (-t • E22 (k := k)) i j
    | Sum.inr _ => 0

/-- The lower-left block for a simple coupled `E_{11}` family. -/
def coupledR11 (t : k) : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => (t • E11 (k := k)) i j

/-- The upper-right block for a simple coupled `E_{22}` family. -/
def coupledQ22 (t : k) : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl i => (-t • E11 (k := k)) i j
    | Sum.inr _ => 0

/-- The lower-left block for a simple coupled `E_{22}` family. -/
def coupledR22 (t : k) : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => (t • E22 (k := k)) i j

/-- In the native Magma local model, one mixed root subgroup becomes the Lean `E22`
coupled family. -/
theorem coupledE22_Block3
    (t : k) :
    Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1 =
      Block3
        (1 : Matrix I I k) 0 (-t • E11 (k := k))
        0 (1 : Matrix I I k) 0
        0 (t • E22 (k := k)) (1 : Matrix I I k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;> fin_cases i <;> fin_cases j <;>
      simp [Block3, coupledQ22, coupledR22, E11, E22]
  · rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [Block3, coupledQ22, coupledR22, E11, E22]
  · fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [Block3, coupledQ22, coupledR22, E11, E22]
  · fin_cases i <;> fin_cases j <;>
      simp [Block3, coupledQ22, coupledR22, E11, E22]

/-- The second transported Magma mixed root subgroup in Lean coordinates. -/
def magmaMixedQ2 (t : k) : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl i => (-t • E12 (k := k)) i j
    | Sum.inr _ => 0

/-- The lower-left block of the second transported Magma mixed root subgroup. -/
def magmaMixedR2 (t : k) : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => (-t • E12 (k := k)) i j

/-- In Lean coordinates, the second transported Magma mixed root subgroup has an explicit
`3 × 3` block form. -/
theorem magmaMixed2_Block3
    (t : k) :
    Matrix.fromBlocks 1 (magmaMixedQ2 (k := k) t) (magmaMixedR2 (k := k) t) 1 =
      Block3
        (1 : Matrix I I k) 0 (-t • E12 (k := k))
        0 (1 : Matrix I I k) 0
        0 (-t • E12 (k := k)) (1 : Matrix I I k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;> fin_cases i <;> fin_cases j <;>
      simp [Block3, magmaMixedQ2, magmaMixedR2, E12]
  · rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [Block3, magmaMixedQ2, magmaMixedR2, E12]
  · fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [Block3, magmaMixedQ2, magmaMixedR2, E12]
  · fin_cases i <;> fin_cases j <;>
      simp [Block3, magmaMixedQ2, magmaMixedR2, E12]

/-- In the current Lean basis, the second transported Magma mixed root subgroup is exactly
the existing `E12` coupled family with parameter `-t`. -/
theorem magmaMixed2_eq_coupledE12
    (t : k) :
    Matrix.fromBlocks 1 (magmaMixedQ2 (k := k) t) (magmaMixedR2 (k := k) t) 1 =
      Matrix.fromBlocks 1 (coupledQ12 (k := k) (-t)) (coupledR12 (k := k) (-t)) 1 := by
  rfl

/-- The second transported Magma mixed root subgroup fixes the mixed one-point pair
pointwise in Lean coordinates. -/
theorem magmaMixed2_pointwise
    (t : k) :
    FixesPairBivector (Matrix.fromBlocks 1 (magmaMixedQ2 (k := k) t) (magmaMixedR2 (k := k) t) 1) := by
  apply coupled_pointwise_family
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [magmaMixedQ2, E12, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [magmaMixedQ2, magmaMixedR2, E12, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [magmaMixedQ2, magmaMixedR2, E12, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [magmaMixedR2, E12, N4.J, N4.onePointRep₁, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [magmaMixedR2, E12, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [magmaMixedR2, E12, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]

/-- The upper-right block attached to a general `2 × 2` matrix in the surviving
mixed coupled cell. -/
def coupledQFrom (X : Matrix I I k) : Matrix W I k :=
  fun i j =>
    match i with
    | Sum.inl i => (N4.J (k := k) * Xᵀ * N4.J (k := k)) i j
    | Sum.inr _ => 0

/-- The lower-left block attached to a general `2 × 2` matrix in the surviving
mixed coupled cell. -/
def coupledRFrom (X : Matrix I I k) : Matrix I W k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr j => X i j

/-- The top `2 × 2` off-diagonal block attached to a general mixed `2 × 2` parameter. -/
def coupledUpperFrom (X : Matrix I I k) : Matrix I I k :=
  !![
    -(X 1 1), X 0 1;
    X 1 0, -(X 0 0)
  ]

/-- The central repeated-support correction attached to a general mixed `2 × 2` parameter. -/
def coupledCenterFrom (X : Matrix I I k) : Matrix I I k :=
  !![
    X 0 1 * X 1 0, X 0 1 * X 1 1;
    X 0 0 * X 1 0, -(X 0 0 * X 1 1)
  ]

/-- The trace-zero central layer on the upper-middle `2 × 2` block in the Magma
unipotent model. -/
def centralBlock (M : Matrix I I k) : Matrix V V k :=
  Block3
    (1 : Matrix I I k) M 0
    0 (1 : Matrix I I k) 0
    0 0 (1 : Matrix I I k)

/-- The explicit `4`-parameter unipotent mixed family matching the quotient
`U_{2,1}/Z_{2,1} ≅ G_a^4`. -/
def coupledKernelFrom (X : Matrix I I k) : Matrix V V k :=
  Block3
    (1 : Matrix I I k) (coupledCenterFrom (k := k) X) (coupledUpperFrom (k := k) X)
    0 (1 : Matrix I I k) 0
    0 X (1 : Matrix I I k)

/-- On the surviving mixed coupled cell, the pointwise stabilizer condition is exactly
the rank-drop equation `det X = 0`. -/
theorem coupledFrom_pointwise_iff_det_zero
    (X : Matrix I I k) :
    FixesPairBivector
      (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) ↔
      X.det = 0 := by
  constructor
  · intro hfix
    rcases
        (fixesPair_identityOffDiagonal_iff
          (k := k)
          (Q := coupledQFrom (k := k) X)
          (R := coupledRFrom (k := k) X)).1 hfix with
      ⟨h11, _, _, _, _, _⟩
    have hentry := congrArg (fun M => M (Sum.inl (0 : N4.I)) (Sum.inl (1 : N4.I))) h11
    have hentry' : -(X 0 1 * X 1 0) + X 1 1 * X 0 0 = 0 := by
      simpa [coupledQFrom, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
        dotProduct, Fin.sum_univ_two] using hentry
    simpa [Matrix.det_fin_two, sub_eq_add_neg, add_comm, add_left_comm, add_assoc, mul_comm,
      mul_left_comm, mul_assoc] using hentry'
  · intro hX
    have hdet' : X 0 0 * X 1 1 - X 0 1 * X 1 0 = 0 := by
      simpa [Matrix.det_fin_two] using hX
    have hdetSwap : X 0 1 * X 1 0 + -(X 0 0 * X 1 1) = 0 := by
      have hdetNeg : -(X 0 0 * X 1 1 - X 0 1 * X 1 0) = 0 := by
        simpa [hdet']
      simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc, mul_comm, mul_left_comm,
        mul_assoc] using hdetNeg
    apply coupled_pointwise_family
    ·
      ext i j
      rcases i with i | i <;> rcases j with j | j
      ·
        fin_cases i <;> fin_cases j
        ·
          simp [coupledQFrom, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
            dotProduct, Fin.sum_univ_two]
          ring
        ·
          simpa [coupledQFrom, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
            dotProduct, Fin.sum_univ_two, sub_eq_add_neg, add_comm, add_left_comm, add_assoc,
            mul_comm, mul_left_comm, mul_assoc] using hdet'
        ·
          simpa [coupledQFrom, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
            dotProduct, Fin.sum_univ_two, sub_eq_add_neg, add_comm, add_left_comm, add_assoc,
            mul_comm, mul_left_comm, mul_assoc] using hdetSwap
        ·
          simp [coupledQFrom, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
            dotProduct, Fin.sum_univ_two]
          ring
      ·
        simp [coupledQFrom, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
          Fin.sum_univ_two]
      ·
        simp [coupledQFrom, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
          Fin.sum_univ_two]
      ·
        simp [coupledQFrom, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
          Fin.sum_univ_two]
    ·
      ext i j
      rcases i with i | i
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledQFrom, coupledRFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two]
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledQFrom, coupledRFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two]
    ·
      ext i j
      rcases j with j | j
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledQFrom, coupledRFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two]
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledQFrom, coupledRFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two]
    ·
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [coupledRFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks, Matrix.transpose_apply,
          Matrix.mul_apply, Fin.sum_univ_two]
    ·
      ext i j
      rcases j with j | j
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledRFrom, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledRFrom, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
    ·
      ext i j
      rcases i with i | i
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledRFrom, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
      ·
        fin_cases i <;> fin_cases j <;>
          simp [coupledRFrom, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks,
            Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]

/-- The intrinsic mixed coupled cell has vanishing Schur correction because its
upper-right and lower-left supports lie in complementary summands of `W`. -/
theorem coupledRFrom_mul_coupledQFrom
    (X : Matrix I I k) :
    coupledRFrom (k := k) X * coupledQFrom (k := k) X = 0 := by
  ext i j
  simp [coupledRFrom, coupledQFrom, Matrix.mul_apply, Fintype.sum_sum_type]

/-- Every element of the intrinsic mixed coupled cell has determinant `1`. -/
theorem coupledFrom_det
    (X : Matrix I I k) :
    Matrix.det (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) = 1 := by
  rw [Matrix.det_fromBlocks_one₁₁, coupledRFrom_mul_coupledQFrom]
  simp

/-- Rank-one `2 × 2` matrices give a clean four-parameter pointwise family on the
intrinsic mixed coupled cell. -/
theorem coupledOuter_pointwise
    (u v : I → k) :
    let X : Matrix I I k := fun i j => u i * v j
    FixesPairBivector
      (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
  intro X
  apply (coupledFrom_pointwise_iff_det_zero (k := k) X).2
  simp [X, Matrix.det_fin_two]
  ring

/-- The rank-one mixed coupled family remains pointwise after multiplication by a
pointwise Levi element. -/
theorem coupledOuterLevi_product_pointwise
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let X : Matrix I I k := fun i j => u i * v j
    FixesPairBivector
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H)) := by
  intro X
  have hU : FixesPairBivector
      (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
    simpa [X] using coupledOuter_pointwise (k := k) u v
  have hL := pointwise_levi_family (k := k) A B H hA hB hH
  constructor
  ·
    calc
      ActBivector rep₁
          ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
            (Block3 A B 0 0 A 0 0 0 H)) =
          ActBivector
            (ActBivector rep₁ (Block3 A B 0 0 A 0 0 0 H))
            (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
              rw [actBivector_mul]
      _ = ActBivector rep₁ (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
            rw [hL.1]
      _ = rep₁ := hU.1
  ·
    calc
      ActBivector rep₂
          ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
            (Block3 A B 0 0 A 0 0 0 H)) =
          ActBivector
            (ActBivector rep₂ (Block3 A B 0 0 A 0 0 0 H))
            (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
              rw [actBivector_mul]
      _ = ActBivector rep₂ (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
            rw [hL.2]
      _ = rep₂ := hU.2

/-- Multiplying a general intrinsic mixed coupled element by a pointwise Levi element
does not change its determinant. -/
theorem coupledFromLevi_product_det
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hH : H.det = 1) :
    Matrix.det
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H)) = 1 := by
  rw [Matrix.det_mul, coupledFrom_det]
  simpa using pointwise_levi_det (k := k) A B H hA hH

/-- The rank-one mixed coupled family also has determinant `1` after multiplication by a
pointwise Levi element. -/
theorem coupledOuterLevi_product_det
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hH : H.det = 1) :
    let X : Matrix I I k := fun i j => u i * v j
    Matrix.det
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H)) = 1 := by
  intro X
  simpa [X] using coupledFromLevi_product_det (k := k) X A B H hA hH

/-- A `2 × 2` matrix over a field has determinant zero exactly when it factors as an
outer product of two vectors. -/
theorem exists_outer_of_det_zero
    (X : Matrix I I k)
    (hX : X.det = 0) :
    ∃ u v : I → k, X = fun i j => u i * v j := by
  have hdetEq : X 0 0 * X 1 1 = X 0 1 * X 1 0 := by
    exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hX)
  by_cases h00 : X 0 0 = 0
  · have hprod : X 0 1 * X 1 0 = 0 := by
      have h0 : (0 : k) = X 0 1 * X 1 0 := by
        simpa [h00] using hdetEq
      simpa using h0.symm
    rcases mul_eq_zero.mp hprod with h01 | h10
    · refine ⟨![0, 1], ![X 1 0, X 1 1], ?_⟩
      ext i j
      fin_cases i <;> fin_cases j <;> simp [h00, h01]
    · refine ⟨![X 0 1, X 1 1], ![0, 1], ?_⟩
      ext i j
      fin_cases i <;> fin_cases j <;> simp [h00, h10]
  · let u : I → k := ![X 0 0, X 1 0]
    let v : I → k := ![1, X 0 1 * (X 0 0)⁻¹]
    have h11 : X 1 1 = X 1 0 * (X 0 1 * (X 0 0)⁻¹) := by
      field_simp [h00]
      simpa [mul_comm, mul_left_comm, mul_assoc] using hdetEq
    refine ⟨u, v, ?_⟩
    ext i j
    fin_cases i <;> fin_cases j
    · simp [u, v]
    · calc
        X 0 1 = X 0 1 * (X 0 0 * (X 0 0)⁻¹) := by simp [h00]
        _ = u 0 * v 1 := by
          simp [u, v]
          ring
    · simp [u, v]
    · simp [u, v, h11]

/-- On the intrinsic mixed coupled cell, pointwise fixing is equivalent to belonging to
the explicit rank-one family `X = u vᵀ`. -/
theorem coupledFrom_pointwise_iff_exists_outer
    (X : Matrix I I k) :
    FixesPairBivector
      (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  constructor
  · intro hfix
    exact exists_outer_of_det_zero
      (k := k)
      X
      ((coupledFrom_pointwise_iff_det_zero (k := k) X).1 hfix)
  · rintro ⟨u, v, rfl⟩
    simpa using coupledOuter_pointwise (k := k) u v

/-- An arbitrary identity-plus-off-diagonal pointwise fixer on the mixed row is
already intrinsic: it comes from a surviving coupled block `X` with `det X = 0`. -/
theorem fixesPair_identityOffDiagonal_iff_exists_coupledFrom_det_zero
    (Q : Matrix W I k)
    (R : Matrix I W k) :
    FixesPairBivector (Matrix.fromBlocks 1 Q R 1) ↔
      ∃ X : Matrix I I k,
        Q = coupledQFrom (k := k) X ∧
        R = coupledRFrom (k := k) X ∧
        X.det = 0 := by
  constructor
  · intro hfix
    rcases
      (fixesPair_identityOffDiagonal_iff
        (k := k)
        (Q := Q)
        (R := R)).1 hfix with
      ⟨h11, h12, h21, h22, h₂left, h₂right⟩
    let X : Matrix I I k := fun i j => R i (Sum.inr j)
    have hR_left0 (i : I) : R i (Sum.inl (0 : I)) = 0 := by
      have hEq := congrArg (fun M => M i (Sum.inl (1 : I))) h₂left
      simpa [N4.onePointRep₂, N4.ω12, N4.J, Matrix.fromBlocks, Matrix.mul_apply,
        Fin.sum_univ_two] using hEq
    have hR_left1 (i : I) : R i (Sum.inl (1 : I)) = 0 := by
      have hEq := congrArg (fun M => M i (Sum.inl (0 : I))) h₂left
      simpa [N4.onePointRep₂, N4.ω12, N4.J, Matrix.fromBlocks, Matrix.mul_apply,
        Fin.sum_univ_two] using hEq
    have hR : R = coupledRFrom (k := k) X := by
      ext i j
      rcases j with j | j
      · fin_cases j
        · simpa [coupledRFrom] using hR_left0 i
        · simpa [coupledRFrom] using hR_left1 i
      · rfl
    have hQ : Q = coupledQFrom (k := k) X := by
      ext i j
      rcases i with i | i
      · fin_cases i <;> fin_cases j
        ·
          have hEq := congrArg (fun M => M (Sum.inl (0 : I)) (1 : I)) h12
          have hcq :
              coupledQFrom (k := k) X (Sum.inl (0 : I)) (0 : I) =
                -R (1 : I) (Sum.inr (1 : I)) := by
            simp [coupledQFrom, X, N4.J, Matrix.transpose_apply, Matrix.mul_apply,
              Matrix.vecMul, dotProduct, Fin.sum_univ_two]
          have hgoal :
              Q (Sum.inl (0 : I)) (0 : I) =
                coupledQFrom (k := k) X (Sum.inl (0 : I)) (0 : I) := by
            rw [hcq, eq_neg_iff_add_eq_zero]
            simpa [N4.onePointRep₁, N4.J, Matrix.fromBlocks, hR_left0, hR_left1,
              Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
              Fin.sum_univ_two, add_comm, add_left_comm, add_assoc] using hEq
          simpa using hgoal
        ·
          have hEq := congrArg (fun M => M (Sum.inl (0 : I)) (0 : I)) h12
          have hcq :
              coupledQFrom (k := k) X (Sum.inl (0 : I)) (1 : I) =
                R (0 : I) (Sum.inr (1 : I)) := by
            simp [coupledQFrom, X, N4.J, Matrix.transpose_apply, Matrix.mul_apply,
              Matrix.vecMul, dotProduct, Fin.sum_univ_two]
          have hgoal :
              Q (Sum.inl (0 : I)) (1 : I) =
                coupledQFrom (k := k) X (Sum.inl (0 : I)) (1 : I) := by
            have hEq' := congrArg Neg.neg hEq
            rw [hcq]
            apply sub_eq_zero.mp
            simpa [sub_eq_add_neg, N4.onePointRep₁, N4.J, Matrix.fromBlocks, hR_left0, hR_left1,
              Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
              Fin.sum_univ_two, add_comm, add_left_comm, add_assoc] using hEq'
          simpa using hgoal
        ·
          have hEq := congrArg (fun M => M (Sum.inl (1 : I)) (1 : I)) h12
          have hcq :
              coupledQFrom (k := k) X (Sum.inl (1 : I)) (0 : I) =
                R (1 : I) (Sum.inr (0 : I)) := by
            simp [coupledQFrom, X, N4.J, Matrix.transpose_apply, Matrix.mul_apply,
              Matrix.vecMul, dotProduct, Fin.sum_univ_two]
          have hgoal :
              Q (Sum.inl (1 : I)) (0 : I) =
                coupledQFrom (k := k) X (Sum.inl (1 : I)) (0 : I) := by
            rw [hcq]
            apply sub_eq_zero.mp
            simpa [sub_eq_add_neg, N4.onePointRep₁, N4.J, Matrix.fromBlocks, hR_left0, hR_left1,
              Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
              Fin.sum_univ_two, add_comm, add_left_comm, add_assoc] using hEq
          simpa using hgoal
        ·
          have hEq := congrArg (fun M => M (Sum.inl (1 : I)) (0 : I)) h12
          have hcq :
              coupledQFrom (k := k) X (Sum.inl (1 : I)) (1 : I) =
                -R (0 : I) (Sum.inr (0 : I)) := by
            simp [coupledQFrom, X, N4.J, Matrix.transpose_apply, Matrix.mul_apply,
              Matrix.vecMul, dotProduct, Fin.sum_univ_two]
          have hgoal :
              Q (Sum.inl (1 : I)) (1 : I) =
                coupledQFrom (k := k) X (Sum.inl (1 : I)) (1 : I) := by
            rw [hcq, eq_neg_iff_add_eq_zero]
            have hEq' := congrArg Neg.neg hEq
            simpa [N4.onePointRep₁, N4.J, Matrix.fromBlocks, hR_left0, hR_left1,
              Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
              Fin.sum_univ_two, add_comm, add_left_comm, add_assoc] using hEq'
          simpa using hgoal
      · fin_cases i <;> fin_cases j
        ·
          have hEq := congrArg (fun M => M (Sum.inr (0 : I)) (1 : I)) h12
          simpa [X, coupledQFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            hR_left0, hR_left1,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two] using hEq
        ·
          have hEq := congrArg (fun M => M (Sum.inr (0 : I)) (0 : I)) h12
          simpa [X, coupledQFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            hR_left0, hR_left1,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two] using hEq
        ·
          have hEq := congrArg (fun M => M (Sum.inr (1 : I)) (1 : I)) h12
          simpa [X, coupledQFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            hR_left0, hR_left1,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two] using hEq
        ·
          have hEq := congrArg (fun M => M (Sum.inr (1 : I)) (0 : I)) h12
          simpa [X, coupledQFrom, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
            hR_left0, hR_left1,
            Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
            Fin.sum_univ_two] using hEq
    have hfixX :
        FixesPairBivector
          (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
      simpa [hQ, hR] using hfix
    refine ⟨X, hQ, hR, ?_⟩
    exact (coupledFrom_pointwise_iff_det_zero (k := k) X).1 hfixX
  · rintro ⟨X, rfl, rfl, hX⟩
    exact (coupledFrom_pointwise_iff_det_zero (k := k) X).2 hX

/-- A simple one-parameter coupled family on the mixed one-point orbit. -/
theorem coupledE12_pointwise
    (t : k) :
    FixesPairBivector (Matrix.fromBlocks 1 (coupledQ12 (k := k) t) (coupledR12 (k := k) t) 1) := by
  apply coupled_pointwise_family
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ12, E12, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ12, coupledR12, E12, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ12, coupledR12, E12, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR12, E12, N4.J, N4.onePointRep₁, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR12, E12, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR12, E12, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]

/-- A second simple one-parameter coupled family on the mixed one-point orbit. -/
theorem coupledE21_pointwise
    (t : k) :
    FixesPairBivector (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1) := by
  apply coupled_pointwise_family
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ21, E21, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ21, coupledR21, E21, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ21, coupledR21, E21, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR21, E21, N4.J, N4.onePointRep₁, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR21, E21, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    rcases i with i | i
    · fin_cases i <;> fin_cases j <;>
        simp [coupledR21, E21, N4.onePointRep₂, N4.ω12, N4.J, Matrix.fromBlocks,
          Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
    · fin_cases i <;> fin_cases j <;>
        simp [coupledR21, E21, N4.onePointRep₂, N4.ω12, N4.J, Matrix.fromBlocks,
          Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]

/-- A third simple one-parameter coupled family on the mixed one-point orbit. -/
theorem coupledE11_pointwise
    (t : k) :
    FixesPairBivector (Matrix.fromBlocks 1 (coupledQ11 (k := k) t) (coupledR11 (k := k) t) 1) := by
  apply coupled_pointwise_family
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ11, E11, E22, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ11, coupledR11, E11, E22, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ11, coupledR11, E11, E22, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR11, E11, N4.J, N4.onePointRep₁, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR11, E11, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR11, E11, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]

/-- A fourth simple one-parameter coupled family on the mixed one-point orbit. -/
theorem coupledE22_pointwise
    (t : k) :
    FixesPairBivector (Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1) := by
  apply coupled_pointwise_family
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ22, E11, E22, N4.J, Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ22, coupledR22, E11, E22, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledQ22, coupledR22, E11, E22, N4.J, N4.onePointRep₁, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR22, E22, N4.J, N4.onePointRep₁, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR22, E22, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  · ext i j
    fin_cases i <;> fin_cases j <;>
      simp [coupledR22, E22, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks, Matrix.transpose_apply,
        Matrix.mul_apply, Fin.sum_univ_two]

/-- The coupled `R12` and `Q12` blocks have vanishing Schur correction. -/
theorem coupledR12_mul_coupledQ12
    (t : k) :
    coupledR12 (k := k) t * coupledQ12 (k := k) t = 0 := by
  ext i j
  simp [Matrix.mul_apply, coupledR12, coupledQ12, Fintype.sum_sum_type]

/-- The coupled `R21` and `Q21` blocks have vanishing Schur correction. -/
theorem coupledR21_mul_coupledQ21
    (t : k) :
    coupledR21 (k := k) t * coupledQ21 (k := k) t = 0 := by
  ext i j
  simp [Matrix.mul_apply, coupledR21, coupledQ21, Fintype.sum_sum_type]

/-- The coupled `R11` and `Q11` blocks have vanishing Schur correction. -/
theorem coupledR11_mul_coupledQ11
    (t : k) :
    coupledR11 (k := k) t * coupledQ11 (k := k) t = 0 := by
  ext i j
  simp [Matrix.mul_apply, coupledR11, coupledQ11, Fintype.sum_sum_type]

/-- The coupled `R22` and `Q22` blocks have vanishing Schur correction. -/
theorem coupledR22_mul_coupledQ22
    (t : k) :
    coupledR22 (k := k) t * coupledQ22 (k := k) t = 0 := by
  ext i j
  simp [Matrix.mul_apply, coupledR22, coupledQ22, Fintype.sum_sum_type]

/-- Each basic explicit mixed coupled generator has determinant `1`. -/
theorem coupledE12_det
    (t : k) :
    Matrix.det (Matrix.fromBlocks 1 (coupledQ12 (k := k) t) (coupledR12 (k := k) t) 1) = 1 := by
  rw [Matrix.det_fromBlocks_one₁₁, coupledR12_mul_coupledQ12]
  simp

/-- Each basic explicit mixed coupled generator has determinant `1`. -/
theorem coupledE21_det
    (t : k) :
    Matrix.det (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1) = 1 := by
  rw [Matrix.det_fromBlocks_one₁₁, coupledR21_mul_coupledQ21]
  simp

/-- Each basic explicit mixed coupled generator has determinant `1`. -/
theorem coupledE11_det
    (t : k) :
    Matrix.det (Matrix.fromBlocks 1 (coupledQ11 (k := k) t) (coupledR11 (k := k) t) 1) = 1 := by
  rw [Matrix.det_fromBlocks_one₁₁, coupledR11_mul_coupledQ11]
  simp

/-- Each basic explicit mixed coupled generator has determinant `1`. -/
theorem coupledE22_det
    (t : k) :
    Matrix.det (Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1) = 1 := by
  rw [Matrix.det_fromBlocks_one₁₁, coupledR22_mul_coupledQ22]
  simp

/-- Setting the `E12` parameter to zero gives the identity matrix. -/
theorem coupledE12_zero_eq_one :
    Matrix.fromBlocks 1 (coupledQ12 (k := k) (0 : k)) (coupledR12 (k := k) (0 : k)) 1 =
      (1 : Matrix V V k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;> simp [Matrix.one_apply]
  · rcases i with i | i <;> simp [coupledQ12]
  · rcases j with j | j <;> simp [coupledR12]
  · simp [Matrix.one_apply]

/-- Setting the `E21` parameter to zero gives the identity matrix. -/
theorem coupledE21_zero_eq_one :
    Matrix.fromBlocks 1 (coupledQ21 (k := k) (0 : k)) (coupledR21 (k := k) (0 : k)) 1 =
      (1 : Matrix V V k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;> simp [Matrix.one_apply]
  · rcases i with i | i <;> simp [coupledQ21]
  · rcases j with j | j <;> simp [coupledR21]
  · simp [Matrix.one_apply]

/-- Setting the `E11` parameter to zero gives the identity matrix. -/
theorem coupledE11_zero_eq_one :
    Matrix.fromBlocks 1 (coupledQ11 (k := k) (0 : k)) (coupledR11 (k := k) (0 : k)) 1 =
      (1 : Matrix V V k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;> simp [Matrix.one_apply]
  · rcases i with i | i <;> simp [coupledQ11]
  · rcases j with j | j <;> simp [coupledR11]
  · simp [Matrix.one_apply]

/-- Setting the `E22` parameter to zero gives the identity matrix. -/
theorem coupledE22_zero_eq_one :
    Matrix.fromBlocks 1 (coupledQ22 (k := k) (0 : k)) (coupledR22 (k := k) (0 : k)) 1 =
      (1 : Matrix V V k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;> simp [Matrix.one_apply]
  · rcases i with i | i <;> simp [coupledQ22]
  · rcases j with j | j <;> simp [coupledR22]
  · simp [Matrix.one_apply]

/-- The older hard-coded `crossUnipotent₂` direction is a congruence-side generator,
but in the bivector model it does not fix the mixed one-point pair unless the parameter
vanishes. -/
theorem cross₂_not_pointwise
    {t : k}
    (ht : t ≠ 0) :
    ¬ FixesPairBivector (Matrix.fromBlocks 1 (t • crossQ₂ (k := k)) (t • crossR₂ (k := k)) 1) := by
  intro hfix
  rcases
      (fixesPair_identityOffDiagonal_iff
        (k := k)
        (Q := t • crossQ₂ (k := k))
        (R := t • crossR₂ (k := k))).1 hfix with
    ⟨_, _, _, _, h₂left, _⟩
  have hentry := congrArg (fun M => M (0 : I) (Sum.inl (0 : N4.I))) h₂left
  simp [crossR₂, E12, N4.onePointRep₂, N4.ω12, Matrix.fromBlocks,
    Matrix.mul_apply, Fin.sum_univ_two, N4.J] at hentry
  exact ht hentry

/-- The older hard-coded `crossUnipotent₁` direction is also a congruence-side
generator, but in the bivector model it does not fix the mixed one-point pair unless
its parameter vanishes. -/
theorem cross₁_not_pointwise
    {t : k}
    (ht : t ≠ 0) :
    ¬ FixesPairBivector (Matrix.fromBlocks 1 (t • crossQ₁ (k := k)) (t • crossR₁ (k := k)) 1) := by
  intro hfix
  rcases
      (fixesPair_identityOffDiagonal_iff
        (k := k)
        (Q := t • crossQ₁ (k := k))
        (R := t • crossR₁ (k := k))).1 hfix with
    ⟨_, h12, _, _, _, _⟩
  have hentry := congrArg (fun M => M (Sum.inl (0 : N4.I)) (1 : I)) h12
  simp [crossQ₁, crossR₁, E21, E22, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
    Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two] at hentry
  exact ht hentry

/-- The subgroup generated by the two explicit coupled one-parameter families still fixes
the mixed one-point pair pointwise. -/
theorem coupledE12E21_product_pointwise
    (s t : k) :
    FixesPairBivector
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) := by
  have h12 := coupledE12_pointwise (k := k) s
  have h21 := coupledE21_pointwise (k := k) t
  constructor
  · rw [actBivector_mul]
    rw [h21.1]
    exact h12.1
  · rw [actBivector_mul]
    rw [h21.2]
    exact h12.2

/-- The two diagonal explicit coupled families also generate a concrete pointwise
subgroup on the mixed one-point orbit. -/
theorem coupledE11E22_product_pointwise
    (s t : k) :
    FixesPairBivector
      ((Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1)) := by
  have h11 := coupledE11_pointwise (k := k) s
  have h22 := coupledE22_pointwise (k := k) t
  constructor
  · rw [actBivector_mul]
    rw [h22.1]
    exact h11.1
  · rw [actBivector_mul]
    rw [h22.2]
    exact h11.2

/-- The four explicit coupled one-parameter families combine to give a concrete
four-parameter pointwise subgroup on the mixed one-point orbit. -/
theorem coupledExplicit_product_pointwise
    (s₁ t₁ s₂ t₂ : k) :
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    FixesPairBivector ((g12 * g21) * (g11 * g22)) := by
  intro g12 g21 g11 g22
  have hA := coupledE12E21_product_pointwise (k := k) s₁ t₁
  have hB := coupledE11E22_product_pointwise (k := k) s₂ t₂
  constructor
  · rw [actBivector_mul]
    rw [hB.1]
    exact hA.1
  · rw [actBivector_mul]
    rw [hB.2]
    exact hA.2

/-- The two off-diagonal explicit mixed coupled families have determinant `1`. -/
theorem coupledE12E21_product_det
    (s t : k) :
    Matrix.det
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) = 1 := by
  rw [Matrix.det_mul, coupledE12_det, coupledE21_det]
  simp

/-- The two diagonal explicit mixed coupled families have determinant `1`. -/
theorem coupledE11E22_product_det
    (s t : k) :
    Matrix.det
      ((Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1)) = 1 := by
  rw [Matrix.det_mul, coupledE11_det, coupledE22_det]
  simp

/-- The full explicit four-parameter mixed coupled family has determinant `1`. -/
theorem coupledExplicit_product_det
    (s₁ t₁ s₂ t₂ : k) :
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    Matrix.det ((g12 * g21) * (g11 * g22)) = 1 := by
  intro g12 g21 g11 g22
  rw [Matrix.det_mul]
  rw [show Matrix.det (g12 * g21) = 1 by simpa [g12, g21] using coupledE12E21_product_det (k := k) s₁ t₁]
  rw [show Matrix.det (g11 * g22) = 1 by simpa [g11, g22] using coupledE11E22_product_det (k := k) s₂ t₂]
  simp

/-- Combining the explicit coupled unipotent families with the Levi family gives a
concrete pointwise subgroup on the mixed one-point orbit. -/
theorem coupledLevi_product_pointwise
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let gU : Matrix V V k :=
      (Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
      (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)
    let gL : Matrix V V k := Block3 A B 0 0 A 0 0 0 H
    FixesPairBivector (gU * gL) := by
  intro gU gL
  have hU := coupledE12E21_product_pointwise (k := k) s t
  have hL := pointwise_levi_family (k := k) A B H hA hB hH
  constructor
  ·
    calc
      ActBivector rep₁ (gU * gL) =
          ActBivector (ActBivector rep₁ gL) gU := by
            rw [actBivector_mul]
      _ = ActBivector rep₁ gU := by rw [hL.1]
      _ = rep₁ := hU.1
  ·
    calc
      ActBivector rep₂ (gU * gL) =
          ActBivector (ActBivector rep₂ gL) gU := by
            rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by rw [hL.2]
      _ = rep₂ := hU.2

/-- Combining the diagonal explicit coupled subgroup with the Levi family gives a
concrete pointwise subgroup on the mixed one-point orbit. -/
theorem coupledDiagonalLevi_product_pointwise
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    FixesPairBivector gU := by
  intro g11 g22 gU
  have hU := coupledE11E22_product_pointwise (k := k) s t
  have hL := pointwise_levi_family (k := k) A B H hA hB hH
  constructor
  ·
    calc
      ActBivector rep₁ gU =
          ActBivector (ActBivector rep₁ (Block3 A B 0 0 A 0 0 0 H)) (g11 * g22) := by
            rw [actBivector_mul]
      _ = ActBivector rep₁ (g11 * g22) := by rw [hL.1]
      _ = rep₁ := hU.1
  ·
    calc
      ActBivector rep₂ gU =
          ActBivector (ActBivector rep₂ (Block3 A B 0 0 A 0 0 0 H)) (g11 * g22) := by
            rw [actBivector_mul]
      _ = ActBivector rep₂ (g11 * g22) := by rw [hL.2]
      _ = rep₂ := hU.2

/-- Combining the full explicit four-parameter coupled subgroup with the Levi family
gives a concrete pointwise subgroup on the mixed one-point orbit. -/
theorem coupledExplicitLevi_product_pointwise
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    FixesPairBivector gU := by
  intro g12 g21 g11 g22 gU
  have hU := coupledExplicit_product_pointwise (k := k) s₁ t₁ s₂ t₂
  have hL := pointwise_levi_family (k := k) A B H hA hB hH
  constructor
  ·
    calc
      ActBivector rep₁ gU =
          ActBivector (ActBivector rep₁ (Block3 A B 0 0 A 0 0 0 H)) ((g12 * g21) * (g11 * g22)) := by
            rw [actBivector_mul]
      _ = ActBivector rep₁ ((g12 * g21) * (g11 * g22)) := by rw [hL.1]
      _ = rep₁ := hU.1
  ·
    calc
      ActBivector rep₂ gU =
          ActBivector (ActBivector rep₂ (Block3 A B 0 0 A 0 0 0 H)) ((g12 * g21) * (g11 * g22)) := by
            rw [actBivector_mul]
      _ = ActBivector rep₂ ((g12 * g21) * (g11 * g22)) := by rw [hL.2]
      _ = rep₂ := hU.2

/-- The explicit product of the four basic mixed root families has a closed `3 × 3`
block form depending on an arbitrary `2 × 2` parameter matrix. This is the concrete
`G_a^4` quotient family suggested by the Magma local model. -/
theorem coupledKernelFrom_eq_product
    (X : Matrix I I k) :
    coupledKernelFrom (k := k) X =
      ((Matrix.fromBlocks
            1
            (coupledQ12 (k := k) (X 0 1))
            (coupledR12 (k := k) (X 0 1))
            1) *
        (Matrix.fromBlocks
            1
            (coupledQ21 (k := k) (X 1 0))
            (coupledR21 (k := k) (X 1 0))
            1)) *
        ((Matrix.fromBlocks
              1
              (coupledQ11 (k := k) (X 0 0))
              (coupledR11 (k := k) (X 0 0))
              1) *
          (Matrix.fromBlocks
              1
              (coupledQ22 (k := k) (X 1 1))
              (coupledR22 (k := k) (X 1 1))
              1)) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    rcases i with i | i <;> rcases j with j | j <;> fin_cases i <;> fin_cases j <;>
      simp [coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, coupledQ12, coupledR12,
        coupledQ21, coupledR21, coupledQ11, coupledR11, coupledQ22, coupledR22, Block3,
        E11, E12, E21, E22, N4.J, Matrix.fromBlocks_multiply, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_two] <;> ring
  ·
    rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, coupledQ12, coupledR12,
        coupledQ21, coupledR21, coupledQ11, coupledR11, coupledQ22, coupledR22, Block3,
        E11, E12, E21, E22, N4.J, Matrix.fromBlocks_multiply, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_two] <;> ring
  ·
    fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, coupledQ12, coupledR12,
        coupledQ21, coupledR21, coupledQ11, coupledR11, coupledQ22, coupledR22, Block3,
        E11, E12, E21, E22, N4.J, Matrix.fromBlocks_multiply, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_two] <;> ring
  ·
    fin_cases i <;> fin_cases j <;>
      simp [coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, coupledQ12, coupledR12,
        coupledQ21, coupledR21, coupledQ11, coupledR11, coupledQ22, coupledR22, Block3,
        E11, E12, E21, E22, N4.J, Matrix.fromBlocks_multiply, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_two] <;> ring

/-- The closed `4`-parameter mixed quotient family fixes the mixed one-point pair
pointwise. -/
theorem coupledKernelFrom_pointwise
    (X : Matrix I I k) :
    FixesPairBivector (coupledKernelFrom (k := k) X) := by
  have hprod :=
    coupledExplicit_product_pointwise
      (k := k)
      (s₁ := X 0 1)
      (t₁ := X 1 0)
      (s₂ := X 0 0)
      (t₂ := X 1 1)
  simpa [coupledKernelFrom_eq_product (k := k) X] using hprod

/-- The closed `4`-parameter mixed quotient family has determinant `1`. -/
theorem coupledKernelFrom_det
    (X : Matrix I I k) :
    Matrix.det (coupledKernelFrom (k := k) X) = 1 := by
  have hprod :=
    coupledExplicit_product_det
      (k := k)
      (s₁ := X 0 1)
      (t₁ := X 1 0)
      (s₂ := X 0 0)
      (t₂ := X 1 1)
  simpa [coupledKernelFrom_eq_product (k := k) X] using hprod

/-- The central upper-middle layer is just the embedded `n = 4` repeated-support upper
shear. -/
theorem centralBlock_eq
    (M : Matrix I I k) :
    centralBlock (k := k) M =
      Matrix.fromBlocks (N4.onePointUpperShear (k := k) M) 0 0 (1 : Matrix I I k) := by
  rw [centralBlock, Block3_diag_eq]
  rfl

/-- The central upper-middle layer fixes the mixed pair exactly when its trace
vanishes. -/
theorem centralBlock_pointwise_iff
    (M : Matrix I I k) :
    FixesPairBivector (centralBlock (k := k) M) ↔
      M 0 0 + M 1 1 = 0 := by
  have hrep₁ :
      ActBivector (rep₁ (k := k)) (centralBlock (k := k) M) =
        rep₁ (k := k) + (M 0 0 + M 1 1) • rep₂ (k := k) := by
    calc
      ActBivector (rep₁ (k := k)) (centralBlock (k := k) M) =
          Matrix.fromBlocks
            (N4.ActBivector (N4.onePointRep₁ (k := k)) (N4.onePointUpperShear (k := k) M))
            0
            0
            ((1 : Matrix I I k) * N4.J * (1 : Matrix I I k)ᵀ) := by
              simpa [rep₁, centralBlock_eq (k := k) M] using
                act_blockDiagonal
                  (k := k)
                  (ΩW := N4.onePointRep₁ (k := k))
                  (ΩI := N4.J (k := k))
                  (H := N4.onePointUpperShear (k := k) M)
                  (E := (1 : Matrix I I k))
      _ =
          Matrix.fromBlocks
            (N4.onePointRep₁ (k := k) + (M 0 0 + M 1 1) • N4.onePointRep₂ (k := k))
            0
            0
            (N4.J (k := k)) := by
              rcases N4Summary.onePoint_upperShear_action (k := k) M with ⟨h₁, _⟩
              simp [h₁]
      _ = rep₁ (k := k) + (M 0 0 + M 1 1) • rep₂ (k := k) := by
              ext i j
              cases i <;> cases j <;>
                simp [rep₁, rep₂, N4.onePointRep₁, N4.onePointRep₂, N4.ω12,
                  Matrix.fromBlocks, Matrix.add_apply]
  have hrep₂ :
      ActBivector (rep₂ (k := k)) (centralBlock (k := k) M) = rep₂ (k := k) := by
    calc
      ActBivector (rep₂ (k := k)) (centralBlock (k := k) M) =
          Matrix.fromBlocks
            (N4.ActBivector (N4.onePointRep₂ (k := k)) (N4.onePointUpperShear (k := k) M))
            0
            0
            ((0 : Matrix I I k)) := by
              simpa [rep₂, centralBlock_eq (k := k) M] using
                act_blockDiagonal
                  (k := k)
                  (ΩW := N4.onePointRep₂ (k := k))
                  (ΩI := (0 : Matrix I I k))
                  (H := N4.onePointUpperShear (k := k) M)
                  (E := (1 : Matrix I I k))
      _ = rep₂ (k := k) := by
            rcases N4Summary.onePoint_upperShear_action (k := k) M with ⟨_, h₂⟩
            simpa [rep₂, N4.onePointRep₂, N4.ω12] using
              congrArg (fun X => Matrix.fromBlocks X 0 0 (0 : Matrix I I k)) h₂
  constructor
  · intro hfix
    have hEq : rep₁ (k := k) + (M 0 0 + M 1 1) • rep₂ (k := k) = rep₁ (k := k) := by
      simpa [hrep₁] using hfix.1
    have hEntry := congrArg
      (fun Ω : Matrix V V k =>
        Ω (Sum.inl (Sum.inl (0 : I))) (Sum.inl (Sum.inl (1 : I)))) hEq
    simpa [rep₁, rep₂, N4.onePointRep₂, N4.ω12, N4.J, Matrix.fromBlocks, Matrix.add_apply] using hEntry
  · intro htrace
    exact ⟨by simpa [htrace] using hrep₁, hrep₂⟩

/-- The trace-zero central layer always has determinant `1`. -/
theorem centralBlock_det
    (M : Matrix I I k) :
    Matrix.det (centralBlock (k := k) M) = 1 := by
  rw [centralBlock_eq, Matrix.det_fromBlocks_zero₂₁]
  simp [N4Summary.onePoint_upperShear_det]

/-- The full Magma unipotent cell factors as its trace-zero central layer times the
concrete `G_a^4` quotient family. -/
theorem Block3_unipotent_eq_centralBlock_mul_coupledKernelFrom
    (M X : Matrix I I k) :
    Block3
      (1 : Matrix I I k) M (coupledUpperFrom (k := k) X)
      0 (1 : Matrix I I k) 0
      0 X (1 : Matrix I I k) =
      centralBlock (k := k) (M - coupledCenterFrom (k := k) X) *
        coupledKernelFrom (k := k) X := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    rcases i with i | i <;> rcases j with j | j <;> fin_cases i <;> fin_cases j <;>
      simp [centralBlock, coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring
  ·
    rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [centralBlock, coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring
  ·
    fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [centralBlock, coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring
  ·
    fin_cases i <;> fin_cases j <;>
      simp [centralBlock, coupledKernelFrom, coupledCenterFrom, coupledUpperFrom, Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring

set_option maxHeartbeats 2000000 in
/-- On the full Magma `3 × 3` unipotent cell, pointwise fixing is equivalent to the
transported mixed-root upper-right block together with the single trace condition on the
central upper-middle block. -/
theorem fixesPair_Block3_unipotent_iff
    (M U X : Matrix I I k) :
    FixesPairBivector
      (Block3
        (1 : Matrix I I k) M U
        0 (1 : Matrix I I k) 0
        0 X (1 : Matrix I I k)) ↔
      U = coupledUpperFrom (k := k) X ∧
        M 0 0 + M 1 1 = -X.det := by
  constructor
  · intro hfix
    have hU01Eq := congrArg
      (fun Ω : Matrix V V k =>
        Ω (Sum.inl (Sum.inl (0 : I))) (Sum.inr (0 : I))) hfix.1
    have hU01raw : X 0 1 - U 0 1 = 0 := by
      simpa [ActBivector, rep₁, Block3, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two, sub_eq_add_neg,
        add_comm, add_left_comm, add_assoc] using hU01Eq
    have hU01 : U 0 1 = X 0 1 := by
      exact (sub_eq_zero.mp hU01raw).symm
    have hU00Eq := congrArg
      (fun Ω : Matrix V V k =>
        Ω (Sum.inl (Sum.inl (0 : I))) (Sum.inr (1 : I))) hfix.1
    have hU00raw : U 0 0 + X 1 1 = 0 := by
      simpa [ActBivector, rep₁, Block3, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two, sub_eq_add_neg,
        add_comm, add_left_comm, add_assoc] using hU00Eq
    have hU00 : U 0 0 = -(X 1 1) := by
      exact eq_neg_iff_add_eq_zero.mpr hU00raw
    have hU11Eq := congrArg
      (fun Ω : Matrix V V k =>
        Ω (Sum.inr (0 : I)) (Sum.inl (Sum.inl (1 : I)))) hfix.1
    have hU11raw : X 0 0 + U 1 1 = 0 := by
      simpa [ActBivector, rep₁, Block3, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two, sub_eq_add_neg,
        add_comm, add_left_comm, add_assoc] using hU11Eq
    have hU11 : U 1 1 = -(X 0 0) := by
      exact (neg_eq_iff_add_eq_zero.mpr hU11raw).symm
    have hU10Eq := congrArg
      (fun Ω : Matrix V V k =>
        Ω (Sum.inr (1 : I)) (Sum.inl (Sum.inl (1 : I)))) hfix.1
    have hU10raw : X 1 0 - U 1 0 = 0 := by
      simpa [ActBivector, rep₁, Block3, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
        Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two, sub_eq_add_neg,
        add_comm, add_left_comm, add_assoc] using hU10Eq
    have hU10 : U 1 0 = X 1 0 := by
      exact (sub_eq_zero.mp hU10raw).symm
    have hU : U = coupledUpperFrom (k := k) X := by
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [coupledUpperFrom, hU00, hU01, hU10, hU11]
    have hTraceEq := congrArg
      (fun Ω : Matrix V V k =>
        Ω (Sum.inl (Sum.inl (0 : I))) (Sum.inl (Sum.inl (1 : I)))) hfix.1
    have hTrace0 : M 0 0 + M 1 1 + X.det = 0 := by
      have hRaw :
          M 0 0 + M 1 1 + (U 0 0 * U 1 1 - U 0 1 * U 1 0) = 0 := by
        simpa [ActBivector, rep₁, Block3, N4.onePointRep₁, N4.J, Matrix.fromBlocks,
          Matrix.transpose_apply, Matrix.mul_apply, Fin.sum_univ_two, sub_eq_add_neg,
          add_comm, add_left_comm, add_assoc, mul_comm, mul_left_comm, mul_assoc] using hTraceEq
      calc
        M 0 0 + M 1 1 + X.det
            = M 0 0 + M 1 1 + (U 0 0 * U 1 1 - U 0 1 * U 1 0) := by
                rw [Matrix.det_fin_two, hU00, hU01, hU10, hU11]
                ring
        _ = 0 := hRaw
    exact ⟨hU, eq_neg_iff_add_eq_zero.mpr hTrace0⟩
  · rintro ⟨hU, hTrace⟩
    let Z : Matrix I I k := M - coupledCenterFrom (k := k) X
    have hZtrace : Z 0 0 + Z 1 1 = 0 := by
      have hsum : M 0 0 + M 1 1 + X.det = 0 := by
        exact eq_neg_iff_add_eq_zero.mp hTrace
      have hrewrite :
          Z 0 0 + Z 1 1 = M 0 0 + M 1 1 + X.det := by
        simp [Z, coupledCenterFrom, Matrix.det_fin_two, sub_eq_add_neg]
        ring
      rw [hrewrite, hsum]
    have hZfix : FixesPairBivector (centralBlock (k := k) Z) := by
      exact (centralBlock_pointwise_iff (k := k) Z).2 hZtrace
    have hKfix : FixesPairBivector (coupledKernelFrom (k := k) X) :=
      coupledKernelFrom_pointwise (k := k) X
    have hprod : FixesPairBivector (centralBlock (k := k) Z * coupledKernelFrom (k := k) X) := by
      constructor
      · calc
          ActBivector (rep₁ (k := k))
              (centralBlock (k := k) Z * coupledKernelFrom (k := k) X) =
              ActBivector
                (ActBivector (rep₁ (k := k)) (coupledKernelFrom (k := k) X))
                (centralBlock (k := k) Z) := by
                  rw [actBivector_mul]
          _ = ActBivector (rep₁ (k := k)) (centralBlock (k := k) Z) := by rw [hKfix.1]
          _ = rep₁ (k := k) := hZfix.1
      · calc
          ActBivector (rep₂ (k := k))
              (centralBlock (k := k) Z * coupledKernelFrom (k := k) X) =
              ActBivector
                (ActBivector (rep₂ (k := k)) (coupledKernelFrom (k := k) X))
                (centralBlock (k := k) Z) := by
                  rw [actBivector_mul]
          _ = ActBivector (rep₂ (k := k)) (centralBlock (k := k) Z) := by rw [hKfix.2]
          _ = rep₂ (k := k) := hZfix.2
    have hEq :
        Block3
          (1 : Matrix I I k) M U
          0 (1 : Matrix I I k) 0
          0 X (1 : Matrix I I k) =
          centralBlock (k := k) Z * coupledKernelFrom (k := k) X := by
      simpa [Z, hU] using
        (Block3_unipotent_eq_centralBlock_mul_coupledKernelFrom (k := k) (M := M) (X := X))
    exact hEq ▸ hprod

/-- Equivalently, the full Magma unipotent cell fixes the mixed pair exactly when it
splits as a trace-zero central factor times the concrete `G_a^4` quotient family. -/
theorem fixesPair_Block3_unipotent_iff_exists_centralBlock_mul_coupledKernelFrom
    (M U X : Matrix I I k) :
    FixesPairBivector
      (Block3
        (1 : Matrix I I k) M U
        0 (1 : Matrix I I k) 0
        0 X (1 : Matrix I I k)) ↔
      ∃ Z : Matrix I I k,
        Z 0 0 + Z 1 1 = 0 ∧
        Block3
          (1 : Matrix I I k) M U
          0 (1 : Matrix I I k) 0
          0 X (1 : Matrix I I k) =
          centralBlock (k := k) Z * coupledKernelFrom (k := k) X := by
  constructor
  · intro hfix
    rcases (fixesPair_Block3_unipotent_iff (k := k) (M := M) (U := U) (X := X)).1 hfix with
      ⟨hU, hTrace⟩
    let Z : Matrix I I k := M - coupledCenterFrom (k := k) X
    refine ⟨Z, ?_, ?_⟩
    · have hsum : M 0 0 + M 1 1 + X.det = 0 := by
        exact eq_neg_iff_add_eq_zero.mp hTrace
      have hrewrite :
          Z 0 0 + Z 1 1 = M 0 0 + M 1 1 + X.det := by
        simp [Z, coupledCenterFrom, Matrix.det_fin_two, sub_eq_add_neg]
        ring
      rw [hrewrite, hsum]
    · simpa [Z, hU] using
        (Block3_unipotent_eq_centralBlock_mul_coupledKernelFrom (k := k) (M := M) (X := X))
  · rintro ⟨Z, hZ, hEq⟩
    have hZfix : FixesPairBivector (centralBlock (k := k) Z) := by
      exact (centralBlock_pointwise_iff (k := k) Z).2 hZ
    have hKfix : FixesPairBivector (coupledKernelFrom (k := k) X) :=
      coupledKernelFrom_pointwise (k := k) X
    have hprod : FixesPairBivector (centralBlock (k := k) Z * coupledKernelFrom (k := k) X) := by
      constructor
      · calc
          ActBivector (rep₁ (k := k))
              (centralBlock (k := k) Z * coupledKernelFrom (k := k) X) =
              ActBivector
                (ActBivector (rep₁ (k := k)) (coupledKernelFrom (k := k) X))
                (centralBlock (k := k) Z) := by
                  rw [actBivector_mul]
          _ = ActBivector (rep₁ (k := k)) (centralBlock (k := k) Z) := by rw [hKfix.1]
          _ = rep₁ (k := k) := hZfix.1
      · calc
          ActBivector (rep₂ (k := k))
              (centralBlock (k := k) Z * coupledKernelFrom (k := k) X) =
              ActBivector
                (ActBivector (rep₂ (k := k)) (coupledKernelFrom (k := k) X))
                (centralBlock (k := k) Z) := by
                  rw [actBivector_mul]
          _ = ActBivector (rep₂ (k := k)) (centralBlock (k := k) Z) := by rw [hKfix.2]
          _ = rep₂ (k := k) := hZfix.2
    simpa [hEq] using hprod

/-- Inside the natural coupled-times-Levi cell for the mixed one-point orbit, pointwise
fixing is equivalent to the explicit coupled block equations. The Levi factor already
fixes the pair pointwise, so the only remaining conditions are those on the coupled
off-diagonal blocks. -/
theorem fixesPair_coupledLevi_iff
    (Q : Matrix W I k)
    (R : Matrix I W k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    FixesPairBivector
      ((Matrix.fromBlocks 1 Q R 1) * (Block3 A B 0 0 A 0 0 0 H)) ↔
      Q * N4.J * Qᵀ = (0 : Matrix W W k) ∧
        N4.onePointRep₁ * Rᵀ + Q * N4.J = (0 : Matrix W I k) ∧
        R * N4.onePointRep₁ + N4.J * Qᵀ = (0 : Matrix I W k) ∧
        R * N4.onePointRep₁ * Rᵀ = (0 : Matrix I I k) ∧
        R * N4.onePointRep₂ = (0 : Matrix I W k) ∧
        N4.onePointRep₂ * Rᵀ = (0 : Matrix W I k) := by
  let gU : Matrix V V k := Matrix.fromBlocks 1 Q R 1
  let gL : Matrix V V k := Block3 A B 0 0 A 0 0 0 H
  have hL := pointwise_levi_family (k := k) A B H hA hB hH
  constructor
  · intro h
    have hU : FixesPairBivector gU := by
      constructor
      · calc
          ActBivector rep₁ gU = ActBivector (ActBivector rep₁ gL) gU := by
            rw [hL.1]
          _ = ActBivector rep₁ (gU * gL) := by
            rw [actBivector_mul]
          _ = rep₁ := h.1
      · calc
          ActBivector rep₂ gU = ActBivector (ActBivector rep₂ gL) gU := by
            rw [hL.2]
          _ = ActBivector rep₂ (gU * gL) := by
            rw [actBivector_mul]
          _ = rep₂ := h.2
    simpa [gU] using (fixesPair_identityOffDiagonal_iff (k := k) Q R).1 hU
  · intro hQR
    have hU : FixesPairBivector gU := by
      simpa [gU] using (fixesPair_identityOffDiagonal_iff (k := k) Q R).2 hQR
    constructor
    · calc
        ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
          rw [actBivector_mul]
        _ = ActBivector rep₁ gU := by rw [hL.1]
        _ = rep₁ := hU.1
    · calc
        ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
          rw [actBivector_mul]
        _ = ActBivector rep₂ gU := by rw [hL.2]
        _ = rep₂ := hU.2

/-- On the natural coupled-from-times-Levi cell, pointwise fixing is exactly the
rank-drop condition on the attached `2 × 2` coupled block. -/
theorem fixesPair_coupledFromLevi_iff_det_zero
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    FixesPairBivector
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H)) ↔
      X.det = 0 := by
  constructor
  · intro hfix
    have hQR :=
      (fixesPair_coupledLevi_iff
        (k := k)
        (Q := coupledQFrom (k := k) X)
        (R := coupledRFrom (k := k) X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 hfix
    have hU :
        FixesPairBivector
          (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
      exact
        (fixesPair_identityOffDiagonal_iff
          (k := k)
          (Q := coupledQFrom (k := k) X)
          (R := coupledRFrom (k := k) X)).2 hQR
    exact (coupledFrom_pointwise_iff_det_zero (k := k) X).1 hU
  · intro hX
    have hU :
        FixesPairBivector
          (Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) := by
      exact (coupledFrom_pointwise_iff_det_zero (k := k) X).2 hX
    have hQR :=
      (fixesPair_identityOffDiagonal_iff
        (k := k)
        (Q := coupledQFrom (k := k) X)
        (R := coupledRFrom (k := k) X)).1 hU
    exact
      (fixesPair_coupledLevi_iff
        (k := k)
        (Q := coupledQFrom (k := k) X)
        (R := coupledRFrom (k := k) X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 hQR

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to the
coupled block lying in the explicit rank-one family `X = u vᵀ`. -/
theorem fixesPair_coupledFromLevi_iff_exists_outer
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    FixesPairBivector
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H)) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  constructor
  · intro hfix
    exact exists_outer_of_det_zero
      (k := k)
      X
      ((fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 hfix)
  · rintro ⟨u, v, rfl⟩
    simpa using coupledOuterLevi_product_pointwise (k := k) u v A B H hA hB hH

/-- The explicit mixed pointwise subgroup has determinant `1`. -/
theorem coupledLevi_product_det
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let gU : Matrix V V k :=
      (Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
      (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)
    let gL : Matrix V V k := Block3 A B 0 0 A 0 0 0 H
    Matrix.det (gU * gL) = 1 := by
  intro gU gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using coupledE12E21_product_det (k := k) s t
  have hLdet : Matrix.det gL = 1 := by
    simpa [gL] using pointwise_levi_det (k := k) A B H hA hH
  rw [Matrix.det_mul, hUdet, hLdet]
  simp

/-- The diagonal explicit mixed pointwise subgroup has determinant `1`. -/
theorem coupledDiagonalLevi_product_det
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    Matrix.det gU = 1 := by
  intro g11 g22 gU
  have hUdet : Matrix.det (g11 * g22) = 1 := by
    simpa [g11, g22] using coupledE11E22_product_det (k := k) s t
  have hLdet : Matrix.det (Block3 A B 0 0 A 0 0 0 H) = 1 := by
    simpa using pointwise_levi_det (k := k) A B H hA hH
  rw [show gU = (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H) by rfl, Matrix.det_mul, hUdet, hLdet]
  simp

/-- The full explicit four-parameter mixed pointwise subgroup has determinant `1`. -/
theorem coupledExplicitLevi_product_det
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    Matrix.det gU = 1 := by
  intro g12 g21 g11 g22 gU
  have hUdet : Matrix.det ((g12 * g21) * (g11 * g22)) = 1 := by
    simpa [g12, g21, g11, g22] using coupledExplicit_product_det (k := k) s₁ t₁ s₂ t₂
  have hLdet : Matrix.det (Block3 A B 0 0 A 0 0 0 H) = 1 := by
    simpa using pointwise_levi_det (k := k) A B H hA hH
  rw [show gU = ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H) by rfl, Matrix.det_mul, hUdet, hLdet]
  simp

/-- The determinant on the mixed one-point line factors as the square of `a^3`. -/
theorem det_in_basis
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) =
      (a * a * a) * (a * a * a) := by
  rw [basis_linearCombination (k := k) (a := a) (b := b)]
  rw [Matrix.det_fromBlocks_zero₂₁]
  have htop := N4Summary.onePoint_det_in_basis (k := k) (a := a) (b := b)
  have hcard : Fintype.card I = 2 := by
    simp [I, N4.I]
  have hbot : Matrix.det (a • (N4.J (k := k))) = a * a := by
    simp [Matrix.det_smul, hcard, N4.J_det, pow_two]
  rw [htop, hbot]
  ring

/-- Rank drop on the mixed one-point line occurs exactly at the repeated-support
projective point. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • (rep₁ (k := k)) + b • (rep₂ (k := k))) = 0 ↔
      a = 0 := by
  rw [det_in_basis (k := k) (a := a) (b := b)]
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

/-- Embedding the `n = 4` one-point Borel action and a compatible lower-right scaling
gives a quotient lift for the mixed one-point orbit. -/
theorem borel_lift_action
    (a b : k)
    (ha : a ≠ 0) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix I I k := !![a, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks H 0 0 E
    ActBivector rep₁ g = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ g = (a * a) • rep₂ := by
  dsimp
  let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
  let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
  let E0 : Matrix I I k := !![a, 0; 0, 1]
  let g0 : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
  have htop := N4Summary.onePoint_borel_lift_action (k := k) (a := a) (b := b) ha
  have hbot : E0 * N4.J * E0ᵀ =
      a • N4.J := by
    dsimp [E0]
    rw [N4.mul_J_transpose_mul]
    simp [N4.J, Matrix.det_fin_two]
  have hblock1 :
      ActBivector rep₁ g0 =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₁ (k := k)) H0)
          0
          0
          (E0 * N4.J * E0ᵀ) := by
    simpa [g0, H0, E0, rep₁] using
      act_blockDiagonal
        (k := k)
        (ΩW := N4.onePointRep₁ (k := k))
        (ΩI := N4.J (k := k))
        (H := H0)
        (E := E0)
  have hblock2 :
      ActBivector rep₂ g0 =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₂ (k := k)) H0)
          0
          0
          (0 : Matrix I I k) := by
    simpa [g0, H0, E0, rep₂] using
      act_blockDiagonal
        (k := k)
        (ΩW := N4.onePointRep₂ (k := k))
        (ΩI := (0 : Matrix I I k))
        (H := H0)
        (E := E0)
  constructor
  · calc
      ActBivector rep₁ g0 =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₁ (k := k)) H0)
          0
          0
          (E0 * N4.J * E0ᵀ) := by
            simpa using hblock1
      _ =
        Matrix.fromBlocks
          (a • (N4.onePointRep₁ (k := k)) + b • (N4.onePointRep₂ (k := k)))
          0
          0
          (a • N4.J (k := k)) := by
            rw [htop.1, hbot]
      _ = a • (rep₁ (k := k)) + b • (rep₂ (k := k)) := by
            simpa using embed_linearCombination (k := k) (a := a) (b := b)
  · calc
      ActBivector rep₂ g0 =
        Matrix.fromBlocks
          (N4.ActBivector (N4.onePointRep₂ (k := k)) H0)
          0
          0
          (0 : Matrix I I k) := by
            simpa using hblock2
      _ =
        Matrix.fromBlocks ((a * a) • (N4.onePointRep₂ (k := k))) 0 0 (0 : Matrix I I k) := by
          simpa [H0, B0] using congrArg (fun X => Matrix.fromBlocks X 0 0 (0 : Matrix I I k)) htop.2
      _ = (a * a) • (rep₂ (k := k)) := by
          simpa using embed_second (k := k) (c := a * a)

/-- The determinant of the standard mixed one-point lift is `a^3`. -/
theorem borel_lift_det
    (a b : k) :
    let B : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B
    let E : Matrix I I k := !![a, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks H 0 0 E
    Matrix.det g = a * a * a := by
  dsimp
  rw [Matrix.det_fromBlocks_zero₂₁, N4Summary.onePoint_borel_lift_det]
  simp [Matrix.det_fin_two]

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing together with the
standard Borel lift gives the expected quotient action. -/
theorem coupledFromLevi_borel_product_lift_action
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (hX : X.det = 0)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
  intro gU B0 H0 E0 gL
  have hU : FixesPairBivector gU := by
    simpa [gU] using
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 hX
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • rep₁ + b • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector rep₁ gU + b • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ + b • rep₂ := by
        rw [hU.1, hU.2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector ((a * a) • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.2
      _ = (a * a) • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
      _ = (a * a) • rep₂ := by
        rw [hU.2]

/-- On the intrinsic coupled-from-times-Levi cell, the standard Borel product gives the
expected quotient action exactly when the coupled block satisfies `det X = 0`. -/
theorem coupledFromLevi_borel_product_lift_action_iff_det_zero
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    (ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ↔
      X.det = 0 := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hAct
    have hL := borel_lift_action (k := k) (a := a) (b := b) ha
    have hmul2 :
        ActBivector (ActBivector rep₂ gL) gU =
          ActBivector rep₂ (gU * gL) := by
      simpa using
        (actBivector_mul (k := k) (Ω := rep₂ (k := k)) (g := gU) (h := gL)).symm
    have hmul1 :
        ActBivector (ActBivector rep₁ gL) gU =
          ActBivector rep₁ (gU * gL) := by
      simpa using
        (actBivector_mul (k := k) (Ω := rep₁ (k := k)) (g := gU) (h := gL)).symm
    have haa : a * a ≠ 0 := mul_ne_zero ha ha
    have hrep2smul : (a * a) • ActBivector rep₂ gU = (a * a) • rep₂ := by
      calc
        (a * a) • ActBivector rep₂ gU =
            ActBivector ((a * a) • rep₂) gU := by
              simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
        _ = ActBivector (ActBivector rep₂ gL) gU := by
              rw [hL.2]
        _ = ActBivector rep₂ (gU * gL) := hmul2
        _ = (a * a) • rep₂ := hAct.2
    have hrep2 : ActBivector rep₂ gU = rep₂ := by
      exact (smul_right_injective _ haa) hrep2smul
    have hrep1lin :
        a • ActBivector rep₁ gU + b • ActBivector rep₂ gU =
          a • rep₁ + b • rep₂ := by
      calc
        a • ActBivector rep₁ gU + b • ActBivector rep₂ gU =
            ActBivector (a • rep₁ + b • rep₂) gU := by
              simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
        _ = ActBivector (ActBivector rep₁ gL) gU := by
              rw [hL.1]
        _ = ActBivector rep₁ (gU * gL) := hmul1
        _ = a • rep₁ + b • rep₂ := hAct.1
    have hrep1lin' :
        a • ActBivector rep₁ gU + b • rep₂ =
          a • rep₁ + b • rep₂ := by
      simpa [hrep2] using hrep1lin
    have hrep1smul : a • ActBivector rep₁ gU = a • rep₁ := by
      exact add_right_cancel hrep1lin'
    have hrep1 : ActBivector rep₁ gU = rep₁ := by
      exact (smul_right_injective _ ha) hrep1smul
    exact
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 ⟨hrep1, hrep2⟩
  · intro hX
    exact
      coupledFromLevi_borel_product_lift_action
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        hX
        a
        b
        ha

/-- On the intrinsic coupled-from-times-Levi cell, the standard Borel product gives the
expected quotient action exactly when the coupled block lies in the explicit
rank-one family `X = u vᵀ`. -/
theorem coupledFromLevi_borel_product_lift_action_iff_exists_outer
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    (ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hAct
    exact exists_outer_of_det_zero
      (k := k)
      X
      ((coupledFromLevi_borel_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hAct)
  · rintro ⟨u, v, rfl⟩
    have hX : Matrix.det (fun i j => u i * v j : Matrix I I k) = 0 := by
      simp [Matrix.det_fin_two]
      ring_nf
    simpa using
      coupledFromLevi_borel_product_lift_action
        (k := k)
        (X := fun i j => u i * v j)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        hX
        a
        b
        ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to getting
the usual quotient action after left multiplication by the standard mixed Borel lift. -/
theorem fixesPair_coupledFromLevi_iff_left_borel_action
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    FixesPairBivector gU ↔
      (ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hfix
    have hX :=
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 hfix
    exact
      (coupledFromLevi_borel_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hAct
    have hX :=
      (coupledFromLevi_borel_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hAct
    exact
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 hX

/-- Left-multiplying the standard mixed Borel lift by an intrinsic
coupled-from-times-Levi element does not change the determinant. -/
theorem coupledFromLevi_borel_product_lift_det
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a := by
  intro gU B0 H0 E0 gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using coupledFromLevi_product_det (k := k) X A B H hA hH
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL, E0, H0, B0] using borel_lift_det (k := k) (a := a) (b := b)
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- On the concrete rank-one intrinsic coupled-from-times-Levi family, left-multiplying
the standard Borel lift gives the expected quotient action. -/
theorem coupledOuterLevi_borel_product_lift_action
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix I I k := fun i j => u i * v j
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
  intro X gU B0 H0 E0 gL
  have hX : X.det = 0 := by
    simp [X, Matrix.det_fin_two]
    ring
  simpa [X] using
    coupledFromLevi_borel_product_lift_action
      (k := k) X A B H hA hB hH hX a b ha

/-- On the concrete rank-one intrinsic coupled-from-times-Levi family, left-multiplying
the standard Borel lift keeps determinant `a^3`. -/
theorem coupledOuterLevi_borel_product_lift_det
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix I I k := fun i j => u i * v j
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a := by
  intro X gU B0 H0 E0 gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [X, gU] using coupledOuterLevi_product_det (k := k) u v A B H hA hH
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL, E0, H0, B0] using borel_lift_det (k := k) (a := a) (b := b)
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- The concrete rank-one intrinsic coupled-from-times-Levi family gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem coupledOuterLevi_left_borel_package
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix I I k := fun i j => u i * v j
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro X gU B0 H0 E0 gL
  exact
    ⟨(coupledOuterLevi_borel_product_lift_action
        (k := k)
        (u := u)
        (v := v)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (coupledOuterLevi_borel_product_lift_det
        (k := k)
        (u := u)
        (v := v)
        (A := A)
        (B := B)
        (H := H)
        hA
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- The explicit pointwise subgroup on the mixed one-point orbit combines with the
standard Borel lift exactly as expected. -/
theorem coupledLevi_borel_product_lift_action
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) *
      (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
  intro gU B0 H0 E0 gL
  have hU := coupledLevi_product_pointwise (k := k) s t A B H hA hB hH
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • rep₁ + b • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector rep₁ gU + b • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ + b • rep₂ := by
        rw [hU.1, hU.2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector ((a * a) • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.2
      _ = (a * a) • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
      _ = (a * a) • rep₂ := by
        rw [hU.2]

/-- The diagonal explicit mixed coupled subgroup combines with the standard Borel
lift exactly as expected. -/
theorem coupledDiagonalLevi_borel_product_lift_action
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
  intro g11 g22 gU B0 H0 E0 gL
  have hU := coupledDiagonalLevi_product_pointwise (k := k) s t A B H hA hB hH
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • rep₁ + b • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector rep₁ gU + b • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ + b • rep₂ := by
        rw [hU.1, hU.2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector ((a * a) • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.2
      _ = (a * a) • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
      _ = (a * a) • rep₂ := by
        rw [hU.2]

/-- The full explicit four-parameter mixed pointwise subgroup combines with the
standard Borel lift exactly as expected. -/
theorem coupledExplicitLevi_borel_product_lift_action
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
  intro g12 g21 g11 g22 gU B0 H0 E0 gL
  have hU := coupledExplicitLevi_product_pointwise (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  · calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • rep₁ + b • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector rep₁ gU + b • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • rep₁ + b • rep₂ := by
        rw [hU.1, hU.2]
  · calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector ((a * a) • rep₂) gU := by
        simpa [gL, E0, H0, B0] using congrArg (fun M => ActBivector M gU) hL.2
      _ = (a * a) • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
      _ = (a * a) • rep₂ := by
        rw [hU.2]

/-- Left-multiplying the mixed Borel lift by the explicit mixed pointwise subgroup does
not change the determinant. -/
theorem coupledLevi_borel_product_lift_det
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) *
      (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a := by
  intro gU B0 H0 E0 gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using coupledLevi_product_det (k := k) s t A B H hA hB hH
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL, E0, H0, B0] using borel_lift_det (k := k) (a := a) (b := b)
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- The explicit off-diagonal mixed coupled-plus-Levi subgroup gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem coupledLevi_left_borel_package
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) *
      (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gU B0 H0 E0 gL
  exact
    ⟨(coupledLevi_borel_product_lift_action
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (coupledLevi_borel_product_lift_det
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- Left-multiplying the mixed Borel lift by the diagonal explicit mixed pointwise
subgroup does not change the determinant. -/
theorem coupledDiagonalLevi_borel_product_lift_det
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a := by
  intro g11 g22 gU B0 H0 E0 gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU, g11, g22] using coupledDiagonalLevi_product_det (k := k) s t A B H hA hB hH
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL, E0, H0, B0] using borel_lift_det (k := k) (a := a) (b := b)
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- The diagonal mixed coupled-plus-Levi subgroup gives the expected left Borel
quotient action together with determinant `a^3`. -/
theorem coupledDiagonalLevi_left_borel_package
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro g11 g22 gU B0 H0 E0 gL
  exact
    ⟨(coupledDiagonalLevi_borel_product_lift_action
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (coupledDiagonalLevi_borel_product_lift_det
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- Left-multiplying the mixed Borel lift by the full explicit four-parameter pointwise
subgroup does not change the determinant. -/
theorem coupledExplicitLevi_borel_product_lift_det
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    Matrix.det (gU * gL) = a * a * a := by
  intro g12 g21 g11 g22 gU B0 H0 E0 gL
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU, g12, g21, g11, g22] using
      coupledExplicitLevi_product_det (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL, E0, H0, B0] using borel_lift_det (k := k) (a := a) (b := b)
  rw [Matrix.det_mul, hUdet, hLdet]
  ring

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem coupledExplicitLevi_left_borel_package
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gL g12 g21 g11 g22 gU
  exact
    ⟨(coupledExplicitLevi_borel_product_lift_action
        (k := k)
        (s₁ := s₁)
        (t₁ := t₁)
        (s₂ := s₂)
        (t₂ := t₂)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (coupledExplicitLevi_borel_product_lift_det
        (k := k)
        (s₁ := s₁)
        (t₁ := t₁)
        (s₂ := s₂)
        (t₂ := t₂)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- Right-multiplying the standard Borel lift by the explicit mixed pointwise subgroup
does not change the quotient action on the chosen basis. -/
theorem borel_coupledLevi_right_product_lift_action
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) *
      (Block3 A B 0 0 A 0 0 0 H)
    ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
  intro gL gU
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  have hU := coupledLevi_product_pointwise (k := k) s t A B H hA hB hH
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
      _ = (a * a) • rep₂ := by
        simpa [gL] using hL.2

/-- Right-multiplying the standard Borel lift by the diagonal explicit mixed coupled
subgroup does not change the quotient action on the chosen basis. -/
theorem borel_coupledDiagonalLevi_right_product_lift_action
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
  intro gL g11 g22 gU
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  have hU := coupledDiagonalLevi_product_pointwise (k := k) s t A B H hA hB hH
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
      _ = (a * a) • rep₂ := by
        simpa [gL] using hL.2

/-- Right-multiplying the standard Borel lift by the full explicit four-parameter mixed
pointwise subgroup does not change the quotient action on the chosen basis. -/
theorem borel_coupledExplicitLevi_right_product_lift_action
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
  intro gL g12 g21 g11 g22 gU
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  have hU := coupledExplicitLevi_product_pointwise (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH
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
      _ = (a * a) • rep₂ := by
        simpa [gL] using hL.2

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected
quotient action on both the left and right of the standard Borel lift. -/
theorem coupledExplicitLevi_left_and_right_borel_action
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
      (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂)) := by
  intro gL g12 g21 g11 g22 gU
  constructor
  · simpa [gL, g12, g21, g11, g22, gU] using
      coupledExplicitLevi_borel_product_lift_action
        (k := k)
        (s₁ := s₁)
        (t₁ := t₁)
        (s₂ := s₂)
        (t₂ := t₂)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha
  · simpa [gL, g12, g21, g11, g22, gU] using
      borel_coupledExplicitLevi_right_product_lift_action
        (k := k)
        (s₁ := s₁)
        (t₁ := t₁)
        (s₂ := s₂)
        (t₂ := t₂)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha

/-- Right-multiplying the mixed Borel lift by the explicit mixed pointwise subgroup does
not change the determinant. -/
theorem borel_coupledLevi_right_product_lift_det
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) *
      (Block3 A B 0 0 A 0 0 0 H)
    Matrix.det (gL * gU) = a * a * a := by
  intro gL gU
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using borel_lift_det (k := k) (a := a) (b := b)
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using coupledLevi_product_det (k := k) s t A B H hA hB hH
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

/-- The explicit off-diagonal mixed coupled-plus-Levi subgroup gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem coupledLevi_right_borel_package
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) *
      (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL gU
  exact
    ⟨(borel_coupledLevi_right_product_lift_action
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (borel_coupledLevi_right_product_lift_det
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- Right-multiplying the mixed Borel lift by the diagonal explicit mixed coupled
subgroup does not change the determinant. -/
theorem borel_coupledDiagonalLevi_right_product_lift_det
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    Matrix.det (gL * gU) = a * a * a := by
  intro gL g11 g22 gU
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using borel_lift_det (k := k) (a := a) (b := b)
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU, g11, g22] using coupledDiagonalLevi_product_det (k := k) s t A B H hA hB hH
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

/-- The diagonal mixed coupled-plus-Levi subgroup gives the expected right Borel
quotient action together with determinant `a^3`. -/
theorem coupledDiagonalLevi_right_borel_package
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL g11 g22 gU
  exact
    ⟨(borel_coupledDiagonalLevi_right_product_lift_action
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (borel_coupledDiagonalLevi_right_product_lift_det
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- The explicit off-diagonal mixed coupled-plus-Levi subgroup gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem coupledLevi_left_and_right_borel_package
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1) *
       (Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1)) *
      (Block3 A B 0 0 A 0 0 0 H)
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) := by
  intro gL gU
  exact
    ⟨(coupledLevi_left_borel_package
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (coupledLevi_right_borel_package
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- The diagonal mixed coupled-plus-Levi subgroup gives the expected quotient action
on both sides of the standard Borel lift, and both products have determinant `a^3`. -/
theorem coupledDiagonalLevi_left_and_right_borel_package
    (s t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := (g11 * g22) * (Block3 A B 0 0 A 0 0 0 H)
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) := by
  intro gL g11 g22 gU
  exact
    ⟨(coupledDiagonalLevi_left_borel_package
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (coupledDiagonalLevi_right_borel_package
        (k := k)
        (s := s)
        (t := t)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- The single mixed `E12` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem coupledE12Levi_left_and_right_borel_package
    (s : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1
    let gU : Matrix V V k := g12 * (Block3 A B 0 0 A 0 0 0 H)
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) := by
  intro gL g12 gU
  simpa [gL, gU, g12, coupledE21_zero_eq_one (k := k), Matrix.mul_assoc] using
    (coupledLevi_left_and_right_borel_package
      (k := k) (s := s) (t := 0) A B H hA hB hH a b ha)

/-- The single mixed `E12` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem coupledE12Levi_left_borel_package
    (s : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1
    let gU : Matrix V V k := g12 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gL g12 gU
  exact (coupledE12Levi_left_and_right_borel_package
    (k := k) s A B H hA hB hH a b ha).1

/-- The single mixed `E12` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem coupledE12Levi_right_borel_package
    (s : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s) (coupledR12 (k := k) s) 1
    let gU : Matrix V V k := g12 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL g12 gU
  exact (coupledE12Levi_left_and_right_borel_package
    (k := k) s A B H hA hB hH a b ha).2

/-- The single mixed `E21` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem coupledE21Levi_left_and_right_borel_package
    (t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1
    let gU : Matrix V V k := g21 * (Block3 A B 0 0 A 0 0 0 H)
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) := by
  intro gL g21 gU
  simpa [gL, gU, g21, coupledE12_zero_eq_one (k := k), Matrix.mul_assoc] using
    (coupledLevi_left_and_right_borel_package
      (k := k) (s := 0) (t := t) A B H hA hB hH a b ha)

/-- The single mixed `E21` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem coupledE21Levi_left_borel_package
    (t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1
    let gU : Matrix V V k := g21 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gL g21 gU
  exact (coupledE21Levi_left_and_right_borel_package
    (k := k) t A B H hA hB hH a b ha).1

/-- The single mixed `E21` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem coupledE21Levi_right_borel_package
    (t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t) (coupledR21 (k := k) t) 1
    let gU : Matrix V V k := g21 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL g21 gU
  exact (coupledE21Levi_left_and_right_borel_package
    (k := k) t A B H hA hB hH a b ha).2

/-- The single mixed `E11` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem coupledE11Levi_left_and_right_borel_package
    (s : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let gU : Matrix V V k := g11 * (Block3 A B 0 0 A 0 0 0 H)
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) := by
  intro gL g11 gU
  simpa [gL, gU, g11, coupledE22_zero_eq_one (k := k), Matrix.mul_assoc] using
    (coupledDiagonalLevi_left_and_right_borel_package
      (k := k) (s := s) (t := 0) A B H hA hB hH a b ha)

/-- The single mixed `E11` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem coupledE11Levi_left_borel_package
    (s : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let gU : Matrix V V k := g11 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gL g11 gU
  exact (coupledE11Levi_left_and_right_borel_package
    (k := k) s A B H hA hB hH a b ha).1

/-- The single mixed `E11` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem coupledE11Levi_right_borel_package
    (s : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s) (coupledR11 (k := k) s) 1
    let gU : Matrix V V k := g11 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL g11 gU
  exact (coupledE11Levi_left_and_right_borel_package
    (k := k) s A B H hA hB hH a b ha).2

/-- The single mixed `E22` generator, combined with the Levi family, gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem coupledE22Levi_left_and_right_borel_package
    (t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := g22 * (Block3 A B 0 0 A 0 0 0 H)
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) := by
  intro gL g22 gU
  simpa [gL, gU, g22, coupledE11_zero_eq_one (k := k), Matrix.mul_assoc] using
    (coupledDiagonalLevi_left_and_right_borel_package
      (k := k) (s := 0) (t := t) A B H hA hB hH a b ha)

/-- The single mixed `E22` generator, combined with the Levi family, gives the expected
left Borel quotient action together with determinant `a^3`. -/
theorem coupledE22Levi_left_borel_package
    (t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := g22 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gL g22 gU
  exact (coupledE22Levi_left_and_right_borel_package
    (k := k) t A B H hA hB hH a b ha).1

/-- The single mixed `E22` generator, combined with the Levi family, gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem coupledE22Levi_right_borel_package
    (t : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t) (coupledR22 (k := k) t) 1
    let gU : Matrix V V k := g22 * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL g22 gU
  exact (coupledE22Levi_left_and_right_borel_package
    (k := k) t A B H hA hB hH a b ha).2

/-- Right-multiplying the mixed Borel lift by the full explicit four-parameter pointwise
subgroup does not change the determinant. -/
theorem borel_coupledExplicitLevi_right_product_lift_det
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    Matrix.det (gL * gU) = a * a * a := by
  intro gL g12 g21 g11 g22 gU
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using borel_lift_det (k := k) (a := a) (b := b)
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU, g12, g21, g11, g22] using
      coupledExplicitLevi_product_det (k := k) s₁ t₁ s₂ t₂ A B H hA hB hH
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem coupledExplicitLevi_right_borel_package
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL g12 g21 g11 g22 gU
  exact
    ⟨(borel_coupledExplicitLevi_right_product_lift_action
        (k := k)
        (s₁ := s₁)
        (t₁ := t₁)
        (s₂ := s₂)
        (t₂ := t₂)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (borel_coupledExplicitLevi_right_product_lift_det
        (k := k)
        (s₁ := s₁)
        (t₁ := t₁)
        (s₂ := s₂)
        (t₂ := t₂)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- The full explicit four-parameter mixed pointwise subgroup gives the expected quotient
action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem coupledExplicitLevi_left_and_right_borel_package
    (s₁ t₁ s₂ t₂ : k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let g12 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ12 (k := k) s₁) (coupledR12 (k := k) s₁) 1
    let g21 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ21 (k := k) t₁) (coupledR21 (k := k) t₁) 1
    let g11 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ11 (k := k) s₂) (coupledR11 (k := k) s₂) 1
    let g22 : Matrix V V k :=
      Matrix.fromBlocks 1 (coupledQ22 (k := k) t₂) (coupledR22 (k := k) t₂) 1
    let gU : Matrix V V k := ((g12 * g21) * (g11 * g22)) * (Block3 A B 0 0 A 0 0 0 H)
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL g12 g21 g11 g22 gU
  constructor
  · exact
      ⟨(coupledExplicitLevi_borel_product_lift_action
          (k := k)
          (s₁ := s₁)
          (t₁ := t₁)
          (s₂ := s₂)
          (t₂ := t₂)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          (a := a)
          (b := b)
          ha),
        (coupledExplicitLevi_borel_product_lift_det
          (k := k)
          (s₁ := s₁)
          (t₁ := t₁)
          (s₂ := s₂)
          (t₂ := t₂)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          (a := a)
          (b := b)
          ha)⟩
  · exact
      ⟨(borel_coupledExplicitLevi_right_product_lift_action
          (k := k)
          (s₁ := s₁)
          (t₁ := t₁)
          (s₂ := s₂)
          (t₂ := t₂)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          (a := a)
          (b := b)
          ha),
        (borel_coupledExplicitLevi_right_product_lift_det
          (k := k)
          (s₁ := s₁)
          (t₁ := t₁)
          (s₂ := s₂)
          (t₂ := t₂)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          (a := a)
          (b := b)
          ha)⟩

/-- Right-multiplying the standard Borel lift by an intrinsic
coupled-from-times-Levi element with `det X = 0` does not change the quotient action. -/
theorem borel_coupledFromLevi_right_product_lift_action
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (hX : X.det = 0)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
  intro gL gU
  have hL := borel_lift_action (k := k) (a := a) (b := b) ha
  have hU : FixesPairBivector gU := by
    simpa [gU] using
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 hX
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
      _ = (a * a) • rep₂ := by
        simpa [gL] using hL.2

/-- On the intrinsic coupled-from-times-Levi cell, right-multiplying the standard Borel
lift gives the expected quotient action exactly when the coupled block satisfies
`det X = 0`. -/
theorem borel_coupledFromLevi_right_product_lift_action_iff_det_zero
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ↔
      X.det = 0 := by
  intro gL gU
  constructor
  · intro hAct
    have hL := borel_lift_action (k := k) (a := a) (b := b) ha
    have hdet : Matrix.det gL = a * a * a := by
      simpa [gL] using borel_lift_det (k := k) (a := a) (b := b)
    have hunit : IsUnit (Matrix.det gL) := by
      rw [hdet]
      exact isUnit_iff_ne_zero.mpr (mul_ne_zero (mul_ne_zero ha ha) ha)
    have hinj := actBivector_injective_of_isUnit_det (k := k) gL hunit
    have hmul1 :
        ActBivector (ActBivector rep₁ gU) gL =
          ActBivector rep₁ (gL * gU) := by
      simpa using
        (actBivector_mul (k := k) (Ω := rep₁ (k := k)) (g := gL) (h := gU)).symm
    have hmul2 :
        ActBivector (ActBivector rep₂ gU) gL =
          ActBivector rep₂ (gL * gU) := by
      simpa using
        (actBivector_mul (k := k) (Ω := rep₂ (k := k)) (g := gL) (h := gU)).symm
    have hrep1 : ActBivector rep₁ gU = rep₁ := by
      apply hinj
      calc
        ActBivector (ActBivector rep₁ gU) gL = ActBivector rep₁ (gL * gU) := hmul1
        _ = a • rep₁ + b • rep₂ := hAct.1
        _ = ActBivector rep₁ gL := by
          simpa [gL] using hL.1.symm
    have hrep2 : ActBivector rep₂ gU = rep₂ := by
      apply hinj
      calc
        ActBivector (ActBivector rep₂ gU) gL = ActBivector rep₂ (gL * gU) := hmul2
        _ = (a * a) • rep₂ := hAct.2
        _ = ActBivector rep₂ gL := by
          simpa [gL] using hL.2.symm
    exact
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 ⟨hrep1, hrep2⟩
  · intro hX
    exact
      borel_coupledFromLevi_right_product_lift_action
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        hX
        a
        b
        ha

/-- On the intrinsic coupled-from-times-Levi cell, right-multiplying the standard Borel
lift gives the expected quotient action exactly when the coupled block lies in the
explicit rank-one family `X = u vᵀ`. -/
theorem borel_coupledFromLevi_right_product_lift_action_iff_exists_outer
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  intro gL gU
  constructor
  · intro hAct
    exact exists_outer_of_det_zero
      (k := k)
      X
      ((borel_coupledFromLevi_right_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hAct)
  · rintro ⟨u, v, rfl⟩
    have hX : Matrix.det (fun i j => u i * v j : Matrix I I k) = 0 := by
      simp [Matrix.det_fin_two]
      ring_nf
    simpa using
      borel_coupledFromLevi_right_product_lift_action
        (k := k)
        (X := fun i j => u i * v j)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        hX
        a
        b
        ha

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to getting
the usual quotient action after right multiplication by the standard mixed Borel lift. -/
theorem fixesPair_coupledFromLevi_iff_right_borel_action
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    FixesPairBivector gU ↔
      (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂) := by
  intro gL gU
  constructor
  · intro hfix
    have hX :=
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 hfix
    exact
      (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hAct
    have hX :=
      (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hAct
    exact
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 hX

/-- Right-multiplying the standard mixed Borel lift by an intrinsic
coupled-from-times-Levi element does not change the determinant. -/
theorem coupledFromLevi_left_borel_action_iff_right_borel_action
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ↔
      (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂) := by
  intro gL gU
  constructor
  · intro hAct
    exact
      (fixesPair_coupledFromLevi_iff_right_borel_action
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 <|
        (fixesPair_coupledFromLevi_iff_left_borel_action
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).2 hAct
  · intro hAct
    exact
      (fixesPair_coupledFromLevi_iff_left_borel_action
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 <|
        (fixesPair_coupledFromLevi_iff_right_borel_action
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).2 hAct

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to
simultaneously getting the usual left and right quotient actions. -/
theorem fixesPair_coupledFromLevi_iff_left_and_right_borel_action
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    FixesPairBivector gU ↔
      ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂)) := by
  intro gL gU
  constructor
  · intro hfix
    exact
      ⟨(fixesPair_coupledFromLevi_iff_left_borel_action
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hfix,
        (fixesPair_coupledFromLevi_iff_right_borel_action
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hfix⟩
  · rintro ⟨hLeft, _hRight⟩
    exact
      (fixesPair_coupledFromLevi_iff_left_borel_action
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hLeft

/-- On the intrinsic coupled-from-times-Levi cell, simultaneously getting the usual left
and right quotient actions is equivalent to the coupled block satisfying `det X = 0`. -/
theorem coupledFromLevi_left_and_right_borel_action_iff_det_zero
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
      (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂)) ↔
      X.det = 0 := by
  intro gL gU
  constructor
  · rintro ⟨hLeft, _hRight⟩
    exact
      (coupledFromLevi_borel_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hLeft
  · intro hX
    have hLeft :
        ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
      exact
        (coupledFromLevi_borel_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).2 hX
    have hRight :
        ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
      exact
        (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).2 hX
    exact ⟨hLeft, hRight⟩

/-- On the intrinsic coupled-from-times-Levi cell, simultaneously getting the usual left
and right quotient actions is equivalent to the coupled block being rank one. -/
theorem coupledFromLevi_left_and_right_borel_action_iff_exists_outer
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
      (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂)) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  intro gL gU
  constructor
  · rintro ⟨hLeft, _hRight⟩
    exact
      (coupledFromLevi_borel_product_lift_action_iff_exists_outer
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hLeft
  · rintro ⟨u, v, huv⟩
    have hLeft :
        ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂ := by
      exact
        (coupledFromLevi_borel_product_lift_action_iff_exists_outer
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).2 ⟨u, v, huv⟩
    have hRight :
        ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
      exact
        (borel_coupledFromLevi_right_product_lift_action_iff_exists_outer
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).2 ⟨u, v, huv⟩
    exact ⟨hLeft, hRight⟩

/-- Right-multiplying the standard mixed Borel lift by an intrinsic
coupled-from-times-Levi element does not change the determinant. -/
theorem borel_coupledFromLevi_right_product_lift_det
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    Matrix.det (gL * gU) = a * a * a := by
  intro gL gU
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using borel_lift_det (k := k) (a := a) (b := b)
  have hUdet : Matrix.det gU = 1 := by
    simpa [gU] using coupledFromLevi_product_det (k := k) X A B H hA hH
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient action
with the matching determinant condition is equivalent to `det X = 0`. -/
theorem coupledFromLevi_left_borel_package_iff_det_zero
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      X.det = 0 := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hPkg
    exact
      (coupledFromLevi_borel_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hPkg.1
  · intro hX
    exact
      ⟨(coupledFromLevi_borel_product_lift_action
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          hX
          a
          b
          ha),
        (coupledFromLevi_borel_product_lift_det
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hH
          a
          b
          ha)⟩

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient action
with the matching determinant condition is equivalent to explicit rank-one
factorization `X = u vᵀ`. -/
theorem coupledFromLevi_left_borel_package_iff_exists_outer
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hPkg
    exact
      exists_outer_of_det_zero
        (k := k)
        X
        ((coupledFromLevi_left_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hPkg)
  · rintro ⟨u, v, rfl⟩
    have hX : Matrix.det (fun i j => u i * v j : Matrix I I k) = 0 := by
      simp [Matrix.det_fin_two]
      ring_nf
    exact
      (coupledFromLevi_left_borel_package_iff_det_zero
        (k := k)
        (X := fun i j => u i * v j)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient action
with the matching determinant condition is equivalent to `det X = 0`. -/
theorem coupledFromLevi_right_borel_package_iff_det_zero
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) ↔
      X.det = 0 := by
  intro gL gU
  constructor
  · intro hPkg
    exact
      (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hPkg.1
  · intro hX
    exact
      ⟨(borel_coupledFromLevi_right_product_lift_action
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          hX
          a
          b
          ha),
        (borel_coupledFromLevi_right_product_lift_det
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hH
          a
          b
          ha)⟩

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient action
with the matching determinant condition is equivalent to explicit rank-one
factorization `X = u vᵀ`. -/
theorem coupledFromLevi_right_borel_package_iff_exists_outer
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  intro gL gU
  constructor
  · intro hPkg
    exact
      exists_outer_of_det_zero
        (k := k)
        X
        ((coupledFromLevi_right_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hPkg)
  · rintro ⟨u, v, rfl⟩
    have hX : Matrix.det (fun i j => u i * v j : Matrix I I k) = 0 := by
      simp [Matrix.det_fin_two]
      ring_nf
    exact
      (coupledFromLevi_right_borel_package_iff_det_zero
        (k := k)
        (X := fun i j => u i * v j)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX

/-- On the intrinsic coupled-from-times-Levi cell, the usual left and right Borel
products give the expected quotient action with determinant `a^3` on both sides
exactly when `det X = 0`. -/
theorem coupledFromLevi_left_and_right_borel_package_iff_det_zero
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) ↔
      X.det = 0 := by
  intro gL gU
  constructor
  · intro hPkg
    exact
      (coupledFromLevi_left_and_right_borel_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 ⟨hPkg.1.1, hPkg.2.1⟩
  · intro hX
    constructor
    · exact
        ⟨(coupledFromLevi_borel_product_lift_action
            (k := k)
            (X := X)
            (A := A)
            (B := B)
            (H := H)
            hA
            hB
            hH
            hX
            a
            b
            ha),
          (coupledFromLevi_borel_product_lift_det
            (k := k)
            (X := X)
            (A := A)
            (B := B)
            (H := H)
            hA
            hH
            a
            b
            ha)⟩
    · exact
        ⟨(borel_coupledFromLevi_right_product_lift_action
            (k := k)
            (X := X)
            (A := A)
            (B := B)
            (H := H)
            hA
            hB
            hH
            hX
            a
            b
            ha),
          (borel_coupledFromLevi_right_product_lift_det
            (k := k)
            (X := X)
            (A := A)
            (B := B)
            (H := H)
            hA
            hH
            a
            b
            ha)⟩

/-- On the intrinsic coupled-from-times-Levi cell, the full two-sided Borel
action-and-determinant package is equivalent to explicit rank-one factorization
`X = u vᵀ`. -/
theorem coupledFromLevi_left_and_right_borel_package_iff_exists_outer
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) ↔
      ∃ u v : I → k, X = fun i j => u i * v j := by
  intro gL gU
  constructor
  · intro hPkg
    exact
      exists_outer_of_det_zero
        (k := k)
        X
        ((coupledFromLevi_left_and_right_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hPkg)
  · rintro ⟨u, v, rfl⟩
    have hX : Matrix.det (fun i j => u i * v j : Matrix I I k) = 0 := by
      simp [Matrix.det_fin_two]
      ring_nf
    exact
      (coupledFromLevi_left_and_right_borel_package_iff_det_zero
        (k := k)
        (X := fun i j => u i * v j)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX

/-- On the intrinsic coupled-from-times-Levi cell, the usual left one-sided
Borel action-and-determinant package is equivalent to the full two-sided
package. -/
theorem coupledFromLevi_left_borel_package_iff_left_and_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) := by
  intro gL gU
  constructor
  · intro hLeftPkg
    have hX :
        X.det = 0 := by
      exact
        (coupledFromLevi_left_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hLeftPkg
    exact
      (coupledFromLevi_left_and_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hPkg
    exact hPkg.1

/-- On the intrinsic coupled-from-times-Levi cell, the usual right one-sided
Borel action-and-determinant package is equivalent to the full two-sided
package. -/
theorem coupledFromLevi_right_borel_package_iff_left_and_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
      Matrix.det (gL * gU) = a * a * a) ↔
      (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) := by
  intro gL gU
  constructor
  · intro hRightPkg
    have hX :
        X.det = 0 := by
      exact
        (coupledFromLevi_right_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hRightPkg
    exact
      (coupledFromLevi_left_and_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hPkg
    exact hPkg.2

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to the
full two-sided Borel action-and-determinant package. -/
theorem fixesPair_coupledFromLevi_iff_left_and_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    FixesPairBivector gU ↔
      (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) := by
  intro gL gU
  constructor
  · intro hFix
    exact
      (coupledFromLevi_left_and_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 <|
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 hFix
  · intro hPkg
    exact
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 <|
      (coupledFromLevi_left_and_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hPkg

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to
the usual left Borel quotient action together with its determinant condition. -/
theorem fixesPair_coupledFromLevi_iff_left_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    FixesPairBivector gU ↔
      ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hFix
    exact
      (coupledFromLevi_left_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 <|
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 hFix
  · intro hPkg
    exact
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 <|
      (coupledFromLevi_left_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hPkg

/-- On the intrinsic coupled-from-times-Levi cell, pointwise fixing is equivalent to
the usual right Borel quotient action together with its determinant condition. -/
theorem fixesPair_coupledFromLevi_iff_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    FixesPairBivector gU ↔
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL gU
  constructor
  · intro hFix
    exact
      (coupledFromLevi_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 <|
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).1 hFix
  · intro hPkg
    exact
      (fixesPair_coupledFromLevi_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH).2 <|
      (coupledFromLevi_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).1 hPkg

/-- On the intrinsic coupled-from-times-Levi cell, the usual left and right Borel
action-and-determinant packages are equivalent. -/
theorem coupledFromLevi_left_borel_package_iff_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
        ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
      Matrix.det (gU * gL) = a * a * a) ↔
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hLeft
    have hX :
        X.det = 0 := by
      exact
        (coupledFromLevi_left_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hLeft
    exact
      (coupledFromLevi_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hRight
    have hX :
        X.det = 0 := by
      exact
        (coupledFromLevi_right_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hRight
    exact
      (coupledFromLevi_left_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX

/-- On the intrinsic coupled-from-times-Levi cell, getting the usual left Borel
quotient action is equivalent to getting the full two-sided Borel
action-and-determinant package. -/
theorem coupledFromLevi_left_borel_action_iff_left_and_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ↔
      (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) := by
  intro gL gU
  constructor
  · intro hLeft
    have hX :
        X.det = 0 := by
      exact
        (coupledFromLevi_borel_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hLeft
    exact
      (coupledFromLevi_left_and_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hPkg
    exact hPkg.1.1

/-- On the intrinsic coupled-from-times-Levi cell, getting the usual right Borel
quotient action is equivalent to getting the full two-sided Borel
action-and-determinant package. -/
theorem coupledFromLevi_right_borel_action_iff_left_and_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ↔
      (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
          Matrix.det (gU * gL) = a * a * a) ∧
        ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
            ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
          Matrix.det (gL * gU) = a * a * a)) := by
  intro gL gU
  constructor
  · intro hRight
    have hX :
        X.det = 0 := by
      exact
        (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hRight
    exact
      (coupledFromLevi_left_and_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hPkg
    exact hPkg.2.1

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient
action is equivalent to the matching one-sided Borel action-and-determinant
package. -/
theorem coupledFromLevi_left_borel_action_iff_left_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    let B0 : Matrix N4.I N4.I k := !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
    let H0 : Matrix W W k := N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k) B0
    let E0 : Matrix I I k := !![a, 0; 0, 1]
    let gL : Matrix V V k := Matrix.fromBlocks H0 0 0 E0
    (ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ↔
      ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gU B0 H0 E0 gL
  constructor
  · intro hLeft
    have hX :
        X.det = 0 := by
      exact
        (coupledFromLevi_borel_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hLeft
    exact
      (coupledFromLevi_left_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hPkg
    exact hPkg.1

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient
action is equivalent to the matching one-sided Borel action-and-determinant
package. -/
theorem coupledFromLevi_right_borel_action_iff_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ↔
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL gU
  constructor
  · intro hRight
    have hX :
        X.det = 0 := by
      exact
        (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hRight
    exact
      (coupledFromLevi_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hPkg
    exact hPkg.1

/-- On the intrinsic coupled-from-times-Levi cell, the usual left Borel quotient
action is equivalent to the opposite one-sided Borel
action-and-determinant package. -/
theorem coupledFromLevi_left_borel_action_iff_right_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ↔
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro gL gU
  constructor
  · intro hLeft
    have hX : X.det = 0 := by
      exact
        (coupledFromLevi_borel_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hLeft
    exact
      (coupledFromLevi_right_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hRightPkg
    have hX : X.det = 0 := by
      exact
        (coupledFromLevi_right_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hRightPkg
    exact
      (coupledFromLevi_borel_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX

/-- On the intrinsic coupled-from-times-Levi cell, the usual right Borel quotient
action is equivalent to the opposite one-sided Borel
action-and-determinant package. -/
theorem coupledFromLevi_right_borel_action_iff_left_borel_package
    (X : Matrix I I k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ↔
      ((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) := by
  intro gL gU
  constructor
  · intro hRight
    have hX : X.det = 0 := by
      exact
        (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hRight
    exact
      (coupledFromLevi_left_borel_package_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX
  · intro hLeftPkg
    have hX : X.det = 0 := by
      exact
        (coupledFromLevi_left_borel_package_iff_det_zero
          (k := k)
          (X := X)
          (A := A)
          (B := B)
          (H := H)
          hA
          hB
          hH
          a
          b
          ha).1 hLeftPkg
    exact
      (borel_coupledFromLevi_right_product_lift_action_iff_det_zero
        (k := k)
        (X := X)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        a
        b
        ha).2 hX

/-- On the concrete rank-one intrinsic coupled-from-times-Levi family, right-multiplying
the standard Borel lift gives the expected quotient action. -/
theorem borel_coupledOuterLevi_right_product_lift_action
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix I I k := fun i j => u i * v j
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ (gL * gU) = (a * a) • rep₂ := by
  intro X gL gU
  have hX : X.det = 0 := by
    simp [X, Matrix.det_fin_two]
    ring
  simpa [X] using
    borel_coupledFromLevi_right_product_lift_action
      (k := k) X A B H hA hB hH hX a b ha

/-- On the concrete rank-one intrinsic coupled-from-times-Levi family, right-multiplying
the standard Borel lift keeps determinant `a^3`. -/
theorem borel_coupledOuterLevi_right_product_lift_det
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix I I k := fun i j => u i * v j
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    Matrix.det (gL * gU) = a * a * a := by
  intro X gL gU
  have hLdet : Matrix.det gL = a * a * a := by
    simpa [gL] using borel_lift_det (k := k) (a := a) (b := b)
  have hUdet : Matrix.det gU = 1 := by
    simpa [X, gU] using coupledOuterLevi_product_det (k := k) u v A B H hA hH
  rw [Matrix.det_mul, hLdet, hUdet]
  ring

/-- The concrete rank-one intrinsic coupled-from-times-Levi family gives the expected
right Borel quotient action together with determinant `a^3`. -/
theorem coupledOuterLevi_right_borel_package
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix I I k := fun i j => u i * v j
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a) := by
  intro X gL gU
  exact
    ⟨(borel_coupledOuterLevi_right_product_lift_action
        (k := k)
        (u := u)
        (v := v)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (borel_coupledOuterLevi_right_product_lift_det
        (k := k)
        (u := u)
        (v := v)
        (A := A)
        (B := B)
        (H := H)
        hA
        hH
        (a := a)
        (b := b)
        ha)⟩

/-- The concrete rank-one intrinsic coupled-from-times-Levi family gives the expected
quotient action on both sides of the standard Borel lift, and both products have
determinant `a^3`. -/
theorem coupledOuterLevi_left_and_right_borel_package
    (u v : I → k)
    (A B H : Matrix I I k)
    (hA : A.det = 1)
    (hB :
      A * N4.J * Bᵀ + B * N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let X : Matrix I I k := fun i j => u i * v j
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N4.onePointScale (k := k) a * N4.onePointUpperShear (k := k)
          (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix N4.I N4.I k))
        0
        0
        (!![a, 0; 0, 1] : Matrix I I k)
    let gU : Matrix V V k :=
      ((Matrix.fromBlocks 1 (coupledQFrom (k := k) X) (coupledRFrom (k := k) X) 1) *
        (Block3 A B 0 0 A 0 0 0 H))
    (((ActBivector rep₁ (gU * gL) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gU * gL) = (a * a) • rep₂) ∧
        Matrix.det (gU * gL) = a * a * a) ∧
      ((ActBivector rep₁ (gL * gU) = a • rep₁ + b • rep₂ ∧
          ActBivector rep₂ (gL * gU) = (a * a) • rep₂) ∧
        Matrix.det (gL * gU) = a * a * a)) := by
  intro X gL gU
  exact
    ⟨(coupledOuterLevi_left_borel_package
        (k := k)
        (u := u)
        (v := v)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha),
      (coupledOuterLevi_right_borel_package
        (k := k)
        (u := u)
        (v := v)
        (A := A)
        (B := B)
        (H := H)
        hA
        hB
        hH
        (a := a)
        (b := b)
        ha)⟩

end N6MixedOnePoint
end Wedge2Formalization
