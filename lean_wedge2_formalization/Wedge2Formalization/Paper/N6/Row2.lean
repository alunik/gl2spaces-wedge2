import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N4
import Wedge2Formalization.Paper.N5
import Wedge2Formalization.N4PaperSummary
import Wedge2Formalization.N6Summary
import Wedge2Formalization.N6WeightedTwoPoint
import Wedge2Formalization.N6DoublePureSingular
import Mathlib.Algebra.Polynomial.Roots
import Mathlib.LinearAlgebra.Matrix.Rank

namespace Wedge2Formalization
namespace Paper
namespace N6

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 6`, row 2.
Representative `⟨e₁∧e₃ + e₂∧e₄ + e₅∧e₆, e₁∧e₄⟩`.
Divisor `3[a]`.
Claimed stabilizer:
`K_L = U_7 \rtimes (SL_2(k) \times SL_2(k))`, exact quotient family `Q_L = B`.
-/
namespace Row2

def rep₁ :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 2:
`e₁∧e₃ + e₂∧e₄ + e₅∧e₆`. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))
    0
    0
    (Wedge2Formalization.N4.J (k := k))

/-- Literal second basis vector from Appendix A, row 2:
`e₁∧e₄`. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))
    0
    0
    0

/-- Explicit basis change sending the literal paper representative to the internal
working representative pair. -/
def paperChange :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))
    0
    0
    (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)

/-- Transport of the literal paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6MixedOnePoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  calc
    Wedge2Formalization.N6MixedOnePoint.ActBivector
        (paperRep₁ (k := k)) (paperChange (k := k)) =
      Matrix.fromBlocks
        (Wedge2Formalization.N4.ActBivector
          (Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))
          (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k)))
        0
        0
        ((1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k) *
          Wedge2Formalization.N4.J (k := k) *
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)ᵀ) := by
          simpa [paperRep₁, paperChange] using
            (Wedge2Formalization.N6MixedOnePoint.act_blockDiagonal
              (k := k)
              (ΩW := Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))
              (ΩI := Wedge2Formalization.N4.J (k := k))
              (H := Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))
              (E := (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
                Wedge2Formalization.N6MixedOnePoint.I k)))
    _ =
      Matrix.fromBlocks
        (Wedge2Formalization.N4.onePointRep₁ (k := k))
        0
        0
        (Wedge2Formalization.N4.J (k := k)) := by
          simpa [Wedge2Formalization.Paper.N4.Row1.rep₁] using
            congrArg
              (fun Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k =>
                Matrix.fromBlocks
                  Ω
                  (0 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N6MixedOnePoint.I k)
                  (0 : Matrix Wedge2Formalization.N6MixedOnePoint.I Wedge2Formalization.N4.V k)
                  (Wedge2Formalization.N4.J (k := k)))
              (Wedge2Formalization.Paper.N4.Row1.paperRep₁_transport (k := k))
    _ = rep₁ (k := k) := by
          simp [rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁]

/-- Transport of the literal paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6MixedOnePoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  calc
    Wedge2Formalization.N6MixedOnePoint.ActBivector
        (paperRep₂ (k := k)) (paperChange (k := k)) =
      Matrix.fromBlocks
        (Wedge2Formalization.N4.ActBivector
          (Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))
          (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k)))
        0
        0
        0 := by
          simpa [paperRep₂, paperChange] using
            (Wedge2Formalization.N6MixedOnePoint.act_blockDiagonal
              (k := k)
              (ΩW := Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))
              (ΩI := (0 : Matrix Wedge2Formalization.N6MixedOnePoint.I
                Wedge2Formalization.N6MixedOnePoint.I k))
              (H := Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))
              (E := (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
                Wedge2Formalization.N6MixedOnePoint.I k)))
    _ =
      Matrix.fromBlocks
        (Wedge2Formalization.N4.onePointRep₂ (k := k))
        0
        0
        0 := by
          simpa [Wedge2Formalization.Paper.N4.Row1.rep₂] using
            congrArg
              (fun Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k =>
                Matrix.fromBlocks
                  Ω
                  (0 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N6MixedOnePoint.I k)
                  (0 : Matrix Wedge2Formalization.N6MixedOnePoint.I Wedge2Formalization.N4.V k)
                  (0 : Matrix Wedge2Formalization.N6MixedOnePoint.I
                    Wedge2Formalization.N6MixedOnePoint.I k))
              (Wedge2Formalization.Paper.N4.Row1.paperRep₂_transport (k := k))
    _ = rep₂ (k := k) := by
          simp [rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂]

def U
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X

def Levi
    (A B H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Wedge2Formalization.N6MixedOnePoint.Block3 A B 0 0 A 0 0 0 H

/-- The displayed pointwise kernel family on the mixed `3[a]+[b]` row, written in the
paper form `U_7 ⋊ (SL₂(k) × SL₂(k))`. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6MixedOnePoint.V
        Wedge2Formalization.N6MixedOnePoint.V k) :=
  { g | ∃ X A B H, A.det = 1 ∧
      A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
      H.det = 1 ∧
      g = U (k := k) X * Levi (k := k) A B H }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

def lift
    (a b : k) :
    Matrix Wedge2Formalization.N6MixedOnePoint.V
      Wedge2Formalization.N6MixedOnePoint.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N4.onePointScale (k := k) a *
      Wedge2Formalization.N4.onePointUpperShear (k := k)
        (!![b * (a⁻¹ * a⁻¹), 0; 0, 0] : Matrix Wedge2Formalization.N4.I
          Wedge2Formalization.N4.I k))
    0
    0
    (!![a, 0; 0, 1] : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)

theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • rep₁ (k := k) + b • rep₂ (k := k)) = 0 ↔
      a = 0 :=
  Wedge2Formalization.N6Summary.mixedOnePoint_det_zero_iff
    (k := k) (a := a) (b := b)

theorem U_pointwise
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
      (U (k := k) X) := by
  simpa [U] using
    Wedge2Formalization.N6Summary.mixedOnePoint_coupledKernelFrom_pointwise
      (k := k) X

theorem Levi_pointwise
    (A B H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
      (Levi (k := k) A B H) := by
  simpa [Levi] using
    Wedge2Formalization.N6Summary.mixedOnePoint_pointwise_levi_family
      (k := k) (A := A) (B := B) (H := H) hA hB hH

theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6MixedOnePoint.V
        Wedge2Formalization.N6MixedOnePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ X A B H, A.det = 1 ∧
        A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
        H.det = 1 ∧
        g = U (k := k) X * Levi (k := k) A B H := by
  rfl

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_7 \rtimes (SL_2(k) \times SL_2(k))`. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6MixedOnePoint.V
        Wedge2Formalization.N6MixedOnePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ X A B H, A.det = 1 ∧
        A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
        H.det = 1 ∧
        g = U (k := k) X * Levi (k := k) A B H :=
  mem_K_iff (k := k) (g := g)

theorem pointwise_of_mem_K
    (g :
      Matrix Wedge2Formalization.N6MixedOnePoint.V
        Wedge2Formalization.N6MixedOnePoint.V k)
    (hg : g ∈ K (k := k)) :
    Wedge2Formalization.N6MixedOnePoint.FixesPairBivector g := by
  rcases (mem_K_iff (k := k) (g := g)).1 hg with
    ⟨X, A, B, H, hA, hB, hH, rfl⟩
  have hU := U_pointwise (k := k) X
  have hL := Levi_pointwise (k := k) A B H hA hB hH
  constructor
  ·
    calc
      Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k))
          (U (k := k) X * Levi (k := k) A B H) =
        Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.N6MixedOnePoint.ActBivector
            (Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k))
            (Levi (k := k) A B H))
          (U (k := k) X) := by
            rw [Wedge2Formalization.N6MixedOnePoint.actBivector_mul]
      _ =
        Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k))
          (U (k := k) X) := by
            rw [hL.1]
      _ = Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k) := hU.1
  ·
    calc
      Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k))
          (U (k := k) X * Levi (k := k) A B H) =
        Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.N6MixedOnePoint.ActBivector
            (Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k))
            (Levi (k := k) A B H))
          (U (k := k) X) := by
            rw [Wedge2Formalization.N6MixedOnePoint.actBivector_mul]
      _ =
        Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k))
          (U (k := k) X) := by
            rw [hL.2]
      _ = Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k) := hU.2

/-- Exact pointwise criterion on the intrinsic coupled-times-Levi cell for the mixed
repeated-support row. -/
theorem cell_pointwise_iff_exists_outer
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (A B H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * Wedge2Formalization.N4.J * Bᵀ +
        B * Wedge2Formalization.N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
      ((Matrix.fromBlocks
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.W
            Wedge2Formalization.N6MixedOnePoint.W k)
          (Wedge2Formalization.N6MixedOnePoint.coupledQFrom (k := k) X)
          (Wedge2Formalization.N6MixedOnePoint.coupledRFrom (k := k) X)
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)) *
        Levi (k := k) A B H) ↔
      ∃ u v : Wedge2Formalization.N6MixedOnePoint.I → k,
        X = fun i j => u i * v j :=
  Wedge2Formalization.N6Summary.mixedOnePoint_coupledFromLevi_iff_exists_outer
    (k := k) (X := X) (A := A) (B := B) (H := H) hA hB hH

private theorem coupledFrom_Block3_eq
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Matrix.fromBlocks
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.W
          Wedge2Formalization.N6MixedOnePoint.W k)
        (Wedge2Formalization.N6MixedOnePoint.coupledQFrom (k := k) X)
        (Wedge2Formalization.N6MixedOnePoint.coupledRFrom (k := k) X)
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
          Wedge2Formalization.N6MixedOnePoint.I k) =
      Wedge2Formalization.N6MixedOnePoint.Block3
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
          Wedge2Formalization.N6MixedOnePoint.I k)
        0
        (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)
        0
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
          Wedge2Formalization.N6MixedOnePoint.I k)
        0
        0
        X
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
          Wedge2Formalization.N6MixedOnePoint.I k) := by
  have hQ :
      Wedge2Formalization.N4.J (k := k) * Xᵀ * Wedge2Formalization.N4.J (k := k) =
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
        Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Matrix.transpose_apply] <;> ring
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    rcases i with i | i <;> rcases j with j | j <;>
      fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.coupledQFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledRFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, hQ,
        Wedge2Formalization.N4.J, Matrix.fromBlocks]
  ·
    rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.coupledQFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledRFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, hQ,
        Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Matrix.vecMul,
        dotProduct, Fin.sum_univ_two, Matrix.transpose_apply] <;> ring
  ·
    fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.coupledQFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledRFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, hQ,
        Wedge2Formalization.N4.J, Matrix.fromBlocks]
  ·
    fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.coupledQFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledRFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, hQ,
        Wedge2Formalization.N4.J, Matrix.fromBlocks]

private theorem centralBlock_mul_coupledKernelFrom_comm
    (Z X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Wedge2Formalization.N6MixedOnePoint.centralBlock (k := k) Z *
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X =
      Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X *
        Wedge2Formalization.N6MixedOnePoint.centralBlock (k := k) Z := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    rcases i with i | i <;> rcases j with j | j <;>
      fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring
  ·
    rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring
  ·
    fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring
  ·
    fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3,
        Matrix.mul_apply, Fin.sum_univ_two] <;> ring

private theorem centralBlock_mul_Levi
    (Z A B H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Wedge2Formalization.N6MixedOnePoint.centralBlock (k := k) Z *
        Levi (k := k) A B H =
      Levi (k := k) A (B + Z * A) H := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  ·
    rcases i with i | i <;> rcases j with j | j <;>
      fin_cases i <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.add_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  ·
    rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.add_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  ·
    fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.add_apply,
        Matrix.mul_apply, Fin.sum_univ_two]
  ·
    fin_cases i <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.centralBlock,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.add_apply,
        Matrix.mul_apply, Fin.sum_univ_two]

private theorem add_mul_preserves_levi_relation
    (Z A B : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * Wedge2Formalization.N4.J * Bᵀ +
        B * Wedge2Formalization.N4.J * Aᵀ = 0)
    (hZ : Z 0 0 + Z 1 1 = 0) :
    A * Wedge2Formalization.N4.J * (B + Z * A)ᵀ +
        (B + Z * A) * Wedge2Formalization.N4.J * Aᵀ = 0 := by
  have hAJAT : A * Wedge2Formalization.N4.J * Aᵀ = Wedge2Formalization.N4.J := by
    simpa [Wedge2Formalization.N4.J, hA] using
      Wedge2Formalization.N4.mul_J_transpose_mul (k := k) (A := A)
  have hZJ : Wedge2Formalization.N4.J * Zᵀ + Z * Wedge2Formalization.N4.J = 0 := by
    have hZ' : Z 1 1 = -Z 0 0 := by
      exact (eq_neg_iff_add_eq_zero).2 (by simpa [add_comm] using hZ)
    have hZ'' : Z 0 0 = -Z 1 1 := by
      exact (eq_neg_iff_add_eq_zero).2 hZ
    ext i j
    fin_cases i <;> fin_cases j
    ·
      simp [Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
        Wedge2Formalization.N4.J, Matrix.transpose_apply]
    ·
      simp [Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
        Wedge2Formalization.N4.J, Matrix.transpose_apply, hZ'']
    ·
      simp [Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
        Wedge2Formalization.N4.J, Matrix.transpose_apply, hZ']
    ·
      simp [Matrix.mul_apply, Matrix.vecMul, dotProduct, Fin.sum_univ_two,
        Wedge2Formalization.N4.J, Matrix.transpose_apply]
  calc
    A * Wedge2Formalization.N4.J * (B + Z * A)ᵀ +
        (B + Z * A) * Wedge2Formalization.N4.J * Aᵀ
        =
      (A * Wedge2Formalization.N4.J * Bᵀ +
        B * Wedge2Formalization.N4.J * Aᵀ) +
      ((A * Wedge2Formalization.N4.J * Aᵀ) * Zᵀ +
        Z * (A * Wedge2Formalization.N4.J * Aᵀ)) := by
          simp [Matrix.transpose_add, Matrix.transpose_mul, Matrix.mul_add, Matrix.add_mul,
            Matrix.mul_assoc, add_assoc, add_left_comm, add_comm]
    _ = (A * Wedge2Formalization.N4.J * Bᵀ +
          B * Wedge2Formalization.N4.J * Aᵀ) +
        (Wedge2Formalization.N4.J * Zᵀ + Z * Wedge2Formalization.N4.J) := by
          simp [hAJAT]
    _ = 0 := by rw [hB, hZJ]; simp

/-- Any pointwise element already lying on the intrinsic coupled-from-times-Levi cell belongs
to the displayed paper kernel family `K`. -/
theorem intrinsic_cell_mem_K
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (A B H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (hB :
      A * Wedge2Formalization.N4.J * Bᵀ +
        B * Wedge2Formalization.N4.J * Aᵀ = 0)
    (hH : H.det = 1)
    (houter : ∃ u v : Wedge2Formalization.N6MixedOnePoint.I → k,
      X = fun i j => u i * v j) :
    ((Matrix.fromBlocks
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.W
          Wedge2Formalization.N6MixedOnePoint.W k)
        (Wedge2Formalization.N6MixedOnePoint.coupledQFrom (k := k) X)
        (Wedge2Formalization.N6MixedOnePoint.coupledRFrom (k := k) X)
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
          Wedge2Formalization.N6MixedOnePoint.I k)) *
      Levi (k := k) A B H) ∈ K (k := k) := by
  have hfixBlock :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)
          0
          (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)
          0
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)
          0
          0
          X
        (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
          Wedge2Formalization.N6MixedOnePoint.I k)) := by
    have hdet0 : X.det = 0 := by
      rcases houter with ⟨u, v, rfl⟩
      simp [Matrix.det_fin_two, mul_assoc, mul_left_comm, mul_comm, sub_eq_add_neg,
        add_comm, add_left_comm, add_assoc]
    exact
      (Wedge2Formalization.N6MixedOnePoint.fixesPair_Block3_unipotent_iff
        (k := k)
        (M := 0)
        (U := Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)
        (X := X)).2
        (by
          constructor
          · rfl
          · simpa [hdet0] using hdet0)
  rcases
    (Wedge2Formalization.N6MixedOnePoint.fixesPair_Block3_unipotent_iff_exists_centralBlock_mul_coupledKernelFrom
      (k := k)
      (M := 0)
      (U := Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)
      (X := X)).1 hfixBlock with ⟨Z, hZ, hEq⟩
  refine (mem_K_iff (k := k) (g := _)).2 ?_
  refine ⟨X, A, B + Z * A, H, hA, ?_, hH, ?_⟩
  · exact add_mul_preserves_levi_relation (k := k) Z A B hA hB hZ
  ·
    calc
      ((Matrix.fromBlocks
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.W
            Wedge2Formalization.N6MixedOnePoint.W k)
          (Wedge2Formalization.N6MixedOnePoint.coupledQFrom (k := k) X)
          (Wedge2Formalization.N6MixedOnePoint.coupledRFrom (k := k) X)
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)) *
        Levi (k := k) A B H)
          =
        (Wedge2Formalization.N6MixedOnePoint.Block3
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)
          0
          (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)
          0
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)
          0
          0
          X
          (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
            Wedge2Formalization.N6MixedOnePoint.I k)) *
        Levi (k := k) A B H := by
            rw [coupledFrom_Block3_eq (k := k) (X := X)]
      _ =
        (Wedge2Formalization.N6MixedOnePoint.centralBlock (k := k) Z *
          Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X) *
        Levi (k := k) A B H := by
            rw [hEq]
      _ =
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X *
          (Wedge2Formalization.N6MixedOnePoint.centralBlock (k := k) Z *
            Levi (k := k) A B H) := by
              rw [← Matrix.mul_assoc,
                centralBlock_mul_coupledKernelFrom_comm (k := k) (Z := Z) (X := X),
                Matrix.mul_assoc]
      _ =
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X *
          Levi (k := k) A (B + Z * A) H := by
              rw [centralBlock_mul_Levi (k := k) (Z := Z) (A := A) (B := B) (H := H)]

private theorem coupledUpperFrom_mul_J
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
        Wedge2Formalization.N4.J (k := k) =
      -(Wedge2Formalization.N4.J (k := k) * Xᵀ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
      Wedge2Formalization.N4.J, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, Matrix.transpose_apply]

private theorem coupledCenterFrom_J_relation
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Wedge2Formalization.N4.J (k := k) *
        (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X)ᵀ +
      Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X *
        Wedge2Formalization.N4.J (k := k) =
      - (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
          Wedge2Formalization.N4.J (k := k) *
          (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)ᵀ) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    (simp [Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
      Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
      Wedge2Formalization.N4.J, Matrix.cons_val', Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_fin_one, vecHead, vecTail]
      ; ring_nf)

private theorem zero_of_mul_J_transpose_eq_zero_left
    (A X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (h : A * Wedge2Formalization.N4.J (k := k) * Xᵀ = 0) :
    X = 0 := by
  have hunit : IsUnit ((A * Wedge2Formalization.N4.J (k := k)).det) := by
    rw [Matrix.det_mul, hA, Wedge2Formalization.N4.J_det, one_mul]
    exact isUnit_one
  have hXt : Xᵀ = 0 := by
    calc
      Xᵀ = (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
              Wedge2Formalization.N6MixedOnePoint.I k) * Xᵀ := by simp
      _ =
          ((A * Wedge2Formalization.N4.J (k := k))⁻¹ *
              (A * Wedge2Formalization.N4.J (k := k))) * Xᵀ := by
            rw [Matrix.nonsing_inv_mul _ hunit]
      _ =
          (A * Wedge2Formalization.N4.J (k := k))⁻¹ *
            ((A * Wedge2Formalization.N4.J (k := k)) * Xᵀ) := by
              simp [Matrix.mul_assoc]
      _ = 0 := by
            have htmp := congrArg
              (fun M =>
                (A * Wedge2Formalization.N4.J (k := k))⁻¹ * M) h
            simpa [Matrix.mul_assoc] using htmp
  simpa [transpose_eq_zero] using hXt

private theorem zero_of_mul_J_transpose_eq_zero_right
    (X H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hH : H.det = 1)
    (h : X * Wedge2Formalization.N4.J (k := k) * Hᵀ = 0) :
    X = 0 := by
  have hunit : IsUnit ((Wedge2Formalization.N4.J (k := k) * Hᵀ).det) := by
    rw [Matrix.det_mul, Matrix.det_transpose, hH, Wedge2Formalization.N4.J_det, one_mul]
    exact isUnit_one
  calc
    X = X * (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k) := by simp
    _ =
        X * ((Wedge2Formalization.N4.J (k := k) * Hᵀ) *
          (Wedge2Formalization.N4.J (k := k) * Hᵀ)⁻¹) := by
            rw [Matrix.mul_nonsing_inv _ hunit]
    _ =
        (X * (Wedge2Formalization.N4.J (k := k) * Hᵀ)) *
          (Wedge2Formalization.N4.J (k := k) * Hᵀ)⁻¹ := by
            simp [Matrix.mul_assoc]
    _ = 0 := by
          have htmp := congrArg
            (fun M =>
              M * (Wedge2Formalization.N4.J (k := k) * Hᵀ)⁻¹) h
          simpa [Matrix.mul_assoc] using htmp

private theorem eq_of_mul_J_transpose_eq_J
    (A D : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hA : A.det = 1)
    (h : A * Wedge2Formalization.N4.J (k := k) * Dᵀ =
      Wedge2Formalization.N4.J (k := k)) :
    D = A := by
  have hAA :
      A * Wedge2Formalization.N4.J (k := k) * Aᵀ =
        Wedge2Formalization.N4.J (k := k) := by
    rw [Wedge2Formalization.N4.mul_J_transpose_mul, hA]
    simp [Wedge2Formalization.N4.J]
  have hsub :
      A * Wedge2Formalization.N4.J (k := k) * (D - A)ᵀ = 0 := by
    calc
      A * Wedge2Formalization.N4.J (k := k) * (D - A)ᵀ
          =
        A * Wedge2Formalization.N4.J (k := k) * Dᵀ -
          A * Wedge2Formalization.N4.J (k := k) * Aᵀ := by
            simp [Matrix.transpose_sub, Matrix.mul_assoc, Matrix.mul_sub]
      _ = 0 := by rw [h, hAA]; simp
  have hzero := zero_of_mul_J_transpose_eq_zero_left
    (k := k) A (D - A) hA hsub
  exact sub_eq_zero.mp hzero

private theorem coupledKernelFrom_mul_Levi
    (X A B H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X *
        Levi (k := k) A B H =
      Wedge2Formalization.N6MixedOnePoint.Block3
        A
        (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A + B)
        (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H)
        0
        A
        0
        0
        (X * A)
        H := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · rcases i with i | i <;> rcases j with j | j <;>
      fin_cases i <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.mul_apply, Matrix.vecMul,
        dotProduct, Fin.sum_univ_two, Matrix.add_apply, add_assoc, add_left_comm, add_comm] <;> ring_nf
  · rcases i with i | i <;> fin_cases i <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.mul_apply, Matrix.vecMul,
        dotProduct, Fin.sum_univ_two, Matrix.add_apply, add_assoc, add_left_comm, add_comm] <;> ring_nf
  · fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.mul_apply, Matrix.vecMul,
        dotProduct, Fin.sum_univ_two, Matrix.add_apply, add_assoc, add_left_comm, add_comm] <;> ring_nf
  · fin_cases i <;> fin_cases j <;>
      simp [Levi, Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom,
        Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom,
        Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.mul_apply, Matrix.vecMul,
        dotProduct, Fin.sum_univ_two, Matrix.add_apply, add_assoc, add_left_comm, add_comm] <;> ring_nf

private theorem row2_Afix
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H)) :
    A * Wedge2Formalization.N4.J (k := k) * Aᵀ =
      Wedge2Formalization.N4.J (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inl 0))) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inl 1))) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inl (Sum.inl 0))) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inl (Sum.inl 1))) hfix'.2

private theorem row2_Ceq
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H)) :
    A * Wedge2Formalization.N4.J (k := k) * Cᵀ = 0 := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inr 0))) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inr 1))) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inl (Sum.inr 0))) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inl (Sum.inr 1))) hfix'.2

private theorem row2_Peq
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H)) :
    A * Wedge2Formalization.N4.J (k := k) * Pᵀ = 0 := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inr 0)) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inr 1)) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inr 0)) hfix'.2
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₂, Wedge2Formalization.N6MixedOnePoint.rep₂,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Wedge2Formalization.N6MixedOnePoint.Block3,
      Matrix.fromBlocks, Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul,
      dotProduct, Fin.sum_univ_two] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inr 1)) hfix'.2

private theorem row2_Hfix
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H))
    (hC0 : C = 0) (hP0 : P = 0) :
    H * Wedge2Formalization.N4.J (k := k) * Hᵀ =
      Wedge2Formalization.N4.J (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inr 0) (Sum.inr 0)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inr 0) (Sum.inr 1)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inr 1) (Sum.inr 0)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inr 1) (Sum.inr 1)) hfix'.1

private theorem row2_Veq
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H))
    (hC0 : C = 0) (hP0 : P = 0) :
    V0 * Wedge2Formalization.N4.J (k := k) * Hᵀ = 0 := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inr 0)) (Sum.inr 0)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inr 0)) (Sum.inr 1)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inr 1)) (Sum.inr 0)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inr 1)) (Sum.inr 1)) hfix'.1

private theorem row2_Drel
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H))
    (hC0 : C = 0) (hP0 : P = 0) (hV0 : V0 = 0) :
    A * Wedge2Formalization.N4.J (k := k) * Dᵀ =
      Wedge2Formalization.N4.J (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0, hV0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inr 0))) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0, hV0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inr 1))) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0, hV0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inl (Sum.inr 0))) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0, hV0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inl (Sum.inr 1))) hfix'.1

private theorem row2_UYeq
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H))
    (hP0 : P = 0) :
    A * Wedge2Formalization.N4.J (k := k) * Yᵀ +
      U0 * Wedge2Formalization.N4.J (k := k) * Hᵀ = 0 := by
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inr 0)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inr 1)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inr 0)) hfix'.1
  ·
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hP0] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inr 1)) hfix'.1

set_option maxHeartbeats 20000000 in
private theorem row2_Top
    (A B U0 C D V0 P Y H : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (hfix' :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H))
    (hC0 : C = 0) (hP0 : P = 0) (hD : D = A) :
    A * Wedge2Formalization.N4.J (k := k) * Bᵀ +
      B * Wedge2Formalization.N4.J (k := k) * Aᵀ +
      U0 * Wedge2Formalization.N4.J (k := k) * U0ᵀ = 0 := by
  have h01 :
      (A * Wedge2Formalization.N4.J (k := k) * Bᵀ +
          B * Wedge2Formalization.N4.J (k := k) * Aᵀ +
          U0 * Wedge2Formalization.N4.J (k := k) * U0ᵀ) 0 1 = 0 := by
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0, hD,
      add_assoc, add_left_comm, add_comm, Matrix.cons_val', Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_fin_one, vecHead, vecTail] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inl 1))) hfix'.1
  have h10 :
      (A * Wedge2Formalization.N4.J (k := k) * Bᵀ +
          B * Wedge2Formalization.N4.J (k := k) * Aᵀ +
          U0 * Wedge2Formalization.N4.J (k := k) * U0ᵀ) 1 0 = 0 := by
    simpa [Wedge2Formalization.N6MixedOnePoint.ActBivector,
      rep₁, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.J,
      Wedge2Formalization.N6MixedOnePoint.Block3, Matrix.fromBlocks,
      Matrix.transpose_apply, Matrix.mul_apply, Matrix.vecMul, dotProduct,
      Fin.sum_univ_two, hC0, hP0, hD,
      add_assoc, add_left_comm, add_comm, Matrix.cons_val', Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.cons_val_fin_one, vecHead, vecTail] using
      congrArg (fun Ω => Ω (Sum.inl (Sum.inl 1)) (Sum.inl (Sum.inl 0))) hfix'.1
  ext i j
  fin_cases i <;> fin_cases j
  ·
    simp [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_two,
      Wedge2Formalization.N4.J]
    ring
  ·
    exact h01
  ·
    exact h10
  ·
    simp [Matrix.add_apply, Matrix.mul_apply, Fin.sum_univ_two,
      Wedge2Formalization.N4.J]
    ring

set_option maxHeartbeats 20000000 in
set_option maxRecDepth 4096 in
theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6MixedOnePoint.FixesPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hfix
    let A : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inl (Sum.inl i)) (Sum.inl (Sum.inl j))
    let B : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inl (Sum.inl i)) (Sum.inl (Sum.inr j))
    let U0 : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inl (Sum.inl i)) (Sum.inr j)
    let C : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inl (Sum.inr i)) (Sum.inl (Sum.inl j))
    let D : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inl (Sum.inr i)) (Sum.inl (Sum.inr j))
    let V0 : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inl (Sum.inr i)) (Sum.inr j)
    let P : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inr i) (Sum.inl (Sum.inl j))
    let Y : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inr i) (Sum.inl (Sum.inr j))
    let H : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      fun i j => g (Sum.inr i) (Sum.inr j)
    have hgBlock :
        g =
          Wedge2Formalization.N6MixedOnePoint.Block3
            A B U0 C D V0 P Y H := by
      ext i j
      rcases i with i | i <;> rcases j with j | j
      · rcases i with i | i <;> rcases j with j | j <;> rfl
      · rcases i with i | i <;> rfl
      · rcases j with j | j <;> rfl
      · rfl
    have hfix' :
        Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
          (Wedge2Formalization.N6MixedOnePoint.Block3
            A B U0 C D V0 P Y H) := by
      simpa [hgBlock] using hfix
    have hAfix := row2_Afix (k := k) A B U0 C D V0 P Y H hfix'
    have hA : A.det = 1 := by
      rw [Wedge2Formalization.N4.mul_J_transpose_mul] at hAfix
      have h01 := congrArg (fun M => M 0 1) hAfix
      simpa [Wedge2Formalization.N4.J] using h01
    have hCeq := row2_Ceq (k := k) A B U0 C D V0 P Y H hfix'
    have hC0 : C = 0 := zero_of_mul_J_transpose_eq_zero_left (k := k) A C hA hCeq
    have hPeq := row2_Peq (k := k) A B U0 C D V0 P Y H hfix'
    have hP0 : P = 0 := zero_of_mul_J_transpose_eq_zero_left (k := k) A P hA hPeq
    have hHfix := row2_Hfix (k := k) A B U0 C D V0 P Y H hfix' hC0 hP0
    have hH : H.det = 1 := by
      rw [Wedge2Formalization.N4.mul_J_transpose_mul] at hHfix
      have h01 := congrArg (fun M => M 0 1) hHfix
      simpa [Wedge2Formalization.N4.J] using h01
    have hVeq := row2_Veq (k := k) A B U0 C D V0 P Y H hfix' hC0 hP0
    have hV0 : V0 = 0 := zero_of_mul_J_transpose_eq_zero_right (k := k) V0 H hH hVeq
    have hDrel := row2_Drel (k := k) A B U0 C D V0 P Y H hfix' hC0 hP0 hV0
    have hD : D = A := eq_of_mul_J_transpose_eq_J (k := k) A D hA hDrel
    have hUYeq := row2_UYeq (k := k) A B U0 C D V0 P Y H hfix' hP0
    have hTop := row2_Top (k := k) A B U0 C D V0 P Y H hfix' hC0 hP0 hD
    have hAunit : IsUnit (A.det) := by
      rw [hA]
      exact isUnit_one
    let X : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k := Y * A⁻¹
    have hXA : X * A = Y := by
      dsimp [X]
      calc
        Y * A⁻¹ * A = Y * (A⁻¹ * A) := by simp [Matrix.mul_assoc]
        _ = Y * (1 : Matrix Wedge2Formalization.N6MixedOnePoint.I
              Wedge2Formalization.N6MixedOnePoint.I k) := by
              rw [Matrix.nonsing_inv_mul _ hAunit]
        _ = Y := by simp
    have hAJXt :
        A * Wedge2Formalization.N4.J (k := k) * Yᵀ =
          Wedge2Formalization.N4.J (k := k) * Xᵀ := by
      calc
        A * Wedge2Formalization.N4.J (k := k) * Yᵀ
            =
          A * Wedge2Formalization.N4.J (k := k) * (X * A)ᵀ := by
              rw [hXA]
        _ =
          A * Wedge2Formalization.N4.J (k := k) * (Aᵀ * Xᵀ) := by
              simp [Matrix.transpose_mul]
        _ =
          (A * Wedge2Formalization.N4.J (k := k) * Aᵀ) * Xᵀ := by
              simp [Matrix.mul_assoc]
        _ = Wedge2Formalization.N4.J (k := k) * Xᵀ := by
              rw [hAfix]
    have hUJ :
        U0 * (Wedge2Formalization.N4.J (k := k) * Hᵀ) =
          -(Wedge2Formalization.N4.J (k := k) * Xᵀ) := by
      have htmp :
          U0 * (Wedge2Formalization.N4.J (k := k) * Hᵀ) =
            -(A * Wedge2Formalization.N4.J (k := k) * Yᵀ) := by
        exact eq_neg_of_add_eq_zero_left (by simpa [Matrix.mul_assoc, add_comm] using hUYeq)
      simpa [Matrix.mul_assoc, hAJXt] using htmp
    have hunitJH :
        IsUnit ((Wedge2Formalization.N4.J (k := k) * Hᵀ).det) := by
      rw [Matrix.det_mul, Matrix.det_transpose, hH, Wedge2Formalization.N4.J_det, one_mul]
      exact isUnit_one
    have hUtarget :
        (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H) *
            (Wedge2Formalization.N4.J (k := k) * Hᵀ) =
          -(Wedge2Formalization.N4.J (k := k) * Xᵀ) := by
      calc
        (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H) *
            (Wedge2Formalization.N4.J (k := k) * Hᵀ)
            =
          Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
            (H * Wedge2Formalization.N4.J (k := k) * Hᵀ) := by
              simp [Matrix.mul_assoc]
        _ =
          Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
            Wedge2Formalization.N4.J (k := k) := by
              rw [hHfix]
        _ = -(Wedge2Formalization.N4.J (k := k) * Xᵀ) := by
              rw [coupledUpperFrom_mul_J (k := k) (X := X)]
    letI := Matrix.invertibleOfIsUnitDet
      (Wedge2Formalization.N4.J (k := k) * Hᵀ) hunitJH
    have hUshape :
        U0 =
          Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H := by
      have hsame :
          U0 * (Wedge2Formalization.N4.J (k := k) * Hᵀ) =
            (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H) *
              (Wedge2Formalization.N4.J (k := k) * Hᵀ) := by
        rw [hUJ, hUtarget]
      exact
        (Matrix.mul_left_inj_of_invertible
          (A := Wedge2Formalization.N4.J (k := k) * Hᵀ)).mp hsame
    have hUJU :
        U0 * Wedge2Formalization.N4.J (k := k) * U0ᵀ =
          Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
            Wedge2Formalization.N4.J (k := k) *
            (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)ᵀ := by
      calc
        U0 * Wedge2Formalization.N4.J (k := k) * U0ᵀ
            =
          (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H) *
            Wedge2Formalization.N4.J (k := k) *
            (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H)ᵀ := by
              rw [hUshape]
        _ =
          Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
            (H * Wedge2Formalization.N4.J (k := k) * Hᵀ) *
            (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)ᵀ := by
              simp [Matrix.mul_assoc, Matrix.transpose_mul]
        _ =
          Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
            Wedge2Formalization.N4.J (k := k) *
            (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)ᵀ := by
              simpa [Matrix.mul_assoc] using congrArg
                (fun M =>
                  Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X *
                    M *
                    (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X)ᵀ)
                hHfix
    let B0 : Matrix Wedge2Formalization.N6MixedOnePoint.I
        Wedge2Formalization.N6MixedOnePoint.I k :=
      B - Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A
    have hB0 :
        A * Wedge2Formalization.N4.J (k := k) * B0ᵀ +
          B0 * Wedge2Formalization.N4.J (k := k) * Aᵀ = 0 := by
      have hcenter :
          A * Wedge2Formalization.N4.J (k := k) *
              (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A)ᵀ +
            (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A) *
              Wedge2Formalization.N4.J (k := k) * Aᵀ =
              Wedge2Formalization.N4.J (k := k) *
                (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X)ᵀ +
              Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X *
                Wedge2Formalization.N4.J (k := k) := by
        calc
          A * Wedge2Formalization.N4.J (k := k) *
              (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A)ᵀ +
            (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A) *
              Wedge2Formalization.N4.J (k := k) * Aᵀ
              =
            (A * Wedge2Formalization.N4.J (k := k) * Aᵀ) *
                (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X)ᵀ +
              Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X *
                (A * Wedge2Formalization.N4.J (k := k) * Aᵀ) := by
                  simp [Matrix.transpose_mul, Matrix.mul_assoc]
          _ =
            Wedge2Formalization.N4.J (k := k) *
                (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X)ᵀ +
              Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X *
                Wedge2Formalization.N4.J (k := k) := by
                  rw [hAfix]
      calc
        A * Wedge2Formalization.N4.J (k := k) * B0ᵀ +
            B0 * Wedge2Formalization.N4.J (k := k) * Aᵀ
            =
          (A * Wedge2Formalization.N4.J (k := k) * Bᵀ +
              B * Wedge2Formalization.N4.J (k := k) * Aᵀ) -
            (A * Wedge2Formalization.N4.J (k := k) *
                (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A)ᵀ +
              (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A) *
                Wedge2Formalization.N4.J (k := k) * Aᵀ) := by
                  dsimp [B0]
                  simp [Matrix.transpose_sub, Matrix.transpose_mul, Matrix.mul_assoc,
                    Matrix.mul_add, Matrix.add_mul, Matrix.mul_sub, Matrix.sub_mul,
                    sub_eq_add_neg, add_assoc, add_left_comm, add_comm]
        _ =
          (A * Wedge2Formalization.N4.J (k := k) * Bᵀ +
              B * Wedge2Formalization.N4.J (k := k) * Aᵀ) -
            (Wedge2Formalization.N4.J (k := k) *
                (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X)ᵀ +
              Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X *
                Wedge2Formalization.N4.J (k := k)) := by
                  rw [hcenter]
        _ =
          -(U0 * Wedge2Formalization.N4.J (k := k) * U0ᵀ) -
            (Wedge2Formalization.N4.J (k := k) *
                (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X)ᵀ +
              Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X *
                Wedge2Formalization.N4.J (k := k)) := by
                  rw [eq_neg_iff_add_eq_zero.mpr hTop]
        _ = 0 := by
              rw [coupledCenterFrom_J_relation (k := k) (X := X), hUJU]
              simp
    have hBlockShape :
        g =
          Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 0 A 0 0 Y H := by
      calc
        g =
          Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 C D V0 P Y H := hgBlock
        _ =
          Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 0 A 0 0 Y H := by
            rw [hC0, hD, hV0, hP0]
    have hKernelProd :
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X *
            Levi (k := k) A B0 H =
          Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 0 A 0 0 Y H := by
      calc
        Wedge2Formalization.N6MixedOnePoint.coupledKernelFrom (k := k) X *
            Levi (k := k) A B0 H
            =
          Wedge2Formalization.N6MixedOnePoint.Block3
            A
            (Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A + B0)
            (Wedge2Formalization.N6MixedOnePoint.coupledUpperFrom (k := k) X * H)
            0
            A
            0
            0
            (X * A)
            H := coupledKernelFrom_mul_Levi (k := k) (X := X) (A := A) (B := B0) (H := H)
        _ =
          Wedge2Formalization.N6MixedOnePoint.Block3 A B U0 0 A 0 0 Y H := by
            have hBraw :
                Wedge2Formalization.N6MixedOnePoint.coupledCenterFrom (k := k) X * A + B0 = B := by
              dsimp [B0]
              simp [sub_eq_add_neg, add_assoc, add_left_comm, add_comm]
            rw [hBraw, hUshape, hXA]
    refine (mem_K_iff (k := k) (g := g)).2 ?_
    refine ⟨X, A, B0, H, hA, hB0, hH, ?_⟩
    exact hBlockShape.trans hKernelProd.symm
  · intro hg
    exact pointwise_of_mem_K (k := k) (g := g) hg

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6MixedOnePoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.mixedOnePoint_borel_lift_action
        (k := k) (a := a) (b := b) ha).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.N4PaperSummary.row1_borelCoeff] using
      (Wedge2Formalization.N6Summary.mixedOnePoint_borel_lift_action
        (k := k) (a := a) (b := b) ha).2

private theorem rep_pair_independent
    {a b : k}
    (h :
      a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h45 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h
  have ha : a = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N6MixedOnePoint.rep₂, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h45
  have h01 := congrArg (fun M => M (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inr 0))) h
  have hb : b = 0 := by
    simpa [ha, rep₁, rep₂, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N6MixedOnePoint.rep₂, Wedge2Formalization.N4.onePointRep₁,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using
      (congrArg (fun M => M (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inl 1))) h)
  exact ⟨ha, hb⟩

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  exact one_ne_zero (rep_pair_independent (k := k) hcomb).1

private theorem rep₁_det_ne_zero : Matrix.det (rep₁ (k := k)) ≠ 0 := by
  intro hdet
  have h' : (1 : k) = 0 :=
    (det_zero_iff (k := k) (a := (1 : k)) (b := (0 : k))).1 (by simpa using hdet)
  exact one_ne_zero h'

theorem det_ne_zero_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6MixedOnePoint.V
        Wedge2Formalization.N6MixedOnePoint.V k)
    (hfix :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector g) :
    Matrix.det g ≠ 0 := by
  intro hg0
  have hdet := congrArg Matrix.det hfix.1
  rw [Wedge2Formalization.N6MixedOnePoint.ActBivector, Matrix.det_mul, Matrix.det_mul,
    Matrix.det_transpose] at hdet
  have hrep0 : Matrix.det (rep₁ (k := k)) = 0 := by
    simpa [hg0] using hdet.symm
  exact rep₁_det_ne_zero (k := k) hrep0

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6MixedOnePoint.V
        Wedge2Formalization.N6MixedOnePoint.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N6MixedOnePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6MixedOnePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h2' :
      Wedge2Formalization.N6MixedOnePoint.ActBivector
          (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    exact hM.2
  have h20 :
      M 1 0 = 0 := by
    have hrep₂_zero : Matrix.det (rep₂ (k := k)) = 0 := by
      have h :
          Matrix.det ((0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k)) = 0 :=
        (det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 rfl
      simpa using h
    have hrep₂det :
        Matrix.det
          (Wedge2Formalization.N6MixedOnePoint.ActBivector
            (rep₂ (k := k)) g) = 0 := by
      rw [Wedge2Formalization.N6MixedOnePoint.ActBivector, Matrix.det_mul, Matrix.det_mul,
        Matrix.det_transpose, hrep₂_zero]
      simp
    exact (det_zero_iff (k := k) (a := M 1 0) (b := M 1 1)).1 (h2'.symm ▸ hrep₂det)
  refine ⟨M, ?_, hM⟩
  constructor
  · simpa [h20]
  ·
    change Matrix.det M ≠ 0
    by_contra hdet
    have hdet_formula : M 0 0 * M 1 1 = 0 := by
      simpa [Matrix.det_fin_two, h20] using hdet
    rcases mul_eq_zero.mp hdet_formula with h00 | h11
    · have hsing :
          Matrix.det
            (Wedge2Formalization.N6MixedOnePoint.ActBivector
              (rep₁ (k := k)) g) = 0 := by
        rcases hM with ⟨h1, _⟩
        rw [h1]
        have hrep₂_zero : Matrix.det (rep₂ (k := k)) = 0 := by
          have h :
              Matrix.det ((0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k)) = 0 :=
            (det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 rfl
          simpa using h
        simpa [h00, h20, hrep₂_zero]
      have hnonsing :
          Matrix.det
            (Wedge2Formalization.N6MixedOnePoint.ActBivector
              (rep₁ (k := k)) g) ≠ 0 := by
        rw [Wedge2Formalization.N6MixedOnePoint.ActBivector, Matrix.det_mul, Matrix.det_mul,
          Matrix.det_transpose]
        simp [rep₁_det_ne_zero (k := k), hg]
      exact hnonsing hsing
    · have hzero :
          Wedge2Formalization.N6MixedOnePoint.ActBivector
              (rep₂ (k := k)) g = 0 := by
        rcases hM with ⟨_, h2⟩
        simpa [ActsOnOrderedPair, h20, h11] using h2
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
      have hrep₂nonzero : rep₂ (k := k) ≠ 0 := by
        intro hzero'
        have hcomb : (0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k) = 0 := by
          simpa using hzero'
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).2
      exact hrep₂nonzero horig

end Row2

end N6
end Paper
end Wedge2Formalization
