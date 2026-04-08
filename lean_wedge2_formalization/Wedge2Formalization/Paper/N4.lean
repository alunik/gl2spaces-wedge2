import Wedge2Formalization.Paper.Core
import Wedge2Formalization.N4PaperTheorems

namespace Wedge2Formalization
namespace Paper
namespace N4

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 4`, row 1.
Paper representative `⟨e₁∧e₃ + e₂∧e₄, e₁∧e₄⟩`.
Divisor `2[a]`.
Claimed stabilizer:
`K_L = U_3 ⋊ SL₂(k)`, `Q_L = B`.
The internal pair `rep₁`, `rep₂` is the working model used by the established
`N4` one-point development. The literal paper representative is exposed as
`paperRep₁`, `paperRep₂`, together with the explicit basis change
`paperChange` transporting it to the internal pair.
-/
namespace Row1

/-- Internal first basis vector used by the existing one-point development. -/
def rep₁ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4.onePointRep₁ (k := k)

/-- Internal second basis vector used by the existing one-point development. -/
def rep₂ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4.onePointRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 1:
`e₁∧e₃ + e₂∧e₄`. -/
def paperRep₁ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Matrix.fromBlocks
    0
    (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (-(1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))
    0

/-- Literal second basis vector from Appendix A, row 1:
`e₁∧e₄`. -/
def paperRep₂ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  let E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k := !![(0 : k), 1; 0, 0]
  Matrix.fromBlocks 0 E (-Eᵀ) 0

/-- Explicit basis change sending the paper representative to the internal
working representative pair. -/
def paperChange : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Matrix.fromBlocks
    (!![(1 : k), 0; 0, 0] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (!![(0 : k), 0; 0, 1] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (!![(0 : k), 1; 0, 0] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (!![(0 : k), 0; 1, 0] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)

/-- Transport of the literal paper `\omega_1` to the internal working
representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N4.ActBivector (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N4.onePointRep₁,
      Wedge2Formalization.N4.ActBivector, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose,
      Matrix.mul_apply, Fin.sum_univ_two]

/-- Transport of the literal paper `\omega_2` to the internal working
representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N4.ActBivector (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N4.onePointRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ActBivector, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.fromBlocks_multiply, Matrix.fromBlocks_transpose,
      Matrix.mul_apply, Fin.sum_univ_two]

/-- Exact pointwise kernel `K_L`. -/
def K : Set (Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :=
  Wedge2Formalization.N4PaperTheorems.Row1.K (k := k)

/-- Displayed unipotent family `U_3`. -/
def U (x y z : k) : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4PaperTheorems.Row1.U (k := k) x y z

/-- Displayed Levi family `SL₂(k)`. -/
def Levi
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4PaperTheorems.Row1.L (k := k) A

/-- Exact coefficient-side quotient family on the ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.N4PaperTheorems.Row1.Qproj (k := k)

/-- Chosen explicit Borel section used for quotient lifts. -/
def lift (a b : k) : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4PaperTheorems.Row1.lift (k := k) a b

/-- Exact pointwise stabilizer statement for row 1. -/
theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N4.FixesOnePointPairBivector (k := k))
      (K (k := k)) :=
  Wedge2Formalization.N4PaperTheorems.Row1.pointwise_stabilizer (k := k)

/-- Raw matrix-shape description of the pointwise kernel. -/
theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :
    g ∈ K (k := k) ↔
      ∃ A B : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        A.det = 1 ∧
        A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
        g = Matrix.fromBlocks A B 0 A := by
  rfl

private theorem shape_eq_U_mul_Levi
    (A B : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1)
    (hrel : A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0) :
    let x := B 0 0 * A 1 1 - B 0 1 * A 1 0
    let y := A 0 0 * B 0 1 - A 0 1 * B 0 0
    let z := B 1 0 * A 1 1 - B 1 1 * A 1 0
    Matrix.fromBlocks A B 0 A = U (k := k) x y z * Levi (k := k) A := by
  have htrace :
      A 0 0 * B 1 1 - A 0 1 * B 1 0 + (B 0 0 * A 1 1 - B 0 1 * A 1 0) = 0 := by
    have h01 := congrArg (fun M => M 0 1) hrel
    simpa [Wedge2Formalization.N4.J, Matrix.mul_apply, Fin.sum_univ_two, add_comm, add_left_comm,
      add_assoc, sub_eq_add_neg] using h01
  have hx :
      B 0 0 * A 1 1 - B 0 1 * A 1 0 = -(A 0 0 * B 1 1 - A 0 1 * B 1 0) := by
    apply (eq_neg_iff_add_eq_zero).2
    simpa [add_comm, add_left_comm, add_assoc] using htrace
  have hdet : A 0 0 * A 1 1 - A 0 1 * A 1 0 = 1 := by
    simpa [Matrix.det_fin_two] using hA
  have hnegx : -(B 0 0 * A 1 1 - B 0 1 * A 1 0) = A 0 0 * B 1 1 - A 0 1 * B 1 0 := by
    simpa [hx]
  let M : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k :=
    !![B 0 0 * A 1 1 - B 0 1 * A 1 0, A 0 0 * B 0 1 - A 0 1 * B 0 0;
      B 1 0 * A 1 1 - B 1 1 * A 1 0, -(B 0 0 * A 1 1 - B 0 1 * A 1 0)]
  have hMA : M * A = B := by
    ext i j
    fin_cases i <;> fin_cases j
    · simp [Matrix.mul_apply, Fin.sum_univ_two, M]
      calc
        (B 0 0 * A 1 1 - B 0 1 * A 1 0) * A 0 0 +
            (A 0 0 * B 0 1 - A 0 1 * B 0 0) * A 1 0 =
          B 0 0 * (A 0 0 * A 1 1 - A 0 1 * A 1 0) := by ring
        _ = B 0 0 := by rw [hdet]; ring
    · simp [Matrix.mul_apply, Fin.sum_univ_two, M]
      calc
        (B 0 0 * A 1 1 - B 0 1 * A 1 0) * A 0 1 +
            (A 0 0 * B 0 1 - A 0 1 * B 0 0) * A 1 1 =
          B 0 1 * (A 0 0 * A 1 1 - A 0 1 * A 1 0) := by ring
        _ = B 0 1 := by rw [hdet]; ring
    · simp [Matrix.mul_apply, Fin.sum_univ_two, M]
      calc
        (B 1 0 * A 1 1 - B 1 1 * A 1 0) * A 0 0 +
            (B 0 1 * A 1 0 - B 0 0 * A 1 1) * A 1 0 =
          (B 1 0 * A 1 1 - B 1 1 * A 1 0) * A 0 0 +
            (-(B 0 0 * A 1 1 - B 0 1 * A 1 0)) * A 1 0 := by ring
        _ =
          (B 1 0 * A 1 1 - B 1 1 * A 1 0) * A 0 0 +
            (A 0 0 * B 1 1 - A 0 1 * B 1 0) * A 1 0 := by rw [hnegx]
        _ = B 1 0 * (A 0 0 * A 1 1 - A 0 1 * A 1 0) := by ring
        _ = B 1 0 := by rw [hdet]; ring
    · simp [Matrix.mul_apply, Fin.sum_univ_two, M]
      calc
        (B 1 0 * A 1 1 - B 1 1 * A 1 0) * A 0 1 +
            (B 0 1 * A 1 0 - B 0 0 * A 1 1) * A 1 1 =
          (B 1 0 * A 1 1 - B 1 1 * A 1 0) * A 0 1 +
            (-(B 0 0 * A 1 1 - B 0 1 * A 1 0)) * A 1 1 := by ring
        _ =
          (B 1 0 * A 1 1 - B 1 1 * A 1 0) * A 0 1 +
            (A 0 0 * B 1 1 - A 0 1 * B 1 0) * A 1 1 := by rw [hnegx]
        _ = B 1 1 * (A 0 0 * A 1 1 - A 0 1 * A 1 0) := by ring
        _ = B 1 1 := by rw [hdet]; ring
  dsimp [U, Levi, Wedge2Formalization.N4PaperTheorems.Row1.U, Wedge2Formalization.N4PaperTheorems.Row1.L,
    Wedge2Formalization.N4PaperSummary.row1_U3, Wedge2Formalization.N4PaperSummary.row1_levi,
    Wedge2Formalization.N4.onePointUpperShear, M]
  rw [Matrix.fromBlocks_multiply, hMA]
  simp

/-- Rank drop on the repeated-support pencil occurs exactly at the repeated point. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • rep₁ (k := k) + b • rep₂ (k := k)) = 0 ↔
      a = 0 :=
  Wedge2Formalization.N4Summary.onePoint_det_zero_iff (k := k) (a := a) (b := b)

/-- The displayed unipotent family lies in the pointwise stabilizer. -/
theorem U_pointwise
    (x y z : k) :
    Wedge2Formalization.N4.FixesOnePointPairBivector (U (k := k) x y z) :=
  Wedge2Formalization.N4PaperTheorems.Row1.U_pointwise (k := k) x y z

/-- The displayed Levi family lies in the pointwise stabilizer. -/
theorem Levi_pointwise
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1) :
    Wedge2Formalization.N4.FixesOnePointPairBivector (Levi (k := k) A) :=
  Wedge2Formalization.N4PaperTheorems.Row1.L_pointwise (k := k) A hA

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        ∃ x y z : k,
          A.det = 1 ∧
          g = U (k := k) x y z * Levi (k := k) A := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) g).1 hg with ⟨A, B, hA, hrel, rfl⟩
    refine ⟨A, B 0 0 * A 1 1 - B 0 1 * A 1 0, A 0 0 * B 0 1 - A 0 1 * B 0 0,
      B 1 0 * A 1 1 - B 1 1 * A 1 0, hA, ?_⟩
    simpa using shape_eq_U_mul_Levi (k := k) A B hA hrel
  · rintro ⟨A, x, y, z, hA, hg⟩
    have hfixU := U_pointwise (k := k) x y z
    have hfixL := Levi_pointwise (k := k) A hA
    have h1 :
        Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.onePointRep₁ (k := k))
            (U (k := k) x y z * Levi (k := k) A) =
          Wedge2Formalization.N4.onePointRep₁ (k := k) := by
      calc
        Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.onePointRep₁ (k := k))
            (U (k := k) x y z * Levi (k := k) A)
            =
          Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.ActBivector
              (Wedge2Formalization.N4.onePointRep₁ (k := k))
              (Levi (k := k) A))
            (U (k := k) x y z) := by
              rw [Wedge2Formalization.N4Summary.actBivector_mul]
        _ =
          Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.onePointRep₁ (k := k))
            (U (k := k) x y z) := by
              simpa [Wedge2Formalization.N4.FixesBivector] using
                congrArg (fun M => Wedge2Formalization.N4.ActBivector M (U (k := k) x y z)) hfixL.1
        _ = Wedge2Formalization.N4.onePointRep₁ (k := k) := hfixU.1
    have h2 :
        Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.onePointRep₂ (k := k))
            (U (k := k) x y z * Levi (k := k) A) =
          Wedge2Formalization.N4.onePointRep₂ (k := k) := by
      calc
        Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.onePointRep₂ (k := k))
            (U (k := k) x y z * Levi (k := k) A)
            =
          Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.ActBivector
              (Wedge2Formalization.N4.onePointRep₂ (k := k))
              (Levi (k := k) A))
            (U (k := k) x y z) := by
              rw [Wedge2Formalization.N4Summary.actBivector_mul]
        _ =
          Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.N4.onePointRep₂ (k := k))
            (U (k := k) x y z) := by
              simpa [Wedge2Formalization.N4.FixesBivector] using
                congrArg (fun M => Wedge2Formalization.N4.ActBivector M (U (k := k) x y z)) hfixL.2
        _ = Wedge2Formalization.N4.onePointRep₂ (k := k) := hfixU.2
    apply (pointwise_stabilizer (k := k) (g := g)).1
    rw [hg]
    exact ⟨h1, h2⟩

/-- Table-facing kernel statement for Appendix A, row 1:
`K_L = U_3 \rtimes SL_2(k)`. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        ∃ x y z : k,
          A.det = 1 ∧
          g = U (k := k) x y z * Levi (k := k) A :=
  mem_K_iff (k := k) (g := g)

/-- The chosen Borel lift acts on the ordered basis by its coefficient matrix. -/
theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N4.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) :=
  Wedge2Formalization.N4PaperTheorems.Row1.quotient_action (k := k) a b ha

/-- Every invertible setwise stabilizer induces a coefficient matrix in the projective Borel. -/
theorem quotient_image
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres : Wedge2Formalization.N4.PreservesOnePointSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N4.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M :=
  Wedge2Formalization.N4PaperTheorems.Row1.quotient_image_projective (k := k) g hg hpres

end Row1

/-! Appendix A, `n = 4`, row 2.
Representative `⟨e₁∧e₂ + e₃∧e₄, e₃∧e₄⟩`.
Divisor `[a] + [b]`.
Claimed stabilizer:
`K_L = SL₂(k) × SL₂(k)`, `Q_L = N_T`.
-/
namespace Row2

/-- First basis vector of the paper representative. -/
def rep₁ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4.splitRep₁ (k := k)

/-- Second basis vector of the paper representative. -/
def rep₂ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4.splitRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 2. -/
def paperRep₁ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 2. -/
def paperRep₂ : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  rep₂ (k := k)

/-- The row-2 paper representative already agrees with the internal working pair. -/
def paperChange : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N4.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N4.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N4.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N4.ActBivector]

/-- Exact pointwise kernel `K_L`. -/
def K : Set (Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :=
  Wedge2Formalization.N4PaperTheorems.Row2.K (k := k)

/-- Displayed Levi family `SL₂(k) × SL₂(k)`. -/
def Levi
    (A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4PaperTheorems.Row2.L (k := k) A D

/-- Displayed coefficient-side normalizer of the split torus. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.N4PaperTheorems.Row2.Q (k := k)

/-- Chosen torus lift. -/
def torusLift (u v : k) : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4PaperTheorems.Row2.torusLift (k := k) u v

/-- Chosen swap-coset lift. -/
def swapLift (u v : k) : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
  Wedge2Formalization.N4PaperTheorems.Row2.swapLift (k := k) u v

/-- Exact pointwise stabilizer statement for row 2. -/
theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N4.FixesSplitPairBivector (k := k))
      (K (k := k)) :=
  Wedge2Formalization.N4PaperTheorems.Row2.pointwise_stabilizer (k := k)

/-- Raw matrix-shape description of the pointwise kernel. -/
theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :
    g ∈ K (k := k) ↔
      ∃ A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        A.det = 1 ∧
        D.det = 1 ∧
        g = Matrix.fromBlocks A 0 0 D := by
  rfl

/-- Rank drop on the split-support pencil occurs exactly at the two support points. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • rep₁ (k := k) + b • rep₂ (k := k)) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  Wedge2Formalization.N4Summary.split_det_zero_iff (k := k) (a := a) (b := b)

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :
    g ∈ K (k := k) ↔
      ∃ A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        A.det = 1 ∧
        D.det = 1 ∧
        g = Levi (k := k) A D := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) g).1 hg with ⟨A, D, hA, hD, rfl⟩
    exact ⟨A, D, hA, hD, rfl⟩
  · rintro ⟨A, D, hA, hD, rfl⟩
    exact (mem_K_shape_iff (k := k) (Levi (k := k) A D)).2 ⟨A, D, hA, hD, rfl⟩

/-- Table-facing kernel statement for Appendix A, row 2:
`K_L = SL_2(k) \times SL_2(k)`. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k) :
    g ∈ K (k := k) ↔
      ∃ A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        A.det = 1 ∧
        D.det = 1 ∧
        g = Levi (k := k) A D :=
  mem_K_iff (k := k) (g := g)

/-- The torus lift acts on the ordered basis by the displayed coefficient matrix. -/
theorem torus_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N4.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (torusLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_torusCoeff (k := k) u v) :=
  Wedge2Formalization.N4PaperTheorems.Row2.torus_action (k := k) u v

/-- The swap-coset lift acts on the ordered basis by the displayed coefficient matrix. -/
theorem swap_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N4.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (swapLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_swapCoeff (k := k) u v) :=
  Wedge2Formalization.N4PaperTheorems.Row2.swap_action (k := k) u v

/-- Every invertible setwise stabilizer induces a coefficient matrix in the split-torus
normalizer. -/
theorem quotient_image
    (g : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres : Wedge2Formalization.N4.PreservesSplitSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N4.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M :=
  Wedge2Formalization.N4PaperTheorems.Row2.quotient_image_projective (k := k) g hg hpres

end Row2

/-! Appendix A, `n = 4`, row 3.
Paper representative `⟨e₂∧e₄, e₃∧e₄⟩`.
Divisor `0`.
Claimed stabilizer:
`K_L = U_5 \rtimes (G_m(k) \times G_m(k))`, `Q_L = GL₂(k)`.
The internal pair `rep₁`, `rep₂` is the working pure-singular model. The
literal paper representative is exposed as `paperRep₁`, `paperRep₂`,
together with the explicit transport `paperChange`.
-/
namespace Row3

/-- Internal first basis vector used by the existing pure-singular development. -/
def rep₁ : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  Wedge2Formalization.N4PureSingular.ω12 (k := k)

/-- Internal second basis vector used by the existing pure-singular development. -/
def rep₂ : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  Wedge2Formalization.N4PureSingular.ω13 (k := k)

/-- Literal first basis vector from Appendix A, row 3:
`e₂∧e₄`. -/
def paperRep₁ : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  !![(0 : k), 0, 0, 0;
     0, 0, 0, 1;
     0, 0, 0, 0;
     0, -1, 0, 0]

/-- Literal second basis vector from Appendix A, row 3:
`e₃∧e₄`. -/
def paperRep₂ : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  !![(0 : k), 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 1;
     0, 0, -1, 0]

/-- Explicit basis change sending the paper pure-singular representative to the
internal pair `ω12, ω13`. -/
def paperChange : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  !![(0 : k), 0, 0, -1;
     0, 1, 0, 0;
     0, 0, 1, 0;
     1, 0, 0, 0]

/-- Transport of the literal paper `\omega_1` to the internal working
representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N4PureSingular.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N4PureSingular.ω12,
      Wedge2Formalization.N4PureSingular.ActBivector, Matrix.mul_apply, Fin.sum_univ_four]

/-- Transport of the literal paper `\omega_2` to the internal working
representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N4PureSingular.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N4PureSingular.ω13,
      Wedge2Formalization.N4PureSingular.ActBivector, Matrix.mul_apply, Fin.sum_univ_four]

/-- Exact pointwise kernel `K_L`. -/
def K : Set (Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) :=
  Wedge2Formalization.N4PaperTheorems.Row3.K (k := k)

/-- Paper-facing `GL_4` unipotent radical `U_5`. -/
def U
    (x y z u v : k) :
    Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  !![(1 : k), x, y, z;
      0, 1, 0, u;
      0, 0, 1, v;
      0, 0, 0, 1]

/-- Paper-facing torus `G_m(k) × G_m(k)`. -/
def Levi (a t : k) :
    Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  !![a, 0, 0, 0;
      0, a⁻¹, 0, 0;
      0, 0, a⁻¹, 0;
      0, 0, 0, t]

/-- The paper-facing `GL_4` stabilizer family, i.e. the raw pointwise kernel together
with the missing nonzero torus coordinate. -/
def KGL : Set (Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) :=
  { g | g ∈ K (k := k) ∧ g 3 3 ≠ 0 }

/-- Exact coefficient-side quotient family `GL₂(k)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.N4PaperTheorems.Row3.Q (k := k)

/-- Chosen `GL₂(k)` lift. -/
def lift
    (α β γ δ : k) :
    Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k :=
  Wedge2Formalization.N4PaperTheorems.Row3.lift (k := k) α β γ δ

/-- Exact pointwise stabilizer statement for row 3. -/
theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N4PureSingular.FixesPureSingularPair (k := k))
      (K (k := k)) :=
  Wedge2Formalization.N4PaperTheorems.Row3.pointwise_stabilizer (k := k)

/-- Exact pointwise stabilizer statement on the actual `GL_4` paper kernel. -/
theorem pointwise_stabilizer_GL :
    ExactFamily
      (fun g =>
        Wedge2Formalization.N4PureSingular.FixesPureSingularPair (k := k) g ∧ g 3 3 ≠ 0)
      (KGL (k := k)) := by
  intro g
  constructor
  · rintro ⟨hfix, h33⟩
    exact ⟨(pointwise_stabilizer (k := k) (g := g)).1 hfix, h33⟩
  · rintro ⟨hg, h33⟩
    exact ⟨(pointwise_stabilizer (k := k) (g := g)).2 hg, h33⟩

/-- Raw matrix-shape description of the pointwise kernel. -/
theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) :
    g ∈ K (k := k) ↔
      g 0 0 ≠ 0 ∧
        g =
          Wedge2Formalization.N4PureSingular.pureSingularShape (k := k)
            (g 0 0) (g 0 1) (g 0 2) (g 0 3) (g 1 3) (g 2 3) (g 3 3) := by
  rfl

/-- Product form of the paper-facing `U_5` and `G_m \times G_m` families. -/
private theorem U_mul_Levi_eq_shape
    (x y z u v a t : k) :
    U (k := k) x y z u v * Levi (k := k) a t =
      Wedge2Formalization.N4PureSingular.pureSingularShape
        (k := k) a (x * a⁻¹) (y * a⁻¹) (z * t) (u * t) (v * t) t := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [U, Levi, Wedge2Formalization.N4PureSingular.pureSingularShape,
      Matrix.mul_apply, Fin.sum_univ_four] <;>
    ring

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) :
    g ∈ K (k := k) ↔
      ∃ a b c d e f t : k,
        a ≠ 0 ∧
        g =
          Wedge2Formalization.N4PureSingular.pureSingularShape
            (k := k) a b c d e f t := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) g).1 hg with ⟨ha, hshape⟩
    exact ⟨g 0 0, g 0 1, g 0 2, g 0 3, g 1 3, g 2 3, g 3 3, ha, hshape⟩
  · rintro ⟨a, b, c, d, e, f, t, ha, rfl⟩
    exact (mem_K_shape_iff (k := k)
      (Wedge2Formalization.N4PureSingular.pureSingularShape
        (k := k) a b c d e f t)).2 ⟨ha, rfl⟩

/-- The displayed unipotent radical lies in the pointwise stabilizer. -/
theorem U_pointwise
    (x y z u v : k) :
    Wedge2Formalization.N4PureSingular.FixesPureSingularPair (k := k)
      (U (k := k) x y z u v) := by
  simpa [U, Wedge2Formalization.N4PureSingular.pureSingularShape]
    using Wedge2Formalization.N4PureSingular.pureSingularShape_fixes
      (k := k) 1 x y z u v 1 one_ne_zero

/-- The displayed torus lies in the pointwise stabilizer. -/
theorem Levi_pointwise
    (a t : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N4PureSingular.FixesPureSingularPair (k := k)
      (Levi (k := k) a t) := by
  simpa [Levi, Wedge2Formalization.N4PureSingular.pureSingularShape]
    using Wedge2Formalization.N4PureSingular.pureSingularShape_fixes
      (k := k) a 0 0 0 0 0 t ha

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_5 \rtimes (G_m(k) \times G_m(k))`, now expressed on the actual
`GL_4` kernel family. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k) :
    g ∈ KGL (k := k) ↔
      ∃ a t x y z u v : k,
        a ≠ 0 ∧
        t ≠ 0 ∧
        g = U (k := k) x y z u v * Levi (k := k) a t := by
  constructor
  · rintro ⟨hg, h33⟩
    rcases (mem_K_iff (k := k) (g := g)).1 hg with ⟨a, b, c, d, e, f, t, ha, rfl⟩
    have ht : t ≠ 0 := by
      simpa [Wedge2Formalization.N4PureSingular.pureSingularShape] using h33
    refine ⟨a, t, b * a, c * a, d * t⁻¹, e * t⁻¹, f * t⁻¹, ha, ht, ?_⟩
    simpa [ha, ht, mul_assoc] using
      (U_mul_Levi_eq_shape (k := k) (b * a) (c * a) (d * t⁻¹) (e * t⁻¹) (f * t⁻¹) a t).symm
  · rintro ⟨a, t, x, y, z, u, v, ha, ht, rfl⟩
    refine ⟨?_, ?_⟩
    · exact
        (mem_K_iff (k := k)
          (g := U (k := k) x y z u v * Levi (k := k) a t)).2
          ⟨a, x * a⁻¹, y * a⁻¹, z * t, u * t, v * t, t, ha,
            by simpa using U_mul_Levi_eq_shape (k := k) x y z u v a t⟩
    · simp [U, Levi, ht]

/-- The chosen `GL₂(k)` lift acts on the ordered basis by the displayed coefficient matrix. -/
theorem quotient_action
    (α β γ δ : k)
    (hΔ : α * δ - β * γ ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N4PureSingular.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) α β γ δ)
      (Wedge2Formalization.N4PaperSummary.row3_GL2Coeff (k := k) α β γ δ) := by
  exact Wedge2Formalization.N4PaperTheorems.Row3.quotient_action (k := k) α β γ δ hΔ

/-- Every invertible setwise stabilizer induces an invertible coefficient matrix on the
ordered pure singular basis. -/
theorem quotient_image
    (g : Matrix Wedge2Formalization.N4PureSingular.I Wedge2Formalization.N4PureSingular.I k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N4PureSingular.PreservesPureSingularSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N4PureSingular.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M :=
  Wedge2Formalization.N4PaperTheorems.Row3.quotient_image_projective (k := k) g hg hpres

end Row3

end N4
end Paper
end Wedge2Formalization
