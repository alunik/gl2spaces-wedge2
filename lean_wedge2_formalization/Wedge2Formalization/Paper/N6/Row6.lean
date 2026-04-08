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

/-! Appendix A, `n = 6`, row 6.
Representative `⟨e₃∧e₅ + e₄∧e₆, e₃∧e₆⟩`.
Divisor `2[a]`.
Claimed stabilizer:
`K_L = U_{11} \rtimes (GL_2(k) \times SL_2(k))`, exact quotient family `Q_L = B`.
-/
namespace Row6

private theorem act_embedded_toBlocks₂₂_N6
    (Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    Matrix.toBlocks₂₂
        (Wedge2Formalization.N6.ActBivector (Matrix.fromBlocks 0 0 0 Ω) g) =
      Wedge2Formalization.N4.ActBivector Ω g.toBlocks₂₂ := by
  rw [← Matrix.fromBlocks_toBlocks g]
  simp [Wedge2Formalization.N6.ActBivector, Wedge2Formalization.N4.ActBivector,
    Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

def rep₁ : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)

def rep₂ : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 6:
`e₃∧e₅ + e₄∧e₆`. -/
def paperRep₁ : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))

/-- Literal second basis vector from Appendix A, row 6:
`e₃∧e₆`. -/
def paperRep₂ : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))

/-- Explicit basis change sending the literal paper representative to the internal
working representative pair. -/
def paperChange : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    0
    0
    (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))

/-- Transport of the literal paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6.ActBivector (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  calc
    Wedge2Formalization.N6.ActBivector (paperRep₁ (k := k)) (paperChange (k := k)) =
        Matrix.fromBlocks 0 0 0
          (Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))
            (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))) := by
          simpa [paperRep₁, paperChange] using
            (Wedge2Formalization.N6.act_embedded_fromBlocks_zeroUpperRight
              (k := k)
              (Ω := Wedge2Formalization.Paper.N4.Row1.paperRep₁ (k := k))
              (A := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
              (C := (0 : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k))
              (D := Wedge2Formalization.Paper.N4.Row1.paperChange (k := k)))
    _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.N4.onePointRep₁ (k := k)) := by
          simpa [Wedge2Formalization.Paper.N4.Row1.rep₁] using
            congrArg
              (fun Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k =>
                Matrix.fromBlocks
                  (0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
                  (0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N4.V k)
                  (0 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N6.I k)
                  Ω)
              (Wedge2Formalization.Paper.N4.Row1.paperRep₁_transport (k := k))
    _ = rep₁ (k := k) := by
          simp [rep₁, Wedge2Formalization.N6OnePoint.rad2OnePointRep₁]

/-- Transport of the literal paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6.ActBivector (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  calc
    Wedge2Formalization.N6.ActBivector (paperRep₂ (k := k)) (paperChange (k := k)) =
        Matrix.fromBlocks 0 0 0
          (Wedge2Formalization.N4.ActBivector
            (Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))
            (Wedge2Formalization.Paper.N4.Row1.paperChange (k := k))) := by
          simpa [paperRep₂, paperChange] using
            (Wedge2Formalization.N6.act_embedded_fromBlocks_zeroUpperRight
              (k := k)
              (Ω := Wedge2Formalization.Paper.N4.Row1.paperRep₂ (k := k))
              (A := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
              (C := (0 : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k))
              (D := Wedge2Formalization.Paper.N4.Row1.paperChange (k := k)))
    _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
          simpa [Wedge2Formalization.Paper.N4.Row1.rep₂] using
            congrArg
              (fun Ω : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k =>
                Matrix.fromBlocks
                  (0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
                  (0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N4.V k)
                  (0 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N6.I k)
                  Ω)
              (Wedge2Formalization.Paper.N4.Row1.paperRep₂_transport (k := k))
    _ = rep₂ (k := k) := by
          simp [rep₂, Wedge2Formalization.N6OnePoint.rad2OnePointRep₂]

/-- Exact pointwise kernel family for the radical repeated-support row. -/
def K : Set (Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :=
  { g |
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ A B : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
            g = Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A B 0 A) }

/-- The arbitrary top `2×2` block on the radical quotient. -/
def TopBlock
    (A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks
    A0
    0
    0
    (1 : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)

/-- The radical unipotent family together with the embedded repeated-support unipotent
factor. -/
def U
    (C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k)
    (x y z : k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    0
    C
    (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z)

/-- The embedded repeated-support `SL₂(k)` Levi factor. -/
def CoreLevi
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    0
    0
    (Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)

/-- Exact coefficient-side projective Borel family. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

/-- Chosen Borel lift on the repeated-support quotient. -/
def lift
    (a b : k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  let B : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k :=
    !![b * (a⁻¹ * a⁻¹), 0; 0, 0]
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    0
    0
    (Wedge2Formalization.N4.onePointScale (k := k) a *
      Wedge2Formalization.N4.onePointUpperShear (k := k) B)

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6OnePoint.FixesRad2OnePointPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N6Summary.rad2OnePoint_pointwise_bivector_iff
        (k := k)
        (A0 := g.toBlocks₁₁)
        (R := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (A := g.toBlocks₂₂.toBlocks₁₁)
        (B₁ := g.toBlocks₂₂.toBlocks₁₂)
        (C₁ := g.toBlocks₂₂.toBlocks₂₁)
        (D := g.toBlocks₂₂.toBlocks₂₂)).1
        (by
          have hEq22 :
              Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
                g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂ = g.toBlocks₂₂ := by
            simpa using Matrix.fromBlocks_toBlocks g.toBlocks₂₂
          have hEq :
              Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁
                (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
                  g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂) = g := by
            rw [hEq22]
            exact Matrix.fromBlocks_toBlocks g
          simpa [hEq] using hg) with
      ⟨hR, hA, hC, hD, hrel⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, g.toBlocks₂₂.toBlocks₁₁,
      g.toBlocks₂₂.toBlocks₁₂, hA, ?_, ?_⟩
    · simpa [hD]
    · calc
        g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
              exact (Matrix.fromBlocks_toBlocks g).symm
        _ =
            Matrix.fromBlocks
              g.toBlocks₁₁
              0
              g.toBlocks₂₁
              (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
                g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂) := by
              rw [hR]
              exact congrArg
                (fun X => Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ X)
                (Matrix.fromBlocks_toBlocks g.toBlocks₂₂).symm
        _ =
            Matrix.fromBlocks
              g.toBlocks₁₁
              0
              g.toBlocks₂₁
              (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
                0 g.toBlocks₂₂.toBlocks₁₁) := by
              rw [hC, hD]
  · rintro ⟨A0, C, A, B, hA, hrel, rfl⟩
    exact
      (Wedge2Formalization.N6Summary.rad2OnePoint_pointwise_bivector_iff
        (k := k) (A0 := A0) (R := 0) (C := C) (A := A) (B₁ := B) (C₁ := 0) (D := A)).2
        ⟨rfl, hA, rfl, rfl, hrel⟩

theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ A B : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            A * Wedge2Formalization.N4.J * Bᵀ + B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
            g = Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A B 0 A) := by
  rfl

/-- The displayed top block, radical unipotent family, and embedded repeated-support
Levi multiply to the raw pointwise block shape. -/
private theorem topBlock_mul_U
    (A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    (C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k)
    (x y z : k) :
    TopBlock (k := k) A0 * U (k := k) C x y z =
      Matrix.fromBlocks
        A0
        0
        C
        (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z) := by
  rw [TopBlock, U, Matrix.fromBlocks_multiply]
  simp

private theorem embedded_mul_coreLevi
    (A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    (C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k)
    (H : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix.fromBlocks A0 0 C H * CoreLevi (k := k) A =
      Matrix.fromBlocks
        A0
        0
        C
        (H * Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) := by
  rw [CoreLevi, Matrix.fromBlocks_multiply]
  simp

private theorem topBlock_mul_U_mul_coreLevi
    (A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    (C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k)
    (x y z : k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    TopBlock (k := k) A0 * U (k := k) C x y z * CoreLevi (k := k) A =
      Matrix.fromBlocks
        A0
        0
        C
        (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
          Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) := by
  calc
    TopBlock (k := k) A0 * U (k := k) C x y z * CoreLevi (k := k) A =
        Matrix.fromBlocks
          A0
          0
          C
          (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z) *
        CoreLevi (k := k) A := by
          rw [topBlock_mul_U]
    _ =
        Matrix.fromBlocks
          A0
          0
          C
          ((Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z) *
            Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) := by
          rw [embedded_mul_coreLevi]

/-- The embedded repeated-support core, i.e. the `U_3 \rtimes SL_2(k)` factor
from Appendix A, `n = 4`, row 1. -/
def OnePointCore :
    Set
      (Matrix Wedge2Formalization.N4.V
        Wedge2Formalization.N4.V k) :=
  Wedge2Formalization.Paper.N4.Row1.K (k := k)

/-- Secondary kernel-membership statement through the embedded `n = 4`, row 1 core. -/
theorem mem_K_core_iff
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ H : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k,
            H ∈ OnePointCore (k := k) ∧
            g = Matrix.fromBlocks A0 0 C H := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨A0, C, A, B, hA, hrel, rfl⟩
    refine ⟨A0, C, Matrix.fromBlocks A B 0 A, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks A B 0 A)).2
        ⟨A, B, hA, hrel, rfl⟩
  · rintro ⟨A0, C, H, hH, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := H)).1 hH with
      ⟨A, B, hA, hrel, rfl⟩
    exact
      (mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A B 0 A))).2
        ⟨A0, C, A, B, hA, hrel, rfl⟩

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ x y z : k,
            ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
              A.det = 1 ∧
              g = TopBlock (k := k) A0 * U (k := k) C x y z * CoreLevi (k := k) A := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨A0, C, A, B, hA, hrel, hEq⟩
    have hcore :
        Matrix.fromBlocks A B 0 A ∈ Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := Matrix.fromBlocks A B 0 A)).2
        ⟨A, B, hA, hrel, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row1.mem_K_iff
        (k := k)
        (g := Matrix.fromBlocks A B 0 A)).1 hcore with
      ⟨A', x, y, z, hA', hcoreEq⟩
    refine ⟨A0, C, x, y, z, A', hA', ?_⟩
    calc
      g = Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A B 0 A) := hEq
      _ =
        Matrix.fromBlocks
          A0
          0
          C
          (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
            Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A') := by
              rw [hcoreEq]
      _ = TopBlock (k := k) A0 * U (k := k) C x y z * CoreLevi (k := k) A' := by
              rw [topBlock_mul_U_mul_coreLevi]
  · rintro ⟨A0, C, x, y, z, A, hA, hEq⟩
    have hcore :
        Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
            Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A ∈
          Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
      (Wedge2Formalization.Paper.N4.Row1.mem_K_iff
        (k := k)
        (g := Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
          Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)).2
        ⟨A, x, y, z, hA, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row1.mem_K_shape_iff
        (k := k)
        (g := Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
          Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)).1 hcore with
      ⟨A', B, hA', hrel, hcoreEq⟩
    apply (mem_K_shape_iff (k := k) (g := g)).2
    refine ⟨A0, C, A', B, hA', hrel, ?_⟩
    rw [hEq, topBlock_mul_U_mul_coreLevi, hcoreEq]

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_{11} \rtimes (GL_2(k) \times SL_2(k))`, written in the displayed
`TopBlock * U * CoreLevi` factorization. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ x y z : k,
            ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
              A.det = 1 ∧
              g = TopBlock (k := k) A0 * U (k := k) C x y z * CoreLevi (k := k) A :=
  mem_K_iff (k := k) (g := g)

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.rad2OnePoint_borel_lift_action
        (k := k)
        (A0 := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
        (a := a) (b := b) ha (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.N4PaperSummary.row1_borelCoeff] using
      (Wedge2Formalization.N6Summary.rad2OnePoint_borel_lift_action
        (k := k)
        (A0 := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
        (a := a) (b := b) ha (C := 0)).2

/-- Rank drop on the repeated-support quotient pencil occurs exactly at the repeated point. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det
      (a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
        b • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 ↔
      a = 0 :=
  Wedge2Formalization.N4Summary.onePoint_det_zero_iff (k := k) (a := a) (b := b)

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          b • (Wedge2Formalization.N4.onePointRep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6OnePoint.rad2OnePointRep₁,
      Wedge2Formalization.N6OnePoint.rad2OnePointRep₂] using h'
  have ha :
      a = 0 := by
    have hdet :
        Matrix.det
          (a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
            b • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 := by
      rw [h'']
      simpa using (Matrix.det_zero (n := Wedge2Formalization.N4.V) (R := k))
    exact
      (Wedge2Formalization.N4Summary.onePoint_det_zero_iff (k := k) (a := a) (b := b)).1 hdet
  have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h''
  have hb : b = 0 := by
    simpa [ha, Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.J,
      Matrix.add_apply, Matrix.smul_apply] using h01
  exact ⟨ha, hb⟩

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  rcases rep_pair_independent (k := k) hcomb with ⟨_, h1⟩
  exact one_ne_zero h1

private theorem lowerRight_action
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k)
    {α β γ δ : k}
    (h1 :
      Wedge2Formalization.N6.ActBivector
        (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) g
        =
      α • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
        β • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k)))
    (h2 :
      Wedge2Formalization.N6.ActBivector
        (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k)) g
        =
      γ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
        δ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k))) :
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.onePointRep₁ (k := k)) g.toBlocks₂₂ =
        α • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.onePointRep₂ (k := k)) ∧
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.onePointRep₂ (k := k)) g.toBlocks₂₂ =
        γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
  have h1' := congrArg Matrix.toBlocks₂₂ h1
  have h2' := congrArg Matrix.toBlocks₂₂ h2
  have hrhs1 :
      Matrix.toBlocks₂₂
        (α • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
          β • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k))) =
        α • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
    ext i j
    change
      (α • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
        β • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (α • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
        β • (Wedge2Formalization.N4.onePointRep₂ (k := k))) i j
    simp [Wedge2Formalization.N6OnePoint.rad2OnePointRep₁, Wedge2Formalization.N6OnePoint.rad2OnePointRep₂,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
  have hrhs2 :
      Matrix.toBlocks₂₂
        (γ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k))) =
        γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
    ext i j
    change
      (γ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
        δ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
        δ • (Wedge2Formalization.N4.onePointRep₂ (k := k))) i j
    simp [Wedge2Formalization.N6OnePoint.rad2OnePointRep₁, Wedge2Formalization.N6OnePoint.rad2OnePointRep₂,
      Wedge2Formalization.N4.onePointRep₁, Wedge2Formalization.N4.onePointRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
  constructor
  · calc
      Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.onePointRep₁ (k := k)) g.toBlocks₂₂
          =
        Matrix.toBlocks₂₂
          (Wedge2Formalization.N6.ActBivector
            (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) g) := by
              symm
              simpa [Wedge2Formalization.N6OnePoint.rad2OnePointRep₁] using
                act_embedded_toBlocks₂₂_N6 (k := k)
                  (Ω := Wedge2Formalization.N4.onePointRep₁ (k := k))
                  (g := g)
      _ = Matrix.toBlocks₂₂
            (α • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
              β • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k))) := by
                simpa using h1'
      _ =
          α • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
            β • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := hrhs1
  · calc
      Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.onePointRep₂ (k := k)) g.toBlocks₂₂
          =
        Matrix.toBlocks₂₂
          (Wedge2Formalization.N6.ActBivector
            (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k)) g) := by
              symm
              simpa [Wedge2Formalization.N6OnePoint.rad2OnePointRep₂] using
                act_embedded_toBlocks₂₂_N6 (k := k)
                  (Ω := Wedge2Formalization.N4.onePointRep₂ (k := k))
                  (g := g)
      _ = Matrix.toBlocks₂₂
            (γ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₁ (k := k)) +
              δ • (Wedge2Formalization.N6OnePoint.rad2OnePointRep₂ (k := k))) := by
                simpa using h2'
      _ =
          γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
            δ • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := hrhs2

theorem quotient_image
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N6OnePoint.PreservesRad2OnePointSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N6.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N6.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  rcases lowerRight_action (k := k) (g := g) h1 h2 with ⟨hD1, hD2⟩
  have hrep₂_det :
      Matrix.det (Wedge2Formalization.N4.onePointRep₂ (k := k)) = 0 := by
    rw [Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Matrix.det_fromBlocks_zero₂₁]
    simpa using (Matrix.det_zero (n := Wedge2Formalization.N4.I) (R := k))
  have hrep₂_det_zero :
      Matrix.det
        (γ • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 := by
    rw [← hD2, Wedge2Formalization.N4Summary.actBivector_det, hrep₂_det]
    simp
  have hγ : γ = 0 :=
    Wedge2Formalization.N4Summary.onePoint_action_case_of_rankDrop
      (k := k) (γ := γ) (δ := δ) hrep₂_det_zero
  have hδ : δ ≠ 0 := by
    intro hδ0
    have hzero :
        Wedge2Formalization.N6.ActBivector (rep₂ (k := k)) g = 0 := by
      simpa [hγ, hδ0] using h2'
    have horig :=
      (actBivector_eq_zero_iff_of_det_ne_zero
        (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
    exact rep₂_ne_zero (k := k) horig
  have hα : α ≠ 0 := by
    intro hα0
    have hzero :
        Wedge2Formalization.N6.ActBivector
          (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g = 0 := by
      calc
        Wedge2Formalization.N6.ActBivector
            (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g
            =
          δ • Wedge2Formalization.N6.ActBivector (rep₁ (k := k)) g +
            (-β) • Wedge2Formalization.N6.ActBivector (rep₂ (k := k)) g := by
              simp [Wedge2Formalization.N6.ActBivector,
                Matrix.mul_add, Matrix.add_mul]
        _ = 0 := by
              rw [h1', h2']
              ext i j
              simp [rep₁, rep₂, hα0, hγ]
              ring
    have horig :=
      (actBivector_eq_zero_iff_of_det_ne_zero
        (Ω := δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) (g := g) hg).1 hzero
    have hcomb_ne :
        δ • rep₁ (k := k) + (-β) • rep₂ (k := k) ≠ 0 := by
      intro hcomb
      rcases rep_pair_independent (k := k) hcomb with ⟨hδ', hβ'⟩
      exact hδ hδ'
    exact hcomb_ne horig
  refine ⟨!![α, β; 0, δ], ?_, ?_⟩
  · constructor
    · simp
    · simp [Matrix.det_fin_two, hα, hδ]
  · constructor
    · simpa [ActsOnOrderedPair] using h1
    · simpa [ActsOnOrderedPair, hγ] using h2

end Row6

end N6
end Paper
end Wedge2Formalization
