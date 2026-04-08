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

/-! Appendix A, `n = 6`, row 10.
Representative `⟨e₂∧e₅ + e₃∧e₆, e₃∧e₅ + e₄∧e₆⟩`.
Divisor `0`.
Claimed stabilizer:
`K_L = U_8 \rtimes (G_m(k) \times G_m(k))`, exact quotient family
`Q_L = GL₂(k)`.
-/
namespace Row10

def rep₁ :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  Wedge2Formalization.N6PureSingularLong.rep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  Wedge2Formalization.N6PureSingularLong.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 10:
`e₂∧e₅ + e₃∧e₆`. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k))

/-- Literal second basis vector from Appendix A, row 10:
`e₃∧e₅ + e₄∧e₆`. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k))

/-- The paper representative already agrees with the internal working representative. -/
def paperChange :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  1

/-- Transport of the literal paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6PureSingularLong.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6PureSingularLong.rep₁,
    Wedge2Formalization.Paper.N5.Row1.rep₁, Wedge2Formalization.N6PureSingularLong.ActBivector]

/-- Transport of the literal paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6PureSingularLong.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6PureSingularLong.rep₂,
    Wedge2Formalization.Paper.N5.Row1.rep₂, Wedge2Formalization.N6PureSingularLong.ActBivector]

/-- The displayed unipotent family on the radical pure singular length-two row. -/
def U
    (C : Matrix Wedge2Formalization.N6PureSingularLong.W
      Wedge2Formalization.N6PureSingularLong.I k)
    (a b c d : k) :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) (1 : k))
    0
    C
    (Matrix.fromBlocks
      (1 : Matrix Wedge2Formalization.N5PureSingularLong.W
        Wedge2Formalization.N5PureSingularLong.W k)
      0
      (Wedge2Formalization.N5PureSingularLong.lowerHankel (k := k) a b c d)
      (1 : Matrix Wedge2Formalization.N5PureSingularLong.I
        Wedge2Formalization.N5PureSingularLong.I k))

/-- The displayed Levi/core family on the radical pure singular length-two row. -/
def Levi
    (u t : k) :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) u)
    0
    0
    (Matrix.fromBlocks
      ((t⁻¹) • (1 : Matrix Wedge2Formalization.N5PureSingularLong.W
        Wedge2Formalization.N5PureSingularLong.W k))
      0
      0
      (t • (1 : Matrix Wedge2Formalization.N5PureSingularLong.I
        Wedge2Formalization.N5PureSingularLong.I k)))

/-- Exact pointwise kernel family for the radical pure singular length-two row. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6PureSingularLong.V
        Wedge2Formalization.N6PureSingularLong.V k) :=
  { g |
      ∃ A0 : Matrix Wedge2Formalization.N6PureSingularLong.I
          Wedge2Formalization.N6PureSingularLong.I k,
        ∃ C0 : Matrix Wedge2Formalization.N6PureSingularLong.W
            Wedge2Formalization.N6PureSingularLong.I k,
          ∃ a u v w z : k,
            a ≠ 0 ∧
              g =
                Matrix.fromBlocks
                  A0
                  0
                  C0
                  (Wedge2Formalization.N5PureSingularLong.pointwiseShape
                    (k := k) a u v w z) }

/-- Exact coefficient-side quotient family `GL₂(k)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

/-- The coefficient matrix attached to the determinant-scaled `GL₂(k)` lift. -/
def coeff
    (α β γ δ : k) :
    Matrix (Fin 2) (Fin 2) k :=
  !![(Wedge2Formalization.N5PureSingularLong.Delta α β γ δ) * α,
     (Wedge2Formalization.N5PureSingularLong.Delta α β γ δ) * γ;
     (Wedge2Formalization.N5PureSingularLong.Delta α β γ δ) * β,
     (Wedge2Formalization.N5PureSingularLong.Delta α β γ δ) * δ]

/-- Chosen determinant-scaled `GL₂(k)` lift. -/
def lift
    (α β γ δ : k) :
    Matrix Wedge2Formalization.N6PureSingularLong.V
      Wedge2Formalization.N6PureSingularLong.V k :=
  let P : Matrix Wedge2Formalization.N5PureSingularLong.W
      Wedge2Formalization.N5PureSingularLong.W k :=
    Wedge2Formalization.N5PureSingularLong.sym2Raw (k := k) α β γ δ
  let Q : Matrix Wedge2Formalization.N5PureSingularLong.I
      Wedge2Formalization.N5PureSingularLong.I k :=
    Wedge2Formalization.N5PureSingularLong.adjointT (k := k) α β γ δ
  let h : Matrix Wedge2Formalization.N6PureSingularLong.W
      Wedge2Formalization.N6PureSingularLong.W k := Matrix.fromBlocks P 0 0 Q
  Matrix.fromBlocks
    (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) (1 : k))
    0
    0
    h

private theorem core_product_eq_shape
    (a u v w z : k) :
    Wedge2Formalization.Paper.N5.Row1.U (k := k) u v w z *
        Wedge2Formalization.Paper.N5.Row1.Levi (k := k) a =
      Wedge2Formalization.N5PureSingularLong.pointwiseShape
        (k := k) a (a * u) (a * v) (a * w) (a * z) := by
  rw [Wedge2Formalization.Paper.N5.Row1.U, Wedge2Formalization.Paper.N5.Row1.Levi,
    Wedge2Formalization.N5PureSingularLong.pointwiseShape, Matrix.fromBlocks_multiply]
  ext i j
  cases i <;> rename_i i <;> cases j <;> rename_i j <;>
    fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N5PureSingularLong.lowerHankel,
        Matrix.fromBlocks, Matrix.smul_apply]

private def topBasis : Wedge2Formalization.N6PureSingularLong.V → k
  | Sum.inl 0 => 1
  | _ => 0

private theorem rep₁_mulVec_topBasis :
    rep₁ (k := k) *ᵥ topBasis (k := k) = 0 := by
  ext i
  cases i <;>
    simp [rep₁, Wedge2Formalization.N6PureSingularLong.rep₁, topBasis, Matrix.mulVec,
      Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type, Fin.sum_univ_one,
      Fin.sum_univ_three, Fin.sum_univ_two]

private theorem rep₂_mulVec_topBasis :
    rep₂ (k := k) *ᵥ topBasis (k := k) = 0 := by
  ext i
  cases i <;>
    simp [rep₂, Wedge2Formalization.N6PureSingularLong.rep₂, topBasis, Matrix.mulVec,
      Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type, Fin.sum_univ_one,
      Fin.sum_univ_three, Fin.sum_univ_two]

private theorem core_commonRadical_zero
    (x : Wedge2Formalization.N5PureSingularLong.V → k)
    (h1 : Wedge2Formalization.N5PureSingularLong.rep₁ (k := k) *ᵥ x = 0)
    (h2 : Wedge2Formalization.N5PureSingularLong.rep₂ (k := k) *ᵥ x = 0) :
    x = 0 := by
  funext i
  cases i with
  | inl i =>
      fin_cases i
      ·
        have hcoord := congrArg (fun v => v (Sum.inr 0)) h1
        simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
          Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
      ·
        have hcoord := congrArg (fun v => v (Sum.inr 1)) h1
        simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
          Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
      ·
        have hcoord := congrArg (fun v => v (Sum.inr 1)) h2
        simpa [Wedge2Formalization.N5PureSingularLong.rep₂,
          Wedge2Formalization.N5PureSingularLong.mulY,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
          Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
  | inr i =>
      fin_cases i
      ·
        have hcoord := congrArg (fun v => v (Sum.inl 0)) h1
        simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
          Fin.sum_univ_three, Fin.sum_univ_two] using hcoord
      ·
        have hcoord := congrArg (fun v => v (Sum.inl 1)) h1
        simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct, Fintype.sum_sum_type,
          Fin.sum_univ_three, Fin.sum_univ_two] using hcoord

private theorem upperRight_zero_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6PureSingularLong.V
        Wedge2Formalization.N6PureSingularLong.V k)
    (hg : Matrix.det g ≠ 0)
    (hfix : Wedge2Formalization.N6PureSingularLong.FixesPairBivector (k := k) g) :
    g.toBlocks₁₂ = 0 := by
  rcases hfix with ⟨h₁, h₂⟩
  have h₁' :
      Wedge2Formalization.N6PureSingularLong.ActBivector (rep₁ (k := k)) g =
        rep₁ (k := k) := by
    simpa [rep₁] using h₁
  have h₂' :
      Wedge2Formalization.N6PureSingularLong.ActBivector (rep₂ (k := k)) g =
        rep₂ (k := k) := by
    simpa [rep₂] using h₂
  have hpre₁ :
      Wedge2Formalization.N6PureSingularLong.rep₁ (k := k) *ᵥ
          (gᵀ *ᵥ topBasis (k := k)) = 0 := by
    have hzero :
        Wedge2Formalization.N6PureSingularLong.ActBivector (rep₁ (k := k)) g *ᵥ
            topBasis (k := k) = 0 := by
      rw [h₁']
      exact rep₁_mulVec_topBasis (k := k)
    have hzero' :
        g *ᵥ
            (Wedge2Formalization.N6PureSingularLong.rep₁ (k := k) *ᵥ
              (gᵀ *ᵥ topBasis (k := k))) = 0 := by
      simpa [Wedge2Formalization.N6PureSingularLong.ActBivector, Matrix.mulVec_mulVec,
        Matrix.mul_assoc] using hzero
    exact Matrix.eq_zero_of_mulVec_eq_zero hg hzero'
  have hpre₂ :
      Wedge2Formalization.N6PureSingularLong.rep₂ (k := k) *ᵥ
          (gᵀ *ᵥ topBasis (k := k)) = 0 := by
    have hzero :
        Wedge2Formalization.N6PureSingularLong.ActBivector (rep₂ (k := k)) g *ᵥ
            topBasis (k := k) = 0 := by
      rw [h₂']
      exact rep₂_mulVec_topBasis (k := k)
    have hzero' :
        g *ᵥ
            (Wedge2Formalization.N6PureSingularLong.rep₂ (k := k) *ᵥ
              (gᵀ *ᵥ topBasis (k := k))) = 0 := by
      simpa [Wedge2Formalization.N6PureSingularLong.ActBivector, Matrix.mulVec_mulVec,
        Matrix.mul_assoc] using hzero
    exact Matrix.eq_zero_of_mulVec_eq_zero hg hzero'
  let x :
      Wedge2Formalization.N5PureSingularLong.V → k :=
    fun j => (gᵀ *ᵥ topBasis (k := k)) (Sum.inr j)
  have hx₁ :
      Wedge2Formalization.N5PureSingularLong.rep₁ (k := k) *ᵥ x = 0 := by
    ext i
    have hcoord := congrArg (fun v => v (Sum.inr i)) hpre₁
    simpa [x, rep₁, Wedge2Formalization.N6PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₁, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_one, Fin.sum_univ_three,
      Fin.sum_univ_two] using hcoord
  have hx₂ :
      Wedge2Formalization.N5PureSingularLong.rep₂ (k := k) *ᵥ x = 0 := by
    ext i
    have hcoord := congrArg (fun v => v (Sum.inr i)) hpre₂
    simpa [x, rep₂, Wedge2Formalization.N6PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.rep₂, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_one, Fin.sum_univ_three,
      Fin.sum_univ_two] using hcoord
  have hx : x = 0 := core_commonRadical_zero (k := k) x hx₁ hx₂
  ext i j
  fin_cases i
  have hcoord := congrArg (fun v => v j) hx
  simpa [x, Matrix.toBlocks₁₂, Matrix.mulVec, Matrix.transpose_apply, topBasis,
    dotProduct, Fintype.sum_sum_type, Fin.sum_univ_one] using hcoord

theorem pointwise_stabilizer
    (g :
      Matrix Wedge2Formalization.N6PureSingularLong.V
        Wedge2Formalization.N6PureSingularLong.V k)
    (hg : Matrix.det g ≠ 0) :
    Wedge2Formalization.N6PureSingularLong.FixesPairBivector (k := k) g ↔
      g ∈ K (k := k) := by
  constructor
  · intro hfix
    have hB : g.toBlocks₁₂ = 0 := upperRight_zero_of_pointwise (k := k) g hg hfix
    have hEq :
        Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ = g := by
      calc
        Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂
            = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
                simpa [hB]
        _ = g := Matrix.fromBlocks_toBlocks g
    have hfix' :
        Wedge2Formalization.N6PureSingularLong.FixesPairBivector
          (Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂) := by
      simpa [hEq] using hfix
    have hDfix :
        Wedge2Formalization.N5PureSingularLong.FixesPairBivector (k := k) g.toBlocks₂₂ := by
      exact
        (Wedge2Formalization.N6PureSingularLong.fixesPair_fromBlocks_zeroUpperRight_iff
          (k := k) (A := g.toBlocks₁₁) (C := g.toBlocks₂₁) (D := g.toBlocks₂₂)).1 hfix'
    have hdetBlock :
        Matrix.det g.toBlocks₁₁ * Matrix.det g.toBlocks₂₂ ≠ 0 := by
      have hdet' :
          Matrix.det (Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂) ≠ 0 := by
        simpa [hEq] using hg
      simpa [Matrix.det_fromBlocks_zero₁₂] using hdet'
    have hDdet : Matrix.det g.toBlocks₂₂ ≠ 0 := (mul_ne_zero_iff.mp hdetBlock).2
    have hDmem :
        g.toBlocks₂₂ ∈ Wedge2Formalization.Paper.N5.Row1.K (k := k) := by
      exact
        (Wedge2Formalization.Paper.N5.Row1.pointwise_stabilizer
          (k := k) g.toBlocks₂₂ hDdet).1 hDfix
    rcases
      (Wedge2Formalization.Paper.N5.Row1.mem_K_iff
        (k := k) g.toBlocks₂₂).1 hDmem with
      ⟨a, u, v, w, z, ha, hD⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, a, a * u, a * v, a * w, a * z, ha, ?_⟩
    have hDshape :
        g.toBlocks₂₂ =
          Wedge2Formalization.N5PureSingularLong.pointwiseShape
            (k := k) a (a * u) (a * v) (a * w) (a * z) := by
      rw [hD, core_product_eq_shape]
    calc
      g =
          Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := hEq.symm
      _ =
          Matrix.fromBlocks
            g.toBlocks₁₁
            0
            g.toBlocks₂₁
            (Wedge2Formalization.N5PureSingularLong.pointwiseShape
              (k := k) a (a * u) (a * v) (a * w) (a * z)) := by
            exact congrArg
              (fun X =>
                Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ X)
              hDshape
  · rintro ⟨A0, C0, a, u, v, w, z, ha, rfl⟩
    exact
      (Wedge2Formalization.N6Summary.pureSingularLong_nested_iff_shape
        (k := k)
        (A0 := A0)
        (B0 := 0)
        (C0 := C0)
        (A1 := Wedge2Formalization.N5PureSingularLong.pointwiseShape
          (k := k) a u v w z |>.toBlocks₁₁)
        (C1 := Wedge2Formalization.N5PureSingularLong.pointwiseShape
          (k := k) a u v w z |>.toBlocks₂₁)
        (D1 := Wedge2Formalization.N5PureSingularLong.pointwiseShape
          (k := k) a u v w z |>.toBlocks₂₂)).2
        ⟨rfl, by
          simpa [Wedge2Formalization.N5PureSingularLong.pointwiseShape] using ha, by
          simp [Wedge2Formalization.N5PureSingularLong.pointwiseShape,
            Wedge2Formalization.N5PureSingularLong.lowerHankel]⟩

theorem mem_K_shape_iff
    (g :
      Matrix Wedge2Formalization.N6PureSingularLong.V
        Wedge2Formalization.N6PureSingularLong.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6PureSingularLong.I
          Wedge2Formalization.N6PureSingularLong.I k,
        ∃ C0 : Matrix Wedge2Formalization.N6PureSingularLong.W
            Wedge2Formalization.N6PureSingularLong.I k,
          ∃ a u v w z : k,
            a ≠ 0 ∧
              g =
                Matrix.fromBlocks
                  A0
                  0
                  C0
                  (Wedge2Formalization.N5PureSingularLong.pointwiseShape
                    (k := k) a u v w z) := by
  rfl

/-- The embedded pure singular length-two core on the bottom-right `5`-space. -/
def Core :
    Set
      (Matrix Wedge2Formalization.N6PureSingularLong.W
        Wedge2Formalization.N6PureSingularLong.W k) :=
  Wedge2Formalization.Paper.N5.Row1.K (k := k)

/-- The embedded pure singular length-two core, i.e. the `U_4 \rtimes G_m(k)`
factor from Appendix A, `n = 5`, row 1. -/
def PureSingularLongCore :
    Set
      (Matrix Wedge2Formalization.N6PureSingularLong.W
        Wedge2Formalization.N6PureSingularLong.W k) :=
  Core (k := k)

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6PureSingularLong.V
        Wedge2Formalization.N6PureSingularLong.V k) :
    g ∈ K (k := k) ↔
      ∃ u0 : k,
        ∃ C0 : Matrix Wedge2Formalization.N6PureSingularLong.W
            Wedge2Formalization.N6PureSingularLong.I k,
          ∃ H : Matrix Wedge2Formalization.N6PureSingularLong.W
              Wedge2Formalization.N6PureSingularLong.W k,
            H ∈ Core (k := k) ∧
              g =
                Matrix.fromBlocks
                  (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) u0)
                  0
                  C0
                  H := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨A0, C0, a, u, v, w, z, ha, rfl⟩
    have hA0 :
        Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) (A0 0 0) = A0 := by
      ext i j
      fin_cases i
      fin_cases j
      simp [Wedge2Formalization.N6PureSingularLong.scalarBlock]
    refine ⟨A0 0 0, C0,
      Wedge2Formalization.N5PureSingularLong.pointwiseShape (k := k) a u v w z, ?_, ?_⟩
    ·
      refine
        (Wedge2Formalization.Paper.N5.Row1.mem_K_iff
          (k := k)
          (g := Wedge2Formalization.N5PureSingularLong.pointwiseShape
            (k := k) a u v w z)).2 ?_
      refine ⟨a, a⁻¹ * u, a⁻¹ * v, a⁻¹ * w, a⁻¹ * z, ha, ?_⟩
      rw [core_product_eq_shape (k := k) (a := a) (u := a⁻¹ * u) (v := a⁻¹ * v)
        (w := a⁻¹ * w) (z := a⁻¹ * z)]
      have hu : a * (a⁻¹ * u) = u := by field_simp [ha]
      have hv : a * (a⁻¹ * v) = v := by field_simp [ha]
      have hw : a * (a⁻¹ * w) = w := by field_simp [ha]
      have hz : a * (a⁻¹ * z) = z := by field_simp [ha]
      simpa [hu, hv, hw, hz]
    ·
      exact congrArg
        (fun X =>
          Matrix.fromBlocks
            X
            0
            C0
            (Wedge2Formalization.N5PureSingularLong.pointwiseShape (k := k) a u v w z))
        hA0.symm
  · rintro ⟨u0, C0, H, hH, hEq⟩
    apply (mem_K_shape_iff (k := k) (g := g)).2
    rcases
      (Wedge2Formalization.Paper.N5.Row1.mem_K_iff
        (k := k)
        (g := H)).1 hH with ⟨a, u, v, w, z, ha, rfl⟩
    refine ⟨Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) u0, C0, a, a * u, a * v,
      a * w, a * z,
      ha, ?_⟩
    calc
      g =
          Matrix.fromBlocks
            (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) u0)
            0
            C0
            (Wedge2Formalization.Paper.N5.Row1.U (k := k) u v w z *
              Wedge2Formalization.Paper.N5.Row1.Levi (k := k) a) := hEq
      _ =
          Matrix.fromBlocks
            (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) u0)
            0
            C0
            (Wedge2Formalization.N5PureSingularLong.pointwiseShape
              (k := k) a (a * u) (a * v) (a * w) (a * z)) := by
                exact congrArg
                  (fun X =>
                    Matrix.fromBlocks
                      (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) u0)
                      0
                      C0
                      X)
                  (core_product_eq_shape (k := k) (a := a) (u := u) (v := v) (w := w) (z := z))

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_8 \rtimes (G_m(k) \times G_m(k))`, written in the explicit local
coordinates of this row rather than through a nested `n = 5` core theorem. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6PureSingularLong.V
        Wedge2Formalization.N6PureSingularLong.V k) :
    g ∈ K (k := k) ↔
      ∃ u0 : k,
        ∃ C0 : Matrix Wedge2Formalization.N6PureSingularLong.W
            Wedge2Formalization.N6PureSingularLong.I k,
          ∃ a u v w z : k,
            a ≠ 0 ∧
              g =
                Matrix.fromBlocks
                  (Wedge2Formalization.N6PureSingularLong.scalarBlock (k := k) u0)
                  0
                  C0
                  (Wedge2Formalization.N5PureSingularLong.pointwiseShape
                    (k := k) a u v w z) := by
  constructor
  · intro hg
    rcases (mem_K_iff (k := k) (g := g)).1 hg with ⟨u0, C0, H, hH, hEq⟩
    rcases
      (Wedge2Formalization.Paper.N5.Row1.mem_K_iff
        (k := k)
        (g := H)).1 hH with
      ⟨a, u, v, w, z, ha, rfl⟩
    refine ⟨u0, C0, a, a * u, a * v, a * w, a * z, ha, ?_⟩
    rw [hEq]
    congr 4
    simpa [ha] using
      (Wedge2Formalization.Paper.N5.Row1.shape_eq_product
        (k := k)
        (a := a)
        (u := a * u)
        (v := a * v)
        (w := a * w)
        (z := a * z)
        ha).symm
  · rintro ⟨u0, C0, a, u, v, w, z, ha, hEq⟩
    apply (mem_K_iff (k := k) (g := g)).2
    refine ⟨u0, C0, Wedge2Formalization.N5PureSingularLong.pointwiseShape (k := k) a u v w z, ?_, hEq⟩
    exact
      (Wedge2Formalization.Paper.N5.Row1.mem_K_iff
        (k := k)
        (g := Wedge2Formalization.N5PureSingularLong.pointwiseShape
          (k := k) a u v w z)).2
        ⟨a, a⁻¹ * u, a⁻¹ * v, a⁻¹ * w, a⁻¹ * z, ha, by
          simpa [ha] using
            (Wedge2Formalization.Paper.N5.Row1.shape_eq_product
              (k := k)
              (a := a)
              (u := u)
              (v := v)
              (w := w)
              (z := z)
              ha)⟩

theorem U_pointwise
    (C : Matrix Wedge2Formalization.N6PureSingularLong.W
      Wedge2Formalization.N6PureSingularLong.I k)
    (a b c d : k) :
    Wedge2Formalization.N6PureSingularLong.FixesPairBivector
      (U (k := k) C a b c d) := by
  simpa [U] using
    Wedge2Formalization.N6Summary.pureSingularLong_lowerHankel_pointwise_family
      (k := k)
      (u := (1 : k))
      (C := C)
      (a := a)
      (b := b)
      (c := c)
      (d := d)

theorem Levi_pointwise
    (u t : k)
    (ht : t ≠ 0) :
    Wedge2Formalization.N6PureSingularLong.FixesPairBivector
      (Levi (k := k) u t) := by
  simpa [Levi] using
    Wedge2Formalization.N6Summary.pureSingularLong_pointwise_scale_family
      (k := k)
      (u := u)
      (t := t)
      ht
      (C := 0)

theorem quotient_action
    (α β γ δ : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6PureSingularLong.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) α β γ δ)
      (coeff (k := k) α β γ δ) := by
  simpa [ActsOnOrderedPair, rep₁, rep₂, lift, coeff,
    Wedge2Formalization.N5PureSingularLong.Delta, smul_smul, mul_assoc] using
    Wedge2Formalization.N6Summary.pureSingularLong_GL2_lift_action
      (k := k)
      (u := (1 : k))
      (α := α)
      (β := β)
      (γ := γ)
      (δ := δ)
      (C := 0)

private theorem rep_pair_independent
    {a b : k}
    (h :
      a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • (Wedge2Formalization.N5PureSingularLong.rep₁ (k := k)) +
          b • (Wedge2Formalization.N5PureSingularLong.rep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6PureSingularLong.rep₁,
      Wedge2Formalization.N6PureSingularLong.rep₂] using h'
  have h03 := congrArg (fun M => M (Sum.inl 0) (Sum.inr 0)) h''
  have ha : a = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Wedge2Formalization.N5PureSingularLong.mulY,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h03
  have h13 := congrArg (fun M => M (Sum.inl 1) (Sum.inr 0)) h''
  have hb : b = 0 := by
    simpa [ha, Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Wedge2Formalization.N5PureSingularLong.mulY,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h13
  exact ⟨ha, hb⟩

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have h' :
      Wedge2Formalization.N5PureSingularLong.rep₁ (k := k) = 0 := by
    simpa [rep₁, Wedge2Formalization.N6PureSingularLong.rep₁] using
      congrArg Matrix.toBlocks₂₂ hzero
  have h03 := congrArg (fun M => M (Sum.inl 0) (Sum.inr 0)) h'
  simpa [rep₁, Wedge2Formalization.N6PureSingularLong.rep₁,
    Wedge2Formalization.N5PureSingularLong.rep₁,
    Wedge2Formalization.N5PureSingularLong.mulX, Matrix.fromBlocks] using h03

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h' :
      Wedge2Formalization.N5PureSingularLong.rep₂ (k := k) = 0 := by
    simpa [rep₂, Wedge2Formalization.N6PureSingularLong.rep₂] using
      congrArg Matrix.toBlocks₂₂ hzero
  have h13 := congrArg (fun M => M (Sum.inl 1) (Sum.inr 0)) h'
  simpa [rep₂, Wedge2Formalization.N6PureSingularLong.rep₂,
    Wedge2Formalization.N5PureSingularLong.rep₂,
    Wedge2Formalization.N5PureSingularLong.mulY, Matrix.fromBlocks] using h13

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6PureSingularLong.V
        Wedge2Formalization.N6PureSingularLong.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N6PureSingularLong.PreservesSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6PureSingularLong.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N6PureSingularLong.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N6PureSingularLong.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  have hdetcoeff : Matrix.det (!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k) ≠ 0 := by
    intro hdet
    by_cases hleft : α = 0 ∧ γ = 0
    · by_cases hright : β = 0 ∧ δ = 0
      · have hz :
            Wedge2Formalization.N6PureSingularLong.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [hleft.1, hleft.2, hright.1] using h1'
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := rep₁ (k := k)) (g := g) hg).1 hz
        exact rep₁_ne_zero (k := k) horig
      · have hzero :
            Wedge2Formalization.N6PureSingularLong.ActBivector
                (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N6PureSingularLong.ActBivector
                (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g
                =
              δ • Wedge2Formalization.N6PureSingularLong.ActBivector (rep₁ (k := k)) g +
                (-β) • Wedge2Formalization.N6PureSingularLong.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N6PureSingularLong.ActBivector,
                    Matrix.mul_add, Matrix.add_mul]
            _ = 0 := by
                  rw [h1', h2']
                  ext i j
                  simp [hleft.1, hleft.2, Matrix.add_apply, Matrix.smul_apply]
                  ring
        have hcomb_zero :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) (g := g) hg).1 hzero
        have hcomb_ne :
            δ • rep₁ (k := k) + (-β) • rep₂ (k := k) ≠ 0 := by
          intro hcomb
          rcases rep_pair_independent (k := k) hcomb with ⟨hδ, hβ'⟩
          exact hright ⟨by simpa using neg_eq_zero.mp hβ', hδ⟩
        exact hcomb_ne hcomb_zero
    · have hzero :
          Wedge2Formalization.N6PureSingularLong.ActBivector
              (γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) g = 0 := by
        have hdet' : α * δ - β * γ = 0 := by
          simpa [Matrix.det_fin_two] using hdet
        have hcoeff : γ * β = α * δ := by
          simpa [mul_comm] using (sub_eq_zero.mp hdet').symm
        calc
          Wedge2Formalization.N6PureSingularLong.ActBivector
              (γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) g
              =
            γ • Wedge2Formalization.N6PureSingularLong.ActBivector (rep₁ (k := k)) g +
              (-α) • Wedge2Formalization.N6PureSingularLong.ActBivector (rep₂ (k := k)) g := by
                simp [Wedge2Formalization.N6PureSingularLong.ActBivector,
                  Matrix.mul_add, Matrix.add_mul]
          _ = 0 := by
                rw [h1', h2']
                ext i j
                simp [Matrix.add_apply, Matrix.smul_apply]
                have hcoeff_zero :
                    γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ = 0 := by
                  calc
                    γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ
                        = (γ * β - α * δ) * rep₂ (k := k) i j := by ring
                    _ = 0 := by simp [hcoeff]
                calc
                  γ * (α * rep₁ (k := k) i j) + γ * (β * rep₂ (k := k) i j) +
                      (-(α * (γ * rep₁ (k := k) i j)) + -(α * (δ * rep₂ (k := k) i j)))
                      = γ * β * rep₂ (k := k) i j - α * rep₂ (k := k) i j * δ := by ring
                  _ = 0 := hcoeff_zero
      have hcomb_zero :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k) (Ω := γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) (g := g) hg).1 hzero
      have hcomb_ne :
          γ • rep₁ (k := k) + (-α) • rep₂ (k := k) ≠ 0 := by
        intro hcomb
        rcases rep_pair_independent (k := k) hcomb with ⟨hγ, hα'⟩
        exact hleft ⟨by simpa using neg_eq_zero.mp hα', hγ⟩
      exact hcomb_ne hcomb_zero
  refine ⟨!![α, β; γ, δ], hdetcoeff, ?_⟩
  exact ⟨h1', h2'⟩

end Row10

end N6
end Paper
end Wedge2Formalization
