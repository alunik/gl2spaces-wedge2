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

/-! Appendix A, `n = 6`, row 9.
Representative `⟨e₄∧e₆, e₅∧e₆⟩`.
Divisor `0`.
Claimed stabilizer:
`K_L = U_{11} \rtimes (GL_3(k) \times G_m(k))`, exact quotient family
`Q_L = GL_2(k)`.
-/
namespace Row9

def rep₁ :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  Wedge2Formalization.N6PureSingular3.rep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  Wedge2Formalization.N6PureSingular3.rep₂ (k := k)

/-- The paper representative used for the first basis vector in Appendix A, row 9. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  rep₁ (k := k)

/-- The paper representative used for the second basis vector in Appendix A, row 9. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  rep₂ (k := k)

/-- The row-9 paper representative already agrees with the internal working one. -/
def paperChange :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6PureSingular3.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6PureSingular3.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6PureSingular3.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6PureSingular3.ActBivector]

/-- Exact pointwise kernel for the radical pure singular `3 + 3` row. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6PureSingular3.V
        Wedge2Formalization.N6PureSingular3.V k) :=
  { g |
      ∃ A : Matrix Wedge2Formalization.N6PureSingular3.I
          Wedge2Formalization.N6PureSingular3.I k,
        ∃ C : Matrix Wedge2Formalization.N6PureSingular3.W
            Wedge2Formalization.N6PureSingular3.I k,
          ∃ a x y : k,
            a ≠ 0 ∧
              g =
                Matrix.fromBlocks
                  A
                  0
                  C
                  (Wedge2Formalization.N3PureSingular.pureSingularShape
                    (k := k) a x y) }

/-- The pure singular core factor on the bottom-right `3`-space. -/
def Core :
    Set
      (Matrix Wedge2Formalization.N6PureSingular3.W
        Wedge2Formalization.N6PureSingular3.W k) :=
  { H |
      ∃ a x y : k,
        a ≠ 0 ∧
          H = Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y }

/-- The embedded pure singular `3 \times 3` core block on the bottom-right
`3`-space, written in the local coordinates of the pure singular model. -/
def PureSingularCore :
    Set
      (Matrix Wedge2Formalization.N6PureSingular3.W
        Wedge2Formalization.N6PureSingular3.W k) :=
  Core (k := k)

/-- A displayed pointwise family on the radical pure singular `3 + 3` row. -/
def U
    (C : Matrix Wedge2Formalization.N6PureSingular3.W
      Wedge2Formalization.N6PureSingular3.I k)
    (x y : k) :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6PureSingular3.I
      Wedge2Formalization.N6PureSingular3.I k)
    0
    C
    (Wedge2Formalization.N3PureSingular.pointwiseUnipotent (k := k) x y)

/-- The displayed Levi/core family on the radical pure singular `3 + 3` row. -/
def Levi
    (A : Matrix Wedge2Formalization.N6PureSingular3.I
      Wedge2Formalization.N6PureSingular3.I k)
    (a : k) :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  Matrix.fromBlocks
    A
    0
    0
    (Wedge2Formalization.N3PureSingular.pointwiseScale (k := k) a)

/-- Exact coefficient-side quotient family `GL₂(k)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

/-- Chosen `GL₂(k)` lift. -/
def lift
    (α β γ δ : k) :
    Matrix Wedge2Formalization.N6PureSingular3.V
      Wedge2Formalization.N6PureSingular3.V k :=
  let h : Matrix Wedge2Formalization.N6PureSingular3.W
      Wedge2Formalization.N6PureSingular3.W k :=
    Wedge2Formalization.N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6PureSingular3.I
      Wedge2Formalization.N6PureSingular3.I k)
    0
    0
    h

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6PureSingular3.FixesPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N6Summary.pureSingular3_pointwise_bivector_iff
        (k := k)
        (A := g.toBlocks₁₁)
        (B := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hg) with
      ⟨hB, hD00, hshape⟩
    let D := g.toBlocks₂₂
    have hshape' :
        D =
          (Wedge2Formalization.N3PureSingular.pureSingularShape
            (k := k) (D 0 0) (D 2 0) (D 2 1)) := by
      simpa [D] using hshape
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, D 0 0, D 2 0, D 2 1, ?_, ?_⟩
    · simpa [D] using hD00
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ D := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ =
          Matrix.fromBlocks
            g.toBlocks₁₁
            0
            g.toBlocks₂₁
            D := by
            rw [hB]
      _ =
          Matrix.fromBlocks
            g.toBlocks₁₁
            0
            g.toBlocks₂₁
            (Wedge2Formalization.N3PureSingular.pureSingularShape
              (k := k) (D 0 0) (D 2 0) (D 2 1)) := by
            exact congrArg
              (fun X =>
                Matrix.fromBlocks
                  g.toBlocks₁₁
                  0
                  g.toBlocks₂₁
                  X) hshape'
  · rintro ⟨A, C, a, x, y, ha, rfl⟩
    exact
      (Wedge2Formalization.N6Summary.pureSingular3_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := 0)
        (C := C)
        (D := Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)).2
        ⟨rfl, ha, rfl⟩

theorem U_pointwise
    (C : Matrix Wedge2Formalization.N6PureSingular3.W
      Wedge2Formalization.N6PureSingular3.I k)
    (x y : k) :
    Wedge2Formalization.N6PureSingular3.FixesPairBivector (U (k := k) C x y) := by
  simpa [U] using
    Wedge2Formalization.N6Summary.pureSingular3_pointwise_unipotent_family
      (k := k)
      (A0 := (1 : Matrix Wedge2Formalization.N6PureSingular3.I
        Wedge2Formalization.N6PureSingular3.I k))
      (C := C)
      (x := x)
      (y := y)

theorem Levi_pointwise
    (A : Matrix Wedge2Formalization.N6PureSingular3.I
      Wedge2Formalization.N6PureSingular3.I k)
    (a : k)
    (ha : a ≠ 0) :
    Wedge2Formalization.N6PureSingular3.FixesPairBivector (Levi (k := k) A a) := by
  simpa [Levi] using
    Wedge2Formalization.N6Summary.pureSingular3_pointwise_scale_family
      (k := k)
      (A0 := A)
      (a := a)
      ha
      (C := 0)

theorem mem_K_shape_iff
    (g :
      Matrix Wedge2Formalization.N6PureSingular3.V
        Wedge2Formalization.N6PureSingular3.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6PureSingular3.I
          Wedge2Formalization.N6PureSingular3.I k,
        ∃ C : Matrix Wedge2Formalization.N6PureSingular3.W
            Wedge2Formalization.N6PureSingular3.I k,
          ∃ a x y : k,
            a ≠ 0 ∧
              g =
                Matrix.fromBlocks
                  A
                  0
                  C
                  (Wedge2Formalization.N3PureSingular.pureSingularShape
                    (k := k) a x y) := by
  rfl

theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6PureSingular3.V
        Wedge2Formalization.N6PureSingular3.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6PureSingular3.I
          Wedge2Formalization.N6PureSingular3.I k,
        ∃ C : Matrix Wedge2Formalization.N6PureSingular3.W
            Wedge2Formalization.N6PureSingular3.I k,
          ∃ H : Matrix Wedge2Formalization.N6PureSingular3.W
              Wedge2Formalization.N6PureSingular3.W k,
            H ∈ Core (k := k) ∧
              g = Matrix.fromBlocks A 0 C H := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with ⟨A, C, a, x, y, ha, hEq⟩
    refine ⟨A, C, Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y,
      ?_, hEq⟩
    exact ⟨a, x, y, ha, rfl⟩
  · rintro ⟨A, C, H, hH, hEq⟩
    apply (mem_K_shape_iff (k := k) (g := g)).2
    rcases hH with ⟨a, x, y, ha, rfl⟩
    exact ⟨A, C, a, x, y, ha, hEq⟩

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_{11} \rtimes (GL_3(k) \times G_m(k))`, written in the explicit local
coordinates of this row rather than through a nested pure-singular core block. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6PureSingular3.V
        Wedge2Formalization.N6PureSingular3.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6PureSingular3.I
          Wedge2Formalization.N6PureSingular3.I k,
        ∃ C : Matrix Wedge2Formalization.N6PureSingular3.W
            Wedge2Formalization.N6PureSingular3.I k,
          ∃ a x y : k,
            a ≠ 0 ∧
              g =
                Matrix.fromBlocks
                  A
                  0
                  C
                  (Wedge2Formalization.N3PureSingular.pureSingularShape
                    (k := k) a x y) :=
  mem_K_shape_iff (k := k) (g := g)

theorem quotient_action
    (α β γ δ : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6PureSingular3.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) α β γ δ)
      (!![α, β; γ, δ]) := by
  simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
    Wedge2Formalization.N6Summary.pureSingular3_GL2_lift_action
      (k := k)
      (A0 := (1 : Matrix Wedge2Formalization.N6PureSingular3.I
        Wedge2Formalization.N6PureSingular3.I k))
      (α := α)
      (β := β)
      (γ := γ)
      (δ := δ)
      (C := 0)

private theorem rep_pair_independent_row9
    {a b : k}
    (h :
      a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h02 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 2)) h
  have ha : a = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6PureSingular3.rep₁,
      Wedge2Formalization.N6PureSingular3.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h02
  have h12 := congrArg (fun M => M (Sum.inr 1) (Sum.inr 2)) h
  have hb : b = 0 := by
    simpa [ha, rep₁, rep₂, Wedge2Formalization.N6PureSingular3.rep₁,
      Wedge2Formalization.N6PureSingular3.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h12
  exact ⟨ha, hb⟩

private theorem rep₂_ne_zero_row9 : rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h02 := congrArg (fun M => M (Sum.inr 1) (Sum.inr 2)) hzero
  simpa [rep₂, Wedge2Formalization.N6PureSingular3.rep₂,
    Wedge2Formalization.N3PureSingular.ω23, Matrix.fromBlocks] using h02

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6PureSingular3.V
        Wedge2Formalization.N6PureSingular3.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N6PureSingular3.PreservesSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6PureSingular3.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have h1' :
      Wedge2Formalization.N6PureSingular3.ActBivector (rep₁ (k := k)) g =
        α • rep₁ (k := k) + β • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h1
  have h2' :
      Wedge2Formalization.N6PureSingular3.ActBivector (rep₂ (k := k)) g =
        γ • rep₁ (k := k) + δ • rep₂ (k := k) := by
    simpa [rep₁, rep₂] using h2
  have hdetcoeff : Matrix.det (!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k) ≠ 0 := by
    intro hdet
    by_cases hleft : α = 0 ∧ γ = 0
    · by_cases hright : β = 0 ∧ δ = 0
      · have hz :
          Wedge2Formalization.N6PureSingular3.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [hleft.1, hleft.2, hright.1] using h1'
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := rep₁ (k := k)) (g := g) hg).1 hz
        have hrep₁nz : rep₁ (k := k) ≠ 0 := by
          intro hzero
          have h02 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 2)) hzero
          simpa [rep₁, Wedge2Formalization.N6PureSingular3.rep₁,
            Wedge2Formalization.N3PureSingular.ω13, Matrix.fromBlocks] using h02
        exact hrep₁nz horig
      · have hzero :
          Wedge2Formalization.N6PureSingular3.ActBivector
              (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N6PureSingular3.ActBivector
                (δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) g
                =
              δ • Wedge2Formalization.N6PureSingular3.ActBivector (rep₁ (k := k)) g +
                (-β) • Wedge2Formalization.N6PureSingular3.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N6PureSingular3.ActBivector,
                    Matrix.mul_add, Matrix.add_mul]
            _ = 0 := by
                  rw [h1', h2']
                  ext i j <;> simp [hleft.1, hleft.2, Matrix.add_apply, Matrix.smul_apply] <;> ring
        have hcomb_zero :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := δ • rep₁ (k := k) + (-β) • rep₂ (k := k)) (g := g) hg).1 hzero
        have hcomb_ne :
            δ • rep₁ (k := k) + (-β) • rep₂ (k := k) ≠ 0 := by
          intro hcomb
          rcases rep_pair_independent_row9 (k := k) hcomb with ⟨hδ, hβ'⟩
          exact hright ⟨by simpa using neg_eq_zero.mp hβ', hδ⟩
        exact hcomb_ne hcomb_zero
    · have hzero :
        Wedge2Formalization.N6PureSingular3.ActBivector
            (γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) g = 0 := by
        have hdet' : α * δ - β * γ = 0 := by
          simpa [Matrix.det_fin_two] using hdet
        have hcoeff : γ * β = α * δ := by
          simpa [mul_comm] using (sub_eq_zero.mp hdet').symm
        calc
          Wedge2Formalization.N6PureSingular3.ActBivector
              (γ • rep₁ (k := k) + (-α) • rep₂ (k := k)) g
              =
            γ • Wedge2Formalization.N6PureSingular3.ActBivector (rep₁ (k := k)) g +
              (-α) • Wedge2Formalization.N6PureSingular3.ActBivector (rep₂ (k := k)) g := by
                simp [Wedge2Formalization.N6PureSingular3.ActBivector,
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
        rcases rep_pair_independent_row9 (k := k) hcomb with ⟨hγ, hα'⟩
        exact hleft ⟨by simpa using neg_eq_zero.mp hα', hγ⟩
      exact hcomb_ne hcomb_zero
  refine ⟨!![α, β; γ, δ], hdetcoeff, ?_⟩
  exact ⟨h1', h2'⟩

end Row9

end N6
end Paper
end Wedge2Formalization
