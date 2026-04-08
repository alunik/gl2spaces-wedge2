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

/-! Appendix A, `n = 6`, row 7.
Representative `⟨e₃∧e₄ + e₅∧e₆, e₅∧e₆⟩`.
Divisor `[a]+[b]`.
Claimed stabilizer:
`K_L = U_8 \rtimes (GL_2(k) \times SL_2(k) \times SL_2(k))`, exact quotient
family `Q_L = N_T`.
-/
namespace Row7

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
  Wedge2Formalization.N6.rad2SplitRep₁ (k := k)

def rep₂ : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Wedge2Formalization.N6.rad2SplitRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 7:
`e₃∧e₄ + e₅∧e₆`. -/
def paperRep₁ : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row2.rep₁ (k := k))

/-- Literal second basis vector from Appendix A, row 7:
`e₅∧e₆`. -/
def paperRep₂ : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N4.Row2.rep₂ (k := k))

/-- The paper representative already agrees with the internal working representative. -/
def paperChange : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  1

/-- Transport of the literal paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6.rad2SplitRep₁,
    Wedge2Formalization.Paper.N4.Row2.rep₁, Wedge2Formalization.N4.splitRep₁,
    Wedge2Formalization.N6.ActBivector]

/-- Transport of the literal paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6.rad2SplitRep₂,
    Wedge2Formalization.Paper.N4.Row2.rep₂, Wedge2Formalization.N4.splitRep₂,
    Wedge2Formalization.N6.ActBivector]

/-- Exact pointwise kernel family for the radical split row. -/
def K : Set (Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :=
  { g |
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g =
              Matrix.fromBlocks
                A0
                0
                C
                (Matrix.fromBlocks A 0 0 E) }

/-- The displayed radical unipotent family on the radical split row. -/
def U
    (C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    0
    C
    (Matrix.fromBlocks
      (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
      0
      0
      (1 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))

/-- The displayed Levi/core family on the radical split row. -/
def Levi
    (A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    (A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks
    A0
    0
    0
    (Matrix.fromBlocks A 0 0 E)

/-- Exact coefficient-side split normalizer family in the ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N4.Row2.Qproj (k := k)

/-- Chosen torus lift. -/
def torusLift
    (u v : k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    0
    0
    (Matrix.fromBlocks
      (!![u, 0; 0, 1] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
      0 0
      (!![v, 0; 0, 1] : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k))

/-- Chosen swap-coset lift. -/
def swapLift
    (u v : k) :
    Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k :=
  let A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k := !![u, 0; 0, 1]
  let E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k := !![v, 0; 0, 1]
  Matrix.fromBlocks
    (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    0
    0
    (Matrix.fromBlocks A 0 0 E * Wedge2Formalization.N4.splitSwap (k := k))

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6.FixesRad2SplitPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N6Summary.rad2Split_pointwise_bivector_iff
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
      ⟨hR, hA, hB, hC, hD⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, g.toBlocks₂₂.toBlocks₁₁, g.toBlocks₂₂.toBlocks₂₂,
      hA, hD, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := by
            rw [hR]
      _ =
          Matrix.fromBlocks
            g.toBlocks₁₁
            0
            g.toBlocks₂₁
            (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ 0 0 g.toBlocks₂₂.toBlocks₂₂) := by
            have hEq22 :
                Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ 0 0 g.toBlocks₂₂.toBlocks₂₂ =
                  g.toBlocks₂₂ := by
              calc
                Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ 0 0 g.toBlocks₂₂.toBlocks₂₂
                    =
                  Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₁₂
                    g.toBlocks₂₂.toBlocks₂₁ g.toBlocks₂₂.toBlocks₂₂ := by
                      rw [hB, hC]
                _ = g.toBlocks₂₂ := by
                      exact Matrix.fromBlocks_toBlocks g.toBlocks₂₂
            rw [hEq22]
  · rintro ⟨A0, C, A, E, hA, hE, rfl⟩
    exact
      (Wedge2Formalization.N6Summary.rad2Split_pointwise_bivector_iff
        (k := k) (A0 := A0) (R := 0) (C := C) (A := A) (B₁ := 0) (C₁ := 0) (D := E)).2
        ⟨rfl, hA, rfl, rfl, hE⟩

theorem U_pointwise
    (C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k) :
    Wedge2Formalization.N6.FixesRad2SplitPairBivector (U (k := k) C) := by
  exact
    (Wedge2Formalization.N6Summary.rad2Split_pointwise_bivector_iff
      (k := k)
      (A0 := 1)
      (R := 0)
      (C := C)
      (A := 1)
      (B₁ := 0)
      (C₁ := 0)
      (D := 1)).2
      ⟨rfl, by simp, rfl, rfl, by simp⟩

theorem Levi_pointwise
    (A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    (A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1)
    (hE : E.det = 1) :
    Wedge2Formalization.N6.FixesRad2SplitPairBivector (Levi (k := k) A0 A E) := by
  exact
    (Wedge2Formalization.N6Summary.rad2Split_pointwise_bivector_iff
      (k := k)
      (A0 := A0)
      (R := 0)
      (C := 0)
      (A := A)
      (B₁ := 0)
      (C₁ := 0)
      (D := E)).2
      ⟨rfl, hA, rfl, rfl, hE⟩

theorem mem_K_shape_iff
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g =
              Matrix.fromBlocks
                A0
                0
                C
                (Matrix.fromBlocks A 0 0 E) := by
  rfl

/-- The displayed Levi family and radical family multiply to the raw pointwise block
shape. -/
private theorem Levi_mul_U_eq_shape
    (A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k)
    (C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k)
    (A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Levi (k := k) A0 A E * U (k := k) C =
      Matrix.fromBlocks
        A0
        0
        (Matrix.fromBlocks A 0 0 E * C)
        (Matrix.fromBlocks A 0 0 E) := by
  rw [U, Levi, Matrix.fromBlocks_multiply]
  ext i j
  fin_cases i <;> fin_cases j <;> simp

/-- Public kernel-membership statement in the displayed paper notation. -/
theorem mem_K_iff
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g = Levi (k := k) A0 A E * U (k := k) C := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with ⟨A0, C, A, E, hA, hE, hEq⟩
    let H : Matrix Wedge2Formalization.N4.V Wedge2Formalization.N4.V k :=
      Matrix.fromBlocks A 0 0 E
    have hdetH : Matrix.det H = 1 := by
      dsimp [H]
      rw [Matrix.det_fromBlocks_zero₂₁]
      simp [hA, hE]
    have hunitH : IsUnit (Matrix.det H) := by rw [hdetH]; exact isUnit_one
    have hHC : H * (H⁻¹ * C) = C := by
      simpa using Matrix.mul_nonsing_inv_cancel_left H C hunitH
    refine ⟨A0, H⁻¹ * C, A, E, hA, hE, ?_⟩
    rw [Levi_mul_U_eq_shape (k := k) (A0 := A0) (C := H⁻¹ * C) (A := A) (E := E)]
    simpa [H, hHC] using hEq
  · rintro ⟨A0, C, A, E, hA, hE, hEq⟩
    apply (mem_K_shape_iff (k := k) (g := g)).2
    refine ⟨A0, Matrix.fromBlocks A 0 0 E * C, A, E, hA, hE, ?_⟩
    rw [hEq, Levi_mul_U_eq_shape (k := k) (A0 := A0) (C := C) (A := A) (E := E)]

/-- Table-facing kernel statement for the Appendix A row
`K_L = U_8 \rtimes (GL_2(k) \times SL_2(k) \times SL_2(k))`. -/
theorem mem_K_table_iff
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k,
        ∃ C : Matrix Wedge2Formalization.N6.W Wedge2Formalization.N6.I k,
          ∃ A E : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g = Levi (k := k) A0 A E * U (k := k) C :=
  mem_K_iff (k := k) (g := g)

theorem torus_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (torusLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_torusCoeff (k := k) u v) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, torusLift] using
      (Wedge2Formalization.N6Summary.rad2Split_torus_lift_action
        (k := k) (A0 := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
        (u := u) (v := v) (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, torusLift,
      Wedge2Formalization.N4PaperSummary.row2_torusCoeff] using
      (Wedge2Formalization.N6Summary.rad2Split_torus_lift_action
        (k := k) (A0 := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
        (u := u) (v := v) (C := 0)).2

theorem swap_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (swapLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_swapCoeff (k := k) u v) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, swapLift] using
      (Wedge2Formalization.N6Summary.rad2Split_swapCoset_lift_action
        (k := k) (A0 := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
        (u := u) (v := v) (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, swapLift,
      Wedge2Formalization.N4PaperSummary.row2_swapCoeff] using
      (Wedge2Formalization.N6Summary.rad2Split_swapCoset_lift_action
        (k := k) (A0 := (1 : Matrix Wedge2Formalization.N6.I Wedge2Formalization.N6.I k))
        (u := u) (v := v) (C := 0)).2

/-- Rank drop on the split-support quotient pencil occurs exactly at the two support points. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det
      (a • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        b • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  Wedge2Formalization.N4Summary.split_det_zero_iff (k := k) (a := a) (b := b)

private theorem rep_pair_linearCombination_ne_zero
    {a b : k}
    (h : a ≠ 0 ∨ b ≠ 0) :
    a • rep₁ (k := k) + b • rep₂ (k := k) ≠ 0 := by
  intro hzero
  have h' := congrArg Matrix.toBlocks₂₂ hzero
  have h'' :
      a • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        b • (Wedge2Formalization.N4.splitRep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6.rad2SplitRep₁, Wedge2Formalization.N6.rad2SplitRep₂] using h'
  have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h''
  have ha : a = 0 := by
    simpa [Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
      Matrix.add_apply, Matrix.smul_apply] using h01
  have h23 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h''
  have hb : b = 0 := by
    simpa [ha, Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
      Matrix.add_apply, Matrix.smul_apply] using h23
  exact h.elim (fun ha' => ha' ha) (fun hb' => hb' hb)

private theorem rep₂_ne_zero : rep₂ (k := k) ≠ 0 := by
  have h :=
    rep_pair_linearCombination_ne_zero (k := k) (a := (0 : k)) (b := (1 : k)) (Or.inr one_ne_zero)
  simpa using h

private theorem lowerRight_action
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k)
    {α β γ δ : k}
    (h1 :
      Wedge2Formalization.N6.ActBivector
        (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) g
        =
      α • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
        β • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k)))
    (h2 :
      Wedge2Formalization.N6.ActBivector
        (Wedge2Formalization.N6.rad2SplitRep₂ (k := k)) g
        =
      γ • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
        δ • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k))) :
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₁ (k := k)) g.toBlocks₂₂ =
        α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.splitRep₂ (k := k)) ∧
    Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₂ (k := k)) g.toBlocks₂₂ =
        γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
  have h1' := congrArg Matrix.toBlocks₂₂ h1
  have h2' := congrArg Matrix.toBlocks₂₂ h2
  have hrhs1 :
      Matrix.toBlocks₂₂
        (α • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
          β • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k))) =
        α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          β • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    ext i j
    change
      (α • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
        β • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        β • (Wedge2Formalization.N4.splitRep₂ (k := k))) i j
    simp [Wedge2Formalization.N6.rad2SplitRep₁, Wedge2Formalization.N6.rad2SplitRep₂,
      Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
    ring
  have hrhs2 :
      Matrix.toBlocks₂₂
        (γ • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k))) =
        γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    ext i j
    change
      (γ • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
        δ • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k))) (Sum.inr i) (Sum.inr j)
        =
      (γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        δ • (Wedge2Formalization.N4.splitRep₂ (k := k))) i j
    simp [Wedge2Formalization.N6.rad2SplitRep₁, Wedge2Formalization.N6.rad2SplitRep₂,
      Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
    ring
  constructor
  · calc
      Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₁ (k := k)) g.toBlocks₂₂
          =
        Matrix.toBlocks₂₂
          (Wedge2Formalization.N6.ActBivector
            (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) g) := by
              symm
              simpa [Wedge2Formalization.N6.rad2SplitRep₁] using
                act_embedded_toBlocks₂₂_N6
                  (k := k)
                  (Ω := Wedge2Formalization.N4.splitRep₁ (k := k))
                  (g := g)
      _ = Matrix.toBlocks₂₂
            (α • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
              β • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k))) := by
                simpa using h1'
      _ =
          α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            β • (Wedge2Formalization.N4.splitRep₂ (k := k)) := hrhs1
  · calc
      Wedge2Formalization.N4.ActBivector (Wedge2Formalization.N4.splitRep₂ (k := k)) g.toBlocks₂₂
          =
        Matrix.toBlocks₂₂
          (Wedge2Formalization.N6.ActBivector
            (Wedge2Formalization.N6.rad2SplitRep₂ (k := k)) g) := by
              symm
              simpa [Wedge2Formalization.N6.rad2SplitRep₂] using
                act_embedded_toBlocks₂₂_N6
                  (k := k)
                  (Ω := Wedge2Formalization.N4.splitRep₂ (k := k))
                  (g := g)
      _ = Matrix.toBlocks₂₂
            (γ • (Wedge2Formalization.N6.rad2SplitRep₁ (k := k)) +
              δ • (Wedge2Formalization.N6.rad2SplitRep₂ (k := k))) := by
                simpa using h2'
      _ =
          γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            δ • (Wedge2Formalization.N4.splitRep₂ (k := k)) := hrhs2

theorem quotient_image
    (g : Matrix Wedge2Formalization.N6.V Wedge2Formalization.N6.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      Wedge2Formalization.N6.PreservesRad2SplitSubspaceBivector (k := k) g) :
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
  have hrep2_det :
      Matrix.det (Wedge2Formalization.N4.splitRep₂ (k := k)) = 0 := by
    simpa using
      (Wedge2Formalization.N4Summary.split_det_zero_iff
        (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl)
  have hrep2_det_zero :
      Matrix.det
        (γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          δ • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    rw [← hD2, Wedge2Formalization.N4Summary.actBivector_det, hrep2_det]
    simp
  have hω12_det :
      Matrix.det
        ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
          (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    have h :=
      (Wedge2Formalization.N4Summary.split_det_zero_iff
        (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (by ring))
    simpa [sub_eq_add_neg] using h
  have hω12_action :
      Wedge2Formalization.N4.ActBivector
        ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
          (Wedge2Formalization.N4.splitRep₂ (k := k))) g.toBlocks₂₂
        =
      (α - γ) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
        (β - δ) • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    calc
      Wedge2Formalization.N4.ActBivector
          ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
            (Wedge2Formalization.N4.splitRep₂ (k := k))) g.toBlocks₂₂
          =
        Wedge2Formalization.N4.ActBivector
          ((Wedge2Formalization.N4.splitRep₁ (k := k)) +
            (-1 : k) • (Wedge2Formalization.N4.splitRep₂ (k := k))) g.toBlocks₂₂ := by
              simp [sub_eq_add_neg]
      _ =
          Wedge2Formalization.N4.ActBivector
              (Wedge2Formalization.N4.splitRep₁ (k := k)) g.toBlocks₂₂ +
            (-1 : k) •
              Wedge2Formalization.N4.ActBivector
                (Wedge2Formalization.N4.splitRep₂ (k := k)) g.toBlocks₂₂ := by
                  rw [Wedge2Formalization.N4Summary.actBivector_add,
                    Wedge2Formalization.N4Summary.actBivector_smul]
      _ =
          (α • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            β • (Wedge2Formalization.N4.splitRep₂ (k := k))) +
          (-1 : k) •
            (γ • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
              δ • (Wedge2Formalization.N4.splitRep₂ (k := k))) := by
                rw [hD1, hD2]
      _ =
          (α - γ) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
            (β - δ) • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
              ext i j
              simp [sub_eq_add_neg]
              ring
  have hω12_det_zero :
      Matrix.det
        ((α - γ) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          (β - δ) • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    rw [← hω12_action, Wedge2Formalization.N4Summary.actBivector_det, hω12_det]
    simp
  have hcoeff_det_ne : α * δ - β * γ ≠ 0 := by
    intro hdet
    by_cases hbd : β = 0 ∧ δ = 0
    · by_cases hca : γ = 0 ∧ α = 0
      · have hrep2_zero :
          Wedge2Formalization.N6.ActBivector (rep₂ (k := k)) g = 0 := by
          simpa [hca.1, hbd.2] using h2'
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (Ω := rep₂ (k := k)) (g := g) hg).1 hrep2_zero
        exact rep₂_ne_zero (k := k) horig
      · have hnonzero : γ ≠ 0 ∨ α ≠ 0 := by
          by_cases hγ : γ = 0
          · right
            intro hα
            exact hca ⟨hγ, hα⟩
          · exact Or.inl hγ
        have hzero2 :
            Wedge2Formalization.N6.ActBivector
              (γ • rep₁ (k := k) - α • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N6.ActBivector
                (γ • rep₁ (k := k) - α • rep₂ (k := k)) g
                =
              γ • Wedge2Formalization.N6.ActBivector (rep₁ (k := k)) g +
                (-α) • Wedge2Formalization.N6.ActBivector (rep₂ (k := k)) g := by
                    simp [sub_eq_add_neg, Wedge2Formalization.N6.ActBivector,
                      Matrix.mul_add, Matrix.add_mul]
            _ =
                γ • (α • rep₁ (k := k) + β • rep₂ (k := k)) +
                  (-α) • (γ • rep₁ (k := k) + δ • rep₂ (k := k)) := by
                    rw [h1', h2']
            _ = 0 := by
                  ext i j
                  simp [hbd.1, hbd.2, sub_eq_add_neg]
                  ring
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (Ω := γ • rep₁ (k := k) - α • rep₂ (k := k)) (g := g) hg).1 hzero2
        have hcomb_ne :
            γ • rep₁ (k := k) - α • rep₂ (k := k) ≠ 0 := by
          have hcomb_ne' :
              γ ≠ 0 ∨ -α ≠ 0 := by
            rcases hnonzero with hγ | hα
            · exact Or.inl hγ
            · exact Or.inr (neg_ne_zero.mpr hα)
          have := rep_pair_linearCombination_ne_zero (k := k) (a := γ) (b := -α) hcomb_ne'
          simpa [sub_eq_add_neg] using this
        exact hcomb_ne horig
    · have hnonzero : β ≠ 0 ∨ δ ≠ 0 := by
        by_cases hβ : β = 0
        · right
          intro hδ
          exact hbd ⟨hβ, hδ⟩
        · exact Or.inl hβ
      have hzero1 :
          Wedge2Formalization.N6.ActBivector
            (δ • rep₁ (k := k) - β • rep₂ (k := k)) g = 0 := by
        calc
          Wedge2Formalization.N6.ActBivector
              (δ • rep₁ (k := k) - β • rep₂ (k := k)) g
              =
            δ • Wedge2Formalization.N6.ActBivector (rep₁ (k := k)) g +
              (-β) • Wedge2Formalization.N6.ActBivector (rep₂ (k := k)) g := by
                simp [sub_eq_add_neg, Wedge2Formalization.N6.ActBivector,
                  Matrix.mul_add, Matrix.add_mul]
          _ =
              δ • (α • rep₁ (k := k) + β • rep₂ (k := k)) +
                (-β) • (γ • rep₁ (k := k) + δ • rep₂ (k := k)) := by
                  rw [h1', h2']
          _ = (α * δ - β * γ) • rep₁ (k := k) := by
                ext i j
                simp [sub_eq_add_neg]
                ring
          _ = 0 := by
                simp [hdet]
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (Ω := δ • rep₁ (k := k) - β • rep₂ (k := k)) (g := g) hg).1 hzero1
      have hcomb_ne :
          δ • rep₁ (k := k) - β • rep₂ (k := k) ≠ 0 := by
        have hcomb_ne' :
            δ ≠ 0 ∨ -β ≠ 0 := by
          rcases hnonzero with hβ | hδ
          · exact Or.inr (neg_ne_zero.mpr hβ)
          · exact Or.inl hδ
        have := rep_pair_linearCombination_ne_zero (k := k) (a := δ) (b := -β) hcomb_ne'
        simpa [sub_eq_add_neg] using this
      exact hcomb_ne horig
  have hM_mem :
      ((!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k)) ∈ Qproj (k := k) := by
    rcases
      Wedge2Formalization.N4Summary.split_action_cases_of_rankDrop
        (k := k)
        (α := α)
        (β := β)
        (γ := γ)
        (δ := δ)
        hcoeff_det_ne
        hrep2_det_zero
        hω12_det_zero with
      hdiag | hswap
    · rcases hdiag with ⟨hγ, hsum⟩
      left
      have hbeta : β = δ - α := by
        calc
          β = α + β - α := by ring
          _ = δ - α := by rw [hsum]
      refine ⟨α, δ, ?_⟩
      ext i j
      fin_cases i <;> fin_cases j
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, hbeta]
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, hbeta]
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, hγ]
      · simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff]
    · rcases hswap with ⟨hαγ, hγδ⟩
      right
      refine ⟨α, α + β, ?_⟩
      have hδ : δ = -α := by
        calc
          δ = γ + δ - γ := by ring
          _ = 0 - γ := by rw [hγδ]
          _ = -γ := by ring
          _ = -α := by rw [hαγ]
      ext i j
      fin_cases i <;> fin_cases j
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff]
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff]
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff, hαγ]
      · simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff, hδ]
  refine ⟨!![α, β; γ, δ], hM_mem, ?_⟩
  constructor
  · simpa [ActsOnOrderedPair] using h1
  · simpa [ActsOnOrderedPair] using h2

end Row7

end N6
end Paper
end Wedge2Formalization
