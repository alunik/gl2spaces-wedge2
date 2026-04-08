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

/-! Appendix A, `n = 6`, row 5.
Representative `⟨e₁∧e₂ + e₃∧e₄, e₃∧e₄ + e₅∧e₆⟩`.
Divisor `[a]+[b]+[c]`.
Claimed stabilizer:
`K_L = SL_2(k) \times SL_2(k) \times SL_2(k)`, `Q_L =` the preimage of `S_3`
in `GL_2(k)`.
-/
namespace Row5

def rep₁ : Matrix Wedge2Formalization.N6ThreePoint.V Wedge2Formalization.N6ThreePoint.V k :=
  Wedge2Formalization.N6ThreePoint.rep₁ (k := k)

def rep₂ : Matrix Wedge2Formalization.N6ThreePoint.V Wedge2Formalization.N6ThreePoint.V k :=
  Wedge2Formalization.N6ThreePoint.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 5. -/
def paperRep₁ : Matrix Wedge2Formalization.N6ThreePoint.V Wedge2Formalization.N6ThreePoint.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 5. -/
def paperRep₂ : Matrix Wedge2Formalization.N6ThreePoint.V Wedge2Formalization.N6ThreePoint.V k :=
  rep₂ (k := k)

/-- The row-5 paper representative already agrees with the internal working pair. -/
def paperChange : Matrix Wedge2Formalization.N6ThreePoint.V Wedge2Formalization.N6ThreePoint.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6ThreePoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6ThreePoint.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6ThreePoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6ThreePoint.ActBivector]

/-- The displayed block-diagonal pointwise family on the three-point row. -/
def Levi
    (A : Matrix Wedge2Formalization.N6ThreePoint.W Wedge2Formalization.N6ThreePoint.W k)
    (D : Matrix Wedge2Formalization.N6ThreePoint.I Wedge2Formalization.N6ThreePoint.I k) :
    Matrix Wedge2Formalization.N6ThreePoint.V Wedge2Formalization.N6ThreePoint.V k :=
  Matrix.fromBlocks A 0 0 D

/-- Exact pointwise kernel family for the three-point row. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6ThreePoint.V
        Wedge2Formalization.N6ThreePoint.V k) :=
  { g |
      ∃ A : Matrix Wedge2Formalization.N6ThreePoint.W
          Wedge2Formalization.N6ThreePoint.W k,
        ∃ D : Matrix Wedge2Formalization.N6ThreePoint.I
            Wedge2Formalization.N6ThreePoint.I k,
          Wedge2Formalization.N4.FixesSplitPairBivector A ∧
            D.det = 1 ∧
            g = Levi (k := k) A D }

/-- The embedded split-support core, i.e. the `SL_2(k) \times SL_2(k)` factor
coming from Appendix A, `n = 4`, row 2. -/
def SplitCore :
    Set
      (Matrix Wedge2Formalization.N6ThreePoint.W
        Wedge2Formalization.N6ThreePoint.W k) :=
  { A | Wedge2Formalization.N4.FixesSplitPairBivector A }

/-- The three distinguished support conditions on the three-point projective line. -/
def supportCondition (a b : k) : Prop :=
  a = 0 ∨ a + b = 0 ∨ b = 0

/-- Exact coefficient-side quotient family preserving the three support points in the
chosen ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M |
      Matrix.det M ≠ 0 ∧
        supportCondition (k := k) (M 0 0) (M 0 1) ∧
        supportCondition (k := k) (M 1 0) (M 1 1) ∧
        supportCondition (k := k) (M 0 0 - M 1 0) (M 0 1 - M 1 1) }

/-- The nontrivial quotient permutation coming from the split swap on the `4`-block. -/
def swapCoeff : Matrix (Fin 2) (Fin 2) k := !![(1 : k), 0; 1, -1]

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6ThreePoint.FixesPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N6Summary.threePoint_pointwise_bivector_iff
        (k := k) (A := g.toBlocks₁₁) (B := g.toBlocks₁₂) (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hg) with
      ⟨hB, hC, hA, hD⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₂, hA, hD, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Levi (k := k) g.toBlocks₁₁ g.toBlocks₂₂ := by
            simp [Levi, hB, hC]
  · rintro ⟨A, D, hA, hD, rfl⟩
    exact
      (Wedge2Formalization.N6Summary.threePoint_pointwise_bivector_iff
        (k := k) (A := A) (B := 0) (C := 0) (D := D)).2 ⟨rfl, rfl, hA, hD⟩

theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6ThreePoint.V
        Wedge2Formalization.N6ThreePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6ThreePoint.W
          Wedge2Formalization.N6ThreePoint.W k,
        ∃ D : Matrix Wedge2Formalization.N6ThreePoint.I
            Wedge2Formalization.N6ThreePoint.I k,
          Wedge2Formalization.N4.FixesSplitPairBivector A ∧
            D.det = 1 ∧
            g = Levi (k := k) A D := by
  rfl

/-- Table-facing kernel statement for the Appendix A row
`K_L = SL_2(k) \times SL_2(k) \times SL_2(k)`, written directly as a
block-diagonal triple-`SL_2` family. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6ThreePoint.V
        Wedge2Formalization.N6ThreePoint.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N6ThreePoint.W
          Wedge2Formalization.N6ThreePoint.W k,
        ∃ D : Matrix Wedge2Formalization.N6ThreePoint.I
            Wedge2Formalization.N6ThreePoint.I k,
          Wedge2Formalization.N4.FixesSplitPairBivector A ∧
            D.det = 1 ∧
            g = Levi (k := k) A D :=
  mem_K_iff (k := k) (g := g)

theorem Levi_pointwise
    (A : Matrix Wedge2Formalization.N6ThreePoint.W
      Wedge2Formalization.N6ThreePoint.W k)
    (D : Matrix Wedge2Formalization.N6ThreePoint.I
      Wedge2Formalization.N6ThreePoint.I k)
    (hA : Wedge2Formalization.N4.FixesSplitPairBivector A)
    (hD : D.det = 1) :
    Wedge2Formalization.N6ThreePoint.FixesPairBivector (Levi (k := k) A D) := by
  simpa [Levi] using
    Wedge2Formalization.N6Summary.threePoint_pointwise_family
      (k := k) (A := A) (D := D) hA hD

theorem swap_action :
    let E : Matrix Wedge2Formalization.N6ThreePoint.I
        Wedge2Formalization.N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
    let g : Matrix Wedge2Formalization.N6ThreePoint.V
        Wedge2Formalization.N6ThreePoint.V k :=
      Matrix.fromBlocks (Wedge2Formalization.N4.splitSwap (k := k)) 0 0 E
    ActsOnOrderedPair
      (Wedge2Formalization.N6ThreePoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      g
      (swapCoeff (k := k)) := by
  intro E g
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, swapCoeff] using
      (Wedge2Formalization.N6Summary.threePoint_swap_lift_action (k := k)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, swapCoeff, sub_eq_add_neg] using
      (Wedge2Formalization.N6Summary.threePoint_swap_lift_action (k := k)).2

/-- Rank drop on the three-point line occurs exactly at the three distinguished support
points. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • rep₁ (k := k) + b • rep₂ (k := k)) = 0 ↔
      supportCondition (k := k) a b := by
  simpa [supportCondition] using
    Wedge2Formalization.N6Summary.threePoint_det_zero_iff (k := k) (a := a) (b := b)

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h45 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h
  have hb : b = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6ThreePoint.rep₁,
      Wedge2Formalization.N6ThreePoint.rep₂, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h45
  have h01 :=
      congrArg
        (fun M => M (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inl 1)))
        h
  have ha : a = 0 := by
    simpa [hb, rep₁, rep₂, Wedge2Formalization.N6ThreePoint.rep₁,
      Wedge2Formalization.N6ThreePoint.rep₂, Wedge2Formalization.N4.splitRep₁,
      Wedge2Formalization.N4.splitRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J, Matrix.fromBlocks,
      Matrix.add_apply, Matrix.smul_apply] using h01
  exact ⟨ha, hb⟩

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  rcases rep_pair_independent (k := k) hcomb with ⟨h1, _⟩
  exact one_ne_zero h1

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6ThreePoint.V
        Wedge2Formalization.N6ThreePoint.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N6ThreePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6ThreePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1' :
      Wedge2Formalization.N6ThreePoint.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2' :
      Wedge2Formalization.N6ThreePoint.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have h12 :
      Wedge2Formalization.N6ThreePoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g =
        (M 0 0 - M 1 0) • rep₁ (k := k) +
          (M 0 1 - M 1 1) • rep₂ (k := k) := by
    calc
      Wedge2Formalization.N6ThreePoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g
          =
        Wedge2Formalization.N6ThreePoint.ActBivector
          (rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) g := by
            simp [sub_eq_add_neg]
      _ =
        Wedge2Formalization.N6ThreePoint.ActBivector (rep₁ (k := k)) g +
          (-1 : k) • Wedge2Formalization.N6ThreePoint.ActBivector (rep₂ (k := k)) g := by
            rw [Wedge2Formalization.N6ThreePoint.actBivector_add,
              Wedge2Formalization.N6ThreePoint.actBivector_smul]
      _ =
        (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
          (-1 : k) • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
              rw [h1', h2']
      _ =
        (M 0 0 - M 1 0) • rep₁ (k := k) +
          (M 0 1 - M 1 1) • rep₂ (k := k) := by
            ext i j
            simp [sub_eq_add_neg]
            ring
  have hdetM : Matrix.det M ≠ 0 := by
    by_contra hdet
    by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
    · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
      · have hzero :
            Wedge2Formalization.N6ThreePoint.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using hM.1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := rep₁ (k := k)) (g := g) hg).1 hzero
        exact rep₁_ne_zero (k := k) horig
      · have hzero :
            Wedge2Formalization.N6ThreePoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N6ThreePoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) •
                  Wedge2Formalization.N6ThreePoint.ActBivector (rep₁ (k := k)) g +
                M 0 1 •
                  Wedge2Formalization.N6ThreePoint.ActBivector (rep₂ (k := k)) g := by
                    simp [Wedge2Formalization.N6ThreePoint.ActBivector, Matrix.mul_add,
                      Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [hM.1, hM.2]
            _ = 0 := by
                  ext i j
                  simp [hcol0.1, hcol0.2, Matrix.add_apply, Matrix.smul_apply]
                  ring_nf
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb_ne :
            (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k) ≠ 0 := by
          intro hcomb
          rcases rep_pair_independent (k := k) hcomb with ⟨h11, h01⟩
          exact hcol1 ⟨h01, by simpa using neg_eq_zero.mp h11⟩
        exact hcomb_ne horig
    · have hzero :
          Wedge2Formalization.N6ThreePoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
        have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
          have : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
            simpa [Matrix.det_fin_two] using hdet
          exact sub_eq_zero.mp this
        calc
          Wedge2Formalization.N6ThreePoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
              =
            (-(M 1 0)) •
                Wedge2Formalization.N6ThreePoint.ActBivector (rep₁ (k := k)) g +
              M 0 0 •
                Wedge2Formalization.N6ThreePoint.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N6ThreePoint.ActBivector, Matrix.mul_add,
                    Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
          _ =
            (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [hM.1, hM.2]
          _ = 0 := by
                ext i j
                have hdet'' :
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j) =
                      M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by
                  calc
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j)
                        = (M 0 0 * M 1 1) * rep₂ (k := k) i j := by ring
                    _ = (M 0 1 * M 1 0) * rep₂ (k := k) i j := by rw [hdet']
                    _ = M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by ring
                simp [Matrix.add_apply, Matrix.smul_apply]
                rw [hdet'']
                ring_nf
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k)
          (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
          (g := g)
          hg).1 hzero
      have hcomb_ne :
          (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k) ≠ 0 := by
        intro hcomb
        rcases rep_pair_independent (k := k) hcomb with ⟨h10, h00⟩
        exact hcol0 ⟨h00, by simpa using neg_eq_zero.mp h10⟩
      exact hcomb_ne horig
  have hrep1_zero :
      Matrix.det (rep₁ (k := k)) = 0 := by
    simpa [supportCondition] using
      ((det_zero_iff (k := k) (a := (1 : k)) (b := (0 : k))).2 (Or.inr (Or.inr rfl)))
  have hrep2_zero :
      Matrix.det (rep₂ (k := k)) = 0 := by
    simpa [supportCondition] using
      ((det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl))
  have hrep12_zero :
      Matrix.det (rep₁ (k := k) - rep₂ (k := k)) = 0 := by
    have h :
        Matrix.det ((1 : k) • rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) = 0 := by
      simpa [supportCondition, sub_eq_add_neg] using
        ((det_zero_iff (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (Or.inl (by ring))))
    simpa [sub_eq_add_neg] using h
  have hsupp1 :
      supportCondition (k := k) (M 0 0) (M 0 1) := by
    have hdet1 :
        Matrix.det (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) = 0 := by
      rw [← h1', Wedge2Formalization.N6ThreePoint.ActBivector, Matrix.det_mul,
        Matrix.det_mul, Matrix.det_transpose, hrep1_zero]
      simp
    exact (det_zero_iff (k := k) (a := M 0 0) (b := M 0 1)).1 hdet1
  have hsupp2 :
      supportCondition (k := k) (M 1 0) (M 1 1) := by
    have hdet2 :
        Matrix.det (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) = 0 := by
      rw [← h2', Wedge2Formalization.N6ThreePoint.ActBivector, Matrix.det_mul,
        Matrix.det_mul, Matrix.det_transpose, hrep2_zero]
      simp
    exact (det_zero_iff (k := k) (a := M 1 0) (b := M 1 1)).1 hdet2
  have hsupp12 :
      supportCondition (k := k) (M 0 0 - M 1 0) (M 0 1 - M 1 1) := by
    have hdet12 :
        Matrix.det ((M 0 0 - M 1 0) • rep₁ (k := k) + (M 0 1 - M 1 1) • rep₂ (k := k)) = 0 := by
      rw [← h12, Wedge2Formalization.N6ThreePoint.ActBivector, Matrix.det_mul,
        Matrix.det_mul, Matrix.det_transpose, hrep12_zero]
      simp
    exact (det_zero_iff (k := k) (a := M 0 0 - M 1 0) (b := M 0 1 - M 1 1)).1 hdet12
  exact ⟨M, ⟨hdetM, hsupp1, hsupp2, hsupp12⟩, hM⟩

end Row5

end N6
end Paper
end Wedge2Formalization
