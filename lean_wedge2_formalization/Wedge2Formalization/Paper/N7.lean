import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N4
import Wedge2Formalization.Paper.N6
import Wedge2Formalization.Paper.N7.Row1
import Wedge2Formalization.Paper.N7.Row2
import Wedge2Formalization.Paper.N7.Row3
import Wedge2Formalization.Paper.N7.Row5
import Wedge2Formalization.Paper.N7.Row6
import Wedge2Formalization.Paper.N7.Row7
import Wedge2Formalization.Paper.N7.Row8
import Wedge2Formalization.Paper.N7.Row9
import Wedge2Formalization.Paper.N7.Row14
import Wedge2Formalization.N7
import Wedge2Formalization.N7OnePoint
import Wedge2Formalization.N7PureSingular
import Wedge2Formalization.N7ThreePoint
import Wedge2Formalization.N7WeightedTwoPoint
import Wedge2Formalization.N7MixedOnePoint
import Wedge2Formalization.N7Summary
import Mathlib.LinearAlgebra.Matrix.Rank

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

private theorem act_embedded_toBlocks₂₂_N7
    (Ω : Matrix Wedge2Formalization.N7.W Wedge2Formalization.N7.W k)
    (g : Matrix Wedge2Formalization.N7.V Wedge2Formalization.N7.V k) :
    Matrix.toBlocks₂₂
        (Wedge2Formalization.N7.ActBivector (Matrix.fromBlocks 0 0 0 Ω) g) =
      Wedge2Formalization.N4.ActBivector Ω g.toBlocks₂₂ := by
  rw [← Matrix.fromBlocks_toBlocks g]
  simp [Wedge2Formalization.N7.ActBivector, Wedge2Formalization.N4.ActBivector,
    Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

private theorem embeddedAct_fromBlocks
    {I W : Type*}
    [Fintype I] [DecidableEq I] [Fintype W] [DecidableEq W]
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    Matrix.fromBlocks A B C D * Matrix.fromBlocks 0 0 0 Ω *
        (Matrix.fromBlocks A B C D)ᵀ =
      Matrix.fromBlocks
        (B * Ω * Bᵀ)
        (B * Ω * Dᵀ)
        (D * Ω * Bᵀ)
        (D * Ω * Dᵀ) := by
  simp [Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply, Matrix.mul_assoc]

private theorem embeddedAct_fromBlocks_zeroUpperRight
    {I W : Type*}
    [Fintype I] [DecidableEq I] [Fintype W] [DecidableEq W]
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    Matrix.fromBlocks A 0 C D * Matrix.fromBlocks 0 0 0 Ω *
        (Matrix.fromBlocks A 0 C D)ᵀ =
      Matrix.fromBlocks
        0
        0
        0
        (D * Ω * Dᵀ) := by
  simpa using embeddedAct_fromBlocks (k := k) Ω A (0 : Matrix I W k) C D

private theorem upperRight_zero_of_fixing_embedded
    {I W : Type*}
    [Fintype I] [DecidableEq I] [Fintype W] [DecidableEq W]
    (Ω : Matrix W W k)
    (hΩ : Matrix.det Ω ≠ 0)
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hfix :
      Matrix.fromBlocks A B C D * Matrix.fromBlocks 0 0 0 Ω *
          (Matrix.fromBlocks A B C D)ᵀ =
        Matrix.fromBlocks 0 0 0 Ω) :
    B = 0 := by
  have htop :
      B * Ω * Dᵀ = 0 := by
    have h' := congrArg Matrix.toBlocks₁₂ hfix
    simpa [embeddedAct_fromBlocks] using h'
  have hbot :
      D * Ω * Dᵀ = Ω := by
    have h' := congrArg Matrix.toBlocks₂₂ hfix
    simpa [embeddedAct_fromBlocks] using h'
  have hdetD : Matrix.det D ≠ 0 := by
    intro hD
    have hdet := congrArg Matrix.det hbot
    rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose, hD, zero_mul, zero_mul] at hdet
    exact hΩ hdet.symm
  have hΩD :
      Matrix.det (Ω * Dᵀ) ≠ 0 := by
    rw [Matrix.det_mul, Matrix.det_transpose]
    exact mul_ne_zero hΩ hdetD
  have hunit : IsUnit (Matrix.det (Ω * Dᵀ)) :=
    isUnit_iff_ne_zero.mpr hΩD
  calc
    B = B * (Ω * Dᵀ) * (Ω * Dᵀ)⁻¹ := by
          symm
          simpa [Matrix.mul_assoc] using
            (Matrix.mul_nonsing_inv_cancel_right
              (A := Ω * Dᵀ) (B := B) hunit)
    _ = 0 := by
          rw [← Matrix.mul_assoc, htop]
          simp

/-!
This module is the public `Paper/` entry point for `n = 7`.
Together with the split row submodules under `Paper/N7/`, it exposes the
Appendix A row namespaces
`Row1` through `Row16` in the same paper-facing style as `Paper/N4`,
`Paper/N5`, and `Paper/N6`.

Rows `1`, `2`, `3`, `5`, `6`, `7`, `8`, `9`, and `14` are imported from the
split submodules `Paper/N7/Row*.lean`. The remaining public rows are defined in
this file.
-/

/-! Appendix A, `n = 7`, row 10.
Representative `S_1^3 + J_{a,2}`.
Divisor `2[a]`.
Claimed stabilizer:
`K_L = \Ga^{12} \rtimes (\GL_3(k) \times \SL_2(k[t]/(t^2)))`, exact quotient
family `Q_L = B`.
-/
namespace Row10

abbrev I := Wedge2Formalization.N7.I
abbrev W := Wedge2Formalization.N7.W
abbrev V := Wedge2Formalization.N7.V

def rep₁ : Matrix V V k :=
  Wedge2Formalization.N7OnePoint.rad3OnePointRep₁ (k := k)

def rep₂ : Matrix V V k :=
  Wedge2Formalization.N7OnePoint.rad3OnePointRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 10. -/
def paperRep₁ : Matrix V V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 10. -/
def paperRep₂ : Matrix V V k :=
  rep₂ (k := k)

/-- The row-10 paper representative already agrees with the internal working one. -/
def paperChange : Matrix V V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N7.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N7.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N7.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N7.ActBivector]

def U
    (C : Matrix W I k)
    (x y z : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (1 : Matrix I I k)
    0
    C
    (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z)

def Levi
    (A0 : Matrix I I k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix V V k :=
  Matrix.fromBlocks
    A0
    0
    0
    (Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)

/-- Exact pointwise kernel family for the three-dimensional radical extension of the
repeated-support `n = 4` row. -/
def OnePointCore : Set (Matrix W W k) :=
  Wedge2Formalization.Paper.N4.Row1.K (k := k)

/-- Exact pointwise kernel family for the three-dimensional radical extension of the
repeated-support `n = 4` row. -/
def K : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ OnePointCore (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N4.Row1.Qproj (k := k)

def lift
    (a b : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (1 : Matrix I I k)
    0
    0
    (Wedge2Formalization.Paper.N4.Row1.lift (k := k) a b)

private theorem fixes_lowerRight
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hfix :
      Wedge2Formalization.N7OnePoint.FixesRad3OnePointPairBivector
        (Matrix.fromBlocks A B C D)) :
    Wedge2Formalization.N4.FixesOnePointPairBivector (k := k) D :=
  (Wedge2Formalization.N7OnePoint.fixesRad3OnePointPair_fromBlocks_iff
    (k := k) (A := A) (B := B) (C := C) (D := D)).1 hfix |>.2

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N7OnePoint.FixesRad3OnePointPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N7OnePoint.fixesRad3OnePointPair_fromBlocks_iff
        (k := k)
        (A := g.toBlocks₁₁)
        (B := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hg) with
      ⟨hB0, hDfix⟩
    have hDmem :
        g.toBlocks₂₂ ∈ Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
      (Wedge2Formalization.Paper.N4.Row1.pointwise_stabilizer
        (k := k) (g := g.toBlocks₂₂)).1 hDfix
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, g.toBlocks₂₂, hDmem, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := by
            simp [hB0]
  · rintro ⟨A0, C, D, hD, rfl⟩
    have hDfix :
        Wedge2Formalization.N4.FixesOnePointPairBivector (k := k) D :=
      (Wedge2Formalization.Paper.N4.Row1.pointwise_stabilizer
        (k := k) (g := D)).2 hD
    exact
      (Wedge2Formalization.N7OnePoint.fixesRad3OnePointPair_fromBlocks_zeroUpperRight_iff
        (k := k) (A := A0) (C := C) (D := D)).2 hDfix

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ A : Matrix Wedge2Formalization.N4.I
              Wedge2Formalization.N4.I k,
            ∃ x y z : k,
              A.det = 1 ∧
                g =
                  Matrix.fromBlocks
                    A0
                    0
                    C
                    (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
                      Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) := by
  constructor
  · rintro ⟨A0, C, D, hD, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row1.mem_K_iff
        (k := k) (g := D)).1 hD with
      ⟨A, x, y, z, hA, hshape⟩
    refine ⟨A0, C, A, x, y, z, hA, ?_⟩
    simpa [hshape]
  · rintro ⟨A0, C, A, x, y, z, hA, rfl⟩
    refine ⟨A0, C, _, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N4.Row1.mem_K_iff
        (k := k)
        (g :=
          Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
            Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)).2
        ⟨A, x, y, z, hA, rfl⟩

/-- Table-facing kernel statement for Appendix A, row 10. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ A : Matrix Wedge2Formalization.N4.I
              Wedge2Formalization.N4.I k,
            ∃ x y z : k,
              A.det = 1 ∧
                g =
                  Matrix.fromBlocks
                    A0
                    0
                    C
                    (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z *
                      Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) }

/-- Table-facing kernel statement for Appendix A, row 10. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) :=
  mem_K_iff (k := k) (g := g)

theorem U_pointwise
    (C : Matrix W I k)
    (x y z : k) :
    Wedge2Formalization.N7OnePoint.FixesRad3OnePointPairBivector
      (U (k := k) C x y z) := by
  have hDfix :
      Wedge2Formalization.N4.FixesOnePointPairBivector
        (k := k) (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z) :=
    Wedge2Formalization.Paper.N4.Row1.U_pointwise (k := k) x y z
  have hDmem :
      Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z ∈
        Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
    (Wedge2Formalization.Paper.N4.Row1.pointwise_stabilizer
      (k := k)
      (g := Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z)).1 hDfix
  apply (pointwise_stabilizer (k := k) (g := U (k := k) C x y z)).2
  refine ⟨1, C, Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z, hDmem, rfl⟩

theorem Levi_pointwise
    (A0 : Matrix I I k)
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1) :
    Wedge2Formalization.N7OnePoint.FixesRad3OnePointPairBivector
      (Levi (k := k) A0 A) := by
  have hDfix :
      Wedge2Formalization.N4.FixesOnePointPairBivector
        (k := k) (Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A) :=
    Wedge2Formalization.Paper.N4.Row1.Levi_pointwise (k := k) A hA
  have hDmem :
      Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A ∈
        Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
    (Wedge2Formalization.Paper.N4.Row1.pointwise_stabilizer
      (k := k)
      (g := Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)).1 hDfix
  apply (pointwise_stabilizer (k := k) (g := Levi (k := k) A0 A)).2
  refine ⟨A0, 0, Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A, hDmem, ?_⟩
  simp [Levi]

theorem det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂
          (a • rep₁ (k := k) + b • rep₂ (k := k))) = 0 ↔
      a = 0 := by
  simpa [rep₁, rep₂, Wedge2Formalization.N7OnePoint.rad3OnePointRep₁,
    Wedge2Formalization.N7OnePoint.rad3OnePointRep₂] using
    Wedge2Formalization.N4Summary.onePoint_det_zero_iff (k := k) (a := a) (b := b)

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N7.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N7OnePoint.rad3OnePoint_borel_lift_action
        (k := k)
        (A0 := (1 : Matrix I I k))
        (a := a)
        (b := b)
        ha
        (C := 0)).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.N4PaperSummary.row1_borelCoeff] using
      (Wedge2Formalization.N7OnePoint.rad3OnePoint_borel_lift_action
        (k := k)
        (A0 := (1 : Matrix I I k))
        (a := a)
        (b := b)
        ha
        (C := 0)).2

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          b • (Wedge2Formalization.N4.onePointRep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7OnePoint.rad3OnePointRep₁,
      Wedge2Formalization.N7OnePoint.rad3OnePointRep₂] using h'
  have ha : a = 0 := by
    have hdet :
        Matrix.det
          (a • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
            b • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 := by
      have hdet' := congrArg Matrix.det h''
      have hzero :
          Matrix.det
            (0 :
              Matrix Wedge2Formalization.N4.V
                Wedge2Formalization.N4.V k) = 0 := by
        simpa using
          (Matrix.det_zero (n := Wedge2Formalization.N4.V) (R := k))
      exact hdet'.trans hzero
    exact
      Wedge2Formalization.N4Summary.onePoint_action_case_of_rankDrop
        (k := k) (γ := a) (δ := b) hdet
  have hb : b = 0 := by
    have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h''
    simpa [ha, Wedge2Formalization.N4.onePointRep₁,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h01
  exact ⟨ha, hb⟩

theorem quotient_image
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N7.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N7.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have hdetM : Matrix.det M ≠ 0 := by
    by_contra hdet
    by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
    · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
      · have hzero :
            Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using h1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := rep₁ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
          simpa using horig
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).1
      · have hzero :
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) • Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g +
                M 0 1 • Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7.ActBivector, Matrix.mul_add, Matrix.add_mul,
                    Matrix.smul_mul, Matrix.mul_smul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
            _ = 0 := by
                  ext i j
                  simp [hcol0.1, hcol0.2, Matrix.add_apply, Matrix.smul_apply]
                  ring_nf
        have horig :
            (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k) = 0 :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        rcases rep_pair_independent (k := k) horig with ⟨h11, h01⟩
        exact hcol1 ⟨h01, by simpa using neg_eq_zero.mp h11⟩
    · exact by
        have hzero :
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
          have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
            exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hdet)
          calc
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
                =
              (-(M 1 0)) • Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g +
                M 0 0 • Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7.ActBivector, Matrix.mul_add, Matrix.add_mul]
            _ =
              (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
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
        have horig :
            (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k) = 0 :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        rcases rep_pair_independent (k := k) horig with ⟨h10, h00⟩
        exact hcol0 ⟨h00, by simpa using neg_eq_zero.mp h10⟩
  have h2blocks :
      Matrix.toBlocks₂₂
          (Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g) =
        M 1 0 • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          M 1 1 • (Wedge2Formalization.N4.onePointRep₂ (k := k)) := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7OnePoint.rad3OnePointRep₁,
      Wedge2Formalization.N7OnePoint.rad3OnePointRep₂] using
      congrArg Matrix.toBlocks₂₂ h2
  have hrep₂det :
      Matrix.det
        (M 1 0 • (Wedge2Formalization.N4.onePointRep₁ (k := k)) +
          M 1 1 • (Wedge2Formalization.N4.onePointRep₂ (k := k))) = 0 := by
    have hbase :
        Matrix.det (Wedge2Formalization.N4.onePointRep₂ (k := k)) = 0 := by
      simpa [Wedge2Formalization.Paper.N4.Row1.rep₂] using
        ((Wedge2Formalization.Paper.N4.Row1.det_zero_iff
          (k := k) (a := (0 : k)) (b := (1 : k))).2 rfl)
    rw [← h2blocks]
    rw [rep₂, Wedge2Formalization.N7OnePoint.rad3OnePointRep₂, act_embedded_toBlocks₂₂_N7,
      Wedge2Formalization.N4.ActBivector, Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
    simp [hbase]
  have h20 : M 1 0 = 0 := by
    exact
      Wedge2Formalization.N4Summary.onePoint_action_case_of_rankDrop
        (k := k) (γ := M 1 0) (δ := M 1 1) hrep₂det
  refine ⟨M, ?_, hM⟩
  exact ⟨h20, by simpa [Matrix.det_fin_two, h20] using hdetM⟩

end Row10

/-! Appendix A, `n = 7`, row 11.
Representative `S_1^3 + J_{a,1} + J_{b,1}`.
Divisor `[a] + [b]`.
Claimed stabilizer:
`K_L = \Ga^{12} \rtimes (\GL_3(k) \times (\SL_2(k) \times \SL_2(k)))`, exact
quotient family `Q_L = N_T`.
-/
namespace Row11

abbrev I := Wedge2Formalization.N7.I
abbrev W := Wedge2Formalization.N7.W
abbrev V := Wedge2Formalization.N7.V

def rep₁ : Matrix V V k :=
  Wedge2Formalization.N7.rad3SplitRep₁ (k := k)

def rep₂ : Matrix V V k :=
  Wedge2Formalization.N7.rad3SplitRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 11. -/
def paperRep₁ : Matrix V V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 11. -/
def paperRep₂ : Matrix V V k :=
  rep₂ (k := k)

/-- The row-11 paper representative already agrees with the internal working one. -/
def paperChange : Matrix V V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N7.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N7.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N7.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N7.ActBivector]

def Levi
    (A0 : Matrix I I k)
    (A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k) :
    Matrix V V k :=
  Matrix.fromBlocks
    A0
    0
    0
    (Matrix.fromBlocks A 0 0 D)

def SplitPairCore : Set (Matrix W W k) :=
  Wedge2Formalization.Paper.N4.Row2.K (k := k)

def K : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D0 : Matrix W W k,
            D0 ∈ SplitPairCore (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D0 }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N4.Row2.Qproj (k := k)

def torusLift
    (u v : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (1 : Matrix I I k)
    0
    0
    (Wedge2Formalization.Paper.N4.Row2.torusLift (k := k) u v)

def swapLift
    (u v : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (1 : Matrix I I k)
    0
    0
    (Wedge2Formalization.Paper.N4.Row2.swapLift (k := k) u v)

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N7.FixesRad3SplitPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N7.fixesRad3SplitPair_fromBlocks_iff
        (k := k)
        (A := g.toBlocks₁₁)
        (B := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hg) with
      ⟨hB0, hDfix⟩
    have hDmem :
        g.toBlocks₂₂ ∈ Wedge2Formalization.Paper.N4.Row2.K (k := k) :=
      (Wedge2Formalization.Paper.N4.Row2.pointwise_stabilizer
        (k := k) (g := g.toBlocks₂₂)).1 hDfix
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, g.toBlocks₂₂, hDmem, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := by
            simp [hB0]
  · rintro ⟨A0, C, D0, hD0, rfl⟩
    have hDfix :
        Wedge2Formalization.N4.FixesSplitPairBivector (k := k) D0 :=
      (Wedge2Formalization.Paper.N4.Row2.pointwise_stabilizer
        (k := k) (g := D0)).2 hD0
    exact
      (Wedge2Formalization.N7.fixesRad3SplitPair_fromBlocks_zeroUpperRight_iff
        (k := k) (A := A0) (C := C) (D := D0)).2 hDfix

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
              D.det = 1 ∧
              g = Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D) := by
  constructor
  · rintro ⟨A0, C, D0, hD0, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N4.Row2.mem_K_iff
        (k := k) (g := D0)).1 hD0 with
      ⟨A, D, hA, hD, hshape⟩
    refine ⟨A0, C, A, D, hA, hD, ?_⟩
    simpa [Wedge2Formalization.Paper.N4.Row2.Levi] using congrArg (Matrix.fromBlocks A0 0 C) hshape
  · rintro ⟨A0, C, A, D, hA, hD, rfl⟩
    refine ⟨A0, C, Matrix.fromBlocks A 0 0 D, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N4.Row2.mem_K_iff
        (k := k) (g := Matrix.fromBlocks A 0 0 D)).2 ⟨A, D, hA, hD, rfl⟩

/-- Table-facing kernel statement for Appendix A, row 11. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
            A.det = 1 ∧
              D.det = 1 ∧
              g = Matrix.fromBlocks A0 0 C (Matrix.fromBlocks A 0 0 D) }

/-- Table-facing kernel statement for Appendix A, row 11. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) :=
  mem_K_iff (k := k) (g := g)

theorem Levi_pointwise
    (A0 : Matrix I I k)
    (A D : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1)
    (hD : D.det = 1) :
    Wedge2Formalization.N7.FixesRad3SplitPairBivector
      (Levi (k := k) A0 A D) := by
  apply (pointwise_stabilizer (k := k) (g := Levi (k := k) A0 A D)).2
  refine ⟨A0, 0, Matrix.fromBlocks A 0 0 D, ?_, ?_⟩
  · exact
      (Wedge2Formalization.Paper.N4.Row2.mem_K_iff
        (k := k) (g := Matrix.fromBlocks A 0 0 D)).2 ⟨A, D, hA, hD, rfl⟩
  · simp [Levi]

theorem det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂
          (a • rep₁ (k := k) + b • rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 := by
  simpa [rep₁, rep₂, Wedge2Formalization.N7.rad3SplitRep₁,
    Wedge2Formalization.N7.rad3SplitRep₂] using
    Wedge2Formalization.N4Summary.split_det_zero_iff (k := k) (a := a) (b := b)

theorem torus_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N7.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (torusLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_torusCoeff (k := k) u v) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, torusLift] using
      (Wedge2Formalization.N7.rad3Split_torus_lift_action
        (k := k)
        (A0 := (1 : Matrix I I k))
        (u := u)
        (v := v)
        (C := 0)).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, torusLift,
      Wedge2Formalization.N4PaperSummary.row2_torusCoeff] using
      (Wedge2Formalization.N7.rad3Split_torus_lift_action
        (k := k)
        (A0 := (1 : Matrix I I k))
        (u := u)
        (v := v)
        (C := 0)).2

theorem swap_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N7.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (swapLift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_swapCoeff (k := k) u v) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, swapLift] using
      (Wedge2Formalization.N7.rad3Split_swapCoset_lift_action
        (k := k)
        (A0 := (1 : Matrix I I k))
        (u := u)
        (v := v)
        (C := 0)).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, swapLift,
      Wedge2Formalization.N4PaperSummary.row2_swapCoeff] using
      (Wedge2Formalization.N7.rad3Split_swapCoset_lift_action
        (k := k)
        (A0 := (1 : Matrix I I k))
        (u := u)
        (v := v)
        (C := 0)).2

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          b • (Wedge2Formalization.N4.splitRep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7.rad3SplitRep₁,
      Wedge2Formalization.N7.rad3SplitRep₂] using h'
  have ha : a = 0 := by
    have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h''
    simpa [Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h01
  have hb : b = 0 := by
    have h45 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h''
    simpa [ha, Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.splitRep₂,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h45
  exact ⟨ha, hb⟩

theorem quotient_image
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N7.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N7.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have h12 :
      Wedge2Formalization.N7.ActBivector (rep₁ (k := k) - rep₂ (k := k)) g =
        (M 0 0 - M 1 0) • rep₁ (k := k) + (M 0 1 - M 1 1) • rep₂ (k := k) := by
    calc
      Wedge2Formalization.N7.ActBivector (rep₁ (k := k) - rep₂ (k := k)) g
          =
        Wedge2Formalization.N7.ActBivector (rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) g := by
            simp [sub_eq_add_neg]
      _ =
        Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g +
          (-1 : k) • Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g := by
            simp [Wedge2Formalization.N7.ActBivector, Matrix.mul_add, Matrix.add_mul,
              Matrix.smul_mul, Matrix.mul_smul]
      _ =
        (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
          (-1 : k) • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
              rw [h1, h2]
      _ = (M 0 0 - M 1 0) • rep₁ (k := k) + (M 0 1 - M 1 1) • rep₂ (k := k) := by
            ext i j
            simp [Matrix.add_apply, Matrix.smul_apply, sub_eq_add_neg]
            ring
  have hdetM : Matrix.det M ≠ 0 := by
    by_contra hdet
    by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
    · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
      · have hzero :
            Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using h1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := rep₁ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
          simpa using horig
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).1
      · have hzero :
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) • Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g +
                M 0 1 • Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7.ActBivector, Matrix.mul_add, Matrix.add_mul,
                    Matrix.smul_mul, Matrix.mul_smul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
            _ = 0 := by
                  ext i j
                  simp [hcol0.1, hcol0.2, Matrix.add_apply, Matrix.smul_apply]
                  ring_nf
        have horig :
            (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k) = 0 :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        rcases rep_pair_independent (k := k) horig with ⟨h11, h01⟩
        exact hcol1 ⟨h01, by simpa using neg_eq_zero.mp h11⟩
    · exact by
        have hzero :
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
          have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
            exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hdet)
          calc
            Wedge2Formalization.N7.ActBivector
                ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
                =
              (-(M 1 0)) • Wedge2Formalization.N7.ActBivector (rep₁ (k := k)) g +
                M 0 0 • Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7.ActBivector, Matrix.mul_add, Matrix.add_mul]
            _ =
              (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
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
        have horig :
            (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k) = 0 :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        rcases rep_pair_independent (k := k) horig with ⟨h10, h00⟩
        exact hcol0 ⟨h00, by simpa using neg_eq_zero.mp h10⟩
  have h2blocks :
      Matrix.toBlocks₂₂
          (Wedge2Formalization.N7.ActBivector (rep₂ (k := k)) g) =
        M 1 0 • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          M 1 1 • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7.rad3SplitRep₁,
      Wedge2Formalization.N7.rad3SplitRep₂] using
      congrArg Matrix.toBlocks₂₂ h2
  have h12blocks :
      Matrix.toBlocks₂₂
          (Wedge2Formalization.N7.ActBivector (rep₁ (k := k) - rep₂ (k := k)) g) =
        (M 0 0 - M 1 0) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          (M 0 1 - M 1 1) • (Wedge2Formalization.N4.splitRep₂ (k := k)) := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7.rad3SplitRep₁,
      Wedge2Formalization.N7.rad3SplitRep₂] using
      congrArg Matrix.toBlocks₂₂ h12
  have hrep2_zero :
      Matrix.det
        (M 1 0 • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          M 1 1 • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    have hbase :
        Matrix.det (Wedge2Formalization.N4.splitRep₂ (k := k)) = 0 := by
      simpa [Wedge2Formalization.Paper.N4.Row2.rep₂] using
        ((Wedge2Formalization.Paper.N4.Row2.det_zero_iff
          (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl))
    rw [← h2blocks]
    rw [rep₂, Wedge2Formalization.N7.rad3SplitRep₂, act_embedded_toBlocks₂₂_N7,
      Wedge2Formalization.N4.ActBivector, Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
    simp [hbase]
  have hrep12_zero :
      Matrix.det
        ((M 0 0 - M 1 0) • (Wedge2Formalization.N4.splitRep₁ (k := k)) +
          (M 0 1 - M 1 1) • (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
    have hbase :
        Matrix.det
          ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
            (Wedge2Formalization.N4.splitRep₂ (k := k))) = 0 := by
      simpa [sub_eq_add_neg] using
        ((Wedge2Formalization.Paper.N4.Row2.det_zero_iff
          (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (by ring)))
    rw [← h12blocks]
    have hsub :
        rep₁ (k := k) - rep₂ (k := k) =
          Matrix.fromBlocks 0 0 0
            ((Wedge2Formalization.N4.splitRep₁ (k := k)) -
              (Wedge2Formalization.N4.splitRep₂ (k := k))) := by
      ext i j
      cases i <;> cases j <;>
        simp [rep₁, rep₂, Wedge2Formalization.N7.rad3SplitRep₁,
          Wedge2Formalization.N7.rad3SplitRep₂, sub_eq_add_neg,
          Matrix.fromBlocks, Matrix.add_apply]
    rw [hsub, act_embedded_toBlocks₂₂_N7, Wedge2Formalization.N4.ActBivector,
      Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
    simp [hbase]
  have hcase :=
    Wedge2Formalization.N4Summary.split_action_cases_of_rankDrop
      (k := k)
      (α := M 0 0)
      (β := M 0 1)
      (γ := M 1 0)
      (δ := M 1 1)
      (hinv := by simpa [Matrix.det_fin_two] using hdetM)
      (hrep2 := hrep2_zero)
      (hω12 := hrep12_zero)
  refine ⟨M, ?_, hM⟩
  rcases hcase with ⟨hγ0, hrel⟩ | ⟨hαγ, hsum⟩
  · left
    refine ⟨M 0 0, M 1 1, ?_⟩
    have h01 : M 0 1 = M 1 1 - M 0 0 := by
      calc
        M 0 1 = M 0 0 + M 0 1 - M 0 0 := by ring
        _ = M 1 1 - M 0 0 := by rw [hrel]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N4PaperSummary.row2_torusCoeff, hγ0, h01]
  · right
    refine ⟨M 0 0, M 0 1 + M 0 0, ?_⟩
    have h11 : M 1 1 = -M 1 0 := by
      exact eq_neg_of_add_eq_zero_left (by simpa [add_comm] using hsum)
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.N4PaperSummary.row2_swapCoeff, hαγ, h11]

end Row11

/-! Appendix A, `n = 7`, row 12.
Representative `S_1 + J_{a,3}`.
Divisor `3[a]`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times \SL_2(k[t]/(t^3)))`, exact quotient
family `Q_L = B`.
-/
namespace Row12

abbrev I := Fin 1
abbrev W := Wedge2Formalization.N6OnePointLong.V
abbrev V := Sum I W

def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

def scalarBlock (a : k) : Matrix I I k :=
  !![a]

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row1.paperRep₁ (k := k))

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row1.paperRep₂ (k := k))

def FixesPairBivector (g : Matrix V V k) : Prop :=
  ActBivector (rep₁ (k := k)) g = rep₁ (k := k) ∧
    ActBivector (rep₂ (k := k)) g = rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 12. -/
def paperRep₁ : Matrix V V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 12. -/
def paperRep₂ : Matrix V V k :=
  rep₂ (k := k)

/-- The row-12 paper representative already agrees with the internal working pair. -/
def paperChange : Matrix V V k :=
  1

theorem paperRep₁_transport :
    ActBivector (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [ActBivector, paperRep₁, paperChange, rep₁]

theorem paperRep₂_transport :
    ActBivector (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [ActBivector, paperRep₂, paperChange, rep₂]

def U
    (C : Matrix W I k)
    (x₁ x₂ y₁ y₂ z₁ z₂ : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (1 : Matrix I I k)
    0
    C
    (Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂)

def Levi
    (a0 a b c d : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (scalarBlock (k := k) a0)
    0
    0
    (Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d)

/-- Public name for the embedded Appendix A, `n = 6`, row 1 core
`U_6 \rtimes \SL_2(k)` on the lower-right `6`-space. -/
def OnePointLongCore : Set (Matrix W W k) :=
  Wedge2Formalization.Paper.N6.Row1.K (k := k)

def K : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ OnePointLongCore (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N6.Row1.Qproj (k := k)

def lift
    (a b : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (1 : Matrix I I k)
    0
    0
    (Wedge2Formalization.Paper.N6.Row1.lift (k := k) a b)

private theorem rep₁_det_ne_zero :
    Matrix.det (Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k)) ≠ 0 := by
  intro hdet
  have hzero :
      Matrix.det ((1 : k) • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
        (0 : k) • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) = 0 := by
    simpa using hdet
  have : (1 : k) = 0 :=
    (Wedge2Formalization.Paper.N6.Row1.det_zero_iff
      (k := k) (a := (1 : k)) (b := (0 : k))).1 hzero
  exact one_ne_zero this

private theorem act_embedded_toBlocks₂₂
    (Ω : Matrix W W k)
    (g : Matrix V V k) :
    Matrix.toBlocks₂₂
        (ActBivector (Matrix.fromBlocks 0 0 0 Ω) g) =
      g.toBlocks₂₂ * Ω * g.toBlocks₂₂ᵀ := by
  rw [← Matrix.fromBlocks_toBlocks g]
  simp [ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

private theorem embed_linearCombination
    (a b : k) :
    Matrix.fromBlocks
        (0 : Matrix I I k)
        (0 : Matrix I W k)
        (0 : Matrix W I k)
        (a • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
          b • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) =
      a • rep₁ (k := k) + b • rep₂ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j
  · fin_cases i
    fin_cases j
    simp [rep₁, rep₂, Wedge2Formalization.Paper.N6.Row1.rep₁,
      Wedge2Formalization.Paper.N6.Row1.rep₂, Wedge2Formalization.Paper.N6.Row1.paperRep₁,
      Wedge2Formalization.Paper.N6.Row1.paperRep₂, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Matrix.fromBlocks, Matrix.add_apply]
  · simp [rep₁, rep₂, Wedge2Formalization.Paper.N6.Row1.rep₁,
      Wedge2Formalization.Paper.N6.Row1.rep₂, Wedge2Formalization.Paper.N6.Row1.paperRep₁,
      Wedge2Formalization.Paper.N6.Row1.paperRep₂, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Matrix.fromBlocks, Matrix.add_apply]
  · simp [rep₁, rep₂, Wedge2Formalization.Paper.N6.Row1.rep₁,
      Wedge2Formalization.Paper.N6.Row1.rep₂, Wedge2Formalization.Paper.N6.Row1.paperRep₁,
      Wedge2Formalization.Paper.N6.Row1.paperRep₂, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Matrix.fromBlocks, Matrix.add_apply]
  · simp [rep₁, rep₂, Wedge2Formalization.Paper.N6.Row1.rep₁,
      Wedge2Formalization.Paper.N6.Row1.rep₂, Wedge2Formalization.Paper.N6.Row1.paperRep₁,
      Wedge2Formalization.Paper.N6.Row1.paperRep₂, Wedge2Formalization.N6OnePointLong.rep₁,
      Wedge2Formalization.N6OnePointLong.rep₂, Matrix.fromBlocks, Matrix.add_apply]

theorem pointwise_stabilizer :
    ExactFamily
      (FixesPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · rintro ⟨h1, h2⟩
    have h1' :
        Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ *
            Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k)) *
            (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
          Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k)) := by
      simpa [ActBivector, rep₁, Matrix.fromBlocks_toBlocks] using h1
    have hB0 :
        g.toBlocks₁₂ = 0 := by
      exact
        upperRight_zero_of_fixing_embedded
          (k := k)
          (Ω := Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k))
          (hΩ := rep₁_det_ne_zero (k := k))
          (A := g.toBlocks₁₁)
          (B := g.toBlocks₁₂)
          (C := g.toBlocks₂₁)
          (D := g.toBlocks₂₂)
          h1'
    have hDfix1 :
        Wedge2Formalization.N6OnePointLong.ActBivector
            (Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k))
            g.toBlocks₂₂ =
          Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) := by
      have h' := congrArg Matrix.toBlocks₂₂ h1
      have h'' := congrArg Matrix.toBlocks₂₂ h1
      simpa [rep₁, hB0, act_embedded_toBlocks₂₂,
        Wedge2Formalization.N6OnePointLong.ActBivector] using h''
    have hDfix2 :
        Wedge2Formalization.N6OnePointLong.ActBivector
            (Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k))
            g.toBlocks₂₂ =
          Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k) := by
      have h'' := congrArg Matrix.toBlocks₂₂ h2
      simpa [rep₂, hB0, act_embedded_toBlocks₂₂,
        Wedge2Formalization.N6OnePointLong.ActBivector] using h''
    have hDmem :
        g.toBlocks₂₂ ∈ Wedge2Formalization.Paper.N6.Row1.K (k := k) :=
      (Wedge2Formalization.Paper.N6.Row1.pointwise_stabilizer
        (k := k) (g := g.toBlocks₂₂)).1 ⟨hDfix1, hDfix2⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, g.toBlocks₂₂, hDmem, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := by
            simp [hB0]
  · rintro ⟨A0, C, D, hD, rfl⟩
    have hDfix :
        Wedge2Formalization.N6OnePointLong.FixesPairBivector (k := k) D :=
      (Wedge2Formalization.Paper.N6.Row1.pointwise_stabilizer
        (k := k) (g := D)).2 hD
    constructor
    · calc
        ActBivector (rep₁ (k := k)) (Matrix.fromBlocks A0 0 C D)
            =
          Matrix.fromBlocks 0 0 0
            (Wedge2Formalization.N6OnePointLong.ActBivector
              (Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k)) D) := by
                simpa [ActBivector, rep₁] using
                  (embeddedAct_fromBlocks_zeroUpperRight
                    (k := k)
                    (Ω := Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k))
                    (A := A0)
                    (C := C)
                    (D := D))
        _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k)) := by
              simpa [Wedge2Formalization.N6OnePointLong.FixesPairBivector] using hDfix.1
        _ = rep₁ (k := k) := rfl
    · calc
        ActBivector (rep₂ (k := k)) (Matrix.fromBlocks A0 0 C D)
            =
          Matrix.fromBlocks 0 0 0
            (Wedge2Formalization.N6OnePointLong.ActBivector
              (Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) D) := by
                simpa [ActBivector, rep₂] using
                  (embeddedAct_fromBlocks_zeroUpperRight
                    (k := k)
                    (Ω := Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k))
                    (A := A0)
                    (C := C)
                    (D := D))
        _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) := by
              simpa [Wedge2Formalization.N6OnePointLong.FixesPairBivector] using hDfix.2
        _ = rep₂ (k := k) := rfl

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ Wedge2Formalization.Paper.N6.Row1.K (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D := by
  rfl

theorem mem_K_shape_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ a0 x₁ x₂ y₁ y₂ z₁ z₂ a b c d : k,
        ∃ C : Matrix W I k,
          a * d - b * c = 1 ∧
          g =
            Matrix.fromBlocks
              (scalarBlock (k := k) a0)
              0
              C
              (Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ *
                Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d) := by
  constructor
  · rintro ⟨A0, C, D, hD, rfl⟩
    rcases
      (Wedge2Formalization.Paper.N6.Row1.mem_K_table_iff
        (k := k) (g := D)).1 hD with
      ⟨a, b, c, d, x₁, x₂, y₁, y₂, z₁, z₂, hdet, hshape⟩
    refine ⟨A0 0 0, x₁, x₂, y₁, y₂, z₁, z₂, a, b, c, d, C, hdet, ?_⟩
    rw [hshape]
    ext i j
    rcases i with i | i <;> rcases j with j | j
    · fin_cases i
      fin_cases j
      simp [scalarBlock, Matrix.fromBlocks]
    · simp [scalarBlock, Matrix.fromBlocks]
    · simp [scalarBlock, Matrix.fromBlocks]
    · simp [scalarBlock, Matrix.fromBlocks]
  · rintro ⟨a0, x₁, x₂, y₁, y₂, z₁, z₂, a, b, c, d, C, hdet, rfl⟩
    refine ⟨scalarBlock (k := k) a0, C,
      Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ *
        Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d, ?_, rfl⟩
    exact
      (Wedge2Formalization.Paper.N6.Row1.mem_K_table_iff
        (k := k)
        (g := Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ *
          Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d)).2
        ⟨a, b, c, d, x₁, x₂, y₁, y₂, z₁, z₂, hdet, rfl⟩

/-- Table-facing kernel statement for Appendix A, row 12, using the named embedded
`U_6 \rtimes \SL_2(k)` core on the lower-right block. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ a0 x₁ x₂ y₁ y₂ z₁ z₂ a b c d : k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ OnePointLongCore (k := k) ∧
            a * d - b * c = 1 ∧
            D =
              Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ *
                Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d ∧
            g = Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C D }

/-- Table-facing kernel statement for Appendix A, row 12, using the named embedded
`U_6 \rtimes \SL_2(k)` core on the lower-right block. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) := by
  constructor
  · rintro h
    rcases (mem_K_shape_iff (k := k) (g := g)).1 h with
      ⟨a0, x₁, x₂, y₁, y₂, z₁, z₂, a, b, c, d, C, hdet, hshape⟩
    refine ⟨a0, x₁, x₂, y₁, y₂, z₁, z₂, a, b, c, d, C,
      Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ *
        Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d,
      ?_, hdet, rfl, hshape⟩
    exact
      (Wedge2Formalization.Paper.N6.Row1.mem_K_table_iff
        (k := k)
        (g := Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ *
          Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d)).2
        ⟨a, b, c, d, x₁, x₂, y₁, y₂, z₁, z₂, hdet, rfl⟩
  · rintro ⟨a0, x₁, x₂, y₁, y₂, z₁, z₂, a, b, c, d, C, D, hD, hdet, hshapeD, hshape⟩
    have hshape' :
        g =
          Matrix.fromBlocks
            (scalarBlock (k := k) a0)
            0
            C
            (Wedge2Formalization.Paper.N6.Row1.U (k := k) x₁ x₂ y₁ y₂ z₁ z₂ *
              Wedge2Formalization.Paper.N6.Row1.Levi (k := k) a b c d) := by
      simpa [hshapeD] using hshape
    exact (mem_K_shape_iff (k := k) (g := g)).2
      ⟨a0, x₁, x₂, y₁, y₂, z₁, z₂, a, b, c, d, C, hdet, hshape'⟩

theorem det_zero_iff
    (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • rep₁ (k := k) + b • rep₂ (k := k))) = 0 ↔
      a = 0 := by
  simpa [rep₁, rep₂] using
    Wedge2Formalization.Paper.N6.Row1.det_zero_iff (k := k) (a := a) (b := b)

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (!![a, b; 0, 1] : Matrix (Fin 2) (Fin 2) k) := by
  constructor
  · calc
      ActBivector (rep₁ (k := k)) (lift (k := k) a b)
          =
        Matrix.fromBlocks 0 0 0
          (Wedge2Formalization.N6OnePointLong.ActBivector
            (Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k))
            (Wedge2Formalization.Paper.N6.Row1.lift (k := k) a b)) := by
              simpa [ActBivector, rep₁, lift] using
                (embeddedAct_fromBlocks_zeroUpperRight
                  (k := k)
                  (Ω := Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k))
                  (A := (1 : Matrix I I k))
                  (C := 0)
                  (D := Wedge2Formalization.Paper.N6.Row1.lift (k := k) a b))
      _ =
        Matrix.fromBlocks 0 0 0
          (a • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
            b • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) := by
              simpa [ActsOnOrderedPair] using
                (Wedge2Formalization.Paper.N6.Row1.quotient_action
                  (k := k) (a := a) (b := b) ha).1
      _ = a • rep₁ (k := k) + b • rep₂ (k := k) := by
            simpa using embed_linearCombination (k := k) (a := a) (b := b)
  · calc
      ActBivector (rep₂ (k := k)) (lift (k := k) a b)
          =
        Matrix.fromBlocks 0 0 0
          (Wedge2Formalization.N6OnePointLong.ActBivector
            (Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k))
            (Wedge2Formalization.Paper.N6.Row1.lift (k := k) a b)) := by
              simpa [ActBivector, rep₂, lift] using
                (embeddedAct_fromBlocks_zeroUpperRight
                  (k := k)
                  (Ω := Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k))
                  (A := (1 : Matrix I I k))
                  (C := 0)
                  (D := Wedge2Formalization.Paper.N6.Row1.lift (k := k) a b))
      _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) := by
            simpa [ActsOnOrderedPair] using
              (Wedge2Formalization.Paper.N6.Row1.quotient_action
                (k := k) (a := a) (b := b) ha).2
      _ = (0 : k) • rep₁ (k := k) + (1 : k) • rep₂ (k := k) := by
            simpa using (embed_linearCombination (k := k) (a := (0 : k)) (b := (1 : k))).symm

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
          b • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k) = 0 := by
    simpa [rep₁, rep₂] using h'
  have ha : a = 0 := by
    have hdet :
        Matrix.det
          (a • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
            b • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) = 0 := by
      have hdet' :
          Matrix.det
              (a • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
                b • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) =
            Matrix.det (0 : Matrix W W k) :=
        congrArg Matrix.det h''
      have hdet0 : Matrix.det (0 : Matrix W W k) = 0 := by
        let _ : Nonempty W := inferInstance
        simpa using (Matrix.det_zero (n := W) (R := k))
      exact Eq.trans hdet' hdet0
    exact
      (Wedge2Formalization.Paper.N6.Row1.det_zero_iff
        (k := k) (a := a) (b := b)).1 hdet
  have hb : b = 0 := by
    have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inr 1)) h''
    simpa [ha, Wedge2Formalization.Paper.N6.Row1.rep₁, Wedge2Formalization.Paper.N6.Row1.rep₂,
      Wedge2Formalization.N6OnePointLong.rep₁, Wedge2Formalization.N6OnePointLong.rep₂,
      Wedge2Formalization.N6OnePointLong.N, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h01
  exact ⟨ha, hb⟩

theorem quotient_image
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have hdetM : Matrix.det M ≠ 0 := by
    by_contra hdet
    by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
    · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
      · have hzero : ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [hcol0.1, hcol1.1] using h1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := rep₁ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
          simpa using horig
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).1
      · have hzero :
            ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k))
                g = 0 := by
          calc
            ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k))
                g
                =
              (-(M 1 1)) • ActBivector (rep₁ (k := k)) g +
                M 0 1 • ActBivector (rep₂ (k := k)) g := by
                  simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
            _ = 0 := by
                  ext i j
                  simp [hcol0.1, hcol0.2, Matrix.add_apply, Matrix.smul_apply]
                  ring
        have horig :
            (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k) = 0 :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        rcases rep_pair_independent (k := k) horig with ⟨h11, h01⟩
        exact hcol1 ⟨h01, neg_eq_zero.mp h11⟩
    · have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
        exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hdet)
      have hzero :
          ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
              g = 0 := by
        calc
          ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
              g
              =
            (-(M 1 0)) • ActBivector (rep₁ (k := k)) g +
              M 0 0 • ActBivector (rep₂ (k := k)) g := by
                simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
          _ =
            (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [h1, h2]
          _ = 0 := by
                ext i j
                simp [Matrix.add_apply, Matrix.smul_apply]
                have hcoeff : M 0 0 * M 1 1 = M 1 0 * M 0 1 := by
                  simpa [mul_comm] using hdet'
                calc
                  -(M 1 0 * (M 0 0 * rep₁ (k := k) i j)) +
                        -(M 1 0 * (M 0 1 * rep₂ (k := k) i j)) +
                        (M 0 0 * (M 1 0 * rep₁ (k := k) i j) +
                          M 0 0 * (M 1 1 * rep₂ (k := k) i j))
                        = (M 0 0 * M 1 1 - M 1 0 * M 0 1) * rep₂ (k := k) i j := by
                            ring
                    _ = 0 := by rw [hcoeff]; ring
      have horig :
          (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k) = 0 :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k)
          (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
          (g := g)
          hg).1 hzero
      rcases rep_pair_independent (k := k) horig with ⟨h10, h00⟩
      exact hcol0 ⟨h00, neg_eq_zero.mp h10⟩
  have h2blocks :
      Matrix.toBlocks₂₂ (ActBivector (rep₂ (k := k)) g) =
        M 1 0 • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
          M 1 1 • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k) := by
    simpa [rep₁, rep₂] using congrArg Matrix.toBlocks₂₂ h2
  have hact₂det :
      Matrix.det (Matrix.toBlocks₂₂ (ActBivector (rep₂ (k := k)) g)) = 0 := by
    have hbase :
        Matrix.det (Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) = 0 := by
      simpa [Wedge2Formalization.Paper.N6.Row1.rep₂] using
        ((Wedge2Formalization.Paper.N6.Row1.det_zero_iff
          (k := k) (a := (0 : k)) (b := (1 : k))).2 rfl)
    have hpaperbase :
        Matrix.det (Wedge2Formalization.Paper.N6.Row1.paperRep₂ (k := k)) = 0 := by
      simpa [Wedge2Formalization.Paper.N6.Row1.paperRep₂] using hbase
    rw [rep₂, act_embedded_toBlocks₂₂]
    calc
      Matrix.det
          (g.toBlocks₂₂ *
            Wedge2Formalization.Paper.N6.Row1.paperRep₂ (k := k) *
            g.toBlocks₂₂ᵀ) =
          Matrix.det g.toBlocks₂₂ *
          Matrix.det (Wedge2Formalization.Paper.N6.Row1.paperRep₂ (k := k)) *
          Matrix.det g.toBlocks₂₂ᵀ := by
            rw [Matrix.det_mul, Matrix.det_mul]
      _ = 0 := by
            rw [hpaperbase]
            ring
  have hrep₂det :
      Matrix.det
        (M 1 0 • Wedge2Formalization.Paper.N6.Row1.rep₁ (k := k) +
          M 1 1 • Wedge2Formalization.Paper.N6.Row1.rep₂ (k := k)) = 0 := by
    rwa [h2blocks] at hact₂det
  have h20 : M 1 0 = 0 := by
    exact
      (Wedge2Formalization.Paper.N6.Row1.det_zero_iff
        (k := k) (a := M 1 0) (b := M 1 1)).1 hrep₂det
  refine ⟨M, ?_, hM⟩
  simpa [Qproj, Wedge2Formalization.Paper.N6.Row1.Qproj] using
    (show M 1 0 = 0 ∧ Matrix.det M ≠ 0 from ⟨h20, hdetM⟩)

end Row12

/-! Appendix A, `n = 7`, row 13.
Representative `S_1 + J_{a,2} + J_{a,1}`.
Divisor `3[a]`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times (U_{2,1}(k) \rtimes
(\SL_2(k) \times \SL_2(k))))`, exact quotient family `Q_L = B`.
-/
namespace Row13

abbrev I := Wedge2Formalization.N7MixedOnePoint.I
abbrev W := Wedge2Formalization.N7MixedOnePoint.W
abbrev V := Wedge2Formalization.N7MixedOnePoint.V

def rep₁ : Matrix V V k :=
  Wedge2Formalization.N7MixedOnePoint.radMixedRep₁ (k := k)

def rep₂ : Matrix V V k :=
  Wedge2Formalization.N7MixedOnePoint.radMixedRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 13:
`e₂∧e₄ + e₃∧e₅ + e₆∧e₇`. -/
def paperRep₁ : Matrix V V k :=
  Matrix.fromBlocks
    0
    0
    0
    (Wedge2Formalization.Paper.N6.Row2.paperRep₁ (k := k))

/-- Literal second basis vector from Appendix A, row 13:
`e₂∧e₅`. -/
def paperRep₂ : Matrix V V k :=
  Matrix.fromBlocks
    0
    0
    0
    (Wedge2Formalization.Paper.N6.Row2.paperRep₂ (k := k))

/-- Explicit basis change sending the literal paper representative to the internal
working representative pair. -/
def paperChange : Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k))
    0
    0
    (Wedge2Formalization.Paper.N6.Row2.paperChange (k := k))

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N7MixedOnePoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  calc
    Wedge2Formalization.N7MixedOnePoint.ActBivector
        (paperRep₁ (k := k)) (paperChange (k := k)) =
      Matrix.fromBlocks
        0
        0
        0
        (Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.Paper.N6.Row2.paperRep₁ (k := k))
          (Wedge2Formalization.Paper.N6.Row2.paperChange (k := k))) := by
            simpa [paperRep₁, paperChange] using
              (Wedge2Formalization.N7MixedOnePoint.act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := Wedge2Formalization.Paper.N6.Row2.paperRep₁ (k := k))
                (a := Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k))
                (C := 0)
                (D := Wedge2Formalization.Paper.N6.Row2.paperChange (k := k)))
    _ =
      Matrix.fromBlocks
        0
        0
        0
        (Wedge2Formalization.Paper.N6.Row2.rep₁ (k := k)) := by
          simpa [Wedge2Formalization.Paper.N6.Row2.rep₁] using
            congrArg
              (fun Ω :
                Matrix Wedge2Formalization.N6MixedOnePoint.V
                  Wedge2Formalization.N6MixedOnePoint.V k =>
                Matrix.fromBlocks
                  (0 : Matrix I I k)
                  (0 : Matrix I W k)
                  (0 : Matrix W I k)
                  Ω)
              (Wedge2Formalization.Paper.N6.Row2.paperRep₁_transport (k := k))
    _ = rep₁ (k := k) := by
          simpa [rep₁, Wedge2Formalization.N7MixedOnePoint.radMixedRep₁,
            Wedge2Formalization.Paper.N6.Row2.rep₁,
            Wedge2Formalization.N6MixedOnePoint.rep₁]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N7MixedOnePoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  calc
    Wedge2Formalization.N7MixedOnePoint.ActBivector
        (paperRep₂ (k := k)) (paperChange (k := k)) =
      Matrix.fromBlocks
        0
        0
        0
        (Wedge2Formalization.N6MixedOnePoint.ActBivector
          (Wedge2Formalization.Paper.N6.Row2.paperRep₂ (k := k))
          (Wedge2Formalization.Paper.N6.Row2.paperChange (k := k))) := by
            simpa [paperRep₂, paperChange] using
              (Wedge2Formalization.N7MixedOnePoint.act_embedded_fromBlocks_zeroUpperRight
                (k := k)
                (Ω := Wedge2Formalization.Paper.N6.Row2.paperRep₂ (k := k))
                (a := Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k))
                (C := 0)
                (D := Wedge2Formalization.Paper.N6.Row2.paperChange (k := k)))
    _ =
      Matrix.fromBlocks
        0
        0
        0
        (Wedge2Formalization.Paper.N6.Row2.rep₂ (k := k)) := by
          simpa [Wedge2Formalization.Paper.N6.Row2.rep₂] using
            congrArg
              (fun Ω :
                Matrix Wedge2Formalization.N6MixedOnePoint.V
                  Wedge2Formalization.N6MixedOnePoint.V k =>
                Matrix.fromBlocks
                  (0 : Matrix I I k)
                  (0 : Matrix I W k)
                  (0 : Matrix W I k)
                  Ω)
              (Wedge2Formalization.Paper.N6.Row2.paperRep₂_transport (k := k))
    _ = rep₂ (k := k) := by
          simpa [rep₂, Wedge2Formalization.N7MixedOnePoint.radMixedRep₂,
            Wedge2Formalization.Paper.N6.Row2.rep₂,
            Wedge2Formalization.N6MixedOnePoint.rep₂]

/-- The displayed unipotent radical `\Ga^6 \rtimes U_{2,1}(k)`. -/
def U
    (C : Matrix W I k)
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k))
    0
    C
    (Wedge2Formalization.Paper.N6.Row2.U (k := k) X)

/-- The displayed Levi `\GL_1(k) \times (\SL_2(k) \times \SL_2(k))`. -/
def Levi
    (a0 : k)
    (A B H : Matrix Wedge2Formalization.N4.I
      Wedge2Formalization.N4.I k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
    0
    0
    (Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H)

/-- Exact pointwise kernel family for the mixed one-point radical extension. -/
def K : Set (Matrix V V k) :=
  { g |
      ∃ a0 : k,
        ∃ C : Matrix W I k,
          ∃ X : Matrix Wedge2Formalization.N6MixedOnePoint.I
              Wedge2Formalization.N6MixedOnePoint.I k,
            ∃ A B H : Matrix Wedge2Formalization.N4.I
                Wedge2Formalization.N4.I k,
              A.det = 1 ∧
                A * Wedge2Formalization.N4.J * Bᵀ +
                  B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
                H.det = 1 ∧
                g =
                  Matrix.fromBlocks
                    (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
                    0
                    C
                    (Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
                      Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) }

/-- The actual `GL_7` paper kernel family, obtained by adding the missing nonzero
`\GL_1(k)` coordinate to the raw pointwise kernel. -/
def KGL : Set (Matrix V V k) :=
  { g | g ∈ K (k := k) ∧ g (Sum.inl (0 : I)) (Sum.inl (0 : I)) ≠ 0 }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N6.Row2.Qproj (k := k)

def lift
    (a b : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k))
    0
    0
    (Wedge2Formalization.Paper.N6.Row2.lift (k := k) a b)

private theorem scalarBlock_eq
    (A : Matrix I I k) :
    Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (A 0 0) = A := by
  ext i j
  fin_cases i <;> fin_cases j
  simp [Wedge2Formalization.N7MixedOnePoint.scalarBlock]

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N7MixedOnePoint.FixesRadMixedPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N7Summary.radMixed_pointwise_bivector_iff
        (k := k)
        (a0 := g.toBlocks₁₁)
        (R := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hg) with
      ⟨hR0, hDfix⟩
    rcases
      (Wedge2Formalization.Paper.N6.Row2.mem_K_iff
        (k := k) (g := g.toBlocks₂₂)).1
        ((Wedge2Formalization.Paper.N6.Row2.pointwise_stabilizer
          (k := k) (g := g.toBlocks₂₂)).1 hDfix) with
      ⟨X, A, B, H, hA, hB, hH, hshape⟩
    refine ⟨g.toBlocks₁₁ 0 0, g.toBlocks₂₁, X, A, B, H, hA, hB, hH, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ g.toBlocks₂₂ := by
            simp [hR0]
      _ =
        Matrix.fromBlocks
          (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (g.toBlocks₁₁ 0 0))
          0
          g.toBlocks₂₁
          g.toBlocks₂₂ := by
            simp [scalarBlock_eq (k := k) g.toBlocks₁₁]
      _ =
        Matrix.fromBlocks
          (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (g.toBlocks₁₁ 0 0))
          0
          g.toBlocks₂₁
          (Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) := by
            rw [hshape]
  · rintro ⟨a0, C, X, A, B, H, hA, hB, hH, rfl⟩
    have hDmem :
        Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H ∈
          Wedge2Formalization.Paper.N6.Row2.K (k := k) :=
      (Wedge2Formalization.Paper.N6.Row2.mem_K_iff
        (k := k)
        (g :=
          Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H)).2
        ⟨X, A, B, H, hA, hB, hH, rfl⟩
    have hDfix :
        Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
          (Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) :=
      (Wedge2Formalization.Paper.N6.Row2.pointwise_stabilizer
        (k := k)
        (g :=
          Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H)).2
        hDmem
    exact
      (Wedge2Formalization.N7MixedOnePoint.fixesRadMixedPair_fromBlocks_zeroUpperRight_iff
        (k := k)
        (a := Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
        (C := C)
        (D :=
          Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H)).2
        hDfix

/-- Exact pointwise stabilizer statement on the actual `GL_7` paper kernel. -/
theorem pointwise_stabilizer_GL :
    ExactFamily
      (fun g =>
        Wedge2Formalization.N7MixedOnePoint.FixesRadMixedPairBivector (k := k) g ∧
          g (Sum.inl (0 : I)) (Sum.inl (0 : I)) ≠ 0)
      (KGL (k := k)) := by
  intro g
  constructor
  · rintro ⟨hfix, h00⟩
    exact ⟨(pointwise_stabilizer (k := k) (g := g)).1 hfix, h00⟩
  · rintro ⟨hg, h00⟩
    exact ⟨(pointwise_stabilizer (k := k) (g := g)).2 hg, h00⟩

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ a0 : k,
        ∃ C : Matrix W I k,
          ∃ X : Matrix Wedge2Formalization.N6MixedOnePoint.I
              Wedge2Formalization.N6MixedOnePoint.I k,
            ∃ A B H : Matrix Wedge2Formalization.N4.I
                Wedge2Formalization.N4.I k,
              A.det = 1 ∧
                A * Wedge2Formalization.N4.J * Bᵀ +
                  B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
                H.det = 1 ∧
                g =
                  Matrix.fromBlocks
                    (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
                    0
                    C
                    (Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
                      Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) := by
  rfl

theorem U_pointwise
    (C : Matrix W I k)
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k) :
    Wedge2Formalization.N7MixedOnePoint.FixesRadMixedPairBivector
      (U (k := k) C X) := by
  have hDfix :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.Paper.N6.Row2.U (k := k) X) :=
    Wedge2Formalization.Paper.N6.Row2.U_pointwise (k := k) X
  exact
    (Wedge2Formalization.N7MixedOnePoint.fixesRadMixedPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (a := Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k))
      (C := C)
      (D := Wedge2Formalization.Paper.N6.Row2.U (k := k) X)).2
      hDfix

theorem Levi_pointwise
    (a0 : k)
    (A B H : Matrix Wedge2Formalization.N4.I
      Wedge2Formalization.N4.I k)
    (hA : A.det = 1)
    (hB :
      A * Wedge2Formalization.N4.J * Bᵀ +
        B * Wedge2Formalization.N4.J * Aᵀ = 0)
    (hH : H.det = 1) :
    Wedge2Formalization.N7MixedOnePoint.FixesRadMixedPairBivector
      (Levi (k := k) a0 A B H) := by
  have hDfix :
      Wedge2Formalization.N6MixedOnePoint.FixesPairBivector
        (Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) :=
    Wedge2Formalization.Paper.N6.Row2.Levi_pointwise (k := k) A B H hA hB hH
  exact
    (Wedge2Formalization.N7MixedOnePoint.fixesRadMixedPair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (a := Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
      (C := 0)
      (D := Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H)).2
      hDfix

private theorem scalarBlock_one_mul
    (a0 : k) :
    Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k) *
        Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0 =
      Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0 := by
  ext i j
  fin_cases i
  fin_cases j
  simp [Wedge2Formalization.N7MixedOnePoint.scalarBlock, Matrix.mul_apply,
    Fin.sum_univ_one]

private theorem mul_scalarBlock_eq_smul
    (C : Matrix W I k)
    (a0 : k) :
    C * Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0 = a0 • C := by
  ext i j
  fin_cases j
  simp [Wedge2Formalization.N7MixedOnePoint.scalarBlock, Matrix.mul_apply,
    Fin.sum_univ_one, smul_eq_mul, mul_comm]

private theorem U_mul_Levi_raw
    (C : Matrix W I k)
    (X : Matrix Wedge2Formalization.N6MixedOnePoint.I
      Wedge2Formalization.N6MixedOnePoint.I k)
    (a0 : k)
    (A B H : Matrix Wedge2Formalization.N4.I
      Wedge2Formalization.N4.I k) :
    U (k := k) C X * Levi (k := k) a0 A B H =
      Matrix.fromBlocks
        (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
        0
        (a0 • C)
        (Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
          Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) := by
  calc
    U (k := k) C X * Levi (k := k) a0 A B H =
        Matrix.fromBlocks
          ((Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k)) *
            (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0) +
            (0 : Matrix I I k) * 0)
          ((Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) (1 : k)) *
              (0 : Matrix I W k) +
            (0 : Matrix I W k) *
              (Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H))
          (C * Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0 +
            Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
              (0 : Matrix W I k))
          (C * (0 : Matrix I W k) +
            Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
              Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) := by
          simp [U, Levi, Matrix.fromBlocks_multiply]
    _ =
        Matrix.fromBlocks
          (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
          0
          (a0 • C)
          (Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) := by
          simp [scalarBlock_one_mul, mul_scalarBlock_eq_smul, zero_mul, mul_zero, add_zero]

/-- Table-facing kernel statement for Appendix A, row 13, expressed on the actual
`GL_7` paper kernel. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ a0 : k,
        ∃ C : Matrix W I k,
          ∃ X : Matrix Wedge2Formalization.N6MixedOnePoint.I
              Wedge2Formalization.N6MixedOnePoint.I k,
            ∃ A B H : Matrix Wedge2Formalization.N4.I
                Wedge2Formalization.N4.I k,
              a0 ≠ 0 ∧
                A.det = 1 ∧
                A * Wedge2Formalization.N4.J * Bᵀ +
                  B * Wedge2Formalization.N4.J * Aᵀ = 0 ∧
                H.det = 1 ∧
                g = U (k := k) C X * Levi (k := k) a0 A B H }

/-- Table-facing kernel statement for Appendix A, row 13, expressed on the actual
`GL_7` paper kernel. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ KGL (k := k) ↔
      g ∈ TableCell (k := k) := by
  constructor
  · rintro ⟨hg, h00⟩
    rcases (mem_K_iff (k := k) (g := g)).1 hg with
      ⟨a0, C0, X, A, B, H, hA, hB, hH, hshape⟩
    have ha0 : a0 ≠ 0 := by
      simpa [hshape, Wedge2Formalization.N7MixedOnePoint.scalarBlock] using h00
    refine ⟨a0, a0⁻¹ • C0, X, A, B, H, ha0, hA, hB, hH, ?_⟩
    calc
      g =
        Matrix.fromBlocks
          (Wedge2Formalization.N7MixedOnePoint.scalarBlock (k := k) a0)
          0
          C0
          (Wedge2Formalization.Paper.N6.Row2.U (k := k) X *
            Wedge2Formalization.Paper.N6.Row2.Levi (k := k) A B H) := hshape
      _ =
        U (k := k) (a0⁻¹ • C0) X * Levi (k := k) a0 A B H := by
          symm
          simpa [ha0] using
            (U_mul_Levi_raw
              (k := k)
              (C := a0⁻¹ • C0)
              (X := X)
              (a0 := a0)
              (A := A)
              (B := B)
              (H := H))
  · rintro ⟨a0, C, X, A, B, H, ha0, hA, hB, hH, rfl⟩
    refine ⟨?_, ?_⟩
    · exact
        (mem_K_iff (k := k)
          (g := U (k := k) C X * Levi (k := k) a0 A B H)).2
          ⟨a0, a0 • C, X, A, B, H, hA, hB, hH,
            by simpa using
              (U_mul_Levi_raw
                (k := k)
                (C := C)
                (X := X)
                (a0 := a0)
                (A := A)
                (B := B)
                (H := H))⟩
    · rw [U_mul_Levi_raw (k := k) (C := C) (X := X) (a0 := a0) (A := A) (B := B) (H := H)]
      simp [ha0, Wedge2Formalization.N7MixedOnePoint.scalarBlock]

theorem det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂
          (a • rep₁ (k := k) + b • rep₂ (k := k))) = 0 ↔
      a = 0 :=
  Wedge2Formalization.N7Summary.radMixed_lowerRight_det_zero_iff
    (k := k) (a := a) (b := b)

theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N7MixedOnePoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a b)
      (Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b) := by
  constructor
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.Paper.N6.Row2.lift] using
      (Wedge2Formalization.N7Summary.radMixed_borel_lift_action
        (k := k)
        (a0 := (1 : k))
        (a := a)
        (b := b)
        ha
        (C := 0)).1
  ·
    simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.Paper.N6.Row2.lift,
      Wedge2Formalization.N4PaperSummary.row1_borelCoeff] using
      (Wedge2Formalization.N7Summary.radMixed_borel_lift_action
        (k := k)
        (a0 := (1 : k))
        (a := a)
        (b := b)
        ha
        (C := 0)).2

private theorem act_embedded_toBlocks₂₂_mixed
    (Ω : Matrix W W k)
    (g : Matrix V V k) :
    Matrix.toBlocks₂₂
        (Wedge2Formalization.N7MixedOnePoint.ActBivector
          (Matrix.fromBlocks 0 0 0 Ω) g) =
      Wedge2Formalization.N6MixedOnePoint.ActBivector Ω g.toBlocks₂₂ := by
  rw [← Matrix.fromBlocks_toBlocks g]
  simp [Wedge2Formalization.N7MixedOnePoint.ActBivector,
    Wedge2Formalization.N6MixedOnePoint.ActBivector,
    Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' :
      a • (Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k)) +
          b • (Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7MixedOnePoint.radMixedRep₁,
      Wedge2Formalization.N7MixedOnePoint.radMixedRep₂] using
      congrArg Matrix.toBlocks₂₂ h
  have ha : a = 0 := by
    simpa [Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N6MixedOnePoint.rep₂, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using
      (congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h')
  have hb : b = 0 := by
    simpa [ha, Wedge2Formalization.N6MixedOnePoint.rep₁,
      Wedge2Formalization.N6MixedOnePoint.rep₂, Wedge2Formalization.N4.onePointRep₁,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using
      (congrArg (fun M => M (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inl 1))) h')
  exact ⟨ha, hb⟩

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  exact one_ne_zero (rep_pair_independent (k := k) hcomb).1

theorem quotient_image
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N7MixedOnePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N7MixedOnePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have hdetM : Matrix.det M ≠ 0 := by
    by_contra hdet
    by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
    · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
      · have hzero :
            Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using hM.1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := rep₁ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
          simpa using horig
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).1
      · have hzero :
            Wedge2Formalization.N7MixedOnePoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N7MixedOnePoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) •
                  Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₁ (k := k)) g +
                M 0 1 •
                  Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₂ (k := k)) g := by
                    simp [Wedge2Formalization.N7MixedOnePoint.ActBivector, Matrix.mul_add,
                      Matrix.add_mul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
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
        rcases rep_pair_independent (k := k) horig with ⟨h11, h01⟩
        exact hcol1 ⟨h01, by simpa using neg_eq_zero.mp h11⟩
    · have hzero :
          Wedge2Formalization.N7MixedOnePoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
        have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
          exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hdet)
        calc
          Wedge2Formalization.N7MixedOnePoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
              =
            (-(M 1 0)) •
                Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₁ (k := k)) g +
              M 0 0 •
                Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7MixedOnePoint.ActBivector, Matrix.mul_add,
                    Matrix.add_mul]
          _ =
            (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [h1, h2]
          _ = 0 := by
                ext i j
                simp [Matrix.add_apply, Matrix.smul_apply]
                have hdet'' :
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j) =
                      M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by
                  calc
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j)
                        = (M 0 0 * M 1 1) * rep₂ (k := k) i j := by ring
                    _ = (M 0 1 * M 1 0) * rep₂ (k := k) i j := by rw [hdet']
                    _ = M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by ring
                rw [hdet'']
                ring_nf
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k)
          (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
          (g := g)
          hg).1 hzero
      rcases rep_pair_independent (k := k) horig with ⟨h10, h00⟩
      exact hcol0 ⟨h00, by simpa using neg_eq_zero.mp h10⟩
  have h2blocks :
      Matrix.toBlocks₂₂
          (Wedge2Formalization.N7MixedOnePoint.ActBivector (rep₂ (k := k)) g) =
        M 1 0 • (Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k)) +
          M 1 1 • (Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k)) := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7MixedOnePoint.radMixedRep₁,
      Wedge2Formalization.N7MixedOnePoint.radMixedRep₂] using
      congrArg Matrix.toBlocks₂₂ h2
  have hrep₂det :
      Matrix.det
        (M 1 0 • (Wedge2Formalization.N6MixedOnePoint.rep₁ (k := k)) +
          M 1 1 • (Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k))) = 0 := by
    have hbase :
        Matrix.det (Wedge2Formalization.N6MixedOnePoint.rep₂ (k := k)) = 0 := by
      simpa [Wedge2Formalization.Paper.N6.Row2.rep₂] using
        ((Wedge2Formalization.Paper.N6.Row2.det_zero_iff
          (k := k) (a := (0 : k)) (b := (1 : k))).2 rfl)
    rw [← h2blocks]
    rw [rep₂, Wedge2Formalization.N7MixedOnePoint.radMixedRep₂,
      act_embedded_toBlocks₂₂_mixed, Wedge2Formalization.N6MixedOnePoint.ActBivector,
      Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose]
    simp [hbase]
  have h20 : M 1 0 = 0 := by
    exact
      (Wedge2Formalization.Paper.N6.Row2.det_zero_iff
        (k := k) (a := M 1 0) (b := M 1 1)).1
        (by simpa [h2blocks] using hrep₂det)
  refine ⟨M, ?_, hM⟩
  exact ⟨h20, by simpa [Matrix.det_fin_two, h20] using hdetM⟩

end Row13

/-! Appendix A, `n = 7`, row 4.
Representative `S_2 + S_1^4`.
Divisor `0`.
Claimed stabilizer:
`K_L = \Ga^{12} \rtimes (\GL_4(k) \times (\Ga^2 \rtimes \Gm(k)))`, exact
quotient family `Q_L = GL_2(k)`.
-/
namespace Row4

abbrev I := Wedge2Formalization.N7PureSingular.I
abbrev W := Wedge2Formalization.N7PureSingular.W
abbrev V := Wedge2Formalization.N7PureSingular.V

def rep₁ : Matrix V V k :=
  Wedge2Formalization.N7PureSingular.rad3PureSingularRep₁ (k := k)

def rep₂ : Matrix V V k :=
  Wedge2Formalization.N7PureSingular.rad3PureSingularRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 4:
`e₁ \wedge e₃`. -/
def paperRep₁ : Matrix V V k :=
  let Ω : Matrix I I k := !![(0 : k), 0, 1; 0, 0, 0; -1, 0, 0]
  Matrix.fromBlocks Ω 0 0 0

/-- Literal second basis vector from Appendix A, row 4:
`e₂ \wedge e₃`. -/
def paperRep₂ : Matrix V V k :=
  let Ω : Matrix I I k := !![(0 : k), 0, 0; 0, 0, 1; 0, -1, 0]
  Matrix.fromBlocks Ω 0 0 0

/-- Explicit basis change sending the literal paper pair to the internal working one. -/
def paperChange : Matrix V V k :=
  let A : Matrix I I k := 0
  let B : Matrix I W k := !![(1 : k), 0, 0, 0;
                              0, 1, 0, 0;
                              0, 0, 1, 0]
  let C : Matrix W I k := !![(0 : k), 0, 1;
                              -1, 0, 0;
                              0, -1, 0;
                              0, 0, 0]
  let D : Matrix W W k := !![(0 : k), 0, 0, 0;
                              0, 0, 0, 0;
                              0, 0, 0, 0;
                              0, 0, 0, 1]
  Matrix.fromBlocks A B C D

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N7PureSingular.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N7PureSingular.rad3PureSingularRep₁,
      Wedge2Formalization.N7PureSingular.ActBivector, Wedge2Formalization.N4PureSingular.ω12,
      Matrix.fromBlocks,
      Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_four]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N7PureSingular.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  ext i j
  rcases i with i | i <;> rcases j with j | j <;>
    fin_cases i <;> fin_cases j <;>
    simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N7PureSingular.rad3PureSingularRep₂,
      Wedge2Formalization.N7PureSingular.ActBivector, Wedge2Formalization.N4PureSingular.ω13,
      Matrix.fromBlocks,
      Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_four]

def K : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₂ =
        Wedge2Formalization.N7PureSingular.upperRightLast (k := k)
          (g.toBlocks₁₂ 0 3) (g.toBlocks₁₂ 1 3) (g.toBlocks₁₂ 2 3) ∧
      g.toBlocks₂₂ 0 0 ≠ 0 ∧
      g.toBlocks₂₂ =
        Wedge2Formalization.N4PureSingular.pureSingularShape
          (k := k)
          (g.toBlocks₂₂ 0 0)
          (g.toBlocks₂₂ 0 1)
          (g.toBlocks₂₂ 0 2)
          (g.toBlocks₂₂ 0 3)
          (g.toBlocks₂₂ 1 3)
          (g.toBlocks₂₂ 2 3)
          (g.toBlocks₂₂ 3 3) }

/-- Public name for the embedded `U_5 \rtimes (G_m(k) \times G_m(k))` core on the
lower-right `4`-space in Appendix A, row 4. -/
def PureSingularCore : Set (Matrix W W k) :=
  { D |
      D 0 0 ≠ 0 ∧
        D =
          Wedge2Formalization.N4PureSingular.pureSingularShape
            (k := k)
            (D 0 0)
            (D 0 1)
            (D 0 2)
            (D 0 3)
            (D 1 3)
            (D 2 3)
            (D 3 3) }

/-- Public name for the upper-right radical block in Appendix A, row 4. -/
def UpperRightRadical : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₂ =
        Wedge2Formalization.N7PureSingular.upperRightLast (k := k)
          (g.toBlocks₁₂ 0 3) (g.toBlocks₁₂ 1 3) (g.toBlocks₁₂ 2 3) }

/-- Public name for the explicit Appendix A row-4 kernel cell. -/
def TableCell : Set (Matrix V V k) :=
  { g | g ∈ UpperRightRadical (k := k) ∧ g.toBlocks₂₂ ∈ PureSingularCore (k := k) }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

def lift
    (α β γ δ : k) :
    Matrix V V k :=
  Matrix.fromBlocks
    (1 : Matrix I I k)
    0
    0
    (Wedge2Formalization.N4PureSingular.pureSingularSetwiseShape
      (k := k) (1 : k) 0 0 0 α γ β δ 0 0 1)

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N7PureSingular.FixesRad3PureSingularPairBivector (k := k))
      (K (k := k)) := by
  intro g
  have hgBlocks :
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
    simpa using (Matrix.fromBlocks_toBlocks g).symm
  constructor
  · intro hg
    have hfixBlocks :
        Wedge2Formalization.N7PureSingular.FixesRad3PureSingularPairBivector (k := k)
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) := by
      rw [← hgBlocks]
      exact hg
    rcases
      (Wedge2Formalization.N7Summary.rad3PureSingular_pointwise_bivector_iff
        (k := k)
        (A0 := g.toBlocks₁₁)
        (R := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1
        hfixBlocks with
      ⟨hR, h00, hD⟩
    exact ⟨hR, h00, hD⟩
  · rintro ⟨hR, h00, hD⟩
    have hfix :
        Wedge2Formalization.N7PureSingular.FixesRad3PureSingularPairBivector (k := k)
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) :=
      (Wedge2Formalization.N7Summary.rad3PureSingular_pointwise_bivector_iff
        (k := k)
        (A0 := g.toBlocks₁₁)
        (R := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).2 ⟨hR, h00, hD⟩
    rw [hgBlocks]
    exact hfix

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g.toBlocks₁₂ =
        Wedge2Formalization.N7PureSingular.upperRightLast (k := k)
          (g.toBlocks₁₂ 0 3) (g.toBlocks₁₂ 1 3) (g.toBlocks₁₂ 2 3) ∧
      g.toBlocks₂₂ 0 0 ≠ 0 ∧
      g.toBlocks₂₂ =
        Wedge2Formalization.N4PureSingular.pureSingularShape
          (k := k)
          (g.toBlocks₂₂ 0 0)
          (g.toBlocks₂₂ 0 1)
          (g.toBlocks₂₂ 0 2)
          (g.toBlocks₂₂ 0 3)
          (g.toBlocks₂₂ 1 3)
          (g.toBlocks₂₂ 2 3)
          (g.toBlocks₂₂ 3 3) := by
  rfl

/-- Table-facing kernel statement for Appendix A, row 4. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) := by
  constructor
  · rintro ⟨hR, h00, hD⟩
    exact ⟨hR, h00, hD⟩
  · rintro ⟨hR, hD⟩
    exact ⟨hR, hD.1, hD.2⟩

theorem quotient_action
    (α β γ δ : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N7PureSingular.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) α β γ δ)
      (!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N7Summary.rad3PureSingular_GL2_lift_action
        (k := k)
        (α := α)
        (β := β)
        (γ := γ)
        (δ := δ)
        (A0 := (1 : Matrix I I k))
        (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N7Summary.rad3PureSingular_GL2_lift_action
        (k := k)
        (α := α)
        (β := β)
        (γ := γ)
        (δ := δ)
        (A0 := (1 : Matrix I I k))
        (C := 0)).2

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h' := congrArg Matrix.toBlocks₂₂ h
  have h'' :
      a • (Wedge2Formalization.N4PureSingular.ω12 (k := k)) +
          b • (Wedge2Formalization.N4PureSingular.ω13 (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7PureSingular.rad3PureSingularRep₁,
      Wedge2Formalization.N7PureSingular.rad3PureSingularRep₂] using h'
  have ha : a = 0 := by
    have h01 := congrArg (fun M => M 0 1) h''
    simpa [Wedge2Formalization.N4PureSingular.ω12,
      Wedge2Formalization.N4PureSingular.ω13, Matrix.add_apply, Matrix.smul_apply] using h01
  have hb : b = 0 := by
    have h02 := congrArg (fun M => M 0 2) h''
    simpa [ha, Wedge2Formalization.N4PureSingular.ω12,
      Wedge2Formalization.N4PureSingular.ω13, Matrix.add_apply, Matrix.smul_apply] using h02
  exact ⟨ha, hb⟩

theorem quotient_image
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N7PureSingular.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N7PureSingular.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      Wedge2Formalization.N7PureSingular.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      Wedge2Formalization.N7PureSingular.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have hdetM : Matrix.det M ≠ 0 := by
    by_contra hdet
    by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
    · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
      · have hzero :
            Wedge2Formalization.N7PureSingular.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using h1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := rep₁ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
          simpa using horig
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).1
      · have hzero :
            Wedge2Formalization.N7PureSingular.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N7PureSingular.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) • Wedge2Formalization.N7PureSingular.ActBivector (rep₁ (k := k)) g +
                M 0 1 • Wedge2Formalization.N7PureSingular.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7PureSingular.ActBivector, Matrix.mul_add, Matrix.add_mul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
            _ = 0 := by
                  ext i j
                  simp [hcol0.1, hcol0.2, Matrix.add_apply, Matrix.smul_apply]
                  ring_nf
        have horig :
            (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k) = 0 :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        rcases rep_pair_independent (k := k) horig with ⟨h11, h01⟩
        exact hcol1 ⟨h01, by simpa using neg_eq_zero.mp h11⟩
    · exact by
        have hzero :
            Wedge2Formalization.N7PureSingular.ActBivector
                ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
          have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
            exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hdet)
          calc
            Wedge2Formalization.N7PureSingular.ActBivector
                ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
                =
              (-(M 1 0)) • Wedge2Formalization.N7PureSingular.ActBivector (rep₁ (k := k)) g +
                M 0 0 • Wedge2Formalization.N7PureSingular.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7PureSingular.ActBivector, Matrix.mul_add, Matrix.add_mul]
            _ =
              (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
            _ = 0 := by
                  ext i j
                  simp [Matrix.add_apply, Matrix.smul_apply]
                  have hdet'' :
                      M 0 0 * (M 1 1 * rep₂ (k := k) i j) =
                        M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by
                    calc
                      M 0 0 * (M 1 1 * rep₂ (k := k) i j)
                          = (M 0 0 * M 1 1) * rep₂ (k := k) i j := by ring
                      _ = (M 0 1 * M 1 0) * rep₂ (k := k) i j := by rw [hdet']
                      _ = M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by ring
                  rw [hdet'']
                  ring_nf
        have horig :
            (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k) = 0 :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
            (g := g)
            hg).1 hzero
        rcases rep_pair_independent (k := k) horig with ⟨h10, h00⟩
        exact hcol0 ⟨h00, by simpa using neg_eq_zero.mp h10⟩
  exact ⟨M, hdetM, hM⟩

end Row4

/-! Appendix A, `n = 7`, row 15.
Representative `S_1 + 2J_{a,1} + J_{b,1}`.
Divisor `[a] + [b]`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times (\Sp_4(k) \times \SL_2(k)))`, exact
quotient family `Q_L = T`.
-/
namespace Row15

abbrev I := Wedge2Formalization.N7WeightedTwoPoint.I
abbrev W := Wedge2Formalization.N7WeightedTwoPoint.W
abbrev V := Wedge2Formalization.N7WeightedTwoPoint.V

def rep₁ : Matrix V V k :=
  Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₁ (k := k)

def rep₂ : Matrix V V k :=
  Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 15. -/
def paperRep₁ : Matrix V V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 15. -/
def paperRep₂ : Matrix V V k :=
  rep₂ (k := k)

/-- The row-15 paper representative already agrees with the internal working one. -/
def paperChange : Matrix V V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N7WeightedTwoPoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N7WeightedTwoPoint.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N7WeightedTwoPoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N7WeightedTwoPoint.ActBivector]

/-- The displayed kernel family: free radical scalar and column, together with the
weighted two-point pointwise kernel on the `6`-dimensional quotient. -/
def Levi
    (a0 : Matrix I I k)
    (C : Matrix W I k)
    (A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
      Wedge2Formalization.N6WeightedTwoPoint.W k)
    (D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
      Wedge2Formalization.N6WeightedTwoPoint.I k) :
    Matrix V V k :=
  Matrix.fromBlocks a0 0 C (Matrix.fromBlocks A 0 0 D)

/-- Public name for the embedded `\Sp_4(k)` factor in Appendix A, row 15. -/
def Sp4Core :
    Set
      (Matrix Wedge2Formalization.N6WeightedTwoPoint.W
        Wedge2Formalization.N6WeightedTwoPoint.W k) :=
  Wedge2Formalization.Paper.N6.Row4.Sp4Core (k := k)

def K : Set (Matrix V V k) :=
  { g |
      ∃ a0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
              Wedge2Formalization.N6WeightedTwoPoint.W k,
            ∃ D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
                Wedge2Formalization.N6WeightedTwoPoint.I k,
              A ∈ Sp4Core (k := k) ∧
                D.det = 1 ∧
                g = Levi (k := k) a0 C A D }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N6.Row4.Qproj (k := k)

def lift
    (u v : k) :
    Matrix V V k :=
  let A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k := !![u, 0; 0, 1]
  let H : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
      Wedge2Formalization.N6WeightedTwoPoint.W k := Matrix.fromBlocks A 0 0 A
  let E : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
      Wedge2Formalization.N6WeightedTwoPoint.I k := !![v, 0; 0, 1]
  let h : Matrix W I k := 0
  Matrix.fromBlocks
    (Wedge2Formalization.N7WeightedTwoPoint.scalarBlock (k := k) (1 : k))
    0
    h
    (Matrix.fromBlocks H 0 0 E)

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N7WeightedTwoPoint.FixesRadWeightedPairBivector (k := k))
      (K (k := k)) := by
  intro g
  have hgBlocks :
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
    simpa using (Matrix.fromBlocks_toBlocks g).symm
  constructor
  · intro hg
    have hfixBlocks :
        Wedge2Formalization.N7WeightedTwoPoint.FixesRadWeightedPairBivector (k := k)
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) := by
      rw [← hgBlocks]
      exact hg
    rcases
      (Wedge2Formalization.N7WeightedTwoPoint.fixesRadWeightedPair_fromBlocks_iff
        (k := k)
        (a := g.toBlocks₁₁)
        (B := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1 hfixBlocks with
      ⟨hR, hDfix⟩
    rcases
      (Wedge2Formalization.Paper.N6.Row4.pointwise_stabilizer
        (k := k)
        (g := g.toBlocks₂₂)).1 hDfix with
      ⟨A, D, hA, hD, hshape⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, A, D, hA, hD, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁ (Matrix.fromBlocks A 0 0 D) := by
            simpa [Wedge2Formalization.Paper.N6.Row4.Levi, hR] using hshape
      _ = Levi (k := k) g.toBlocks₁₁ g.toBlocks₂₁ A D := by
            rfl
  · rintro ⟨a0, C, A, D, hA, hD, rfl⟩
    have h6 :
        Wedge2Formalization.N6WeightedTwoPoint.FixesPairBivector
          (Matrix.fromBlocks A 0 0 D) := by
      exact
        (Wedge2Formalization.N6Summary.weightedTwoPoint_pointwise_bivector_iff
          (k := k) (A := A) (B := 0) (C := 0) (D := D)).2 ⟨rfl, rfl, hA, hD⟩
    exact
      (Wedge2Formalization.N7WeightedTwoPoint.fixesRadWeightedPair_fromBlocks_iff
        (k := k)
        (a := a0)
        (B := 0)
        (C := C)
        (D := Matrix.fromBlocks A 0 0 D)).2 ⟨rfl, h6⟩

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ a0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
              Wedge2Formalization.N6WeightedTwoPoint.W k,
            ∃ D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
                Wedge2Formalization.N6WeightedTwoPoint.I k,
              A ∈ Sp4Core (k := k) ∧
                D.det = 1 ∧
                g = Levi (k := k) a0 C A D := by
  rfl

/-- Table-facing kernel statement for Appendix A, row 15. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ a0 : Matrix I I k,
        ∃ C : Matrix W I k,
        ∃ A : Matrix Wedge2Formalization.N6WeightedTwoPoint.W
            Wedge2Formalization.N6WeightedTwoPoint.W k,
          ∃ D : Matrix Wedge2Formalization.N6WeightedTwoPoint.I
              Wedge2Formalization.N6WeightedTwoPoint.I k,
              A ∈ Sp4Core (k := k) ∧
                D.det = 1 ∧
                g = Levi (k := k) a0 C A D }

/-- Table-facing kernel statement for Appendix A, row 15. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) :=
  mem_K_iff (k := k) (g := g)

theorem torus_action
    (u v : k) :
    ActsOnOrderedPair
      (Wedge2Formalization.N7WeightedTwoPoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) u v)
      (Wedge2Formalization.N4PaperSummary.row2_torusCoeff (k := k) u v) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N7Summary.radWeighted_torus_lift_action
        (k := k)
        (a0 := (1 : k))
        (u := u)
        (v := v)
        (C := 0)).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift,
      Wedge2Formalization.N4PaperSummary.row2_torusCoeff] using
      (Wedge2Formalization.N7Summary.radWeighted_torus_lift_action
        (k := k)
        (a0 := (1 : k))
        (u := u)
        (v := v)
        (C := 0)).2

/-- Rank drop on the weighted-support line is detected by the lower-right determinant. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂ (a • rep₁ (k := k) + b • rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  Wedge2Formalization.N7Summary.radWeighted_lowerRight_det_zero_iff
    (k := k) (a := a) (b := b)

private theorem act_embedded_toBlocks₂₂_weighted
    (Ω : Matrix Wedge2Formalization.N7WeightedTwoPoint.W
      Wedge2Formalization.N7WeightedTwoPoint.W k)
    (g : Matrix V V k) :
    Matrix.toBlocks₂₂
        (Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (Matrix.fromBlocks 0 0 0 Ω) g) =
      Wedge2Formalization.N6WeightedTwoPoint.ActBivector Ω g.toBlocks₂₂ := by
  rw [← Matrix.fromBlocks_toBlocks g]
  simp [Wedge2Formalization.N7WeightedTwoPoint.ActBivector,
    Wedge2Formalization.N6WeightedTwoPoint.ActBivector, Matrix.fromBlocks_transpose,
    Matrix.fromBlocks_multiply]

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have hlower :
      a • (Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)) +
          b • (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)) = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₁,
      Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂] using
      congrArg Matrix.toBlocks₂₂ h
  have htop := congrArg
      (fun M =>
        M (Sum.inl (Sum.inl (0 : Wedge2Formalization.N4.I)))
          (Sum.inl (Sum.inl (1 : Wedge2Formalization.N4.I)))) hlower
  have ha : a = 0 := by
    have htop' := htop
    simp [Wedge2Formalization.N6WeightedTwoPoint.rep₁,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, Wedge2Formalization.N4.splitRep₁,
      Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] at htop'
    exact htop'
  have hbot := congrArg
      (fun M =>
        M (Sum.inr (0 : Wedge2Formalization.N6WeightedTwoPoint.I))
          (Sum.inr (1 : Wedge2Formalization.N6WeightedTwoPoint.I))) hlower
  have hb : b = 0 := by
    simpa [ha, Wedge2Formalization.N6WeightedTwoPoint.rep₁,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂,
      Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.J, Matrix.fromBlocks,
      Matrix.add_apply, Matrix.smul_apply] using hbot
  exact ⟨ha, hb⟩

private def rep₂_left :
    Matrix V Wedge2Formalization.N4.I k :=
  fun i j =>
    match i with
    | Sum.inr (Sum.inr i) => if i = j then 1 else 0
    | _ => 0

private def rep₂_right :
    Matrix Wedge2Formalization.N4.I V k :=
  fun i j =>
    match j with
    | Sum.inr (Sum.inr j) => Wedge2Formalization.N4.J (k := k) i j
    | _ => 0

private theorem rep₂_factor :
    rep₂ (k := k) = rep₂_left (k := k) * rep₂_right (k := k) := by
  ext i j
  rcases i with _ | (_ | i) <;> rcases j with _ | (_ | j)
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
      Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
      Fin.sum_univ_two]
  · fin_cases i <;> fin_cases j <;>
      simp [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
        Wedge2Formalization.N6WeightedTwoPoint.rep₂, rep₂_left, rep₂_right,
        Wedge2Formalization.N4.J, Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type,
        Fin.sum_univ_two]

private theorem rep₂_rank_le_two :
    (rep₂ (k := k)).rank ≤ 2 := by
  rw [rep₂_factor (k := k)]
  exact (Matrix.rank_mul_le_left (rep₂_left (k := k)) (rep₂_right (k := k))).trans <| by
    simpa using (Matrix.rank_le_card_width (rep₂_left (k := k)))

private theorem rep₂_ne_zero :
    rep₂ (k := k) ≠ 0 := by
  intro hzero
  have hcoord := congrArg (fun M => M (Sum.inr (Sum.inr 0)) (Sum.inr (Sum.inr 1))) hzero
  simpa [rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
    Wedge2Formalization.N6WeightedTwoPoint.rep₂, Wedge2Formalization.N4.J,
    Matrix.fromBlocks] using hcoord

private def topProj :
    Matrix Wedge2Formalization.N4.V V k :=
  fun i j =>
    match j with
    | Sum.inl _ => 0
    | Sum.inr (Sum.inl j) => if i = j then 1 else 0
    | Sum.inr (Sum.inr _) => 0

private def topIncl :
    Matrix V Wedge2Formalization.N4.V k :=
  fun i j =>
    match i with
    | Sum.inl _ => 0
    | Sum.inr (Sum.inl i) => if i = j then 1 else 0
    | Sum.inr (Sum.inr _) => 0

private theorem topProj_mul_topRep_mul_topIncl :
    topProj (k := k) * (rep₁ (k := k) - rep₂ (k := k)) * topIncl (k := k) =
      Wedge2Formalization.N4.splitRep₁ (k := k) := by
  ext i j
  cases i <;> cases j <;>
    simp [topProj, topIncl, rep₁, rep₂, sub_eq_add_neg,
      Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₁,
      Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
      Wedge2Formalization.N6WeightedTwoPoint.rep₁,
      Wedge2Formalization.N6WeightedTwoPoint.rep₂,
      Wedge2Formalization.N4.splitRep₁, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.mul_apply, Fintype.sum_sum_type, Fin.sum_univ_two]

private theorem topCombo_rank_ge_four
    (a : k)
    (ha : a ≠ 0) :
    4 ≤ (a • (rep₁ (k := k) - rep₂ (k := k))).rank := by
  have hUnit : IsUnit (Matrix.det (a • Wedge2Formalization.N4.splitRep₁ (k := k))) := by
    refine isUnit_iff_ne_zero.mpr ?_
    rw [Matrix.det_smul, Wedge2Formalization.N6WeightedTwoPoint.det_splitRep₁]
    simp [ha]
  calc
    4 = (a • Wedge2Formalization.N4.splitRep₁ (k := k)).rank := by
      simpa [Wedge2Formalization.N4.V] using
        (Matrix.rank_of_isUnit
          (a • Wedge2Formalization.N4.splitRep₁ (k := k))
          ((a • Wedge2Formalization.N4.splitRep₁ (k := k)).isUnit_iff_isUnit_det.mpr hUnit)).symm
    _ = (topProj (k := k) * (a • (rep₁ (k := k) - rep₂ (k := k))) * topIncl (k := k)).rank := by
      have hEq :
          topProj (k := k) * (a • (rep₁ (k := k) - rep₂ (k := k))) * topIncl (k := k) =
            a • Wedge2Formalization.N4.splitRep₁ (k := k) := by
        calc
          topProj (k := k) * (a • (rep₁ (k := k) - rep₂ (k := k))) * topIncl (k := k)
              = a • (topProj (k := k) * (rep₁ (k := k) - rep₂ (k := k)) * topIncl (k := k)) := by
                  simp [Matrix.mul_assoc]
          _ = a • Wedge2Formalization.N4.splitRep₁ (k := k) := by
                rw [topProj_mul_topRep_mul_topIncl (k := k)]
      rw [hEq]
    _ ≤ (a • (rep₁ (k := k) - rep₂ (k := k))).rank := by
      exact (Matrix.rank_mul_le_left _ _).trans <| Matrix.rank_mul_le_right _ _

private theorem act_rank_eq
    (Ω : Matrix V V k)
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0) :
    (Wedge2Formalization.N7WeightedTwoPoint.ActBivector Ω g).rank = Ω.rank := by
  have hg_unit : IsUnit (Matrix.det g) := isUnit_iff_ne_zero.mpr hg
  have hgt_unit : IsUnit (Matrix.det gᵀ) := by
    simpa [Matrix.det_transpose] using hg_unit
  calc
    (Wedge2Formalization.N7WeightedTwoPoint.ActBivector Ω g).rank
        = (g * Ω * gᵀ).rank := by rfl
    _ = (g * Ω).rank := by
          simpa [Wedge2Formalization.N7WeightedTwoPoint.ActBivector, Matrix.mul_assoc] using
            (Matrix.rank_mul_eq_left_of_isUnit_det (A := gᵀ) (B := g * Ω) hgt_unit)
    _ = Ω.rank := by
          simpa using (Matrix.rank_mul_eq_right_of_isUnit_det (A := g) (B := Ω) hg_unit)

private theorem combo_rank_le_two_of_leftCoeff_zero
    (c : k) :
    ((0 : k) • rep₁ (k := k) + c • rep₂ (k := k)).rank ≤ 2 := by
  rw [zero_smul, zero_add]
  have hfactor :
      c • rep₂ (k := k) = (c • rep₂_left (k := k)) * rep₂_right (k := k) := by
    rw [rep₂_factor (k := k)]
    ext i j
    simp [Matrix.mul_apply, Fin.sum_univ_two]
    ring
  rw [hfactor]
  exact (Matrix.rank_mul_le_left (c • rep₂_left (k := k)) (rep₂_right (k := k))).trans <| by
    simpa using (Matrix.rank_le_card_width (c • rep₂_left (k := k)))

theorem quotient_image
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N7WeightedTwoPoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N7WeightedTwoPoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have h12 :
      Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g =
        (M 0 0 - M 1 0) • rep₁ (k := k) +
          (M 0 1 - M 1 1) • rep₂ (k := k) := by
    calc
      Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g
          =
        Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) g := by
            simp [sub_eq_add_neg]
      _ =
        Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₁ (k := k)) g +
          (-1 : k) • Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g := by
            simp [Wedge2Formalization.N7WeightedTwoPoint.ActBivector, Matrix.mul_add, Matrix.add_mul]
      _ =
        (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
          (-1 : k) • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
            rw [h1, h2]
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
            Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using hM.1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := rep₁ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
          simpa using horig
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).1
      · have hzero :
            Wedge2Formalization.N7WeightedTwoPoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N7WeightedTwoPoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) • Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₁ (k := k)) g +
                M 0 1 • Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7WeightedTwoPoint.ActBivector, Matrix.mul_add,
                    Matrix.add_mul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
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
        rcases rep_pair_independent (k := k) horig with ⟨h11, h01⟩
        exact hcol1 ⟨h01, by simpa using neg_eq_zero.mp h11⟩
    · have hzero :
          Wedge2Formalization.N7WeightedTwoPoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
        have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
          exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hdet)
        calc
          Wedge2Formalization.N7WeightedTwoPoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
              =
            (-(M 1 0)) • Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₁ (k := k)) g +
              M 0 0 • Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g := by
                simp [Wedge2Formalization.N7WeightedTwoPoint.ActBivector, Matrix.mul_add,
                  Matrix.add_mul]
          _ =
            (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [h1, h2]
          _ = 0 := by
                ext i j
                simp [Matrix.add_apply, Matrix.smul_apply]
                have hdet'' :
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j) =
                      M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by
                  calc
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j)
                        = (M 0 0 * M 1 1) * rep₂ (k := k) i j := by ring
                    _ = (M 0 1 * M 1 0) * rep₂ (k := k) i j := by rw [hdet']
                    _ = M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by ring
                rw [hdet'']
                ring_nf
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k)
          (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
          (g := g)
          hg).1 hzero
      rcases rep_pair_independent (k := k) horig with ⟨h10, h00⟩
      exact hcol0 ⟨h00, by simpa using neg_eq_zero.mp h10⟩
  have hsupp2 :
      M 1 0 = 0 ∨ M 1 0 + M 1 1 = 0 := by
    have h2' :
        Wedge2Formalization.N6WeightedTwoPoint.ActBivector
            (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)) g.toBlocks₂₂ =
          M 1 0 • (Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)) +
            M 1 1 • (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)) := by
      simpa [rep₁, rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₁,
        Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂,
        act_embedded_toBlocks₂₂_weighted] using congrArg Matrix.toBlocks₂₂ h2
    have hrep2_zero :
        Matrix.det (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)) = 0 := by
      simpa using
        ((Wedge2Formalization.Paper.N6.Row4.det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2
          (Or.inl rfl))
    have hdet2 :
        Matrix.det
            (Wedge2Formalization.N6WeightedTwoPoint.ActBivector
              (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)) g.toBlocks₂₂) = 0 := by
      rw [Wedge2Formalization.N6WeightedTwoPoint.ActBivector, Matrix.det_mul, Matrix.det_mul,
        Matrix.det_transpose, hrep2_zero]
      simp
    exact
      (Wedge2Formalization.Paper.N6.Row4.det_zero_iff
        (k := k) (a := M 1 0) (b := M 1 1)).1 (by simpa [h2'] using hdet2)
  have hγ : M 1 0 = 0 := by
    rcases hsupp2 with hγ0 | hsum
    · exact hγ0
    · by_contra hγ0
      have hrewrite :
          M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) =
            M 1 0 • (rep₁ (k := k) - rep₂ (k := k)) := by
        have hδeq : M 1 1 = -M 1 0 := by
          apply eq_neg_of_add_eq_zero_left
          simpa [add_comm] using hsum
        ext i j
        simp [hδeq, sub_eq_add_neg]
        ring
      have hlarge :
          4 ≤ (Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g).rank := by
        rw [h2, hrewrite]
        exact topCombo_rank_ge_four (k := k) (M 1 0) hγ0
      have hsmall :
          (Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g).rank ≤ 2 := by
        rw [act_rank_eq (k := k) (Ω := rep₂ (k := k)) (g := g) hg]
        exact rep₂_rank_le_two (k := k)
      omega
  have hδ : M 1 1 ≠ 0 := by
    intro hδ0
    have hzero :
        Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g = 0 := by
      simpa [hγ, hδ0] using h2
    have horig :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) (Ω := rep₂ (k := k)) (g := g) hg).1 hzero
    exact rep₂_ne_zero (k := k) horig
  have htop_action :
      Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g =
        M 0 0 • rep₁ (k := k) + (M 0 1 - M 1 1) • rep₂ (k := k) := by
    calc
      Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g
          =
        Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) g := by
            simp [sub_eq_add_neg]
      _ =
        Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₁ (k := k)) g +
          (-1 : k) • Wedge2Formalization.N7WeightedTwoPoint.ActBivector (rep₂ (k := k)) g := by
            simp [Wedge2Formalization.N7WeightedTwoPoint.ActBivector, Matrix.mul_add, Matrix.add_mul]
      _ = (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
            (-1 : k) • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
              rw [h1, h2]
      _ = M 0 0 • rep₁ (k := k) + (M 0 1 - M 1 1) • rep₂ (k := k) := by
            ext i j
            simp [hγ, sub_eq_add_neg]
            ring
  have hα : M 0 0 ≠ 0 := by
    intro hα0
    have hlarge :
        4 ≤ (Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g).rank := by
      rw [act_rank_eq (k := k) (Ω := rep₁ (k := k) - rep₂ (k := k)) (g := g) hg]
      simpa using topCombo_rank_ge_four (k := k) (1 : k) one_ne_zero
    have hsmall :
        (Wedge2Formalization.N7WeightedTwoPoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g).rank ≤ 2 := by
      rw [htop_action, hα0]
      simpa using combo_rank_le_two_of_leftCoeff_zero (k := k) (M 0 1 - M 1 1)
    omega
  have hrep12_zero :
      Matrix.det ((Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)) -
        (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k))) = 0 := by
    simpa [sub_eq_add_neg] using
      ((Wedge2Formalization.Paper.N6.Row4.det_zero_iff
        (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (by ring)))
  have htop_action₂₂ :
      Wedge2Formalization.N6WeightedTwoPoint.ActBivector
          ((Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)) -
            (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)))
          g.toBlocks₂₂ =
        M 0 0 • (Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)) +
          (M 0 1 - M 1 1) • (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k)) := by
    have htop_blocks := congrArg Matrix.toBlocks₂₂ htop_action
    have hsub :
        rep₁ (k := k) - rep₂ (k := k) =
          Matrix.fromBlocks 0 0 0
            ((Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)) -
              (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k))) := by
      ext i j
      cases i <;> cases j <;>
        simp [rep₁, rep₂, sub_eq_add_neg, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₁,
          Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂, Matrix.fromBlocks,
          Matrix.add_apply, Matrix.smul_apply]
    rw [hsub] at htop_blocks
    simpa [rep₁, rep₂, Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₁,
      Wedge2Formalization.N7WeightedTwoPoint.radWeightedRep₂, hsub,
      act_embedded_toBlocks₂₂_weighted, Matrix.toBlocks_fromBlocks₂₂] using htop_blocks
  have halpha_plus : M 0 0 + M 0 1 = M 1 1 := by
    have hcombo_zero :
        Matrix.det
          (M 0 0 • (Wedge2Formalization.N6WeightedTwoPoint.rep₁ (k := k)) +
            (M 0 1 - M 1 1) • (Wedge2Formalization.N6WeightedTwoPoint.rep₂ (k := k))) = 0 := by
      rw [← htop_action₂₂, Wedge2Formalization.N6WeightedTwoPoint.ActBivector, Matrix.det_mul,
        Matrix.det_mul, Matrix.det_transpose, hrep12_zero]
      simp
    rcases (Wedge2Formalization.Paper.N6.Row4.det_zero_iff
      (k := k) (a := M 0 0) (b := M 0 1 - M 1 1)).1 hcombo_zero with hα0 | hsum0
    · exact False.elim (hα hα0)
    · calc
        M 0 0 + M 0 1 = M 0 0 + (M 0 1 - M 1 1) + M 1 1 := by ring
        _ = 0 + M 1 1 := by rw [hsum0]
        _ = M 1 1 := by ring
  refine ⟨M, ?_, hM⟩
  constructor
  · simpa [Wedge2Formalization.Paper.N6.Row4.Qproj, hγ]
  · constructor
    · simpa [hγ] using halpha_plus
    · exact hdetM

end Row15

/-! Appendix A, `n = 7`, row 16.
Representative `S_1 + J_{a,1} + J_{b,1} + J_{c,1}`.
Divisor `[a] + [b] + [c]`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times (\SL_2(k) \times \SL_2(k) \times
\SL_2(k)))`, exact quotient family `Q_L =` the preimage of `S_3` in `GL_2(k)`.
-/
namespace Row16

abbrev I := Wedge2Formalization.N7ThreePoint.I
abbrev W := Wedge2Formalization.N7ThreePoint.W
abbrev V := Wedge2Formalization.N7ThreePoint.V

def rep₁ : Matrix V V k :=
  Wedge2Formalization.N7ThreePoint.radThreePointRep₁ (k := k)

def rep₂ : Matrix V V k :=
  Wedge2Formalization.N7ThreePoint.radThreePointRep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 16. -/
def paperRep₁ : Matrix V V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 16. -/
def paperRep₂ : Matrix V V k :=
  rep₂ (k := k)

/-- The row-16 paper representative already agrees with the internal working one. -/
def paperChange : Matrix V V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N7ThreePoint.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N7ThreePoint.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N7ThreePoint.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N7ThreePoint.ActBivector]

def Levi
    (a0 : Matrix I I k)
    (C0 : Matrix W I k)
    (A : Matrix Wedge2Formalization.N6ThreePoint.W
      Wedge2Formalization.N6ThreePoint.W k)
    (D : Matrix Wedge2Formalization.N6ThreePoint.I
      Wedge2Formalization.N6ThreePoint.I k) :
    Matrix V V k :=
  Matrix.fromBlocks a0 0 C0 (Matrix.fromBlocks A 0 0 D)

/-- Public name for the embedded `\SL_2(k) \times \SL_2(k)` factor in
Appendix A, row 16. -/
def SplitCore :
    Set
      (Matrix Wedge2Formalization.N6ThreePoint.W
        Wedge2Formalization.N6ThreePoint.W k) :=
  Wedge2Formalization.Paper.N6.Row5.SplitCore (k := k)

def K : Set (Matrix V V k) :=
  { g |
      ∃ a0 : Matrix I I k,
        ∃ C0 : Matrix W I k,
          ∃ A : Matrix Wedge2Formalization.N6ThreePoint.W
              Wedge2Formalization.N6ThreePoint.W k,
            ∃ D : Matrix Wedge2Formalization.N6ThreePoint.I
                Wedge2Formalization.N6ThreePoint.I k,
              A ∈ SplitCore (k := k) ∧
                D.det = 1 ∧
                g = Levi (k := k) a0 C0 A D }

def supportCondition (a b : k) : Prop :=
  Wedge2Formalization.Paper.N6.Row5.supportCondition (k := k) a b

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N6.Row5.Qproj (k := k)

def swapCoeff : Matrix (Fin 2) (Fin 2) k :=
  Wedge2Formalization.Paper.N6.Row5.swapCoeff (k := k)

def swapLift : Matrix V V k :=
  let E : Matrix Wedge2Formalization.N6ThreePoint.I
      Wedge2Formalization.N6ThreePoint.I k := !![(1 : k), 0; 0, (-1 : k)]
  let g6 : Matrix W W k :=
    Matrix.fromBlocks (Wedge2Formalization.N4.splitSwap (k := k)) 0 0 E
  Matrix.fromBlocks
    (Wedge2Formalization.N7WeightedTwoPoint.scalarBlock (k := k) (1 : k))
    0
    0
    g6

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N7ThreePoint.FixesRadThreePointPairBivector (k := k))
      (K (k := k)) := by
  intro g
  have hgBlocks :
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
    simpa using (Matrix.fromBlocks_toBlocks g).symm
  constructor
  · intro hg
    have hfixBlocks :
        Wedge2Formalization.N7ThreePoint.FixesRadThreePointPairBivector (k := k)
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) := by
      rw [← hgBlocks]
      exact hg
    rcases
      (Wedge2Formalization.N7Summary.radThreePoint_pointwise_bivector_iff
        (k := k)
        (a := g.toBlocks₁₁)
        (R := g.toBlocks₁₂)
        (C0 := g.toBlocks₂₁)
        (A := g.toBlocks₂₂.toBlocks₁₁)
        (B := g.toBlocks₂₂.toBlocks₁₂)
        (C := g.toBlocks₂₂.toBlocks₂₁)
        (D := g.toBlocks₂₂.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hfixBlocks) with
      ⟨hR, hB, hC, hA, hD⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₁, g.toBlocks₂₂.toBlocks₁₁, g.toBlocks₂₂.toBlocks₂₂, hA, hD, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 g.toBlocks₂₁
            (Matrix.fromBlocks g.toBlocks₂₂.toBlocks₁₁ 0 0 g.toBlocks₂₂.toBlocks₂₂) := by
            rw [(Matrix.fromBlocks_toBlocks g.toBlocks₂₂).symm]
            simp [hR, hB, hC]
      _ = Levi (k := k) g.toBlocks₁₁ g.toBlocks₂₁ g.toBlocks₂₂.toBlocks₁₁ g.toBlocks₂₂.toBlocks₂₂ := by
            rfl
  · rintro ⟨a0, C0, A, D, hA, hD, rfl⟩
    have h6 :
        Wedge2Formalization.N6ThreePoint.FixesPairBivector (Matrix.fromBlocks A 0 0 D) := by
      exact
        (Wedge2Formalization.N6Summary.threePoint_pointwise_bivector_iff
          (k := k) (A := A) (B := 0) (C := 0) (D := D)).2 ⟨rfl, rfl, hA, hD⟩
    exact
      (Wedge2Formalization.N7ThreePoint.fixesRadThreePointPair_fromBlocks_iff
        (k := k)
        (a := a0)
        (R := 0)
        (C0 := C0)
        (A := A)
        (B := 0)
        (C := 0)
        (D := D)).2 ⟨rfl, rfl, rfl, hA, hD⟩

theorem mem_K_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ a0 : Matrix I I k,
        ∃ C0 : Matrix W I k,
          ∃ A : Matrix Wedge2Formalization.N6ThreePoint.W
              Wedge2Formalization.N6ThreePoint.W k,
            ∃ D : Matrix Wedge2Formalization.N6ThreePoint.I
                Wedge2Formalization.N6ThreePoint.I k,
              A ∈ SplitCore (k := k) ∧
                D.det = 1 ∧
                g = Levi (k := k) a0 C0 A D := by
  rfl

/-- Table-facing kernel statement for Appendix A, row 16. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ a0 : Matrix I I k,
        ∃ C0 : Matrix W I k,
        ∃ A : Matrix Wedge2Formalization.N6ThreePoint.W
            Wedge2Formalization.N6ThreePoint.W k,
          ∃ D : Matrix Wedge2Formalization.N6ThreePoint.I
              Wedge2Formalization.N6ThreePoint.I k,
              A ∈ SplitCore (k := k) ∧
                D.det = 1 ∧
                g = Levi (k := k) a0 C0 A D }

/-- Table-facing kernel statement for Appendix A, row 16. -/
theorem mem_K_table_iff
    (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) :=
  mem_K_iff (k := k) (g := g)

theorem swap_action :
    ActsOnOrderedPair
      (Wedge2Formalization.N7ThreePoint.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (swapLift (k := k))
      (swapCoeff (k := k)) := by
  let g6 : Matrix W W k :=
    Matrix.fromBlocks (Wedge2Formalization.N4.splitSwap (k := k)) 0 0
      (!![(1 : k), 0; 0, (-1 : k)] : Matrix Wedge2Formalization.N6ThreePoint.I
        Wedge2Formalization.N6ThreePoint.I k)
  have hswap := Wedge2Formalization.N6Summary.threePoint_swap_lift_action (k := k)
  constructor
  ·
    calc
      Wedge2Formalization.N7ThreePoint.ActBivector (rep₁ (k := k)) (swapLift (k := k))
          = Matrix.fromBlocks 0 0 0
              (Wedge2Formalization.N6ThreePoint.ActBivector
                (Wedge2Formalization.N6ThreePoint.rep₁ (k := k)) g6) := by
                simpa [rep₁, swapLift, g6,
                  Wedge2Formalization.N7ThreePoint.radThreePointRep₁] using
                  (Wedge2Formalization.N7ThreePoint.act_embedded_fromBlocks_zeroUpperRight
                    (k := k)
                    (Ω := Wedge2Formalization.N6ThreePoint.rep₁ (k := k))
                    (a := Wedge2Formalization.N7WeightedTwoPoint.scalarBlock (k := k) (1 : k))
                    (C := 0)
                    (D := g6))
      _ = Matrix.fromBlocks 0 0 0 (Wedge2Formalization.N6ThreePoint.rep₁ (k := k)) := by
            simpa [g6] using hswap.1
      _ = Matrix.fromBlocks 0 0 0
            ((swapCoeff (k := k) 0 0) • Wedge2Formalization.N6ThreePoint.rep₁ (k := k) +
              (swapCoeff (k := k) 0 1) • Wedge2Formalization.N6ThreePoint.rep₂ (k := k)) := by
            simp [swapCoeff, Wedge2Formalization.Paper.N6.Row5.swapCoeff]
      _ = (swapCoeff (k := k) 0 0) • rep₁ (k := k) + (swapCoeff (k := k) 0 1) • rep₂ (k := k) := by
            simpa [rep₁, rep₂] using
              (Wedge2Formalization.N7ThreePoint.embed_linearCombination
                (k := k)
                (α := swapCoeff (k := k) 0 0)
                (β := swapCoeff (k := k) 0 1))
  ·
    calc
      Wedge2Formalization.N7ThreePoint.ActBivector (rep₂ (k := k)) (swapLift (k := k))
          = Matrix.fromBlocks 0 0 0
              (Wedge2Formalization.N6ThreePoint.ActBivector
                (Wedge2Formalization.N6ThreePoint.rep₂ (k := k)) g6) := by
                simpa [rep₂, swapLift, g6,
                  Wedge2Formalization.N7ThreePoint.radThreePointRep₂] using
                  (Wedge2Formalization.N7ThreePoint.act_embedded_fromBlocks_zeroUpperRight
                    (k := k)
                    (Ω := Wedge2Formalization.N6ThreePoint.rep₂ (k := k))
                    (a := Wedge2Formalization.N7WeightedTwoPoint.scalarBlock (k := k) (1 : k))
                    (C := 0)
                    (D := g6))
      _ = Matrix.fromBlocks 0 0 0
            (Wedge2Formalization.N6ThreePoint.rep₁ (k := k) + (-1 : k) •
              Wedge2Formalization.N6ThreePoint.rep₂ (k := k)) := by
            simpa [g6, sub_eq_add_neg] using hswap.2
      _ = Matrix.fromBlocks 0 0 0
            ((swapCoeff (k := k) 1 0) • Wedge2Formalization.N6ThreePoint.rep₁ (k := k) +
              (swapCoeff (k := k) 1 1) • Wedge2Formalization.N6ThreePoint.rep₂ (k := k)) := by
            simp [swapCoeff, Wedge2Formalization.Paper.N6.Row5.swapCoeff, sub_eq_add_neg]
      _ = (swapCoeff (k := k) 1 0) • rep₁ (k := k) + (swapCoeff (k := k) 1 1) • rep₂ (k := k) := by
            simpa [rep₁, rep₂] using
              (Wedge2Formalization.N7ThreePoint.embed_linearCombination
                (k := k)
                (α := swapCoeff (k := k) 1 0)
                (β := swapCoeff (k := k) 1 1))

/-- Rank drop on the three-point line is detected by the lower-right determinant. -/
theorem det_zero_iff
    (a b : k) :
    Matrix.det
        (Matrix.toBlocks₂₂ (a • rep₁ (k := k) + b • rep₂ (k := k))) = 0 ↔
      supportCondition (k := k) a b := by
  simpa [supportCondition] using
    Wedge2Formalization.N7Summary.radThreePoint_lowerRight_det_zero_iff
      (k := k) (a := a) (b := b)

private theorem act_embedded_toBlocks₂₂_threePoint
    (Ω : Matrix Wedge2Formalization.N7ThreePoint.W
      Wedge2Formalization.N7ThreePoint.W k)
    (g : Matrix V V k) :
    Matrix.toBlocks₂₂
        (Wedge2Formalization.N7ThreePoint.ActBivector
          (Matrix.fromBlocks 0 0 0 Ω) g) =
      Wedge2Formalization.N6ThreePoint.ActBivector Ω g.toBlocks₂₂ := by
  rw [← Matrix.fromBlocks_toBlocks g]
  simp [Wedge2Formalization.N7ThreePoint.ActBivector,
    Wedge2Formalization.N6ThreePoint.ActBivector, Matrix.fromBlocks_transpose,
    Matrix.fromBlocks_multiply]

private theorem rep_pair_independent
    {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have hbot := congrArg (fun M => M (Sum.inr (Sum.inr 0)) (Sum.inr (Sum.inr 1))) h
  have hb : b = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N7ThreePoint.radThreePointRep₁,
      Wedge2Formalization.N7ThreePoint.radThreePointRep₂, Wedge2Formalization.N6ThreePoint.rep₁,
      Wedge2Formalization.N6ThreePoint.rep₂, Wedge2Formalization.N4.splitRep₁,
      Wedge2Formalization.N4.splitRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J, Matrix.fromBlocks,
      Matrix.add_apply, Matrix.smul_apply] using hbot
  have htop :=
    congrArg
      (fun M =>
        M (Sum.inr (Sum.inl (Sum.inl 0))) (Sum.inr (Sum.inl (Sum.inl 1))))
      h
  have ha : a = 0 := by
    simpa [hb, rep₁, rep₂, Wedge2Formalization.N7ThreePoint.radThreePointRep₁,
      Wedge2Formalization.N7ThreePoint.radThreePointRep₂, Wedge2Formalization.N6ThreePoint.rep₁,
      Wedge2Formalization.N6ThreePoint.rep₂, Wedge2Formalization.N4.splitRep₁,
      Wedge2Formalization.N4.splitRep₂, Wedge2Formalization.N4.ω12,
      Wedge2Formalization.N4.ω34, Wedge2Formalization.N4.J, Matrix.fromBlocks,
      Matrix.add_apply, Matrix.smul_apply] using htop
  exact ⟨ha, hb⟩

theorem quotient_image
    (g : Matrix V V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N7ThreePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N7ThreePoint.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1 :
      Wedge2Formalization.N7ThreePoint.ActBivector (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.1
  have h2 :
      Wedge2Formalization.N7ThreePoint.ActBivector (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    simpa [ActsOnOrderedPair] using hM.2
  have h12 :
      Wedge2Formalization.N7ThreePoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g =
        (M 0 0 - M 1 0) • rep₁ (k := k) +
          (M 0 1 - M 1 1) • rep₂ (k := k) := by
    calc
      Wedge2Formalization.N7ThreePoint.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g
          =
        Wedge2Formalization.N7ThreePoint.ActBivector
          (rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) g := by
            simp [sub_eq_add_neg]
      _ =
        Wedge2Formalization.N7ThreePoint.ActBivector (rep₁ (k := k)) g +
          (-1 : k) • Wedge2Formalization.N7ThreePoint.ActBivector (rep₂ (k := k)) g := by
            simp [Wedge2Formalization.N7ThreePoint.ActBivector, Matrix.mul_add,
              Matrix.add_mul]
      _ =
        (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
          (-1 : k) • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
            rw [h1, h2]
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
            Wedge2Formalization.N7ThreePoint.ActBivector (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using hM.1
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := rep₁ (k := k))
            (g := g)
            hg).1 hzero
        have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
          simpa using horig
        exact one_ne_zero (rep_pair_independent (k := k) hcomb).1
      · have hzero :
            Wedge2Formalization.N7ThreePoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N7ThreePoint.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) • Wedge2Formalization.N7ThreePoint.ActBivector (rep₁ (k := k)) g +
                M 0 1 • Wedge2Formalization.N7ThreePoint.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N7ThreePoint.ActBivector, Matrix.mul_add,
                    Matrix.add_mul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1, h2]
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
        rcases rep_pair_independent (k := k) horig with ⟨h11, h01⟩
        exact hcol1 ⟨h01, by simpa using neg_eq_zero.mp h11⟩
    · have hzero :
          Wedge2Formalization.N7ThreePoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
        have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
          exact sub_eq_zero.mp (by simpa [Matrix.det_fin_two] using hdet)
        calc
          Wedge2Formalization.N7ThreePoint.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
              =
            (-(M 1 0)) • Wedge2Formalization.N7ThreePoint.ActBivector (rep₁ (k := k)) g +
              M 0 0 • Wedge2Formalization.N7ThreePoint.ActBivector (rep₂ (k := k)) g := by
                simp [Wedge2Formalization.N7ThreePoint.ActBivector, Matrix.mul_add,
                  Matrix.add_mul]
          _ =
            (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [h1, h2]
          _ = 0 := by
                ext i j
                simp [Matrix.add_apply, Matrix.smul_apply]
                have hdet'' :
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j) =
                      M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by
                  calc
                    M 0 0 * (M 1 1 * rep₂ (k := k) i j)
                        = (M 0 0 * M 1 1) * rep₂ (k := k) i j := by ring
                    _ = (M 0 1 * M 1 0) * rep₂ (k := k) i j := by rw [hdet']
                    _ = M 1 0 * (M 0 1 * rep₂ (k := k) i j) := by ring
                rw [hdet'']
                ring_nf
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k)
          (Ω := (-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k))
          (g := g)
          hg).1 hzero
      rcases rep_pair_independent (k := k) horig with ⟨h10, h00⟩
      exact hcol0 ⟨h00, by simpa using neg_eq_zero.mp h10⟩
  have hrep1_zero :
      Matrix.det (Wedge2Formalization.N6ThreePoint.rep₁ (k := k)) = 0 := by
    simpa [supportCondition] using
      ((det_zero_iff (k := k) (a := (1 : k)) (b := (0 : k))).2 (Or.inr (Or.inr rfl)))
  have hrep2_zero :
      Matrix.det (Wedge2Formalization.N6ThreePoint.rep₂ (k := k)) = 0 := by
    simpa [supportCondition] using
      ((det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl))
  have hrep12_zero :
      Matrix.det ((Wedge2Formalization.N6ThreePoint.rep₁ (k := k)) -
        (Wedge2Formalization.N6ThreePoint.rep₂ (k := k))) = 0 := by
    have h :
        Matrix.det ((1 : k) • rep₁ (k := k) + (-1 : k) • rep₂ (k := k)).toBlocks₂₂ = 0 := by
      simpa [supportCondition, sub_eq_add_neg] using
        ((det_zero_iff (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (Or.inl (by ring))))
    simpa [rep₁, rep₂, sub_eq_add_neg,
      Wedge2Formalization.N7ThreePoint.radThreePointRep₁,
      Wedge2Formalization.N7ThreePoint.radThreePointRep₂] using h
  have hsupp1 :
      supportCondition (k := k) (M 0 0) (M 0 1) := by
    have hdet1 :
        Matrix.det
            (Matrix.toBlocks₂₂
              (Wedge2Formalization.N7ThreePoint.ActBivector (rep₁ (k := k)) g)) = 0 := by
      rw [rep₁, Wedge2Formalization.N7ThreePoint.radThreePointRep₁,
        act_embedded_toBlocks₂₂_threePoint, Wedge2Formalization.N6ThreePoint.ActBivector,
        Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose, hrep1_zero]
      simp
    exact (det_zero_iff (k := k) (a := M 0 0) (b := M 0 1)).1 (by simpa [h1] using hdet1)
  have hsupp2 :
      supportCondition (k := k) (M 1 0) (M 1 1) := by
    have hdet2 :
        Matrix.det
            (Matrix.toBlocks₂₂
              (Wedge2Formalization.N7ThreePoint.ActBivector (rep₂ (k := k)) g)) = 0 := by
      rw [rep₂, Wedge2Formalization.N7ThreePoint.radThreePointRep₂,
        act_embedded_toBlocks₂₂_threePoint, Wedge2Formalization.N6ThreePoint.ActBivector,
        Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose, hrep2_zero]
      simp
    exact (det_zero_iff (k := k) (a := M 1 0) (b := M 1 1)).1 (by simpa [h2] using hdet2)
  have hsupp12 :
      supportCondition (k := k) (M 0 0 - M 1 0) (M 0 1 - M 1 1) := by
    have hdet12 :
        Matrix.det
            (Matrix.toBlocks₂₂
              (Wedge2Formalization.N7ThreePoint.ActBivector
                (rep₁ (k := k) - rep₂ (k := k)) g)) = 0 := by
      have hsub :
          rep₁ (k := k) - rep₂ (k := k) =
            Matrix.fromBlocks 0 0 0
              ((Wedge2Formalization.N6ThreePoint.rep₁ (k := k)) -
                (Wedge2Formalization.N6ThreePoint.rep₂ (k := k))) := by
        ext i j
        cases i <;> cases j <;>
          simp [rep₁, rep₂, sub_eq_add_neg, Wedge2Formalization.N7ThreePoint.radThreePointRep₁,
            Wedge2Formalization.N7ThreePoint.radThreePointRep₂, Matrix.fromBlocks,
            Matrix.add_apply, Matrix.smul_apply]
      rw [hsub, act_embedded_toBlocks₂₂_threePoint, Wedge2Formalization.N6ThreePoint.ActBivector,
        Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose, hrep12_zero]
      simp
    exact (det_zero_iff
      (k := k) (a := M 0 0 - M 1 0) (b := M 0 1 - M 1 1)).1 (by simpa [h12] using hdet12)
  exact ⟨M, ⟨hdetM, hsupp1, hsupp2, hsupp12⟩, hM⟩

end Row16

end N7
end Paper
end Wedge2Formalization
