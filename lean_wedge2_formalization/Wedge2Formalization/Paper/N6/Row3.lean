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

/-! Appendix A, `n = 6`, row 3.
Representative `⟨e₁∧e₃ + e₂∧e₄ + e₅∧e₆, e₁∧e₄ + e₅∧e₆⟩`.
Divisor `2[a]+[b]`.
Claimed stabilizer:
`K_L = U_3 ⋊ (SL₂(k) × SL₂(k))`, exact projective quotient family
`Q_L = T`.
-/
namespace Row3

open Polynomial

def rep₁ :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  Wedge2Formalization.N6OnePointPlusSimple.rep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  Wedge2Formalization.N6OnePointPlusSimple.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 3. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 3. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  rep₂ (k := k)

/-- The row-3 paper representative already agrees with the internal working pair. -/
def paperChange :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6OnePointPlusSimple.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6OnePointPlusSimple.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6OnePointPlusSimple.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6OnePointPlusSimple.ActBivector]

/-- Exact pointwise kernel for the direct-sum `2[a]+[b]` row. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6OnePointPlusSimple.V
        Wedge2Formalization.N6OnePointPlusSimple.V k) :=
  { g |
      ∃ A : Matrix Wedge2Formalization.N6OnePointPlusSimple.W
          Wedge2Formalization.N6OnePointPlusSimple.W k,
        ∃ E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
            Wedge2Formalization.N6OnePointPlusSimple.I k,
          Wedge2Formalization.N4.FixesOnePointPairBivector (k := k) A ∧
            E.det = 1 ∧
            g = Matrix.fromBlocks A 0 0 E }

/-- The displayed unipotent family `U_3` on the top repeated-support block. -/
def U
    (x y z : k) :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z)
    0
    0
    (1 : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
      Wedge2Formalization.N6OnePointPlusSimple.I k)

/-- The displayed Levi family `SL₂(k) × SL₂(k)`. -/
def Levi
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
      Wedge2Formalization.N6OnePointPlusSimple.I k) :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A)
    0
    0
    E

/-- Exact coefficient-side torus family in the ordered basis `(rep₁, rep₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ M 0 0 + M 0 1 = M 1 1 ∧ Matrix.det M ≠ 0 }

/-- Chosen torus lift. -/
def lift
    (a : k) :
    Matrix Wedge2Formalization.N6OnePointPlusSimple.V
      Wedge2Formalization.N6OnePointPlusSimple.V k :=
  let B : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k :=
    !![(a * a - a) * (a⁻¹ * a⁻¹), 0; 0, 0]
  let H : Matrix Wedge2Formalization.N6OnePointPlusSimple.W
      Wedge2Formalization.N6OnePointPlusSimple.W k :=
    Wedge2Formalization.N4.onePointScale (k := k) a *
      Wedge2Formalization.N4.onePointUpperShear (k := k) B
  let E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
      Wedge2Formalization.N6OnePointPlusSimple.I k := !![a * a, 0; 0, 1]
  Matrix.fromBlocks H 0 0 E

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector (k := k))
      (K (k := k)) := by
  intro g
  constructor
  · intro hg
    rcases
      (Wedge2Formalization.N6Summary.onePointPlusSimple_pointwise_bivector_iff
        (k := k)
        (A := g.toBlocks₁₁)
        (B := g.toBlocks₁₂)
        (C := g.toBlocks₂₁)
        (D := g.toBlocks₂₂)).1
        (by simpa [Matrix.fromBlocks_toBlocks] using hg) with
      ⟨hB, hC, hA, hD⟩
    refine ⟨g.toBlocks₁₁, g.toBlocks₂₂, hA, hD, ?_⟩
    calc
      g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
            exact (Matrix.fromBlocks_toBlocks g).symm
      _ = Matrix.fromBlocks g.toBlocks₁₁ 0 0 g.toBlocks₂₂ := by simp [hB, hC]
  · rintro ⟨A, E, hA, hE, rfl⟩
    exact
      (Wedge2Formalization.N6Summary.onePointPlusSimple_pointwise_bivector_iff
        (k := k)
        (A := A)
        (B := 0)
        (C := 0)
        (D := E)).2 ⟨rfl, rfl, hA, hE⟩

theorem det_zero_iff
    (a b : k) :
    Matrix.det (a • rep₁ (k := k) + b • rep₂ (k := k)) = 0 ↔
      a = 0 ∨ a + b = 0 :=
  Wedge2Formalization.N6Summary.onePointPlusSimple_det_zero_iff
    (k := k) (a := a) (b := b)

private theorem block_diag_eq_U_mul_Levi
    (A : Matrix Wedge2Formalization.N6OnePointPlusSimple.W
      Wedge2Formalization.N6OnePointPlusSimple.W k)
    (E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
      Wedge2Formalization.N6OnePointPlusSimple.I k)
    (hA :
      Wedge2Formalization.N4.FixesOnePointPairBivector (k := k) A) :
    ∃ A0 : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
      ∃ x y z : k,
        A0.det = 1 ∧
        Matrix.fromBlocks A 0 0 E =
          U (k := k) x y z * Levi (k := k) A0 E := by
  have hAK :
      A ∈ Wedge2Formalization.Paper.N4.Row1.K (k := k) :=
    (Wedge2Formalization.Paper.N4.Row1.pointwise_stabilizer
      (k := k)
      (g := A)).1 hA
  rcases
    (Wedge2Formalization.Paper.N4.Row1.mem_K_iff
      (k := k)
      (g := A)).1 hAK with
    ⟨A0, x, y, z, hA0, hEq⟩
  refine ⟨A0, x, y, z, hA0, ?_⟩
  rw [hEq]
  simp [U, Levi, Matrix.fromBlocks_multiply]

theorem U_pointwise
    (x y z : k) :
    Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector (U (k := k) x y z) := by
  apply (pointwise_stabilizer (k := k) (g := U (k := k) x y z)).2
  refine ⟨Wedge2Formalization.Paper.N4.Row1.U (k := k) x y z,
    (1 : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
      Wedge2Formalization.N6OnePointPlusSimple.I k), ?_, ?_, rfl⟩
  · exact Wedge2Formalization.Paper.N4.Row1.U_pointwise (k := k) x y z
  · simp

theorem Levi_pointwise
    (A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k)
    (hA : A.det = 1)
    (E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
      Wedge2Formalization.N6OnePointPlusSimple.I k)
    (hE : E.det = 1) :
    Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector (Levi (k := k) A E) := by
  apply (pointwise_stabilizer (k := k) (g := Levi (k := k) A E)).2
  refine ⟨Wedge2Formalization.Paper.N4.Row1.Levi (k := k) A, E, ?_, hE, rfl⟩
  exact Wedge2Formalization.Paper.N4.Row1.Levi_pointwise (k := k) A hA

theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6OnePointPlusSimple.V
        Wedge2Formalization.N6OnePointPlusSimple.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        ∃ x y z : k,
          ∃ E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
              Wedge2Formalization.N6OnePointPlusSimple.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g = U (k := k) x y z * Levi (k := k) A E := by
  constructor
  · rintro ⟨A, E, hA, hE, rfl⟩
    rcases block_diag_eq_U_mul_Levi (k := k) A E hA with ⟨A0, x, y, z, hA0, hEq⟩
    exact ⟨A0, x, y, z, E, hA0, hE, hEq⟩
  · rintro ⟨A, x, y, z, E, hA, hE, rfl⟩
    have hfixU := U_pointwise (k := k) x y z
    have hfixL := Levi_pointwise (k := k) A hA E hE
    have h1 :
        Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (rep₁ (k := k))
            (U (k := k) x y z * Levi (k := k) A E) =
          rep₁ (k := k) := by
      calc
        Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (rep₁ (k := k))
            (U (k := k) x y z * Levi (k := k) A E)
            =
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (Wedge2Formalization.N6OnePointPlusSimple.ActBivector
              (rep₁ (k := k))
              (Levi (k := k) A E))
            (U (k := k) x y z) := by
              rw [Wedge2Formalization.N6OnePointPlusSimple.actBivector_mul]
        _ =
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (rep₁ (k := k))
            (U (k := k) x y z) := by
              simpa [Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector] using
                congrArg
                  (fun M =>
                    Wedge2Formalization.N6OnePointPlusSimple.ActBivector M
                      (U (k := k) x y z))
                  hfixL.1
        _ = rep₁ (k := k) := hfixU.1
    have h2 :
        Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (rep₂ (k := k))
            (U (k := k) x y z * Levi (k := k) A E) =
          rep₂ (k := k) := by
      calc
        Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (rep₂ (k := k))
            (U (k := k) x y z * Levi (k := k) A E)
            =
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (Wedge2Formalization.N6OnePointPlusSimple.ActBivector
              (rep₂ (k := k))
              (Levi (k := k) A E))
            (U (k := k) x y z) := by
              rw [Wedge2Formalization.N6OnePointPlusSimple.actBivector_mul]
        _ =
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (rep₂ (k := k))
            (U (k := k) x y z) := by
              simpa [Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector] using
                congrArg
                  (fun M =>
                    Wedge2Formalization.N6OnePointPlusSimple.ActBivector M
                      (U (k := k) x y z))
                  hfixL.2
        _ = rep₂ (k := k) := hfixU.2
    exact (pointwise_stabilizer (k := k) (g := U (k := k) x y z * Levi (k := k) A E)).1 ⟨h1, h2⟩

/-- Table-facing kernel statement for Appendix A, row 3:
`K_L = U_3 \rtimes (SL_2(k) \times SL_2(k))`. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6OnePointPlusSimple.V
        Wedge2Formalization.N6OnePointPlusSimple.V k) :
    g ∈ K (k := k) ↔
      ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
        ∃ x y z : k,
          ∃ E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
              Wedge2Formalization.N6OnePointPlusSimple.I k,
            A.det = 1 ∧
            E.det = 1 ∧
            g = U (k := k) x y z * Levi (k := k) A E :=
  mem_K_iff (k := k) (g := g)

theorem coeff_mem_Qproj
    (a : k)
    (ha : a ≠ 0) :
    !![a, a * a - a; 0, a * a] ∈ Qproj (k := k) := by
  constructor
  · simp
  · constructor
    · change a + (a * a - a) = a * a
      ring
    · simpa [Matrix.det_fin_two] using mul_ne_zero ha (mul_ne_zero ha ha)

theorem quotient_action
    (a : k)
    (ha : a ≠ 0) :
    ActsOnOrderedPair
      (Wedge2Formalization.N6OnePointPlusSimple.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      (lift (k := k) a)
      (!![a, a * a - a; 0, a * a]) := by
  constructor
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.onePointPlusSimple_torus_lift_action
        (k := k) (a := a) ha).1
  · simpa [ActsOnOrderedPair, rep₁, rep₂, lift] using
      (Wedge2Formalization.N6Summary.onePointPlusSimple_torus_lift_action
        (k := k) (a := a) ha).2

theorem rep_pair_independent
    {a b : k}
    (h :
      a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h01 := congrArg
    (fun M => M (Sum.inl (Sum.inl 0)) (Sum.inl (Sum.inl 1))) h
  have hb : b = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6OnePointPlusSimple.rep₁,
      Wedge2Formalization.N6OnePointPlusSimple.rep₂, Wedge2Formalization.N4.onePointRep₁,
      Wedge2Formalization.N4.onePointRep₂, Wedge2Formalization.N4.ω12, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h01
  have h45 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h
  have ha : a = 0 := by
    simpa [hb, rep₁, rep₂, Wedge2Formalization.N6OnePointPlusSimple.rep₁,
      Wedge2Formalization.N6OnePointPlusSimple.rep₂, Wedge2Formalization.N4.J,
      Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply] using h45
  exact ⟨ha, hb⟩

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  exact one_ne_zero (rep_pair_independent (k := k) hcomb).1

private theorem rep_pair_linearCombination_ne_zero
    {a b : k}
    (h : a ≠ 0 ∨ b ≠ 0) :
    a • rep₁ (k := k) + b • rep₂ (k := k) ≠ 0 := by
  intro hzero
  rcases rep_pair_independent (k := k) hzero with ⟨ha, hb⟩
  exact h.elim (fun ha' => ha' ha) (fun hb' => hb' hb)

private theorem rep₁_det_ne_zero : Matrix.det (rep₁ (k := k)) ≠ 0 := by
  intro hdet
  rcases (det_zero_iff (k := k) (a := (1 : k)) (b := (0 : k))).1 (by simpa using hdet) with
    h1 | h1
  · exact one_ne_zero h1
  · have : (1 : k) = 0 := by simpa using h1
    exact one_ne_zero this

theorem rep₂_det_zero : Matrix.det (rep₂ (k := k)) = 0 := by
  simpa using
    ((det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl))

private theorem topRep_det_zero :
    Matrix.det (rep₁ (k := k) - rep₂ (k := k)) = 0 := by
  simpa [sub_eq_add_neg] using
    ((det_zero_iff (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (by ring)))

theorem quotient_image
    [Infinite k]
    (g :
      Matrix Wedge2Formalization.N6OnePointPlusSimple.V
        Wedge2Formalization.N6OnePointPlusSimple.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N6OnePointPlusSimple.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6OnePointPlusSimple.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  have h1' :
      Wedge2Formalization.N6OnePointPlusSimple.ActBivector
          (rep₁ (k := k)) g =
        M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k) := by
    exact hM.1
  have h2' :
      Wedge2Formalization.N6OnePointPlusSimple.ActBivector
          (rep₂ (k := k)) g =
        M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k) := by
    exact hM.2
  have hdetM : Matrix.det M ≠ 0 := by
    by_contra hdet
    by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
    · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
      · have hzero :
            Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                (rep₁ (k := k)) g = 0 := by
          simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using h1'
        have horig :=
          (actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := rep₁ (k := k)) (g := g) hg).1 hzero
        exact rep₁_ne_zero (k := k) horig
      · have hzero :
            Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
          calc
            Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
                =
              (-(M 1 1)) •
                  Wedge2Formalization.N6OnePointPlusSimple.ActBivector (rep₁ (k := k)) g +
                M 0 1 •
                  Wedge2Formalization.N6OnePointPlusSimple.ActBivector (rep₂ (k := k)) g := by
                    simp [Wedge2Formalization.N6OnePointPlusSimple.ActBivector, Matrix.mul_add,
                      Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
            _ =
              (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1', h2']
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
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
        have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
          have : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
            simpa [Matrix.det_fin_two] using hdet
          exact sub_eq_zero.mp this
        calc
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
              =
            (-(M 1 0)) •
                Wedge2Formalization.N6OnePointPlusSimple.ActBivector (rep₁ (k := k)) g +
              M 0 0 •
                Wedge2Formalization.N6OnePointPlusSimple.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N6OnePointPlusSimple.ActBivector, Matrix.mul_add,
                    Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
          _ =
            (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [h1', h2']
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
  have hrep₂det :
      Matrix.det
        (Wedge2Formalization.N6OnePointPlusSimple.ActBivector
          (rep₂ (k := k)) g) = 0 := by
    rw [Wedge2Formalization.N6OnePointPlusSimple.ActBivector, Matrix.det_mul, Matrix.det_mul,
      Matrix.det_transpose, rep₂_det_zero (k := k)]
    simp
  have h20_or :
      M 1 0 = 0 ∨ M 1 0 + M 1 1 = 0 :=
    (det_zero_iff (k := k) (a := M 1 0) (b := M 1 1)).1 (h2'.symm ▸ hrep₂det)
  have h20 : M 1 0 = 0 := by
    rcases h20_or with h20 | hsum
    · exact h20
    · by_contra h20
      have hrow0sum_ne : M 0 0 + M 0 1 ≠ 0 := by
        intro hsum0
        have hdet0 : Matrix.det M = 0 := by
          rw [Matrix.det_fin_two]
          have h11eq : M 1 1 = -M 1 0 :=
            eq_neg_of_add_eq_zero_left (by simpa [add_comm] using hsum)
          have h01eq : M 0 1 = -M 0 0 :=
            eq_neg_of_add_eq_zero_left (by simpa [add_comm] using hsum0)
          rw [h11eq, h01eq]
          ring
        exact hdetM hdet0
      let P : Polynomial k :=
        ((((Polynomial.C (M 0 0) * Polynomial.X + Polynomial.C (M 1 0)) ^ 2) *
          (Polynomial.C (M 0 0 + M 0 1) * Polynomial.X +
            Polynomial.C (M 1 0 + M 1 1))) ^ 2)
      let Q : Polynomial k :=
        Polynomial.C (Matrix.det g * Matrix.det g) *
          ((((Polynomial.X ^ 2) * (Polynomial.X + 1)) ^ 2))
      have hpoly_fun : ∀ t : k, Polynomial.eval t P = Polynomial.eval t Q := by
        intro t
        have hcomb :
            Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                (t • rep₁ (k := k) + rep₂ (k := k)) g =
              (t * M 0 0 + M 1 0) • rep₁ (k := k) +
                (t * M 0 1 + M 1 1) • rep₂ (k := k) := by
          calc
            Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                (t • rep₁ (k := k) + rep₂ (k := k)) g
                =
              t • Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                    (rep₁ (k := k)) g +
                Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                    (rep₂ (k := k)) g := by
                      simp [Wedge2Formalization.N6OnePointPlusSimple.ActBivector,
                        Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
            _ =
              t • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
                (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                  rw [h1', h2']
            _ =
              (t * M 0 0 + M 1 0) • rep₁ (k := k) +
                (t * M 0 1 + M 1 1) • rep₂ (k := k) := by
                  ext i j
                  simp [Matrix.add_apply, Matrix.smul_apply]
                  ring
        calc
          Polynomial.eval t P
              =
            (((t * M 0 0 + M 1 0) * (t * M 0 0 + M 1 0)) *
                ((t * M 0 0 + M 1 0) + (t * M 0 1 + M 1 1))) *
              (((t * M 0 0 + M 1 0) * (t * M 0 0 + M 1 0)) *
                ((t * M 0 0 + M 1 0) + (t * M 0 1 + M 1 1))) := by
                  simp [P]
                  ring
          _ =
            Matrix.det
              ((t * M 0 0 + M 1 0) • rep₁ (k := k) +
                (t * M 0 1 + M 1 1) • rep₂ (k := k)) := by
                  symm
                  simpa [rep₁, rep₂] using
                    (Wedge2Formalization.N6OnePointPlusSimple.det_in_basis
                      (k := k) (a := t * M 0 0 + M 1 0) (b := t * M 0 1 + M 1 1))
          _ =
            Matrix.det
              (Wedge2Formalization.N6OnePointPlusSimple.ActBivector
                (t • rep₁ (k := k) + rep₂ (k := k)) g) := by rw [hcomb]
          _ =
            Matrix.det g * Matrix.det g *
              Matrix.det (t • rep₁ (k := k) + rep₂ (k := k)) := by
                rw [Wedge2Formalization.N6OnePointPlusSimple.ActBivector, Matrix.det_mul,
                  Matrix.det_mul, Matrix.det_transpose]
                ring
          _ =
            Polynomial.eval t Q := by
                have hdet_basis_t :
                    Matrix.det (t • rep₁ (k := k) + rep₂ (k := k)) =
                      (((t * t) * (t + 1)) * ((t * t) * (t + 1))) := by
                  simpa [rep₁, rep₂] using
                    (Wedge2Formalization.N6OnePointPlusSimple.det_in_basis
                      (k := k) (a := t) (b := (1 : k)))
                calc
                  Matrix.det g * Matrix.det g * Matrix.det (t • rep₁ (k := k) + rep₂ (k := k))
                      =
                    Matrix.det g * Matrix.det g *
                      (((t * t) * (t + 1)) * ((t * t) * (t + 1))) := by
                        rw [hdet_basis_t]
                _ = Polynomial.eval t Q := by
                      simp [Q, hg]
                      ring
      have hpoly : P = Q := Polynomial.funext hpoly_fun
      have hcoeff2 := congrArg (fun f : Polynomial k => f.coeff 2) hpoly
      have hPcoeff2 :
          P.coeff 2 =
            M 1 0 * M 1 0 * (M 1 0 * M 1 0) *
              ((M 0 0 + M 0 1) * (M 0 0 + M 0 1)) := by
        let A : Polynomial k := Polynomial.C (M 0 0) * Polynomial.X + Polynomial.C (M 1 0)
        let c : k := M 0 0 + M 0 1
        have hcoeff_base :
            ((((A ^ 2) * (Polynomial.C c * Polynomial.X)) ^ 2 : Polynomial k).coeff 2) =
              M 1 0 * M 1 0 * (M 1 0 * M 1 0) * (c * c) := by
          have hPshape :
              (((A ^ 2) * (Polynomial.C c * Polynomial.X)) ^ 2 : Polynomial k) =
                (Polynomial.C (c * c) * A ^ 4) * Polynomial.X ^ 2 := by
            simp [pow_two, mul_assoc]
            ring
          rw [hPshape, Polynomial.coeff_mul_X_pow]
          simp [A, c, pow_two, Polynomial.coeff_zero_eq_eval_zero]
          ring
        dsimp [A, c] at hcoeff_base
        dsimp [P]
        simp [hsum]
        simpa [Polynomial.C_add] using hcoeff_base
      have hQcoeff2 : Q.coeff 2 = 0 := by
        have hQshape :
            Q =
              (Polynomial.C (Matrix.det g * Matrix.det g) * (Polynomial.X + 1) ^ 2) *
                Polynomial.X ^ 4 := by
          dsimp [Q]
          simp [pow_two, mul_assoc]
          ring_nf
          simp
        rw [hQshape, Polynomial.coeff_mul_X_pow']
        simp
      have hcoeff2' :
          M 1 0 * M 1 0 * (M 1 0 * M 1 0) *
            ((M 0 0 + M 0 1) * (M 0 0 + M 0 1)) = 0 := by
        have hcoeff2' : P.coeff 2 = Q.coeff 2 := by
          simpa using hcoeff2
        rw [hPcoeff2, hQcoeff2] at hcoeff2'
        exact hcoeff2'
      have hγ4_ne : M 1 0 * M 1 0 * (M 1 0 * M 1 0) ≠ 0 := by
        exact mul_ne_zero (mul_ne_zero h20 h20) (mul_ne_zero h20 h20)
      have hsum0sq : (M 0 0 + M 0 1) * (M 0 0 + M 0 1) = 0 := by
        exact (mul_eq_zero.mp hcoeff2').resolve_left hγ4_ne
      have hsum0 : M 0 0 + M 0 1 = 0 := by
        exact mul_self_eq_zero.mp (by simpa using hsum0sq)
      exact hrow0sum_ne hsum0
  have h11 : M 1 1 ≠ 0 := by
    intro h11
    have hdet0 : Matrix.det M = 0 := by
      rw [Matrix.det_fin_two, h20, h11]
      ring
    exact hdetM hdet0
  have htop :
      Wedge2Formalization.N6OnePointPlusSimple.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g =
        M 0 0 • rep₁ (k := k) + (M 0 1 - M 1 1) • rep₂ (k := k) := by
    calc
      Wedge2Formalization.N6OnePointPlusSimple.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g
          =
        Wedge2Formalization.N6OnePointPlusSimple.ActBivector
          (rep₁ (k := k) + (-1 : k) • rep₂ (k := k)) g := by
            simp [sub_eq_add_neg]
      _ =
        Wedge2Formalization.N6OnePointPlusSimple.ActBivector
            (rep₁ (k := k)) g +
          (-1 : k) •
            Wedge2Formalization.N6OnePointPlusSimple.ActBivector
              (rep₂ (k := k)) g := by
                simp [Wedge2Formalization.N6OnePointPlusSimple.ActBivector, Matrix.mul_add,
                  Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ =
        (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
          (-1 : k) • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
            rw [h1', h2']
      _ =
        M 0 0 • rep₁ (k := k) + (M 0 1 - M 1 1) • rep₂ (k := k) := by
          ext i j
          simp [h20, Matrix.add_apply, Matrix.smul_apply, sub_eq_add_neg]
          ring
  have htopdet :
      Matrix.det
        (Wedge2Formalization.N6OnePointPlusSimple.ActBivector
          (rep₁ (k := k) - rep₂ (k := k)) g) = 0 := by
    rw [Wedge2Formalization.N6OnePointPlusSimple.ActBivector, Matrix.det_mul, Matrix.det_mul,
      Matrix.det_transpose, topRep_det_zero (k := k)]
    simp
  have hrel_or :
      M 0 0 = 0 ∨ M 0 0 + (M 0 1 - M 1 1) = 0 :=
    (det_zero_iff (k := k) (a := M 0 0) (b := M 0 1 - M 1 1)).1 (htop.symm ▸ htopdet)
  have h00 : M 0 0 ≠ 0 := by
    intro h00zero
    have hdet0 : Matrix.det M = 0 := by
      rw [Matrix.det_fin_two, h20, h00zero]
      ring
    exact hdetM hdet0
  have hrel : M 0 0 + M 0 1 = M 1 1 := by
    rcases hrel_or with h00zero | hrelzero
    · exact False.elim (h00 h00zero)
    · have hrelzero' : (M 0 0 + M 0 1) - M 1 1 = 0 := by
          simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hrelzero
      exact sub_eq_zero.mp hrelzero'
  refine ⟨M, ?_, hM⟩
  constructor
  · simpa [h20]
  · constructor
    · simpa using hrel
    · simpa [Matrix.det_fin_two, h20] using hdetM

end Row3

end N6
end Paper
end Wedge2Formalization
