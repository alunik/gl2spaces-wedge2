import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N6.Row3

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 14.
Representative `S_1 + J_{a,2} + J_{b,1}`.
Divisor `2[a]+[b]`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times (\SL_2(k) \times \SL_2(k)))`, exact quotient
family `Q_L = T`.

This row is a radical-1 extension of N6 Row 3 (OnePointPlusSimple).
-/
namespace Row14

abbrev I := Fin 1
abbrev W := Wedge2Formalization.N6OnePointPlusSimple.V
abbrev V := Sum I W

def scalarBlock (a : k) : Matrix I I k := !![a]

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k))

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k))

def paperRep₁ : Matrix V V k := rep₁ (k := k)
def paperRep₂ : Matrix V V k := rep₂ (k := k)
def paperChange : Matrix V V k := 1

theorem paperRep₁_transport :
    (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (paperRep₁ (k := k)) (paperChange (k := k)) = rep₁ (k := k) := by
  simp [paperRep₁, paperChange]

theorem paperRep₂_transport :
    (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (paperRep₂ (k := k)) (paperChange (k := k)) = rep₂ (k := k) := by
  simp [paperRep₂, paperChange]

/-- The reps are linearly independent (from the W-block N6.Row3 independence). -/
theorem rep_pair_independent {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) : a = 0 ∧ b = 0 := by
  have h22 := congrArg toBlocks₂₂ h
  simp only [rep₁, rep₂, fromBlocks_smul, fromBlocks_add, toBlocks₂₂, fromBlocks,
    Matrix.of_apply, Sum.elim_inr, Matrix.zero_apply, Matrix.smul_apply,
    Matrix.add_apply, smul_zero, add_zero] at h22
  exact Wedge2Formalization.Paper.N6.Row3.rep_pair_independent (k := k) h22

/-- Public name for the embedded Appendix A, `n = 6`, row 3 core
`U_3 \rtimes (\SL_2(k) \times \SL_2(k))` on the lower-right `6`-space. -/
def OnePointPlusSimpleCore : Set (Matrix W W k) :=
  Wedge2Formalization.Paper.N6.Row3.K (k := k)

def K : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ OnePointPlusSimpleCore (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N6.Row3.Qproj (k := k)

def lift (a : k) : Matrix V V k :=
  Matrix.fromBlocks (1 : Matrix I I k) 0 0
    (Wedge2Formalization.Paper.N6.Row3.lift (k := k) a)

private theorem n6_rep₁_det_ne_zero :
    Matrix.det (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k)) ≠ 0 := by
  intro hdet
  rcases (Wedge2Formalization.Paper.N6.Row3.det_zero_iff
    (k := k) (a := 1) (b := 0)).1 (by simpa using hdet) with h | h
  · exact one_ne_zero h
  · simp at h

/-- The exact pointwise kernel family for the radical-1 extension. -/
theorem pointwise_stabilizer :
    ∀ g : Matrix V V k,
      (g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
        g * rep₂ (k := k) * gᵀ = rep₂ (k := k)) ↔
      g ∈ K (k := k) := by
  intro g
  set A0 := g.toBlocks₁₁
  set B := g.toBlocks₁₂
  set C := g.toBlocks₂₁
  set D := g.toBlocks₂₂
  have hgEq : g = Matrix.fromBlocks A0 B C D := (Matrix.fromBlocks_toBlocks g).symm
  constructor
  · -- Forward: fixing ⟹ membership in K
    rintro ⟨h1, h2⟩
    -- Compute the (1,2) block of g * rep₁ * gᵀ
    have h1_blocks : fromBlocks A0 B C D * fromBlocks 0 0 0
          (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k)) *
          (fromBlocks A0 B C D)ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k)) := by
      rw [← hgEq]; exact h1
    -- Extract (1,2) block: B * Ω₁ * Dᵀ = 0
    have htop : B * Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) * Dᵀ = 0 := by
      have := congrArg Matrix.toBlocks₁₂ h1_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- Extract (2,2) block: D * Ω₁ * Dᵀ = Ω₁
    have hbot : D * Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) * Dᵀ =
        Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) := by
      have := congrArg Matrix.toBlocks₂₂ h1_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- det(D) ≠ 0 from D * Ω * Dᵀ = Ω and det(Ω) ≠ 0
    have hdetD : Matrix.det D ≠ 0 := by
      intro hD
      have := congrArg Matrix.det hbot
      rw [Matrix.det_mul, Matrix.det_mul, Matrix.det_transpose, hD,
        zero_mul, zero_mul] at this
      exact n6_rep₁_det_ne_zero (k := k) this.symm
    -- B = 0 from B * Ω * Dᵀ = 0 with Ω * Dᵀ invertible
    have hB0 : B = 0 := by
      have hΩD_det : Matrix.det
          (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) * Dᵀ) ≠ 0 := by
        rw [Matrix.det_mul, Matrix.det_transpose]
        exact mul_ne_zero (n6_rep₁_det_ne_zero (k := k)) hdetD
      have hunit := isUnit_iff_ne_zero.mpr hΩD_det
      calc B = B * (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) * Dᵀ) *
              (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) * Dᵀ)⁻¹ := by
                symm; simpa [Matrix.mul_assoc] using
                  Matrix.mul_nonsing_inv_cancel_right
                    (A := Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) * Dᵀ) (B := B) hunit
        _ = 0 := by rw [← Matrix.mul_assoc, htop]; simp
    -- D fixes the N6 pair
    have hDfix : Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector (k := k) D := by
      constructor
      · -- D * rep₁ * Dᵀ = rep₁
        have := congrArg Matrix.toBlocks₂₂ h1
        simpa [rep₁, hgEq, hB0, fromBlocks_multiply, fromBlocks_transpose,
          Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector,
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector] using this
      · -- D * rep₂ * Dᵀ = rep₂
        have := congrArg Matrix.toBlocks₂₂ h2
        simpa [rep₂, hgEq, hB0, fromBlocks_multiply, fromBlocks_transpose,
          Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector,
          Wedge2Formalization.N6OnePointPlusSimple.ActBivector] using this
    -- D ∈ N6.Row3.K
    have hDmem : D ∈ Wedge2Formalization.Paper.N6.Row3.K (k := k) :=
      (Wedge2Formalization.Paper.N6.Row3.pointwise_stabilizer (k := k) D).1 hDfix
    exact ⟨A0, C, D, hDmem, by rw [hgEq]; simp [hB0]⟩
  · -- Backward: membership in K ⟹ fixing
    rintro ⟨A0', C', D', hD', rfl⟩
    have hDfix : Wedge2Formalization.N6OnePointPlusSimple.FixesPairBivector (k := k) D' :=
      (Wedge2Formalization.Paper.N6.Row3.pointwise_stabilizer (k := k) D').2 hD'
    constructor
    · -- fromBlocks A0' 0 C' D' * rep₁ * (fromBlocks A0' 0 C' D')ᵀ = rep₁
      show fromBlocks A0' 0 C' D' * fromBlocks 0 0 0
          (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.1
    · show fromBlocks A0' 0 C' D' * fromBlocks 0 0 0
          (Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.2

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ Wedge2Formalization.Paper.N6.Row3.K (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D := Iff.rfl

/-- Public name for the explicit Appendix A row-14 kernel cell. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ a0 : k,
        ∃ A : Matrix Wedge2Formalization.N4.I Wedge2Formalization.N4.I k,
          ∃ x y z : k,
            ∃ E : Matrix Wedge2Formalization.N6OnePointPlusSimple.I
                Wedge2Formalization.N6OnePointPlusSimple.I k,
              ∃ C : Matrix W I k,
                A.det = 1 ∧ E.det = 1 ∧
                g = Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C
                    (Wedge2Formalization.Paper.N6.Row3.U (k := k) x y z *
                      Wedge2Formalization.Paper.N6.Row3.Levi (k := k) A E) }

theorem mem_K_table_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) := by
  constructor
  · rintro ⟨A0, C, D, hD, rfl⟩
    rcases (Wedge2Formalization.Paper.N6.Row3.mem_K_table_iff (k := k) D).1 hD with
      ⟨A, x, y, z, E, hA, hE, hshape⟩
    exact ⟨A0 0 0, A, x, y, z, E, C, hA, hE, by
      rw [hshape]; ext i j
      rcases i with i | i <;> rcases j with j | j
      · fin_cases i; fin_cases j; simp [scalarBlock, fromBlocks]
      · simp [fromBlocks]
      · simp [fromBlocks]
      · simp [fromBlocks]⟩
  · rintro ⟨a0, A, x, y, z, E, C, hA, hE, rfl⟩
    exact ⟨scalarBlock (k := k) a0, C, _, (Wedge2Formalization.Paper.N6.Row3.mem_K_table_iff
      (k := k) _).2 ⟨A, x, y, z, E, hA, hE, rfl⟩, rfl⟩

theorem det_zero_iff (a b : k) :
    Matrix.det (Matrix.toBlocks₂₂ (a • rep₁ (k := k) + b • rep₂ (k := k))) = 0 ↔
      a = 0 ∨ a + b = 0 := by
  show Matrix.det (a • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
      b • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k)) = 0 ↔ a = 0 ∨ a + b = 0
  exact Wedge2Formalization.Paper.N6.Row3.det_zero_iff (k := k) a b

private theorem embed_act
    (D Ω : Matrix W W k) :
    fromBlocks (1 : Matrix I I k) 0 0 D *
      fromBlocks (0 : Matrix I I k) 0 0 Ω *
      (fromBlocks (1 : Matrix I I k) 0 0 D)ᵀ =
    fromBlocks 0 0 0 (D * Ω * Dᵀ) := by
  simp [fromBlocks_multiply, fromBlocks_transpose]

private theorem embed_lincomb
    (α β : k) (Ω₁ Ω₂ : Matrix W W k) :
    α • fromBlocks (0 : Matrix I I k) 0 0 Ω₁ +
      β • fromBlocks (0 : Matrix I I k) 0 0 Ω₂ =
    fromBlocks 0 0 0 (α • Ω₁ + β • Ω₂) := by
  simp [fromBlocks_smul, fromBlocks_add]

theorem quotient_action (a : k) (ha : a ≠ 0) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) a)
      (!![a, a * a - a; 0, a * a]) := by
  have h6 := Wedge2Formalization.Paper.N6.Row3.quotient_action (k := k) a ha
  simp only [ActsOnOrderedPair, Wedge2Formalization.N6OnePointPlusSimple.ActBivector] at h6
  refine ⟨?_, ?_⟩
  · show lift (k := k) a * rep₁ (k := k) * (lift (k := k) a)ᵀ =
        a • rep₁ (k := k) + (a * a - a) • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h6.1
  · show lift (k := k) a * rep₂ (k := k) * (lift (k := k) a)ᵀ =
        (0 : k) • rep₁ (k := k) + (a * a) • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h6.2

/-- Every invertible setwise stabilizer induces a coefficient matrix in the projective quotient.
The proof follows the same determinant-pencil argument as N6.Row3.quotient_image:
1. det(M) ≠ 0 from linear independence of rep₁, rep₂ + g invertible
2. M 1 0 = 0 from W-block determinant pencil structure (det_zero_iff)
3. M 0 0 + M 0 1 = M 1 1 from polynomial identity over infinite k
Since rep = fromBlocks 0 0 0 Ω_W, the W-block pencil determinant is inherited
from the parent level (N6.Row3). -/
theorem quotient_image [Infinite k]
    (g : Matrix V V k) (hg : Matrix.det g ≠ 0)
    (hpres : PreservesSpanPair
        (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
        (rep₁ (k := k)) (rep₂ (k := k)) g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
        (rep₁ (k := k)) (rep₂ (k := k)) g M := by
  rcases hpres with ⟨M, hact1, hact2⟩
  refine ⟨M, ?_, hact1, hact2⟩
  simp only [Qproj, Wedge2Formalization.Paper.N6.Row3.Qproj, Set.mem_setOf_eq]
  -- Independence of rep₁, rep₂ (W-block delegates to parent)
  have rep_indep : ∀ a b : k, a • rep₁ (k := k) + b • rep₂ (k := k) = 0 → a = 0 ∧ b = 0 := by
    intro a b h
    have h22 := congrArg toBlocks₂₂ h
    simp only [rep₁, rep₂, fromBlocks_smul, fromBlocks_add, toBlocks₂₂, fromBlocks,
      Matrix.of_apply, Sum.elim_inr, Matrix.zero_apply, Matrix.smul_apply,
      Matrix.add_apply, smul_zero, add_zero] at h22
    exact Wedge2Formalization.Paper.N6.Row3.rep_pair_independent (k := k) h22
  -- det_zero_iff on W-block pencil: det(a • N6rep₁ + b • N6rep₂) = 0 ↔ a = 0 ∨ a + b = 0
  -- From hact2, the (2,2) block gives:
  --   toBlocks₂₂ (M 1 0 • rep₁ + M 1 1 • rep₂) = M 1 0 • N6rep₁ + M 1 1 • N6rep₂
  -- Argument for M 1 0 = 0:
  --   det(g * rep₂ * gᵀ) has zero W-block det (since rep₂ has zero W-block det at parent)
  --   ⟹ det_zero_iff gives M 1 0 = 0 ∨ M 1 0 + M 1 1 = 0
  --   Polynomial identity over infinite k eliminates the second case
  -- Argument for M 0 0 + M 0 1 = M 1 1:
  --   Similarly from rep₁ pencil structure
  -- Argument for det(M) ≠ 0:
  --   If det M = 0, both images are proportional, collapsing the span.
  --   Contradicts invertibility of g + rep independence.
  -- See N6.Row3.quotient_image for the full 100-line proof of the analogous statement.
  -- Unfold ActsOnOrderedPair
  simp only [ActsOnOrderedPair] at hact1 hact2
  -- Step 1: Extract W-block action. Since rep = fromBlocks 0 0 0 Ω_W,
  -- toBlocks₂₂ of the action gives D * Ω_W * Dᵀ = M-coeff combination.
  -- Extract pencil det = 0 from hact2: taking toBlocks₂₂ of both sides,
  -- LHS = D * N6rep₂ * Dᵀ, RHS = M 1 0 • N6rep₁ + M 1 1 • N6rep₂
  -- det(LHS) = det(D)² * det(N6rep₂) = 0 since det(N6rep₂) = 0
  have h22_act2 : Matrix.det (M 1 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
      M 1 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k)) = 0 := by
    -- Take det of toBlocks₂₂ of both sides of hact2.
    -- LHS: det(toBlocks₂₂(g * rep₂ * gᵀ)) = det(D * N6rep₂ * Dᵀ) = det(D)² * det(N6rep₂) = 0
    -- RHS: det(toBlocks₂₂(M10 • rep₁ + M11 • rep₂)) = det(M10 • N6rep₁ + M11 • N6rep₂)
    -- Equating gives the result.
    -- Compute det of toBlocks₂₂ of the LHS of hact2
    have key : toBlocks₂₂ (g * rep₂ (k := k) * gᵀ) =
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) * g.toBlocks₂₂ᵀ := by
      ext i j; simp [rep₂, toBlocks₂₂, fromBlocks, mul_apply, Fintype.sum_sum_type, transpose_apply]
    -- Compute det of toBlocks₂₂ of the RHS of hact2
    have rhs_eq : toBlocks₂₂ (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) =
        M 1 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
          M 1 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := by
      ext i j; simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, Matrix.smul_apply, Matrix.add_apply]
    -- Equate dets: det(LHS) = det(RHS) where LHS = D*Ω₂*Dᵀ, RHS = pencil
    have h1 : det (toBlocks₂₂ (g * rep₂ (k := k) * gᵀ)) = 0 := by
      rw [key, det_mul, det_mul, det_transpose,
        Wedge2Formalization.Paper.N6.Row3.rep₂_det_zero (k := k)]; simp
    rw [congrArg (fun X => det (toBlocks₂₂ X)) hact2, rhs_eq] at h1
    exact h1
  -- Step 2: det(D * N6rep₂ * Dᵀ) = det(D)² * det(N6rep₂) = 0
  have h10_cases :=
    (Wedge2Formalization.Paper.N6.Row3.det_zero_iff (k := k) (M 1 0) (M 1 1)).1 h22_act2
  -- Step 3: Similarly for rep₁ - rep₂ pencil (for the sum condition)
  -- det(D*(N6rep₁-N6rep₂)*Dᵀ) = 0 since det(N6rep₁-N6rep₂) = 0
  -- This gives: det((M00-M10)•N6rep₁ + (M01-M11)•N6rep₂) = 0
  -- Combined with h10_cases and det(M) ≠ 0, forces Borel conditions.
  -- det(M) ≠ 0: same argument as Row 3
  have hdetM : Matrix.det M ≠ 0 := by
    intro hdet0
    have hM_dep : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
      simpa [Matrix.det_fin_two] using hdet0
    have hcomb_zero : (M 1 1) • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k) = 0 :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
        calc g * ((M 1 1) • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k)) * gᵀ
            = (M 1 1) • (g * rep₁ (k := k) * gᵀ) +
              (-(M 0 1)) • (g * rep₂ (k := k) * gᵀ) := by
              simp [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul]
          _ = 0 := by rw [hact1, hact2]; ext i j
                      simp only [Matrix.add_apply, Matrix.smul_apply, Matrix.zero_apply, smul_eq_mul, neg_mul]
                      linear_combination (rep₁ (k := k) i j) * hM_dep)
    rcases rep_indep _ _ hcomb_zero with ⟨h11, h01neg⟩
    have h01 : M 0 1 = 0 := neg_eq_zero.mp h01neg
    have hcomb2 : (-(M 1 0)) • rep₁ (k := k) + (M 0 0) • rep₂ (k := k) = 0 :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
        calc g * ((-(M 1 0)) • rep₁ (k := k) + (M 0 0) • rep₂ (k := k)) * gᵀ
            = (-(M 1 0)) • (g * rep₁ (k := k) * gᵀ) +
              (M 0 0) • (g * rep₂ (k := k) * gᵀ) := by
              simp [Matrix.mul_add, Matrix.add_mul, Matrix.mul_smul, Matrix.smul_mul]
          _ = 0 := by rw [hact1, hact2]; ext i j
                      simp only [Matrix.add_apply, Matrix.smul_apply, Matrix.zero_apply, smul_eq_mul, neg_mul]
                      linear_combination (rep₂ (k := k) i j) * hM_dep)
    have h10 : M 1 0 = 0 := neg_eq_zero.mp (rep_indep _ _ hcomb2).1
    have hrep2_zero : rep₂ (k := k) = 0 := by
      have : g * rep₂ (k := k) * gᵀ = 0 := by
        have := hact2; rw [h10, h11] at this; simpa using this
      exact (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 this
    exact absurd ((rep_indep 0 1 (by simpa using hrep2_zero)).2) one_ne_zero
  -- W-block extraction helper
  have w_block_lhs (Ω : Matrix W W k) :
      toBlocks₂₂ (g * fromBlocks (0 : Matrix I I k) 0 0 Ω * gᵀ) =
      g.toBlocks₂₂ * Ω * g.toBlocks₂₂ᵀ := by
    ext i j; simp [toBlocks₂₂, fromBlocks, mul_apply, Fintype.sum_sum_type, transpose_apply]
  have w_block_rhs (a b : k) :
      toBlocks₂₂ (a • rep₁ (k := k) + b • rep₂ (k := k)) =
      a • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        b • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := by
    ext i j; simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, smul_apply, add_apply]
  -- Extract W-block ActBivector
  have hD_act (Ω : Matrix W W k) (a b : k)
      (hact : g * fromBlocks (0 : Matrix I I k) 0 0 Ω * gᵀ =
        a • rep₁ (k := k) + b • rep₂ (k := k)) :
      g.toBlocks₂₂ * Ω * g.toBlocks₂₂ᵀ =
      a • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        b • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := by
    rw [← w_block_lhs, ← w_block_rhs]; exact congrArg toBlocks₂₂ hact
  -- D preserves the N6 span with coefficient matrix M
  have hD_act1 := hD_act _ _ _ hact1
  have hD_act2 := hD_act _ _ _ hact2
  -- D satisfies ActBivector
  have hD1 : Wedge2Formalization.N6OnePointPlusSimple.ActBivector
      (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k)) g.toBlocks₂₂ =
    M 0 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
      M 0 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := by
    simp only [Wedge2Formalization.N6OnePointPlusSimple.ActBivector]; exact hD_act1
  have hD2 : Wedge2Formalization.N6OnePointPlusSimple.ActBivector
      (Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k)) g.toBlocks₂₂ =
    M 1 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
      M 1 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := by
    simp only [Wedge2Formalization.N6OnePointPlusSimple.ActBivector]; exact hD_act2
  -- det(D) ≠ 0: if det D = 0, det(D*Ω*Dᵀ) = 0 for all Ω,
  -- so M00•r1+M01•r2 and M10•r1+M11•r2 both have det 0.
  -- By det_zero_iff + det M ≠ 0, this gives contradiction.
  have hdetD : Matrix.det g.toBlocks₂₂ ≠ 0 := by
    intro hD0
    -- Both pencil values have det 0
    have p1 : det (M 0 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        M 0 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k)) = 0 := by
      rw [← hD_act1, det_mul, det_mul, det_transpose, hD0]; simp
    have p2 := h22_act2
    -- det_zero_iff on both
    have c1 := (Wedge2Formalization.Paper.N6.Row3.det_zero_iff (k := k) (M 0 0) (M 0 1)).1 p1
    -- c1 : M00 = 0 ∨ M00+M01 = 0
    -- h10_cases : M10 = 0 ∨ M10+M11 = 0
    -- det M = M00*M11 - M01*M10
    -- Enumerate all 4 cases and show det M = 0 in each
    rcases c1 with h00 | h0001 <;> rcases h10_cases with h10 | h1011
    · -- M00=0, M10=0: det = 0*M11 - M01*0 = 0
      exact hdetM (by simp [det_fin_two, h00, h10])
    · -- M00=0, M10+M11=0: evaluate pencil at t=1
      -- g*(rep₁+rep₂)*gᵀ = (M00+M10)•rep₁ + (M01+M11)•rep₂ = M10•rep₁ + (M01-M10)•rep₂
      have h_sum : g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ =
          (M 0 0 + M 1 0) • rep₁ (k := k) + (M 0 1 + M 1 1) • rep₂ (k := k) := by
        calc g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ
            = g * rep₁ (k := k) * gᵀ + g * rep₂ (k := k) * gᵀ := by
              simp [mul_add, add_mul]
          _ = _ := by rw [hact1, hact2]; ext i j; simp [add_apply, smul_apply]; ring
      -- (2,2) block det:
      have h22_sum : det (toBlocks₂₂ (g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ)) = 0 := by
        rw [show toBlocks₂₂ (g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ) =
          g.toBlocks₂₂ * (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
            Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k)) * g.toBlocks₂₂ᵀ from by
          ext i j; simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, mul_apply,
            Fintype.sum_sum_type, transpose_apply, add_apply]]
        simp [det_mul, det_transpose, hD0]
      rw [h_sum] at h22_sum
      have h22_rhs : toBlocks₂₂ ((M 0 0 + M 1 0) • rep₁ (k := k) + (M 0 1 + M 1 1) • rep₂ (k := k)) =
        (M 0 0 + M 1 0) • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
          (M 0 1 + M 1 1) • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := by
        ext i j; simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, smul_apply, add_apply]
      rw [h22_rhs] at h22_sum
      -- det_zero_iff on (M00+M10, M01+M11) = (M10, M01-M10):
      rcases (Wedge2Formalization.Paper.N6.Row3.det_zero_iff (k := k)
        (M 0 0 + M 1 0) (M 0 1 + M 1 1)).1 h22_sum with h_a | h_ab
      · -- M00+M10 = 0, with M00=0: M10 = 0. det M = 0.
        have hM10 : M 1 0 = 0 := by linear_combination h_a - h00
        exact hdetM (by simp [det_fin_two, h00, hM10])
      · -- (M00+M10)+(M01+M11) = 0. With M00=0, M10+M11=0: M01 = 0. det M = 0.
        have hM01 : M 0 1 = 0 := by linear_combination h_ab - h1011 - h00
        exact hdetM (by simp [det_fin_two, h00, hM01])
    · -- M00+M01=0, M10=0: evaluate sum pencil at t=1
      -- Same h_sum and h22_sum as above
      have h_sum : g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ =
          (M 0 0 + M 1 0) • rep₁ (k := k) + (M 0 1 + M 1 1) • rep₂ (k := k) := by
        calc g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ
            = g * rep₁ (k := k) * gᵀ + g * rep₂ (k := k) * gᵀ := by
              simp [mul_add, add_mul]
          _ = _ := by rw [hact1, hact2]; ext i j; simp [add_apply, smul_apply]; ring
      have h22_sum : det (toBlocks₂₂ (g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ)) = 0 := by
        rw [show toBlocks₂₂ (g * (rep₁ (k := k) + rep₂ (k := k)) * gᵀ) =
          g.toBlocks₂₂ * (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
            Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k)) * g.toBlocks₂₂ᵀ from by
          ext i j; simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, mul_apply,
            Fintype.sum_sum_type, transpose_apply, add_apply]]
        simp [det_mul, det_transpose, hD0]
      rw [h_sum] at h22_sum
      have h22_rhs : toBlocks₂₂ ((M 0 0 + M 1 0) • rep₁ (k := k) + (M 0 1 + M 1 1) • rep₂ (k := k)) =
        (M 0 0 + M 1 0) • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
          (M 0 1 + M 1 1) • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := by
        ext i j; simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, smul_apply, add_apply]
      rw [h22_rhs] at h22_sum
      rcases (Wedge2Formalization.Paper.N6.Row3.det_zero_iff (k := k)
        (M 0 0 + M 1 0) (M 0 1 + M 1 1)).1 h22_sum with h_a | h_ab
      · -- M00+M10 = 0, with M10=0: M00 = 0. Then M01 = -M00 = 0. det M = 0.
        have hM00 : M 0 0 = 0 := by linear_combination h_a - h10
        exact hdetM (by
          have hM01 : M 0 1 = 0 := by linear_combination h0001 - hM00
          simp [det_fin_two, hM00, h10])
      · -- (M00+M10)+(M01+M11) = 0. With M10=0: M00+M01+M11 = 0.
        -- With M00+M01=0: M11 = 0. det M = M00*0 - (-M00)*0 = 0.
        have hM11 : M 1 1 = 0 := by linear_combination h_ab - h0001 - h10
        exact hdetM (by rw [det_fin_two, h10, hM11]; ring)
    · -- M00+M01=0, M10+M11=0: det = M00*(-M10)-(-M00)*M10 = -M00*M10+M00*M10 = 0
      exact hdetM (by
        have h01 : M 0 1 = -(M 0 0) := by linear_combination h0001
        have h11 : M 1 1 = -(M 1 0) := by linear_combination h1011
        rw [det_fin_two, h01, h11]; ring)
  -- Now D has det ≠ 0 and PreservesSpanPair. Delegate to N6.Row3.quotient_image.
  have hD_pres : PreservesSpanPair
      (Wedge2Formalization.N6OnePointPlusSimple.ActBivector (k := k))
      (Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k))
      (Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k))
      g.toBlocks₂₂ := ⟨M, hD1, hD2⟩
  rcases Wedge2Formalization.Paper.N6.Row3.quotient_image (k := k)
    g.toBlocks₂₂ hdetD hD_pres with ⟨M', hM'mem, hM'act1, hM'act2⟩
  -- Uniqueness M' = M: from rep independence at the N6 level
  have n6_indep := @Wedge2Formalization.Paper.N6.Row3.rep_pair_independent k _
  have hM'_eq : M' = M := by
    -- hD1 and hM'act1 both give the action on N6.rep₁.
    -- Difference: (M'00-M00)•rep₁ + (M'01-M01)•rep₂ = 0.
    have eq1 : hD1.symm.trans hM'act1 = hD1.symm.trans hM'act1 := rfl
    -- Actually just use the equalities directly:
    -- hD1 : ActBivector rep₁ D = M00•rep₁ + M01•rep₂
    -- hM'act1 : ActBivector rep₁ D = M'00•rep₁ + M'01•rep₂
    -- So M00•rep₁ + M01•rep₂ = M'00•rep₁ + M'01•rep₂
    have heq1 : M 0 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        M 0 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) =
      M' 0 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        M' 0 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := hD1.symm.trans hM'act1
    have heq2 : M 1 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        M 1 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) =
      M' 1 0 • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        M' 1 1 • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) := hD2.symm.trans hM'act2
    -- Difference = 0
    have d1 : (M 0 0 - M' 0 0) • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        (M 0 1 - M' 0 1) • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) = 0 := by
      have := sub_eq_zero.mpr heq1
      rw [show _ - _ = (M 0 0 - M' 0 0) • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        (M 0 1 - M' 0 1) • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) from by
        simp [sub_smul]; abel] at this; exact this
    have d2 : (M 1 0 - M' 1 0) • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        (M 1 1 - M' 1 1) • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) = 0 := by
      have := sub_eq_zero.mpr heq2
      rw [show _ - _ = (M 1 0 - M' 1 0) • Wedge2Formalization.Paper.N6.Row3.rep₁ (k := k) +
        (M 1 1 - M' 1 1) • Wedge2Formalization.Paper.N6.Row3.rep₂ (k := k) from by
        simp [sub_smul]; abel] at this; exact this
    rcases n6_indep d1 with ⟨e00, e01⟩
    rcases n6_indep d2 with ⟨e10, e11⟩
    ext i j; fin_cases i <;> fin_cases j <;>
      [exact (sub_eq_zero.mp e00).symm; exact (sub_eq_zero.mp e01).symm;
       exact (sub_eq_zero.mp e10).symm; exact (sub_eq_zero.mp e11).symm]
  rw [← hM'_eq]
  exact ⟨hM'mem.1, hM'mem.2.1, hM'mem.2.2⟩

end Row14

end N7
end Paper
end Wedge2Formalization
