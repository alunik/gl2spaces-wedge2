import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N5
import Mathlib.Tactic.LinearCombination
import Mathlib.LinearAlgebra.Matrix.ToLinearEquiv

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 6.
Representative `S₂ + S₁² + J_{a,1}`.
Divisor `[a]`.
Claimed stabilizer:
`K_L = \Ga^{10} \rtimes (\Gm(k) \times (U_4 \rtimes (\Gm(k) \times \SL_2(k))))`,
exact quotient family `Q_L = B` (Borel).

This row is a radical-2 extension of N5 Row 3 (SimplePoint).
The Borel condition M 1 0 = 0 is witnessed by the quotient_image theorem.
-/
namespace Row6

abbrev I := Fin 2
abbrev W := Wedge2Formalization.N5SimplePoint.V
abbrev V := Sum I W

def scalarBlock (a b : k) : Matrix I I k := !![a, b; 0, a]

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k))

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k))

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

/-- The reps are linearly independent (from the W-block N5.Row3 independence). -/
theorem rep_pair_independent {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) : a = 0 ∧ b = 0 := by
  have ha : a = 0 := by
    have := congrArg (fun M => M (Sum.inr (Sum.inr 0)) (Sum.inr (Sum.inr 1))) h
    simp [rep₁, rep₂, Wedge2Formalization.Paper.N5.Row3.rep₁,
      Wedge2Formalization.Paper.N5.Row3.rep₂,
      Wedge2Formalization.N5SimplePoint.rep₁,
      Wedge2Formalization.N5SimplePoint.rep₂,
      Wedge2Formalization.N4.J, fromBlocks] at this
    exact this
  have hb : b = 0 := by
    have := congrArg (fun M => M (Sum.inr (Sum.inl 1)) (Sum.inr (Sum.inl 2))) h
    simp [ha, rep₁, rep₂, Wedge2Formalization.Paper.N5.Row3.rep₁,
      Wedge2Formalization.Paper.N5.Row3.rep₂,
      Wedge2Formalization.N5SimplePoint.rep₁,
      Wedge2Formalization.N5SimplePoint.rep₂,
      Wedge2Formalization.N3PureSingular.ω23, fromBlocks] at this
    exact this
  exact ⟨ha, hb⟩

/-- N5-level rep independence, derived from the N7-level independence by W-block extraction. -/
private theorem n5_rep_pair_independent {a b : k}
    (h : a • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
         b • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) = 0) : a = 0 ∧ b = 0 := by
  have h7 : a • rep₁ (k := k) + b • rep₂ (k := k) = 0 := by
    ext i j; rcases i with i | i <;> rcases j with j | j
    · simp [rep₁, rep₂, fromBlocks]
    · simp [rep₁, rep₂, fromBlocks]
    · simp [rep₁, rep₂, fromBlocks]
    · have := congrArg (fun M => M i j) h
      simp [rep₁, rep₂, fromBlocks, smul_apply, add_apply] at this ⊢
      exact this
  exact rep_pair_independent (k := k) h7

/-- Public name for the embedded Appendix A, `n = 5`, row 3 core
`U_4 \rtimes (G_m(k) \times \SL_2(k))` on the lower-right `5`-space. -/
def SimplePointCore : Set (Matrix W W k) :=
  Wedge2Formalization.Paper.N5.Row3.K (k := k)

/-- The kernel: parametric form with D in the embedded row-3 core. -/
def K : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ SimplePointCore (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D }

/-- The exact quotient family: Borel (upper triangular with nonzero determinant). -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

def lift (a b : k) : Matrix V V k :=
  Matrix.fromBlocks (1 : Matrix I I k) 0 0
    (Wedge2Formalization.Paper.N5.Row3.lift (k := k) a b)

/-- Any vector killed by both N5 reps is zero. This is the entry-level computation:
ker(ω₁₃) ∩ ker(ω₂₃) = {0} and J is invertible. -/
private theorem n5_mulVec_zero_of_reps (v : W → k)
    (h1 : Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) *ᵥ v = 0)
    (h2 : Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) *ᵥ v = 0) :
    v = 0 := by
  -- Convert mulVec to matrix product: use a diagonal-1 column vector approach.
  -- Instead, convert to matrix entry equations via congrFun + Pi.zero_apply.
  have entry1 : ∀ i, (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) *ᵥ v) i = 0 :=
    fun i => by rw [h1]; rfl
  have entry2 : ∀ i, (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) *ᵥ v) i = 0 :=
    fun i => by rw [h2]; rfl
  -- Unfold mulVec at specific indices to get v = 0.
  -- rep₁ = fromBlocks ω₁₃ 0 0 J where ω₁₃ = !![0,0,1;0,0,0;-1,0,0] and J = !![0,1;-1,0]
  -- rep₂ = fromBlocks ω₂₃ 0 0 0 where ω₂₃ = !![0,0,0;0,0,1;0,-1,0]
  ext (j | j)
  · fin_cases j
    · -- v(inl 0) = 0 from rep₁ at row (inl 2): -v(inl 0) = 0
      have := entry1 (Sum.inl 2)
      simp only [mulVec, dotProduct, Fintype.sum_sum_type,
        Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Fin.sum_univ_three, Fin.sum_univ_two, Pi.zero_apply] at this
      simp at this; exact this
    · -- v(inl 1) = 0 from rep₂ at row (inl 2): -v(inl 1) = 0
      have := entry2 (Sum.inl 2)
      simp only [mulVec, dotProduct, Fintype.sum_sum_type,
        Wedge2Formalization.Paper.N5.Row3.rep₂,
        Wedge2Formalization.N5SimplePoint.rep₂,
        Wedge2Formalization.N3PureSingular.ω23,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Fin.sum_univ_three, Fin.sum_univ_two, Pi.zero_apply] at this
      simp at this; exact this
    · -- v(inl 2) = 0 from rep₁ at row (inl 0): v(inl 2) = 0
      have := entry1 (Sum.inl 0)
      simp only [mulVec, dotProduct, Fintype.sum_sum_type,
        Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Fin.sum_univ_three, Fin.sum_univ_two, Pi.zero_apply] at this
      simp at this; exact this
  · fin_cases j
    · -- v(inr 0) = 0 from rep₁ at row (inr 1): -v(inr 0) = 0
      have := entry1 (Sum.inr 1)
      simp only [mulVec, dotProduct, Fintype.sum_sum_type,
        Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Fin.sum_univ_three, Fin.sum_univ_two, Pi.zero_apply] at this
      simp at this; exact this
    · -- v(inr 1) = 0 from rep₁ at row (inr 0): v(inr 1) = 0
      have := entry1 (Sum.inr 0)
      simp only [mulVec, dotProduct, Fintype.sum_sum_type,
        Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Fin.sum_univ_three, Fin.sum_univ_two, Pi.zero_apply] at this
      simp at this; exact this

/-- B = 0 from B * rep₁ = 0 ∧ B * rep₂ = 0 at the N5 level. -/
private theorem B_eq_zero_of_mul_reps_zero
    (B : Matrix I W k)
    (h1 : B * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) = 0)
    (h2 : B * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) = 0) :
    B = 0 := by
  ext i j
  fin_cases i <;> rcases j with (j | j)
  -- B(0, inl j):
  · fin_cases j
    · -- B(0, inl 0) = 0 from (B * rep₁)(0, inl 2) = 0
      have := congrArg (fun M => M 0 (Sum.inl 2)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
    · -- B(0, inl 1) = 0 from (B * rep₂)(0, inl 2) = 0
      have := congrArg (fun M => M 0 (Sum.inl 2)) h2
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₂,
        Wedge2Formalization.N5SimplePoint.rep₂,
        Wedge2Formalization.N3PureSingular.ω23,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
    · -- B(0, inl 2) = 0 from (B * rep₁)(0, inl 0) = 0
      have := congrArg (fun M => M 0 (Sum.inl 0)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
  -- B(0, inr j):
  · fin_cases j
    · -- B(0, inr 0) = 0 from (B * rep₁)(0, inr 1) = 0
      have := congrArg (fun M => M 0 (Sum.inr 1)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
    · -- B(0, inr 1) = 0 from (B * rep₁)(0, inr 0) = 0
      have := congrArg (fun M => M 0 (Sum.inr 0)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
  -- B(1, inl j):
  · fin_cases j
    · have := congrArg (fun M => M 1 (Sum.inl 2)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
    · have := congrArg (fun M => M 1 (Sum.inl 2)) h2
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₂,
        Wedge2Formalization.N5SimplePoint.rep₂,
        Wedge2Formalization.N3PureSingular.ω23,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
    · have := congrArg (fun M => M 1 (Sum.inl 0)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
  -- B(1, inr j):
  · fin_cases j
    · have := congrArg (fun M => M 1 (Sum.inr 1)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this
    · have := congrArg (fun M => M 1 (Sum.inr 0)) h1
      simp only [Wedge2Formalization.Paper.N5.Row3.rep₁,
        Wedge2Formalization.N5SimplePoint.rep₁,
        Wedge2Formalization.N3PureSingular.ω13, Wedge2Formalization.N4.J,
        mul_apply, Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two,
        fromBlocks, of_apply, Sum.elim_inl, Sum.elim_inr, zero_apply,
        cons_val', cons_val_zero, cons_val_one, head_cons, head_fin_const,
        Pi.zero_apply] at this
      simp at this; exact this

/-- det(D) ≠ 0 when D preserves the N5 rep span with an invertible coefficient matrix.
Uses the vector kernel argument: ker(Dᵀ) ⊆ ker(rep₁) ∩ ker(rep₂) = {0}. -/
private theorem det_ne_zero_of_preserves_span
    (D : Matrix W W k)
    (M : Matrix (Fin 2) (Fin 2) k) (hdetM : Matrix.det M ≠ 0)
    (hD1 : D * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) * Dᵀ =
      M 0 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        M 0 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k))
    (hD2 : D * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) * Dᵀ =
      M 1 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        M 1 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) :
    Matrix.det D ≠ 0 := by
  intro hD0
  -- From det D = 0, get a nonzero vector v with Dᵀ *ᵥ v = 0
  have hDt0 : Matrix.det Dᵀ = 0 := by rwa [Matrix.det_transpose]
  obtain ⟨v, hv_ne, hDtv⟩ := Matrix.exists_mulVec_eq_zero_iff.mpr hDt0
  -- From D * rep_i * Dᵀ = pencil_i, right-multiply both sides by v (as mulVec)
  -- LHS: (D * rep_i * Dᵀ) *ᵥ v = D *ᵥ (rep_i *ᵥ (Dᵀ *ᵥ v)) = D *ᵥ (rep_i *ᵥ 0) = 0
  have hpv1 : (M 0 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
      M 0 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) *ᵥ v = 0 := by
    rw [← hD1, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, hDtv,
      Matrix.mulVec_zero, Matrix.mulVec_zero]
  have hpv2 : (M 1 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
      M 1 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) *ᵥ v = 0 := by
    rw [← hD2, ← Matrix.mulVec_mulVec, ← Matrix.mulVec_mulVec, hDtv,
      Matrix.mulVec_zero, Matrix.mulVec_zero]
  -- Expand: M_i0 • (rep₁ *ᵥ v) + M_i1 • (rep₂ *ᵥ v) = 0
  -- Since det M ≠ 0, we can solve: rep₁ *ᵥ v = 0 and rep₂ *ᵥ v = 0.
  -- Use: det(M) • (rep₁ *ᵥ v) = M11 • (pencil₁ *ᵥ v) - M01 • (pencil₂ *ᵥ v) = 0
  -- Use linear algebra: from the two pencil equations and det M ≠ 0,
  -- derive that rep₁ *ᵥ v = 0 and rep₂ *ᵥ v = 0.
  -- Expand pencil equations into components
  have e1 := hpv1
  have e2 := hpv2
  rw [Matrix.add_mulVec, Matrix.smul_mulVec, Matrix.smul_mulVec] at e1 e2
  -- rep₁ *ᵥ v = 0
  have hsmul1 : Matrix.det M • (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) *ᵥ v) = 0 := by
    -- M11•e1 - M01•e2 = (M00*M11-M01*M10) • (rep₁ *ᵥ v) = det(M) • (rep₁ *ᵥ v)
    ext i
    have h1i := congrFun e1 i
    have h2i := congrFun e2 i
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at h1i h2i ⊢
    -- Goal: det M * (rep₁ *ᵥ v) i = 0
    -- h1i: M 0 0 * (rep₁ *ᵥ v) i + M 0 1 * (rep₂ *ᵥ v) i = 0
    -- h2i: M 1 0 * (rep₁ *ᵥ v) i + M 1 1 * (rep₂ *ᵥ v) i = 0
    -- linear_combination: M11 * h1i - M01 * h2i
    rw [show Matrix.det M = M 0 0 * M 1 1 - M 0 1 * M 1 0 by
      simp [Matrix.det_fin_two]]
    linear_combination M 1 1 * h1i - M 0 1 * h2i
  -- rep₂ *ᵥ v = 0
  have hsmul2 : Matrix.det M • (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) *ᵥ v) = 0 := by
    ext i
    have h1i := congrFun e1 i
    have h2i := congrFun e2 i
    simp only [Pi.add_apply, Pi.smul_apply, smul_eq_mul, Pi.zero_apply] at h1i h2i ⊢
    rw [show Matrix.det M = M 0 0 * M 1 1 - M 0 1 * M 1 0 by
      simp [Matrix.det_fin_two]]
    linear_combination -(M 1 0) * h1i + M 0 0 * h2i
  -- Since det M ≠ 0, cancel to get rep_i *ᵥ v = 0
  have hrv1 : Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) *ᵥ v = 0 := by
    rcases (smul_eq_zero.mp hsmul1) with h | h
    · exact absurd h hdetM
    · exact h
  have hrv2 : Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) *ᵥ v = 0 := by
    rcases (smul_eq_zero.mp hsmul2) with h | h
    · exact absurd h hdetM
    · exact h
  -- v = 0 from entry extraction
  exact hv_ne (n5_mulVec_zero_of_reps (k := k) v hrv1 hrv2)

/-- The exact pointwise kernel family for the radical-2 extension. -/
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
    have h1_blocks : fromBlocks A0 B C D * fromBlocks 0 0 0
          (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k)) *
          (fromBlocks A0 B C D)ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k)) := by
      rw [← hgEq]; exact h1
    have h2_blocks : fromBlocks A0 B C D * fromBlocks 0 0 0
          (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) *
          (fromBlocks A0 B C D)ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) := by
      rw [← hgEq]; exact h2
    -- Extract (1,2) block: B * Ω_i * Dᵀ = 0
    have htop1 : B * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) * Dᵀ = 0 := by
      have := congrArg Matrix.toBlocks₁₂ h1_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    have htop2 : B * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) * Dᵀ = 0 := by
      have := congrArg Matrix.toBlocks₁₂ h2_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- Extract (2,2) block: D * Ω_i * Dᵀ = Ω_i
    have hbot1 : D * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) * Dᵀ =
        Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) := by
      have := congrArg Matrix.toBlocks₂₂ h1_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    have hbot2 : D * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) * Dᵀ =
        Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) := by
      have := congrArg Matrix.toBlocks₂₂ h2_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- D fixes the N5 pair, so D ∈ N5.Row3.K
    have hDfix : Wedge2Formalization.N5SimplePoint.FixesPairBivector (k := k) D :=
      ⟨hbot1, hbot2⟩
    have hDmem : D ∈ Wedge2Formalization.Paper.N5.Row3.K (k := k) :=
      (Wedge2Formalization.Paper.N5.Row3.pointwise_stabilizer (k := k) D).1 hDfix
    -- det(D) ≠ 0: D fixes the pair with coefficient matrix I₂ (identity),
    -- which has det 1 ≠ 0.
    have hdetD : Matrix.det D ≠ 0 :=
      det_ne_zero_of_preserves_span (k := k) D (1 : Matrix (Fin 2) (Fin 2) k)
        (by simp [Matrix.det_fin_two])
        (by simp [hbot1])
        (by simp [hbot2])
    -- B = 0: from B * rep_i * Dᵀ = 0 with Dᵀ invertible → B * rep_i = 0
    have hB0 : B = 0 := by
      have hunit : IsUnit (Matrix.det Dᵀ) := by
        rw [Matrix.det_transpose]; exact isUnit_iff_ne_zero.mpr hdetD
      have h1_cancel : B * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) = 0 := by
        calc B * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k)
            = B * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) * Dᵀ * Dᵀ⁻¹ :=
              (Matrix.mul_nonsing_inv_cancel_right Dᵀ
                (B * Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k)) hunit).symm
          _ = 0 := by rw [htop1]; simp
      have h2_cancel : B * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) = 0 := by
        calc B * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)
            = B * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) * Dᵀ * Dᵀ⁻¹ :=
              (Matrix.mul_nonsing_inv_cancel_right Dᵀ
                (B * Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) hunit).symm
          _ = 0 := by rw [htop2]; simp
      exact B_eq_zero_of_mul_reps_zero (k := k) B h1_cancel h2_cancel
    exact ⟨A0, C, D, hDmem, by rw [hgEq]; simp [hB0]⟩
  · -- Backward: membership in K ⟹ fixing
    rintro ⟨A0', C', D', hD', rfl⟩
    have hDfix : Wedge2Formalization.N5SimplePoint.FixesPairBivector (k := k) D' :=
      (Wedge2Formalization.Paper.N5.Row3.pointwise_stabilizer (k := k) D').2 hD'
    constructor
    · show fromBlocks A0' 0 C' D' * fromBlocks 0 0 0
          (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.1
    · show fromBlocks A0' 0 C' D' * fromBlocks 0 0 0
          (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.2

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ Wedge2Formalization.Paper.N5.Row3.K (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D := Iff.rfl

/-- Public name for the explicit Appendix A row-6 kernel cell. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ a x y u v : k,
          ∃ E : Matrix Wedge2Formalization.N5SimplePoint.I
              Wedge2Formalization.N5SimplePoint.I k,
            ∃ C : Matrix W I k,
              a ≠ 0 ∧ E.det = 1 ∧
              g = Matrix.fromBlocks A0 0 C
                  (Matrix.fromBlocks
                    (Wedge2Formalization.N3PureSingular.pureSingularShape (k := k) a x y)
                    (!![(0 : k), 0; 0, 0; u, v])
                    (!![a * (u * E 0 1 - v * E 0 0), 0, 0;
                       a * (u * E 1 1 - v * E 1 0), 0, 0])
                    E) }

theorem mem_K_table_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) := by
  constructor
  · rintro ⟨A0, C, D, hD, rfl⟩
    rcases (Wedge2Formalization.Paper.N5.Row3.mem_K_table_iff (k := k) D).1 hD with
      ⟨a, x, y, u, v, E, ha, hE, hshape⟩
    exact ⟨A0, a, x, y, u, v, E, C, ha, hE, by rw [hshape]⟩
  · rintro ⟨A0, a, x, y, u, v, E, C, ha, hE, rfl⟩
    exact ⟨A0, C, _, (Wedge2Formalization.Paper.N5.Row3.mem_K_table_iff
      (k := k) _).2 ⟨a, x, y, u, v, E, ha, hE, rfl⟩, rfl⟩

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

theorem quotient_action (a b : k) (ha : a ≠ 0) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) a b)
      (!![a, b; 0, a * a]) := by
  have h5 := Wedge2Formalization.Paper.N5.Row3.quotient_action (k := k) a b ha
  simp only [ActsOnOrderedPair, Wedge2Formalization.N5SimplePoint.ActBivector] at h5
  have hcoeff : Wedge2Formalization.N4PaperSummary.row1_borelCoeff (k := k) a b = !![a, b; 0, a * a] := by
    rfl
  rw [hcoeff] at h5
  refine ⟨?_, ?_⟩
  · show lift (k := k) a b * rep₁ (k := k) * (lift (k := k) a b)ᵀ =
        a • rep₁ (k := k) + b • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h5.1
  · show lift (k := k) a b * rep₂ (k := k) * (lift (k := k) a b)ᵀ =
        (0 : k) • rep₁ (k := k) + (a * a) • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h5.2

/-- Every invertible setwise stabilizer induces a coefficient matrix in the Borel.
The proof delegates to N5.Row3.quotient_image via W-block extraction. -/
theorem quotient_image
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
  simp only [Qproj, Set.mem_setOf_eq]
  -- Unfold ActsOnOrderedPair
  simp only [ActsOnOrderedPair] at hact1 hact2
  -- Step 1: det(M) ≠ 0 (standard adjoint kernel argument)
  have rep_indep : ∀ a b : k, a • rep₁ (k := k) + b • rep₂ (k := k) = 0 → a = 0 ∧ b = 0 :=
    fun a b h => rep_pair_independent (k := k) h
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
  -- W-block extraction helpers
  have w_block_lhs (Ω : Matrix W W k) :
      toBlocks₂₂ (g * fromBlocks (0 : Matrix I I k) 0 0 Ω * gᵀ) =
      g.toBlocks₂₂ * Ω * g.toBlocks₂₂ᵀ := by
    ext i j; simp [toBlocks₂₂, fromBlocks, mul_apply, Fintype.sum_sum_type, transpose_apply]
  have w_block_rhs (a b : k) :
      toBlocks₂₂ (a • rep₁ (k := k) + b • rep₂ (k := k)) =
      a • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        b • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) := by
    ext i j; simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, smul_apply, add_apply]
  -- Extract W-block action
  have hD_act (Ω : Matrix W W k) (a b : k)
      (hact : g * fromBlocks (0 : Matrix I I k) 0 0 Ω * gᵀ =
        a • rep₁ (k := k) + b • rep₂ (k := k)) :
      g.toBlocks₂₂ * Ω * g.toBlocks₂₂ᵀ =
      a • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        b • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) := by
    rw [← w_block_lhs, ← w_block_rhs]; exact congrArg toBlocks₂₂ hact
  have hD_act1 := hD_act _ _ _ hact1
  have hD_act2 := hD_act _ _ _ hact2
  -- Step 2: det(D) ≠ 0 via vector kernel argument
  have hdetD : Matrix.det g.toBlocks₂₂ ≠ 0 :=
    det_ne_zero_of_preserves_span (k := k) g.toBlocks₂₂ M hdetM hD_act1 hD_act2
  -- Step 3: D satisfies N5 ActBivector / PreservesSubspaceBivector
  have hD1 : Wedge2Formalization.N5SimplePoint.ActBivector
      (Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k)) g.toBlocks₂₂ =
    M 0 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
      M 0 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) := by
    simp only [Wedge2Formalization.N5SimplePoint.ActBivector]; exact hD_act1
  have hD2 : Wedge2Formalization.N5SimplePoint.ActBivector
      (Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k)) g.toBlocks₂₂ =
    M 1 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
      M 1 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) := by
    simp only [Wedge2Formalization.N5SimplePoint.ActBivector]; exact hD_act2
  -- Step 4: Delegate to N5.Row3.quotient_image to get M' ∈ Qproj(N5.Row3)
  have hD_pres : Wedge2Formalization.N5SimplePoint.PreservesSubspaceBivector (k := k)
      g.toBlocks₂₂ :=
    ⟨M 0 0, M 0 1, M 1 0, M 1 1, hD1, hD2⟩
  rcases Wedge2Formalization.Paper.N5.Row3.quotient_image (k := k)
    g.toBlocks₂₂ hdetD hD_pres with ⟨M', hM'mem, hM'act1, hM'act2⟩
  -- Step 5: Uniqueness M' = M from N5-level rep independence
  have hM'_eq : M' = M := by
    have heq1 : M 0 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        M 0 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) =
      M' 0 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        M' 0 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) := hD1.symm.trans hM'act1
    have heq2 : M 1 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        M 1 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) =
      M' 1 0 • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        M' 1 1 • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) := hD2.symm.trans hM'act2
    have d1 : (M 0 0 - M' 0 0) • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        (M 0 1 - M' 0 1) • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) = 0 := by
      have := sub_eq_zero.mpr heq1
      rw [show _ - _ = (M 0 0 - M' 0 0) • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        (M 0 1 - M' 0 1) • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) from by
        simp [sub_smul]; abel] at this; exact this
    have d2 : (M 1 0 - M' 1 0) • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        (M 1 1 - M' 1 1) • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) = 0 := by
      have := sub_eq_zero.mpr heq2
      rw [show _ - _ = (M 1 0 - M' 1 0) • Wedge2Formalization.Paper.N5.Row3.rep₁ (k := k) +
        (M 1 1 - M' 1 1) • Wedge2Formalization.Paper.N5.Row3.rep₂ (k := k) from by
        simp [sub_smul]; abel] at this; exact this
    rcases n5_rep_pair_independent (k := k) d1 with ⟨e00, e01⟩
    rcases n5_rep_pair_independent (k := k) d2 with ⟨e10, e11⟩
    ext i j; fin_cases i <;> fin_cases j <;>
      [exact (sub_eq_zero.mp e00).symm; exact (sub_eq_zero.mp e01).symm;
       exact (sub_eq_zero.mp e10).symm; exact (sub_eq_zero.mp e11).symm]
  -- Transfer Qproj membership: M' ∈ N5.Row3.Qproj = { M | M 1 0 = 0 ∧ det M ≠ 0 }
  rw [← hM'_eq]
  exact ⟨hM'mem.1, hM'mem.2⟩

end Row6

end N7
end Paper
end Wedge2Formalization
