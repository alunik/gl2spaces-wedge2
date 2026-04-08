import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N5
import Mathlib.Tactic.LinearCombination

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 2.
Representative `S₃ + S₁²`.
Divisor `0`.
Claimed stabilizer:
`K_L = \Ga^{10} \rtimes (\Gm(k) \times (U_4 \rtimes \Gm(k)))`, exact quotient
family `Q_L = \GL_2(k)`.

This row is a radical-2 extension of N5 Row 1 (PureSingularLong).
-/
namespace Row2

abbrev I := Fin 2
abbrev W := Wedge2Formalization.N5PureSingularLong.V
abbrev V := Sum I W

def scalarBlock (a b : k) : Matrix I I k := !![a, b; 0, a]

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k))

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k))

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

/-- The reps are linearly independent (from the W-block N5.Row1 independence). -/
theorem rep_pair_independent {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) : a = 0 ∧ b = 0 := by
  -- Extract entry (inr(inl 0), inr(inr 0)): a * mulX(0,0) + b * mulY(0,0) = a = 0
  have ha : a = 0 := by
    have := congrArg (fun M => M (Sum.inr (Sum.inl 0)) (Sum.inr (Sum.inr 0))) h
    simp [rep₁, rep₂, Wedge2Formalization.Paper.N5.Row1.rep₁,
      Wedge2Formalization.Paper.N5.Row1.rep₂,
      Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Wedge2Formalization.N5PureSingularLong.mulY,
      fromBlocks] at this
    exact this
  -- Extract entry (inr(inl 1), inr(inr 0)): 0 + b * mulY(1,0) = b = 0
  have hb : b = 0 := by
    have := congrArg (fun M => M (Sum.inr (Sum.inl 1)) (Sum.inr (Sum.inr 0))) h
    simp [ha, rep₁, rep₂, Wedge2Formalization.Paper.N5.Row1.rep₁,
      Wedge2Formalization.Paper.N5.Row1.rep₂,
      Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Wedge2Formalization.N5PureSingularLong.mulY,
      fromBlocks] at this
    exact this
  exact ⟨ha, hb⟩

/-- Public name for the embedded Appendix A, `n = 5`, row 1 core
`U_4 \rtimes G_m(k)` on the lower-right `5`-space. -/
def LongPureCore : Set (Matrix W W k) :=
  Wedge2Formalization.Paper.N5.Row1.K (k := k)

/-- The kernel: parametric form following Row 3/Row 14 pattern. -/
def K : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ LongPureCore (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

def lift (α β γ δ : k) : Matrix V V k :=
  Matrix.fromBlocks (1 : Matrix I I k) 0 0
    (Wedge2Formalization.Paper.N5.Row1.lift (k := k) α β γ δ)

/-- The common right-kernel of the N5PureSingularLong rep pair is trivial.
rep₁ kernel = span{e₃}, rep₂ kernel = span{e₁}, so their intersection is {0}. -/
private theorem commonKernel_zero
    (v : Wedge2Formalization.N5PureSingularLong.V → k)
    (hv₁ : Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) *ᵥ v = 0)
    (hv₂ : Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) *ᵥ v = 0) :
    v = 0 := by
  funext i
  cases i with
  | inl i =>
      fin_cases i
      · -- v (inl 0) = 0: from rep₁ entry (inr 0): -(v(inl 0)) = 0
        have := congrArg (fun w => w (Sum.inr 0)) hv₁
        simp [Wedge2Formalization.Paper.N5.Row1.rep₁,
          Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct,
          Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] at this
        exact this
      · -- v (inl 1) = 0: from rep₁ entry (inr 1): -(v(inl 1)) = 0
        have := congrArg (fun w => w (Sum.inr 1)) hv₁
        simp [Wedge2Formalization.Paper.N5.Row1.rep₁,
          Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct,
          Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] at this
        exact this
      · -- v (inl 2) = 0: from rep₂ entry (inr 1): -(v(inl 2)) = 0
        have := congrArg (fun w => w (Sum.inr 1)) hv₂
        simp [Wedge2Formalization.Paper.N5.Row1.rep₂,
          Wedge2Formalization.N5PureSingularLong.rep₂,
          Wedge2Formalization.N5PureSingularLong.mulY,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct,
          Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] at this
        exact this
  | inr i =>
      fin_cases i
      · -- v (inr 0) = 0: from rep₁ entry (inl 0): v(inr 0) = 0
        have := congrArg (fun w => w (Sum.inl 0)) hv₁
        simp [Wedge2Formalization.Paper.N5.Row1.rep₁,
          Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct,
          Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] at this
        exact this
      · -- v (inr 1) = 0: from rep₁ entry (inl 1): v(inr 1) = 0
        have := congrArg (fun w => w (Sum.inl 1)) hv₁
        simp [Wedge2Formalization.Paper.N5.Row1.rep₁,
          Wedge2Formalization.N5PureSingularLong.rep₁,
          Wedge2Formalization.N5PureSingularLong.mulX,
          Matrix.mulVec, Matrix.fromBlocks, dotProduct,
          Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two] at this
        exact this

/-- Any matrix fixing the N5PureSingularLong pair has nonzero determinant.
Proof: if det = 0, there's v ≠ 0 in common kernel of both reps, contradicting triviality. -/
private theorem det_ne_zero_of_fixing_pair
    (D : Matrix W W k)
    (hDfix : Wedge2Formalization.N5PureSingularLong.FixesPairBivector (k := k) D) :
    Matrix.det D ≠ 0 := by
  intro hD0
  obtain ⟨v, hvne, hv0⟩ :=
    (Matrix.exists_mulVec_eq_zero_iff (M := Dᵀ)).2
      (by simpa [Matrix.det_transpose] using hD0)
  have hv₁ : Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) *ᵥ v = 0 := by
    have hzero :
        Wedge2Formalization.N5PureSingularLong.ActBivector
          (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) D *ᵥ v = 0 := by
      calc
        Wedge2Formalization.N5PureSingularLong.ActBivector
            (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) D *ᵥ v
            = D *ᵥ (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) *ᵥ (Dᵀ *ᵥ v)) := by
              simp [Wedge2Formalization.N5PureSingularLong.ActBivector,
                Matrix.mulVec_mulVec, Matrix.mul_assoc]
        _ = 0 := by simp [hv0]
    have hfix₁ : Wedge2Formalization.N5PureSingularLong.ActBivector
        (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) D =
      Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) := by
      simpa [Wedge2Formalization.Paper.N5.Row1.rep₁] using hDfix.1
    rw [hfix₁] at hzero
    exact hzero
  have hv₂ : Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) *ᵥ v = 0 := by
    have hzero :
        Wedge2Formalization.N5PureSingularLong.ActBivector
          (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) D *ᵥ v = 0 := by
      calc
        Wedge2Formalization.N5PureSingularLong.ActBivector
            (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) D *ᵥ v
            = D *ᵥ (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) *ᵥ (Dᵀ *ᵥ v)) := by
              simp [Wedge2Formalization.N5PureSingularLong.ActBivector,
                Matrix.mulVec_mulVec, Matrix.mul_assoc]
        _ = 0 := by simp [hv0]
    have hfix₂ : Wedge2Formalization.N5PureSingularLong.ActBivector
        (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) D =
      Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) := by
      simpa [Wedge2Formalization.Paper.N5.Row1.rep₂] using hDfix.2
    rw [hfix₂] at hzero
    exact hzero
  exact hvne (commonKernel_zero (k := k) v hv₁ hv₂)

/-- For the N5PureSingularLong orbit, the upper-right block vanishes when the pair is
fixed. The argument uses:
1. D fixes pair -> det(D) ≠ 0 (from common kernel triviality)
2. B * rep_i * Dᵀ = 0 with Dᵀ invertible -> B * rep_i = 0
3. B in left-kernel of both reps -> B = 0 (by entry extraction) -/
private theorem upperRight_zero_of_fixing_pair
    (B : Matrix I W k) (D : Matrix W W k)
    (hDfix : Wedge2Formalization.N5PureSingularLong.FixesPairBivector (k := k) D)
    (htop₁ : B * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * Dᵀ = 0)
    (htop₂ : B * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * Dᵀ = 0) :
    B = 0 := by
  have hDdet : Matrix.det D ≠ 0 := det_ne_zero_of_fixing_pair (k := k) D hDfix
  have hDetTr : Matrix.det Dᵀ ≠ 0 := by rwa [Matrix.det_transpose]
  have hDTu : IsUnit (Matrix.det Dᵀ) := IsUnit.mk0 _ hDetTr
  have hDTinv : Dᵀ * (Dᵀ)⁻¹ = (1 : Matrix W W k) := Matrix.mul_nonsing_inv _ hDTu
  have cancel : ∀ (X : Matrix I W k), X * Dᵀ = 0 → X = 0 := by
    intro X hX
    have : X = X * Dᵀ * (Dᵀ)⁻¹ := by
      conv_lhs => rw [show X = X * (1 : Matrix W W k) from (Matrix.mul_one X).symm]
      rw [← hDTinv, ← Matrix.mul_assoc]
    rw [this, hX, Matrix.zero_mul]
  have hBR1 : B * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) = 0 := cancel _ htop₁
  have hBR2 : B * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) = 0 := cancel _ htop₂
  -- B is 2×5. We show each entry is zero via entry extraction from B * rep = 0.
  -- From B * rep₁ = 0 at (i, inr 0): B(i, inl 0) = 0
  -- From B * rep₁ = 0 at (i, inr 1): B(i, inl 1) = 0
  -- From B * rep₁ = 0 at (i, inl 0): B(i, inr 0) = 0
  -- From B * rep₁ = 0 at (i, inl 1): B(i, inr 1) = 0 (from entry (i, inl 2) = 0 trivially)
  -- From B * rep₂ = 0 at (i, inr 1): B(i, inl 2) = 0
  have simp_defs := @fun (i : I) (j : W) =>
    show (B * Wedge2Formalization.N5PureSingularLong.rep₁ (k := k)) i j = 0 from
      congr_fun (congr_fun hBR1 i) j
  have simp_defs2 := @fun (i : I) (j : W) =>
    show (B * Wedge2Formalization.N5PureSingularLong.rep₂ (k := k)) i j = 0 from
      congr_fun (congr_fun hBR2 i) j
  -- Extract each entry B(i,j) = 0 from B * rep = 0.
  -- Use simpa to combine the simp on hypothesis and goal.
  have hB_0_l0 : B 0 (Sum.inl 0) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two] using simp_defs 0 (Sum.inr 0)
  have hB_0_l1 : B 0 (Sum.inl 1) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two] using simp_defs 0 (Sum.inr 1)
  have hB_0_l2 : B 0 (Sum.inl 2) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.mulY,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two] using simp_defs2 0 (Sum.inr 1)
  have hB_0_r0 : B 0 (Sum.inr 0) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two,
      Matrix.transpose_apply, Matrix.neg_apply] using simp_defs 0 (Sum.inl 0)
  have hB_0_r1 : B 0 (Sum.inr 1) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two,
      Matrix.transpose_apply, Matrix.neg_apply] using simp_defs 0 (Sum.inl 1)
  have hB_1_l0 : B 1 (Sum.inl 0) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two] using simp_defs 1 (Sum.inr 0)
  have hB_1_l1 : B 1 (Sum.inl 1) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two] using simp_defs 1 (Sum.inr 1)
  have hB_1_l2 : B 1 (Sum.inl 2) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.mulY,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two] using simp_defs2 1 (Sum.inr 1)
  have hB_1_r0 : B 1 (Sum.inr 0) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two,
      Matrix.transpose_apply, Matrix.neg_apply] using simp_defs 1 (Sum.inl 0)
  have hB_1_r1 : B 1 (Sum.inr 1) = 0 := by
    simpa [Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Matrix.mul_apply, fromBlocks, Fintype.sum_sum_type,
      Fin.sum_univ_three, Fin.sum_univ_two,
      Matrix.transpose_apply, Matrix.neg_apply] using simp_defs 1 (Sum.inl 1)
  ext i j
  fin_cases i <;> rcases j with j | j <;> fin_cases j <;>
    first | exact hB_0_l0 | exact hB_0_l1 | exact hB_0_l2
          | exact hB_0_r0 | exact hB_0_r1 | exact hB_1_l0
          | exact hB_1_l1 | exact hB_1_l2 | exact hB_1_r0 | exact hB_1_r1

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
  have hgEq : g = fromBlocks A0 B C D := (fromBlocks_toBlocks g).symm
  constructor
  · rintro ⟨h1, h2⟩
    have h1_blocks : fromBlocks A0 B C D *
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) *
        (fromBlocks A0 B C D)ᵀ =
      fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) := by
      rw [← hgEq]; exact h1
    have h2_blocks : fromBlocks A0 B C D *
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) *
        (fromBlocks A0 B C D)ᵀ =
      fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) := by
      rw [← hgEq]; exact h2
    -- Extract (1,2) blocks
    have htop₁ : B * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * Dᵀ = 0 := by
      have := congrArg toBlocks₁₂ h1_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    have htop₂ : B * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * Dᵀ = 0 := by
      have := congrArg toBlocks₁₂ h2_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- Extract (2,2) blocks: D fixes N5 pair
    have hDfix : Wedge2Formalization.N5PureSingularLong.FixesPairBivector (k := k) D := by
      constructor
      · have := congrArg toBlocks₂₂ h1_blocks
        simpa [fromBlocks_multiply, fromBlocks_transpose] using this
      · have := congrArg toBlocks₂₂ h2_blocks
        simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- B = 0
    have hB0 : B = 0 := upperRight_zero_of_fixing_pair (k := k) B D hDfix htop₁ htop₂
    -- D ∈ K: use det(D) ≠ 0 and pointwise_stabilizer
    have hDdet : Matrix.det D ≠ 0 := det_ne_zero_of_fixing_pair (k := k) D hDfix
    have hDmem : D ∈ Wedge2Formalization.Paper.N5.Row1.K (k := k) :=
      (Wedge2Formalization.Paper.N5.Row1.pointwise_stabilizer (k := k) D hDdet).1 hDfix
    exact ⟨A0, C, D, hDmem, by rw [hgEq]; simp [hB0]⟩
  · rintro ⟨A0', C', D', hD', rfl⟩
    -- D' ∈ K means D' = U * Levi a with a ≠ 0, so det(D') ≠ 0
    have hD'det : Matrix.det D' ≠ 0 := by
      rcases (Wedge2Formalization.Paper.N5.Row1.mem_K_iff (k := k) D').1 hD' with
        ⟨a, u, v, w, z, ha, hshape⟩
      rw [hshape]
      simp [Wedge2Formalization.Paper.N5.Row1.U,
        Wedge2Formalization.Paper.N5.Row1.Levi,
        Wedge2Formalization.N5PureSingularLong.lowerHankel,
        Matrix.det_mul, Matrix.det_one, ha, inv_ne_zero ha]
    have hDfix : Wedge2Formalization.N5PureSingularLong.FixesPairBivector (k := k) D' :=
      (Wedge2Formalization.Paper.N5.Row1.pointwise_stabilizer (k := k) D' hD'det).2 hD'
    constructor
    · show fromBlocks A0' 0 C' D' *
          fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.1
    · show fromBlocks A0' 0 C' D' *
          fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.2

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ Wedge2Formalization.Paper.N5.Row1.K (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D := Iff.rfl

/-- Public name for the explicit Appendix A row-2 kernel cell. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ a u v w z : k,
        ∃ A0 : Matrix I I k,
          ∃ C : Matrix W I k,
            a ≠ 0 ∧
            g = Matrix.fromBlocks A0 0 C
                (Wedge2Formalization.Paper.N5.Row1.U (k := k) u v w z *
                  Wedge2Formalization.Paper.N5.Row1.Levi (k := k) a) }

theorem mem_K_table_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) := by
  constructor
  · rintro ⟨A0, C, D, hD, rfl⟩
    rcases (Wedge2Formalization.Paper.N5.Row1.mem_K_table_iff (k := k) D).1 hD with
      ⟨a, u, v, w, z, ha, hshape⟩
    exact ⟨a, u, v, w, z, A0, C, ha, by rw [hshape]⟩
  · rintro ⟨a, u, v, w, z, A0, C, ha, rfl⟩
    exact ⟨A0, C, _,
      (Wedge2Formalization.Paper.N5.Row1.mem_K_table_iff (k := k) _).2
        ⟨a, u, v, w, z, ha, rfl⟩, rfl⟩

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

theorem quotient_action (α β γ δ : k) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) α β γ δ)
      (Wedge2Formalization.Paper.N5.Row1.coeff (k := k) α β γ δ) := by
  -- N5.Row1 quotient action
  have h5 := Wedge2Formalization.Paper.N5.Row1.quotient_action (k := k) α β γ δ
  simp only [ActsOnOrderedPair] at h5
  refine ⟨?_, ?_⟩
  · change lift (k := k) α β γ δ * rep₁ (k := k) * (lift (k := k) α β γ δ)ᵀ =
        (Wedge2Formalization.Paper.N5.Row1.coeff (k := k) α β γ δ) 0 0 •
          rep₁ (k := k) +
        (Wedge2Formalization.Paper.N5.Row1.coeff (k := k) α β γ δ) 0 1 •
          rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h5.1
  · change lift (k := k) α β γ δ * rep₂ (k := k) * (lift (k := k) α β γ δ)ᵀ =
        (Wedge2Formalization.Paper.N5.Row1.coeff (k := k) α β γ δ) 1 0 •
          rep₁ (k := k) +
        (Wedge2Formalization.Paper.N5.Row1.coeff (k := k) α β γ δ) 1 1 •
          rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h5.2

/-- Every invertible setwise stabilizer induces a coefficient matrix in the projective GL₂.
The proof uses rep independence and the actBivector_eq_zero_iff_of_det_ne_zero argument. -/
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
  dsimp only at hact1 hact2
  refine ⟨M, ?_, hact1, hact2⟩
  -- Qproj = { M | det M ≠ 0 }
  simp only [Qproj, Set.mem_setOf_eq]
  intro hdet0
  have hM_dep : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
    simpa [Matrix.det_fin_two] using hdet0
  have rep_indep : ∀ a b : k, a • rep₁ (k := k) + b • rep₂ (k := k) = 0 → a = 0 ∧ b = 0 :=
    fun a b h => rep_pair_independent (k := k) h
  -- If det M = 0, the combination (M 1 1)•rep₁ + (-(M 0 1))•rep₂ maps to 0
  have hcomb_zero : (M 1 1) • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k) = 0 :=
    (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
      calc g * ((M 1 1) • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k)) * gᵀ
          = (M 1 1) • (g * rep₁ (k := k) * gᵀ) +
            (-(M 0 1)) • (g * rep₂ (k := k) * gᵀ) := by
            simp [Matrix.mul_add, Matrix.add_mul]
        _ = 0 := by
            rw [hact1, hact2]; ext i j
            simp only [Matrix.add_apply, Matrix.smul_apply,
              Matrix.zero_apply, smul_eq_mul, neg_mul]
            linear_combination (rep₁ (k := k) i j) * hM_dep)
  rcases rep_indep _ _ hcomb_zero with ⟨h11, h01neg⟩
  have h01 : M 0 1 = 0 := neg_eq_zero.mp h01neg
  -- Similarly the other combination
  have hcomb2 : (-(M 1 0)) • rep₁ (k := k) + (M 0 0) • rep₂ (k := k) = 0 :=
    (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
      calc g * ((-(M 1 0)) • rep₁ (k := k) + (M 0 0) • rep₂ (k := k)) * gᵀ
          = (-(M 1 0)) • (g * rep₁ (k := k) * gᵀ) +
            (M 0 0) • (g * rep₂ (k := k) * gᵀ) := by
            simp [Matrix.mul_add, Matrix.add_mul]
        _ = 0 := by
            rw [hact1, hact2]; ext i j
            simp only [Matrix.add_apply, Matrix.smul_apply,
              Matrix.zero_apply, smul_eq_mul, neg_mul]
            linear_combination (rep₂ (k := k) i j) * hM_dep)
  have h10 : M 1 0 = 0 := neg_eq_zero.mp (rep_indep _ _ hcomb2).1
  -- All of M except M 0 0 is zero. hact2: g rep₂ gᵀ = 0 • rep₁ + 0 • rep₂ = 0
  have hrep2_zero : rep₂ (k := k) = 0 := by
    have : g * rep₂ (k := k) * gᵀ = 0 := by
      have := hact2; rw [h10, h11] at this; simpa using this
    exact (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 this
  exact absurd ((rep_indep 0 1 (by simpa using hrep2_zero)).2) one_ne_zero

end Row2

end N7
end Paper
end Wedge2Formalization
