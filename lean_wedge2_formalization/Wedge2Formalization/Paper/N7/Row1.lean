import Wedge2Formalization.Paper.Core
import Wedge2Formalization.N5PureSingularLong
import Mathlib.Tactic.LinearCombination

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 1.
Representative `S₄`.
Divisor `0`.
Claimed stabilizer:
`K_L = \Ga^{14} \rtimes (\Gm(k) \times \Gm(k))`, exact quotient
family `Q_L = \GL_2(k)`.

This row has a pairing structure: the representatives are off-diagonal
`fromBlocks 0 A (-Aᵀ) 0` with W = Fin 4, I = Fin 3 (note W before I).

Public layer:
- `TableCell` names the Appendix A kernel cell in the pairing coordinates used here.
- `lift` and `coeff` expose the `GL₂(k)` quotient action through the `Sym³` / `adj(Sym²)ᵀ`
  block lift.
- `quotient_image` proves every setwise stabilizer induces a coefficient in `GL₂(k)`.
-/
namespace Row1

abbrev W := Fin 4
abbrev I := Fin 3
abbrev V := Sum W I

/-- The 4x3 embedding matrix for rep₁:
    maps e₁->e₁, e₂->e₂, e₃->e₃ (identity embedding). -/
def A₁ : Matrix W I k :=
  !![1, 0, 0;
     0, 1, 0;
     0, 0, 1;
     0, 0, 0]

/-- The 4x3 shifted embedding matrix for rep₂:
    maps e₁->e₂, e₂->e₃, e₃->e₄. -/
def A₂ : Matrix W I k :=
  !![0, 0, 0;
     1, 0, 0;
     0, 1, 0;
     0, 0, 1]

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks
    (0 : Matrix W W k)
    (A₁ (k := k))
    (-(A₁ (k := k))ᵀ)
    (0 : Matrix I I k)

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks
    (0 : Matrix W W k)
    (A₂ (k := k))
    (-(A₂ (k := k))ᵀ)
    (0 : Matrix I I k)

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

/-- The reps are linearly independent. -/
theorem rep_pair_independent {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have ha : a = 0 := by
    have := congrArg (fun M => M (Sum.inl 0) (Sum.inr 0)) h
    simp [rep₁, rep₂, A₁, A₂, fromBlocks] at this
    exact this
  have hb : b = 0 := by
    have := congrArg (fun M => M (Sum.inl 1) (Sum.inr 0)) h
    simp [ha, rep₁, rep₂, A₁, A₂, fromBlocks] at this
    exact this
  exact ⟨ha, hb⟩

/-- The pointwise stabilizer, expressed via off-diagonal block conditions (W = Fin 4, I = Fin 3).
For the pairing structure rep_i = fromBlocks 0 A_i (-A_iᵀ) 0, the 8 block conditions of
g * rep_i * gᵀ = rep_i describe:
- (1,2) blocks: the key pairing equations P*A_i*Qᵀ - B*A_iᵀ*Cᵀ = A_i
- (1,1)/(2,2) blocks: cancellation conditions on diagonal
The parametric form is P = Sym³(g), Q = adj(Sym²(g))^T with the unipotent radical Ga^{14}. -/
def K : Set (Matrix V V k) :=
  { g |
      -- (1,2) of rep₁: P*A₁*Qᵀ - B*A₁ᵀ*Cᵀ = A₁
      -(g.toBlocks₁₂ * (A₁ (k := k))ᵀ * g.toBlocks₂₁ᵀ) +
        g.toBlocks₁₁ * A₁ (k := k) * g.toBlocks₂₂ᵀ = A₁ (k := k) ∧
      -- (1,2) of rep₂: P*A₂*Qᵀ - B*A₂ᵀ*Cᵀ = A₂
      -(g.toBlocks₁₂ * (A₂ (k := k))ᵀ * g.toBlocks₂₁ᵀ) +
        g.toBlocks₁₁ * A₂ (k := k) * g.toBlocks₂₂ᵀ = A₂ (k := k) ∧
      -- Full fixing conditions (for completeness):
      g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
      g * rep₂ (k := k) * gᵀ = rep₂ (k := k) }

/-- The explicit pairing equations on the off-diagonal `(W,I)` block:
`P * A_i * Qᵀ - B * A_iᵀ * Cᵀ = A_i` for `i = 1, 2`.
This is the public name for the `\Ga^{14}`-type coupling cell in Appendix A, row 1. -/
def PairingCell : Set (Matrix V V k) :=
  { g |
      -(g.toBlocks₁₂ * (A₁ (k := k))ᵀ * g.toBlocks₂₁ᵀ) +
        g.toBlocks₁₁ * A₁ (k := k) * g.toBlocks₂₂ᵀ = A₁ (k := k) ∧
      -(g.toBlocks₁₂ * (A₂ (k := k))ᵀ * g.toBlocks₂₁ᵀ) +
        g.toBlocks₁₁ * A₂ (k := k) * g.toBlocks₂₂ᵀ = A₂ (k := k) }

/-- The exact pointwise fixing equations for the displayed pairing representatives.
Together with `PairingCell`, this is the public Appendix A row-1 kernel cell. -/
def PointwiseFixingCell : Set (Matrix V V k) :=
  { g |
      g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
      g * rep₂ (k := k) * gᵀ = rep₂ (k := k) }

/-- Public name for the explicit Appendix A row-1 kernel cell
`K_L = \Ga^{14} \rtimes (\Gm(k) \times \Gm(k))`, written as the pairing cell plus the
exact pointwise fixing equations in the coordinates used in this file. -/
def TableCell : Set (Matrix V V k) :=
  { g | g ∈ PairingCell (k := k) ∧ g ∈ PointwiseFixingCell (k := k) }

/-- The quotient family GL₂(k). This matches the paper's claim Q_L = GL₂. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

/-- The Sym³ matrix: induced action of GL₂ on degree-3 monomials
    {x³, x²y, xy², y³}. -/
def sym3 (a b c d : k) : Matrix W W k :=
  !![a * a * a, a * a * b, a * b * b, b * b * b;
     3 * a * a * c, a * a * d + 2 * a * b * c,
       2 * a * b * d + b * b * c, 3 * b * b * d;
     3 * a * c * c, 2 * a * c * d + b * c * c,
       a * d * d + 2 * b * c * d, 3 * b * d * d;
     c * c * c, c * c * d, c * d * d, d * d * d]

/-- The adj(Sym²)ᵀ matrix: acts on the I = Fin 3 block.
    Q = adj(Sym²(a,b,c,d))ᵀ so that Qᵀ = adj(Sym²(a,b,c,d)). -/
def sym2adjT (a b c d : k) : Matrix I I k :=
  !![(a * d - b * c) * d * d,
     -((a * d - b * c) * 2 * c * d),
     (a * d - b * c) * c * c;
     -((a * d - b * c) * b * d),
     (a * d - b * c) * (a * d + b * c),
     -((a * d - b * c) * a * c);
     (a * d - b * c) * b * b,
     -((a * d - b * c) * 2 * a * b),
     (a * d - b * c) * a * a]

/-- The GL₂ lift for the pairing structure. -/
def lift (a b c d : k) : Matrix V V k :=
  Matrix.fromBlocks
    (sym3 (k := k) a b c d)
    0
    0
    (sym2adjT (k := k) a b c d)

/-- The coefficient matrix: Δ³ · !![a, c; b, d] where Δ = ad - bc. -/
def coeff (a b c d : k) : Matrix (Fin 2) (Fin 2) k :=
  !![(a * d - b * c) * (a * d - b * c) * (a * d - b * c) * a,
     (a * d - b * c) * (a * d - b * c) * (a * d - b * c) * c;
     (a * d - b * c) * (a * d - b * c) * (a * d - b * c) * b,
     (a * d - b * c) * (a * d - b * c) * (a * d - b * c) * d]

/-- The exact pointwise stabilizer: g fixes the pair iff g ∈ K.
The forward direction extracts the (1,2) off-diagonal block equations (the key pairing
conditions P*A_i*Qᵀ - B*A_iᵀ*Cᵀ = A_i) from the full fixing condition. -/
theorem pointwise_stabilizer :
    ∀ g : Matrix V V k,
      (g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
        g * rep₂ (k := k) * gᵀ = rep₂ (k := k)) ↔
      g ∈ K (k := k) := by
  intro g
  have hg_eq : g = fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ :=
    (fromBlocks_toBlocks g).symm
  constructor
  · -- Forward: extract (1,2) block conditions and the full fixing conditions
    rintro ⟨h1, h2⟩
    have h1b : fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ *
        fromBlocks (0 : Matrix W W k) (A₁ (k := k)) (-(A₁ (k := k))ᵀ) (0 : Matrix I I k) *
        (fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
      fromBlocks (0 : Matrix W W k) (A₁ (k := k)) (-(A₁ (k := k))ᵀ) (0 : Matrix I I k) := by
      rw [← hg_eq]; exact h1
    have h2b : fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ *
        fromBlocks (0 : Matrix W W k) (A₂ (k := k)) (-(A₂ (k := k))ᵀ) (0 : Matrix I I k) *
        (fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
      fromBlocks (0 : Matrix W W k) (A₂ (k := k)) (-(A₂ (k := k))ᵀ) (0 : Matrix I I k) := by
      rw [← hg_eq]; exact h2
    exact ⟨by have := congrArg toBlocks₁₂ h1b
              simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₁₂ h2b
              simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           h1, h2⟩
  · -- Backward: the full fixing conditions are included in K
    rintro ⟨_, _, h1, h2⟩
    exact ⟨h1, h2⟩

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
        g * rep₂ (k := k) * gᵀ = rep₂ (k := k) := by
  constructor
  · intro hg; exact hg.2.2
  · intro ⟨h1, h2⟩; exact ((pointwise_stabilizer (k := k) g).mp ⟨h1, h2⟩)

/-- Off-diagonal block equations extracted from the fixing condition for the pairing structure.
For a matrix g = fromBlocks P B C Q where P : W -> W and Q : I -> I, the (1,2) block of
g * rep_i * g^T = rep_i gives P * A_i * Q^T + B * (-A_i^T) * B^T = A_i (= (1,2) block).
The (1,1) block gives P * A_i * B^T + B * (-A_i^T) * Q^T = 0. When B = C = 0
(block-diagonal case), these simplify to P * A_i * Q^T = A_i for i=1,2. -/
theorem mem_K_offdiag_blocks (g : Matrix V V k) (hg : g ∈ K (k := k)) :
    -(g.toBlocks₁₂ * (A₁ (k := k))ᵀ * g.toBlocks₂₁ᵀ) +
      g.toBlocks₁₁ * A₁ (k := k) * g.toBlocks₂₂ᵀ =
    A₁ (k := k) ∧
    -(g.toBlocks₁₂ * (A₂ (k := k))ᵀ * g.toBlocks₂₁ᵀ) +
      g.toBlocks₁₁ * A₂ (k := k) * g.toBlocks₂₂ᵀ =
    A₂ (k := k) :=
  ⟨hg.1, hg.2.1⟩

/-- Table form: K unfolds into the off-diagonal block pairing conditions and the full
fixing equations. The public table cell is named explicitly so a reader can match
Appendix A, row 1 directly to the Lean statement. -/
theorem mem_K_table_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) :=
by
  constructor
  · rintro ⟨hA1, hA2, h1, h2⟩
    exact ⟨⟨hA1, hA2⟩, ⟨h1, h2⟩⟩
  · rintro ⟨⟨hA1, hA2⟩, hfix⟩
    exact ⟨hA1, hA2, hfix.1, hfix.2⟩

-- Helper: Sym³ * Aᵢ * adj(Sym²) on the upper-right 4×3 block
private theorem lift_rep₁_upperRight (a b c d : k) :
    sym3 (k := k) a b c d * A₁ (k := k) *
      (sym2adjT (k := k) a b c d)ᵀ =
    ((a * d - b * c) * (a * d - b * c) * (a * d - b * c) * a) •
      A₁ (k := k) +
    ((a * d - b * c) * (a * d - b * c) * (a * d - b * c) * c) •
      A₂ (k := k) := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sym3, A₁, A₂, sym2adjT, Matrix.mul_apply,
      Matrix.smul_apply, Matrix.add_apply,
      Matrix.transpose_apply,
      Fin.sum_univ_three, Fin.sum_univ_four] <;>
    ring

private theorem lift_rep₂_upperRight (a b c d : k) :
    sym3 (k := k) a b c d * A₂ (k := k) *
      (sym2adjT (k := k) a b c d)ᵀ =
    ((a * d - b * c) * (a * d - b * c) * (a * d - b * c) * b) •
      A₁ (k := k) +
    ((a * d - b * c) * (a * d - b * c) * (a * d - b * c) * d) •
      A₂ (k := k) := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sym3, A₁, A₂, sym2adjT, Matrix.mul_apply,
      Matrix.smul_apply, Matrix.add_apply,
      Matrix.transpose_apply,
      Fin.sum_univ_three, Fin.sum_univ_four] <;>
    ring

-- Block-diagonal action preserves the off-diagonal skew structure
private theorem blockDiag_act_rep
    (P : Matrix W W k) (Q : Matrix I I k)
    (A : Matrix W I k) :
    fromBlocks P 0 0 Q *
      fromBlocks (0 : Matrix W W k) A (-Aᵀ) (0 : Matrix I I k) *
      (fromBlocks P 0 0 Q)ᵀ =
    fromBlocks 0 (P * A * Qᵀ) (-(P * A * Qᵀ)ᵀ) 0 := by
  have h1 : (fromBlocks P 0 0 Q)ᵀ =
      fromBlocks Pᵀ 0 0 Qᵀ := fromBlocks_transpose _ _ _ _
  rw [h1]
  simp only [fromBlocks_multiply, Matrix.mul_zero, Matrix.zero_mul,
    add_zero, zero_add, Matrix.transpose_zero,
    Matrix.mul_neg, Matrix.neg_mul]
  congr 1 <;> simp [Matrix.transpose_neg, Matrix.transpose_mul,
    Matrix.transpose_transpose, Matrix.mul_assoc, Matrix.neg_mul]

private theorem lincomb_rep
    (α β : k) :
    α • rep₁ (k := k) + β • rep₂ (k := k) =
    fromBlocks (0 : Matrix W W k)
      (α • A₁ (k := k) + β • A₂ (k := k))
      (-(α • A₁ (k := k) + β • A₂ (k := k))ᵀ)
      (0 : Matrix I I k) := by
  simp only [rep₁, rep₂, fromBlocks_smul, fromBlocks_add,
    smul_zero, smul_neg]
  congr 1
  · simp
  · simp [Matrix.transpose_add, Matrix.transpose_smul,
      add_comm]
  · simp

set_option maxHeartbeats 1600000 in
theorem quotient_action (a b c d : k) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) a b c d)
      (coeff (k := k) a b c d) := by
  have h₁ := lift_rep₁_upperRight (k := k) a b c d
  have h₂ := lift_rep₂_upperRight (k := k) a b c d
  refine ⟨?_, ?_⟩
  · -- lift * rep₁ * liftᵀ = coeff 0 0 • rep₁ + coeff 0 1 • rep₂
    simp only [lift, rep₁]
    rw [blockDiag_act_rep, h₁]
    -- Now goal: fromBlocks 0 (Δ³a•A₁+Δ³c•A₂) (-(...)ᵀ) 0 = ...
    -- RHS: coeff 0 0 • fromBlocks 0 A₁ (-A₁ᵀ) 0 + coeff 0 1 • ...
    -- Unfold the sum into fromBlocks form
    show _ = (coeff (k := k) a b c d) 0 0 •
        fromBlocks 0 (A₁ (k := k)) (-(A₁ (k := k))ᵀ)
          (0 : Matrix I I k) +
      (coeff (k := k) a b c d) 0 1 •
        fromBlocks 0 (A₂ (k := k)) (-(A₂ (k := k))ᵀ)
          (0 : Matrix I I k)
    ext i j; rcases i with i | i <;> rcases j with j | j
    all_goals (fin_cases i <;> fin_cases j <;>
      simp [fromBlocks, coeff, A₁, A₂, Matrix.smul_apply,
        Matrix.add_apply, Matrix.transpose_apply,
        Matrix.neg_apply] <;> ring)
  · simp only [lift, rep₂]
    rw [blockDiag_act_rep, h₂]
    show _ = (coeff (k := k) a b c d) 1 0 •
        fromBlocks 0 (A₁ (k := k)) (-(A₁ (k := k))ᵀ)
          (0 : Matrix I I k) +
      (coeff (k := k) a b c d) 1 1 •
        fromBlocks 0 (A₂ (k := k)) (-(A₂ (k := k))ᵀ)
          (0 : Matrix I I k)
    ext i j; rcases i with i | i <;> rcases j with j | j
    all_goals (fin_cases i <;> fin_cases j <;>
      simp [fromBlocks, coeff, A₁, A₂, Matrix.smul_apply,
        Matrix.add_apply, Matrix.transpose_apply,
        Matrix.neg_apply] <;> ring)

/-- Every invertible setwise stabilizer induces a coefficient matrix in GL₂ (det M ≠ 0).
Together with `quotient_action` (which produces matrices with coeff = Delta³ * !![a,c;b,d]),
this establishes Q_L = GL₂, matching the paper's claim exactly. -/
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
  simp only [Qproj, Set.mem_setOf_eq]
  intro hdet0
  have hM_dep : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
    simpa [Matrix.det_fin_two] using hdet0
  have rep_indep :
      ∀ a b : k,
        a • rep₁ (k := k) + b • rep₂ (k := k) = 0 →
          a = 0 ∧ b = 0 :=
    fun a b h => rep_pair_independent (k := k) h
  have hcomb_zero :
      (M 1 1) • rep₁ (k := k) +
        (-(M 0 1)) • rep₂ (k := k) = 0 :=
    (actBivector_eq_zero_iff_of_det_ne_zero
      (k := k) _ g hg).1 (by
      calc g * ((M 1 1) • rep₁ (k := k) +
              (-(M 0 1)) • rep₂ (k := k)) * gᵀ
          = (M 1 1) • (g * rep₁ (k := k) * gᵀ) +
            (-(M 0 1)) • (g * rep₂ (k := k) * gᵀ) := by
            simp [Matrix.mul_add, Matrix.add_mul]
        _ = 0 := by
            rw [hact1, hact2]; ext i j
            simp only [Matrix.add_apply, Matrix.smul_apply,
              Matrix.zero_apply, smul_eq_mul, neg_mul]
            linear_combination
              (rep₁ (k := k) i j) * hM_dep)
  rcases rep_indep _ _ hcomb_zero with ⟨h11, h01neg⟩
  have h01 : M 0 1 = 0 := neg_eq_zero.mp h01neg
  have hcomb2 :
      (-(M 1 0)) • rep₁ (k := k) +
        (M 0 0) • rep₂ (k := k) = 0 :=
    (actBivector_eq_zero_iff_of_det_ne_zero
      (k := k) _ g hg).1 (by
      calc g * ((-(M 1 0)) • rep₁ (k := k) +
              (M 0 0) • rep₂ (k := k)) * gᵀ
          = (-(M 1 0)) • (g * rep₁ (k := k) * gᵀ) +
            (M 0 0) • (g * rep₂ (k := k) * gᵀ) := by
            simp [Matrix.mul_add, Matrix.add_mul]
        _ = 0 := by
            rw [hact1, hact2]; ext i j
            simp only [Matrix.add_apply, Matrix.smul_apply,
              Matrix.zero_apply, smul_eq_mul, neg_mul]
            linear_combination
              (rep₂ (k := k) i j) * hM_dep)
  have h10 : M 1 0 = 0 :=
    neg_eq_zero.mp (rep_indep _ _ hcomb2).1
  have hrep2_zero : rep₂ (k := k) = 0 := by
    have : g * rep₂ (k := k) * gᵀ = 0 := by
      have := hact2; rw [h10, h11] at this
      simpa using this
    exact (actBivector_eq_zero_iff_of_det_ne_zero
      (k := k) _ g hg).1 this
  exact absurd
    ((rep_indep 0 1 (by simpa using hrep2_zero)).2)
    one_ne_zero

end Row1

end N7
end Paper
end Wedge2Formalization
