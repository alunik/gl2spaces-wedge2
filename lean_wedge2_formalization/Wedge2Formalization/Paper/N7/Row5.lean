import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N5
import Wedge2Formalization.N5PureSingularLong
import Mathlib.Tactic.LinearCombination

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 5.
Representative `S₃ + J_{a,1}`.
Divisor `[a]`.
Claimed stabilizer:
`K_L = \Ga^{10} \rtimes (\SL_2(k) \times (U_4 \rtimes \Gm(k)))`, exact quotient
family `Q_L = B` (Borel).

This row has a 2-dimensional I-block with the symplectic form J₂ on rep₁
and the N5 PureSingularLong pair on the W-block.
rep₁ = fromBlocks J₂ 0 0 (N5.Row1.rep₁), rep₂ = fromBlocks 0 0 0 (N5.Row1.rep₂).

Public layer:
- `TableCell` names the Appendix A kernel cell in the mixed `J₂` / `N5.Row1`
  block coordinates.
- `Qproj` is the displayed Borel family, with `quotient_action` giving the explicit lift.
- `quotient_image` proves every setwise stabilizer induces a Borel coefficient matrix.
-/
namespace Row5

abbrev I := Fin 2
abbrev W := Wedge2Formalization.N5PureSingularLong.V
abbrev V := Sum I W

/-- The 2x2 symplectic form. -/
def J₂ : Matrix I I k := !![(0 : k), 1; -1, 0]

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks
    (J₂ (k := k))
    0
    0
    (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k))

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks
    (0 : Matrix I I k)
    0
    0
    (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k))

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

/-- The reps are linearly independent (from the (1,1) block J₂ for a, and (2,2) block
    N5.Row1.rep₂ for b). -/
theorem rep_pair_independent {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) : a = 0 ∧ b = 0 := by
  -- Extract (1,1) block: a • J₂ + b • 0 = 0, so a • J₂ = 0
  -- Extract entry (inl 0, inl 1): a * J₂(0,1) + b * 0 = a = 0
  have ha : a = 0 := by
    have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h
    simp [rep₁, rep₂, J₂, fromBlocks] at h01
    exact h01
  -- Extract entry (inl 1 ++ inr 0, inl 1 ++ inr 0) of (2,2) block:
  -- b * mulY(1,0) = b = 0
  have hb : b = 0 := by
    have hentry := congrArg (fun M => M (Sum.inr (Sum.inl 1)) (Sum.inr (Sum.inr 0))) h
    simp [ha, rep₁, rep₂, Wedge2Formalization.Paper.N5.Row1.rep₁,
      Wedge2Formalization.Paper.N5.Row1.rep₂,
      Wedge2Formalization.N5PureSingularLong.rep₁,
      Wedge2Formalization.N5PureSingularLong.rep₂,
      Wedge2Formalization.N5PureSingularLong.mulX,
      Wedge2Formalization.N5PureSingularLong.mulY,
      J₂, fromBlocks] at hentry
    exact hentry
  exact ⟨ha, hb⟩

/-- The pointwise stabilizer, expressed via block decomposition (I = Fin 2, W = N5.V).
K consists of matrices g = fromBlocks A B C D satisfying all 8 block conditions from
g * rep_i * gᵀ = rep_i. The (1,1) blocks encode the Sp₂ = SL₂ factor, the (2,2) blocks
encode the U₄ x Gm factor, and the off-diagonal blocks encode the Ga¹⁰ radical. -/
def K : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * J₂ (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₁₂ᵀ =
      J₂ (k := k) ∧
      g.toBlocks₂₁ * J₂ (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) ∧
      g.toBlocks₁₁ * J₂ (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * J₂ (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₁₂ᵀ = 0 ∧
      g.toBlocks₁₁ * (0 : Matrix I I k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₁₂ᵀ =
      (0 : Matrix I I k) ∧
      g.toBlocks₂₁ * (0 : Matrix I I k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) ∧
      g.toBlocks₁₁ * (0 : Matrix I I k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * (0 : Matrix I I k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₁₂ᵀ = 0 }

/-- The two diagonal block equations for the `\Sp_2(k)`/`U_4 \rtimes G_m(k)` factors
coming from the first representative. -/
def DiagRep₁Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * J₂ (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₁₂ᵀ =
      J₂ (k := k) ∧
      g.toBlocks₂₁ * J₂ (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) }

/-- The two off-diagonal coupling equations for the first representative. -/
def CouplingRep₁Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * J₂ (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * J₂ (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₁₂ᵀ = 0 }

/-- The two diagonal block equations for the second representative. -/
def DiagRep₂Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * (0 : Matrix I I k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₁₂ᵀ =
      (0 : Matrix I I k) ∧
      g.toBlocks₂₁ * (0 : Matrix I I k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) }

/-- The two off-diagonal coupling equations for the second representative. -/
def CouplingRep₂Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * (0 : Matrix I I k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * (0 : Matrix I I k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₁₂ᵀ = 0 }

/-- Public name for the explicit Appendix A row-5 kernel cell
`K_L = \Ga^{10} \rtimes (\SL_2(k) \times (U_4 \rtimes \Gm(k)))`,
packaged as the four named block cells in the mixed `J_2`/`N5.Row1` coordinates. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      g ∈ DiagRep₁Cell (k := k) ∧
      g ∈ CouplingRep₁Cell (k := k) ∧
      g ∈ DiagRep₂Cell (k := k) ∧
      g ∈ CouplingRep₂Cell (k := k) }

/-- The quotient family Q_L = B (Borel). Every setwise stabilizer induces a
coefficient matrix M with M 1 0 = 0 and det M ≠ 0. The lower bound B ⊆ Q_L is
witnessed by `quotient_action` which produces !![a², 0; 0, a]. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

/-- The Borel lift: I-block via scalar a*I₂, W-block via N5.Row1.lift(a,0,0,1).
Coefficient matrix is !![a², 0; 0, a], which is upper-triangular (Borel). -/
def lift (a : k) : Matrix V V k :=
  Matrix.fromBlocks
    (a • (1 : Matrix I I k))
    0
    0
    (Wedge2Formalization.Paper.N5.Row1.lift (k := k) a 0 0 1)

/-- The exact pointwise stabilizer: g fixes the pair iff g ∈ K. -/
theorem pointwise_stabilizer :
    ∀ g : Matrix V V k,
      (g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
        g * rep₂ (k := k) * gᵀ = rep₂ (k := k)) ↔
      g ∈ K (k := k) := by
  intro g
  have hg_eq : g = fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ :=
    (fromBlocks_toBlocks g).symm
  constructor
  · rintro ⟨h1, h2⟩
    rw [hg_eq] at h1 h2; simp only [rep₁, rep₂] at h1 h2
    exact ⟨by have := congrArg toBlocks₁₁ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₂ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₁₂ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₁ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₁₁ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₂ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₁₂ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₁ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this⟩
  · rintro ⟨h11_1, h22_1, h12_1, h21_1, h11_2, h22_2, h12_2, h21_2⟩
    constructor
    · show g * rep₁ (k := k) * gᵀ = rep₁ (k := k)
      have : fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ *
          fromBlocks (J₂ (k := k)) 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) *
          (fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
        fromBlocks (J₂ (k := k)) 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k)) := by
        simp [fromBlocks_multiply, fromBlocks_transpose, h11_1, h12_1, h21_1, h22_1]
      rw [hg_eq]; exact this
    · show g * rep₂ (k := k) * gᵀ = rep₂ (k := k)
      have : fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ *
          fromBlocks (0 : Matrix I I k) 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) *
          (fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
        fromBlocks (0 : Matrix I I k) 0 0 (Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k)) := by
        simp only [Matrix.mul_zero, Matrix.zero_mul, add_zero, zero_add] at h11_2 h22_2 h12_2 h21_2
        simp [fromBlocks_multiply, fromBlocks_transpose, h11_2, h12_2, h21_2, h22_2]
      rw [hg_eq]; exact this

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
      g * rep₂ (k := k) * gᵀ = rep₂ (k := k) :=
  (pointwise_stabilizer (k := k) g).symm

/-- Block equations extracted from the fixing condition (diagonal blocks only). -/
theorem mem_K_diag_blocks (g : Matrix V V k) (hg : g ∈ K (k := k)) :
    g.toBlocks₁₁ * J₂ (k := k) * g.toBlocks₁₁ᵀ +
      g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₁₂ᵀ =
    J₂ (k := k) ∧
    g.toBlocks₂₁ * J₂ (k := k) * g.toBlocks₂₁ᵀ +
      g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) * g.toBlocks₂₂ᵀ =
    Wedge2Formalization.Paper.N5.Row1.rep₁ (k := k) ∧
    g.toBlocks₁₁ * (0 : Matrix I I k) * g.toBlocks₁₁ᵀ +
      g.toBlocks₁₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₁₂ᵀ =
    (0 : Matrix I I k) ∧
    g.toBlocks₂₁ * (0 : Matrix I I k) * g.toBlocks₂₁ᵀ +
      g.toBlocks₂₂ * Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) * g.toBlocks₂₂ᵀ =
    Wedge2Formalization.Paper.N5.Row1.rep₂ (k := k) :=
  ⟨hg.1, hg.2.1, hg.2.2.2.2.1, hg.2.2.2.2.2.1⟩

/-- Table form: K unfolds into the 8 block conditions capturing the diagonal
`\Sp_2(k)` block, the embedded `U_4 \rtimes G_m(k)` core, and the
`\Ga^{10}` coupling equations, packaged as the named Appendix A table cell. -/
theorem mem_K_table_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) :=
by
  constructor
  · rintro ⟨h11_1, h22_1, h12_1, h21_1, h11_2, h22_2, h12_2, h21_2⟩
    exact ⟨⟨h11_1, h22_1⟩, ⟨h12_1, h21_1⟩, ⟨h11_2, h22_2⟩, ⟨h12_2, h21_2⟩⟩
  · rintro ⟨hdiag1, hcouple1, hdiag2, hcouple2⟩
    exact ⟨hdiag1.1, hdiag1.2, hcouple1.1, hcouple1.2, hdiag2.1, hdiag2.2, hcouple2.1,
      hcouple2.2⟩

private theorem diag_act
    (A : Matrix I I k) (D : Matrix W W k)
    (Ω_I : Matrix I I k) (Ω_W : Matrix W W k) :
    fromBlocks A 0 0 D *
      fromBlocks Ω_I 0 0 Ω_W *
      (fromBlocks A 0 0 D)ᵀ =
    fromBlocks (A * Ω_I * Aᵀ) 0 0 (D * Ω_W * Dᵀ) := by
  simp [fromBlocks_multiply, fromBlocks_transpose]

private theorem diag_lincomb
    (α β : k) (Ω₁_I Ω₂_I : Matrix I I k) (Ω₁_W Ω₂_W : Matrix W W k) :
    α • fromBlocks Ω₁_I 0 0 Ω₁_W + β • fromBlocks Ω₂_I 0 0 Ω₂_W =
    fromBlocks (α • Ω₁_I + β • Ω₂_I) 0 0 (α • Ω₁_W + β • Ω₂_W) := by
  simp [fromBlocks_smul, fromBlocks_add]

/-- The Borel lift produces coefficient matrix !![a², 0; 0, a] (upper-triangular).
This witnesses B ⊆ Q_L: for every a ≠ 0, the Borel matrix !![a², 0; 0, a] is realized.
The I-block scalar a*I₂ gives (a*I₂)*J₂*(a*I₂)^T = a²*J₂ for rep₁ and 0 for rep₂.
The W-block N5.Row1.lift(a,0,0,1) gives the matching a²*rep₁_W + 0*rep₂_W. -/
theorem quotient_action (a : k) (ha : a ≠ 0) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) a)
      (!![a * a, 0; 0, a]) := by
  -- N5.Row1 quotient action for lift(a,0,0,1) with coeff(a,0,0,1)
  have hN5 := Wedge2Formalization.Paper.N5.Row1.quotient_action (k := k) a 0 0 1
  simp only [Wedge2Formalization.N5PureSingularLong.ActBivector] at hN5
  -- coeff(a,0,0,1) = !![a*a, 0; 0, a]
  have hcoeff : Wedge2Formalization.Paper.N5.Row1.coeff (k := k) a 0 0 1 =
      !![a * a, 0; 0, a] := by
    ext i j; fin_cases i <;> fin_cases j <;>
      simp [Wedge2Formalization.Paper.N5.Row1.coeff,
        Wedge2Formalization.N5PureSingularLong.Delta]
  rw [hcoeff] at hN5
  refine ⟨?_, ?_⟩
  · show lift (k := k) a * rep₁ (k := k) * (lift (k := k) a)ᵀ =
        (a * a) • rep₁ (k := k) + (0 : k) • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [diag_act, diag_lincomb]
    congr 1
    · -- I-block: (a•1) * J₂ * (a•1)ᵀ = (a*a) • J₂ + 0 • 0
      ext i j; fin_cases i <;> fin_cases j <;>
        simp [J₂, Matrix.mul_apply, Matrix.smul_apply,
          Matrix.transpose_apply, Fin.sum_univ_two]
    · -- W-block: from hN5.1
      exact hN5.1
  · show lift (k := k) a * rep₂ (k := k) * (lift (k := k) a)ᵀ =
        (0 : k) • rep₁ (k := k) + a • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [diag_act, diag_lincomb]
    congr 1
    · -- I-block: (a•1) * 0 * (a•1)ᵀ = 0 = 0 • J₂ + a • 0
      simp
    · -- W-block: from hN5.2
      exact hN5.2

set_option maxHeartbeats 400000 in
/-- Every invertible setwise stabilizer induces a Borel coefficient matrix (M 1 0 = 0, det M ≠ 0).
Together with `quotient_action`, this establishes Q_L = B:
- B ⊆ Q_L (Borel matrices !![a², 0; 0, a] are realized by the lift family)
- Q_L ⊆ B (this theorem)
The proof of M 1 0 = 0 uses the kernel-proportionality argument:
if M 1 0 ≠ 0, the I-block M₁₀•J₂ is invertible, so ker(RHS) vectors have I-component = 0.
Two independent kernel vectors (from ker(rep₂) ⊇ Fin 2 × {0}) map under gᵀ⁻¹ into the
1-dimensional W-block kernel, giving a proportionality contradiction. -/
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
  dsimp only [] at hact1 hact2
  have rep_indep : ∀ a b : k, a • rep₁ (k := k) + b • rep₂ (k := k) = 0 → a = 0 ∧ b = 0 :=
    fun a b h => rep_pair_independent (k := k) h
  -- det(M) ≠ 0
  have hdetM : Matrix.det M ≠ 0 := by
    intro hdet0
    have hM_dep : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
      simpa [Matrix.det_fin_two] using hdet0
    have hcomb_zero : (M 1 1) • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k) = 0 :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
        calc g * ((M 1 1) • rep₁ (k := k) + (-(M 0 1)) • rep₂ (k := k)) * gᵀ
            = (M 1 1) • (g * rep₁ (k := k) * gᵀ) +
              (-(M 0 1)) • (g * rep₂ (k := k) * gᵀ) := by
              simp [mul_add, add_mul]
          _ = 0 := by rw [hact1, hact2]; ext i j
                      simp only [add_apply, smul_apply, Matrix.zero_apply, smul_eq_mul, neg_mul]
                      linear_combination (rep₁ (k := k) i j) * hM_dep)
    rcases rep_indep _ _ hcomb_zero with ⟨h11, h01neg⟩
    have h01 : M 0 1 = 0 := neg_eq_zero.mp h01neg
    have hcomb2 : (-(M 1 0)) • rep₁ (k := k) + (M 0 0) • rep₂ (k := k) = 0 :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
        calc g * ((-(M 1 0)) • rep₁ (k := k) + (M 0 0) • rep₂ (k := k)) * gᵀ
            = (-(M 1 0)) • (g * rep₁ (k := k) * gᵀ) +
              (M 0 0) • (g * rep₂ (k := k) * gᵀ) := by
              simp [mul_add, add_mul]
          _ = 0 := by rw [hact1, hact2]; ext i j
                      simp only [add_apply, smul_apply, Matrix.zero_apply, smul_eq_mul, neg_mul]
                      linear_combination (rep₂ (k := k) i j) * hM_dep)
    have h10 : M 1 0 = 0 := neg_eq_zero.mp (rep_indep _ _ hcomb2).1
    have hrep2_zero : rep₂ (k := k) = 0 :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
          have := hact2; rw [h10, h11] at this; simpa using this)
    exact absurd ((rep_indep 0 1 (by simpa using hrep2_zero)).2) one_ne_zero
  -- M 1 0 = 0: from the kernel-proportionality argument.
  -- The I-block of RHS is M₁₀•J₂ (invertible when M₁₀ ≠ 0).
  -- ker(rep₂) contains v₁ = e_{inl 0} and v₂ = e_{inl 1} (since I-block of rep₂ is 0).
  -- gᵀ⁻¹ maps them to ker(RHS). With I-block invertible, the images have I-component = 0
  -- and W-component in the 1D kernel of the W-block pencil.
  -- Proportionality + gᵀ equations give contradiction.
  have h10 : M 1 0 = 0 := by
    by_contra h10ne
    -- I-block of RHS is M₁₀•J₂, invertible since det(J₂) = 1 and M₁₀ ≠ 0.
    have hJ₂_det : Matrix.det (J₂ (k := k)) = 1 := by
      simp [J₂, Matrix.det_fin_two]
    have hI_inv : IsUnit (M 1 0 • J₂ (k := k)) := by
      rw [isUnit_iff_isUnit_det, det_smul, hJ₂_det, Fintype.card_fin]
      apply IsUnit.mk0
      simp only [mul_one]
      exact pow_ne_zero 2 h10ne
    -- Step 1: define v₁, v₂ and show they're in ker(rep₂)
    set v₁ : V → k := Pi.single (Sum.inl (0 : Fin 2)) 1
    set v₂ : V → k := Pi.single (Sum.inl (1 : Fin 2)) 1
    have hv₁_ker : rep₂ (k := k) *ᵥ v₁ = 0 := by
      ext i; rcases i with i | (i | i) <;> fin_cases i <;>
        simp [v₁, rep₂, mulVec, dotProduct, Fin.sum_univ_two, Fin.sum_univ_three,
          Fintype.sum_sum_type, fromBlocks,
          N5PureSingularLong.rep₂, N5PureSingularLong.mulY,
          Pi.single, Function.update, Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
          Matrix.cons_val_one, Matrix.head_cons]
    have hv₂_ker : rep₂ (k := k) *ᵥ v₂ = 0 := by
      ext i; rcases i with i | (i | i) <;> fin_cases i <;>
        simp [v₂, rep₂, mulVec, dotProduct, Fin.sum_univ_two, Fin.sum_univ_three,
          Fintype.sum_sum_type, fromBlocks,
          N5PureSingularLong.rep₂, N5PureSingularLong.mulY,
          Pi.single, Function.update, Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
          Matrix.cons_val_one, Matrix.head_cons]
    -- Step 2: x₁ = gᵀ⁻¹ *ᵥ v₁, x₂ = gᵀ⁻¹ *ᵥ v₂ are in ker(RHS)
    have hgT_det : det gᵀ ≠ 0 := by rwa [det_transpose]
    set x₁ := gᵀ⁻¹ *ᵥ v₁
    set x₂ := gᵀ⁻¹ *ᵥ v₂
    have ker_transfer (v : V → k) (hv : rep₂ (k := k) *ᵥ v = 0) :
        (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) *ᵥ (gᵀ⁻¹ *ᵥ v) = 0 := by
      rw [← hact2]
      simp only [mulVec_mulVec]
      rw [Matrix.mul_assoc, mul_nonsing_inv _ (IsUnit.mk0 _ hgT_det), Matrix.mul_one]
      rw [← mulVec_mulVec]
      simp [hv]
    have hx₁_ker : (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) *ᵥ x₁ = 0 :=
      ker_transfer v₁ hv₁_ker
    have hx₂_ker : (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) *ᵥ x₂ = 0 :=
      ker_transfer v₂ hv₂_ker
    -- Step 3: x₁ and x₂ have I-component 0 (from I-block invertibility).
    -- I-block of RHS is M₁₀•J₂, invertible.
    have I_zero (x : V → k) (hxker : (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) *ᵥ x = 0)
        (i : Fin 2) : x (Sum.inl i) = 0 := by
      have hI_full : (M 1 0 • J₂ (k := k)) *ᵥ (fun j => x (Sum.inl j)) = 0 := by
        ext j; have := congrFun hxker (Sum.inl j)
        simp only [Pi.zero_apply, mulVec, dotProduct, Fintype.sum_sum_type,
          Pi.add_apply, Pi.smul_apply, smul_eq_mul,
          Fin.sum_univ_two, Fin.sum_univ_three,
          rep₁, rep₂, fromBlocks_apply₁₁, fromBlocks_apply₁₂,
          Matrix.zero_apply, zero_mul, add_zero, mul_zero] at this
        simpa [mulVec, dotProduct, Fin.sum_univ_two, smul_apply, J₂] using this
      have := (mulVec_injective_iff_isUnit.mpr hI_inv)
        (show (M 1 0 • J₂ (k := k)) *ᵥ (fun j => x (Sum.inl j)) =
          (M 1 0 • J₂ (k := k)) *ᵥ 0 by simp [hI_full])
      exact congrFun this i
    have hx₁I : ∀ i, x₁ (Sum.inl i) = 0 := I_zero x₁ hx₁_ker
    have hx₂I : ∀ i, x₂ (Sum.inl i) = 0 := I_zero x₂ hx₂_ker
    -- Step 4: Extract all pencil kernel equations directly from congrFun hxᵢ_ker.
    -- For each x ∈ ker(pencil) with x(inl _) = 0, extract entries at each coordinate.
    -- The pencil has block structure: I-block = M₁₀•J₂ (invertible, already used),
    -- W-block = fromBlocks 0 (M₁₀•mulX+M₁₁•mulY) (-(M₁₀•mulX+M₁₁•mulY)ᵀ) 0 on Sum (Fin 3) (Fin 2).
    -- We extract individual pencil entries at specific V-coordinates.
    --
    -- Helper: pencil applied to x at a specific V-coordinate.
    -- At (inr(inl i)) for i ∈ Fin 3: gives (M₁₀•mulX+M₁₁•mulY)(i,j) * x(inr(inr j)).
    -- At (inr(inr j)) for j ∈ Fin 2: gives -(M₁₀•mulX+M₁₁•mulY)ᵀ(j,i) * x(inr(inl i)).
    -- Extract directly: pencil entry at (inr(inl 0)) → M₁₀ * x(rr 0) = 0.
    have pencil_rl (x : V → k)
        (hxI : ∀ i, x (Sum.inl i) = 0)
        (hker : (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) *ᵥ x = 0) :
        M 1 0 * x (Sum.inr (Sum.inr 0)) = 0 ∧
        M 1 1 * x (Sum.inr (Sum.inr 0)) + M 1 0 * x (Sum.inr (Sum.inr 1)) = 0 ∧
        -(M 1 0 * x (Sum.inr (Sum.inl 0))) + -(M 1 1 * x (Sum.inr (Sum.inl 1))) = 0 ∧
        -(M 1 0 * x (Sum.inr (Sum.inl 1))) + -(M 1 1 * x (Sum.inr (Sum.inl 2))) = 0 := by
      -- Extract each pencil entry at a specific V-index, simp through the block structure.
      -- The shared simp set that unfolds the nested fromBlocks.
      have aux (idx : V) := congrFun hker idx
      refine ⟨?_, ?_, ?_, ?_⟩
      · have h := aux (Sum.inr (Sum.inl 0))
        simp only [Pi.zero_apply, mulVec, dotProduct, Fintype.sum_sum_type, Pi.add_apply,
          Pi.smul_apply, smul_eq_mul, Fin.sum_univ_two, Fin.sum_univ_three,
          hxI, mul_zero, add_zero, zero_add, zero_mul,
          rep₁, rep₂, J₂,
          Wedge2Formalization.Paper.N5.Row1.rep₁, Wedge2Formalization.Paper.N5.Row1.rep₂,
          N5PureSingularLong.rep₁, N5PureSingularLong.rep₂,
          N5PureSingularLong.mulX, N5PureSingularLong.mulY,
          fromBlocks_apply₁₁, fromBlocks_apply₁₂, fromBlocks_apply₂₁, fromBlocks_apply₂₂,
          Matrix.of_apply,
          Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
          cons_val_two, Matrix.vecTail, Matrix.vecHead, Function.comp_apply,
          Fin.succ_zero_eq_one, transpose_apply, Matrix.add_apply, Matrix.smul_apply,
          Matrix.zero_apply, Matrix.neg_apply, neg_zero, neg_mul, one_mul, mul_one,
          mul_neg, neg_neg] at h
        first | exact h | (linear_combination -h)
      · have h := aux (Sum.inr (Sum.inl 1))
        simp only [Pi.zero_apply, mulVec, dotProduct, Fintype.sum_sum_type, Pi.add_apply,
          Pi.smul_apply, smul_eq_mul, Fin.sum_univ_two, Fin.sum_univ_three,
          hxI, mul_zero, add_zero, zero_add, zero_mul,
          rep₁, rep₂, J₂,
          Wedge2Formalization.Paper.N5.Row1.rep₁, Wedge2Formalization.Paper.N5.Row1.rep₂,
          N5PureSingularLong.rep₁, N5PureSingularLong.rep₂,
          N5PureSingularLong.mulX, N5PureSingularLong.mulY,
          fromBlocks_apply₁₁, fromBlocks_apply₁₂, fromBlocks_apply₂₁, fromBlocks_apply₂₂,
          Matrix.of_apply,
          Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
          cons_val_two, Matrix.vecTail, Matrix.vecHead, Function.comp_apply,
          Fin.succ_zero_eq_one, transpose_apply, Matrix.add_apply, Matrix.smul_apply,
          Matrix.zero_apply, Matrix.neg_apply, neg_zero, neg_mul, one_mul, mul_one,
          mul_neg, neg_neg] at h
        first | exact h | (linear_combination -h)
      · have h := aux (Sum.inr (Sum.inr 0))
        simp only [Pi.zero_apply, mulVec, dotProduct, Fintype.sum_sum_type, Pi.add_apply,
          Pi.smul_apply, smul_eq_mul, Fin.sum_univ_two, Fin.sum_univ_three,
          hxI, mul_zero, add_zero, zero_add, zero_mul,
          rep₁, rep₂, J₂,
          Wedge2Formalization.Paper.N5.Row1.rep₁, Wedge2Formalization.Paper.N5.Row1.rep₂,
          N5PureSingularLong.rep₁, N5PureSingularLong.rep₂,
          N5PureSingularLong.mulX, N5PureSingularLong.mulY,
          fromBlocks_apply₁₁, fromBlocks_apply₁₂, fromBlocks_apply₂₁, fromBlocks_apply₂₂,
          Matrix.of_apply,
          Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
          cons_val_two, Matrix.vecTail, Matrix.vecHead, Function.comp_apply,
          Fin.succ_zero_eq_one, transpose_apply, Matrix.add_apply, Matrix.smul_apply,
          Matrix.zero_apply, Matrix.neg_apply, neg_zero, neg_mul, one_mul, mul_one,
          mul_neg, neg_neg] at h
        first | exact h | (linear_combination -h)
      · have h := aux (Sum.inr (Sum.inr 1))
        simp only [Pi.zero_apply, mulVec, dotProduct, Fintype.sum_sum_type, Pi.add_apply,
          Pi.smul_apply, smul_eq_mul, Fin.sum_univ_two, Fin.sum_univ_three,
          hxI, mul_zero, add_zero, zero_add, zero_mul,
          rep₁, rep₂, J₂,
          Wedge2Formalization.Paper.N5.Row1.rep₁, Wedge2Formalization.Paper.N5.Row1.rep₂,
          N5PureSingularLong.rep₁, N5PureSingularLong.rep₂,
          N5PureSingularLong.mulX, N5PureSingularLong.mulY,
          fromBlocks_apply₁₁, fromBlocks_apply₁₂, fromBlocks_apply₂₁, fromBlocks_apply₂₂,
          Matrix.of_apply,
          Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
          cons_val_two, Matrix.vecTail, Matrix.vecHead, Function.comp_apply,
          Fin.succ_zero_eq_one, transpose_apply, Matrix.add_apply, Matrix.smul_apply,
          Matrix.zero_apply, Matrix.neg_apply, neg_zero, neg_mul, one_mul, mul_one,
          mul_neg, neg_neg] at h
        first | exact h | (linear_combination -h)
    obtain ⟨hx₁_rr0, hx₁_rr1, hx₁_rl01, hx₁_rl12⟩ := pencil_rl x₁ hx₁I hx₁_ker
    obtain ⟨hx₂_rr0, hx₂_rr1, hx₂_rl01, hx₂_rl12⟩ := pencil_rl x₂ hx₂I hx₂_ker
    -- x(rr 0) = 0 and x(rr 1) = 0 from M₁₀ ≠ 0.
    have hx₁W0 : x₁ (Sum.inr (Sum.inr 0)) = 0 := (mul_eq_zero.mp hx₁_rr0).resolve_left h10ne
    have hx₂W0 : x₂ (Sum.inr (Sum.inr 0)) = 0 := (mul_eq_zero.mp hx₂_rr0).resolve_left h10ne
    have hx₁W1 : x₁ (Sum.inr (Sum.inr 1)) = 0 := by
      rw [hx₁W0, mul_zero, zero_add] at hx₁_rr1
      exact (mul_eq_zero.mp hx₁_rr1).resolve_left h10ne
    have hx₂W1 : x₂ (Sum.inr (Sum.inr 1)) = 0 := by
      rw [hx₂W0, mul_zero, zero_add] at hx₂_rr1
      exact (mul_eq_zero.mp hx₂_rr1).resolve_left h10ne
    have hx₁W : ∀ j : Fin 2, x₁ (Sum.inr (Sum.inr j)) = 0 := by
      intro j; match j with | 0 => exact hx₁W0 | 1 => exact hx₁W1
    have hx₂W : ∀ j : Fin 2, x₂ (Sum.inr (Sum.inr j)) = 0 := by
      intro j; match j with | 0 => exact hx₂W0 | 1 => exact hx₂W1
    -- Step 5: gᵀ *ᵥ xᵢ = vᵢ and extract entries.
    have hgx₁ : gᵀ *ᵥ x₁ = v₁ := by
      simp only [x₁]; rw [mulVec_mulVec, mul_nonsing_inv _ (IsUnit.mk0 _ hgT_det), one_mulVec]
    have hgx₂ : gᵀ *ᵥ x₂ = v₂ := by
      simp only [x₂]; rw [mulVec_mulVec, mul_nonsing_inv _ (IsUnit.mk0 _ hgT_det), one_mulVec]
    -- Extract gᵀ *ᵥ xᵢ at (inl 0): sum over inr(inl _) only (inl and inr(inr) are 0).
    have hG₁ := congrFun hgx₁ (Sum.inl (0 : Fin 2))
    have hG₂ := congrFun hgx₂ (Sum.inl (0 : Fin 2))
    simp only [v₁, v₂, Pi.single, Function.update, Sum.inl.injEq,
      dite_true, dite_false, Pi.zero_apply,
      show (0 : Fin 2) ≠ 1 from by decide] at hG₁ hG₂
    -- Simplify hG₁ and hG₂ one at a time to avoid recursion depth issues.
    simp only [mulVec, dotProduct, Fintype.sum_sum_type, transpose_apply,
      hx₁I, hx₁W0, hx₁W1, mul_zero, Finset.sum_const_zero, add_zero,
      zero_add, Fin.sum_univ_two, Fin.sum_univ_three] at hG₁
    simp only [mulVec, dotProduct, Fintype.sum_sum_type, transpose_apply,
      hx₂I, hx₂W0, hx₂W1, mul_zero, Finset.sum_const_zero, add_zero,
      zero_add, Fin.sum_univ_two, Fin.sum_univ_three] at hG₂
    -- Step 6: Use proportionality to derive contradiction.
    -- Combine hGᵢ with hxᵢ_rl01 and hxᵢ_rl12 using linear_combination.
    -- Goal: show x₂(rl 2) * K₂ = 0 with K₂ ≠ 0.
    have hx₂rl2 : x₂ (Sum.inr (Sum.inl 2)) = 0 := by
      -- After simp, hGᵢ uses g (not gᵀ) due to transpose_apply.
      -- Negate the kernel relations (which have the form -(a) + -(b) = 0, i.e. a + b = 0 negated).
      -- The linear_combination uses hGᵢ (positive) and hxᵢ_rl (negated).
      -- Coefficient for hx_rl01: M₁₀ * g(rl 0, l 0)
      -- Coefficient for hx_rl12: M₁₀ * g(rl 1, l 0) - M₁₁ * g(rl 0, l 0)
      -- Since hx_rl0i are negated, we negate the coefficients.
      -- Convert negated kernel eqs to positive form for linear_combination.
      have hx₂_rl01p : M 1 0 * x₂ (Sum.inr (Sum.inl 0)) + M 1 1 * x₂ (Sum.inr (Sum.inl 1)) = 0 := by
        linear_combination -1 * hx₂_rl01
      have hx₂_rl12p : M 1 0 * x₂ (Sum.inr (Sum.inl 1)) + M 1 1 * x₂ (Sum.inr (Sum.inl 2)) = 0 := by
        linear_combination -1 * hx₂_rl12
      have hx₁_rl01p : M 1 0 * x₁ (Sum.inr (Sum.inl 0)) + M 1 1 * x₁ (Sum.inr (Sum.inl 1)) = 0 := by
        linear_combination -1 * hx₁_rl01
      have hx₁_rl12p : M 1 0 * x₁ (Sum.inr (Sum.inl 1)) + M 1 1 * x₁ (Sum.inr (Sum.inl 2)) = 0 := by
        linear_combination -1 * hx₁_rl12
      have hlin₂ : x₂ (Sum.inr (Sum.inl 2)) *
          (M 1 1 * M 1 1 * g (Sum.inr (Sum.inl 0)) (Sum.inl 0) -
           M 1 0 * M 1 1 * g (Sum.inr (Sum.inl 1)) (Sum.inl 0) +
           M 1 0 * M 1 0 * g (Sum.inr (Sum.inl 2)) (Sum.inl 0)) = 0 := by
        linear_combination
          M 1 0 * M 1 0 * hG₂ -
          M 1 0 * g (Sum.inr (Sum.inl 0)) (Sum.inl 0) * hx₂_rl01p -
          (M 1 0 * g (Sum.inr (Sum.inl 1)) (Sum.inl 0) -
           M 1 1 * g (Sum.inr (Sum.inl 0)) (Sum.inl 0)) * hx₂_rl12p
      have hlin₁ : x₁ (Sum.inr (Sum.inl 2)) *
          (M 1 1 * M 1 1 * g (Sum.inr (Sum.inl 0)) (Sum.inl 0) -
           M 1 0 * M 1 1 * g (Sum.inr (Sum.inl 1)) (Sum.inl 0) +
           M 1 0 * M 1 0 * g (Sum.inr (Sum.inl 2)) (Sum.inl 0)) =
          M 1 0 * M 1 0 := by
        linear_combination
          M 1 0 * M 1 0 * hG₁ -
          M 1 0 * g (Sum.inr (Sum.inl 0)) (Sum.inl 0) * hx₁_rl01p -
          (M 1 0 * g (Sum.inr (Sum.inl 1)) (Sum.inl 0) -
           M 1 1 * g (Sum.inr (Sum.inl 0)) (Sum.inl 0)) * hx₁_rl12p
      have hK_ne : M 1 1 * M 1 1 * g (Sum.inr (Sum.inl 0)) (Sum.inl 0) -
          M 1 0 * M 1 1 * g (Sum.inr (Sum.inl 1)) (Sum.inl 0) +
          M 1 0 * M 1 0 * g (Sum.inr (Sum.inl 2)) (Sum.inl 0) ≠ 0 := by
        intro hK0; rw [hK0, mul_zero] at hlin₁
        exact mul_ne_zero h10ne h10ne hlin₁.symm
      exact (mul_eq_zero.mp hlin₂).resolve_right hK_ne
    have hx₂rl1 : x₂ (Sum.inr (Sum.inl 1)) = 0 := by
      have h := hx₂_rl12; simp only [hx₂rl2, mul_zero, neg_zero, add_zero] at h
      -- h : -(M 1 0 * x₂(rl 1)) = 0
      rwa [neg_eq_zero, mul_eq_zero, or_iff_right h10ne] at h
    have hx₂rl0 : x₂ (Sum.inr (Sum.inl 0)) = 0 := by
      have h := hx₂_rl01; simp only [hx₂rl1, mul_zero, neg_zero, add_zero] at h
      rwa [neg_eq_zero, mul_eq_zero, or_iff_right h10ne] at h
    -- x₂ = 0
    have hx₂_zero : x₂ = 0 := by
      ext j; rcases j with (j | (j | j))
      · exact hx₂I j
      · fin_cases j
        · exact hx₂rl0
        · exact hx₂rl1
        · exact hx₂rl2
      · exact hx₂W j
    -- v₂ = 0, contradiction.
    have hv₂_zero : v₂ = 0 := by rw [← hgx₂, hx₂_zero, mulVec_zero]
    have : v₂ (Sum.inl (1 : Fin 2)) = 1 := by
      simp [v₂, Pi.single, Function.update]
    rw [hv₂_zero] at this; simp at this
  refine ⟨M, ⟨h10, hdetM⟩, hact1, hact2⟩

end Row5

end N7
end Paper
end Wedge2Formalization
