import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N4
import Wedge2Formalization.N3PureSingular
import Mathlib.Tactic.LinearCombination

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 9.
Representative `S₂ + J_{a,1} + J_{b,1}`.
Divisor `[a]+[b]`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times (\SL_2(k) \times \SL_2(k)))`, exact quotient
family `Q_L = N_T` (normalizer of the split torus).

This row is a direct sum of the N3PureSingular pair on a 3-dimensional I-block
and the N4 split-support pair on a 4-dimensional W-block.

Public layer:
- `TableCell` names the Appendix A kernel cell in the mixed `N3PureSingular` /
  split-support block coordinates.
- `Qproj` is the displayed normalizer-of-torus family, with both torus and swap lifts
  exposed publicly.
- `quotient_image` proves every setwise stabilizer induces a coefficient matrix in that
  `N_T` quotient family.
-/
namespace Row9

abbrev I := Wedge2Formalization.N3PureSingular.I  -- Fin 3
abbrev W := Wedge2Formalization.N4.V              -- Fin 2 ⊕ Fin 2
abbrev V := Sum I W

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N3PureSingular.ω13 (k := k))
    0
    0
    (Wedge2Formalization.N4.splitRep₁ (k := k))

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N3PureSingular.ω23 (k := k))
    0
    0
    (Wedge2Formalization.N4.splitRep₂ (k := k))

def paperRep₁ : Matrix V V k := rep₁ (k := k)
def paperRep₂ : Matrix V V k := rep₂ (k := k)
def paperChange : Matrix V V k := 1

/-- The pointwise stabilizer, expressed via block decomposition (I = Fin 3, W = N4.V).
K consists of matrices g = fromBlocks A B C D satisfying all 8 block conditions from
g * rep_i * gᵀ = rep_i. The (1,1) blocks encode the Gm factor, the (2,2) blocks encode
the SL₂ x SL₂ factor, and the off-diagonal blocks encode the Ga⁶ unipotent radical. -/
def K : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ =
      Wedge2Formalization.N3PureSingular.ω13 (k := k) ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.N4.splitRep₁ (k := k) ∧
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ = 0 ∧
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₁₂ᵀ =
      Wedge2Formalization.N3PureSingular.ω23 (k := k) ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.N4.splitRep₂ (k := k) ∧
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₁₂ᵀ = 0 }

/-- The two diagonal block equations for the first representative, encoding the
`\GL_1(k)`/`\SL_2(k) \times \SL_2(k)` part of Appendix A, row 9. -/
def DiagRep₁Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ =
      Wedge2Formalization.N3PureSingular.ω13 (k := k) ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.N4.splitRep₁ (k := k) }

/-- The two coupling equations for the first representative. -/
def CouplingRep₁Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ = 0 }

/-- The two diagonal block equations for the second representative. -/
def DiagRep₂Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₁₂ᵀ =
      Wedge2Formalization.N3PureSingular.ω23 (k := k) ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.N4.splitRep₂ (k := k) }

/-- The two coupling equations for the second representative. -/
def CouplingRep₂Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₁₂ᵀ = 0 }

/-- Public name for the explicit Appendix A row-9 kernel cell
`K_L = \Ga^6 \rtimes (\Gm(k) \times (\SL_2(k) \times \SL_2(k)))`,
packaged as the four named block cells in the mixed `N3PureSingular`/split-support
coordinates. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      g ∈ DiagRep₁Cell (k := k) ∧
      g ∈ CouplingRep₁Cell (k := k) ∧
      g ∈ DiagRep₂Cell (k := k) ∧
      g ∈ CouplingRep₂Cell (k := k) }

/-- The normalizer of the split torus N_T. Every setwise stabilizer induces a
coefficient matrix in this family. The condition M 1 0 * (M 1 0 + M 1 1) = 0 means
the second row sends rep₂ to one of the two rank-drop lines of the N4 split pencil.
Combined with det M ≠ 0, this gives the N_T condition.
- Torus coset: `!![u, v-u; 0, v]` (from the torus lift) has M 1 0 = 0
- Swap coset: `!![u, v-u; u, -u]` (from the swap lift) has M 1 0 + M 1 1 = 0
Both satisfy M 1 0 * (M 1 0 + M 1 1) = 0.
`quotient_action` and `quotient_action_swap` witness the lower bound N_T ⊆ Q_L.
`quotient_image` proves the upper bound Q_L ⊆ N_T via the kernel-proportionality argument. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 * (M 1 0 + M 1 1) = 0 ∧ Matrix.det M ≠ 0 }

/-- The torus lift acts on the I-block via `pureLift u (v-u) 0 v 1` and on the
W-block via `N4.Row2.torusLift u v`, producing coefficient matrix `!![u, v-u; 0, v]`. -/
def lift (u v : k) : Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N3PureSingular.pureLift (k := k) u (v - u) 0 v 1)
    0
    0
    (Wedge2Formalization.Paper.N4.Row2.torusLift (k := k) u v)

/-- The swap-coset lift acts on the I-block via `pureLift u (v-u) u (-u) 1` and on the
W-block via `N4.Row2.swapLift u v`, producing coefficient matrix `!![u, v-u; u, -u]`. -/
def swapLift (u v : k) : Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N3PureSingular.pureLift (k := k) u (v - u) u (-u) 1)
    0
    0
    (Wedge2Formalization.Paper.N4.Row2.swapLift (k := k) u v)

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
          fromBlocks (N3PureSingular.ω13 (k := k)) 0 0 (N4.splitRep₁ (k := k)) *
          (fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
        fromBlocks (N3PureSingular.ω13 (k := k)) 0 0 (N4.splitRep₁ (k := k)) := by
        simp [fromBlocks_multiply, fromBlocks_transpose, h11_1, h12_1, h21_1, h22_1]
      rw [hg_eq]; exact this
    · show g * rep₂ (k := k) * gᵀ = rep₂ (k := k)
      have : fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ *
          fromBlocks (N3PureSingular.ω23 (k := k)) 0 0 (N4.splitRep₂ (k := k)) *
          (fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
        fromBlocks (N3PureSingular.ω23 (k := k)) 0 0 (N4.splitRep₂ (k := k)) := by
        simp [fromBlocks_multiply, fromBlocks_transpose, h11_2, h12_2, h21_2, h22_2]
      rw [hg_eq]; exact this

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
      g * rep₂ (k := k) * gᵀ = rep₂ (k := k) :=
  (pointwise_stabilizer (k := k) g).symm

/-- The (2,2) diagonal block of any fixing element preserves the N4 split pair,
and the (1,1) block preserves the N3PureSingular pair (modulo coupling corrections).
For the full fixing set, the diagonal blocks satisfy:
- `A * ω₁₃ * Aᵀ + B * splitRep₁ * Bᵀ = ω₁₃`  (I-block of rep₁ equation)
- `D * splitRep₁ * Dᵀ + C * ω₁₃ * Cᵀ = splitRep₁`  (W-block of rep₁ equation)
When B = C = 0, these reduce to A fixing ω₁₃/ω₂₃ and D fixing splitRep₁/splitRep₂. -/
theorem mem_K_diag_blocks (g : Matrix V V k) (hg : g ∈ K (k := k)) :
    g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
      g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ =
    Wedge2Formalization.N3PureSingular.ω13 (k := k) ∧
    g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
      g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ =
    Wedge2Formalization.N4.splitRep₁ (k := k) ∧
    g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
      g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₁₂ᵀ =
    Wedge2Formalization.N3PureSingular.ω23 (k := k) ∧
    g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
      g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₂ (k := k) * g.toBlocks₂₂ᵀ =
    Wedge2Formalization.N4.splitRep₂ (k := k) :=
  ⟨hg.1, hg.2.1, hg.2.2.2.2.1, hg.2.2.2.2.2.1⟩

/-- Table form: K unfolds into the 8 block conditions capturing
the named Appendix A row-9 table cell. -/
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

theorem paperRep₁_transport :
    (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (paperRep₁ (k := k)) (paperChange (k := k)) = rep₁ (k := k) := by
  simp [paperRep₁, paperChange]

theorem paperRep₂_transport :
    (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (paperRep₂ (k := k)) (paperChange (k := k)) = rep₂ (k := k) := by
  simp [paperRep₂, paperChange]

/-- The reps are linearly independent (from the (1,1) block N3 structure). -/
theorem rep_pair_independent {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) : a = 0 ∧ b = 0 := by
  have h11 := congrArg toBlocks₁₁ h
  simp only [rep₁, rep₂, fromBlocks_smul, fromBlocks_add, toBlocks₁₁, fromBlocks,
    Matrix.of_apply, Sum.elim_inl, smul_zero, add_zero, Matrix.zero_apply] at h11
  have h02 := congrFun (congrFun h11 0) 2
  have h12 := congrFun (congrFun h11 1) 2
  simp [N3PureSingular.ω13, N3PureSingular.ω23, add_apply, smul_apply, Fin.sum_univ_three,
    Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.head_cons] at h02 h12
  exact ⟨by linear_combination h02, by linear_combination h12⟩

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

/-- The torus lift produces coefficient matrix `!![u, v-u; 0, v]` (torus coset of N_T).
Both blocks act with the same coefficient matrix: I-block via pureLift, W-block via
N4.Row2.torusLift. -/
theorem quotient_action (u v : k) (hu : u ≠ 0) (hv : v ≠ 0) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) u v)
      (!![u, v - u; 0, v]) := by
  -- N3 pure lift actions
  have hN3_1 := Wedge2Formalization.N3PureSingular.pureLift_act_ω13 (k := k) u (v-u) 0 v 1
  have hN3_2 := Wedge2Formalization.N3PureSingular.pureLift_act_ω23 (k := k) u (v-u) 0 v 1
  simp only [Wedge2Formalization.N3PureSingular.ActBivector] at hN3_1 hN3_2
  -- N4 torus action
  have hN4 := Wedge2Formalization.Paper.N4.Row2.torus_action (k := k) u v
  simp only [ActsOnOrderedPair, Wedge2Formalization.N4.ActBivector,
    Wedge2Formalization.N4PaperSummary.row2_torusCoeff,
    Wedge2Formalization.Paper.N4.Row2.rep₁, Wedge2Formalization.Paper.N4.Row2.rep₂] at hN4
  refine ⟨?_, ?_⟩
  · show lift (k := k) u v * rep₁ (k := k) * (lift (k := k) u v)ᵀ =
        u • rep₁ (k := k) + (v - u) • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [diag_act, diag_lincomb]
    congr 1
    · rw [hN3_1]; ring_nf
    · exact hN4.1
  · show lift (k := k) u v * rep₂ (k := k) * (lift (k := k) u v)ᵀ =
        (0 : k) • rep₁ (k := k) + v • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [diag_act, diag_lincomb]
    congr 1
    · rw [hN3_2]; ring_nf
    · exact hN4.2

/-- The swap-coset lift produces coefficient matrix `!![u, v-u; u, -u]` (swap coset of N_T).
Both blocks act with the same coefficient matrix: I-block via pureLift, W-block via
N4.Row2.swapLift. -/
theorem quotient_action_swap (u v : k) (hu : u ≠ 0) (hv : v ≠ 0) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (swapLift (k := k) u v)
      (!![u, v - u; u, -u]) := by
  -- N3 pure lift actions
  have hN3_1 := Wedge2Formalization.N3PureSingular.pureLift_act_ω13 (k := k) u (v-u) u (-u) 1
  have hN3_2 := Wedge2Formalization.N3PureSingular.pureLift_act_ω23 (k := k) u (v-u) u (-u) 1
  simp only [Wedge2Formalization.N3PureSingular.ActBivector] at hN3_1 hN3_2
  -- N4 swap action
  have hN4 := Wedge2Formalization.Paper.N4.Row2.swap_action (k := k) u v
  simp only [ActsOnOrderedPair, Wedge2Formalization.N4.ActBivector,
    Wedge2Formalization.N4PaperSummary.row2_swapCoeff,
    Wedge2Formalization.Paper.N4.Row2.rep₁, Wedge2Formalization.Paper.N4.Row2.rep₂] at hN4
  refine ⟨?_, ?_⟩
  · show swapLift (k := k) u v * rep₁ (k := k) * (swapLift (k := k) u v)ᵀ =
        u • rep₁ (k := k) + (v - u) • rep₂ (k := k)
    simp only [rep₁, rep₂, swapLift]
    rw [diag_act, diag_lincomb]
    congr 1
    · rw [hN3_1]; ring_nf
    · exact hN4.1
  · show swapLift (k := k) u v * rep₂ (k := k) * (swapLift (k := k) u v)ᵀ =
        u • rep₁ (k := k) + (-u) • rep₂ (k := k)
    simp only [rep₁, rep₂, swapLift]
    rw [diag_act, diag_lincomb]
    congr 1
    · rw [hN3_2]; ring_nf
    · exact hN4.2

/-- The torus coefficient matrix belongs to N_T: M 1 0 = 0 so M 1 0 * (M 1 0 + M 1 1) = 0. -/
theorem torusCoeff_mem_Qproj (u v : k) (hu : u ≠ 0) (hv : v ≠ 0) :
    !![u, v - u; (0 : k), v] ∈ Qproj (k := k) := by
  constructor
  · -- M 1 0 * (M 1 0 + M 1 1) = 0 * (0 + v) = 0
    simp [Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.head_cons]
  · show Matrix.det !![u, v - u; (0 : k), v] ≠ 0
    simp only [Matrix.det_fin_two, Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.head_cons, Matrix.head_fin_const, mul_zero, sub_zero]
    exact mul_ne_zero hu hv

/-- The swap coefficient matrix belongs to N_T: M 1 0 + M 1 1 = u + (-u) = 0
so M 1 0 * (M 1 0 + M 1 1) = 0. -/
theorem swapCoeff_mem_Qproj (u v : k) (hu : u ≠ 0) (hv : v ≠ 0) :
    !![u, v - u; u, -u] ∈ Qproj (k := k) := by
  constructor
  · -- M 1 0 * (M 1 0 + M 1 1) = u * (u + (-u)) = u * 0 = 0
    simp [Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.head_cons]
  · show Matrix.det !![u, v - u; u, -u] ≠ 0
    simp only [Matrix.det_fin_two, Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
      Matrix.cons_val_one, Matrix.head_cons, Matrix.head_fin_const]
    intro h
    have hvu : v * u = 0 := by linear_combination -h
    exact mul_ne_zero hv hu hvu

set_option maxHeartbeats 400000 in
/-- Every invertible setwise stabilizer induces a coefficient matrix in N_T
(M 1 0 * (M 1 0 + M 1 1) = 0 ∧ det M ≠ 0). Together with `quotient_action` and
`quotient_action_swap`, this establishes Q_L = N_T.
The proof of M 1 0 * (M 1 0 + M 1 1) = 0 uses the kernel-proportionality argument:
when M₁₀*(M₁₀+M₁₁) ≠ 0, the W-block M₁₀•splitRep₁+M₁₁•splitRep₂ has det = (M₁₀*(M₁₀+M₁₁))² ≠ 0,
so ker(RHS) vectors have W-component = 0. The I-block M₁₀•ω₁₃+M₁₁•ω₂₃ has 1D kernel (M₁₀ ≠ 0).
Two independent kernel vectors from ker(rep₂) (in the first Fin 2 block of the W-space)
are forced into this 1D kernel, giving a proportionality contradiction. -/
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
  -- M 1 0 * (M 1 0 + M 1 1) = 0: from the kernel-proportionality argument.
  -- When M₁₀*(M₁₀+M₁₁) ≠ 0, the W-block det = (M₁₀*(M₁₀+M₁₁))² ≠ 0.
  -- ker(rep₂) contains v₁ = e_{inr(inl 0)} and v₂ = e_{inr(inl 1)}.
  -- gᵀ⁻¹ maps them to ker(RHS). W-block invertible ⟹ W-component = 0.
  -- I-block 1D kernel (M₁₀ ≠ 0) ⟹ proportional ⟹ contradiction.
  have hNT : M 1 0 * (M 1 0 + M 1 1) = 0 := by
    by_contra hNTne
    -- M₁₀ ≠ 0 follows from the product being nonzero
    have h10ne : M 1 0 ≠ 0 := left_ne_zero_of_mul hNTne
    -- W-block det ≠ 0
    have hW_det : det (M 1 0 • N4.splitRep₁ (k := k) + M 1 1 • N4.splitRep₂ (k := k)) ≠ 0 := by
      rw [show det (M 1 0 • N4.splitRep₁ (k := k) + M 1 1 • N4.splitRep₂ (k := k)) =
        (M 1 0 * (M 1 0 + M 1 1)) * (M 1 0 * (M 1 0 + M 1 1)) from
          Wedge2Formalization.N4Summary.split_det_in_basis (k := k) (M 1 0) (M 1 1)]
      exact mul_ne_zero hNTne hNTne
    -- Step 1: define v₁, v₂ and show they're in ker(rep₂)
    -- v₁ = e_{inr(inl 0)}, v₂ = e_{inr(inl 1)} are in ker(rep₂) since
    -- splitRep₂ = ω34 = fromBlocks 0 0 0 J, and these are in the first Fin 2 block.
    set v₁ : V → k := Pi.single (Sum.inr (Sum.inl (0 : Fin 2))) 1
    set v₂ : V → k := Pi.single (Sum.inr (Sum.inl (1 : Fin 2))) 1
    have hv₁_ker : rep₂ (k := k) *ᵥ v₁ = 0 := by
      ext i; rcases i with i | (i | i) <;> fin_cases i <;>
        simp [v₁, rep₂, mulVec, dotProduct, Fin.sum_univ_two, Fin.sum_univ_three,
          Fintype.sum_sum_type, fromBlocks, N3PureSingular.ω23,
          N4.splitRep₂, N4.ω34, N4.J,
          Pi.single, Function.update, Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero,
          Matrix.cons_val_one, Matrix.head_cons]
    have hv₂_ker : rep₂ (k := k) *ᵥ v₂ = 0 := by
      ext i; rcases i with i | (i | i) <;> fin_cases i <;>
        simp [v₂, rep₂, mulVec, dotProduct, Fin.sum_univ_two, Fin.sum_univ_three,
          Fintype.sum_sum_type, fromBlocks, N3PureSingular.ω23,
          N4.splitRep₂, N4.ω34, N4.J,
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
    -- Step 3: x₁ and x₂ have W-component 0 (from W-block invertibility).
    have hW_inv : IsUnit (M 1 0 • N4.splitRep₁ (k := k) + M 1 1 • N4.splitRep₂ (k := k)) := by
      rw [isUnit_iff_isUnit_det]
      apply IsUnit.mk0
      rw [Wedge2Formalization.N4Summary.split_det_in_basis]
      exact mul_ne_zero hNTne hNTne
    -- W-zero: extract the W-block equation from the full pencil equation.
    -- The (2,2) block of RHS gives (M₁₀•splitRep₁+M₁₁•splitRep₂) *ᵥ (x ∘ inr) = something.
    -- But since the cross-blocks of rep₁/rep₂ are 0 and x(inl _) doesn't contribute to the
    -- (inr, _) entries, the W-block equation follows.
    -- We prove this by showing the full pencil applied to x has the W-block equation.
    have W_zero (x : V → k) (hx : (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) *ᵥ x = 0)
        (j : W) : x (Sum.inr j) = 0 := by
      have hW_full : (M 1 0 • N4.splitRep₁ (k := k) + M 1 1 • N4.splitRep₂ (k := k)) *ᵥ
          (fun i => x (Sum.inr i)) = 0 := by
        ext j'
        have h := congrFun hx (Sum.inr j')
        simp only [Pi.zero_apply, mulVec, dotProduct, Fintype.sum_sum_type,
          rep₁, rep₂, fromBlocks_apply₂₁, fromBlocks_apply₂₂,
          add_apply, smul_apply, smul_eq_mul, mul_zero, Finset.sum_const_zero, zero_add,
          Matrix.zero_apply, zero_mul, add_zero] at h
        simpa [mulVec, dotProduct, Fintype.sum_sum_type, add_apply, smul_apply, smul_eq_mul] using h
      have := (mulVec_injective_iff_isUnit.mpr hW_inv)
        (show (M 1 0 • N4.splitRep₁ (k := k) + M 1 1 • N4.splitRep₂ (k := k)) *ᵥ
          (fun i => x (Sum.inr i)) =
        (M 1 0 • N4.splitRep₁ (k := k) + M 1 1 • N4.splitRep₂ (k := k)) *ᵥ 0 by
          simp [hW_full])
      exact congrFun this j
    have hx₁W : ∀ j, x₁ (Sum.inr j) = 0 := W_zero x₁ hx₁_ker
    have hx₂W : ∀ j, x₂ (Sum.inr j) = 0 := W_zero x₂ hx₂_ker
    -- Step 4: I-block kernel analysis.
    -- The I-block pencil is M₁₀•ω₁₃+M₁₁•ω₂₃ on Fin 3.
    -- Entry at (inl 0): M₁₀*x(l2) = 0 ⟹ x(l2) = 0 (M₁₀ ≠ 0).
    -- Entry at (inl 2): -M₁₀*x(l0)-M₁₁*x(l1) = 0.
    have pencil_entry_I (x : V → k) (hxW : ∀ j, x (Sum.inr j) = 0)
        (hker : (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) *ᵥ x = 0)
        (i : Fin 3) :
        M 1 0 * (N3PureSingular.ω13 (k := k) *ᵥ (fun j => x (Sum.inl j))) i +
        M 1 1 * (N3PureSingular.ω23 (k := k) *ᵥ (fun j => x (Sum.inl j))) i = 0 := by
      have h := congrFun hker (Sum.inl i)
      simp only [Pi.zero_apply, mulVec, dotProduct, Fintype.sum_sum_type,
        Pi.add_apply, Pi.smul_apply, smul_eq_mul,
        Fin.sum_univ_three, Fin.sum_univ_two,
        hxW, mul_zero, add_zero, zero_add,
        rep₁, rep₂, fromBlocks_apply₁₁, fromBlocks_apply₁₂,
        Matrix.zero_apply, zero_mul] at h
      convert h using 1
      simp [mulVec, dotProduct, Fin.sum_univ_three]
      ring
    -- Extract I-kernel entries.
    have hx₁l2 : x₁ (Sum.inl 2) = 0 := by
      have h := pencil_entry_I x₁ hx₁W hx₁_ker 0
      simp only [mulVec, dotProduct, Fin.sum_univ_three,
        N3PureSingular.ω13, N3PureSingular.ω23, Matrix.of_apply,
        Matrix.cons_val_zero, Matrix.cons_val_one, cons_val_two,
        Matrix.head_cons, Matrix.vecTail, Matrix.vecHead,
        Function.comp_apply, Fin.succ_zero_eq_one,
        zero_mul, mul_one, mul_zero, add_zero, zero_add, neg_zero, neg_mul, one_mul] at h
      exact (mul_eq_zero.mp h).resolve_left h10ne
    have hx₂l2 : x₂ (Sum.inl 2) = 0 := by
      have h := pencil_entry_I x₂ hx₂W hx₂_ker 0
      simp only [mulVec, dotProduct, Fin.sum_univ_three,
        N3PureSingular.ω13, N3PureSingular.ω23, Matrix.of_apply,
        Matrix.cons_val_zero, Matrix.cons_val_one, cons_val_two,
        Matrix.head_cons, Matrix.vecTail, Matrix.vecHead,
        Function.comp_apply, Fin.succ_zero_eq_one,
        zero_mul, mul_one, mul_zero, add_zero, zero_add, neg_zero, neg_mul, one_mul] at h
      exact (mul_eq_zero.mp h).resolve_left h10ne
    have hIker₁_2 : M 1 0 * x₁ (Sum.inl 0) + M 1 1 * x₁ (Sum.inl 1) = 0 := by
      have h := pencil_entry_I x₁ hx₁W hx₁_ker 2
      simp only [mulVec, dotProduct, Fin.sum_univ_three,
        N3PureSingular.ω13, N3PureSingular.ω23, Matrix.of_apply,
        Matrix.cons_val_zero, Matrix.cons_val_one, cons_val_two,
        Matrix.head_cons, Matrix.vecTail, Matrix.vecHead,
        Function.comp_apply, Fin.succ_zero_eq_one,
        zero_mul, mul_one, mul_zero, add_zero, zero_add, neg_mul, one_mul, neg_zero, mul_neg] at h
      linear_combination -h
    have hIker₂_2 : M 1 0 * x₂ (Sum.inl 0) + M 1 1 * x₂ (Sum.inl 1) = 0 := by
      have h := pencil_entry_I x₂ hx₂W hx₂_ker 2
      simp only [mulVec, dotProduct, Fin.sum_univ_three,
        N3PureSingular.ω13, N3PureSingular.ω23, Matrix.of_apply,
        Matrix.cons_val_zero, Matrix.cons_val_one, cons_val_two,
        Matrix.head_cons, Matrix.vecTail, Matrix.vecHead,
        Function.comp_apply, Fin.succ_zero_eq_one,
        zero_mul, mul_one, mul_zero, add_zero, zero_add, neg_mul, one_mul, neg_zero, mul_neg] at h
      linear_combination -h
    -- Step 5: Extract gᵀ *ᵥ xᵢ = vᵢ at (inr(inl 0)).
    have hgx₁ : gᵀ *ᵥ x₁ = v₁ := by
      simp only [x₁]; rw [mulVec_mulVec, mul_nonsing_inv _ (IsUnit.mk0 _ hgT_det), one_mulVec]
    have hgx₂ : gᵀ *ᵥ x₂ = v₂ := by
      simp only [x₂]; rw [mulVec_mulVec, mul_nonsing_inv _ (IsUnit.mk0 _ hgT_det), one_mulVec]
    have hG₁ := congrFun hgx₁ (Sum.inr (Sum.inl (0 : Fin 2)))
    have hG₂ := congrFun hgx₂ (Sum.inr (Sum.inl (0 : Fin 2)))
    simp only [v₁, v₂, Pi.single, Function.update, Sum.inr.injEq, Sum.inl.injEq,
      dite_true, dite_false, if_true, if_false, Pi.zero_apply,
      show (0 : Fin 2) ≠ 1 from by decide,
      show (1 : Fin 2) ≠ 0 from by decide] at hG₁ hG₂
    simp only [mulVec, dotProduct, Fintype.sum_sum_type, transpose_apply,
      hx₁W, hx₂W, hx₁l2, hx₂l2, mul_zero, Finset.sum_const_zero, add_zero,
      Fin.sum_univ_three, Fin.sum_univ_two] at hG₁ hG₂
    -- Step 6: Factor out via I-kernel relation and derive contradiction.
    -- Same technique as Row 7: linear_combination to get x₂(l1) * K = 0 with K ≠ 0.
    have hx₂l1 : x₂ (Sum.inl 1) = 0 := by
      have hlin₂ : x₂ (Sum.inl 1) *
          (-(M 1 1) * g (Sum.inl 0) (Sum.inr (Sum.inl (0 : Fin 2))) +
           M 1 0 * g (Sum.inl 1) (Sum.inr (Sum.inl (0 : Fin 2)))) = 0 := by
        linear_combination
          M 1 0 * hG₂ - g (Sum.inl 0) (Sum.inr (Sum.inl (0 : Fin 2))) * hIker₂_2
      have hlin₁ : x₁ (Sum.inl 1) *
          (-(M 1 1) * g (Sum.inl 0) (Sum.inr (Sum.inl (0 : Fin 2))) +
           M 1 0 * g (Sum.inl 1) (Sum.inr (Sum.inl (0 : Fin 2)))) = M 1 0 := by
        linear_combination
          M 1 0 * hG₁ - g (Sum.inl 0) (Sum.inr (Sum.inl (0 : Fin 2))) * hIker₁_2
      have hK_ne : -(M 1 1) * g (Sum.inl 0) (Sum.inr (Sum.inl (0 : Fin 2))) +
          M 1 0 * g (Sum.inl 1) (Sum.inr (Sum.inl (0 : Fin 2))) ≠ 0 := by
        intro hK0; rw [hK0, mul_zero] at hlin₁; exact h10ne hlin₁.symm
      exact (mul_eq_zero.mp hlin₂).resolve_right hK_ne
    have hx₂l0 : x₂ (Sum.inl 0) = 0 := by
      have : M 1 0 * x₂ (Sum.inl 0) + M 1 1 * x₂ (Sum.inl 1) = 0 := hIker₂_2
      rw [hx₂l1, mul_zero, add_zero] at this
      exact (mul_eq_zero.mp this).resolve_left h10ne
    have hx₂_zero : x₂ = 0 := by
      ext j; rcases j with (j | j)
      · fin_cases j
        · exact hx₂l0
        · exact hx₂l1
        · exact hx₂l2
      · exact hx₂W j
    have hv₂_zero : v₂ = 0 := by rw [← hgx₂, hx₂_zero, mulVec_zero]
    have : v₂ (Sum.inr (Sum.inl (1 : Fin 2))) = 1 := by
      simp [v₂, Pi.single, Function.update]
    rw [hv₂_zero] at this; simp at this
  refine ⟨M, ⟨hNT, hdetM⟩, hact1, hact2⟩

end Row9

end N7
end Paper
end Wedge2Formalization
