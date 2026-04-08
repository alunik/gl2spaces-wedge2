import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N4
import Wedge2Formalization.N3PureSingular
import Wedge2Formalization.N6WeightedTwoPoint
import Mathlib.Tactic.LinearCombination
import Mathlib.LinearAlgebra.Matrix.Rank

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 8.
Representative `S₂ + 2J_{a,1}`.
Divisor `2[a]`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times \Sp_4(k))`, exact quotient
family `Q_L = B` (Borel).

This row is a direct sum of the N3PureSingular pair on a 3-dimensional I-block
and a symplectic pair on the 4-dimensional W-block, where ω₁ has the N4 split form
(e₁∧e₂ + e₃∧e₄) but ω₂ = 0 on the W-block (yielding Sp₄ as the W-block kernel).

Public layer:
- `TableCell` names the Appendix A kernel cell in the mixed `N3PureSingular` /
  split-symplectic block coordinates.
- `Qproj` is the displayed Borel family, with `quotient_action` giving the explicit lift.
- `quotient_image` proves every setwise stabilizer induces a coefficient matrix in that
  Borel quotient family.
-/
namespace Row8

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
    (0 : Matrix W W k)

def paperRep₁ : Matrix V V k := rep₁ (k := k)
def paperRep₂ : Matrix V V k := rep₂ (k := k)
def paperChange : Matrix V V k := 1

/-- The pointwise stabilizer, expressed via block decomposition (I = Fin 3, W = N4.V).
K consists of matrices g = fromBlocks A B C D satisfying:
- (1,1) block conditions: A*ω₁₃*Aᵀ + B*splitRep₁*Bᵀ = ω₁₃ and A*ω₂₃*Aᵀ = ω₂₃
  (the Gm factor from pureSingularShape, with Sp₄ coupling correction)
- (2,2) block conditions: C*ω₁₃*Cᵀ + D*splitRep₁*Dᵀ = splitRep₁ and C*ω₂₃*Cᵀ = 0
  (the Sp₄(k) factor on the W-block)
- (1,2) coupling conditions: A*ω₁₃*Cᵀ + B*splitRep₁*Dᵀ = 0 and A*ω₂₃*Cᵀ = 0
  (the Ga⁶ unipotent radical) -/
def K : Set (Matrix V V k) :=
  { g |
      -- (1,1) of rep₁
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ =
      Wedge2Formalization.N3PureSingular.ω13 (k := k) ∧
      -- (2,2) of rep₁
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ =
      Wedge2Formalization.N4.splitRep₁ (k := k) ∧
      -- (1,2) of rep₁
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ =
      0 ∧
      -- (2,1) of rep₁
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ =
      0 ∧
      -- (1,1) of rep₂
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₁₂ * (0 : Matrix W W k) * g.toBlocks₁₂ᵀ =
      Wedge2Formalization.N3PureSingular.ω23 (k := k) ∧
      -- (2,2) of rep₂
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * (0 : Matrix W W k) * g.toBlocks₂₂ᵀ =
      (0 : Matrix W W k) ∧
      -- (1,2) of rep₂
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * (0 : Matrix W W k) * g.toBlocks₂₂ᵀ =
      0 ∧
      -- (2,1) of rep₂
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * (0 : Matrix W W k) * g.toBlocks₁₂ᵀ =
      0 }

/-- The two diagonal block equations for the first representative, encoding the
`\GL_1(k)`/`\Sp_4(k)` part of Appendix A, row 8. -/
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
        g.toBlocks₁₂ * (0 : Matrix W W k) * g.toBlocks₁₂ᵀ =
      Wedge2Formalization.N3PureSingular.ω23 (k := k) ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₂₂ * (0 : Matrix W W k) * g.toBlocks₂₂ᵀ =
      (0 : Matrix W W k) }

/-- The two coupling equations for the second representative. -/
def CouplingRep₂Cell : Set (Matrix V V k) :=
  { g |
      g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * (0 : Matrix W W k) * g.toBlocks₂₂ᵀ = 0 ∧
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
        g.toBlocks₂₂ * (0 : Matrix W W k) * g.toBlocks₁₂ᵀ = 0 }

/-- Public name for the explicit Appendix A row-8 kernel cell
`K_L = \Ga^6 \rtimes (\Gm(k) \times \Sp_4(k))`,
packaged as the four named block cells in the mixed `N3PureSingular`/split-symplectic
coordinates. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      g ∈ DiagRep₁Cell (k := k) ∧
      g ∈ CouplingRep₁Cell (k := k) ∧
      g ∈ DiagRep₂Cell (k := k) ∧
      g ∈ CouplingRep₂Cell (k := k) }

/-- The quotient family Q_L = B (Borel). -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

/-- The Borel lift: I-block via pureLift(a, 0, 0, a², 1), W-block via torusLift(a,a). -/
def lift (a : k) : Matrix V V k :=
  Matrix.fromBlocks
    (Wedge2Formalization.N3PureSingular.pureLift (k := k) a 0 0 (a * a) 1)
    0
    0
    (Wedge2Formalization.Paper.N4.Row2.torusLift (k := k) a a)

/-- The exact pointwise stabilizer: g fixes the pair iff g ∈ K.
The forward direction decomposes g into I/W blocks (A, B, C, D) and extracts all
8 block equations from the 2 fixing conditions. The backward direction reconstructs
the full matrix equations from the block conditions using fromBlocks algebra. -/
theorem pointwise_stabilizer :
    ∀ g : Matrix V V k,
      (g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
        g * rep₂ (k := k) * gᵀ = rep₂ (k := k)) ↔
      g ∈ K (k := k) := by
  intro g
  have hg_eq : g = fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ :=
    (fromBlocks_toBlocks g).symm
  constructor
  · -- Forward: decompose g into blocks, extract all 8 block conditions
    rintro ⟨h1, h2⟩
    rw [hg_eq] at h1 h2; simp only [rep₁, rep₂] at h1 h2
    exact ⟨by have := congrArg toBlocks₁₁ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₂ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₁₂ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₁ h1; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₁₁ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₂ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₁₂ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this,
           by have := congrArg toBlocks₂₁ h2; simpa [fromBlocks_multiply, fromBlocks_transpose] using this⟩
  · -- Backward: reconstruct the 2 full fixing equations from 8 block conditions
    rintro ⟨h11_1, h22_1, h12_1, h21_1, h11_2, h22_2, h12_2, h21_2⟩
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
          fromBlocks (N3PureSingular.ω23 (k := k)) 0 0 (0 : Matrix W W k) *
          (fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)ᵀ =
        fromBlocks (N3PureSingular.ω23 (k := k)) 0 0 (0 : Matrix W W k) := by
        simp only [Matrix.mul_zero, Matrix.zero_mul, add_zero, zero_add] at h11_2 h22_2 h12_2 h21_2
        simp [fromBlocks_multiply, fromBlocks_transpose, h11_2, h12_2, h21_2, h22_2]
      rw [hg_eq]; exact this

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g * rep₁ (k := k) * gᵀ = rep₁ (k := k) ∧
      g * rep₂ (k := k) * gᵀ = rep₂ (k := k) :=
  (pointwise_stabilizer (k := k) g).symm

/-- Block equations extracted from the fixing condition. The (1,1) and (2,2) blocks
of g * rep_i * g^T = rep_i give these relations. -/
theorem mem_K_diag_blocks (g : Matrix V V k) (hg : g ∈ K (k := k)) :
    g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₁₁ᵀ +
      g.toBlocks₁₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₁₂ᵀ =
    Wedge2Formalization.N3PureSingular.ω13 (k := k) ∧
    g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
      g.toBlocks₂₂ * Wedge2Formalization.N4.splitRep₁ (k := k) * g.toBlocks₂₂ᵀ =
    Wedge2Formalization.N4.splitRep₁ (k := k) ∧
    g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₁₁ᵀ +
      g.toBlocks₁₂ * (0 : Matrix W W k) * g.toBlocks₁₂ᵀ =
    Wedge2Formalization.N3PureSingular.ω23 (k := k) ∧
    g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
      g.toBlocks₂₂ * (0 : Matrix W W k) * g.toBlocks₂₂ᵀ =
    (0 : Matrix W W k) :=
  ⟨hg.1, hg.2.1, hg.2.2.2.2.1, hg.2.2.2.2.2.1⟩

/-- Table form: K unfolds into the 8 block conditions capturing the group structure
of the named Appendix A row-8 table cell. -/
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

/-- The Borel lift produces coefficient matrix !![a, 0; 0, a*a] (upper-triangular).
This witnesses B ⊆ Q_L: for every a ≠ 0, the Borel matrix !![a, 0; 0, a²] is realized. -/
theorem quotient_action (a : k) (ha : a ≠ 0) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) a)
      (!![a, 0; 0, a * a]) := by
  -- N3 pure lift actions
  have hN3_1 := Wedge2Formalization.N3PureSingular.pureLift_act_ω13 (k := k) a 0 0 (a * a) 1
  have hN3_2 := Wedge2Formalization.N3PureSingular.pureLift_act_ω23 (k := k) a 0 0 (a * a) 1
  simp only [Wedge2Formalization.N3PureSingular.ActBivector] at hN3_1 hN3_2
  -- N4 torus action: torusLift(a,a) acts on (splitRep₁, splitRep₂) with coeff !![a, 0; 0, a]
  have hN4 := Wedge2Formalization.Paper.N4.Row2.torus_action (k := k) a a
  simp only [ActsOnOrderedPair, Wedge2Formalization.N4.ActBivector] at hN4
  refine ⟨?_, ?_⟩
  · show lift (k := k) a * rep₁ (k := k) * (lift (k := k) a)ᵀ =
        a • rep₁ (k := k) + (0 : k) • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [diag_act, diag_lincomb]
    congr 1
    · rw [hN3_1]; ring_nf
    · -- W-block: torusLift(a,a)*splitRep₁*torusLift(a,a)ᵀ = a•splitRep₁ + 0•0 = a•splitRep₁
      simp only [smul_zero, add_zero]
      have := hN4.1
      simp only [Wedge2Formalization.Paper.N4.Row2.rep₁,
        Wedge2Formalization.Paper.N4.Row2.rep₂,
        Wedge2Formalization.N4PaperSummary.row2_torusCoeff] at this
      simpa [sub_self, zero_smul, add_zero] using this
  · show lift (k := k) a * rep₂ (k := k) * (lift (k := k) a)ᵀ =
        (0 : k) • rep₁ (k := k) + (a * a) • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [diag_act, diag_lincomb]
    congr 1
    · rw [hN3_2]; ring_nf
    · simp

-- Helper: row i of A * B is in span of rows of B
private lemma row_mul_mem_span
    {m n o : Type*} [Fintype m] [Fintype n] [Fintype o] [DecidableEq o]
    (A : Matrix m n k) (B : Matrix n o k) (i : m) :
    (A * B) i ∈ Submodule.span k (Set.range (fun j : n => B j)) := by
  have key : (A * B) i = ∑ x : n, A i x • B x := by
    ext j; simp [mul_apply, Finset.sum_apply, Pi.smul_apply, smul_eq_mul]
  rw [key]
  exact Submodule.sum_mem _ (fun x _ =>
    Submodule.smul_mem _ _ (Submodule.subset_span ⟨x, rfl⟩))

/-- The determinant of A * B is zero when A : m x n and B : n x m with card n < card m. -/
private theorem det_mul_eq_zero_of_card_lt
    {m n : Type*} [Fintype m] [Fintype n] [DecidableEq m] [DecidableEq n]
    (A : Matrix m n k) (B : Matrix n m k)
    (h : Fintype.card n < Fintype.card m) :
    det (A * B) = 0 := by
  apply det_eq_zero_of_not_linearIndependent_rows
  intro hli
  set S := Submodule.span k (Set.range (fun j : n => B j))
  have hmem : ∀ i, (A * B) i ∈ S := row_mul_mem_span A B
  have hli_sub : LinearIndependent k (fun i : m => ⟨(A * B) i, hmem i⟩ : m → S) :=
    LinearIndependent.of_comp S.subtype hli
  have hfr : Module.finrank k S ≤ Fintype.card n := finrank_range_le_card (fun j : n => B j)
  have hcard : Fintype.card m ≤ Module.finrank k S := hli_sub.fintype_card_le_finrank
  omega

/-- Every invertible setwise stabilizer induces a coefficient matrix in B (Borel). -/
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
  simp only [ActsOnOrderedPair] at hact1 hact2
  refine ⟨M, ?_, hact1, hact2⟩
  simp only [Qproj, Set.mem_setOf_eq]
  have rep_indep : ∀ a b : k, a • rep₁ (k := k) + b • rep₂ (k := k) = 0 → a = 0 ∧ b = 0 := by
    intro a b h
    have h11 := congrArg toBlocks₁₁ h
    simp only [rep₁, rep₂, fromBlocks_smul, fromBlocks_add, toBlocks₁₁, fromBlocks,
        Matrix.of_apply, Sum.elim_inl, smul_zero, add_zero, Matrix.zero_apply] at h11
    have h02 := congrFun (congrFun h11 0) 2
    have h12 := congrFun (congrFun h11 1) 2
    simp [N3PureSingular.ω13, N3PureSingular.ω23, add_apply, smul_apply, Fin.sum_univ_three,
        Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one,
        Matrix.head_cons] at h02 h12
    exact ⟨by linear_combination h02, by linear_combination h12⟩
  -- Part 1: M 1 0 = 0 via (2,2) block determinant argument.
  have h22_rhs : toBlocks₂₂ (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) =
      M 1 0 • Wedge2Formalization.N4.splitRep₁ (k := k) := by
    ext i j
    simp [rep₁, rep₂, toBlocks₂₂, fromBlocks, smul_apply, add_apply]
  have h22_lhs : toBlocks₂₂ (g * rep₂ (k := k) * gᵀ) =
      g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ := by
    ext i j
    simp [rep₂, toBlocks₂₂, toBlocks₂₁, fromBlocks, mul_apply, Fintype.sum_sum_type,
      transpose_apply]
  have h22_eq : g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ =
      M 1 0 • Wedge2Formalization.N4.splitRep₁ (k := k) := by
    rw [← h22_lhs, ← h22_rhs]; exact congrArg toBlocks₂₂ hact2
  have hdet_lhs : det (g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) *
      g.toBlocks₂₁ᵀ) = 0 := by
    have : g.toBlocks₂₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ =
        g.toBlocks₂₁ * (Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ) :=
      Matrix.mul_assoc _ _ _
    rw [this]
    exact det_mul_eq_zero_of_card_lt g.toBlocks₂₁
      (Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ)
      (by simp [W, Wedge2Formalization.N3PureSingular.I, Wedge2Formalization.N4.V,
          Wedge2Formalization.N4.I])
  have hdet_rhs : det (M 1 0 • Wedge2Formalization.N4.splitRep₁ (k := k)) = (M 1 0) ^ 4 := by
    rw [det_smul]
    simp [Wedge2Formalization.N4.V, Wedge2Formalization.N4.I,
      Wedge2Formalization.N6WeightedTwoPoint.det_splitRep₁]
  have hM10_pow : (M 1 0) ^ 4 = 0 := by rw [← hdet_rhs, ← h22_eq, hdet_lhs]
  have h10 : M 1 0 = 0 := by
    have h4 := hM10_pow
    rw [show (4 : ℕ) = 2 + 2 from rfl, pow_add] at h4
    rcases mul_eq_zero.mp h4 with h2 | h2 <;> {
      rw [sq] at h2
      exact (mul_eq_zero.mp h2).elim id id
    }
  -- Part 2: det M ≠ 0 via standard argument.
  constructor
  · exact h10
  · intro hdet0
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
    have h10' : M 1 0 = 0 := neg_eq_zero.mp (rep_indep _ _ hcomb2).1
    have hrep2_zero : rep₂ (k := k) = 0 :=
      (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 (by
          have := hact2; rw [h10', h11] at this; simpa using this)
    exact absurd ((rep_indep 0 1 (by simpa using hrep2_zero)).2) one_ne_zero

end Row8

end N7
end Paper
end Wedge2Formalization
