import Wedge2Formalization.Paper.Core
import Wedge2Formalization.Paper.N6.Row11

namespace Wedge2Formalization
namespace Paper
namespace N7

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 7`, row 3.
Representative `S_2^2 + S_1`.
Divisor `0`.
Claimed stabilizer:
`K_L = \Ga^6 \rtimes (\Gm(k) \times (U_6 \rtimes \GL_2(k)))`, exact quotient
family `Q_L = \GL_2(k)`.

This row is a radical-1 extension of N6 Row 11 (DoublePureSingular).
Note: N6.Row11 reps have det = 0 (pure singular), so the B = 0 argument requires
the pair-fixing condition rather than a single det argument.
-/
namespace Row3

abbrev I := Fin 1
abbrev W := Wedge2Formalization.N6DoublePureSingular.V
abbrev V := Sum I W

def scalarBlock (a : k) : Matrix I I k := !![a]

def rep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k))

def rep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k))

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

/-- The reps are linearly independent (from the W-block N6.Row11 independence). -/
theorem rep_pair_independent {a b : k}
    (h : a • rep₁ (k := k) + b • rep₂ (k := k) = 0) : a = 0 ∧ b = 0 := by
  have h22 := congrArg toBlocks₂₂ h
  simp only [rep₁, rep₂, fromBlocks_smul, fromBlocks_add, toBlocks₂₂, fromBlocks,
    Matrix.of_apply, Sum.elim_inr, Matrix.zero_apply, Matrix.smul_apply,
    Matrix.add_apply, smul_zero, add_zero] at h22
  exact Wedge2Formalization.Paper.N6.Row11.rep_pair_independent (k := k) h22

/-- Public name for the embedded Appendix A, `n = 6`, row 11 core
`U_6 \rtimes \GL_2(k)` on the lower-right `6`-space. -/
def DoublePureCore : Set (Matrix W W k) :=
  Wedge2Formalization.Paper.N6.Row11.K (k := k)

def K : Set (Matrix V V k) :=
  { g |
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ DoublePureCore (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D }

def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  Wedge2Formalization.Paper.N6.Row11.Qproj (k := k)

-- The GL₂ lift embeds the N6 quotient representative (fromBlocks G 0 0 G) in the W-block.
-- Note: this is the GL₂ representative from the quotient action, NOT the full Levi element.
def lift (α β γ δ : k) : Matrix V V k :=
  let G := Wedge2Formalization.N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
  Matrix.fromBlocks (1 : Matrix I I k) 0 0 (Matrix.fromBlocks G 0 0 G)

/-- For the DoublePureSingular orbit, the upper-right block vanishes when the pair is
fixed. The argument uses:
1. D fixes pair → det(D) ≠ 0 (since the common kernel of rep₁, rep₂ is trivial)
2. B * Ωᵢ * Dᵀ = 0 with Dᵀ invertible → B * Ωᵢ = 0
3. B in left-kernel of both Ω₁ and Ω₂ → B = 0 (images span the full space) -/
private theorem upperRight_zero_of_fixing_pair
    (B : Matrix I W k) (D : Matrix W W k)
    (hDfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector (k := k) D)
    (htop₁ : B * Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k) * Dᵀ = 0)
    (htop₂ : B * Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k) * Dᵀ = 0) :
    B = 0 := by
  -- D fixes the pair, so D ∈ K. From K-table shape, det(D) ≠ 0.
  have hDmem : D ∈ Wedge2Formalization.Paper.N6.Row11.K (k := k) :=
    (Wedge2Formalization.Paper.N6.Row11.pointwise_stabilizer (k := k) D).1 hDfix
  rcases (Wedge2Formalization.Paper.N6.Row11.mem_K_table_iff (k := k) D).1 hDmem with
    ⟨a, b, c, d, p, q, r, s, u, v, hdet, hshape⟩
  have hDdet : Matrix.det D ≠ 0 :=
    Wedge2Formalization.Paper.N6.Row11.det_ne_zero_of_pointwise (k := k) D hDfix
  -- Cancel Dᵀ from B * repᵢ * Dᵀ = 0
  have hDetTr : Matrix.det Dᵀ ≠ 0 := by rwa [Matrix.det_transpose]
  -- Cancel Dᵀ (invertible) from B * repᵢ * Dᵀ = 0
  have hDTu : IsUnit (Matrix.det Dᵀ) := IsUnit.mk0 _ hDetTr
  have hDTinv : Dᵀ * (Dᵀ)⁻¹ = (1 : Matrix W W k) := Matrix.mul_nonsing_inv _ hDTu
  have cancel : ∀ (X : Matrix I W k), X * Dᵀ = 0 → X = 0 := by
    intro X hX
    have : X = X * Dᵀ * (Dᵀ)⁻¹ := by
      conv_lhs => rw [show X = X * (1 : Matrix W W k) from (Matrix.mul_one X).symm]
      rw [← hDTinv, ← Matrix.mul_assoc]
    rw [this, hX, Matrix.zero_mul]
  have hBR1 : B * Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k) = 0 :=
    cancel _ htop₁
  have hBR2 : B * Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k) = 0 :=
    cancel _ htop₂
  -- B is 1×6. The common left-kernel of rep₁ = ω₁₃⊕ω₁₃ and rep₂ = ω₂₃⊕ω₂₃ is trivial.
  -- Unfold N6.Row11.rep to N6DoublePureSingular.rep for simp
  simp only [Wedge2Formalization.Paper.N6.Row11.rep₁, Wedge2Formalization.Paper.N6.Row11.rep₂] at hBR1 hBR2
  have e1 := congr_fun (congr_fun hBR1 0) (Sum.inl 0)
  have e2 := congr_fun (congr_fun hBR1 0) (Sum.inl 2)
  have e3 := congr_fun (congr_fun hBR1 0) (Sum.inr 0)
  have e4 := congr_fun (congr_fun hBR1 0) (Sum.inr 2)
  have e5 := congr_fun (congr_fun hBR2 0) (Sum.inl 2)
  have e6 := congr_fun (congr_fun hBR2 0) (Sum.inr 2)
  simp only [Matrix.zero_apply, Matrix.mul_apply, Fintype.sum_sum_type,
    Wedge2Formalization.N6DoublePureSingular.rep₁,
    Wedge2Formalization.N6DoublePureSingular.rep₂,
    Wedge2Formalization.N3PureSingular.ω13,
    Wedge2Formalization.N3PureSingular.ω23,
    fromBlocks, Fin.sum_univ_three,
    Matrix.of_apply, Matrix.cons_val', Matrix.cons_val_zero, Matrix.cons_val_one,
    Matrix.head_cons, Matrix.head_fin_const, mul_zero, mul_one, mul_neg,
    add_zero, zero_add, neg_eq_zero] at e1 e2 e3 e4 e5 e6
  ext i j; fin_cases i; rcases j with j | j <;> fin_cases j <;> simp_all

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
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k)) *
        (fromBlocks A0 B C D)ᵀ =
      fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k)) := by
      rw [← hgEq]; exact h1
    have h2_blocks : fromBlocks A0 B C D *
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k)) *
        (fromBlocks A0 B C D)ᵀ =
      fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k)) := by
      rw [← hgEq]; exact h2
    -- Extract (1,2) blocks
    have htop₁ : B * Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k) * Dᵀ = 0 := by
      have := congrArg toBlocks₁₂ h1_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    have htop₂ : B * Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k) * Dᵀ = 0 := by
      have := congrArg toBlocks₁₂ h2_blocks
      simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- Extract (2,2) blocks → D fixes N6 pair
    have hDfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector (k := k) D := by
      constructor
      · have := congrArg toBlocks₂₂ h1_blocks
        simpa [fromBlocks_multiply, fromBlocks_transpose] using this
      · have := congrArg toBlocks₂₂ h2_blocks
        simpa [fromBlocks_multiply, fromBlocks_transpose] using this
    -- B = 0 from the modified argument
    have hB0 : B = 0 := upperRight_zero_of_fixing_pair (k := k) B D hDfix htop₁ htop₂
    -- D ∈ K
    have hDmem : D ∈ Wedge2Formalization.Paper.N6.Row11.K (k := k) :=
      (Wedge2Formalization.Paper.N6.Row11.pointwise_stabilizer (k := k) D).1 hDfix
    exact ⟨A0, C, D, hDmem, by rw [hgEq]; simp [hB0]⟩
  · rintro ⟨A0', C', D', hD', rfl⟩
    have hDfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector (k := k) D' :=
      (Wedge2Formalization.Paper.N6.Row11.pointwise_stabilizer (k := k) D').2 hD'
    constructor
    · show fromBlocks A0' 0 C' D' *
          fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₁ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.1
    · show fromBlocks A0' 0 C' D' *
          fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k)) *
          (fromBlocks A0' 0 C' D')ᵀ =
        fromBlocks 0 0 0 (Wedge2Formalization.Paper.N6.Row11.rep₂ (k := k))
      simp [fromBlocks_multiply, fromBlocks_transpose]
      exact hDfix.2

theorem mem_K_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      ∃ A0 : Matrix I I k,
        ∃ C : Matrix W I k,
          ∃ D : Matrix W W k,
            D ∈ Wedge2Formalization.Paper.N6.Row11.K (k := k) ∧
              g = Matrix.fromBlocks A0 0 C D := Iff.rfl

/-- Public name for the explicit Appendix A row-3 kernel cell. -/
def TableCell : Set (Matrix V V k) :=
  { g |
      ∃ a0 a b c d p q r s u v : k,
        ∃ C : Matrix W I k,
          a * d - b * c ≠ 0 ∧
          g = Matrix.fromBlocks (scalarBlock (k := k) a0) 0 C
              (Wedge2Formalization.Paper.N6.Row11.Levi (k := k) a b c d *
                Wedge2Formalization.Paper.N6.Row11.U (k := k) p s q u r v) }

theorem mem_K_table_iff (g : Matrix V V k) :
    g ∈ K (k := k) ↔
      g ∈ TableCell (k := k) := by
  constructor
  · rintro ⟨A0, C, D, hD, rfl⟩
    rcases (Wedge2Formalization.Paper.N6.Row11.mem_K_table_iff (k := k) D).1 hD with
      ⟨a, b, c, d, p, q, r, s, u, v, hdet, hshape⟩
    exact ⟨A0 0 0, a, b, c, d, p, q, r, s, u, v, C, hdet, by
      rw [hshape]; ext i j
      rcases i with i | i <;> rcases j with j | j
      · fin_cases i; fin_cases j; simp [scalarBlock, fromBlocks]
      · simp [fromBlocks]
      · simp [fromBlocks]
      · simp [fromBlocks]⟩
  · rintro ⟨a0, a, b, c, d, p, q, r, s, u, v, C, hdet, rfl⟩
    exact ⟨scalarBlock (k := k) a0, C, _,
      (Wedge2Formalization.Paper.N6.Row11.mem_K_table_iff (k := k) _).2
        ⟨a, b, c, d, p, q, r, s, u, v, hdet, rfl⟩, rfl⟩

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

theorem quotient_action (α β γ δ : k) (hdet : α * δ - β * γ ≠ 0) :
    ActsOnOrderedPair
      (fun Ω g => g * Ω * gᵀ : Matrix V V k → Matrix V V k → Matrix V V k)
      (rep₁ (k := k)) (rep₂ (k := k))
      (lift (k := k) α β γ δ)
      (!![α, β; γ, δ]) := by
  have h6 := Wedge2Formalization.Paper.N6.Row11.quotient_action (k := k) α β γ δ
  simp only [ActsOnOrderedPair, Wedge2Formalization.N6DoublePureSingular.ActBivector] at h6
  refine ⟨?_, ?_⟩
  · show lift (k := k) α β γ δ * rep₁ (k := k) * (lift (k := k) α β γ δ)ᵀ =
        α • rep₁ (k := k) + β • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h6.1
  · show lift (k := k) α β γ δ * rep₂ (k := k) * (lift (k := k) α β γ δ)ᵀ =
        γ • rep₁ (k := k) + δ • rep₂ (k := k)
    simp only [rep₁, rep₂, lift]
    rw [embed_act, embed_lincomb]
    exact congrArg (fromBlocks 0 0 0) h6.2

/-- Every invertible setwise stabilizer induces a coefficient matrix in the projective GL₂.
The proof delegates to the W-block pencil determinant structure:
since rep = fromBlocks 0 0 0 Ω_W, the W-block action inherits the parent N6.Row11
pencil, and det(M) ≠ 0 follows from rep independence + g invertible.
For GL₂ quotient, M is invertible with no further constraints. -/
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
  simp only [Qproj, Wedge2Formalization.Paper.N6.Row11.Qproj, Set.mem_setOf_eq]
  -- Unfold ActsOnOrderedPair to get raw equalities
  simp only [ActsOnOrderedPair] at hact1 hact2
  intro hdet0
  have hM_dep : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
    simpa [Matrix.det_fin_two] using hdet0
  -- Rep independence from parent N6.Row11
  have rep_indep : ∀ a b : k, a • rep₁ (k := k) + b • rep₂ (k := k) = 0 → a = 0 ∧ b = 0 := by
    intro a b h
    have h22 := congrArg toBlocks₂₂ h
    simp only [rep₁, rep₂, fromBlocks_smul, fromBlocks_add, toBlocks₂₂, fromBlocks,
      Matrix.of_apply, Sum.elim_inr, Matrix.zero_apply, Matrix.smul_apply,
      Matrix.add_apply, smul_zero, add_zero] at h22
    exact Wedge2Formalization.Paper.N6.Row11.rep_pair_independent (k := k) h22
  -- If det M = 0, the combination (M 1 1)•rep₁ + (-(M 0 1))•rep₂ maps to 0
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
  -- Similarly the other combination
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
  -- All of M except M 0 0 is zero. hact2: g rep₂ gᵀ = 0 • rep₁ + 0 • rep₂ = 0
  have hrep2_zero : rep₂ (k := k) = 0 := by
    have : g * rep₂ (k := k) * gᵀ = 0 := by
      have := hact2; rw [h10, h11] at this; simpa using this
    exact (actBivector_eq_zero_iff_of_det_ne_zero (k := k) _ g hg).1 this
  exact absurd ((rep_indep 0 1 (by simpa using hrep2_zero)).2) one_ne_zero

end Row3

end N7
end Paper
end Wedge2Formalization
