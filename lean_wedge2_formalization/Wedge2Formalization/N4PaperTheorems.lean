import Wedge2Formalization.N4PaperSummary

namespace Wedge2Formalization
namespace N4PaperTheorems

open Matrix

variable {k : Type*} [Field k]

/-!
# `n = 4` Paper-Facing Theorems

This is the file a reader should open first for the finished geometric `n = 4`
stabilizer calculations.

Each row is presented as a tiny API:

- `K` for the exact pointwise kernel;
- `U` and `L` when the paper displays a unipotent / Levi decomposition;
- `Qproj` for the exact projective quotient family, and `Q` for an explicit section
  when the paper uses one;
- `pointwise_stabilizer` for the exact kernel statement;
- `quotient_action` or its split pieces for the action on the ordered basis.

The lower-level matrix calculations remain available in `N4PaperSummary`,
`N4Summary`, `N4`, and `N4PureSingular`.
-/

namespace Row1

/-- The repeated-support row `J_{a,2}`. -/
def K : Set (Matrix N4.V N4.V k) :=
  N4PaperSummary.row1_pointwiseKernel (k := k)

/-- The displayed unipotent radical `U_3`. -/
def U (x y z : k) : Matrix N4.V N4.V k :=
  N4PaperSummary.row1_U3 (k := k) x y z

/-- The displayed Levi `SL_2(k)` family. -/
def L (A : Matrix N4.I N4.I k) : Matrix N4.V N4.V k :=
  N4PaperSummary.row1_levi (k := k) A

/-- The exact coefficient-side projective Borel family in the ordered basis `(ω₁, ω₂)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | M 1 0 = 0 ∧ Matrix.det M ≠ 0 }

/-- A convenient explicit section inside the projective Borel family. -/
def Q : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | ∃ a b : k, a ≠ 0 ∧ M = N4PaperSummary.row1_borelCoeff (k := k) a b }

/-- The standard lift of the Borel family. -/
def lift (a b : k) : Matrix N4.V N4.V k :=
  N4PaperSummary.row1_borelLift (k := k) a b

/-- Exact paper-facing statement for the pointwise stabilizer. -/
theorem pointwise_stabilizer :
    PaperCore.ExactFamily
      (N4.FixesOnePointPairBivector (k := k))
      (K (k := k)) :=
  N4PaperSummary.row1_pointwise_stabilizer (k := k)

/-- The displayed unipotent radical lies in the pointwise stabilizer. -/
theorem U_pointwise
    (x y z : k) :
    N4.FixesOnePointPairBivector (U (k := k) x y z) :=
  N4PaperSummary.row1_U3_pointwise (k := k) x y z

/-- The displayed Levi family lies in the pointwise stabilizer. -/
theorem L_pointwise
    (A : Matrix N4.I N4.I k)
    (hA : A.det = 1) :
    N4.FixesOnePointPairBivector (L (k := k) A) :=
  N4PaperSummary.row1_levi_pointwise (k := k) A hA

/-- The coefficient matrix attached to the standard lift belongs to the Borel family. -/
theorem coeff_mem_Q
    (a b : k)
    (ha : a ≠ 0) :
    N4PaperSummary.row1_borelCoeff (k := k) a b ∈ Q (k := k) :=
  ⟨a, b, ha, rfl⟩

/-- The coefficient matrix attached to the standard lift lies in the exact projective Borel
family. -/
theorem coeff_mem_Qproj
    (a b : k)
    (ha : a ≠ 0) :
    N4PaperSummary.row1_borelCoeff (k := k) a b ∈ Qproj (k := k) := by
  constructor
  · simp [N4PaperSummary.row1_borelCoeff]
  · simp [N4PaperSummary.row1_borelCoeff, Matrix.det_fin_two, ha]

/-- The quotient action is the displayed Borel action on the ordered basis. -/
theorem quotient_action
    (a b : k)
    (ha : a ≠ 0) :
    PaperCore.ActsOnOrderedPair
      (N4.ActBivector (k := k))
      (N4.onePointRep₁ (k := k))
      (N4.onePointRep₂ (k := k))
      (lift (k := k) a b)
      (N4PaperSummary.row1_borelCoeff (k := k) a b) :=
  N4PaperSummary.row1_borel_lift (k := k) a b ha

/-- The repeated-support second basis vector is nonzero. -/
private theorem rep₂_ne_zero : N4.onePointRep₂ (k := k) ≠ 0 := by
  intro hzero
  have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) hzero
  simpa [N4.onePointRep₂, N4.ω12, N4.J] using h01

/-- Every invertible setwise stabilizer induces an upper-triangular coefficient action; this is
the exact statement that the quotient on `P(L)` is the Borel. -/
theorem quotient_image_projective
    (g : Matrix N4.V N4.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres : N4.PreservesOnePointSubspaceBivector (k := k) g) :
    ∃ M ∈ Qproj (k := k),
      PaperCore.ActsOnOrderedPair
        (N4.ActBivector (k := k))
        (N4.onePointRep₁ (k := k))
        (N4.onePointRep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have hrep₂_det :
      Matrix.det (N4.onePointRep₂ (k := k)) = 0 := by
    rw [N4.onePointRep₂, N4.ω12, Matrix.det_fromBlocks_zero₂₁]
    simpa using (Matrix.det_zero (n := N4.I) (R := k))
  have hrep₂_det_zero :
      Matrix.det (γ • (N4.onePointRep₁ (k := k)) + δ • (N4.onePointRep₂ (k := k))) = 0 := by
    rw [← h2, N4Summary.actBivector_det, hrep₂_det]
    simp
  have hγ : γ = 0 :=
    N4Summary.onePoint_action_case_of_rankDrop
      (k := k) (γ := γ) (δ := δ) hrep₂_det_zero
  have hrep₁_det :
      Matrix.det (N4.onePointRep₁ (k := k)) = 1 := by
    simpa using (N4Summary.onePoint_det_in_basis (k := k) (a := (1 : k)) (b := (0 : k)))
  have hrep₁_det_ne_zero :
      Matrix.det (α • (N4.onePointRep₁ (k := k)) + β • (N4.onePointRep₂ (k := k))) ≠ 0 := by
    rw [← h1, N4Summary.actBivector_det, hrep₁_det]
    simpa [mul_assoc] using (mul_ne_zero hg hg)
  have hα : α ≠ 0 := by
    intro hα0
    exact hrep₁_det_ne_zero ((N4Summary.onePoint_det_zero_iff (k := k) (a := α) (b := β)).2 hα0)
  have hδ : δ ≠ 0 := by
    intro hδ0
    have hzero : N4.ActBivector (N4.onePointRep₂ (k := k)) g = 0 := by
      simpa [hγ, hδ0] using h2
    have hrep₂_zero :=
      (N4Summary.actBivector_eq_zero_iff_of_det_ne_zero
        (k := k) (Ω := N4.onePointRep₂ (k := k)) (g := g) hg).1 hzero
    exact rep₂_ne_zero (k := k) hrep₂_zero
  refine ⟨!![α, β; 0, δ], ?_, ?_⟩
  · constructor
    · simp
    · simp [Matrix.det_fin_two, hα, hδ]
  · constructor
    · simpa [PaperCore.ActsOnOrderedPair] using h1
    · simpa [PaperCore.ActsOnOrderedPair, hγ] using h2

end Row1

namespace Row2

/-- The split-support row `J_{a,1} + J_{b,1}`. -/
def K : Set (Matrix N4.V N4.V k) :=
  N4PaperSummary.row2_pointwiseKernel (k := k)

/-- The displayed product Levi `SL_2(k) × SL_2(k)`. -/
def L (A D : Matrix N4.I N4.I k) : Matrix N4.V N4.V k :=
  Matrix.fromBlocks A 0 0 D

/-- The coefficient-side normalizer of the split torus in the ordered basis. -/
def Q : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M |
      (∃ u v : k, M = N4PaperSummary.row2_torusCoeff (k := k) u v) ∨
      (∃ u v : k, M = N4PaperSummary.row2_swapCoeff (k := k) u v) }

/-- The standard torus lift. -/
def torusLift (u v : k) : Matrix N4.V N4.V k :=
  N4PaperSummary.row2_torusLift (k := k) u v

/-- The standard swap-coset lift. -/
def swapLift (u v : k) : Matrix N4.V N4.V k :=
  N4PaperSummary.row2_swapLift (k := k) u v

/-- Exact paper-facing statement for the pointwise stabilizer. -/
theorem pointwise_stabilizer :
    PaperCore.ExactFamily
      (N4.FixesSplitPairBivector (k := k))
      (K (k := k)) :=
  N4PaperSummary.row2_pointwise_stabilizer (k := k)

/-- The standard torus coefficient matrix belongs to the displayed quotient family. -/
theorem torusCoeff_mem_Q
    (u v : k) :
    N4PaperSummary.row2_torusCoeff (k := k) u v ∈ Q (k := k) :=
  Or.inl ⟨u, v, rfl⟩

/-- The standard swap coefficient matrix belongs to the displayed quotient family. -/
theorem swapCoeff_mem_Q
    (u v : k) :
    N4PaperSummary.row2_swapCoeff (k := k) u v ∈ Q (k := k) :=
  Or.inr ⟨u, v, rfl⟩

/-- The torus part of the quotient action. -/
theorem torus_action
    (u v : k) :
    PaperCore.ActsOnOrderedPair
      (N4.ActBivector (k := k))
      (N4.splitRep₁ (k := k))
      (N4.splitRep₂ (k := k))
      (torusLift (k := k) u v)
      (N4PaperSummary.row2_torusCoeff (k := k) u v) :=
  N4PaperSummary.row2_torus_lift (k := k) u v

/-- The swap-coset part of the quotient action. -/
theorem swap_action
    (u v : k) :
    PaperCore.ActsOnOrderedPair
      (N4.ActBivector (k := k))
      (N4.splitRep₁ (k := k))
      (N4.splitRep₂ (k := k))
      (swapLift (k := k) u v)
      (N4PaperSummary.row2_swapCoeff (k := k) u v) :=
  N4PaperSummary.row2_swap_lift (k := k) u v

/-- The split-support basis vectors are linearly independent. -/
private theorem rep_pair_independent
    {a b : k}
    (h : a • (N4.splitRep₁ (k := k)) + b • (N4.splitRep₂ (k := k)) = 0) :
    a = 0 ∧ b = 0 := by
  have h01 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 1)) h
  have ha : a = 0 := by
    simpa [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, N4.J] using h01
  have h23 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) h
  have hab : a + b = 0 := by
    simpa [N4.splitRep₁, N4.splitRep₂, N4.ω12, N4.ω34, N4.J] using h23
  have hb : b = 0 := by
    simpa [ha] using hab
  exact ⟨ha, hb⟩

/-- A nontrivial linear combination of the split-support basis vectors is nonzero. -/
private theorem rep_pair_linearCombination_ne_zero
    {a b : k}
    (h : a ≠ 0 ∨ b ≠ 0) :
    a • (N4.splitRep₁ (k := k)) + b • (N4.splitRep₂ (k := k)) ≠ 0 := by
  intro hzero
  rcases rep_pair_independent (k := k) hzero with ⟨ha, hb⟩
  rcases h with ha' | hb'
  · exact ha' ha
  · exact hb' hb

/-- Every invertible setwise stabilizer induces a coefficient matrix in the split-torus
normalizer. -/
theorem quotient_image_projective
    (g : Matrix N4.V N4.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres : N4.PreservesSplitSubspaceBivector (k := k) g) :
    ∃ M ∈ Q (k := k),
      PaperCore.ActsOnOrderedPair
        (N4.ActBivector (k := k))
        (N4.splitRep₁ (k := k))
        (N4.splitRep₂ (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have hrep2_det :
      Matrix.det (N4.splitRep₂ (k := k)) = 0 := by
    simpa using
      (N4Summary.split_det_zero_iff (k := k) (a := (0 : k)) (b := (1 : k))).2 (Or.inl rfl)
  have hrep2_det_zero :
      Matrix.det (γ • (N4.splitRep₁ (k := k)) + δ • (N4.splitRep₂ (k := k))) = 0 := by
    rw [← h2, N4Summary.actBivector_det, hrep2_det]
    simp
  have hω12_det :
      Matrix.det ((N4.splitRep₁ (k := k)) - (N4.splitRep₂ (k := k))) = 0 := by
    have h :=
      (N4Summary.split_det_zero_iff (k := k) (a := (1 : k)) (b := (-1 : k))).2 (Or.inr (by ring))
    simpa [sub_eq_add_neg] using h
  have hω12_action :
      N4.ActBivector ((N4.splitRep₁ (k := k)) - (N4.splitRep₂ (k := k))) g =
        (α - γ) • (N4.splitRep₁ (k := k)) + (β - δ) • (N4.splitRep₂ (k := k)) := by
    calc
      N4.ActBivector ((N4.splitRep₁ (k := k)) - (N4.splitRep₂ (k := k))) g
          =
            N4.ActBivector
              ((N4.splitRep₁ (k := k)) + (-1 : k) • (N4.splitRep₂ (k := k))) g := by
                simp [sub_eq_add_neg]
      _ =
          N4.ActBivector (N4.splitRep₁ (k := k)) g +
            (-1 : k) • N4.ActBivector (N4.splitRep₂ (k := k)) g := by
              rw [N4Summary.actBivector_add, N4Summary.actBivector_smul]
      _ =
          (α • (N4.splitRep₁ (k := k)) + β • (N4.splitRep₂ (k := k))) +
            (-1 : k) • (γ • (N4.splitRep₁ (k := k)) + δ • (N4.splitRep₂ (k := k))) := by
              rw [h1, h2]
      _ = (α - γ) • (N4.splitRep₁ (k := k)) + (β - δ) • (N4.splitRep₂ (k := k)) := by
            ext i j
            simp [sub_eq_add_neg]
            ring
  have hω12_det_zero :
      Matrix.det
        ((α - γ) • (N4.splitRep₁ (k := k)) + (β - δ) • (N4.splitRep₂ (k := k))) = 0 := by
    rw [← hω12_action, N4Summary.actBivector_det, hω12_det]
    simp
  have hcoeff_det_ne : α * δ - β * γ ≠ 0 := by
    intro hdet
    by_cases hbd : β = 0 ∧ δ = 0
    · rcases hbd with ⟨hβ, hδ⟩
      have hzero2 :
          N4.ActBivector
            (γ • (N4.splitRep₁ (k := k)) - α • (N4.splitRep₂ (k := k))) g = 0 := by
        calc
          N4.ActBivector
              (γ • (N4.splitRep₁ (k := k)) - α • (N4.splitRep₂ (k := k))) g
              =
                γ • N4.ActBivector (N4.splitRep₁ (k := k)) g +
                  (-α) • N4.ActBivector (N4.splitRep₂ (k := k)) g := by
                    rw [sub_eq_add_neg, N4Summary.actBivector_add, N4Summary.actBivector_smul]
                    rw [show -(α • (N4.splitRep₂ (k := k))) = (-α) • (N4.splitRep₂ (k := k)) by simp]
                    rw [N4Summary.actBivector_smul]
          _ =
              γ • (α • (N4.splitRep₁ (k := k)) + β • (N4.splitRep₂ (k := k))) +
                (-α) • (γ • (N4.splitRep₁ (k := k)) + δ • (N4.splitRep₂ (k := k))) := by
                  rw [h1, h2]
          _ = 0 := by
                ext i j
                simp [hβ, hδ, sub_eq_add_neg]
                ring
      by_cases hca : γ = 0 ∧ α = 0
      · rcases hca with ⟨hγ, hα⟩
        have hrep2_zero : N4.ActBivector (N4.splitRep₂ (k := k)) g = 0 := by
          simpa [hγ, hδ] using h2
        have hrep2_orig_zero :=
          (N4Summary.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := N4.splitRep₂ (k := k)) (g := g) hg).1 hrep2_zero
        have hrep2_ne :
            N4.splitRep₂ (k := k) ≠ 0 := by
          intro hzero
          have h01 := congrArg (fun M => M (Sum.inr 0) (Sum.inr 1)) hzero
          simpa [N4.splitRep₂, N4.ω34, N4.J] using h01
        exact hrep2_ne hrep2_orig_zero
      · have hnonzero : γ ≠ 0 ∨ α ≠ 0 := by
          by_cases hγ : γ = 0
          · right
            intro hα
            exact hca ⟨hγ, hα⟩
          · exact Or.inl hγ
        have hcomb_ne :
            γ • (N4.splitRep₁ (k := k)) - α • (N4.splitRep₂ (k := k)) ≠ 0 := by
          have hcomb_ne' :
              γ • (N4.splitRep₁ (k := k)) + (-α) • (N4.splitRep₂ (k := k)) ≠ 0 := by
            apply rep_pair_linearCombination_ne_zero (k := k)
            rcases hnonzero with hγ | hα
            · exact Or.inl hγ
            · exact Or.inr (neg_ne_zero.mpr hα)
          simpa [sub_eq_add_neg] using hcomb_ne'
        have horig_zero :=
          (N4Summary.actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := γ • (N4.splitRep₁ (k := k)) - α • (N4.splitRep₂ (k := k)))
            (g := g)
            hg).1 hzero2
        exact hcomb_ne horig_zero
    · have hnonzero : β ≠ 0 ∨ δ ≠ 0 := by
        by_cases hβ : β = 0
        · right
          intro hδ
          exact hbd ⟨hβ, hδ⟩
        · exact Or.inl hβ
      have hcomb_ne :
          δ • (N4.splitRep₁ (k := k)) - β • (N4.splitRep₂ (k := k)) ≠ 0 := by
        have hcomb_ne' :
            δ • (N4.splitRep₁ (k := k)) + (-β) • (N4.splitRep₂ (k := k)) ≠ 0 := by
          apply rep_pair_linearCombination_ne_zero (k := k)
          rcases hnonzero with hβ | hδ
          · exact Or.inr (neg_ne_zero.mpr hβ)
          · exact Or.inl hδ
        simpa [sub_eq_add_neg] using hcomb_ne'
      have hzero1 :
          N4.ActBivector
            (δ • (N4.splitRep₁ (k := k)) - β • (N4.splitRep₂ (k := k))) g = 0 := by
        calc
          N4.ActBivector
              (δ • (N4.splitRep₁ (k := k)) - β • (N4.splitRep₂ (k := k))) g
              =
                δ • N4.ActBivector (N4.splitRep₁ (k := k)) g +
                  (-β) • N4.ActBivector (N4.splitRep₂ (k := k)) g := by
                    rw [sub_eq_add_neg, N4Summary.actBivector_add, N4Summary.actBivector_smul]
                    rw [show -(β • (N4.splitRep₂ (k := k))) = (-β) • (N4.splitRep₂ (k := k)) by simp]
                    rw [N4Summary.actBivector_smul]
          _ =
              δ • (α • (N4.splitRep₁ (k := k)) + β • (N4.splitRep₂ (k := k))) +
                (-β) • (γ • (N4.splitRep₁ (k := k)) + δ • (N4.splitRep₂ (k := k))) := by
                  rw [h1, h2]
          _ = (α * δ - β * γ) • (N4.splitRep₁ (k := k)) := by
                ext i j
                simp [sub_eq_add_neg]
                ring
          _ = 0 := by
                simp [hdet]
      have horig_zero :=
        (N4Summary.actBivector_eq_zero_iff_of_det_ne_zero
          (k := k)
          (Ω := δ • (N4.splitRep₁ (k := k)) - β • (N4.splitRep₂ (k := k)))
          (g := g)
          hg).1 hzero1
      exact hcomb_ne horig_zero
  have hM_mem : ((!![α, β; γ, δ] : Matrix (Fin 2) (Fin 2) k)) ∈ Q (k := k) := by
    rcases
      N4Summary.split_action_cases_of_rankDrop
        (k := k)
        (α := α)
        (β := β)
        (γ := γ)
        (δ := δ)
        hcoeff_det_ne
        hrep2_det_zero
        hω12_det_zero with
      hdiag | hswap
    · rcases hdiag with ⟨hγ, hsum⟩
      left
      have hbeta : β = δ - α := by
        calc
          β = α + β - α := by ring
          _ = δ - α := by rw [hsum]
      refine ⟨α, δ, ?_⟩
      ext i j
      fin_cases i <;> fin_cases j
      · simp [N4PaperSummary.row2_torusCoeff, hbeta]
      · simp [N4PaperSummary.row2_torusCoeff, hbeta]
      · simp [N4PaperSummary.row2_torusCoeff, hγ]
      · simp [N4PaperSummary.row2_torusCoeff]
    · rcases hswap with ⟨hαγ, hγδ⟩
      right
      refine ⟨α, α + β, ?_⟩
      have hδ : δ = -α := by
        calc
          δ = γ + δ - γ := by ring
          _ = 0 - γ := by rw [hγδ]
          _ = -γ := by ring
          _ = -α := by rw [hαγ]
      ext i j
      fin_cases i <;> fin_cases j
      · simp [N4PaperSummary.row2_swapCoeff]
      · simp [N4PaperSummary.row2_swapCoeff]
      · simp [N4PaperSummary.row2_swapCoeff, hαγ]
      · simp [N4PaperSummary.row2_swapCoeff, hδ]
  refine ⟨!![α, β; γ, δ], hM_mem, ?_⟩
  constructor
  · simpa [PaperCore.ActsOnOrderedPair] using h1
  · simpa [PaperCore.ActsOnOrderedPair] using h2

end Row2

namespace Row3

/-- The pure singular row `S_1 + S_2`. -/
def K : Set (Matrix N4PureSingular.I N4PureSingular.I k) :=
  N4PaperSummary.row3_pointwiseKernel (k := k)

/-- The displayed unipotent radical `U_2`, viewed in the local `3 × 3` block model. -/
def U (x y : k) : Matrix N3PureSingular.I N3PureSingular.I k :=
  N3PureSingular.pointwiseUnipotent (k := k) x y

/-- The displayed torus `G_m`, viewed in the local `3 × 3` block model. -/
def L (a : k) : Matrix N3PureSingular.I N3PureSingular.I k :=
  N3PureSingular.pointwiseScale (k := k) a

/-- The coefficient-side group `GL_2(k)` in the ordered basis `(ω_{12}, ω_{13})`. -/
def Q : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

/-- The standard `GL_2(k)` lift. -/
def lift (α β γ δ : k) : Matrix N4PureSingular.I N4PureSingular.I k :=
  N4PaperSummary.row3_GL2Lift (k := k) α β γ δ

/-- Exact paper-facing statement for the pointwise stabilizer. -/
theorem pointwise_stabilizer :
    PaperCore.ExactFamily
      (N4PureSingular.FixesPureSingularPair (k := k))
      (K (k := k)) :=
  N4PaperSummary.row3_pointwise_stabilizer (k := k)

/-- The displayed local unipotent family fixes the pure singular pair pointwise. -/
theorem U_pointwise
    (x y : k) :
    N3PureSingular.FixesPairBivector (U (k := k) x y) :=
  N3PureSingular.pointwise_unipotent_family (k := k) x y

/-- The displayed local torus fixes the pure singular pair pointwise. -/
theorem L_pointwise
    (a : k)
    (ha : a ≠ 0) :
    N3PureSingular.FixesPairBivector (L (k := k) a) :=
  N3PureSingular.pointwise_scale_family (k := k) a ha

/-- The coefficient matrix attached to the standard lift lies in `GL_2(k)`. -/
theorem coeff_mem_Q
    (α β γ δ : k)
    (hΔ : α * δ - β * γ ≠ 0) :
    N4PaperSummary.row3_GL2Coeff (k := k) α β γ δ ∈ Q (k := k) := by
  simpa [Q, N4PaperSummary.row3_GL2Coeff, Matrix.det_fin_two] using hΔ

/-- The quotient action is the displayed `GL_2(k)` action on the ordered basis. -/
theorem quotient_action
    (α β γ δ : k)
    (hΔ : α * δ - β * γ ≠ 0) :
    PaperCore.ActsOnOrderedPair
      (N4PureSingular.ActBivector (k := k))
      (N4PureSingular.ω12 (k := k))
      (N4PureSingular.ω13 (k := k))
      (lift (k := k) α β γ δ)
      (N4PaperSummary.row3_GL2Coeff (k := k) α β γ δ) := by
  exact N4PaperSummary.row3_GL2_lift (k := k) α β γ δ

/-- The pure singular basis vectors are linearly independent. -/
private theorem rep_pair_independent
    {a b : k}
    (h : a • (N4PureSingular.ω12 (k := k)) + b • (N4PureSingular.ω13 (k := k)) = 0) :
    a = 0 ∧ b = 0 := by
  have h01 := congrArg (fun M => M 0 1) h
  have ha : a = 0 := by
    simpa [N4PureSingular.ω12, N4PureSingular.ω13] using h01
  have h02 := congrArg (fun M => M 0 2) h
  have hb : b = 0 := by
    simpa [N4PureSingular.ω12, N4PureSingular.ω13] using h02
  exact ⟨ha, hb⟩

/-- A nontrivial linear combination of the pure singular basis vectors is nonzero. -/
private theorem rep_pair_linearCombination_ne_zero
    {a b : k}
    (h : a ≠ 0 ∨ b ≠ 0) :
    a • (N4PureSingular.ω12 (k := k)) + b • (N4PureSingular.ω13 (k := k)) ≠ 0 := by
  intro hzero
  rcases rep_pair_independent (k := k) hzero with ⟨ha, hb⟩
  rcases h with ha' | hb'
  · exact ha' ha
  · exact hb' hb

/-- Every invertible setwise stabilizer induces an invertible coefficient matrix on the
ordered pure singular basis. -/
theorem quotient_image_projective
    (g : Matrix N4PureSingular.I N4PureSingular.I k)
    (hg : Matrix.det g ≠ 0)
    (hpres : N4PureSingular.PreservesPureSingularSubspaceBivector (k := k) g) :
    ∃ M ∈ Q (k := k),
      PaperCore.ActsOnOrderedPair
        (N4PureSingular.ActBivector (k := k))
        (N4PureSingular.ω12 (k := k))
        (N4PureSingular.ω13 (k := k))
        g M := by
  rcases hpres with ⟨α, β, γ, δ, h1, h2⟩
  have hcoeff_det_ne : α * δ - β * γ ≠ 0 := by
    intro hdet
    by_cases hbd : β = 0 ∧ δ = 0
    · rcases hbd with ⟨hβ, hδ⟩
      have hzero2 :
          N4PureSingular.ActBivector
            (γ • (N4PureSingular.ω12 (k := k)) - α • (N4PureSingular.ω13 (k := k))) g = 0 := by
        calc
          N4PureSingular.ActBivector
              (γ • (N4PureSingular.ω12 (k := k)) - α • (N4PureSingular.ω13 (k := k))) g
              =
                γ • N4PureSingular.ActBivector (N4PureSingular.ω12 (k := k)) g +
                  (-α) • N4PureSingular.ActBivector (N4PureSingular.ω13 (k := k)) g := by
                    rw [sub_eq_add_neg, N4Summary.pureSingular_actBivector_add,
                      N4Summary.pureSingular_actBivector_smul]
                    rw [show -(α • (N4PureSingular.ω13 (k := k))) = (-α) • (N4PureSingular.ω13 (k := k)) by simp]
                    rw [N4Summary.pureSingular_actBivector_smul]
          _ =
              γ • (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) +
                (-α) • (γ • (N4PureSingular.ω12 (k := k)) + δ • (N4PureSingular.ω13 (k := k))) := by
                  rw [h1, h2]
          _ = 0 := by
                ext i j
                simp [hβ, hδ, sub_eq_add_neg]
                ring
      by_cases hca : γ = 0 ∧ α = 0
      · rcases hca with ⟨hγ, hα⟩
        have hω13_zero : N4PureSingular.ActBivector (N4PureSingular.ω13 (k := k)) g = 0 := by
          simpa [hγ, hδ] using h2
        have hω13_orig_zero :=
          (N4Summary.pureSingular_actBivector_eq_zero_iff_of_det_ne_zero
            (k := k) (Ω := N4PureSingular.ω13 (k := k)) (g := g) hg).1 hω13_zero
        have hω13_ne : N4PureSingular.ω13 (k := k) ≠ 0 := by
          intro hzero
          have h02 := congrArg (fun M => M 0 2) hzero
          simpa [N4PureSingular.ω13] using h02
        exact hω13_ne hω13_orig_zero
      · have hnonzero : γ ≠ 0 ∨ α ≠ 0 := by
          by_cases hγ : γ = 0
          · right
            intro hα
            exact hca ⟨hγ, hα⟩
          · exact Or.inl hγ
        have hcomb_ne :
            γ • (N4PureSingular.ω12 (k := k)) - α • (N4PureSingular.ω13 (k := k)) ≠ 0 := by
          have hcomb_ne' :
              γ • (N4PureSingular.ω12 (k := k)) + (-α) • (N4PureSingular.ω13 (k := k)) ≠ 0 := by
            apply rep_pair_linearCombination_ne_zero (k := k)
            rcases hnonzero with hγ | hα
            · exact Or.inl hγ
            · exact Or.inr (neg_ne_zero.mpr hα)
          simpa [sub_eq_add_neg] using hcomb_ne'
        have horig_zero :=
          (N4Summary.pureSingular_actBivector_eq_zero_iff_of_det_ne_zero
            (k := k)
            (Ω := γ • (N4PureSingular.ω12 (k := k)) - α • (N4PureSingular.ω13 (k := k)))
            (g := g)
            hg).1 hzero2
        exact hcomb_ne horig_zero
    · have hnonzero : β ≠ 0 ∨ δ ≠ 0 := by
        by_cases hβ : β = 0
        · right
          intro hδ
          exact hbd ⟨hβ, hδ⟩
        · exact Or.inl hβ
      have hcomb_ne :
          δ • (N4PureSingular.ω12 (k := k)) - β • (N4PureSingular.ω13 (k := k)) ≠ 0 := by
        have hcomb_ne' :
            δ • (N4PureSingular.ω12 (k := k)) + (-β) • (N4PureSingular.ω13 (k := k)) ≠ 0 := by
          apply rep_pair_linearCombination_ne_zero (k := k)
          rcases hnonzero with hβ | hδ
          · exact Or.inr (neg_ne_zero.mpr hβ)
          · exact Or.inl hδ
        simpa [sub_eq_add_neg] using hcomb_ne'
      have hzero1 :
          N4PureSingular.ActBivector
            (δ • (N4PureSingular.ω12 (k := k)) - β • (N4PureSingular.ω13 (k := k))) g = 0 := by
        calc
          N4PureSingular.ActBivector
              (δ • (N4PureSingular.ω12 (k := k)) - β • (N4PureSingular.ω13 (k := k))) g
              =
                δ • N4PureSingular.ActBivector (N4PureSingular.ω12 (k := k)) g +
                  (-β) • N4PureSingular.ActBivector (N4PureSingular.ω13 (k := k)) g := by
                    rw [sub_eq_add_neg, N4Summary.pureSingular_actBivector_add,
                      N4Summary.pureSingular_actBivector_smul]
                    rw [show -(β • (N4PureSingular.ω13 (k := k))) = (-β) • (N4PureSingular.ω13 (k := k)) by simp]
                    rw [N4Summary.pureSingular_actBivector_smul]
          _ =
              δ • (α • (N4PureSingular.ω12 (k := k)) + β • (N4PureSingular.ω13 (k := k))) +
                (-β) • (γ • (N4PureSingular.ω12 (k := k)) + δ • (N4PureSingular.ω13 (k := k))) := by
                  rw [h1, h2]
          _ = (α * δ - β * γ) • (N4PureSingular.ω12 (k := k)) := by
                ext i j
                simp [sub_eq_add_neg]
                ring
          _ = 0 := by
                simp [hdet]
      have horig_zero :=
        (N4Summary.pureSingular_actBivector_eq_zero_iff_of_det_ne_zero
          (k := k)
          (Ω := δ • (N4PureSingular.ω12 (k := k)) - β • (N4PureSingular.ω13 (k := k)))
          (g := g)
          hg).1 hzero1
      exact hcomb_ne horig_zero
  refine ⟨!![α, β; γ, δ], ?_, ?_⟩
  · simpa [Q, Matrix.det_fin_two] using hcoeff_det_ne
  · constructor
    · simpa [PaperCore.ActsOnOrderedPair] using h1
    · simpa [PaperCore.ActsOnOrderedPair] using h2

end Row3

end N4PaperTheorems
end Wedge2Formalization
