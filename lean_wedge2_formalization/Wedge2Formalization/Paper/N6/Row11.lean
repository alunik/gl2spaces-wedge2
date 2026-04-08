import Wedge2Formalization.Paper.Core
import Wedge2Formalization.N6Summary
import Wedge2Formalization.N6DoublePureSingular

namespace Wedge2Formalization
namespace Paper
namespace N6

open Matrix

variable {k : Type*} [Field k]

/-! Appendix A, `n = 6`, row 11.
Representative `⟨e₁∧e₃ + e₄∧e₆, e₂∧e₃ + e₅∧e₆⟩`.
Divisor `0`.
Claimed stabilizer:
`K_L = U_6 ⋊ GL₂(k)`, exact quotient family `Q_L = GL_2(k)`.
-/
namespace Row11

def rep₁ :
    Matrix Wedge2Formalization.N6DoublePureSingular.V
      Wedge2Formalization.N6DoublePureSingular.V k :=
  Wedge2Formalization.N6DoublePureSingular.rep₁ (k := k)

def rep₂ :
    Matrix Wedge2Formalization.N6DoublePureSingular.V
      Wedge2Formalization.N6DoublePureSingular.V k :=
  Wedge2Formalization.N6DoublePureSingular.rep₂ (k := k)

/-- Literal first basis vector from Appendix A, row 11. -/
def paperRep₁ :
    Matrix Wedge2Formalization.N6DoublePureSingular.V
      Wedge2Formalization.N6DoublePureSingular.V k :=
  rep₁ (k := k)

/-- Literal second basis vector from Appendix A, row 11. -/
def paperRep₂ :
    Matrix Wedge2Formalization.N6DoublePureSingular.V
      Wedge2Formalization.N6DoublePureSingular.V k :=
  rep₂ (k := k)

/-- The row-11 paper representative already agrees with the internal working pair. -/
def paperChange :
    Matrix Wedge2Formalization.N6DoublePureSingular.V
      Wedge2Formalization.N6DoublePureSingular.V k :=
  1

/-- Transport of the paper `\omega_1` to the internal working representative. -/
theorem paperRep₁_transport :
    Wedge2Formalization.N6DoublePureSingular.ActBivector
      (paperRep₁ (k := k)) (paperChange (k := k)) =
      rep₁ (k := k) := by
  simp [paperRep₁, paperChange, rep₁, Wedge2Formalization.N6DoublePureSingular.ActBivector]

/-- Transport of the paper `\omega_2` to the internal working representative. -/
theorem paperRep₂_transport :
    Wedge2Formalization.N6DoublePureSingular.ActBivector
      (paperRep₂ (k := k)) (paperChange (k := k)) =
      rep₂ (k := k) := by
  simp [paperRep₂, paperChange, rep₂, Wedge2Formalization.N6DoublePureSingular.ActBivector]

/-- The displayed six-parameter unipotent family `U_6`. -/
def U
    (x y z w r s : k) :
    Matrix Wedge2Formalization.N6DoublePureSingular.V
      Wedge2Formalization.N6DoublePureSingular.V k :=
  Wedge2Formalization.N6DoublePureSingular.coupledUnipotent
    (k := k) x y z w r s

/-- The displayed pointwise `GL₂(k)` family from the local `S_2^2` kernel model. -/
def Levi
    (a b c d : k) :
    Matrix Wedge2Formalization.N6DoublePureSingular.V
      Wedge2Formalization.N6DoublePureSingular.V k :=
  Wedge2Formalization.N6DoublePureSingular.magmaPointwiseGL2
    (k := k) a b c d

/-- The displayed kernel cell `U_6 ⋊ GL₂(k)` in the parameterization used by the
internal Magma-inspired local model. -/
def K :
    Set
      (Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) :=
  { g |
      ∃ a b c d p q r s u v : k,
        a * d - b * c ≠ 0 ∧
        g =
          Wedge2Formalization.N6DoublePureSingular.magmaKernelCell
            (k := k) a b c d p q r s u v }

/-- Exact coefficient-side quotient family `GL₂(k)`. -/
def Qproj : Set (Matrix (Fin 2) (Fin 2) k) :=
  { M | Matrix.det M ≠ 0 }

/-- Raw membership statement for the displayed kernel cell. -/
theorem mem_K_shape_iff
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) :
    g ∈ K (k := k) ↔
      ∃ a b c d p q r s u v : k,
        a * d - b * c ≠ 0 ∧
        g =
          Wedge2Formalization.N6DoublePureSingular.magmaKernelCell
            (k := k) a b c d p q r s u v := by
  rfl

/-- Paper-facing factorization of the displayed local kernel cell as `U_6` times the
pointwise `GL₂(k)` family. -/
theorem mem_K_iff
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) :
    g ∈ K (k := k) ↔
      ∃ a b c d p q r s u v : k,
        a * d - b * c ≠ 0 ∧
        g = Levi (k := k) a b c d * U (k := k) p s q u r v := by
  constructor
  · intro hg
    rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
      ⟨a, b, c, d, p, q, r, s, u, v, hΔ, hg⟩
    refine ⟨a, b, c, d, p, q, r, s, u, v, hΔ, ?_⟩
    simpa [Levi, U, Wedge2Formalization.N6DoublePureSingular.magmaKernelCell] using hg
  · rintro ⟨a, b, c, d, p, q, r, s, u, v, hΔ, rfl⟩
    exact
      (mem_K_shape_iff
        (k := k)
        (g := Levi (k := k) a b c d * U (k := k) p s q u r v)).2
        ⟨a, b, c, d, p, q, r, s, u, v, hΔ, by
          simp [Levi, U, Wedge2Formalization.N6DoublePureSingular.magmaKernelCell]⟩

/-- Table-facing kernel statement for Appendix A, row 11:
`K_L = U_6 \rtimes GL_2(k)`. -/
theorem mem_K_table_iff
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) :
    g ∈ K (k := k) ↔
      ∃ a b c d p q r s u v : k,
        a * d - b * c ≠ 0 ∧
        g = Levi (k := k) a b c d * U (k := k) p s q u r v :=
  mem_K_iff (k := k) (g := g)

/-- Exact pointwise criterion on the intrinsic block-diagonal-times-coupled cell for
the doubled pure singular row. -/
theorem cell_pointwise_iff_shape
    (x y z w r s : k)
    (G : Matrix Wedge2Formalization.N6DoublePureSingular.W
      Wedge2Formalization.N6DoublePureSingular.W k)
    (H : Matrix Wedge2Formalization.N6DoublePureSingular.I
      Wedge2Formalization.N6DoublePureSingular.I k) :
    Wedge2Formalization.N6DoublePureSingular.FixesPairBivector
      ((Matrix.fromBlocks G 0 0 H) *
        Wedge2Formalization.N6DoublePureSingular.coupledUnipotent
          (k := k) x y z w r s) ↔
      G 0 0 ≠ 0 ∧
        H 0 0 ≠ 0 ∧
        G =
          Wedge2Formalization.N3PureSingular.pureSingularShape
            (k := k) (G 0 0) (G 2 0) (G 2 1) ∧
        H =
          Wedge2Formalization.N3PureSingular.pureSingularShape
            (k := k) (H 0 0) (H 2 0) (H 2 1) :=
  Wedge2Formalization.N6Summary.doublePureSingular_blockDiagonal_coupled_right_product_iff_shape
    (k := k) (x := x) (y := y) (z := z) (w := w) (r := r) (s := s) (G := G) (H := H)

/-- The displayed unipotent family fixes the pair pointwise. -/
theorem U_pointwise
    (x y z w r s : k) :
    Wedge2Formalization.N6DoublePureSingular.FixesPairBivector (U (k := k) x y z w r s) :=
  Wedge2Formalization.N6Summary.doublePureSingular_coupled_unipotent_family
    (k := k) x y z w r s

/-- The displayed pointwise `GL₂(k)` family fixes the pair pointwise. -/
theorem Levi_pointwise
    (a b c d : k)
    (hΔ : a * d - b * c ≠ 0) :
    Wedge2Formalization.N6DoublePureSingular.FixesPairBivector
      (Levi (k := k) a b c d) :=
  Wedge2Formalization.N6Summary.doublePureSingular_magmaPointwiseGL2_family
    (k := k) a b c d hΔ

/-- The full displayed kernel cell fixes the pair pointwise. -/
theorem kernelCell_pointwise
    (a b c d p q r s u v : k)
    (hΔ : a * d - b * c ≠ 0) :
    Wedge2Formalization.N6DoublePureSingular.FixesPairBivector
      (Wedge2Formalization.N6DoublePureSingular.magmaKernelCell
        (k := k) a b c d p q r s u v) :=
  Wedge2Formalization.N6Summary.doublePureSingular_magmaKernelCell_pointwise
    (k := k) a b c d p q r s u v hΔ

theorem pointwise_of_mem_K
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hg : g ∈ K (k := k)) :
    Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g := by
  rcases (mem_K_shape_iff (k := k) (g := g)).1 hg with
    ⟨a, b, c, d, p, q, r, s, u, v, hΔ, rfl⟩
  exact kernelCell_pointwise (k := k) a b c d p q r s u v hΔ

/-- The simultaneous `GL₂(k)` lift has the expected quotient action on the ordered
basis `(rep₁, rep₂)`. -/
theorem quotient_action
    (α β γ δ : k) :
    let G : Matrix Wedge2Formalization.N6DoublePureSingular.W
        Wedge2Formalization.N6DoublePureSingular.W k :=
      Wedge2Formalization.N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k :=
      Matrix.fromBlocks G 0 0 G
    ActsOnOrderedPair
      (Wedge2Formalization.N6DoublePureSingular.ActBivector (k := k))
      (rep₁ (k := k))
      (rep₂ (k := k))
      g
      (!![α, β; γ, δ]) := by
  intro G g
  simpa [ActsOnOrderedPair, rep₁, rep₂] using
    Wedge2Formalization.N6Summary.doublePureSingular_GL2_lift_action
      (k := k) (α := α) (β := β) (γ := γ) (δ := δ)

theorem rep_pair_independent
    {a b : k}
    (h :
      a • rep₁ (k := k) + b • rep₂ (k := k) = 0) :
    a = 0 ∧ b = 0 := by
  have h02 := congrArg (fun M => M (Sum.inl 0) (Sum.inl 2)) h
  have ha : a = 0 := by
    simpa [rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h02
  have h12 := congrArg (fun M => M (Sum.inl 1) (Sum.inl 2)) h
  have hb : b = 0 := by
    simpa [ha, rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.fromBlocks, Matrix.add_apply,
      Matrix.smul_apply] using h12
  exact ⟨ha, hb⟩

private theorem rep₁_ne_zero : rep₁ (k := k) ≠ 0 := by
  intro hzero
  have hcomb : (1 : k) • rep₁ (k := k) + (0 : k) • rep₂ (k := k) = 0 := by
    simpa using hzero
  rcases rep_pair_independent (k := k) hcomb with ⟨h1, _⟩
  exact one_ne_zero h1

theorem quotient_image
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hg : Matrix.det g ≠ 0)
    (hpres :
      PreservesSpanPair
        (Wedge2Formalization.N6DoublePureSingular.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g) :
    ∃ M ∈ Qproj (k := k),
      ActsOnOrderedPair
        (Wedge2Formalization.N6DoublePureSingular.ActBivector (k := k))
        (rep₁ (k := k))
        (rep₂ (k := k))
        g M := by
  rcases hpres with ⟨M, hM⟩
  refine ⟨M, ?_, hM⟩
  change Matrix.det M ≠ 0
  by_contra hdet
  by_cases hcol0 : M 0 0 = 0 ∧ M 1 0 = 0
  · by_cases hcol1 : M 0 1 = 0 ∧ M 1 1 = 0
    · have hzero :
          Wedge2Formalization.N6DoublePureSingular.ActBivector
              (rep₁ (k := k)) g = 0 := by
        rcases hM with ⟨h1, _⟩
        simpa [ActsOnOrderedPair, hcol0.1, hcol1.1] using h1
      have horig :=
        (actBivector_eq_zero_iff_of_det_ne_zero
          (k := k) (Ω := rep₁ (k := k)) (g := g) hg).1 hzero
      exact rep₁_ne_zero (k := k) horig
    · have hzero :
          Wedge2Formalization.N6DoublePureSingular.ActBivector
              ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g = 0 := by
        rcases hM with ⟨h1, h2⟩
        calc
          Wedge2Formalization.N6DoublePureSingular.ActBivector
              ((-(M 1 1)) • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) g
              =
            (-(M 1 1)) •
                Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g +
              M 0 1 •
                Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N6DoublePureSingular.ActBivector, Matrix.mul_add,
                    Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
          _ =
            (-(M 1 1)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 1 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [h1, h2]
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
        Wedge2Formalization.N6DoublePureSingular.ActBivector
            ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g = 0 := by
        rcases hM with ⟨h1, h2⟩
        have hdet' : M 0 0 * M 1 1 = M 0 1 * M 1 0 := by
          have : M 0 0 * M 1 1 - M 0 1 * M 1 0 = 0 := by
            simpa [Matrix.det_fin_two] using hdet
          exact sub_eq_zero.mp this
        calc
          Wedge2Formalization.N6DoublePureSingular.ActBivector
              ((-(M 1 0)) • rep₁ (k := k) + M 0 0 • rep₂ (k := k)) g
              =
            (-(M 1 0)) •
                Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g +
              M 0 0 •
                Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g := by
                  simp [Wedge2Formalization.N6DoublePureSingular.ActBivector, Matrix.mul_add,
                    Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
          _ =
            (-(M 1 0)) • (M 0 0 • rep₁ (k := k) + M 0 1 • rep₂ (k := k)) +
              M 0 0 • (M 1 0 • rep₁ (k := k) + M 1 1 • rep₂ (k := k)) := by
                rw [h1, h2]
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

private theorem commonKernel_zero
    (x : Wedge2Formalization.N6DoublePureSingular.V → k)
    (h1 : rep₁ (k := k) *ᵥ x = 0)
    (h2 : rep₂ (k := k) *ᵥ x = 0) :
    x = 0 := by
  funext i
  cases i with
  | inl i =>
      fin_cases i
      ·
        have hcoord := congrArg (fun v => v (Sum.inl 2)) h1
        simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
          Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
          dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using hcoord
      ·
        have hcoord := congrArg (fun v => v (Sum.inl 2)) h2
        simpa [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
          Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
          dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using hcoord
      ·
        have hcoord := congrArg (fun v => v (Sum.inl 0)) h1
        simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
          Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
          dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using hcoord
  | inr i =>
      fin_cases i
      ·
        have hcoord := congrArg (fun v => v (Sum.inr 2)) h1
        simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
          Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
          dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using hcoord
      ·
        have hcoord := congrArg (fun v => v (Sum.inr 2)) h2
        simpa [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
          Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
          dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using hcoord
      ·
        have hcoord := congrArg (fun v => v (Sum.inr 0)) h1
        simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
          Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
          dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using hcoord

theorem det_ne_zero_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix :
      Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g) :
    Matrix.det g ≠ 0 := by
  intro hg0
  obtain ⟨v, hvne, hv0⟩ :=
    (Matrix.exists_mulVec_eq_zero_iff (M := gᵀ)).2
      (by simpa [Matrix.det_transpose] using hg0)
  have hv1 :
      rep₁ (k := k) *ᵥ v = 0 := by
    have hzero :
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g *ᵥ v = 0 := by
      calc
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g *ᵥ v
            =
          g *ᵥ (rep₁ (k := k) *ᵥ (gᵀ *ᵥ v)) := by
            simp [Wedge2Formalization.N6DoublePureSingular.ActBivector, Matrix.mulVec_mulVec,
              Matrix.mul_assoc]
        _ = 0 := by simp [hv0]
    have hfix1' :
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g =
          rep₁ (k := k) := by
      simpa [rep₁] using hfix.1
    rw [hfix1'] at hzero
    exact hzero
  have hv2 :
      rep₂ (k := k) *ᵥ v = 0 := by
    have hzero :
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g *ᵥ v = 0 := by
      calc
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g *ᵥ v
            =
          g *ᵥ (rep₂ (k := k) *ᵥ (gᵀ *ᵥ v)) := by
            simp [Wedge2Formalization.N6DoublePureSingular.ActBivector, Matrix.mulVec_mulVec,
              Matrix.mul_assoc]
        _ = 0 := by simp [hv0]
    have hfix2' :
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g =
          rep₂ (k := k) := by
      simpa [rep₂] using hfix.2
    rw [hfix2'] at hzero
    exact hzero
  exact hvne (commonKernel_zero (k := k) v hv1 hv2)

private theorem preserves_rightKernel
    (Ω :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hg : Matrix.det g ≠ 0)
    (hfix : Wedge2Formalization.N6DoublePureSingular.ActBivector Ω g = Ω)
    (v : Wedge2Formalization.N6DoublePureSingular.V → k)
    (hv : Ω *ᵥ v = 0) :
    Ω *ᵥ (gᵀ *ᵥ v) = 0 := by
  have hzero :
      Wedge2Formalization.N6DoublePureSingular.ActBivector Ω g *ᵥ v = 0 := by
    simpa [hfix] using hv
  have hzero' :
      g *ᵥ (Ω *ᵥ (gᵀ *ᵥ v)) = 0 := by
    simpa [Wedge2Formalization.N6DoublePureSingular.ActBivector, Matrix.mulVec_mulVec,
      Matrix.mul_assoc] using hzero
  exact Matrix.eq_zero_of_mulVec_eq_zero hg hzero'

private theorem rep₁_kernel_shape
    (x : Wedge2Formalization.N6DoublePureSingular.V → k)
    (hx : rep₁ (k := k) *ᵥ x = 0) :
    x (Sum.inl 0) = 0 ∧ x (Sum.inl 2) = 0 ∧
      x (Sum.inr 0) = 0 ∧ x (Sum.inr 2) = 0 := by
  have h00 := congrArg (fun v => v (Sum.inl 2)) hx
  have h02 := congrArg (fun v => v (Sum.inl 0)) hx
  have h30 := congrArg (fun v => v (Sum.inr 2)) hx
  have h32 := congrArg (fun v => v (Sum.inr 0)) hx
  refine ⟨?_, ?_, ?_, ?_⟩
  · simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h00
  · simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h02
  · simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h30
  · simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h32

private theorem rep₂_kernel_shape
    (x : Wedge2Formalization.N6DoublePureSingular.V → k)
    (hx : rep₂ (k := k) *ᵥ x = 0) :
    x (Sum.inl 1) = 0 ∧ x (Sum.inl 2) = 0 ∧
      x (Sum.inr 1) = 0 ∧ x (Sum.inr 2) = 0 := by
  have h01 := congrArg (fun v => v (Sum.inl 2)) hx
  have h02 := congrArg (fun v => v (Sum.inl 1)) hx
  have h31 := congrArg (fun v => v (Sum.inr 2)) hx
  have h32 := congrArg (fun v => v (Sum.inr 1)) hx
  refine ⟨?_, ?_, ?_, ?_⟩
  · simpa [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h01
  · simpa [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h02
  · simpa [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h31
  · simpa [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h32

private theorem repSum_kernel_shape
    (x : Wedge2Formalization.N6DoublePureSingular.V → k)
    (hx : (rep₁ (k := k) + rep₂ (k := k)) *ᵥ x = 0) :
    x (Sum.inl 2) = 0 ∧ x (Sum.inr 2) = 0 ∧
      x (Sum.inl 0) + x (Sum.inl 1) = 0 ∧
      x (Sum.inr 0) + x (Sum.inr 1) = 0 := by
  have h02 := congrArg (fun v => v (Sum.inl 0)) hx
  have h32 := congrArg (fun v => v (Sum.inr 0)) hx
  have h01 := congrArg (fun v => v (Sum.inl 2)) hx
  have h31 := congrArg (fun v => v (Sum.inr 2)) hx
  refine ⟨?_, ?_, ?_, ?_⟩
  · simpa [rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h02
  · simpa [rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three] using h32
  ·
    have h01' : -(x (Sum.inl 0) + x (Sum.inl 1)) = 0 := by
      simpa [rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
        Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
        dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, sub_eq_add_neg,
        add_comm, add_left_comm, add_assoc] using h01
    have := congrArg Neg.neg h01'
    simpa using this
  ·
    have h31' : -(x (Sum.inr 0) + x (Sum.inr 1)) = 0 := by
      simpa [rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
        Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
        Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
        dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, sub_eq_add_neg,
        add_comm, add_left_comm, add_assoc] using h31
    have := congrArg Neg.neg h31'
    simpa using this

private theorem transpose_mulVec_single
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (a b : Wedge2Formalization.N6DoublePureSingular.V) :
    (gᵀ *ᵥ Pi.single a (1 : k)) b = g a b := by
  rcases a with a | a <;> rcases b with b | b <;>
    fin_cases a <;> fin_cases b <;>
    simp [Matrix.mulVec, Matrix.transpose_apply, dotProduct, Fintype.sum_sum_type,
      Fin.sum_univ_three]

private theorem transpose_mulVec_single_sub_single
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (a b c : Wedge2Formalization.N6DoublePureSingular.V) :
    (gᵀ *ᵥ (Pi.single a (1 : k) - Pi.single b (1 : k))) c = g a c - g b c := by
  rw [mulVec_sub]
  simp [transpose_mulVec_single, Pi.sub_apply]

private structure Row11BlockEqns
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) : Prop where
  inv00 :
    g.toBlocks₁₁ 0 0 * g.toBlocks₁₁ 2 2 +
      g.toBlocks₁₂ 0 0 * g.toBlocks₁₂ 2 2 = 1
  inv01 :
    g.toBlocks₂₁ 0 0 * g.toBlocks₁₁ 2 2 +
      g.toBlocks₂₂ 0 0 * g.toBlocks₁₂ 2 2 = 0
  inv10 :
    g.toBlocks₁₁ 0 0 * g.toBlocks₂₁ 2 2 +
      g.toBlocks₁₂ 0 0 * g.toBlocks₂₂ 2 2 = 0
  inv11 :
    g.toBlocks₂₁ 0 0 * g.toBlocks₂₁ 2 2 +
      g.toBlocks₂₂ 0 0 * g.toBlocks₂₂ 2 2 = 1
  mid03 :
    g.toBlocks₁₁ 2 2 * g.toBlocks₂₁ 2 0 -
      g.toBlocks₁₁ 2 0 * g.toBlocks₂₁ 2 2 +
      g.toBlocks₁₂ 2 2 * g.toBlocks₂₂ 2 0 -
      g.toBlocks₁₂ 2 0 * g.toBlocks₂₂ 2 2 = 0
  mid14 :
    g.toBlocks₁₁ 2 2 * g.toBlocks₂₁ 2 1 -
      g.toBlocks₁₁ 2 1 * g.toBlocks₂₁ 2 2 +
      g.toBlocks₁₂ 2 2 * g.toBlocks₂₂ 2 1 -
      g.toBlocks₁₂ 2 1 * g.toBlocks₂₂ 2 2 = 0

private theorem ω13_cross_entry22_of_shapes
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hAshape :
      g.toBlocks₁₁ =
        !![g.toBlocks₁₁ 0 0, 0, 0;
            0, g.toBlocks₁₁ 0 0, 0;
            g.toBlocks₁₁ 2 0, g.toBlocks₁₁ 2 1, g.toBlocks₁₁ 2 2])
    (hBshape :
      g.toBlocks₁₂ =
        !![g.toBlocks₁₂ 0 0, 0, 0;
            0, g.toBlocks₁₂ 0 0, 0;
            g.toBlocks₁₂ 2 0, g.toBlocks₁₂ 2 1, g.toBlocks₁₂ 2 2])
    (hCshape :
      g.toBlocks₂₁ =
        !![g.toBlocks₂₁ 0 0, 0, 0;
            0, g.toBlocks₂₁ 0 0, 0;
            g.toBlocks₂₁ 2 0, g.toBlocks₂₁ 2 1, g.toBlocks₂₁ 2 2])
    (hDshape :
      g.toBlocks₂₂ =
        !![g.toBlocks₂₂ 0 0, 0, 0;
            0, g.toBlocks₂₂ 0 0, 0;
            g.toBlocks₂₂ 2 0, g.toBlocks₂₂ 2 1, g.toBlocks₂₂ 2 2]) :
    (g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₂ᵀ) 2 2 =
      g.toBlocks₁₁ 2 0 * g.toBlocks₂₁ 2 2 -
        g.toBlocks₁₁ 2 2 * g.toBlocks₂₁ 2 0 +
        g.toBlocks₁₂ 2 0 * g.toBlocks₂₂ 2 2 -
        g.toBlocks₁₂ 2 2 * g.toBlocks₂₂ 2 0 := by
  have hCshapeT :
      g.toBlocks₂₁ᵀ =
        !![g.toBlocks₂₁ 0 0, 0, g.toBlocks₂₁ 2 0;
            0, g.toBlocks₂₁ 0 0, g.toBlocks₂₁ 2 1;
            0, 0, g.toBlocks₂₁ 2 2] := by
    rw [hCshape]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.transpose_apply]
  have hDshapeT :
      g.toBlocks₂₂ᵀ =
        !![g.toBlocks₂₂ 0 0, 0, g.toBlocks₂₂ 2 0;
            0, g.toBlocks₂₂ 0 0, g.toBlocks₂₂ 2 1;
            0, 0, g.toBlocks₂₂ 2 2] := by
    rw [hDshape]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.transpose_apply]
  rw [hAshape, hBshape, hCshapeT, hDshapeT]
  simp [Wedge2Formalization.N3PureSingular.ω13, Matrix.mul_apply, Fin.sum_univ_three]
  ring

private theorem ω23_cross_entry22_of_shapes
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hAshape :
      g.toBlocks₁₁ =
        !![g.toBlocks₁₁ 0 0, 0, 0;
            0, g.toBlocks₁₁ 0 0, 0;
            g.toBlocks₁₁ 2 0, g.toBlocks₁₁ 2 1, g.toBlocks₁₁ 2 2])
    (hBshape :
      g.toBlocks₁₂ =
        !![g.toBlocks₁₂ 0 0, 0, 0;
            0, g.toBlocks₁₂ 0 0, 0;
            g.toBlocks₁₂ 2 0, g.toBlocks₁₂ 2 1, g.toBlocks₁₂ 2 2])
    (hCshape :
      g.toBlocks₂₁ =
        !![g.toBlocks₂₁ 0 0, 0, 0;
            0, g.toBlocks₂₁ 0 0, 0;
            g.toBlocks₂₁ 2 0, g.toBlocks₂₁ 2 1, g.toBlocks₂₁ 2 2])
    (hDshape :
      g.toBlocks₂₂ =
        !![g.toBlocks₂₂ 0 0, 0, 0;
            0, g.toBlocks₂₂ 0 0, 0;
            g.toBlocks₂₂ 2 0, g.toBlocks₂₂ 2 1, g.toBlocks₂₂ 2 2]) :
    (g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
        g.toBlocks₁₂ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₂ᵀ) 2 2 =
      g.toBlocks₁₁ 2 1 * g.toBlocks₂₁ 2 2 -
        g.toBlocks₁₁ 2 2 * g.toBlocks₂₁ 2 1 +
        g.toBlocks₁₂ 2 1 * g.toBlocks₂₂ 2 2 -
        g.toBlocks₁₂ 2 2 * g.toBlocks₂₂ 2 1 := by
  have hCshapeT :
      g.toBlocks₂₁ᵀ =
        !![g.toBlocks₂₁ 0 0, 0, g.toBlocks₂₁ 2 0;
            0, g.toBlocks₂₁ 0 0, g.toBlocks₂₁ 2 1;
            0, 0, g.toBlocks₂₁ 2 2] := by
    rw [hCshape]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.transpose_apply]
  have hDshapeT :
      g.toBlocks₂₂ᵀ =
        !![g.toBlocks₂₂ 0 0, 0, g.toBlocks₂₂ 2 0;
            0, g.toBlocks₂₂ 0 0, g.toBlocks₂₂ 2 1;
            0, 0, g.toBlocks₂₂ 2 2] := by
    rw [hDshape]
    ext i j
    fin_cases i <;> fin_cases j <;> simp [Matrix.transpose_apply]
  rw [hAshape, hBshape, hCshapeT, hDshapeT]
  simp [Wedge2Formalization.N3PureSingular.ω23, Matrix.mul_apply, Fin.sum_univ_three]
  ring

set_option maxRecDepth 5000 in
private theorem blockEqns_of_shapes
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hAshape :
      g.toBlocks₁₁ =
        !![g.toBlocks₁₁ 0 0, 0, 0;
            0, g.toBlocks₁₁ 0 0, 0;
            g.toBlocks₁₁ 2 0, g.toBlocks₁₁ 2 1, g.toBlocks₁₁ 2 2])
    (hBshape :
      g.toBlocks₁₂ =
        !![g.toBlocks₁₂ 0 0, 0, 0;
            0, g.toBlocks₁₂ 0 0, 0;
            g.toBlocks₁₂ 2 0, g.toBlocks₁₂ 2 1, g.toBlocks₁₂ 2 2])
    (hCshape :
      g.toBlocks₂₁ =
        !![g.toBlocks₂₁ 0 0, 0, 0;
            0, g.toBlocks₂₁ 0 0, 0;
            g.toBlocks₂₁ 2 0, g.toBlocks₂₁ 2 1, g.toBlocks₂₁ 2 2])
    (hDshape :
      g.toBlocks₂₂ =
        !![g.toBlocks₂₂ 0 0, 0, 0;
            0, g.toBlocks₂₂ 0 0, 0;
            g.toBlocks₂₂ 2 0, g.toBlocks₂₂ 2 1, g.toBlocks₂₂ 2 2])
    (hfix1' :
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₁ (k := k))
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) =
        rep₁ (k := k))
    (hfix2' :
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₂ (k := k))
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂) =
        rep₂ (k := k)) :
    Row11BlockEqns (k := k) g := by
  have hact1 :=
    Wedge2Formalization.N6DoublePureSingular.act_diagPair_fromBlocks_local
      (k := k) (Ω := Wedge2Formalization.N3PureSingular.ω13 (k := k))
      g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂
  have hact2 :=
    Wedge2Formalization.N6DoublePureSingular.act_diagPair_fromBlocks_local
      (k := k) (Ω := Wedge2Formalization.N3PureSingular.ω23 (k := k))
      g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂
  refine ⟨?_, ?_, ?_, ?_, ?_, ?_⟩
  ·
    have hcoord := congrArg
      (fun M => (Matrix.toBlocks₁₁ M) 0 2) (hact1.symm.trans hfix1')
    rw [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      hAshape, hBshape, hCshape, hDshape] at hcoord
    simpa [Wedge2Formalization.N3PureSingular.ω13,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using hcoord
  ·
    have hcoord := congrArg
      (fun M => (Matrix.toBlocks₂₁ M) 0 2) (hact1.symm.trans hfix1')
    rw [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      hAshape, hBshape, hCshape, hDshape] at hcoord
    simpa [Wedge2Formalization.N3PureSingular.ω13,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using hcoord
  ·
    have hcoord := congrArg
      (fun M => (Matrix.toBlocks₁₂ M) 0 2) (hact1.symm.trans hfix1')
    rw [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      hAshape, hBshape, hCshape, hDshape] at hcoord
    simpa [Wedge2Formalization.N3PureSingular.ω13,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using hcoord
  ·
    have hcoord := congrArg
      (fun M => (Matrix.toBlocks₂₂ M) 0 2) (hact1.symm.trans hfix1')
    rw [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
      hAshape, hBshape, hCshape, hDshape] at hcoord
    simpa [Wedge2Formalization.N3PureSingular.ω13,
      Matrix.fromBlocks, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_three] using hcoord
  ·
    have hcoord := congrArg
      (fun M => M (Sum.inl 2) (Sum.inr 2)) (hact1.symm.trans hfix1')
    have hentry :
        (g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₁ᵀ +
            g.toBlocks₁₂ * Wedge2Formalization.N3PureSingular.ω13 (k := k) * g.toBlocks₂₂ᵀ) 2 2 = 0 := by
      simpa [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁, Matrix.fromBlocks] using hcoord
    rw [ω13_cross_entry22_of_shapes (k := k) (g := g) hAshape hBshape hCshape hDshape] at hentry
    have hentry' := congrArg Neg.neg hentry
    ring_nf at hentry' ⊢
    exact hentry'
  ·
    have hcoord := congrArg
      (fun M => M (Sum.inl 2) (Sum.inr 2)) (hact2.symm.trans hfix2')
    have hentry :
        (g.toBlocks₁₁ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₁ᵀ +
            g.toBlocks₁₂ * Wedge2Formalization.N3PureSingular.ω23 (k := k) * g.toBlocks₂₂ᵀ) 2 2 = 0 := by
      simpa [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂, Matrix.fromBlocks] using hcoord
    rw [ω23_cross_entry22_of_shapes (k := k) (g := g) hAshape hBshape hCshape hDshape] at hentry
    have hentry' := congrArg Neg.neg hentry
    ring_nf at hentry' ⊢
    exact hentry'

private structure Row11ParamEqns
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) : Prop where
  hΔ :
    g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
      g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0 ≠ 0
  ht0 :
    g.toBlocks₁₁ 2 2 =
      g.toBlocks₂₂ 0 0 /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
  ht1 :
    g.toBlocks₁₂ 2 2 =
      (-g.toBlocks₂₁ 0 0) /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
  hu0 :
    g.toBlocks₂₁ 2 2 =
      (-g.toBlocks₁₂ 0 0) /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
  hu1 :
    g.toBlocks₂₂ 2 2 =
      g.toBlocks₁₁ 0 0 /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
  hqq :
    g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 +
      g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0 =
    g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 0 +
      g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 0
  huu :
    g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 +
      g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1 =
    g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 1 +
      g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 1

private theorem blockEqns_hΔ
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
        g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0 ≠ 0 := by
  let M : Matrix (Fin 2) (Fin 2) k :=
    !![g.toBlocks₁₁ 0 0, g.toBlocks₁₂ 0 0; g.toBlocks₂₁ 0 0, g.toBlocks₂₂ 0 0]
  let N : Matrix (Fin 2) (Fin 2) k :=
    !![g.toBlocks₁₁ 2 2, g.toBlocks₂₁ 2 2; g.toBlocks₁₂ 2 2, g.toBlocks₂₂ 2 2]
  have hMN : M * N = 1 := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simpa [M, N, hEqns.inv00, hEqns.inv01, hEqns.inv10, hEqns.inv11]
  intro hΔ0
  have hdetMN := congrArg Matrix.det hMN
  have hdetM : Matrix.det M = 0 := by
    simp [M, Matrix.det_fin_two, hΔ0]
  have : (0 : k) = 1 := by
    calc
      (0 : k) = Matrix.det M * Matrix.det N := by rw [hdetM, zero_mul]
      _ = 1 := by simpa [Matrix.det_one] using hdetMN
  exact zero_ne_one this

private theorem blockEqns_ht0
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    g.toBlocks₁₁ 2 2 =
      g.toBlocks₂₂ 0 0 /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0) := by
  have hΔ := blockEqns_hΔ (k := k) (g := g) hEqns
  apply (eq_div_iff hΔ).2
  calc
    g.toBlocks₁₁ 2 2 *
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
        =
      g.toBlocks₂₂ 0 0 * (g.toBlocks₁₁ 0 0 * g.toBlocks₁₁ 2 2 +
        g.toBlocks₁₂ 0 0 * g.toBlocks₁₂ 2 2) -
      g.toBlocks₁₂ 0 0 * (g.toBlocks₂₁ 0 0 * g.toBlocks₁₁ 2 2 +
        g.toBlocks₂₂ 0 0 * g.toBlocks₁₂ 2 2) := by ring
    _ = g.toBlocks₂₂ 0 0 := by rw [hEqns.inv00, hEqns.inv01]; ring

private theorem blockEqns_ht1
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    g.toBlocks₁₂ 2 2 =
      (-g.toBlocks₂₁ 0 0) /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0) := by
  have hΔ := blockEqns_hΔ (k := k) (g := g) hEqns
  apply (eq_div_iff hΔ).2
  calc
    g.toBlocks₁₂ 2 2 *
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
        =
      g.toBlocks₁₁ 0 0 * (g.toBlocks₂₁ 0 0 * g.toBlocks₁₁ 2 2 +
        g.toBlocks₂₂ 0 0 * g.toBlocks₁₂ 2 2) -
      g.toBlocks₂₁ 0 0 * (g.toBlocks₁₁ 0 0 * g.toBlocks₁₁ 2 2 +
        g.toBlocks₁₂ 0 0 * g.toBlocks₁₂ 2 2) := by ring
    _ = -g.toBlocks₂₁ 0 0 := by rw [hEqns.inv01, hEqns.inv00]; ring

private theorem blockEqns_hu0
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    g.toBlocks₂₁ 2 2 =
      (-g.toBlocks₁₂ 0 0) /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0) := by
  have hΔ := blockEqns_hΔ (k := k) (g := g) hEqns
  apply (eq_div_iff hΔ).2
  calc
    g.toBlocks₂₁ 2 2 *
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
        =
      g.toBlocks₂₂ 0 0 * (g.toBlocks₁₁ 0 0 * g.toBlocks₂₁ 2 2 +
        g.toBlocks₁₂ 0 0 * g.toBlocks₂₂ 2 2) -
      g.toBlocks₁₂ 0 0 * (g.toBlocks₂₁ 0 0 * g.toBlocks₂₁ 2 2 +
        g.toBlocks₂₂ 0 0 * g.toBlocks₂₂ 2 2) := by ring
    _ = -g.toBlocks₁₂ 0 0 := by rw [hEqns.inv10, hEqns.inv11]; ring

private theorem blockEqns_hu1
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    g.toBlocks₂₂ 2 2 =
      g.toBlocks₁₁ 0 0 /
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0) := by
  have hΔ := blockEqns_hΔ (k := k) (g := g) hEqns
  apply (eq_div_iff hΔ).2
  calc
    g.toBlocks₂₂ 2 2 *
        (g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
          g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0)
        =
      g.toBlocks₁₁ 0 0 * (g.toBlocks₂₁ 0 0 * g.toBlocks₂₁ 2 2 +
        g.toBlocks₂₂ 0 0 * g.toBlocks₂₂ 2 2) -
      g.toBlocks₂₁ 0 0 * (g.toBlocks₁₁ 0 0 * g.toBlocks₂₁ 2 2 +
        g.toBlocks₁₂ 0 0 * g.toBlocks₂₂ 2 2) := by ring
    _ = g.toBlocks₁₁ 0 0 := by rw [hEqns.inv11, hEqns.inv10]; ring

private theorem blockEqns_hqq
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 +
        g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0 =
      g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 0 +
        g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 0 := by
  let Δ : k :=
    g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
      g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0
  have hΔ := blockEqns_hΔ (k := k) (g := g) hEqns
  have ht0 := blockEqns_ht0 (k := k) (g := g) hEqns
  have ht1 := blockEqns_ht1 (k := k) (g := g) hEqns
  have hu0 := blockEqns_hu0 (k := k) (g := g) hEqns
  have hu1 := blockEqns_hu1 (k := k) (g := g) hEqns
  have htmp := hEqns.mid03
  rw [ht0, ht1, hu0, hu1] at htmp
  have hnum :
      g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0 +
          g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 -
          g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 0 -
          g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 0 = 0 := by
    have hexpr :
        g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0 +
            g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 -
            g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 0 -
            g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 0 =
          Δ *
            (g.toBlocks₂₂ 0 0 / Δ * g.toBlocks₂₁ 2 0 -
              g.toBlocks₁₁ 2 0 * (-g.toBlocks₁₂ 0 0 / Δ) +
              (-g.toBlocks₂₁ 0 0 / Δ) * g.toBlocks₂₂ 2 0 -
              g.toBlocks₁₂ 2 0 * (g.toBlocks₁₁ 0 0 / Δ)) := by
      field_simp [Δ, hΔ]
      ring
    calc
      g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0 +
          g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 -
          g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 0 -
          g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 0
          =
        Δ *
          (g.toBlocks₂₂ 0 0 / Δ * g.toBlocks₂₁ 2 0 -
            g.toBlocks₁₁ 2 0 * (-g.toBlocks₁₂ 0 0 / Δ) +
            (-g.toBlocks₂₁ 0 0 / Δ) * g.toBlocks₂₂ 2 0 -
            g.toBlocks₁₂ 2 0 * (g.toBlocks₁₁ 0 0 / Δ)) := hexpr
      _ = 0 := by rw [htmp, mul_zero]
  apply sub_eq_zero.mp
  simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hnum

private theorem blockEqns_huu
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 +
        g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1 =
      g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 1 +
        g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 1 := by
  let Δ : k :=
    g.toBlocks₁₁ 0 0 * g.toBlocks₂₂ 0 0 -
      g.toBlocks₁₂ 0 0 * g.toBlocks₂₁ 0 0
  have hΔ := blockEqns_hΔ (k := k) (g := g) hEqns
  have ht0 := blockEqns_ht0 (k := k) (g := g) hEqns
  have ht1 := blockEqns_ht1 (k := k) (g := g) hEqns
  have hu0 := blockEqns_hu0 (k := k) (g := g) hEqns
  have hu1 := blockEqns_hu1 (k := k) (g := g) hEqns
  have htmp := hEqns.mid14
  rw [ht0, ht1, hu0, hu1] at htmp
  have hnum :
      g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1 +
          g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 -
          g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 1 -
          g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 1 = 0 := by
    have hexpr :
        g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1 +
            g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 -
            g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 1 -
            g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 1 =
          Δ *
            (g.toBlocks₂₂ 0 0 / Δ * g.toBlocks₂₁ 2 1 -
              g.toBlocks₁₁ 2 1 * (-g.toBlocks₁₂ 0 0 / Δ) +
              (-g.toBlocks₂₁ 0 0 / Δ) * g.toBlocks₂₂ 2 1 -
              g.toBlocks₁₂ 2 1 * (g.toBlocks₁₁ 0 0 / Δ)) := by
      field_simp [Δ, hΔ]
      ring
    calc
      g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1 +
          g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 -
          g.toBlocks₂₁ 0 0 * g.toBlocks₂₂ 2 1 -
          g.toBlocks₁₁ 0 0 * g.toBlocks₁₂ 2 1
          =
        Δ *
          (g.toBlocks₂₂ 0 0 / Δ * g.toBlocks₂₁ 2 1 -
            g.toBlocks₁₁ 2 1 * (-g.toBlocks₁₂ 0 0 / Δ) +
            (-g.toBlocks₂₁ 0 0 / Δ) * g.toBlocks₂₂ 2 1 -
            g.toBlocks₁₂ 2 1 * (g.toBlocks₁₁ 0 0 / Δ)) := hexpr
      _ = 0 := by rw [htmp, mul_zero]
  apply sub_eq_zero.mp
  simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hnum

private theorem paramEqns_of_blockEqns
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hEqns : Row11BlockEqns (k := k) g) :
    Row11ParamEqns (k := k) g := by
  exact ⟨blockEqns_hΔ (k := k) (g := g) hEqns,
    blockEqns_ht0 (k := k) (g := g) hEqns,
    blockEqns_ht1 (k := k) (g := g) hEqns,
    blockEqns_hu0 (k := k) (g := g) hEqns,
    blockEqns_hu1 (k := k) (g := g) hEqns,
    blockEqns_hqq (k := k) (g := g) hEqns,
    blockEqns_huu (k := k) (g := g) hEqns⟩

private structure Row11Shapes
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) : Prop where
  Ashape :
    g.toBlocks₁₁ =
      !![g.toBlocks₁₁ 0 0, 0, 0;
          0, g.toBlocks₁₁ 0 0, 0;
          g.toBlocks₁₁ 2 0, g.toBlocks₁₁ 2 1, g.toBlocks₁₁ 2 2]
  Bshape :
    g.toBlocks₁₂ =
      !![g.toBlocks₁₂ 0 0, 0, 0;
          0, g.toBlocks₁₂ 0 0, 0;
          g.toBlocks₁₂ 2 0, g.toBlocks₁₂ 2 1, g.toBlocks₁₂ 2 2]
  Cshape :
    g.toBlocks₂₁ =
      !![g.toBlocks₂₁ 0 0, 0, 0;
          0, g.toBlocks₂₁ 0 0, 0;
          g.toBlocks₂₁ 2 0, g.toBlocks₂₁ 2 1, g.toBlocks₂₁ 2 2]
  Dshape :
    g.toBlocks₂₂ =
      !![g.toBlocks₂₂ 0 0, 0, 0;
          0, g.toBlocks₂₂ 0 0, 0;
          g.toBlocks₂₂ 2 0, g.toBlocks₂₂ 2 1, g.toBlocks₂₂ 2 2]

private structure Row11ZeroEntries
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) : Prop where
  z01 : g (Sum.inl 0) (Sum.inl 1) = 0
  z02 : g (Sum.inl 0) (Sum.inl 2) = 0
  z04 : g (Sum.inl 0) (Sum.inr 1) = 0
  z05 : g (Sum.inl 0) (Sum.inr 2) = 0
  z10 : g (Sum.inl 1) (Sum.inl 0) = 0
  z12 : g (Sum.inl 1) (Sum.inl 2) = 0
  z13 : g (Sum.inl 1) (Sum.inr 0) = 0
  z15 : g (Sum.inl 1) (Sum.inr 2) = 0
  z31 : g (Sum.inr 0) (Sum.inl 1) = 0
  z32 : g (Sum.inr 0) (Sum.inl 2) = 0
  z34 : g (Sum.inr 0) (Sum.inr 1) = 0
  z35 : g (Sum.inr 0) (Sum.inr 2) = 0
  z40 : g (Sum.inr 1) (Sum.inl 0) = 0
  z42 : g (Sum.inr 1) (Sum.inl 2) = 0
  z43 : g (Sum.inr 1) (Sum.inr 0) = 0
  z45 : g (Sum.inr 1) (Sum.inr 2) = 0

private structure Row11DiagEqns
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) : Prop where
  eq0011 : g (Sum.inl 0) (Sum.inl 0) = g (Sum.inl 1) (Sum.inl 1)
  eq0314 : g (Sum.inl 0) (Sum.inr 0) = g (Sum.inl 1) (Sum.inr 1)
  eq3041 : g (Sum.inr 0) (Sum.inl 0) = g (Sum.inr 1) (Sum.inl 1)
  eq3344 : g (Sum.inr 0) (Sum.inr 0) = g (Sum.inr 1) (Sum.inr 1)

private structure Row11EntryEqns
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k) : Prop where
  z01 : g (Sum.inl 0) (Sum.inl 1) = 0
  z02 : g (Sum.inl 0) (Sum.inl 2) = 0
  z04 : g (Sum.inl 0) (Sum.inr 1) = 0
  z05 : g (Sum.inl 0) (Sum.inr 2) = 0
  z10 : g (Sum.inl 1) (Sum.inl 0) = 0
  z12 : g (Sum.inl 1) (Sum.inl 2) = 0
  z13 : g (Sum.inl 1) (Sum.inr 0) = 0
  z15 : g (Sum.inl 1) (Sum.inr 2) = 0
  z31 : g (Sum.inr 0) (Sum.inl 1) = 0
  z32 : g (Sum.inr 0) (Sum.inl 2) = 0
  z34 : g (Sum.inr 0) (Sum.inr 1) = 0
  z35 : g (Sum.inr 0) (Sum.inr 2) = 0
  z40 : g (Sum.inr 1) (Sum.inl 0) = 0
  z42 : g (Sum.inr 1) (Sum.inl 2) = 0
  z43 : g (Sum.inr 1) (Sum.inr 0) = 0
  z45 : g (Sum.inr 1) (Sum.inr 2) = 0
  eq0011 : g (Sum.inl 0) (Sum.inl 0) = g (Sum.inl 1) (Sum.inl 1)
  eq0314 : g (Sum.inl 0) (Sum.inr 0) = g (Sum.inl 1) (Sum.inr 1)
  eq3041 : g (Sum.inr 0) (Sum.inl 0) = g (Sum.inr 1) (Sum.inl 1)
  eq3344 : g (Sum.inr 0) (Sum.inr 0) = g (Sum.inr 1) (Sum.inr 1)

set_option maxHeartbeats 20000000 in
private theorem zeroEntries_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g) :
    Row11ZeroEntries (k := k) g := by
  have hg : Matrix.det g ≠ 0 := det_ne_zero_of_pointwise (k := k) g hfix
  have hk0 :=
    preserves_rightKernel (k := k) (Ω := rep₂ (k := k)) (g := g) hg hfix.2
      (Pi.single (Sum.inl 0) (1 : k)) (by
        ext i
        rcases i with i | i <;> fin_cases i <;>
          simp [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
            Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
            dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three])
  have hk1 :=
    preserves_rightKernel (k := k) (Ω := rep₁ (k := k)) (g := g) hg hfix.1
      (Pi.single (Sum.inl 1) (1 : k)) (by
        ext i
        rcases i with i | i <;> fin_cases i <;>
          simp [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
            Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
            dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three])
  have hk3 :=
    preserves_rightKernel (k := k) (Ω := rep₂ (k := k)) (g := g) hg hfix.2
      (Pi.single (Sum.inr 0) (1 : k)) (by
        ext i
        rcases i with i | i <;> fin_cases i <;>
          simp [rep₂, Wedge2Formalization.N6DoublePureSingular.rep₂,
            Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
            dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three])
  have hk4 :=
    preserves_rightKernel (k := k) (Ω := rep₁ (k := k)) (g := g) hg hfix.1
      (Pi.single (Sum.inr 1) (1 : k)) (by
        ext i
        rcases i with i | i <;> fin_cases i <;>
          simp [rep₁, Wedge2Formalization.N6DoublePureSingular.rep₁,
            Wedge2Formalization.N3PureSingular.ω13, Matrix.mulVec, Matrix.fromBlocks,
            dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three])
  rcases rep₂_kernel_shape (k := k) (x := gᵀ *ᵥ Pi.single (Sum.inl 0) (1 : k)) hk0 with
    ⟨h01, h02, h04, h05⟩
  rcases rep₁_kernel_shape (k := k) (x := gᵀ *ᵥ Pi.single (Sum.inl 1) (1 : k)) hk1 with
    ⟨h10, h12, h13, h15⟩
  rcases rep₂_kernel_shape (k := k) (x := gᵀ *ᵥ Pi.single (Sum.inr 0) (1 : k)) hk3 with
    ⟨h31, h32, h34, h35⟩
  rcases rep₁_kernel_shape (k := k) (x := gᵀ *ᵥ Pi.single (Sum.inr 1) (1 : k)) hk4 with
    ⟨h40, h42, h43, h45⟩
  have h01' : g (Sum.inl 0) (Sum.inl 1) = 0 := by
    simpa [transpose_mulVec_single] using h01
  have h02' : g (Sum.inl 0) (Sum.inl 2) = 0 := by
    simpa [transpose_mulVec_single] using h02
  have h04' : g (Sum.inl 0) (Sum.inr 1) = 0 := by
    simpa [transpose_mulVec_single] using h04
  have h05' : g (Sum.inl 0) (Sum.inr 2) = 0 := by
    simpa [transpose_mulVec_single] using h05
  have h10' : g (Sum.inl 1) (Sum.inl 0) = 0 := by
    simpa [transpose_mulVec_single] using h10
  have h12' : g (Sum.inl 1) (Sum.inl 2) = 0 := by
    simpa [transpose_mulVec_single] using h12
  have h13' : g (Sum.inl 1) (Sum.inr 0) = 0 := by
    simpa [transpose_mulVec_single] using h13
  have h15' : g (Sum.inl 1) (Sum.inr 2) = 0 := by
    simpa [transpose_mulVec_single] using h15
  have h31' : g (Sum.inr 0) (Sum.inl 1) = 0 := by
    simpa [transpose_mulVec_single] using h31
  have h32' : g (Sum.inr 0) (Sum.inl 2) = 0 := by
    simpa [transpose_mulVec_single] using h32
  have h34' : g (Sum.inr 0) (Sum.inr 1) = 0 := by
    simpa [transpose_mulVec_single] using h34
  have h35' : g (Sum.inr 0) (Sum.inr 2) = 0 := by
    simpa [transpose_mulVec_single] using h35
  have h40' : g (Sum.inr 1) (Sum.inl 0) = 0 := by
    simpa [transpose_mulVec_single] using h40
  have h42' : g (Sum.inr 1) (Sum.inl 2) = 0 := by
    simpa [transpose_mulVec_single] using h42
  have h43' : g (Sum.inr 1) (Sum.inr 0) = 0 := by
    simpa [transpose_mulVec_single] using h43
  have h45' : g (Sum.inr 1) (Sum.inr 2) = 0 := by
    simpa [transpose_mulVec_single] using h45
  exact
    ⟨h01', h02', h04', h05', h10', h12', h13', h15',
      h31', h32', h34', h35', h40', h42', h43', h45'⟩

private theorem repSum_mulVec_leftDiff_zero :
    (rep₁ (k := k) + rep₂ (k := k)) *ᵥ
        (Pi.single (Sum.inl 0) (1 : k) - Pi.single (Sum.inl 1) (1 : k)) = 0 := by
  ext i
  rcases i with i | i <;> fin_cases i <;>
    simp [rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, sub_eq_add_neg, add_comm]

private theorem repSum_mulVec_rightDiff_zero :
    (rep₁ (k := k) + rep₂ (k := k)) *ᵥ
        (Pi.single (Sum.inr 0) (1 : k) - Pi.single (Sum.inr 1) (1 : k)) = 0 := by
  ext i
  rcases i with i | i <;> fin_cases i <;>
    simp [rep₁, rep₂, Wedge2Formalization.N6DoublePureSingular.rep₁,
      Wedge2Formalization.N6DoublePureSingular.rep₂, Wedge2Formalization.N3PureSingular.ω13,
      Wedge2Formalization.N3PureSingular.ω23, Matrix.mulVec, Matrix.fromBlocks,
      dotProduct, Fintype.sum_sum_type, Fin.sum_univ_three, sub_eq_add_neg, add_comm]

set_option maxHeartbeats 20000000 in
private theorem leftDiagEqns_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g)
    (hZ : Row11ZeroEntries (k := k) g) :
    g (Sum.inl 0) (Sum.inl 0) = g (Sum.inl 1) (Sum.inl 1) ∧
      g (Sum.inl 0) (Sum.inr 0) = g (Sum.inl 1) (Sum.inr 1) := by
  have hg : Matrix.det g ≠ 0 := det_ne_zero_of_pointwise (k := k) g hfix
  have hfix1' :
      Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g =
        rep₁ (k := k) := by
    simpa [rep₁] using hfix.1
  have hfix2' :
      Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g =
        rep₂ (k := k) := by
    simpa [rep₂] using hfix.2
  have hsum :
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₁ (k := k) + rep₂ (k := k)) g =
        rep₁ (k := k) + rep₂ (k := k) := by
    calc
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₁ (k := k) + rep₂ (k := k)) g
          =
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g +
          Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g := by
            ext i j
            simp [Wedge2Formalization.N6DoublePureSingular.ActBivector, Matrix.mul_add,
              Matrix.add_mul, Matrix.mul_assoc]
    _ = rep₁ (k := k) + rep₂ (k := k) := by rw [hfix1', hfix2']
  have hk01 :=
    preserves_rightKernel (k := k) (Ω := rep₁ (k := k) + rep₂ (k := k)) (g := g) hg hsum
      (Pi.single (Sum.inl 0) (1 : k) - Pi.single (Sum.inl 1) (1 : k))
      (repSum_mulVec_leftDiff_zero (k := k))
  rcases repSum_kernel_shape
      (k := k)
      (x := gᵀ *ᵥ (Pi.single (Sum.inl 0) (1 : k) - Pi.single (Sum.inl 1) (1 : k))) hk01 with
    ⟨_, _, h0011, h0314⟩
  have hx0011 :
      (g (Sum.inl 0) (Sum.inl 0) - g (Sum.inl 1) (Sum.inl 0)) +
        (g (Sum.inl 0) (Sum.inl 1) - g (Sum.inl 1) (Sum.inl 1)) = 0 := by
    simpa [transpose_mulVec_single_sub_single, hZ.z10, hZ.z01] using h0011
  have hx0314 :
      (g (Sum.inl 0) (Sum.inr 0) - g (Sum.inl 1) (Sum.inr 0)) +
        (g (Sum.inl 0) (Sum.inr 1) - g (Sum.inl 1) (Sum.inr 1)) = 0 := by
    simpa [transpose_mulVec_single_sub_single, hZ.z13, hZ.z04] using h0314
  have h0011' :
      g (Sum.inl 0) (Sum.inl 0) = g (Sum.inl 1) (Sum.inl 1) := by
    have hx : g (Sum.inl 0) (Sum.inl 0) - g (Sum.inl 1) (Sum.inl 1) = 0 := by
      simpa [hZ.z01, hZ.z10, sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hx0011
    exact sub_eq_zero.mp hx
  have h0314' :
      g (Sum.inl 0) (Sum.inr 0) = g (Sum.inl 1) (Sum.inr 1) := by
    have hx : g (Sum.inl 0) (Sum.inr 0) - g (Sum.inl 1) (Sum.inr 1) = 0 := by
      simpa [hZ.z04, hZ.z13, sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hx0314
    exact sub_eq_zero.mp hx
  exact ⟨h0011', h0314'⟩

set_option maxHeartbeats 20000000 in
private theorem rightDiagEqns_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g)
    (hZ : Row11ZeroEntries (k := k) g) :
    g (Sum.inr 0) (Sum.inl 0) = g (Sum.inr 1) (Sum.inl 1) ∧
      g (Sum.inr 0) (Sum.inr 0) = g (Sum.inr 1) (Sum.inr 1) := by
  have hg : Matrix.det g ≠ 0 := det_ne_zero_of_pointwise (k := k) g hfix
  have hfix1' :
      Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g =
        rep₁ (k := k) := by
    simpa [rep₁] using hfix.1
  have hfix2' :
      Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g =
        rep₂ (k := k) := by
    simpa [rep₂] using hfix.2
  have hsum :
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₁ (k := k) + rep₂ (k := k)) g =
        rep₁ (k := k) + rep₂ (k := k) := by
    calc
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₁ (k := k) + rep₂ (k := k)) g
          =
        Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₁ (k := k)) g +
          Wedge2Formalization.N6DoublePureSingular.ActBivector (rep₂ (k := k)) g := by
            ext i j
            simp [Wedge2Formalization.N6DoublePureSingular.ActBivector, Matrix.mul_add,
              Matrix.add_mul, Matrix.mul_assoc]
    _ = rep₁ (k := k) + rep₂ (k := k) := by rw [hfix1', hfix2']
  have hk34 :=
    preserves_rightKernel (k := k) (Ω := rep₁ (k := k) + rep₂ (k := k)) (g := g) hg hsum
      (Pi.single (Sum.inr 0) (1 : k) - Pi.single (Sum.inr 1) (1 : k))
      (repSum_mulVec_rightDiff_zero (k := k))
  rcases repSum_kernel_shape
      (k := k)
      (x := gᵀ *ᵥ (Pi.single (Sum.inr 0) (1 : k) - Pi.single (Sum.inr 1) (1 : k))) hk34 with
    ⟨_, _, h3041, h3344⟩
  have hx3041 :
      (g (Sum.inr 0) (Sum.inl 0) - g (Sum.inr 1) (Sum.inl 0)) +
        (g (Sum.inr 0) (Sum.inl 1) - g (Sum.inr 1) (Sum.inl 1)) = 0 := by
    simpa [transpose_mulVec_single_sub_single, hZ.z40, hZ.z31] using h3041
  have hx3344 :
      (g (Sum.inr 0) (Sum.inr 0) - g (Sum.inr 1) (Sum.inr 0)) +
        (g (Sum.inr 0) (Sum.inr 1) - g (Sum.inr 1) (Sum.inr 1)) = 0 := by
    simpa [transpose_mulVec_single_sub_single, hZ.z43, hZ.z34] using h3344
  have h3041' :
      g (Sum.inr 0) (Sum.inl 0) = g (Sum.inr 1) (Sum.inl 1) := by
    have hx : g (Sum.inr 0) (Sum.inl 0) - g (Sum.inr 1) (Sum.inl 1) = 0 := by
      simpa [hZ.z31, hZ.z40, sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hx3041
    exact sub_eq_zero.mp hx
  have h3344' :
      g (Sum.inr 0) (Sum.inr 0) = g (Sum.inr 1) (Sum.inr 1) := by
    have hx : g (Sum.inr 0) (Sum.inr 0) - g (Sum.inr 1) (Sum.inr 1) = 0 := by
      simpa [hZ.z34, hZ.z43, sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hx3344
    exact sub_eq_zero.mp hx
  exact ⟨h3041', h3344'⟩

set_option maxHeartbeats 20000000 in
private theorem diagEqns_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g)
    (hZ : Row11ZeroEntries (k := k) g) :
    Row11DiagEqns (k := k) g := by
  rcases leftDiagEqns_of_pointwise (k := k) (g := g) hfix hZ with ⟨h0011, h0314⟩
  rcases rightDiagEqns_of_pointwise (k := k) (g := g) hfix hZ with ⟨h3041, h3344⟩
  exact ⟨h0011, h0314, h3041, h3344⟩

private theorem entryEqns_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g) :
    Row11EntryEqns (k := k) g := by
  let hZ := zeroEntries_of_pointwise (k := k) (g := g) hfix
  let hD := diagEqns_of_pointwise (k := k) (g := g) hfix hZ
  exact
    ⟨hZ.z01, hZ.z02, hZ.z04, hZ.z05, hZ.z10, hZ.z12, hZ.z13, hZ.z15,
      hZ.z31, hZ.z32, hZ.z34, hZ.z35, hZ.z40, hZ.z42, hZ.z43, hZ.z45,
      hD.eq0011, hD.eq0314, hD.eq3041, hD.eq3344⟩

private theorem shapes_of_entryEqns
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hE : Row11EntryEqns (k := k) g) :
    Row11Shapes (k := k) g := by
  have hAshape :
      g.toBlocks₁₁ =
        !![g.toBlocks₁₁ 0 0, 0, 0;
            0, g.toBlocks₁₁ 0 0, 0;
            g.toBlocks₁₁ 2 0, g.toBlocks₁₁ 2 1, g.toBlocks₁₁ 2 2] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.toBlocks₁₁, hE.z01, hE.z02, hE.z10, hE.z12, hE.eq0011]
  have hBshape :
      g.toBlocks₁₂ =
        !![g.toBlocks₁₂ 0 0, 0, 0;
            0, g.toBlocks₁₂ 0 0, 0;
            g.toBlocks₁₂ 2 0, g.toBlocks₁₂ 2 1, g.toBlocks₁₂ 2 2] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.toBlocks₁₂, hE.z04, hE.z05, hE.z13, hE.z15, hE.eq0314]
  have hCshape :
      g.toBlocks₂₁ =
        !![g.toBlocks₂₁ 0 0, 0, 0;
            0, g.toBlocks₂₁ 0 0, 0;
            g.toBlocks₂₁ 2 0, g.toBlocks₂₁ 2 1, g.toBlocks₂₁ 2 2] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.toBlocks₂₁, hE.z31, hE.z32, hE.z40, hE.z42, hE.eq3041]
  have hDshape :
      g.toBlocks₂₂ =
        !![g.toBlocks₂₂ 0 0, 0, 0;
            0, g.toBlocks₂₂ 0 0, 0;
            g.toBlocks₂₂ 2 0, g.toBlocks₂₂ 2 1, g.toBlocks₂₂ 2 2] := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.toBlocks₂₂, hE.z34, hE.z35, hE.z43, hE.z45, hE.eq3344]
  exact ⟨hAshape, hBshape, hCshape, hDshape⟩

private theorem shapes_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g) :
    Row11Shapes (k := k) g := by
  exact shapes_of_entryEqns (k := k) (g := g) (entryEqns_of_pointwise (k := k) (g := g) hfix)

private theorem det_cancel_left
    (a b c d x y : k)
    (hΔ : a * d - b * c ≠ 0) :
    d / (a * d - b * c) * (a * x + c * y) +
        (-c) / (a * d - b * c) * (b * x + d * y) =
      x := by
  let Δ : k := a * d - b * c
  calc
    d / (a * d - b * c) * (a * x + c * y) +
        (-c) / (a * d - b * c) * (b * x + d * y)
        = (d * (a * x + c * y) + (-c) * (b * x + d * y)) * Δ⁻¹ := by
            simp [Δ, div_eq_mul_inv]
            ring
    _ = ((d * a - c * b) * x) * Δ⁻¹ := by
          ring
    _ = ((a * d - b * c) * x) * Δ⁻¹ := by
          congr 1
          ring
    _ = x * (Δ * Δ⁻¹) := by
          simp [Δ]
          ring
    _ = x := by
          simp [Δ, hΔ]

private theorem det_cancel_right
    (a b c d x y : k)
    (hΔ : a * d - b * c ≠ 0) :
    (-b) / (a * d - b * c) * (a * x + c * y) +
        a / (a * d - b * c) * (b * x + d * y) =
      y := by
  let Δ : k := a * d - b * c
  calc
    (-b) / (a * d - b * c) * (a * x + c * y) +
        a / (a * d - b * c) * (b * x + d * y)
        = ((-b) * (a * x + c * y) + a * (b * x + d * y)) * Δ⁻¹ := by
            simp [Δ, div_eq_mul_inv]
            ring
    _ = ((-(b * c) + a * d) * y) * Δ⁻¹ := by
          ring
    _ = ((a * d - b * c) * y) * Δ⁻¹ := by
          congr 1
          ring
    _ = y * (Δ * Δ⁻¹) := by
          simp [Δ]
          ring
    _ = y := by
          simp [Δ, hΔ]

private theorem pointwiseAC_formula
    (a b c d p q s u : k) :
    Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u =
    !![a, 0, 0;
        0, a, 0;
        d / (a * d - b * c) * p + (-c) / (a * d - b * c) * q,
        d / (a * d - b * c) * s + (-c) / (a * d - b * c) * u,
        d / (a * d - b * c)] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA,
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB,
      Wedge2Formalization.N6DoublePureSingular.scalarTri,
      Wedge2Formalization.N6DoublePureSingular.coupledA,
      Wedge2Formalization.N6DoublePureSingular.coupledC,
      Matrix.mul_apply, dotProduct, Fin.sum_univ_three]

private theorem pointwiseBD_formula
    (a b c d q r u v : k) :
    Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v =
    !![b, 0, 0;
        0, b, 0;
        d / (a * d - b * c) * q + (-c) / (a * d - b * c) * r,
        d / (a * d - b * c) * u + (-c) / (a * d - b * c) * v,
        (-c) / (a * d - b * c)] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA,
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB,
      Wedge2Formalization.N6DoublePureSingular.scalarTri,
      Wedge2Formalization.N6DoublePureSingular.coupledB,
      Wedge2Formalization.N6DoublePureSingular.coupledD,
      Matrix.mul_apply, dotProduct, Fin.sum_univ_three]

private theorem pointwiseCA_formula
    (a b c d p q s u : k) :
    Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u =
    !![c, 0, 0;
        0, c, 0;
        (-b) / (a * d - b * c) * p + a / (a * d - b * c) * q,
        (-b) / (a * d - b * c) * s + a / (a * d - b * c) * u,
        (-b) / (a * d - b * c)] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC,
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD,
      Wedge2Formalization.N6DoublePureSingular.scalarTri,
      Wedge2Formalization.N6DoublePureSingular.coupledA,
      Wedge2Formalization.N6DoublePureSingular.coupledC,
      Matrix.mul_apply, dotProduct, Fin.sum_univ_three]

private theorem pointwiseDD_formula
    (a b c d q r u v : k) :
    Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
        Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v =
    !![d, 0, 0;
        0, d, 0;
        (-b) / (a * d - b * c) * q + a / (a * d - b * c) * r,
        (-b) / (a * d - b * c) * u + a / (a * d - b * c) * v,
        a / (a * d - b * c)] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC,
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD,
      Wedge2Formalization.N6DoublePureSingular.scalarTri,
      Wedge2Formalization.N6DoublePureSingular.coupledB,
      Wedge2Formalization.N6DoublePureSingular.coupledD,
      Matrix.mul_apply, dotProduct, Fin.sum_univ_three]

private theorem Ablock_of_shapes_params
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hShapes : Row11Shapes (k := k) g)
    (hParams : Row11ParamEqns (k := k) g) :
    g.toBlocks₁₁ =
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledA
          (k := k)
          (g.toBlocks₁₁ 0 0 * g.toBlocks₁₁ 2 0 + g.toBlocks₂₁ 0 0 * g.toBlocks₂₁ 2 0)
          (g.toBlocks₁₁ 0 0 * g.toBlocks₁₁ 2 1 + g.toBlocks₂₁ 0 0 * g.toBlocks₂₁ 2 1) +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledC
          (k := k)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1) := by
  let a : k := g.toBlocks₁₁ 0 0
  let b : k := g.toBlocks₁₂ 0 0
  let c : k := g.toBlocks₂₁ 0 0
  let d : k := g.toBlocks₂₂ 0 0
  let p : k := a * g.toBlocks₁₁ 2 0 + c * g.toBlocks₂₁ 2 0
  let q : k := b * g.toBlocks₁₁ 2 0 + d * g.toBlocks₂₁ 2 0
  let s : k := a * g.toBlocks₁₁ 2 1 + c * g.toBlocks₂₁ 2 1
  let u : k := b * g.toBlocks₁₁ 2 1 + d * g.toBlocks₂₁ 2 1
  calc
    g.toBlocks₁₁ = !![g.toBlocks₁₁ 0 0, 0, 0;
        0, g.toBlocks₁₁ 0 0, 0;
        g.toBlocks₁₁ 2 0, g.toBlocks₁₁ 2 1, g.toBlocks₁₁ 2 2] := hShapes.Ashape
    _ = !![a, 0, 0;
        0, a, 0;
        d / (a * d - b * c) * p + (-c) / (a * d - b * c) * q,
        d / (a * d - b * c) * s + (-c) / (a * d - b * c) * u,
        d / (a * d - b * c)] := by
          ext i j
          fin_cases i <;> fin_cases j
          · simp [a]
          · simp
          · simp
          · simp
          · simp [a]
          · simp
          · symm
            simpa [p, q, a, b, c, d] using
              (det_cancel_left
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₁ 2 0) (y := g.toBlocks₂₁ 2 0) hParams.hΔ)
          · symm
            simpa [s, u, a, b, c, d] using
              (det_cancel_left
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₁ 2 1) (y := g.toBlocks₂₁ 2 1) hParams.hΔ)
          · simpa [a, b, c, d] using hParams.ht0
    _ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u := by
          exact (pointwiseAC_formula (k := k) a b c d p q s u).symm

private theorem Bblock_of_shapes_params
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hShapes : Row11Shapes (k := k) g)
    (hParams : Row11ParamEqns (k := k) g) :
    g.toBlocks₁₂ =
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledB
          (k := k)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1) +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledD
          (k := k)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₂ 2 0 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₂ 2 0)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₂ 2 1 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₂ 2 1) := by
  let a : k := g.toBlocks₁₁ 0 0
  let b : k := g.toBlocks₁₂ 0 0
  let c : k := g.toBlocks₂₁ 0 0
  let d : k := g.toBlocks₂₂ 0 0
  let q : k := b * g.toBlocks₁₁ 2 0 + d * g.toBlocks₂₁ 2 0
  let r : k := b * g.toBlocks₁₂ 2 0 + d * g.toBlocks₂₂ 2 0
  let u : k := b * g.toBlocks₁₁ 2 1 + d * g.toBlocks₂₁ 2 1
  let v : k := b * g.toBlocks₁₂ 2 1 + d * g.toBlocks₂₂ 2 1
  have hq0 :
      q = a * g.toBlocks₁₂ 2 0 + c * g.toBlocks₂₂ 2 0 := by
    simpa [a, b, c, d, q] using hParams.hqq
  have hu0 :
      u = a * g.toBlocks₁₂ 2 1 + c * g.toBlocks₂₂ 2 1 := by
    simpa [a, b, c, d, u] using hParams.huu
  calc
    g.toBlocks₁₂ = !![g.toBlocks₁₂ 0 0, 0, 0;
        0, g.toBlocks₁₂ 0 0, 0;
        g.toBlocks₁₂ 2 0, g.toBlocks₁₂ 2 1, g.toBlocks₁₂ 2 2] := hShapes.Bshape
    _ = !![b, 0, 0;
        0, b, 0;
        d / (a * d - b * c) * q + (-c) / (a * d - b * c) * r,
        d / (a * d - b * c) * u + (-c) / (a * d - b * c) * v,
        (-c) / (a * d - b * c)] := by
          ext i j
          fin_cases i <;> fin_cases j
          · simp [b]
          · simp
          · simp
          · simp
          · simp [b]
          · simp
          · rw [hq0]
            symm
            simpa [r, a, b, c, d] using
              (det_cancel_left
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₂ 2 0) (y := g.toBlocks₂₂ 2 0) hParams.hΔ)
          · rw [hu0]
            symm
            simpa [v, a, b, c, d] using
              (det_cancel_left
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₂ 2 1) (y := g.toBlocks₂₂ 2 1) hParams.hΔ)
          · simpa [a, b, c, d] using hParams.ht1
    _ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v := by
          exact (pointwiseBD_formula (k := k) a b c d q r u v).symm

private theorem Cblock_of_shapes_params
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hShapes : Row11Shapes (k := k) g)
    (hParams : Row11ParamEqns (k := k) g) :
    g.toBlocks₂₁ =
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledA
          (k := k)
          (g.toBlocks₁₁ 0 0 * g.toBlocks₁₁ 2 0 + g.toBlocks₂₁ 0 0 * g.toBlocks₂₁ 2 0)
          (g.toBlocks₁₁ 0 0 * g.toBlocks₁₁ 2 1 + g.toBlocks₂₁ 0 0 * g.toBlocks₂₁ 2 1) +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledC
          (k := k)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1) := by
  let a : k := g.toBlocks₁₁ 0 0
  let b : k := g.toBlocks₁₂ 0 0
  let c : k := g.toBlocks₂₁ 0 0
  let d : k := g.toBlocks₂₂ 0 0
  let p : k := a * g.toBlocks₁₁ 2 0 + c * g.toBlocks₂₁ 2 0
  let q : k := b * g.toBlocks₁₁ 2 0 + d * g.toBlocks₂₁ 2 0
  let s : k := a * g.toBlocks₁₁ 2 1 + c * g.toBlocks₂₁ 2 1
  let u : k := b * g.toBlocks₁₁ 2 1 + d * g.toBlocks₂₁ 2 1
  calc
    g.toBlocks₂₁ = !![g.toBlocks₂₁ 0 0, 0, 0;
        0, g.toBlocks₂₁ 0 0, 0;
        g.toBlocks₂₁ 2 0, g.toBlocks₂₁ 2 1, g.toBlocks₂₁ 2 2] := hShapes.Cshape
    _ = !![c, 0, 0;
        0, c, 0;
        (-b) / (a * d - b * c) * p + a / (a * d - b * c) * q,
        (-b) / (a * d - b * c) * s + a / (a * d - b * c) * u,
        (-b) / (a * d - b * c)] := by
          ext i j
          fin_cases i <;> fin_cases j
          · simp [c]
          · simp
          · simp
          · simp
          · simp [c]
          · simp
          · symm
            simpa [p, q, a, b, c, d] using
              (det_cancel_right
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₁ 2 0) (y := g.toBlocks₂₁ 2 0) hParams.hΔ)
          · symm
            simpa [s, u, a, b, c, d] using
              (det_cancel_right
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₁ 2 1) (y := g.toBlocks₂₁ 2 1) hParams.hΔ)
          · simpa [a, b, c, d] using hParams.hu0
    _ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u := by
          exact (pointwiseCA_formula (k := k) a b c d p q s u).symm

private theorem Dblock_of_shapes_params
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hShapes : Row11Shapes (k := k) g)
    (hParams : Row11ParamEqns (k := k) g) :
    g.toBlocks₂₂ =
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledB
          (k := k)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 0 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 0)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₁ 2 1 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₁ 2 1) +
      Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD
          (k := k) (g.toBlocks₁₁ 0 0) (g.toBlocks₁₂ 0 0) (g.toBlocks₂₁ 0 0) (g.toBlocks₂₂ 0 0) *
        Wedge2Formalization.N6DoublePureSingular.coupledD
          (k := k)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₂ 2 0 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₂ 2 0)
          (g.toBlocks₁₂ 0 0 * g.toBlocks₁₂ 2 1 + g.toBlocks₂₂ 0 0 * g.toBlocks₂₂ 2 1) := by
  let a : k := g.toBlocks₁₁ 0 0
  let b : k := g.toBlocks₁₂ 0 0
  let c : k := g.toBlocks₂₁ 0 0
  let d : k := g.toBlocks₂₂ 0 0
  let q : k := b * g.toBlocks₁₁ 2 0 + d * g.toBlocks₂₁ 2 0
  let r : k := b * g.toBlocks₁₂ 2 0 + d * g.toBlocks₂₂ 2 0
  let u : k := b * g.toBlocks₁₁ 2 1 + d * g.toBlocks₂₁ 2 1
  let v : k := b * g.toBlocks₁₂ 2 1 + d * g.toBlocks₂₂ 2 1
  have hq0 :
      q = a * g.toBlocks₁₂ 2 0 + c * g.toBlocks₂₂ 2 0 := by
    simpa [a, b, c, d, q] using hParams.hqq
  have hu0 :
      u = a * g.toBlocks₁₂ 2 1 + c * g.toBlocks₂₂ 2 1 := by
    simpa [a, b, c, d, u] using hParams.huu
  calc
    g.toBlocks₂₂ = !![g.toBlocks₂₂ 0 0, 0, 0;
        0, g.toBlocks₂₂ 0 0, 0;
        g.toBlocks₂₂ 2 0, g.toBlocks₂₂ 2 1, g.toBlocks₂₂ 2 2] := hShapes.Dshape
    _ = !![d, 0, 0;
        0, d, 0;
        (-b) / (a * d - b * c) * q + a / (a * d - b * c) * r,
        (-b) / (a * d - b * c) * u + a / (a * d - b * c) * v,
        a / (a * d - b * c)] := by
          ext i j
          fin_cases i <;> fin_cases j
          · simp [d]
          · simp
          · simp
          · simp
          · simp [d]
          · simp
          · rw [hq0]
            symm
            simpa [r, a, b, c, d] using
              (det_cancel_right
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₂ 2 0) (y := g.toBlocks₂₂ 2 0) hParams.hΔ)
          · rw [hu0]
            symm
            simpa [v, a, b, c, d] using
              (det_cancel_right
                (a := g.toBlocks₁₁ 0 0) (b := g.toBlocks₁₂ 0 0)
                (c := g.toBlocks₂₁ 0 0) (d := g.toBlocks₂₂ 0 0)
                (x := g.toBlocks₁₂ 2 1) (y := g.toBlocks₂₂ 2 1) hParams.hΔ)
          · simpa [a, b, c, d] using hParams.hu1
    _ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v := by
          exact (pointwiseDD_formula (k := k) a b c d q r u v).symm

set_option maxHeartbeats 20000000 in
private theorem mem_K_of_shapes_params
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hShapes : Row11Shapes (k := k) g)
    (hParams : Row11ParamEqns (k := k) g) :
    g ∈ K (k := k) := by
  let a : k := g.toBlocks₁₁ 0 0
  let b : k := g.toBlocks₁₂ 0 0
  let c : k := g.toBlocks₂₁ 0 0
  let d : k := g.toBlocks₂₂ 0 0
  let p : k := a * g.toBlocks₁₁ 2 0 + c * g.toBlocks₂₁ 2 0
  let q : k := b * g.toBlocks₁₁ 2 0 + d * g.toBlocks₂₁ 2 0
  let r : k := b * g.toBlocks₁₂ 2 0 + d * g.toBlocks₂₂ 2 0
  let s : k := a * g.toBlocks₁₁ 2 1 + c * g.toBlocks₂₁ 2 1
  let u : k := b * g.toBlocks₁₁ 2 1 + d * g.toBlocks₂₁ 2 1
  let v : k := b * g.toBlocks₁₂ 2 1 + d * g.toBlocks₂₂ 2 1
  have hq0 :
      q = a * g.toBlocks₁₂ 2 0 + c * g.toBlocks₂₂ 2 0 := by
    simpa [a, b, c, d, q] using hParams.hqq
  have hu0 :
      u = a * g.toBlocks₁₂ 2 1 + c * g.toBlocks₂₂ 2 1 := by
    simpa [a, b, c, d, u] using hParams.huu
  have hAblock :
      g.toBlocks₁₁ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u := by
    simpa [a, b, c, d, p, q, s, u] using
      (Ablock_of_shapes_params (k := k) (g := g) hShapes hParams)
  have hBblock :
      g.toBlocks₁₂ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v := by
    simpa [a, b, c, d, q, r, u, v] using
      (Bblock_of_shapes_params (k := k) (g := g) hShapes hParams)
  have hCblock :
      g.toBlocks₂₁ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u := by
    simpa [a, b, c, d, p, q, s, u] using
      (Cblock_of_shapes_params (k := k) (g := g) hShapes hParams)
  have hDblock :
      g.toBlocks₂₂ =
        Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
          Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
            Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v := by
    simpa [a, b, c, d, q, r, u, v] using
      (Dblock_of_shapes_params (k := k) (g := g) hShapes hParams)
  refine (mem_K_shape_iff (k := k) (g := g)).2 ?_
  refine ⟨a, b, c, d, p, q, r, s, u, v, hParams.hΔ, ?_⟩
  calc
    g = Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂ := by
      exact (Matrix.fromBlocks_toBlocks g).symm
    _ =
        Matrix.fromBlocks
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
            Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u)
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
            Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v)
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledA (k := k) p s +
            Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledC (k := k) q u)
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledB (k := k) q u +
            Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d *
              Wedge2Formalization.N6DoublePureSingular.coupledD (k := k) r v) := by
          rw [hAblock, hBblock, hCblock, hDblock]
    _ =
        Matrix.fromBlocks
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseA (k := k) a b c d)
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseB (k := k) a b c d)
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseC (k := k) a b c d)
          (Wedge2Formalization.N6DoublePureSingular.magmaPointwiseD (k := k) a b c d) *
        Wedge2Formalization.N6DoublePureSingular.coupledUnipotent (k := k) p s q u r v := by
          simp [Wedge2Formalization.N6DoublePureSingular.coupledUnipotent,
            Matrix.fromBlocks_multiply]
    _ =
        Wedge2Formalization.N6DoublePureSingular.magmaKernelCell
          (k := k) a b c d p q r s u v := by
          rfl

private theorem mem_K_of_pointwise
    (g :
      Matrix Wedge2Formalization.N6DoublePureSingular.V
        Wedge2Formalization.N6DoublePureSingular.V k)
    (hfix : Wedge2Formalization.N6DoublePureSingular.FixesPairBivector g) :
    g ∈ K (k := k) := by
  have hShapes := shapes_of_pointwise (k := k) (g := g) hfix
  have hfix1_blocks :
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₁ (k := k))
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)
        =
      rep₁ (k := k) := by
    simpa [Matrix.fromBlocks_toBlocks] using hfix.1
  have hfix2_blocks :
      Wedge2Formalization.N6DoublePureSingular.ActBivector
          (rep₂ (k := k))
          (Matrix.fromBlocks g.toBlocks₁₁ g.toBlocks₁₂ g.toBlocks₂₁ g.toBlocks₂₂)
        =
      rep₂ (k := k) := by
    simpa [Matrix.fromBlocks_toBlocks] using hfix.2
  have hEqns :=
    blockEqns_of_shapes (k := k) (g := g) hShapes.Ashape hShapes.Bshape hShapes.Cshape
      hShapes.Dshape hfix1_blocks hfix2_blocks
  have hParams := paramEqns_of_blockEqns (k := k) (g := g) hEqns
  exact mem_K_of_shapes_params (k := k) (g := g) hShapes hParams

theorem pointwise_stabilizer :
    ExactFamily
      (Wedge2Formalization.N6DoublePureSingular.FixesPairBivector)
      (K (k := k)) := by
  intro g
  constructor
  · exact mem_K_of_pointwise (k := k) (g := g)
  · exact pointwise_of_mem_K (k := k) g

end Row11

end N6
end Paper
end Wedge2Formalization
