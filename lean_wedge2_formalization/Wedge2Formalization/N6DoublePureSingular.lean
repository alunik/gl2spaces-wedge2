import Wedge2Formalization.N3PureSingular
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.Tactic

open Matrix

namespace Wedge2Formalization
namespace N6DoublePureSingular

variable {k : Type*} [Field k]

/-- The first `3`-dimensional pure singular block. -/
abbrev W := N3PureSingular.I

/-- The second `3`-dimensional pure singular block. -/
abbrev I := N3PureSingular.I

/-- The ambient `6`-dimensional space. -/
abbrev V := W ⊕ I

/-- The natural action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The bivector action is multiplicative in the acting matrix. -/
theorem actBivector_mul
    (Ω : Matrix V V k)
    (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  simp [ActBivector, Matrix.mul_assoc]

/-- Exact stabilizer of a bivector. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  ActBivector Ω g = Ω

/-- The direct-sum pure singular representative. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks (N3PureSingular.ω13 (k := k)) 0 0 (N3PureSingular.ω13 (k := k))

/-- The second basis vector of the direct-sum pure singular representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks (N3PureSingular.ω23 (k := k)) 0 0 (N3PureSingular.ω23 (k := k))

/-- Exact stabilizer of the pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rep₁ g ∧ FixesBivector rep₂ g

/-- Setwise stabilizer of the `2`-space in the chosen basis. -/
def PreservesSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector rep₁ g = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ g = γ • rep₁ + δ • rep₂

/-- The `3 × 3` upper-left pointwise factor in the direct-sum `S_2^2` row. -/
def coupledA (a b : k) : Matrix W W k :=
  !![(1 : k), 0, 0;
      0, 1, 0;
      a, b, 1]

/-- The `3 × 3` upper-right coupled factor in the direct-sum `S_2^2` row. -/
def coupledB (z w : k) : Matrix W I k :=
  !![(0 : k), 0, 0;
      0, 0, 0;
      z, w, 0]

/-- The `3 × 3` lower-left coupled factor in the direct-sum `S_2^2` row. -/
def coupledC (z w : k) : Matrix I W k :=
  !![(0 : k), 0, 0;
      0, 0, 0;
      z, w, 0]

/-- The `3 × 3` lower-right pointwise factor in the direct-sum `S_2^2` row. -/
def coupledD (r s : k) : Matrix I I k :=
  !![(1 : k), 0, 0;
      0, 1, 0;
      r, s, 1]

/-- A concrete six-parameter unipotent family in the pointwise stabilizer of the
direct-sum pure singular pair. This is the `Sym_2(k) ⊕ Sym_2(k)` part of the kernel in
the basis used by the paper. -/
def coupledUnipotent (x y z w r s : k) : Matrix V V k :=
  Matrix.fromBlocks
    (coupledA (k := k) x y)
    (coupledB (k := k) z w)
    (coupledC (k := k) z w)
    (coupledD (k := k) r s)

/-- A `3 × 3` matrix that scales the first two basis vectors by `u` and the third by
`t`. This is the basic building block of the transported Magma pointwise `GL₂` family
on the direct-sum pure singular row. -/
def scalarTri (u t : k) : Matrix W W k :=
  !![u, 0, 0;
      0, u, 0;
      0, 0, t]

/-- Conjugating `ω13` by scalar-triangular matrices gives the expected transported
simple bivector with independent top-right and bottom-left coefficients. -/
theorem scalarTri_mul_ω13_mul_transpose
    (u t u' t' : k) :
    scalarTri (k := k) u t * N3PureSingular.ω13 * (scalarTri (k := k) u' t')ᵀ =
      !![(0 : k), 0, u * t';
          0, 0, 0;
          -(t * u'), 0, 0] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [scalarTri, N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three]

/-- Conjugating `ω23` by scalar-triangular matrices gives the expected transported
simple bivector with independent top-right and bottom-left coefficients. -/
theorem scalarTri_mul_ω23_mul_transpose
    (u t u' t' : k) :
    scalarTri (k := k) u t * N3PureSingular.ω23 * (scalarTri (k := k) u' t')ᵀ =
      !![(0 : k), 0, 0;
          0, 0, u * t';
          0, -(t * u'), 0] := by
  ext i j
  fin_cases i <;> fin_cases j <;>
    simp [scalarTri, N3PureSingular.ω23, Matrix.mul_apply, Matrix.transpose_apply,
      Fin.sum_univ_three]

/-- The `3 × 3` upper-left block of the pointwise `GL₂` family coming from the native
Magma `S_2^2` local kernel model, after transport to the current Lean basis. -/
def magmaPointwiseA (a b c d : k) : Matrix W W k :=
  scalarTri (k := k) a (d / (a * d - b * c))

/-- The `3 × 3` upper-right block of the transported Magma pointwise `GL₂` family. -/
def magmaPointwiseB (a b c d : k) : Matrix W I k :=
  scalarTri (k := k) b (-c / (a * d - b * c))

/-- The `3 × 3` lower-left block of the transported Magma pointwise `GL₂` family. -/
def magmaPointwiseC (a b c d : k) : Matrix I W k :=
  scalarTri (k := k) c (-b / (a * d - b * c))

/-- The `3 × 3` lower-right block of the transported Magma pointwise `GL₂` family. -/
def magmaPointwiseD (a b c d : k) : Matrix I I k :=
  scalarTri (k := k) d (a / (a * d - b * c))

/-- The transported Magma pointwise `GL₂` family inside the `S_2^2` local kernel. -/
def magmaPointwiseGL2 (a b c d : k) : Matrix V V k :=
  Matrix.fromBlocks
    (magmaPointwiseA (k := k) a b c d)
    (magmaPointwiseB (k := k) a b c d)
    (magmaPointwiseC (k := k) a b c d)
    (magmaPointwiseD (k := k) a b c d)

/-- The full transported Magma `S_2^2` local kernel cell is the pointwise `GL₂`
family followed by the existing coupled unipotent family. The parameter order matches
the original Magma symmetric matrices
`C₁ = [[p, q], [q, r]]`, `C₂ = [[s, u], [u, v]]`. -/
def magmaKernelCell (a b c d p q r s u v : k) : Matrix V V k :=
  magmaPointwiseGL2 (k := k) a b c d *
    coupledUnipotent (k := k) p s q u r v

-- Local copy of the direct-sum block-action formula, placed here so the Magma
-- pointwise theorem can use it before the general theorem is introduced below.
theorem act_diagPair_fromBlocks_local
    (Ω : Matrix W W k)
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    ActBivector (Matrix.fromBlocks Ω 0 0 Ω) (Matrix.fromBlocks A B C D) =
      Matrix.fromBlocks
        (A * Ω * Aᵀ + B * Ω * Bᵀ)
        (A * Ω * Cᵀ + B * Ω * Dᵀ)
        (C * Ω * Aᵀ + D * Ω * Bᵀ)
        (C * Ω * Cᵀ + D * Ω * Dᵀ) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
      Matrix.mul_assoc]

-- The transported Magma pointwise `GL₂` family fixes the direct-sum pure singular
-- pair pointwise. In the paired Magma basis this is the standard
-- `diag(T, T, T^{-T})` family.
set_option maxHeartbeats 1000000 in
-- The proof expands the transported Magma block formulas entrywise.
theorem magmaPointwiseGL2_family
    (a b c d : k)
    (hΔ : a * d - b * c ≠ 0) :
    FixesPairBivector (magmaPointwiseGL2 (k := k) a b c d) := by
  have hΔ' : d * a - b * c ≠ 0 := by
    simpa [mul_comm] using hΔ
  have hΔ'' : -(c * b) + a * d ≠ 0 := by
    simpa [sub_eq_add_neg, mul_comm, add_comm, add_left_comm, add_assoc] using hΔ
  have hΔ''' : -(b * c) + a * d ≠ 0 := by
    simpa [mul_comm] using hΔ''
  have hAA13 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseA (k := k) a b c d)ᵀ =
        !![(0 : k), 0, a * (d / (a * d - b * c));
            0, 0, 0;
            -((d / (a * d - b * c)) * a), 0, 0] := by
    simpa [magmaPointwiseA] using
      (scalarTri_mul_ω13_mul_transpose (k := k) a (d / (a * d - b * c))
        a (d / (a * d - b * c)))
  have hBB13 :
      magmaPointwiseB (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        !![(0 : k), 0, b * (-c / (a * d - b * c));
            0, 0, 0;
            -((-c / (a * d - b * c)) * b), 0, 0] := by
    simpa [magmaPointwiseB] using
      (scalarTri_mul_ω13_mul_transpose (k := k) b (-c / (a * d - b * c))
        b (-c / (a * d - b * c)))
  have hAC13 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseC (k := k) a b c d)ᵀ =
        !![(0 : k), 0, a * (-b / (a * d - b * c));
            0, 0, 0;
            -((d / (a * d - b * c)) * c), 0, 0] := by
    simpa [magmaPointwiseA, magmaPointwiseC] using
      (scalarTri_mul_ω13_mul_transpose (k := k) a (d / (a * d - b * c))
        c (-b / (a * d - b * c)))
  have hBD13 :
      magmaPointwiseB (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        !![(0 : k), 0, b * (a / (a * d - b * c));
            0, 0, 0;
            -((-c / (a * d - b * c)) * d), 0, 0] := by
    simpa [magmaPointwiseB, magmaPointwiseD] using
      (scalarTri_mul_ω13_mul_transpose (k := k) b (-c / (a * d - b * c))
        d (a / (a * d - b * c)))
  have hCA13 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseA (k := k) a b c d)ᵀ =
        !![(0 : k), 0, c * (d / (a * d - b * c));
            0, 0, 0;
            -((-b / (a * d - b * c)) * a), 0, 0] := by
    simpa [magmaPointwiseA, magmaPointwiseC] using
      (scalarTri_mul_ω13_mul_transpose (k := k) c (-b / (a * d - b * c))
        a (d / (a * d - b * c)))
  have hDB13 :
      magmaPointwiseD (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        !![(0 : k), 0, d * (-c / (a * d - b * c));
            0, 0, 0;
            -((a / (a * d - b * c)) * b), 0, 0] := by
    simpa [magmaPointwiseB, magmaPointwiseD] using
      (scalarTri_mul_ω13_mul_transpose (k := k) d (a / (a * d - b * c))
        b (-c / (a * d - b * c)))
  have hCC13 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseC (k := k) a b c d)ᵀ =
        !![(0 : k), 0, c * (-b / (a * d - b * c));
            0, 0, 0;
            -((-b / (a * d - b * c)) * c), 0, 0] := by
    simpa [magmaPointwiseC] using
      (scalarTri_mul_ω13_mul_transpose (k := k) c (-b / (a * d - b * c))
        c (-b / (a * d - b * c)))
  have hDD13 :
      magmaPointwiseD (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        !![(0 : k), 0, d * (a / (a * d - b * c));
            0, 0, 0;
            -((a / (a * d - b * c)) * d), 0, 0] := by
    simpa [magmaPointwiseD] using
      (scalarTri_mul_ω13_mul_transpose (k := k) d (a / (a * d - b * c))
        d (a / (a * d - b * c)))
  have hAA23 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseA (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, a * (d / (a * d - b * c));
            0, -((d / (a * d - b * c)) * a), 0] := by
    simpa [magmaPointwiseA] using
      (scalarTri_mul_ω23_mul_transpose (k := k) a (d / (a * d - b * c))
        a (d / (a * d - b * c)))
  have hBB23 :
      magmaPointwiseB (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, b * (-c / (a * d - b * c));
            0, -((-c / (a * d - b * c)) * b), 0] := by
    simpa [magmaPointwiseB] using
      (scalarTri_mul_ω23_mul_transpose (k := k) b (-c / (a * d - b * c))
        b (-c / (a * d - b * c)))
  have hAC23 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseC (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, a * (-b / (a * d - b * c));
            0, -((d / (a * d - b * c)) * c), 0] := by
    simpa [magmaPointwiseA, magmaPointwiseC] using
      (scalarTri_mul_ω23_mul_transpose (k := k) a (d / (a * d - b * c))
        c (-b / (a * d - b * c)))
  have hBD23 :
      magmaPointwiseB (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, b * (a / (a * d - b * c));
            0, -((-c / (a * d - b * c)) * d), 0] := by
    simpa [magmaPointwiseB, magmaPointwiseD] using
      (scalarTri_mul_ω23_mul_transpose (k := k) b (-c / (a * d - b * c))
        d (a / (a * d - b * c)))
  have hCA23 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseA (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, c * (d / (a * d - b * c));
            0, -((-b / (a * d - b * c)) * a), 0] := by
    simpa [magmaPointwiseA, magmaPointwiseC] using
      (scalarTri_mul_ω23_mul_transpose (k := k) c (-b / (a * d - b * c))
        a (d / (a * d - b * c)))
  have hDB23 :
      magmaPointwiseD (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, d * (-c / (a * d - b * c));
            0, -((a / (a * d - b * c)) * b), 0] := by
    simpa [magmaPointwiseB, magmaPointwiseD] using
      (scalarTri_mul_ω23_mul_transpose (k := k) d (a / (a * d - b * c))
        b (-c / (a * d - b * c)))
  have hCC23 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseC (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, c * (-b / (a * d - b * c));
            0, -((-b / (a * d - b * c)) * c), 0] := by
    simpa [magmaPointwiseC] using
      (scalarTri_mul_ω23_mul_transpose (k := k) c (-b / (a * d - b * c))
        c (-b / (a * d - b * c)))
  have hDD23 :
      magmaPointwiseD (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        !![(0 : k), 0, 0;
            0, 0, d * (a / (a * d - b * c));
            0, -((a / (a * d - b * c)) * d), 0] := by
    simpa [magmaPointwiseD] using
      (scalarTri_mul_ω23_mul_transpose (k := k) d (a / (a * d - b * c))
        d (a / (a * d - b * c)))
  have h11_13 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseA (k := k) a b c d)ᵀ +
        magmaPointwiseB (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        N3PureSingular.ω13 (k := k) := by
    rw [hAA13, hBB13]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [N3PureSingular.ω13, Matrix.add_apply] <;>
      field_simp [hΔ, hΔ', hΔ''] <;>
      ring_nf
  have h12_13 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseC (k := k) a b c d)ᵀ +
        magmaPointwiseB (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        0 := by
    rw [hAC13, hBD13]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.add_apply] <;>
      field_simp [hΔ, hΔ', hΔ''] <;>
      ring_nf
  have h21_13 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseA (k := k) a b c d)ᵀ +
        magmaPointwiseD (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        0 := by
    rw [hCA13, hDB13]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.add_apply] <;>
      field_simp [hΔ, hΔ', hΔ''] <;>
      ring_nf
  have h22_13 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseC (k := k) a b c d)ᵀ +
        magmaPointwiseD (k := k) a b c d * N3PureSingular.ω13 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        N3PureSingular.ω13 (k := k) := by
    rw [hCC13, hDD13]
    ext i j
    fin_cases i <;> fin_cases j
    · simp [N3PureSingular.ω13, Matrix.add_apply]
    · simp [N3PureSingular.ω13, Matrix.add_apply]
    ·
      have hscalar :
          -(c * (b * (a * d - b * c)⁻¹)) + d * (a * (a * d - b * c)⁻¹) = 1 := by
        calc
          -(c * (b * (a * d - b * c)⁻¹)) + d * (a * (a * d - b * c)⁻¹) =
              (a * d - b * c) * (a * d - b * c)⁻¹ := by ring
          _ = 1 := by exact mul_inv_cancel₀ hΔ
      simpa [N3PureSingular.ω13, Matrix.add_apply, div_eq_mul_inv] using hscalar
    · simp [N3PureSingular.ω13, Matrix.add_apply]
    · simp [N3PureSingular.ω13, Matrix.add_apply]
    · simp [N3PureSingular.ω13, Matrix.add_apply]
    ·
      simp [N3PureSingular.ω13, Matrix.add_apply]
      field_simp [hΔ, hΔ', hΔ''']
      ring_nf
    · simp [N3PureSingular.ω13, Matrix.add_apply]
    · simp [N3PureSingular.ω13, Matrix.add_apply]
  have h11_23 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseA (k := k) a b c d)ᵀ +
        magmaPointwiseB (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        N3PureSingular.ω23 (k := k) := by
    rw [hAA23, hBB23]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [N3PureSingular.ω23, Matrix.add_apply] <;>
      field_simp [hΔ, hΔ', hΔ''] <;>
      ring_nf
  have h12_23 :
      magmaPointwiseA (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseC (k := k) a b c d)ᵀ +
        magmaPointwiseB (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        0 := by
    rw [hAC23, hBD23]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.add_apply] <;>
      field_simp [hΔ, hΔ', hΔ''] <;>
      ring_nf
  have h21_23 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseA (k := k) a b c d)ᵀ +
        magmaPointwiseD (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseB (k := k) a b c d)ᵀ =
        0 := by
    rw [hCA23, hDB23]
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Matrix.add_apply] <;>
      field_simp [hΔ, hΔ', hΔ''] <;>
      ring_nf
  have h22_23 :
      magmaPointwiseC (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseC (k := k) a b c d)ᵀ +
        magmaPointwiseD (k := k) a b c d * N3PureSingular.ω23 *
          (magmaPointwiseD (k := k) a b c d)ᵀ =
        N3PureSingular.ω23 (k := k) := by
    rw [hCC23, hDD23]
    ext i j
    fin_cases i <;> fin_cases j
    · simp [N3PureSingular.ω23, Matrix.add_apply]
    · simp [N3PureSingular.ω23, Matrix.add_apply]
    · simp [N3PureSingular.ω23, Matrix.add_apply]
    · simp [N3PureSingular.ω23, Matrix.add_apply]
    · simp [N3PureSingular.ω23, Matrix.add_apply]
    ·
      have hscalar :
          -(c * (b * (a * d - b * c)⁻¹)) + d * (a * (a * d - b * c)⁻¹) = 1 := by
        calc
          -(c * (b * (a * d - b * c)⁻¹)) + d * (a * (a * d - b * c)⁻¹) =
              (a * d - b * c) * (a * d - b * c)⁻¹ := by ring
          _ = 1 := by exact mul_inv_cancel₀ hΔ
      simpa [N3PureSingular.ω23, Matrix.add_apply, div_eq_mul_inv] using hscalar
    · simp [N3PureSingular.ω23, Matrix.add_apply]
    ·
      simp [N3PureSingular.ω23, Matrix.add_apply]
      field_simp [hΔ, hΔ', hΔ''']
      ring_nf
    · simp [N3PureSingular.ω23, Matrix.add_apply]
  constructor
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₁ (magmaPointwiseGL2 (k := k) a b c d) =
          Matrix.fromBlocks
            (magmaPointwiseA (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseA (k := k) a b c d)ᵀ +
              magmaPointwiseB (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseB (k := k) a b c d)ᵀ)
            (magmaPointwiseA (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseC (k := k) a b c d)ᵀ +
              magmaPointwiseB (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseD (k := k) a b c d)ᵀ)
            (magmaPointwiseC (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseA (k := k) a b c d)ᵀ +
              magmaPointwiseD (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseB (k := k) a b c d)ᵀ)
            (magmaPointwiseC (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseC (k := k) a b c d)ᵀ +
              magmaPointwiseD (k := k) a b c d * N3PureSingular.ω13 *
                (magmaPointwiseD (k := k) a b c d)ᵀ) := by
            simpa [rep₁, magmaPointwiseGL2] using
              (act_diagPair_fromBlocks_local (k := k) (Ω := N3PureSingular.ω13)
                (A := magmaPointwiseA (k := k) a b c d)
                (B := magmaPointwiseB (k := k) a b c d)
                (C := magmaPointwiseC (k := k) a b c d)
                (D := magmaPointwiseD (k := k) a b c d))
      _ = rep₁ := by simp [rep₁, h11_13, h12_13, h21_13, h22_13]
  ·
    rw [FixesBivector]
    calc
      ActBivector rep₂ (magmaPointwiseGL2 (k := k) a b c d) =
          Matrix.fromBlocks
            (magmaPointwiseA (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseA (k := k) a b c d)ᵀ +
              magmaPointwiseB (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseB (k := k) a b c d)ᵀ)
            (magmaPointwiseA (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseC (k := k) a b c d)ᵀ +
              magmaPointwiseB (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseD (k := k) a b c d)ᵀ)
            (magmaPointwiseC (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseA (k := k) a b c d)ᵀ +
              magmaPointwiseD (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseB (k := k) a b c d)ᵀ)
            (magmaPointwiseC (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseC (k := k) a b c d)ᵀ +
              magmaPointwiseD (k := k) a b c d * N3PureSingular.ω23 *
                (magmaPointwiseD (k := k) a b c d)ᵀ) := by
            simpa [rep₂, magmaPointwiseGL2] using
              (act_diagPair_fromBlocks_local (k := k) (Ω := N3PureSingular.ω23)
                (A := magmaPointwiseA (k := k) a b c d)
                (B := magmaPointwiseB (k := k) a b c d)
                (C := magmaPointwiseC (k := k) a b c d)
                (D := magmaPointwiseD (k := k) a b c d))
      _ = rep₂ := by simp [rep₂, h11_23, h12_23, h21_23, h22_23]

/-- Acting on a direct-sum bivector `diag(Ω, Ω)` with a `2 × 2` block matrix gives the
expected block formula. -/
theorem act_diagPair_fromBlocks
    (Ω : Matrix W W k)
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    ActBivector (Matrix.fromBlocks Ω 0 0 Ω) (Matrix.fromBlocks A B C D) =
      Matrix.fromBlocks
        (A * Ω * Aᵀ + B * Ω * Bᵀ)
        (A * Ω * Cᵀ + B * Ω * Dᵀ)
        (C * Ω * Aᵀ + D * Ω * Bᵀ)
        (C * Ω * Cᵀ + D * Ω * Dᵀ) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
      Matrix.mul_add, Matrix.add_mul, Matrix.mul_assoc]

/-- Acting by a block diagonal matrix respects the direct-sum decomposition. -/
theorem act_blockDiagonal
    (G : Matrix W W k)
    (H : Matrix I I k) :
    ActBivector rep₁ (Matrix.fromBlocks G 0 0 H) =
      Matrix.fromBlocks
        (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G)
        0
        0
        (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) H) ∧
    ActBivector rep₂ (Matrix.fromBlocks G 0 0 H) =
      Matrix.fromBlocks
        (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G)
        0
        0
        (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) H) := by
  constructor <;>
  ·
    ext i j
    cases i <;> cases j <;>
      simp [ActBivector, rep₁, rep₂, N3PureSingular.ActBivector,
        Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

/-- For block diagonal matrices, fixing the direct-sum pure singular pair is equivalent
to fixing the pure singular pair on each summand. -/
theorem fixesPair_blockDiagonal_iff
    (G : Matrix W W k)
    (H : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks G 0 0 H) ↔
      N3PureSingular.FixesPairBivector G ∧ N3PureSingular.FixesPairBivector H := by
  constructor
  ·
    rintro ⟨h₁, h₂⟩
    have hact₁ := (act_blockDiagonal (k := k) G H).1
    have hact₂ := (act_blockDiagonal (k := k) G H).2
    have h₁' :
        Matrix.fromBlocks
          (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G)
          0
          0
          (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) H) =
        Matrix.fromBlocks
          (N3PureSingular.ω13 (k := k))
          0
          0
          (N3PureSingular.ω13 (k := k)) := by
      calc
        Matrix.fromBlocks
            (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G)
            0
            0
            (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) H) =
          ActBivector rep₁ (Matrix.fromBlocks G 0 0 H) := by
            symm
            exact hact₁
        _ =
          Matrix.fromBlocks
            (N3PureSingular.ω13 (k := k))
            0
            0
            (N3PureSingular.ω13 (k := k)) := by
              simpa [FixesBivector, rep₁] using h₁
    have h₂' :
        Matrix.fromBlocks
          (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G)
          0
          0
          (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) H) =
        Matrix.fromBlocks
          (N3PureSingular.ω23 (k := k))
          0
          0
          (N3PureSingular.ω23 (k := k)) := by
      calc
        Matrix.fromBlocks
            (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G)
            0
            0
            (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) H) =
          ActBivector rep₂ (Matrix.fromBlocks G 0 0 H) := by
            symm
            exact hact₂
        _ =
          Matrix.fromBlocks
            (N3PureSingular.ω23 (k := k))
            0
            0
            (N3PureSingular.ω23 (k := k)) := by
              simpa [FixesBivector, rep₂] using h₂
    constructor
    ·
      constructor
      ·
        rw [N3PureSingular.FixesBivector]
        have h11 := congrArg Matrix.toBlocks₁₁ h₁'
        simpa using h11
      ·
        rw [N3PureSingular.FixesBivector]
        have h11 := congrArg Matrix.toBlocks₁₁ h₂'
        simpa using h11
    ·
      constructor
      ·
        rw [N3PureSingular.FixesBivector]
        have h22 := congrArg Matrix.toBlocks₂₂ h₁'
        simpa using h22
      ·
        rw [N3PureSingular.FixesBivector]
        have h22 := congrArg Matrix.toBlocks₂₂ h₂'
        simpa using h22
  ·
    rintro ⟨hG, hH⟩
    constructor
    ·
      rw [FixesBivector]
      calc
        ActBivector rep₁ (Matrix.fromBlocks G 0 0 H) =
          Matrix.fromBlocks
            (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G)
            0
            0
            (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) H) :=
              (act_blockDiagonal (k := k) G H).1
        _ = rep₁ := by
              simpa [N3PureSingular.FixesBivector, rep₁] using And.intro hG.1 hH.1
    ·
      rw [FixesBivector]
      calc
        ActBivector rep₂ (Matrix.fromBlocks G 0 0 H) =
          Matrix.fromBlocks
            (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G)
            0
            0
            (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) H) :=
              (act_blockDiagonal (k := k) G H).2
        _ = rep₂ := by
              simpa [N3PureSingular.FixesBivector, rep₂] using And.intro hG.2 hH.2

/-- The block diagonal pointwise stabilizers have exactly the expected pair of pure
singular shapes on the two summands. -/
theorem fixesPair_blockDiagonal_iff_shape
    (G : Matrix W W k)
    (H : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks G 0 0 H) ↔
      G 0 0 ≠ 0 ∧ H 0 0 ≠ 0 ∧
        G = N3PureSingular.pureSingularShape (k := k) (G 0 0) (G 2 0) (G 2 1) ∧
        H = N3PureSingular.pureSingularShape (k := k) (H 0 0) (H 2 0) (H 2 1) := by
  rw [fixesPair_blockDiagonal_iff]
  constructor
  ·
    rintro ⟨hG, hH⟩
    rcases (N3PureSingular.fixesPairBivector_iff_shape (k := k) G).1 hG with ⟨hG00, hGshape⟩
    rcases (N3PureSingular.fixesPairBivector_iff_shape (k := k) H).1 hH with ⟨hH00, hHshape⟩
    exact ⟨hG00, hH00, hGshape, hHshape⟩
  ·
    rintro ⟨hG00, hH00, hGshape, hHshape⟩
    constructor
    ·
      exact (N3PureSingular.fixesPairBivector_iff_shape (k := k) G).2 ⟨hG00, hGshape⟩
    ·
      exact (N3PureSingular.fixesPairBivector_iff_shape (k := k) H).2 ⟨hH00, hHshape⟩

/-- Putting the full pure singular pointwise shape independently on the two summands
gives an explicit block-diagonal pointwise family on the direct-sum `S_2^2` row. -/
theorem pointwise_shape_family
    (a b x y r s : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b r s)) := by
  refine (fixesPair_blockDiagonal_iff_shape (k := k)
    (G := N3PureSingular.pureSingularShape (k := k) a x y)
    (H := N3PureSingular.pureSingularShape (k := k) b r s)).2 ?_
  exact ⟨ha, hb, rfl, rfl⟩

/-- Independent pointwise scalings on the two pure singular summands fix the pair. -/
theorem pointwise_scale_family
    (a b : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        (N3PureSingular.pointwiseScale (k := k) b)) := by
  constructor
  ·
    have hact := (act_blockDiagonal (k := k)
      (N3PureSingular.pointwiseScale (k := k) a)
      (N3PureSingular.pointwiseScale (k := k) b)).1
    have hleft :
        N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k))
          (N3PureSingular.pointwiseScale (k := k) a) =
          N3PureSingular.ω13 (k := k) := by
      exact (N3PureSingular.pointwise_scale_family (k := k) a ha).1
    have hright :
        N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k))
          (N3PureSingular.pointwiseScale (k := k) b) =
          N3PureSingular.ω13 (k := k) := by
      exact (N3PureSingular.pointwise_scale_family (k := k) b hb).1
    rw [FixesBivector]
    calc
      ActBivector rep₁
          (Matrix.fromBlocks
            (N3PureSingular.pointwiseScale (k := k) a)
            0
            0
            (N3PureSingular.pointwiseScale (k := k) b)) =
        Matrix.fromBlocks
          (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k))
            (N3PureSingular.pointwiseScale (k := k) a))
          0
          0
          (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k))
            (N3PureSingular.pointwiseScale (k := k) b)) := hact
      _ = rep₁ := by simpa [rep₁, hleft, hright]
  ·
    have hact := (act_blockDiagonal (k := k)
      (N3PureSingular.pointwiseScale (k := k) a)
      (N3PureSingular.pointwiseScale (k := k) b)).2
    have hleft :
        N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k))
          (N3PureSingular.pointwiseScale (k := k) a) =
          N3PureSingular.ω23 (k := k) := by
      exact (N3PureSingular.pointwise_scale_family (k := k) a ha).2
    have hright :
        N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k))
          (N3PureSingular.pointwiseScale (k := k) b) =
          N3PureSingular.ω23 (k := k) := by
      exact (N3PureSingular.pointwise_scale_family (k := k) b hb).2
    rw [FixesBivector]
    calc
      ActBivector rep₂
          (Matrix.fromBlocks
            (N3PureSingular.pointwiseScale (k := k) a)
            0
            0
            (N3PureSingular.pointwiseScale (k := k) b)) =
        Matrix.fromBlocks
          (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k))
            (N3PureSingular.pointwiseScale (k := k) a))
          0
          0
          (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k))
            (N3PureSingular.pointwiseScale (k := k) b)) := hact
      _ = rep₂ := by simpa [rep₂, hleft, hright]

-- The six-parameter coupled unipotent family fixes the direct-sum pure singular pair
-- pointwise. This proof is a large entrywise verification, so it needs a larger
-- heartbeat budget than the default.
set_option maxHeartbeats 1000000 in
-- The proof is a direct entrywise check of the transported Magma kernel coordinates.
theorem coupled_unipotent_family
    (x y z w r s : k) :
    FixesPairBivector (coupledUnipotent (k := k) x y z w r s) := by
  constructor <;> rw [FixesBivector]
  ·
    ext i j
    rcases i with i | i <;> rcases j with j | j
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₁, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω13, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₁, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω13, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₁, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω13, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₁, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω13, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]
  ·
    ext i j
    rcases i with i | i <;> rcases j with j | j
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₂, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω23, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₂, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω23, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₂, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω23, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]
    ·
      fin_cases i <;> fin_cases j <;>
        simp [ActBivector, rep₂, coupledUnipotent, coupledA, coupledB, coupledC, coupledD,
          N3PureSingular.ω23, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_three]

/-- The full transported Magma `S_2^2` local kernel cell fixes the direct-sum pure
singular pair pointwise. -/
theorem magmaKernelCell_pointwise
    (a b c d p q r s u v : k)
    (hΔ : a * d - b * c ≠ 0) :
    FixesPairBivector (magmaKernelCell (k := k) a b c d p q r s u v) := by
  have hG := magmaPointwiseGL2_family (k := k) a b c d hΔ
  have hU := coupled_unipotent_family (k := k) p s q u r v
  constructor
  ·
    calc
      ActBivector rep₁ (magmaKernelCell (k := k) a b c d p q r s u v) =
          ActBivector
            (ActBivector rep₁ (coupledUnipotent (k := k) p s q u r v))
            (magmaPointwiseGL2 (k := k) a b c d) := by
              rw [magmaKernelCell, actBivector_mul]
      _ = ActBivector rep₁ (magmaPointwiseGL2 (k := k) a b c d) := by
            simpa [FixesBivector] using congrArg (fun M =>
              ActBivector M (magmaPointwiseGL2 (k := k) a b c d)) hU.1
      _ = rep₁ := hG.1
  ·
    calc
      ActBivector rep₂ (magmaKernelCell (k := k) a b c d p q r s u v) =
          ActBivector
            (ActBivector rep₂ (coupledUnipotent (k := k) p s q u r v))
            (magmaPointwiseGL2 (k := k) a b c d) := by
              rw [magmaKernelCell, actBivector_mul]
      _ = ActBivector rep₂ (magmaPointwiseGL2 (k := k) a b c d) := by
            simpa [FixesBivector] using congrArg (fun M =>
              ActBivector M (magmaPointwiseGL2 (k := k) a b c d)) hU.2
      _ = rep₂ := hG.2

/-- The coupled unipotent family combines with the independent pointwise scalings on
the two pure singular summands to give a larger concrete pointwise subgroup on the
direct-sum `S_2^2` orbit. -/
theorem coupled_scale_product_pointwise
    (x y z w r s a b : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gU : Matrix V V k := coupledUnipotent (k := k) x y z w r s
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        (N3PureSingular.pointwiseScale (k := k) b)
    FixesPairBivector (gU * gS) := by
  intro gU gS
  have hU := coupled_unipotent_family (k := k) x y z w r s
  have hS := pointwise_scale_family (k := k) a b ha hb
  constructor
  ·
    calc
      ActBivector rep₁ (gU * gS) = ActBivector (ActBivector rep₁ gS) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gU := by
        simpa [gS] using congrArg (fun M => ActBivector M gU) hS.1
      _ = rep₁ := hU.1
  ·
    calc
      ActBivector rep₂ (gU * gS) = ActBivector (ActBivector rep₂ gS) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by
        simpa [gS] using congrArg (fun M => ActBivector M gU) hS.2
      _ = rep₂ := hU.2

/-- The coupled unipotent family also combines with the full pointwise pure singular
shape on each summand, not only with the scalar torus part. -/
theorem coupled_shape_product_pointwise
    (u v z w r s a b x y p q : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gU : Matrix V V k := coupledUnipotent (k := k) u v z w r s
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    FixesPairBivector (gU * gS) := by
  intro gU gS
  have hU := coupled_unipotent_family (k := k) u v z w r s
  have hS := pointwise_shape_family (k := k) a b x y p q ha hb
  constructor
  ·
    calc
      ActBivector rep₁ (gU * gS) = ActBivector (ActBivector rep₁ gS) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gU := by
        simpa [gS] using congrArg (fun M => ActBivector M gU) hS.1
      _ = rep₁ := hU.1
  ·
    calc
      ActBivector rep₂ (gU * gS) = ActBivector (ActBivector rep₂ gS) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by
        simpa [gS] using congrArg (fun M => ActBivector M gU) hS.2
      _ = rep₂ := hU.2

/-- Right-multiplying a block-diagonal pointwise shape by the coupled unipotent family
does not change whether the direct-sum pure singular pair is fixed. -/
theorem blockDiagonal_coupled_right_product_iff_shape
    (x y z w r s : k)
    (G : Matrix W W k)
    (H : Matrix I I k) :
    FixesPairBivector
        ((Matrix.fromBlocks G 0 0 H) * coupledUnipotent (k := k) x y z w r s) ↔
      G 0 0 ≠ 0 ∧ H 0 0 ≠ 0 ∧
        G = N3PureSingular.pureSingularShape (k := k) (G 0 0) (G 2 0) (G 2 1) ∧
        H = N3PureSingular.pureSingularShape (k := k) (H 0 0) (H 2 0) (H 2 1) := by
  have hU := coupled_unipotent_family (k := k) x y z w r s
  constructor
  · rintro ⟨h₁, h₂⟩
    have h₁' : ActBivector rep₁ (Matrix.fromBlocks G 0 0 H) = rep₁ := by
      have htmp : ActBivector rep₁
          ((Matrix.fromBlocks G 0 0 H) * coupledUnipotent (k := k) x y z w r s) = rep₁ := by
        simpa [FixesBivector] using h₁
      rw [actBivector_mul] at htmp
      rw [hU.1] at htmp
      exact htmp
    have h₂' : ActBivector rep₂ (Matrix.fromBlocks G 0 0 H) = rep₂ := by
      have htmp : ActBivector rep₂
          ((Matrix.fromBlocks G 0 0 H) * coupledUnipotent (k := k) x y z w r s) = rep₂ := by
        simpa [FixesBivector] using h₂
      rw [actBivector_mul] at htmp
      rw [hU.2] at htmp
      exact htmp
    exact
      (fixesPair_blockDiagonal_iff_shape (k := k) (G := G) (H := H)).1
        ⟨h₁', h₂'⟩
  · intro hshape
    have hS :
        FixesPairBivector (Matrix.fromBlocks G 0 0 H) :=
      (fixesPair_blockDiagonal_iff_shape (k := k) (G := G) (H := H)).2 hshape
    constructor
    · calc
        ActBivector rep₁
            ((Matrix.fromBlocks G 0 0 H) * coupledUnipotent (k := k) x y z w r s) =
          ActBivector
            (ActBivector rep₁ (coupledUnipotent (k := k) x y z w r s))
            (Matrix.fromBlocks G 0 0 H) := by
              rw [actBivector_mul]
        _ = ActBivector rep₁ (Matrix.fromBlocks G 0 0 H) := by
              rw [hU.1]
        _ = rep₁ := hS.1
    · calc
        ActBivector rep₂
            ((Matrix.fromBlocks G 0 0 H) * coupledUnipotent (k := k) x y z w r s) =
          ActBivector
            (ActBivector rep₂ (coupledUnipotent (k := k) x y z w r s))
            (Matrix.fromBlocks G 0 0 H) := by
              rw [actBivector_mul]
        _ = ActBivector rep₂ (Matrix.fromBlocks G 0 0 H) := by
              rw [hU.2]
        _ = rep₂ := hS.2

/-- Applying the same `GL₂` lift on both pure singular summands preserves the chosen
basis with the expected coefficients. -/
theorem GL2_lift_action
    (α β γ δ : k) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix V V k := Matrix.fromBlocks G 0 0 G
    ActBivector rep₁ g = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ g = γ • rep₁ + δ • rep₂ := by
  intro G g
  have hact := act_blockDiagonal (k := k) G G
  have hG := N3PureSingular.GL2_lift_action (k := k) α β γ δ
  constructor
  ·
    calc
      ActBivector rep₁ g =
        Matrix.fromBlocks
          (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G)
          0
          0
          (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G) := hact.1
      _ = Matrix.fromBlocks
          (α • N3PureSingular.ω13 + β • N3PureSingular.ω23)
          0
          0
          (α • N3PureSingular.ω13 + β • N3PureSingular.ω23) := by
            simpa [G] using hG.1
      _ = α • rep₁ + β • rep₂ := by
        ext i j
        cases i <;> cases j <;>
          simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]
  ·
    calc
      ActBivector rep₂ g =
        Matrix.fromBlocks
          (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G)
          0
          0
          (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G) := hact.2
      _ = Matrix.fromBlocks
          (γ • N3PureSingular.ω13 + δ • N3PureSingular.ω23)
          0
          0
          (γ • N3PureSingular.ω13 + δ • N3PureSingular.ω23) := by
            simpa [G] using hG.2
      _ = γ • rep₁ + δ • rep₂ := by
        ext i j
        cases i <;> cases j <;>
          simp [rep₁, rep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The full block-diagonal pure singular pointwise shape has determinant `ab`. -/
theorem pointwise_shape_det
    (a b x y p q : k) :
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    Matrix.det gS = a * b := by
  intro gS
  rw [Matrix.det_fromBlocks_zero₂₁,
    N3PureSingular.pureSingularShape_det,
    N3PureSingular.pureSingularShape_det]

/-- The simultaneous `GL₂` lift on the two pure singular summands has determinant
`(αδ-βγ)^2`. -/
theorem GL2_lift_det
    (α β γ δ : k) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let g : Matrix V V k := Matrix.fromBlocks G 0 0 G
    Matrix.det g = (α * δ - β * γ) ^ 2 := by
  intro G g
  rw [Matrix.det_fromBlocks_zero₂₁]
  simpa [G, N3PureSingular.pureLift_det, pow_two]

/-- The coupled unipotent family and the simultaneous `GL₂` lift combine exactly as
expected for the semidirect product description of the direct-sum `S_2^2` row. -/
theorem coupled_GL2_product_lift_action
    (x y z w r s α β γ δ : k) :
    let gU : Matrix V V k := coupledUnipotent (k := k) x y z w r s
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    ActBivector rep₁ (gU * gL) = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ (gU * gL) = γ • rep₁ + δ • rep₂ := by
  intro gU G gL
  have hU := coupled_unipotent_family (k := k) x y z w r s
  have hL := GL2_lift_action (k := k) α β γ δ
  constructor
  ·
    calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (α • rep₁ + β • rep₂) gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.1
      _ = α • ActBivector rep₁ gU + β • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = α • rep₁ + β • rep₂ := by
        rw [hU.1, hU.2]
  ·
    calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (γ • rep₁ + δ • rep₂) gU := by
        simpa [gL] using congrArg (fun M => ActBivector M gU) hL.2
      _ = γ • ActBivector rep₁ gU + δ • ActBivector rep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = γ • rep₁ + δ • rep₂ := by
        rw [hU.1, hU.2]

/-- The full block-diagonal pure singular pointwise shape combines with the simultaneous
`GL₂` lift exactly as expected. -/
theorem shape_GL2_product_lift_action
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    ActBivector rep₁ (gS * gL) = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ (gS * gL) = γ • rep₁ + δ • rep₂ := by
  intro gS G gL
  have hS := pointwise_shape_family (k := k) a b x y p q ha hb
  have hL := GL2_lift_action (k := k) α β γ δ
  constructor
  · calc
      ActBivector rep₁ (gS * gL) = ActBivector (ActBivector rep₁ gL) gS := by
        rw [actBivector_mul]
      _ = ActBivector (α • rep₁ + β • rep₂) gS := by
        simpa [gL] using congrArg (fun M => ActBivector M gS) hL.1
      _ = α • ActBivector rep₁ gS + β • ActBivector rep₂ gS := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = α • rep₁ + β • rep₂ := by
        rw [hS.1, hS.2]
  · calc
      ActBivector rep₂ (gS * gL) = ActBivector (ActBivector rep₂ gL) gS := by
        rw [actBivector_mul]
      _ = ActBivector (γ • rep₁ + δ • rep₂) gS := by
        simpa [gL] using congrArg (fun M => ActBivector M gS) hL.2
      _ = γ • ActBivector rep₁ gS + δ • ActBivector rep₂ gS := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = γ • rep₁ + δ • rep₂ := by
        rw [hS.1, hS.2]

/-- Left-multiplying the simultaneous `GL₂` lift by the full block-diagonal pure
singular pointwise shape changes the determinant only by the pointwise factor `ab`. -/
theorem shape_GL2_product_lift_det
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    Matrix.det (gS * gL) = (a * b) * (α * δ - β * γ) ^ 2 := by
  intro gS G gL
  rw [Matrix.det_mul]
  rw [show Matrix.det gS = a * b by
      simpa [gS] using pointwise_shape_det (k := k) a b x y p q]
  rw [show Matrix.det gL = (α * δ - β * γ) ^ 2 by
      simpa [gL] using GL2_lift_det (k := k) α β γ δ]


/-- Right-multiplying the simultaneous `GL₂` lift by the coupled unipotent family does
not change the quotient action on the chosen basis. -/
theorem GL2_coupled_right_product_lift_action
    (x y z w r s α β γ δ : k) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    let gU : Matrix V V k := coupledUnipotent (k := k) x y z w r s
    ActBivector rep₁ (gL * gU) = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ (gL * gU) = γ • rep₁ + δ • rep₂ := by
  intro G gL gU
  have hL := GL2_lift_action (k := k) α β γ δ
  have hU := coupled_unipotent_family (k := k) x y z w r s
  constructor
  · calc
      ActBivector rep₁ (gL * gU) = ActBivector (ActBivector rep₁ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.1
      _ = α • rep₁ + β • rep₂ := by
        simpa [gL] using hL.1
  · calc
      ActBivector rep₂ (gL * gU) = ActBivector (ActBivector rep₂ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hU.2
      _ = γ • rep₁ + δ • rep₂ := by
        simpa [gL] using hL.2

/-- Right-multiplying the simultaneous `GL₂` lift by the full block-diagonal pure
singular pointwise shape does not change the quotient action. -/
theorem GL2_shape_right_product_lift_action
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    ActBivector rep₁ (gL * gS) = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ (gL * gS) = γ • rep₁ + δ • rep₂ := by
  intro G gL gS
  have hL := GL2_lift_action (k := k) α β γ δ
  have hS := pointwise_shape_family (k := k) a b x y p q ha hb
  constructor
  · calc
      ActBivector rep₁ (gL * gS) = ActBivector (ActBivector rep₁ gS) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gS] using congrArg (fun M => ActBivector M gL) hS.1
      _ = α • rep₁ + β • rep₂ := by
        simpa [gL] using hL.1
  · calc
      ActBivector rep₂ (gL * gS) = ActBivector (ActBivector rep₂ gS) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gS] using congrArg (fun M => ActBivector M gL) hS.2
      _ = γ • rep₁ + δ • rep₂ := by
        simpa [gL] using hL.2

/-- Right-multiplying the simultaneous `GL₂` lift by the full block-diagonal pure
singular pointwise shape changes the determinant only by the pointwise factor `ab`. -/
theorem GL2_shape_right_product_lift_det
    (a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    Matrix.det (gL * gS) = (a * b) * (α * δ - β * γ) ^ 2 := by
  intro G gL gS
  rw [Matrix.det_mul]
  rw [show Matrix.det gL = (α * δ - β * γ) ^ 2 by
      simpa [gL] using GL2_lift_det (k := k) α β γ δ]
  rw [show Matrix.det gS = a * b by
      simpa [gS] using pointwise_shape_det (k := k) a b x y p q]
  ring

/-- The full block-diagonal pure singular pointwise shape together with the coupled
unipotent family combines with the simultaneous `GL₂` lift exactly as expected. -/
theorem coupledShape_GL2_product_lift_action
    (u v z w r s a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let gU : Matrix V V k := coupledUnipotent (k := k) u v z w r s
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let gP : Matrix V V k := gU * gS
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    ActBivector rep₁ (gP * gL) = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ (gP * gL) = γ • rep₁ + δ • rep₂ := by
  intro gU gS gP G gL
  have hP := coupled_shape_product_pointwise (k := k) u v z w r s a b x y p q ha hb
  have hL := GL2_lift_action (k := k) α β γ δ
  constructor
  · calc
      ActBivector rep₁ (gP * gL) = ActBivector (ActBivector rep₁ gL) gP := by
        rw [actBivector_mul]
      _ = ActBivector (α • rep₁ + β • rep₂) gP := by
        simpa [gL] using congrArg (fun M => ActBivector M gP) hL.1
      _ = α • ActBivector rep₁ gP + β • ActBivector rep₂ gP := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = α • rep₁ + β • rep₂ := by
        rw [hP.1, hP.2]
  · calc
      ActBivector rep₂ (gP * gL) = ActBivector (ActBivector rep₂ gL) gP := by
        rw [actBivector_mul]
      _ = ActBivector (γ • rep₁ + δ • rep₂) gP := by
        simpa [gL] using congrArg (fun M => ActBivector M gP) hL.2
      _ = γ • ActBivector rep₁ gP + δ • ActBivector rep₂ gP := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = γ • rep₁ + δ • rep₂ := by
        rw [hP.1, hP.2]

/-- Right-multiplying the simultaneous `GL₂` lift by the combined coupled-and-shape
pointwise subgroup does not change the quotient action. -/
theorem GL2_coupledShape_right_product_lift_action
    (u v z w r s a b x y p q α β γ δ : k)
    (ha : a ≠ 0)
    (hb : b ≠ 0) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) α β γ δ (1 : k)
    let gL : Matrix V V k := Matrix.fromBlocks G 0 0 G
    let gU : Matrix V V k := coupledUnipotent (k := k) u v z w r s
    let gS : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        0
        0
        (N3PureSingular.pureSingularShape (k := k) b p q)
    let gP : Matrix V V k := gU * gS
    ActBivector rep₁ (gL * gP) = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ (gL * gP) = γ • rep₁ + δ • rep₂ := by
  intro G gL gU gS gP
  have hL := GL2_lift_action (k := k) α β γ δ
  have hP := coupled_shape_product_pointwise (k := k) u v z w r s a b x y p q ha hb
  constructor
  · calc
      ActBivector rep₁ (gL * gP) = ActBivector (ActBivector rep₁ gP) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gL := by
        simpa [gP] using congrArg (fun M => ActBivector M gL) hP.1
      _ = α • rep₁ + β • rep₂ := by
        simpa [gL] using hL.1
  · calc
      ActBivector rep₂ (gL * gP) = ActBivector (ActBivector rep₂ gP) gL := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gL := by
        simpa [gP] using congrArg (fun M => ActBivector M gL) hP.2
      _ = γ • rep₁ + δ • rep₂ := by
        simpa [gL] using hL.2

end N6DoublePureSingular
end Wedge2Formalization
