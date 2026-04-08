import Wedge2Formalization.N3PureSingular
import Wedge2Formalization.N4
import Mathlib.LinearAlgebra.Matrix.Block

open Matrix

namespace Wedge2Formalization
namespace N5SimplePoint

variable {k : Type*} [Field k]

/-- The `3`-dimensional pure singular summand. -/
abbrev W := N3PureSingular.I

/-- The `2`-dimensional simple summand. -/
abbrev I := N4.I

/-- The ambient `5`-dimensional space. -/
abbrev V := W ⊕ I

/-- Exact stabilizer of a bivector. -/
def FixesBivector (Ω : Matrix V V k) (g : Matrix V V k) : Prop :=
  g * Ω * gᵀ = Ω

/-- The natural action on bivectors. -/
def ActBivector (Ω : Matrix V V k) (g : Matrix V V k) : Matrix V V k :=
  g * Ω * gᵀ

/-- The bivector action is multiplicative in the acting matrix. -/
theorem actBivector_mul
    (Ω : Matrix V V k)
    (g h : Matrix V V k) :
    ActBivector Ω (g * h) = ActBivector (ActBivector Ω h) g := by
  simp [ActBivector, Matrix.mul_assoc]

/-- The direct-sum representative of divisor type `[a]`. -/
def rep₁ : Matrix V V k :=
  Matrix.fromBlocks (N3PureSingular.ω13 (k := k)) 0 0 (N4.J (k := k))

/-- The second basis vector for the direct-sum representative. -/
def rep₂ : Matrix V V k :=
  Matrix.fromBlocks (N3PureSingular.ω23 (k := k)) 0 0 0

/-- Exact stabilizer of the pair. -/
def FixesPairBivector (g : Matrix V V k) : Prop :=
  FixesBivector rep₁ g ∧ FixesBivector rep₂ g

/-- Setwise stabilizer of the `2`-space in the chosen basis. -/
def PreservesSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector rep₁ g = α • rep₁ + β • rep₂ ∧
      ActBivector rep₂ g = γ • rep₁ + δ • rep₂

/-- The basic `4`-parameter unipotent family in the pointwise kernel. -/
def pointwiseUnipotent (x y p q : k) : Matrix V V k :=
  Matrix.fromBlocks
    (!![(1 : k), 0, 0;
        0, 1, 0;
        x, y, 1])
    (!![(0 : k), 0;
        0, 0;
        p, q])
    (!![-q, 0, 0;
        p, 0, 0])
    (1 : Matrix I I k)

/-- Block diagonal action respects the direct-sum decomposition. -/
theorem act_blockDiagonal
    (G : Matrix W W k)
    (E : Matrix I I k) :
    ActBivector rep₁ (Matrix.fromBlocks G 0 0 E) =
      Matrix.fromBlocks
        (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G)
        0
        0
        (E * N4.J * Eᵀ) ∧
    ActBivector rep₂ (Matrix.fromBlocks G 0 0 E) =
      Matrix.fromBlocks
        (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G)
        0
        0
        0 := by
  constructor <;>
  ·
    ext i j
    cases i <;> cases j <;>
      simp [ActBivector, rep₁, rep₂, N3PureSingular.ActBivector,
        Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply]

/-- The obvious Levi subgroup fixes the pair pointwise. -/
theorem pointwise_levi_family
    (a : k)
    (ha : a ≠ 0)
    (E : Matrix I I k)
    (hE : E.det = 1) :
    FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        E) := by
  constructor
  ·
    have hact := (act_blockDiagonal (k := k)
      (N3PureSingular.pointwiseScale (k := k) a)
      E).1
    have htop :
        N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k))
          (N3PureSingular.pointwiseScale (k := k) a) =
          N3PureSingular.ω13 (k := k) := by
      exact (N3PureSingular.pointwise_scale_family (k := k) a ha).1
    have hbot : E * N4.J * Eᵀ = N4.J := by
      simpa [N4.J, hE] using (N4.mul_J_transpose_mul (k := k) E)
    rw [FixesBivector]
    calc
      ActBivector rep₁ (Matrix.fromBlocks (N3PureSingular.pointwiseScale (k := k) a) 0 0 E) =
          Matrix.fromBlocks (N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k))
            (N3PureSingular.pointwiseScale (k := k) a)) 0 0 (E * N4.J * Eᵀ) := hact
      _ = rep₁ := by simpa [rep₁, htop, hbot]
  ·
    have hact := (act_blockDiagonal (k := k)
      (N3PureSingular.pointwiseScale (k := k) a)
      E).2
    have htop :
        N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k))
          (N3PureSingular.pointwiseScale (k := k) a) =
          N3PureSingular.ω23 (k := k) := by
      exact (N3PureSingular.pointwise_scale_family (k := k) a ha).2
    rw [FixesBivector]
    calc
      ActBivector rep₂ (Matrix.fromBlocks (N3PureSingular.pointwiseScale (k := k) a) 0 0 E) =
          Matrix.fromBlocks (N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k))
            (N3PureSingular.pointwiseScale (k := k) a)) 0 0 0 := hact
      _ = rep₂ := by simpa [rep₂, htop]

/-- The basic `4`-parameter unipotent family fixes the direct-sum `[a]` pair
pointwise. -/
theorem pointwise_unipotent_family
    (x y p q : k) :
    FixesPairBivector (pointwiseUnipotent (k := k) x y p q) := by
  constructor <;> rw [FixesBivector]
  ·
    ext i j
    rcases i with i | i <;> rcases j with j | j
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₁, pointwiseUnipotent, N3PureSingular.ω13, N4.J,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₁, pointwiseUnipotent, N3PureSingular.ω13, N4.J,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₁, pointwiseUnipotent, N3PureSingular.ω13, N4.J,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₁, pointwiseUnipotent, N3PureSingular.ω13, N4.J,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf
  ·
    ext i j
    rcases i with i | i <;> rcases j with j | j
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₂, pointwiseUnipotent, N3PureSingular.ω23,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₂, pointwiseUnipotent, N3PureSingular.ω23,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₂, pointwiseUnipotent, N3PureSingular.ω23,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf
    ·
      fin_cases i <;> fin_cases j <;>
        simp [rep₂, pointwiseUnipotent, N3PureSingular.ω23,
          Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply,
          Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three]
        <;> ring_nf

/-- Combining the pointwise unipotent family with the Levi family gives a concrete
pointwise subgroup of the direct-sum `[a]` stabilizer. -/
theorem pointwise_product_family
    (x y p q a : k)
    (ha : a ≠ 0)
    (E : Matrix I I k)
    (hE : E.det = 1) :
    let gU : Matrix V V k := pointwiseUnipotent (k := k) x y p q
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        E
    FixesPairBivector (gU * gL) := by
  intro gU gL
  have hU := pointwise_unipotent_family (k := k) x y p q
  have hL := pointwise_levi_family (k := k) a ha E hE
  constructor <;> rw [FixesBivector]
  ·
    calc
      ActBivector rep₁ (gU * gL) = ActBivector (ActBivector rep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₁ gU := by
        simpa [FixesBivector] using congrArg (fun M => ActBivector M gU) hL.1
      _ = rep₁ := by simpa [FixesBivector] using hU.1
  ·
    calc
      ActBivector rep₂ (gU * gL) = ActBivector (ActBivector rep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector rep₂ gU := by
        simpa [FixesBivector] using congrArg (fun M => ActBivector M gU) hL.2
      _ = rep₂ := by simpa [FixesBivector] using hU.2

/-- Every matrix of the expected block shape on the direct-sum `[a]` orbit fixes the
pair pointwise. This packages the pointwise unipotent radical together with the
`G_m x SL_2` Levi factor in one explicit formula. -/
theorem pointwise_shape_family
    (a x y u v : k)
    (ha : a ≠ 0)
    (E : Matrix I I k)
    (hE : E.det = 1) :
    let B : Matrix W I k := !![(0 : k), 0; 0, 0; u, v]
    let C : Matrix I W k :=
      !![a * (u * E 0 1 - v * E 0 0), 0, 0;
         a * (u * E 1 1 - v * E 1 0), 0, 0]
    FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        B
        C
        E) := by
  intro B C
  let p : k := u * E 1 1 - v * E 1 0
  let q : k := v * E 0 0 - u * E 0 1
  let xU : k := x * a⁻¹
  let yU : k := y * a⁻¹
  let A : Matrix W W k :=
    !![(1 : k), 0, 0;
       0, 1, 0;
       xU, yU, 1]
  let BU : Matrix W I k :=
    !![(0 : k), 0;
       0, 0;
       p, q]
  let CU : Matrix I W k :=
    !![-q, 0, 0;
       p, 0, 0]
  have hprod :=
    pointwise_product_family
      (k := k)
      (x := xU)
      (y := yU)
      (p := p)
      (q := q)
      (a := a⁻¹)
      (ha := inv_ne_zero ha)
      (E := E)
      hE
  have hdet0 : E 0 0 * E 1 1 - E 0 1 * E 1 0 = 1 := by
    simpa [Matrix.det_fin_two] using hE
  have hdet : E 1 1 * E 0 0 - E 1 0 * E 0 1 = 1 := by
    calc
      E 1 1 * E 0 0 - E 1 0 * E 0 1 = E 0 0 * E 1 1 - E 0 1 * E 1 0 := by ring
      _ = 1 := hdet0
  have hA :
      A * N3PureSingular.pointwiseScale (k := k) a⁻¹ =
        N3PureSingular.pureSingularShape (k := k) a x y := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [A, xU, yU, N3PureSingular.pointwiseScale, N3PureSingular.pureSingularShape,
        Matrix.mul_apply, Fin.sum_univ_three, ha] <;>
      field_simp [ha] <;>
      ring
  have hBU :
      BU * E = B := by
    ext i j
    fin_cases i <;> fin_cases j
    · simp [BU, B, Matrix.mul_apply, Fin.sum_univ_two]
    · simp [BU, B, Matrix.mul_apply, Fin.sum_univ_two]
    · simp [BU, B, Matrix.mul_apply, Fin.sum_univ_two]
    · simp [BU, B, Matrix.mul_apply, Fin.sum_univ_two]
    ·
      calc
        (BU * E) 2 0 = p * E 0 0 + q * E 1 0 := by
          simp [BU, Matrix.mul_apply, Fin.sum_univ_two]
        _ = u := by
          calc
            (u * E 1 1 - v * E 1 0) * E 0 0 + (v * E 0 0 - u * E 0 1) * E 1 0 =
                u * (E 1 1 * E 0 0 - E 1 0 * E 0 1) := by ring
            _ = u := by rw [hdet]; ring
          <;> simp [p, q]
        _ = B 2 0 := by simp [B]
    ·
      calc
        (BU * E) 2 1 = p * E 0 1 + q * E 1 1 := by
          simp [BU, Matrix.mul_apply, Fin.sum_univ_two]
        _ = v := by
          calc
            (u * E 1 1 - v * E 1 0) * E 0 1 + (v * E 0 0 - u * E 0 1) * E 1 1 =
                v * (E 1 1 * E 0 0 - E 1 0 * E 0 1) := by ring
            _ = v := by rw [hdet]; ring
          <;> simp [p, q]
        _ = B 2 1 := by simp [B]
  have hCU :
      CU * N3PureSingular.pointwiseScale (k := k) a⁻¹ = C := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [CU, C, p, q, N3PureSingular.pointwiseScale, Matrix.mul_apply,
        Fin.sum_univ_three, ha] <;>
      field_simp [ha] <;>
      ring
  have hEq :
      pointwiseUnipotent (k := k) xU yU p q *
          Matrix.fromBlocks
            (N3PureSingular.pointwiseScale (k := k) a⁻¹)
            0
            0
            E =
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C
          E := by
    calc
      pointwiseUnipotent (k := k) xU yU p q *
          Matrix.fromBlocks
            (N3PureSingular.pointwiseScale (k := k) a⁻¹)
            0
            0
            E =
        Matrix.fromBlocks
          (A * N3PureSingular.pointwiseScale (k := k) a⁻¹)
          (BU * E)
          (CU * N3PureSingular.pointwiseScale (k := k) a⁻¹)
          E := by
            simp [pointwiseUnipotent, A, BU, CU, Matrix.fromBlocks_multiply]
      _ =
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C
          E := by rw [hA, hBU, hCU]
  simpa [hEq] using hprod

/-- Once the upper-left block is already in the pure singular shape, pointwise fixing
forces the remaining blocks on the direct-sum `[a]` orbit to have the expected form. -/
theorem fixesPair_fromBlocks_pureShape_iff
    (a x y : k)
    (ha : a ≠ 0)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    FixesPairBivector
      (Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        B
        C
        D) ↔
      D.det = 1 ∧
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C
          D =
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          (!![(0 : k), 0; 0, 0; B 2 0, B 2 1])
          (!![a * (B 2 0 * D 0 1 - B 2 1 * D 0 0), 0, 0;
             a * (B 2 0 * D 1 1 - B 2 1 * D 1 0), 0, 0])
          D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    have h₂tr :
        N3PureSingular.pureSingularShape (k := k) a x y * N3PureSingular.ω23 (k := k) * Cᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₂
      simpa [FixesBivector, rep₂, N3PureSingular.pureSingularShape, N3PureSingular.ω23,
        Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have hC02 : C 0 2 = 0 := by
      have hij := congrArg (fun M => M 1 0) h₂tr
      simpa [N3PureSingular.pureSingularShape, N3PureSingular.ω23, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, ha] using hij
    have hC12 : C 1 2 = 0 := by
      have hij := congrArg (fun M => M 1 1) h₂tr
      simpa [N3PureSingular.pureSingularShape, N3PureSingular.ω23, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, ha] using hij
    have hC01cases : a = 0 ∨ C 0 1 = 0 := by
      have hij := congrArg (fun M => M 2 0) h₂tr
      simpa [N3PureSingular.pureSingularShape, N3PureSingular.ω23, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, hC02] using hij
    have hC01 : C 0 1 = 0 := hC01cases.resolve_left ha
    have hC11cases : a = 0 ∨ C 1 1 = 0 := by
      have hij := congrArg (fun M => M 2 1) h₂tr
      simpa [N3PureSingular.pureSingularShape, N3PureSingular.ω23, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, hC12] using hij
    have hC11 : C 1 1 = 0 := hC11cases.resolve_left ha
    have h₁br :
        C * N3PureSingular.ω13 (k := k) * Cᵀ + D * N4.J * Dᵀ = N4.J := by
      have h' := congrArg Matrix.toBlocks₂₂ h₁
      simpa [FixesBivector, rep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using h'
    have hCω : C * N3PureSingular.ω13 (k := k) * Cᵀ = 0 := by
      ext i j
      fin_cases i <;> fin_cases j
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
    have hDJ : D * N4.J * Dᵀ = N4.J := by
      rw [hCω, zero_add] at h₁br
      exact h₁br
    have hdetD : D.det = 1 := by
      rw [N4.mul_J_transpose_mul] at hDJ
      have hij := congrArg (fun M => M 0 1) hDJ
      simpa [N4.J] using hij
    let T : Matrix W I k :=
      N3PureSingular.pureSingularShape (k := k) a x y * N3PureSingular.ω13 (k := k) * Cᵀ
    let U : Matrix W I k := B * N4.J * Dᵀ
    have h₁tr :
        T + U = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₁
      simpa [T, U, FixesBivector, rep₁, N3PureSingular.pureSingularShape, N3PureSingular.ω13,
        Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have hDdet0 : D 0 0 * D 1 1 - D 0 1 * D 1 0 = 1 := by
      simpa [Matrix.det_fin_two] using hdetD
    have hDdet0swap : D 0 0 * D 1 1 - D 1 0 * D 0 1 = 1 := by
      simpa [mul_comm] using hDdet0
    have hDdet0neg : D 0 1 * D 1 0 - D 1 1 * D 0 0 = -1 := by
      calc
        D 0 1 * D 1 0 - D 1 1 * D 0 0 = -(D 0 0 * D 1 1 - D 0 1 * D 1 0) := by ring
        _ = -1 := by rw [hDdet0]
    have hvecB00 : (![0, 0, a] ᵥ* Cᵀ) 0 = 0 := by
      simp [Matrix.vecMul, dotProduct, Matrix.transpose_apply, Fin.sum_univ_three, hC02]
    have hvecB01 : (![0, 0, a] ᵥ* Cᵀ) 1 = 0 := by
      simp [Matrix.vecMul, dotProduct, Matrix.transpose_apply, Fin.sum_univ_three, hC12]
    have hvecZero0 : (![0, 0, 0] ᵥ* Cᵀ) 0 = 0 := by
      simp [Matrix.vecMul, dotProduct, Matrix.transpose_apply, Fin.sum_univ_three]
    have hvecZero1 : (![0, 0, 0] ᵥ* Cᵀ) 1 = 0 := by
      simp [Matrix.vecMul, dotProduct, Matrix.transpose_apply, Fin.sum_univ_three]
    have hvecC00 : (![-a⁻¹, 0, x] ᵥ* Cᵀ) 0 = -(a⁻¹ * C 0 0) := by
      simp [Matrix.vecMul, dotProduct, Matrix.transpose_apply, Fin.sum_univ_three, hC02]
    have hvecC10 : (![-a⁻¹, 0, x] ᵥ* Cᵀ) 1 = -(a⁻¹ * C 1 0) := by
      simp [Matrix.vecMul, dotProduct, Matrix.transpose_apply, Fin.sum_univ_three, hC12]
    have hB00eq : -B 0 1 * D 0 0 + B 0 0 * D 0 1 = 0 := by
      have hij := congrArg (fun M => M 0 0) h₁tr
      simp [T, U, N3PureSingular.pureSingularShape, N3PureSingular.ω13, N4.J, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_two, hC01, hC02, ha] at hij
      rw [hvecB00] at hij
      simpa using hij
    have hB01eq : -B 0 1 * D 1 0 + B 0 0 * D 1 1 = 0 := by
      have hij := congrArg (fun M => M 0 1) h₁tr
      simp [T, U, N3PureSingular.pureSingularShape, N3PureSingular.ω13, N4.J, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_two, hC01, hC02, ha] at hij
      rw [hvecB01] at hij
      simpa using hij
    have hB10eq : -B 1 1 * D 0 0 + B 1 0 * D 0 1 = 0 := by
      have hij := congrArg (fun M => M 1 0) h₁tr
      simp [T, U, N3PureSingular.pureSingularShape, N3PureSingular.ω13, N4.J, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_two, hC11, hC12, ha] at hij
      rw [hvecZero0] at hij
      simpa using hij
    have hB11eq : -B 1 1 * D 1 0 + B 1 0 * D 1 1 = 0 := by
      have hij := congrArg (fun M => M 1 1) h₁tr
      simp [T, U, N3PureSingular.pureSingularShape, N3PureSingular.ω13, N4.J, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_two, hC11, hC12, ha] at hij
      rw [hvecZero1] at hij
      simpa using hij
    have hB01 : B 0 1 = 0 := by
      have htmp : -B 0 1 * (D 0 0 * D 1 1 - D 1 0 * D 0 1) = 0 := by
        calc
          -B 0 1 * (D 0 0 * D 1 1 - D 1 0 * D 0 1)
              = (-B 0 1 * D 0 0 + B 0 0 * D 0 1) * D 1 1 -
                  (-B 0 1 * D 1 0 + B 0 0 * D 1 1) * D 0 1 := by ring
          _ = 0 := by rw [hB00eq, hB01eq]; ring
      rw [hDdet0swap] at htmp
      simpa using htmp
    have hB00 : B 0 0 = 0 := by
      have htmp : B 0 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) = 0 := by
        calc
          B 0 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0)
              = (-B 0 1 * D 0 0 + B 0 0 * D 0 1) * D 1 0 -
                  (-B 0 1 * D 1 0 + B 0 0 * D 1 1) * D 0 0 := by ring
          _ = 0 := by rw [hB00eq, hB01eq]; ring
      have : -B 0 0 = 0 := by
        rw [hDdet0neg] at htmp
        simpa using htmp
      simpa using this
    have hB11 : B 1 1 = 0 := by
      have htmp : -B 1 1 * (D 0 0 * D 1 1 - D 1 0 * D 0 1) = 0 := by
        calc
          -B 1 1 * (D 0 0 * D 1 1 - D 1 0 * D 0 1)
              = (-B 1 1 * D 0 0 + B 1 0 * D 0 1) * D 1 1 -
                  (-B 1 1 * D 1 0 + B 1 0 * D 1 1) * D 0 1 := by ring
          _ = 0 := by rw [hB10eq, hB11eq]; ring
      rw [hDdet0swap] at htmp
      simpa using htmp
    have hB10 : B 1 0 = 0 := by
      have htmp : B 1 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) = 0 := by
        calc
          B 1 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0)
              = (-B 1 1 * D 0 0 + B 1 0 * D 0 1) * D 1 0 -
                  (-B 1 1 * D 1 0 + B 1 0 * D 1 1) * D 0 0 := by ring
          _ = 0 := by rw [hB10eq, hB11eq]; ring
      have : -B 1 0 = 0 := by
        rw [hDdet0neg] at htmp
        simpa using htmp
      simpa using this
    have hC00eq : -(a⁻¹ * C 0 0) + (-B 2 1 * D 0 0 + B 2 0 * D 0 1) = 0 := by
      have hij := congrArg (fun M => M 2 0) h₁tr
      simp [T, U, N3PureSingular.pureSingularShape, N3PureSingular.ω13, N4.J, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_two, hC01, hC02, ha, hB00, hB01] at hij
      rw [hvecC00] at hij
      simpa [add_assoc] using hij
    have hC10eq : -(a⁻¹ * C 1 0) + (-B 2 1 * D 1 0 + B 2 0 * D 1 1) = 0 := by
      have hij := congrArg (fun M => M 2 1) h₁tr
      simp [T, U, N3PureSingular.pureSingularShape, N3PureSingular.ω13, N4.J, Matrix.mul_apply,
        Matrix.transpose_apply, Fin.sum_univ_three, Fin.sum_univ_two, hC11, hC12, ha, hB10, hB11] at hij
      rw [hvecC10] at hij
      simpa [add_assoc] using hij
    have hC00 : C 0 0 = a * (B 2 0 * D 0 1 - B 2 1 * D 0 0) := by
      calc
        C 0 0 = a * (a⁻¹ * C 0 0) := by field_simp [ha]
        _ = a * (B 2 0 * D 0 1 - B 2 1 * D 0 0) := by
          have : a⁻¹ * C 0 0 = B 2 0 * D 0 1 - B 2 1 * D 0 0 := by
            have htmp : -(a⁻¹ * C 0 0) = -(-B 2 1 * D 0 0 + B 2 0 * D 0 1) := by
              exact add_eq_zero_iff_eq_neg.mp hC00eq
            simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using congrArg Neg.neg htmp
          rw [this]
    have hC10 : C 1 0 = a * (B 2 0 * D 1 1 - B 2 1 * D 1 0) := by
      calc
        C 1 0 = a * (a⁻¹ * C 1 0) := by field_simp [ha]
        _ = a * (B 2 0 * D 1 1 - B 2 1 * D 1 0) := by
          have : a⁻¹ * C 1 0 = B 2 0 * D 1 1 - B 2 1 * D 1 0 := by
            have htmp : -(a⁻¹ * C 1 0) = -(-B 2 1 * D 1 0 + B 2 0 * D 1 1) := by
              exact add_eq_zero_iff_eq_neg.mp hC10eq
            simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using congrArg Neg.neg htmp
          rw [this]
    refine ⟨hdetD, ?_⟩
    ext i j
    cases i with
    | inl i =>
        cases j with
        | inl j =>
            fin_cases i <;> fin_cases j <;>
              simp [Matrix.fromBlocks, hC00, hC01, hC02, hC10, hC11, hC12]
        | inr j =>
            fin_cases i <;> fin_cases j <;>
              simp [Matrix.fromBlocks, hB00, hB01, hB10, hB11]
    | inr i =>
        cases j with
        | inl j =>
            fin_cases i <;> fin_cases j <;>
              simp [Matrix.fromBlocks, hC00, hC01, hC02, hC10, hC11, hC12]
        | inr j =>
            fin_cases i <;> fin_cases j <;>
              simp [Matrix.fromBlocks]
  · rintro ⟨hdetD, hshape⟩
    rw [hshape]
    exact
      pointwise_shape_family
        (k := k)
        (a := a)
        (x := x)
        (y := y)
        (u := B 2 0)
        (v := B 2 1)
        ha
        (E := D)
        hdetD

set_option maxHeartbeats 1200000 in
/-- The direct-sum `[a]` pointwise stabilizer has the expected explicit shape. -/
theorem fixesPair_fromBlocks_iff_shape
    (A : Matrix W W k)
    (B : Matrix W I k)
    (C : Matrix I W k)
    (D : Matrix I I k) :
    FixesPairBivector (Matrix.fromBlocks A B C D) ↔
      A 0 0 ≠ 0 ∧
        D.det = 1 ∧
          Matrix.fromBlocks A B C D =
            Matrix.fromBlocks
              (N3PureSingular.pureSingularShape (k := k) (A 0 0) (A 2 0) (A 2 1))
              (!![(0 : k), 0; 0, 0; B 2 0, B 2 1])
              (!![(A 0 0) * (B 2 0 * D 0 1 - B 2 1 * D 0 0), 0, 0;
                 (A 0 0) * (B 2 0 * D 1 1 - B 2 1 * D 1 0), 0, 0])
              D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    have hA23eq :
        N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) A =
          N3PureSingular.ω23 (k := k) := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h₂
      simpa [FixesBivector, rep₂, N3PureSingular.ActBivector,
        Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have h2301 : A 0 1 * A 1 2 - A 0 2 * A 1 1 = 0 := by
      have hij := congrArg (fun M => M 0 1) hA23eq
      simp [N3PureSingular.ActBivector, N3PureSingular.ω23, Matrix.mul_apply,
        Fin.sum_univ_three] at hij
      ring_nf at hij ⊢
      exact hij
    have h2302 : A 0 1 * A 2 2 - A 0 2 * A 2 1 = 0 := by
      have hij := congrArg (fun M => M 0 2) hA23eq
      simp [N3PureSingular.ActBivector, N3PureSingular.ω23, Matrix.mul_apply,
        Fin.sum_univ_three] at hij
      ring_nf at hij ⊢
      exact hij
    have h2312 : A 1 1 * A 2 2 - A 1 2 * A 2 1 = 1 := by
      have hij := congrArg (fun M => M 1 2) hA23eq
      simp [N3PureSingular.ActBivector, N3PureSingular.ω23, Matrix.mul_apply,
        Fin.sum_univ_three] at hij
      ring_nf at hij ⊢
      exact hij
    have hA02mul : A 0 2 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) = 0 := by
      calc
        A 0 2 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) =
            (A 0 2 * A 1 1 - A 0 1 * A 1 2) * A 2 2 -
              (A 0 2 * A 2 1 - A 0 1 * A 2 2) * A 1 2 := by
                ring
        _ = 0 := by
              have h1 : A 0 2 * A 1 1 - A 0 1 * A 1 2 = 0 := by
                simpa using congrArg Neg.neg h2301
              have h2 : A 0 2 * A 2 1 - A 0 1 * A 2 2 = 0 := by
                simpa using congrArg Neg.neg h2302
              rw [h1, h2]
              ring
    have hA02 : A 0 2 = 0 := by
      have : A 0 2 * 1 = 0 := by simpa [h2312] using hA02mul
      simpa using this
    have hA01mul : A 0 1 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) = 0 := by
      calc
        A 0 1 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) =
            (A 0 1 * A 2 2 - A 0 2 * A 2 1) * A 1 1 -
              (A 0 1 * A 1 2 - A 0 2 * A 1 1) * A 2 1 := by
                ring
        _ = 0 := by
              rw [h2302, h2301]
              ring
    have hA01 : A 0 1 = 0 := by
      have : A 0 1 * 1 = 0 := by simpa [h2312] using hA01mul
      simpa using this
    have h₂tr : A * N3PureSingular.ω23 (k := k) * Cᵀ = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₂
      simpa [FixesBivector, rep₂, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have hC02eq : -A 1 2 * C 0 1 + A 1 1 * C 0 2 = 0 := by
      have hij := congrArg (fun M => M 1 0) h₂tr
      simp [N3PureSingular.ω23, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02] at hij
      ring_nf at hij ⊢
      exact hij
    have hC01eq : -A 2 2 * C 0 1 + A 2 1 * C 0 2 = 0 := by
      have hij := congrArg (fun M => M 2 0) h₂tr
      simp [N3PureSingular.ω23, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02] at hij
      ring_nf at hij ⊢
      exact hij
    have hC12eq : -A 1 2 * C 1 1 + A 1 1 * C 1 2 = 0 := by
      have hij := congrArg (fun M => M 1 1) h₂tr
      simp [N3PureSingular.ω23, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02] at hij
      ring_nf at hij ⊢
      exact hij
    have hC11eq : -A 2 2 * C 1 1 + A 2 1 * C 1 2 = 0 := by
      have hij := congrArg (fun M => M 2 1) h₂tr
      simp [N3PureSingular.ω23, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02] at hij
      ring_nf at hij ⊢
      exact hij
    have hC01mul : C 0 1 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) = 0 := by
      calc
        C 0 1 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) =
            (A 2 2 * C 0 1 - A 2 1 * C 0 2) * A 1 1 -
              (A 1 2 * C 0 1 - A 1 1 * C 0 2) * A 2 1 := by
                ring
        _ = 0 := by
              have h1 : A 2 2 * C 0 1 - A 2 1 * C 0 2 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using
                  congrArg Neg.neg hC01eq
              have h2 : A 1 2 * C 0 1 - A 1 1 * C 0 2 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using
                  congrArg Neg.neg hC02eq
              rw [h1, h2]
              ring
    have hC01 : C 0 1 = 0 := by
      have : C 0 1 * 1 = 0 := by simpa [h2312] using hC01mul
      simpa using this
    have hC02mul : C 0 2 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) = 0 := by
      calc
        C 0 2 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) =
            (A 1 1 * C 0 2 - A 1 2 * C 0 1) * A 2 2 -
              (A 2 1 * C 0 2 - A 2 2 * C 0 1) * A 1 2 := by
                ring
        _ = 0 := by
              have h1 : A 1 1 * C 0 2 - A 1 2 * C 0 1 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hC02eq
              have h2 : A 2 1 * C 0 2 - A 2 2 * C 0 1 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hC01eq
              rw [h1, h2]
              ring
    have hC02 : C 0 2 = 0 := by
      have : C 0 2 * 1 = 0 := by simpa [h2312] using hC02mul
      simpa using this
    have hC11mul : C 1 1 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) = 0 := by
      calc
        C 1 1 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) =
            (A 2 2 * C 1 1 - A 2 1 * C 1 2) * A 1 1 -
              (A 1 2 * C 1 1 - A 1 1 * C 1 2) * A 2 1 := by
                ring
        _ = 0 := by
              have h1 : A 2 2 * C 1 1 - A 2 1 * C 1 2 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using
                  congrArg Neg.neg hC11eq
              have h2 : A 1 2 * C 1 1 - A 1 1 * C 1 2 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using
                  congrArg Neg.neg hC12eq
              rw [h1, h2]
              ring
    have hC11 : C 1 1 = 0 := by
      have : C 1 1 * 1 = 0 := by simpa [h2312] using hC11mul
      simpa using this
    have hC12mul : C 1 2 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) = 0 := by
      calc
        C 1 2 * (A 1 1 * A 2 2 - A 1 2 * A 2 1) =
            (A 1 1 * C 1 2 - A 1 2 * C 1 1) * A 2 2 -
              (A 2 1 * C 1 2 - A 2 2 * C 1 1) * A 1 2 := by
                ring
        _ = 0 := by
              have h1 : A 1 1 * C 1 2 - A 1 2 * C 1 1 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hC12eq
              have h2 : A 2 1 * C 1 2 - A 2 2 * C 1 1 = 0 := by
                simpa [sub_eq_add_neg, add_comm, add_left_comm, add_assoc] using hC11eq
              rw [h1, h2]
              ring
    have hC12 : C 1 2 = 0 := by
      have : C 1 2 * 1 = 0 := by simpa [h2312] using hC12mul
      simpa using this
    have hCω : C * N3PureSingular.ω13 (k := k) * Cᵀ = 0 := by
      ext i j
      fin_cases i <;> fin_cases j
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
      ·
        simp [N3PureSingular.ω13, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_two, Fin.sum_univ_three, hC01, hC02, hC11, hC12]
    have h₁br :
        C * N3PureSingular.ω13 (k := k) * Cᵀ + D * N4.J * Dᵀ = N4.J := by
      have h' := congrArg Matrix.toBlocks₂₂ h₁
      simpa [FixesBivector, rep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using h'
    have hDJ : D * N4.J * Dᵀ = N4.J := by
      rw [hCω, zero_add] at h₁br
      exact h₁br
    have hDdet : D.det = 1 := by
      rw [N4.mul_J_transpose_mul] at hDJ
      have hij := congrArg (fun M => M 0 1) hDJ
      simpa [N4.J] using hij
    have hDdet0 : D 0 0 * D 1 1 - D 0 1 * D 1 0 = 1 := by
      simpa [Matrix.det_fin_two] using hDdet
    have hDdet0swap : D 0 0 * D 1 1 - D 1 0 * D 0 1 = 1 := by
      simpa [mul_comm] using hDdet0
    have hDdet0neg : D 0 1 * D 1 0 - D 1 1 * D 0 0 = -1 := by
      calc
        D 0 1 * D 1 0 - D 1 1 * D 0 0 = -(D 0 0 * D 1 1 - D 0 1 * D 1 0) := by ring
        _ = -1 := by rw [hDdet0]
    let T : Matrix W I k := A * N3PureSingular.ω13 (k := k) * Cᵀ
    let U : Matrix W I k := B * N4.J * Dᵀ
    have h₁tr :
        T + U = 0 := by
      ext i j
      have hij := congrArg (fun M => M (Sum.inl i) (Sum.inr j)) h₁
      simpa [T, U, FixesBivector, rep₁, Matrix.fromBlocks_transpose,
        Matrix.fromBlocks_multiply] using hij
    have hB00eq : -B 0 1 * D 0 0 + B 0 0 * D 0 1 = 0 := by
      have hij := congrArg (fun M => M 0 0) h₁tr
      simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_three, Fin.sum_univ_two, hA02, hC01, hC02] at hij
      ring_nf at hij ⊢
      simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
    have hB01eq : -B 0 1 * D 1 0 + B 0 0 * D 1 1 = 0 := by
      have hij := congrArg (fun M => M 0 1) h₁tr
      simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
        Fin.sum_univ_three, Fin.sum_univ_two, hA02, hC11, hC12] at hij
      ring_nf at hij ⊢
      simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
    have hB00mul : B 0 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) = 0 := by
      calc
        B 0 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) =
            (-B 0 1 * D 0 0 + B 0 0 * D 0 1) * D 1 0 -
              (-B 0 1 * D 1 0 + B 0 0 * D 1 1) * D 0 0 := by
                ring_nf
        _ = 0 := by rw [hB00eq, hB01eq]; ring
    have hB00 : B 0 0 = 0 := by
      have : -B 0 0 = 0 := by
        rw [hDdet0neg] at hB00mul
        simpa using hB00mul
      simpa using this
    have hB01D00 : -B 0 1 * D 0 0 = 0 := by simpa [hB00] using hB00eq
    have hB01D10 : -B 0 1 * D 1 0 = 0 := by simpa [hB00] using hB01eq
    have hB01mul : -B 0 1 * (D 0 0 * D 1 1 - D 1 0 * D 0 1) = 0 := by
      calc
        -B 0 1 * (D 0 0 * D 1 1 - D 1 0 * D 0 1) =
            (-B 0 1 * D 0 0) * D 1 1 - (-B 0 1 * D 1 0) * D 0 1 := by ring
        _ = 0 := by rw [hB01D00, hB01D10]; ring
    have hB01 : B 0 1 = 0 := by
      rw [hDdet0swap] at hB01mul
      have : -B 0 1 = 0 := by simpa using hB01mul
      simpa using this
    have hB10formula :
        B 1 0 = A 1 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) := by
      have hEq10 : -A 1 2 * C 0 0 - B 1 1 * D 0 0 + B 1 0 * D 0 1 = 0 := by
        have hij := congrArg (fun M => M 1 0) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02, hC01, hC02] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hEq11 : -A 1 2 * C 1 0 - B 1 1 * D 1 0 + B 1 0 * D 1 1 = 0 := by
        have hij := congrArg (fun M => M 1 1) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02, hC11, hC12] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hTmp :
          -A 1 2 * C 0 0 * D 1 0 + A 1 2 * C 1 0 * D 0 0 +
            B 1 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) = 0 := by
        calc
          -A 1 2 * C 0 0 * D 1 0 + A 1 2 * C 1 0 * D 0 0 +
              B 1 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) =
            (-A 1 2 * C 0 0 - B 1 1 * D 0 0 + B 1 0 * D 0 1) * D 1 0 -
              (-A 1 2 * C 1 0 - B 1 1 * D 1 0 + B 1 0 * D 1 1) * D 0 0 := by ring
          _ = 0 := by rw [hEq10, hEq11]; ring
      have hCoeff :
          D 0 1 * D 1 0 - D 1 1 * D 0 0 = -1 := by
        calc
          D 0 1 * D 1 0 - D 1 1 * D 0 0 = -(D 0 0 * D 1 1 - D 0 1 * D 1 0) := by ring
          _ = -1 := by rw [hDdet0]
      have hTmp' :
          B 1 0 = A 1 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) := by
        have : -B 1 0 = -(A 1 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0)) := by
          rw [hCoeff] at hTmp
          have hTmp0 : -B 1 0 + A 1 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) = 0 := by
            calc
              -B 1 0 + A 1 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) =
                  -B 1 0 + (A 1 2 * C 1 0 * D 0 0 + -(A 1 2 * C 0 0 * D 1 0)) := by ring
              _ = 0 := by
                simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm, mul_assoc,
                  mul_left_comm, mul_comm] using hTmp
          exact add_eq_zero_iff_eq_neg.mp hTmp0
        exact neg_injective this
      exact hTmp'
    have hB11formula :
        B 1 1 = A 1 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) := by
      have hEq10 : -A 1 2 * C 0 0 - B 1 1 * D 0 0 + B 1 0 * D 0 1 = 0 := by
        have hij := congrArg (fun M => M 1 0) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02, hC01, hC02] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hEq11 : -A 1 2 * C 1 0 - B 1 1 * D 1 0 + B 1 0 * D 1 1 = 0 := by
        have hij := congrArg (fun M => M 1 1) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02, hC11, hC12] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hTmp :
          -A 1 2 * C 0 0 * D 1 1 + A 1 2 * C 1 0 * D 0 1 +
            B 1 1 * (D 1 0 * D 0 1 - D 0 0 * D 1 1) = 0 := by
        calc
          -A 1 2 * C 0 0 * D 1 1 + A 1 2 * C 1 0 * D 0 1 +
              B 1 1 * (D 1 0 * D 0 1 - D 0 0 * D 1 1) =
            (-A 1 2 * C 0 0 - B 1 1 * D 0 0 + B 1 0 * D 0 1) * D 1 1 -
              (-A 1 2 * C 1 0 - B 1 1 * D 1 0 + B 1 0 * D 1 1) * D 0 1 := by ring
          _ = 0 := by rw [hEq10, hEq11]; ring
      have hCoeff :
          D 1 0 * D 0 1 - D 0 0 * D 1 1 = -1 := by
        calc
          D 1 0 * D 0 1 - D 0 0 * D 1 1 = -(D 0 0 * D 1 1 - D 0 1 * D 1 0) := by ring
          _ = -1 := by rw [hDdet0]
      have hTmp' :
          B 1 1 = A 1 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) := by
        have : -B 1 1 = -(A 1 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1)) := by
          rw [hCoeff] at hTmp
          have hTmp0 : -B 1 1 + A 1 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) = 0 := by
            calc
              -B 1 1 + A 1 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) =
                  -B 1 1 + (A 1 2 * C 1 0 * D 0 1 + -(A 1 2 * C 0 0 * D 1 1)) := by ring
              _ = 0 := by
                simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm, mul_assoc,
                  mul_left_comm, mul_comm] using hTmp
          exact add_eq_zero_iff_eq_neg.mp hTmp0
        exact neg_injective this
      exact hTmp'
    have hB20formula :
        B 2 0 = A 2 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) := by
      have hEq20 : -A 2 2 * C 0 0 - B 2 1 * D 0 0 + B 2 0 * D 0 1 = 0 := by
        have hij := congrArg (fun M => M 2 0) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02, hC01, hC02] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hEq21 : -A 2 2 * C 1 0 - B 2 1 * D 1 0 + B 2 0 * D 1 1 = 0 := by
        have hij := congrArg (fun M => M 2 1) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hC11, hC12] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hTmp :
          -A 2 2 * C 0 0 * D 1 0 + A 2 2 * C 1 0 * D 0 0 +
            B 2 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) = 0 := by
        calc
          -A 2 2 * C 0 0 * D 1 0 + A 2 2 * C 1 0 * D 0 0 +
              B 2 0 * (D 0 1 * D 1 0 - D 1 1 * D 0 0) =
            (-A 2 2 * C 0 0 - B 2 1 * D 0 0 + B 2 0 * D 0 1) * D 1 0 -
              (-A 2 2 * C 1 0 - B 2 1 * D 1 0 + B 2 0 * D 1 1) * D 0 0 := by ring
          _ = 0 := by rw [hEq20, hEq21]; ring
      have hCoeff :
          D 0 1 * D 1 0 - D 1 1 * D 0 0 = -1 := by
        calc
          D 0 1 * D 1 0 - D 1 1 * D 0 0 = -(D 0 0 * D 1 1 - D 0 1 * D 1 0) := by ring
          _ = -1 := by rw [hDdet0]
      have hTmp' :
          B 2 0 = A 2 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) := by
        have : -B 2 0 = -(A 2 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0)) := by
          rw [hCoeff] at hTmp
          have hTmp0 : A 2 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) + -B 2 0 = 0 := by
            calc
              A 2 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) + -B 2 0 =
                  A 2 2 * C 1 0 * D 0 0 + (-B 2 0 + -(A 2 2 * C 0 0 * D 1 0)) := by ring
              _ = 0 := by
                simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm, mul_assoc,
                  mul_left_comm, mul_comm] using hTmp
          have hEq : A 2 2 * (C 1 0 * D 0 0 - C 0 0 * D 1 0) = B 2 0 := by
            simpa using eq_neg_of_add_eq_zero_left hTmp0
          simpa using congrArg Neg.neg hEq.symm
        exact neg_injective this
      exact hTmp'
    have hB21formula :
        B 2 1 = A 2 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) := by
      have hEq20 : -A 2 2 * C 0 0 - B 2 1 * D 0 0 + B 2 0 * D 0 1 = 0 := by
        have hij := congrArg (fun M => M 2 0) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hA01, hA02, hC01, hC02] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hEq21 : -A 2 2 * C 1 0 - B 2 1 * D 1 0 + B 2 0 * D 1 1 = 0 := by
        have hij := congrArg (fun M => M 2 1) h₁tr
        simp [T, U, N3PureSingular.ω13, N4.J, Matrix.mul_apply, Matrix.transpose_apply,
          Fin.sum_univ_three, Fin.sum_univ_two, hC11, hC12] at hij
        simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm] using hij
      have hTmp :
          -A 2 2 * C 0 0 * D 1 1 + A 2 2 * C 1 0 * D 0 1 +
            B 2 1 * (D 1 0 * D 0 1 - D 0 0 * D 1 1) = 0 := by
        calc
          -A 2 2 * C 0 0 * D 1 1 + A 2 2 * C 1 0 * D 0 1 +
              B 2 1 * (D 1 0 * D 0 1 - D 0 0 * D 1 1) =
            (-A 2 2 * C 0 0 - B 2 1 * D 0 0 + B 2 0 * D 0 1) * D 1 1 -
              (-A 2 2 * C 1 0 - B 2 1 * D 1 0 + B 2 0 * D 1 1) * D 0 1 := by ring
          _ = 0 := by rw [hEq20, hEq21]; ring
      have hCoeff :
          D 1 0 * D 0 1 - D 0 0 * D 1 1 = -1 := by
        calc
          D 1 0 * D 0 1 - D 0 0 * D 1 1 = -(D 0 0 * D 1 1 - D 0 1 * D 1 0) := by ring
          _ = -1 := by rw [hDdet0]
      have hTmp' :
          B 2 1 = A 2 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) := by
        have : -B 2 1 = -(A 2 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1)) := by
          rw [hCoeff] at hTmp
          have hTmp0 : A 2 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) + -B 2 1 = 0 := by
            calc
              A 2 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) + -B 2 1 =
                  A 2 2 * C 1 0 * D 0 1 + (-B 2 1 + -(A 2 2 * C 0 0 * D 1 1)) := by ring
              _ = 0 := by
                simpa [sub_eq_add_neg, add_assoc, add_left_comm, add_comm, mul_assoc,
                  mul_left_comm, mul_comm] using hTmp
          have hEq : A 2 2 * (C 1 0 * D 0 1 - C 0 0 * D 1 1) = B 2 1 := by
            simpa using eq_neg_of_add_eq_zero_left hTmp0
          simpa using congrArg Neg.neg hEq.symm
        exact neg_injective this
      exact hTmp'
    have hBminor : B 1 1 * B 2 0 - B 1 0 * B 2 1 = 0 := by
      rw [hB10formula, hB11formula, hB20formula, hB21formula]
      ring
    have hBJBT : B * N4.J * Bᵀ = (0 : Matrix W W k) := by
      ext i j
      fin_cases i <;> fin_cases j
      · simp [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
          Fin.sum_univ_three, hB00, hB01]
      · simp [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
          Fin.sum_univ_three, hB00, hB01]
      · simp [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
          Fin.sum_univ_three, hB00, hB01]
      · simp [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
          Fin.sum_univ_three, hB00, hB01]
      ·
          simpa [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
            Fin.sum_univ_three, hB00, hB01, sub_eq_add_neg, add_assoc,
            add_left_comm, add_comm, mul_assoc, mul_left_comm, mul_comm] using
            congrArg Neg.neg hBminor
      ·
          simpa [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
            Fin.sum_univ_three, hB00, hB01, sub_eq_add_neg, add_assoc,
            add_left_comm, add_comm, mul_assoc, mul_left_comm, mul_comm] using
            congrArg Neg.neg hBminor
      ·
          simpa [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
            Fin.sum_univ_three, hB00, hB01, sub_eq_add_neg, add_assoc,
            add_left_comm, add_comm, mul_assoc, mul_left_comm, mul_comm]
      ·
          simpa [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
            Fin.sum_univ_three, hB00, hB01, sub_eq_add_neg, add_assoc,
            add_left_comm, add_comm, mul_assoc, mul_left_comm, mul_comm] using hBminor
      ·
          simpa [N4.J, Matrix.mul_apply, Matrix.transpose_apply, Fin.sum_univ_two,
            Fin.sum_univ_three, hB00, hB01, sub_eq_add_neg, add_assoc,
            add_left_comm, add_comm, mul_assoc, mul_left_comm, mul_comm]
    have hA13eq :
        N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) A =
          N3PureSingular.ω13 (k := k) := by
      have hTopLeft :
          A * N3PureSingular.ω13 (k := k) * Aᵀ + B * N4.J * Bᵀ =
            N3PureSingular.ω13 (k := k) := by
        ext i j
        have hij := congrArg (fun M => M (Sum.inl i) (Sum.inl j)) h₁
        simpa [FixesBivector, rep₁, Matrix.fromBlocks_transpose,
          Matrix.fromBlocks_multiply] using hij
      ext i j
      have hij := congrArg (fun M => M i j) hTopLeft
      simpa [N3PureSingular.ActBivector, hBJBT, Matrix.add_apply] using hij
    have hAfix :
        N3PureSingular.FixesPairBivector A := by
      exact ⟨by simpa [N3PureSingular.FixesBivector] using hA13eq,
        by simpa [N3PureSingular.FixesBivector] using hA23eq⟩
    rcases
      (N3PureSingular.fixesPairBivector_iff_shape (k := k) (g := A)).1 hAfix with
      ⟨hAshape0, hAshape⟩
    have hfix' :
        FixesPairBivector
          (Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) (A 0 0) (A 2 0) (A 2 1))
            B
            C
            D) := by
      have h₁' := h₁
      have h₂' := h₂
      rw [hAshape] at h₁' h₂'
      exact ⟨h₁', h₂'⟩
    rcases
      (fixesPair_fromBlocks_pureShape_iff
        (k := k) (a := A 0 0) (x := A 2 0) (y := A 2 1) hAshape0
        (B := B) (C := C) (D := D)).1 hfix' with
      ⟨hDdet', hshape'⟩
    refine ⟨hAshape0, hDdet', ?_⟩
    rw [hAshape]
    exact hshape'
  · rintro ⟨hA00, hDdet, hshape⟩
    rw [hshape]
    exact
      pointwise_shape_family
        (k := k)
        (a := A 0 0)
        (x := A 2 0)
        (y := A 2 1)
        (u := B 2 0)
        (v := B 2 1)
        hA00
        (E := D)
        hDdet

/-- A block diagonal family realising the standard Borel action on the pair. -/
theorem borel_lift_action
    (a b : k)
    (ha : a ≠ 0) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix I I k := !![a, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks G 0 0 E
    ActBivector rep₁ g = a • rep₁ + b • rep₂ ∧
      ActBivector rep₂ g = (a * a) • rep₂ := by
  intro G E g
  have hact := act_blockDiagonal (k := k) G E
  constructor
  ·
    rw [hact.1]
    have htop :
        N3PureSingular.ActBivector (N3PureSingular.ω13 (k := k)) G =
          a • (N3PureSingular.ω13 (k := k)) + b • (N3PureSingular.ω23 (k := k)) := by
      simpa [G] using (N3PureSingular.GL2_lift_action (k := k) a b 0 (a * a)).1
    have hbot : E * N4.J * Eᵀ = a • N4.J := by
      ext i j
      fin_cases i <;> fin_cases j <;>
        simp [E, N4.J, Matrix.mul_apply, Fin.sum_univ_two, Matrix.smul_apply]
        <;> ring_nf
    rw [rep₁, rep₂]
    ext i j
    cases i <;> cases j <;>
      simp [htop, hbot, Matrix.fromBlocks, Matrix.add_apply, Matrix.smul_apply]
  ·
    rw [hact.2]
    have htop :
        N3PureSingular.ActBivector (N3PureSingular.ω23 (k := k)) G =
          (a * a) • (N3PureSingular.ω23 (k := k)) := by
      simpa [G] using (N3PureSingular.GL2_lift_action (k := k) a b 0 (a * a)).2
    rw [rep₂]
    ext i j
    cases i <;> cases j <;>
      simp [htop, Matrix.fromBlocks, Matrix.smul_apply]

/-- The determinant of the standard Borel lift on the direct-sum `[a]` orbit is
`a^4`. -/
theorem borel_lift_det
    (a b : k) :
    let G : Matrix W W k := N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix I I k := !![a, 0; 0, 1]
    let g : Matrix V V k := Matrix.fromBlocks G 0 0 E
    Matrix.det g = (a * a) * (a * a) := by
  dsimp
  rw [Matrix.det_fromBlocks_zero₂₁, N3PureSingular.pureLift_det]
  simp [Matrix.det_fin_two]
  ring

end N5SimplePoint
end Wedge2Formalization
