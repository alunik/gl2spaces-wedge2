import Wedge2Formalization.N5SimplePoint
import Mathlib.LinearAlgebra.Matrix.Block
import Mathlib.LinearAlgebra.Matrix.SchurComplement

open Matrix

namespace Wedge2Formalization
namespace N6SimplePoint

variable {k : Type*} [Field k]

/-- The one-dimensional radical factor. -/
abbrev I := Fin 1

/-- The `5`-dimensional quotient carrying the direct-sum `[a]` orbit. -/
abbrev W := N5SimplePoint.V

/-- The ambient `6`-dimensional space. -/
abbrev V := I ⊕ W

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

/-- The `1 × 1` scalar block in the radical direction. -/
def scalarBlock (a : k) : Matrix I I k := !![a]

/-- The embedded direct-sum `[a]` representative in dimension `6`. -/
def radSimpleRep₁ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N5SimplePoint.rep₁ (k := k))

/-- The second basis vector of the embedded direct-sum `[a]` representative. -/
def radSimpleRep₂ : Matrix V V k :=
  Matrix.fromBlocks 0 0 0 (N5SimplePoint.rep₂ (k := k))

/-- Exact stabilizer of the pair. -/
def FixesRadSimplePairBivector (g : Matrix V V k) : Prop :=
  FixesBivector radSimpleRep₁ g ∧ FixesBivector radSimpleRep₂ g

/-- Setwise stabilizer of the `2`-space in the chosen basis. -/
def PreservesRadSimpleSubspaceBivector (g : Matrix V V k) : Prop :=
  ∃ α β γ δ : k,
    ActBivector radSimpleRep₁ g = α • radSimpleRep₁ + β • radSimpleRep₂ ∧
      ActBivector radSimpleRep₂ g = γ • radSimpleRep₁ + δ • radSimpleRep₂

/-- Acting on an embedded bivector with zero upper-right block only changes the
lower-right block. -/
theorem act_embedded_fromBlocks_zeroUpperRight
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
      Matrix.fromBlocks 0 0 0 (N5SimplePoint.ActBivector Ω D) := by
  ext i j
  cases i <;> cases j <;>
    simp [ActBivector, N5SimplePoint.ActBivector, Matrix.fromBlocks_transpose,
      Matrix.fromBlocks_multiply]

/-- Pointwise fixing of an embedded bivector reduces to the lower-right block when the
upper-right block is zero. -/
theorem fixes_embedded_fromBlocks_zeroUpperRight_iff
    (Ω : Matrix W W k)
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) ↔
      N5SimplePoint.FixesBivector Ω D := by
  constructor
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω at h
    rw [act_embedded_fromBlocks_zeroUpperRight] at h
    have h' := congrArg Matrix.toBlocks₂₂ h
    simpa [N5SimplePoint.FixesBivector, N5SimplePoint.ActBivector] using h'
  · intro h
    change
      ActBivector (Matrix.fromBlocks 0 0 0 Ω) (Matrix.fromBlocks A 0 C D) =
        Matrix.fromBlocks 0 0 0 Ω
    rw [act_embedded_fromBlocks_zeroUpperRight]
    simpa [N5SimplePoint.FixesBivector, N5SimplePoint.ActBivector] using h

/-- Pairwise lower-right reduction for the radical extension of the `[a]` family. -/
theorem fixesRadSimplePair_fromBlocks_zeroUpperRight_iff
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadSimplePairBivector (Matrix.fromBlocks A 0 C D) ↔
      N5SimplePoint.FixesPairBivector D := by
  constructor
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5SimplePoint.rep₁ (k := k)) (A := A) (C := C) (D := D)).1 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5SimplePoint.rep₂ (k := k)) (A := A) (C := C) (D := D)).1 h₂⟩
  · rintro ⟨h₁, h₂⟩
    exact
      ⟨(fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5SimplePoint.rep₁ (k := k)) (A := A) (C := C) (D := D)).2 h₁,
        (fixes_embedded_fromBlocks_zeroUpperRight_iff
          (Ω := N5SimplePoint.rep₂ (k := k)) (A := A) (C := C) (D := D)).2 h₂⟩

/-- Even without assuming the upper-right block vanishes, pointwise fixing forces the
lower-right block to fix the underlying `5`-dimensional direct-sum `[a]` pair. -/
theorem fixesRadSimplePair_lowerRight
    (A : Matrix I I k)
    (B : Matrix I W k)
    (C : Matrix W I k)
    (D : Matrix W W k) :
    FixesRadSimplePairBivector (Matrix.fromBlocks A B C D) →
      N5SimplePoint.FixesPairBivector D := by
  rintro ⟨h₁, h₂⟩
  refine ⟨?_, ?_⟩
  ·
    ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₁
    simpa [FixesBivector, radSimpleRep₁, N5SimplePoint.FixesBivector,
      N5SimplePoint.ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
  ·
    ext i j
    have hij := congrArg (fun M => M (Sum.inr i) (Sum.inr j)) h₂
    simpa [FixesBivector, radSimpleRep₂, N5SimplePoint.FixesBivector,
      N5SimplePoint.ActBivector, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij

/-- Setwise preservation also reduces to the lower-right block when the upper-right block
vanishes. -/
theorem preservesRadSimpleSubspace_fromBlocks_zeroUpperRight
    (A : Matrix I I k)
    (C : Matrix W I k)
    (D : Matrix W W k)
    (hD : N5SimplePoint.PreservesSubspaceBivector D) :
    PreservesRadSimpleSubspaceBivector (Matrix.fromBlocks A 0 C D) := by
  rcases hD with ⟨α, β, γ, δ, h₁, h₂⟩
  refine ⟨α, β, γ, δ, ?_, ?_⟩
  ·
    calc
      ActBivector radSimpleRep₁ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N5SimplePoint.ActBivector (N5SimplePoint.rep₁ (k := k)) D) := by
              simpa [radSimpleRep₁] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N5SimplePoint.rep₁ (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (α • N5SimplePoint.rep₁ + β • N5SimplePoint.rep₂) := by simp [h₁]
      _ = α • radSimpleRep₁ + β • radSimpleRep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [radSimpleRep₁, radSimpleRep₂, Matrix.fromBlocks, Matrix.add_apply]
  ·
    calc
      ActBivector radSimpleRep₂ (Matrix.fromBlocks A 0 C D)
          = Matrix.fromBlocks 0 0 0 (N5SimplePoint.ActBivector (N5SimplePoint.rep₂ (k := k)) D) := by
              simpa [radSimpleRep₂] using
                act_embedded_fromBlocks_zeroUpperRight
                  (Ω := N5SimplePoint.rep₂ (k := k)) (A := A) (C := C) (D := D)
      _ = Matrix.fromBlocks 0 0 0 (γ • N5SimplePoint.rep₁ + δ • N5SimplePoint.rep₂) := by simp [h₂]
      _ = γ • radSimpleRep₁ + δ • radSimpleRep₂ := by
            ext i j
            cases i <;> cases j <;>
              simp [radSimpleRep₁, radSimpleRep₂, Matrix.fromBlocks, Matrix.add_apply]

/-- The obvious pointwise family on the radical extension of the `[a]` orbit. -/
theorem radSimple_pointwise_family
    (u a : k)
    (ha : a ≠ 0)
    (C : Matrix W I k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a)
          0
          0
          E)) := by
  have hD :
      N5SimplePoint.FixesPairBivector
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a)
          0
          0
          E) := by
    exact N5SimplePoint.pointwise_levi_family (k := k) a ha E hE
  exact
    (fixesRadSimplePair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := scalarBlock (k := k) u)
      (C := C)
      (D := Matrix.fromBlocks
        (N3PureSingular.pointwiseScale (k := k) a)
        0
        0
        E)).2 hD

/-- The basic `4`-parameter unipotent family on the `5`-dimensional `[a]` quotient
also extends through the radical direction. -/
theorem radSimple_pointwise_unipotent_family
    (u : k)
    (C : Matrix W I k)
    (x y p q : k) :
    FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)) := by
  have hD :
      N5SimplePoint.FixesPairBivector
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q) := by
    exact N5SimplePoint.pointwise_unipotent_family (k := k) x y p q
  exact
    (fixesRadSimplePair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := scalarBlock (k := k) u)
      (C := C)
      (D := N5SimplePoint.pointwiseUnipotent (k := k) x y p q)).2 hD

/-- Combining the explicit pointwise unipotent family on the `5`-dimensional quotient
with the Levi family gives a concrete pointwise subgroup on the direct-sum `[a]`
orbit in dimension `6`. -/
theorem radSimple_pointwise_product_family
    (u : k)
    (C : Matrix W I k)
    (x y p q a : k)
    (ha : a ≠ 0)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    let gU : Matrix V V k :=
      Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)
    let gL : Matrix V V k :=
      Matrix.fromBlocks
        (scalarBlock (k := k) (1 : k))
        0
        0
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a)
          0
          0
          E)
    FixesRadSimplePairBivector (gU * gL) := by
  intro gU gL
  have hU := radSimple_pointwise_unipotent_family (k := k) (u := u) (C := C)
    (x := x) (y := y) (p := p) (q := q)
  have hL :=
    radSimple_pointwise_family
      (k := k)
      (u := (1 : k))
      (a := a)
      ha
      (C := 0)
      (E := E)
      hE
  constructor <;> rw [FixesBivector]
  ·
    calc
      ActBivector radSimpleRep₁ (gU * gL) = ActBivector (ActBivector radSimpleRep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector radSimpleRep₁ gU := by
        simpa [FixesBivector] using congrArg (fun M => ActBivector M gU) hL.1
      _ = radSimpleRep₁ := by simpa [FixesBivector] using hU.1
  ·
    calc
      ActBivector radSimpleRep₂ (gU * gL) = ActBivector (ActBivector radSimpleRep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector radSimpleRep₂ gU := by
        simpa [FixesBivector] using congrArg (fun M => ActBivector M gU) hL.2
      _ = radSimpleRep₂ := by simpa [FixesBivector] using hU.2

/-- Every block matrix built from the expected `n = 5` pointwise shape on the quotient
and an arbitrary radical column fixes the direct-sum `[a]` pair in dimension `6`
pointwise. -/
theorem radSimple_pointwise_shape_family
    (u : k)
    (C : Matrix W I k)
    (a x y p q : k)
    (ha : a ≠ 0)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE : E.det = 1) :
    let B : Matrix N5SimplePoint.W N5SimplePoint.I k := !![(0 : k), 0; 0, 0; p, q]
    let C0 : Matrix N5SimplePoint.I N5SimplePoint.W k :=
      !![a * (p * E 0 1 - q * E 0 0), 0, 0;
         a * (p * E 1 1 - q * E 1 0), 0, 0]
    FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E)) := by
  intro B C0
  have hD :
      N5SimplePoint.FixesPairBivector
        (Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E) := by
    exact
      N5SimplePoint.pointwise_shape_family
        (k := k)
        (a := a)
        (x := x)
        (y := y)
        (u := p)
        (v := q)
        ha
        (E := E)
        hE
  exact
    (fixesRadSimplePair_fromBlocks_zeroUpperRight_iff
      (k := k)
      (A := scalarBlock (k := k) u)
      (C := C)
      (D := Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        B
        C0
        E)).2 hD

/-- Once the quotient block is already in the pure singular shape, pointwise fixing
forces the remaining quotient blocks to have the expected direct-sum `[a]` form, while
the radical scalar and radical column stay free. -/
theorem fixesRadSimplePair_fromBlocks_nested_pureShape_iff
    (u : k)
    (C : Matrix W I k)
    (a x y : k)
    (B : Matrix N5SimplePoint.W N5SimplePoint.I k)
    (C0 : Matrix N5SimplePoint.I N5SimplePoint.W k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (ha : a ≠ 0) :
    FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E)) ↔
      E.det = 1 ∧
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E =
        Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          (!![(0 : k), 0; 0, 0; B 2 0, B 2 1])
          (!![a * (B 2 0 * E 0 1 - B 2 1 * E 0 0), 0, 0;
             a * (B 2 0 * E 1 1 - B 2 1 * E 1 0), 0, 0])
          E := by
  constructor
  · intro h
    have hD :
        N5SimplePoint.FixesPairBivector
          (Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            B
            C0
            E) := by
      exact
        (fixesRadSimplePair_fromBlocks_zeroUpperRight_iff
          (k := k)
          (A := scalarBlock (k := k) u)
          (C := C)
          (D := Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            B
            C0
            E)).1 h
    exact
      (N5SimplePoint.fixesPair_fromBlocks_pureShape_iff
        (k := k) (a := a) (x := x) (y := y) ha (B := B) (C := C0) (D := E)).1 hD
  · intro h
    exact
      (fixesRadSimplePair_fromBlocks_zeroUpperRight_iff
        (k := k)
        (A := scalarBlock (k := k) u)
        (C := C)
        (D := Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C0
          E)).2 <|
        (N5SimplePoint.fixesPair_fromBlocks_pureShape_iff
          (k := k) (a := a) (x := x) (y := y) ha (B := B) (C := C0) (D := E)).2 h

/-- The explicit nested pointwise shape on the `5`-dimensional quotient has determinant
`a`. -/
theorem nestedPointwiseShape_det
    (a x y u v : k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (ha : a ≠ 0)
    (hE : E.det = 1) :
    Matrix.det
      (Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        (!![(0 : k), 0; 0, 0; u, v])
        (!![a * (u * E 0 1 - v * E 0 0), 0, 0;
           a * (u * E 1 1 - v * E 1 0), 0, 0])
        E) = a := by
  let A : Matrix N5SimplePoint.W N5SimplePoint.W k :=
    N3PureSingular.pureSingularShape (k := k) a x y
  let B : Matrix N5SimplePoint.W N5SimplePoint.I k := !![(0 : k), 0; 0, 0; u, v]
  let C : Matrix N5SimplePoint.I N5SimplePoint.W k :=
    !![a * (u * E 0 1 - v * E 0 0), 0, 0;
       a * (u * E 1 1 - v * E 1 0), 0, 0]
  let Ainv : Matrix N5SimplePoint.W N5SimplePoint.W k :=
    !![a⁻¹, 0, 0;
       0, a⁻¹, 0;
       -x, -y, a]
  have hAinv : A * Ainv = 1 := by
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [A, Ainv, N3PureSingular.pureSingularShape, Matrix.mul_apply, Fin.sum_univ_three, ha] <;>
      field_simp [ha] <;>
      ring
  letI : Invertible A := invertibleOfRightInverse A Ainv hAinv
  have hCAB : C * (⅟A) * B = 0 := by
    change C * Ainv * B = 0
    ext i j
    fin_cases i <;> fin_cases j <;>
      simp [Ainv, B, C, Matrix.mul_apply, Fin.sum_univ_two, Fin.sum_univ_three] <;>
      ring
  have hdetA : Matrix.det A = a := by
    simpa [A, N3PureSingular.pureSingularShape, Matrix.det_fin_three, ha]
  calc
    Matrix.det (Matrix.fromBlocks A B C E) = Matrix.det A * Matrix.det (E - C * (⅟A) * B) := by
      rw [Matrix.det_fromBlocks₁₁ A B C E]
    _ = Matrix.det A * Matrix.det E := by rw [hCAB, sub_zero]
    _ = a := by rw [hdetA, hE]; ring

/-- For the direct-sum `[a]` row in dimension `6`, pointwise fixing forces the
upper-right radical row to vanish, and then the lower-right block has exactly the
expected nested shape coming from the `5`-dimensional quotient. -/
theorem fixesRadSimplePair_fromBlocks_nested_iff_shape
    (u : k)
    (B0 : Matrix I W k)
    (C0 : Matrix W I k)
    (a x y : k)
    (B : Matrix N5SimplePoint.W N5SimplePoint.I k)
    (C : Matrix N5SimplePoint.I N5SimplePoint.W k)
    (E : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (ha : a ≠ 0) :
    FixesRadSimplePairBivector
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        B0
        C0
        (Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C
          E)) ↔
      B0 = 0 ∧
        E.det = 1 ∧
          Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            B
            C
            E =
          Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            (!![(0 : k), 0; 0, 0; B 2 0, B 2 1])
            (!![a * (B 2 0 * E 0 1 - B 2 1 * E 0 0), 0, 0;
               a * (B 2 0 * E 1 1 - B 2 1 * E 1 0), 0, 0])
            E := by
  constructor
  · intro h
    have hD :
        N5SimplePoint.FixesPairBivector
          (Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            B
            C
            E) :=
      fixesRadSimplePair_lowerRight
        (k := k)
        (A := scalarBlock (k := k) u)
        (B := B0)
        (C := C0)
        (D := Matrix.fromBlocks
          (N3PureSingular.pureSingularShape (k := k) a x y)
          B
          C
          E)
        h
    rcases
      (N5SimplePoint.fixesPair_fromBlocks_pureShape_iff
        (k := k) (a := a) (x := x) (y := y) ha (B := B) (C := C) (D := E)).1 hD with
      ⟨hE, hshape⟩
    rcases h with ⟨h₁, h₂⟩
    have hDdet :
        Matrix.det
          (Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            B
            C
            E) = a := by
      rw [hshape]
      simpa using
        nestedPointwiseShape_det
          (k := k)
          (a := a)
          (x := x)
          (y := y)
          (u := B 2 0)
          (v := B 2 1)
          (E := E)
          ha
          hE
    have hunitT :
        IsUnit
          ((Matrix.fromBlocks
            (N3PureSingular.pureSingularShape (k := k) a x y)
            B
            C
            E)ᵀ).det := by
      simpa [Matrix.det_transpose, hDdet] using (isUnit_iff_ne_zero.mpr ha)
    have h₂tr :
        B0 * N5SimplePoint.rep₂ (k := k) *
            (Matrix.fromBlocks
              (N3PureSingular.pureSingularShape (k := k) a x y)
              B
              C
              E)ᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₂
      simpa [FixesBivector, radSimpleRep₂, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have h₂zero : B0 * N5SimplePoint.rep₂ (k := k) = 0 := by
      have htmp := congrArg (fun M => M * ((Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        B
        C
        E)ᵀ)⁻¹) h₂tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hB02 : B0 0 (Sum.inl 2) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 1)) h₂zero
      have hneg : -B0 0 (Sum.inl 2) = 0 := by
        simpa [N5SimplePoint.rep₂, N3PureSingular.ω23, Matrix.mul_apply,
          Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two, N4.J] using hij
      exact neg_eq_zero.mp hneg
    have hB01 : B0 0 (Sum.inl 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 2)) h₂zero
      simpa [N5SimplePoint.rep₂, N3PureSingular.ω23, Matrix.mul_apply,
        Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two, N4.J, hB02] using hij
    have h₁tr :
        B0 * N5SimplePoint.rep₁ (k := k) *
            (Matrix.fromBlocks
              (N3PureSingular.pureSingularShape (k := k) a x y)
              B
              C
              E)ᵀ = 0 := by
      ext i j
      fin_cases i
      have hij := congrArg (fun M => M (Sum.inl 0) (Sum.inr j)) h₁
      simpa [FixesBivector, radSimpleRep₁, Matrix.fromBlocks_transpose, Matrix.fromBlocks_multiply] using hij
    have h₁zero : B0 * N5SimplePoint.rep₁ (k := k) = 0 := by
      have htmp := congrArg (fun M => M * ((Matrix.fromBlocks
        (N3PureSingular.pureSingularShape (k := k) a x y)
        B
        C
        E)ᵀ)⁻¹) h₁tr
      simpa [Matrix.mul_assoc, Matrix.mul_nonsing_inv _ hunitT] using htmp
    have hB00 : B0 0 (Sum.inl 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inl 2)) h₁zero
      simpa [N5SimplePoint.rep₁, N3PureSingular.ω13, Matrix.mul_apply,
        Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two, N4.J, hB02] using hij
    have hB04 : B0 0 (Sum.inr 1) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 0)) h₁zero
      have hneg : -B0 0 (Sum.inr 1) = 0 := by
        simpa [N5SimplePoint.rep₁, N3PureSingular.ω13, Matrix.mul_apply,
          Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two, N4.J, hB02] using hij
      exact neg_eq_zero.mp hneg
    have hB03 : B0 0 (Sum.inr 0) = 0 := by
      have hij := congrArg (fun M => M 0 (Sum.inr 1)) h₁zero
      simpa [N5SimplePoint.rep₁, N3PureSingular.ω13, Matrix.mul_apply,
        Fintype.sum_sum_type, Fin.sum_univ_three, Fin.sum_univ_two, N4.J, hB02, hB04] using hij
    refine ⟨?_, hE, hshape⟩
    ext i j
    fin_cases i
    cases j with
    | inl j =>
        fin_cases j <;> simp [hB00, hB01, hB02]
    | inr j =>
        fin_cases j <;> simp [hB03, hB04]
  · rintro ⟨hB0, hE, hshape⟩
    rw [hB0]
    exact
      (fixesRadSimplePair_fromBlocks_nested_pureShape_iff
        (k := k)
        (u := u)
        (C := C0)
        (a := a)
        (x := x)
        (y := y)
        (B := B)
        (C0 := C)
        (E := E)
        ha).2 ⟨hE, hshape⟩

/-- The standard Borel lift on the `[a]` quotient extends through the radical direction. -/
theorem radSimple_borel_lift_action
    (u a b : k)
    (ha : a ≠ 0)
    (C : Matrix W I k) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks G 0 0 E
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) u) 0 C h
    ActBivector radSimpleRep₁ g = a • radSimpleRep₁ + b • radSimpleRep₂ ∧
      ActBivector radSimpleRep₂ g = (a * a) • radSimpleRep₂ := by
  intro G E h g
  have hD := N5SimplePoint.borel_lift_action (k := k) (a := a) (b := b) ha
  constructor
  ·
    calc
      ActBivector radSimpleRep₁ g = Matrix.fromBlocks 0 0 0 (N5SimplePoint.ActBivector (N5SimplePoint.rep₁ (k := k)) h) := by
        simpa [g, h, radSimpleRep₁] using
          act_embedded_fromBlocks_zeroUpperRight
            (k := k)
            (Ω := N5SimplePoint.rep₁ (k := k))
            (A := scalarBlock (k := k) u)
            (C := C)
            (D := h)
      _ = Matrix.fromBlocks 0 0 0 (a • N5SimplePoint.rep₁ + b • N5SimplePoint.rep₂) := by
        simpa [G, E, h] using hD.1
      _ = a • radSimpleRep₁ + b • radSimpleRep₂ := by
        ext i j
        cases i <;> cases j <;>
          simp [radSimpleRep₁, radSimpleRep₂, Matrix.fromBlocks, Matrix.add_apply]
  ·
    calc
      ActBivector radSimpleRep₂ g = Matrix.fromBlocks 0 0 0 (N5SimplePoint.ActBivector (N5SimplePoint.rep₂ (k := k)) h) := by
        simpa [g, h, radSimpleRep₂] using
          act_embedded_fromBlocks_zeroUpperRight
            (k := k)
            (Ω := N5SimplePoint.rep₂ (k := k))
            (A := scalarBlock (k := k) u)
            (C := C)
            (D := h)
      _ = Matrix.fromBlocks 0 0 0 ((a * a) • N5SimplePoint.rep₂) := by
        simpa [G, E, h] using hD.2
      _ = (a * a) • radSimpleRep₂ := by
        ext i j
        cases i <;> cases j <;>
          simp [radSimpleRep₂, Matrix.fromBlocks, Matrix.smul_apply]

/-- The determinant of the standard Borel lift on the `n = 6` direct-sum `[a]`
orbit is the radical scalar times `a^4`. -/
theorem radSimple_borel_lift_det
    (u a b : k)
    (C : Matrix W I k) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks G 0 0 E
    let g : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) u) 0 C h
    Matrix.det g = u * ((a * a) * (a * a)) := by
  dsimp
  rw [Matrix.det_fromBlocks_zero₁₂, N5SimplePoint.borel_lift_det]
  simp [scalarBlock]

/-- The explicit pointwise subgroup on the direct-sum `[a]` orbit combines with the
standard Borel lift exactly as expected. -/
theorem radSimple_pointwise_borel_product_lift_action
    (u : k)
    (C : Matrix W I k)
    (x y p q a₀ : k)
    (ha₀ : a₀ ≠ 0)
    (E0 : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE0 : E0.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let gU : Matrix V V k :=
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)) *
      (Matrix.fromBlocks
        (scalarBlock (k := k) (1 : k))
        0
        0
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a₀)
          0
          0
          E0))
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks G 0 0 E
    let gL : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) 1) 0 0 h
    ActBivector radSimpleRep₁ (gU * gL) = a • radSimpleRep₁ + b • radSimpleRep₂ ∧
      ActBivector radSimpleRep₂ (gU * gL) = (a * a) • radSimpleRep₂ := by
  intro gU G E h gL
  have hU :=
    radSimple_pointwise_product_family
      (k := k)
      (u := u)
      (C := C)
      (x := x)
      (y := y)
      (p := p)
      (q := q)
      (a := a₀)
      ha₀
      (E := E0)
      hE0
  have hL := radSimple_borel_lift_action (k := k) (u := 1) (a := a) (b := b) ha (C := 0)
  have hUg1 : ActBivector radSimpleRep₁ gU = radSimpleRep₁ := by
    simpa [gU] using hU.1
  have hUg2 : ActBivector radSimpleRep₂ gU = radSimpleRep₂ := by
    simpa [gU] using hU.2
  constructor
  · calc
      ActBivector radSimpleRep₁ (gU * gL) = ActBivector (ActBivector radSimpleRep₁ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector (a • radSimpleRep₁ + b • radSimpleRep₂) gU := by
        simpa [gL, h, G, E] using congrArg (fun M => ActBivector M gU) hL.1
      _ = a • ActBivector radSimpleRep₁ gU + b • ActBivector radSimpleRep₂ gU := by
        simp [ActBivector, Matrix.mul_add, Matrix.add_mul, Matrix.smul_mul, Matrix.mul_smul]
      _ = a • radSimpleRep₁ + b • radSimpleRep₂ := by
        rw [hUg1, hUg2]
  · calc
      ActBivector radSimpleRep₂ (gU * gL) = ActBivector (ActBivector radSimpleRep₂ gL) gU := by
        rw [actBivector_mul]
      _ = ActBivector ((a * a) • radSimpleRep₂) gU := by
        simpa [gL, h, G, E] using congrArg (fun M => ActBivector M gU) hL.2
      _ = (a * a) • ActBivector radSimpleRep₂ gU := by
        simp [ActBivector, Matrix.smul_mul, Matrix.mul_smul]
      _ = (a * a) • radSimpleRep₂ := by
        rw [hUg2]

/-- Right-multiplying the standard Borel lift by the explicit pointwise subgroup on the
direct-sum `[a]` orbit does not change the quotient action. -/
theorem radSimple_borel_pointwise_right_product_lift_action
    (u : k)
    (C : Matrix W I k)
    (x y p q a₀ : k)
    (ha₀ : a₀ ≠ 0)
    (E0 : Matrix N5SimplePoint.I N5SimplePoint.I k)
    (hE0 : E0.det = 1)
    (a b : k)
    (ha : a ≠ 0) :
    let G : Matrix N5SimplePoint.W N5SimplePoint.W k :=
      N3PureSingular.pureLift (k := k) a b 0 (a * a) (1 : k)
    let E : Matrix N5SimplePoint.I N5SimplePoint.I k := !![a, 0; 0, 1]
    let h : Matrix W W k := Matrix.fromBlocks G 0 0 E
    let gL : Matrix V V k := Matrix.fromBlocks (scalarBlock (k := k) 1) 0 0 h
    let gU : Matrix V V k :=
      (Matrix.fromBlocks
        (scalarBlock (k := k) u)
        0
        C
        (N5SimplePoint.pointwiseUnipotent (k := k) x y p q)) *
      (Matrix.fromBlocks
        (scalarBlock (k := k) (1 : k))
        0
        0
        (Matrix.fromBlocks
          (N3PureSingular.pointwiseScale (k := k) a₀)
          0
          0
          E0))
    ActBivector radSimpleRep₁ (gL * gU) = a • radSimpleRep₁ + b • radSimpleRep₂ ∧
      ActBivector radSimpleRep₂ (gL * gU) = (a * a) • radSimpleRep₂ := by
  intro G E h gL gU
  have hL := radSimple_borel_lift_action (k := k) (u := 1) (a := a) (b := b) ha (C := 0)
  have hU :=
    radSimple_pointwise_product_family
      (k := k)
      (u := u)
      (C := C)
      (x := x)
      (y := y)
      (p := p)
      (q := q)
      (a := a₀)
      ha₀
      (E := E0)
      hE0
  have hUg1 : ActBivector radSimpleRep₁ gU = radSimpleRep₁ := by
    simpa [gU] using hU.1
  have hUg2 : ActBivector radSimpleRep₂ gU = radSimpleRep₂ := by
    simpa [gU] using hU.2
  constructor
  · calc
      ActBivector radSimpleRep₁ (gL * gU) = ActBivector (ActBivector radSimpleRep₁ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector radSimpleRep₁ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hUg1
      _ = a • radSimpleRep₁ + b • radSimpleRep₂ := by
        simpa [gL, h, G, E] using hL.1
  · calc
      ActBivector radSimpleRep₂ (gL * gU) = ActBivector (ActBivector radSimpleRep₂ gU) gL := by
        rw [actBivector_mul]
      _ = ActBivector radSimpleRep₂ gL := by
        simpa [gU] using congrArg (fun M => ActBivector M gL) hUg2
      _ = (a * a) • radSimpleRep₂ := by
        simpa [gL, h, G, E] using hL.2

end N6SimplePoint
end Wedge2Formalization
