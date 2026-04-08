import Mathlib.Data.Matrix.Basic

namespace Wedge2Formalization
namespace PaperCore

open Matrix

variable {V : Type*} [Fintype V] [DecidableEq V]
variable {k : Type*} [Field k]

/-- The literal pointwise stabilizer attached to a fixing predicate. -/
def PointwiseStabilizer
    (Fixes : Matrix V V k → Prop) : Set (Matrix V V k) :=
  { g | Fixes g }

/-- A named family is exact for a fixing predicate if it captures precisely the
pointwise stabilizer. -/
def ExactFamily
    (Fixes : Matrix V V k → Prop)
    (K : Set (Matrix V V k)) : Prop :=
  ∀ g, Fixes g ↔ g ∈ K

/-- The induced action on the ordered basis `(Ω₁, Ω₂)` is given by the coefficient
matrix `M`. This is the paper-facing notion used for quotient-action statements. -/
def ActsOnOrderedPair
    (Act : Matrix V V k → Matrix V V k → Matrix V V k)
    (Ω₁ Ω₂ : Matrix V V k)
    (g : Matrix V V k)
    (M : Matrix (Fin 2) (Fin 2) k) : Prop :=
  Act Ω₁ g = M 0 0 • Ω₁ + M 0 1 • Ω₂ ∧
    Act Ω₂ g = M 1 0 • Ω₁ + M 1 1 • Ω₂

/-- Public setwise stabilizer notion for the span of an ordered pair. -/
def PreservesSpanPair
    (Act : Matrix V V k → Matrix V V k → Matrix V V k)
    (Ω₁ Ω₂ : Matrix V V k)
    (g : Matrix V V k) : Prop :=
  ∃ M : Matrix (Fin 2) (Fin 2) k, ActsOnOrderedPair Act Ω₁ Ω₂ g M

end PaperCore
end Wedge2Formalization
