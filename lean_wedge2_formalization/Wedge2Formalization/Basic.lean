/-!
This project explores a Lean formalization of the stabilizer calculations for
2-dimensional subspaces of `∧² V`.

The intended strategy is deliberately elementary. We do not formalize the full
weak-congruence classification. Instead we start from an explicit candidate list
of representatives and then prove the stabilizer formulas directly by matrix
calculations.

The first concrete experiment is the `n = 4` split-support representative from
the paper. It is treated in `Wedge2Formalization.N4`.
-/
