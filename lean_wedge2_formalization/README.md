# Wedge2Formalization

Lean formalization of the geometric stabilizer calculations for `2`-planes in
`∧² V`, organized to match the geometric tables and Appendix A crosswalk of the
paper.

## Paper-facing entry points

For referee-style browsing, start in the public `Paper/` layer:

- `Wedge2Formalization/Paper/N4.lean`
- `Wedge2Formalization/Paper/N5.lean`
- `Wedge2Formalization/Paper/N6.lean`
- `Wedge2Formalization/Paper/N7.lean`

Every Appendix A row also has a dedicated file-level entry point under
`Wedge2Formalization/Paper/N*/Row*.lean`. For `n = 6` and part of `n = 7`
these are the main implementation files; for `n = 4`, `n = 5`, and the
remaining `n = 7` rows they are thin public wrappers designed for direct
referee navigation and doc browsing.

## Public API

Each public row namespace is intended to match one Appendix A table row and to
expose the same small reviewer-facing interface:

- `rep₁`, `rep₂`
- `paperRep₁`, `paperRep₂`, `paperChange`
- `paperRep₁_transport`, `paperRep₂_transport`
- `K`, `Qproj`
- `pointwise_stabilizer`
- `mem_K_iff`
- `mem_K_table_iff`
- `quotient_action` or the relevant quotient-action theorem
- `quotient_image`

When a divisor or rank-drop calculation is part of the argument, the public row
also exposes `det_zero_iff`.

The intention is that a reader can pick a table line from Appendix A and check
the corresponding Lean namespace without reading the internal construction
files first.

## Quick referee check

To verify one row from the paper:

1. Find the namespace for that row in the map below, or open the matching file
   `Wedge2Formalization/Paper/Nn/Rowm.lean`.
2. Read `paperRep₁`, `paperRep₂`, `paperRep₁_transport`, and
   `paperRep₂_transport` to identify the literal paper representative with the
   internal working model.
3. Read `mem_K_table_iff` for the kernel statement and `quotient_image` for the
   quotient statement.

This is the public interface used by the paper crosswalk.

Rendered documentation is built from the same repository by the GitHub Actions
workflow in
[`lean_wedge2_formalization/.github/workflows/lean_action_ci.yml`](/Users/aluna/Skew%20Symmetric/lean_wedge2_formalization/.github/workflows/lean_action_ci.yml),
so the Lean and Magma companion materials stay in one place.

## Appendix A namespace map

### `n = 4`

- row 1 (`J_{a,2}`): `Wedge2Formalization.Paper.N4.Row1`
- row 2 (`J_{a,1} + J_{b,1}`): `Wedge2Formalization.Paper.N4.Row2`
- row 3 (`S_1 + S_2`): `Wedge2Formalization.Paper.N4.Row3`

### `n = 5`

- row 1 (`S_3`): `Wedge2Formalization.Paper.N5.Row1`
- row 2 (`S_2 + S_1^2`): `Wedge2Formalization.Paper.N5.Row2`
- row 3 (`S_2 + J_{a,1}`): `Wedge2Formalization.Paper.N5.Row3`
- row 4 (`S_1 + J_{a,2}`): `Wedge2Formalization.Paper.N5.Row4`
- row 5 (`S_1 + J_{a,1} + J_{b,1}`): `Wedge2Formalization.Paper.N5.Row5`

### `n = 6`

- row 1 (`J_{a,3}`): `Wedge2Formalization.Paper.N6.Row1`
- row 2 (`J_{a,2} + J_{a,1}`): `Wedge2Formalization.Paper.N6.Row2`
- row 3 (`J_{a,2} + J_{b,1}`): `Wedge2Formalization.Paper.N6.Row3`
- row 4 (`2J_{a,1} + J_{b,1}`): `Wedge2Formalization.Paper.N6.Row4`
- row 5 (`J_{a,1} + J_{b,1} + J_{c,1}`): `Wedge2Formalization.Paper.N6.Row5`
- row 6 (`S_1^2 + J_{a,2}`): `Wedge2Formalization.Paper.N6.Row6`
- row 7 (`S_1^2 + J_{a,1} + J_{b,1}`): `Wedge2Formalization.Paper.N6.Row7`
- row 8 (`S_1 + S_2 + J_{a,1}`): `Wedge2Formalization.Paper.N6.Row8`
- row 9 (`S_1^3 + S_2`): `Wedge2Formalization.Paper.N6.Row9`
- row 10 (`S_1 + S_3`): `Wedge2Formalization.Paper.N6.Row10`
- row 11 (`S_2^2`): `Wedge2Formalization.Paper.N6.Row11`

### `n = 7`

- row 1 (`S_4`): `Wedge2Formalization.Paper.N7.Row1`
- row 2 (`S_3 + S_1^2`): `Wedge2Formalization.Paper.N7.Row2`
- row 3 (`S_2^2 + S_1`): `Wedge2Formalization.Paper.N7.Row3`
- row 4 (`S_2 + S_1^4`): `Wedge2Formalization.Paper.N7.Row4`
- row 5 (`S_3 + J_{a,1}`): `Wedge2Formalization.Paper.N7.Row5`
- row 6 (`S_2 + S_1^2 + J_{a,1}`): `Wedge2Formalization.Paper.N7.Row6`
- row 7 (`S_2 + J_{a,2}`): `Wedge2Formalization.Paper.N7.Row7`
- row 8 (`S_2 + 2J_{a,1}`): `Wedge2Formalization.Paper.N7.Row8`
- row 9 (`S_2 + J_{a,1} + J_{b,1}`): `Wedge2Formalization.Paper.N7.Row9`
- row 10 (`S_1^3 + J_{a,2}`): `Wedge2Formalization.Paper.N7.Row10`
- row 11 (`S_1^3 + J_{a,1} + J_{b,1}`): `Wedge2Formalization.Paper.N7.Row11`
- row 12 (`S_1 + J_{a,3}`): `Wedge2Formalization.Paper.N7.Row12`
- row 13 (`S_1 + J_{a,2} + J_{a,1}`): `Wedge2Formalization.Paper.N7.Row13`
- row 14 (`S_1 + J_{a,2} + J_{b,1}`): `Wedge2Formalization.Paper.N7.Row14`
- row 15 (`S_1 + 2J_{a,1} + J_{b,1}`): `Wedge2Formalization.Paper.N7.Row15`
- row 16 (`S_1 + J_{a,1} + J_{b,1} + J_{c,1}`): `Wedge2Formalization.Paper.N7.Row16`

## Build

The whole library builds with:

```sh
lake build Wedge2Formalization
```

The most useful public targets are:

```sh
lake build Wedge2Formalization.Paper.N4
lake build Wedge2Formalization.Paper.N5
lake build Wedge2Formalization.Paper.N6
lake build Wedge2Formalization.Paper.N7
```

## Internal layer

The lower-level files in `Wedge2Formalization/` keep the block-matrix,
coordinate, and determinant computations. The `Paper/` layer is only a
translation layer from those internal proofs to the language of Appendix A.

The formalization deliberately stays concrete: stabilizers and quotients are
presented as explicit matrix families rather than abstract algebraic-group
objects.
