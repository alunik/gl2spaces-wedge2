# `gl2spaces-wedge2`

This repository is the Magma companion package for the geometry-first paper on
the action of `GL(V)` on `Gr(2, \wedge^2 V)`.

Its purpose is to keep the orbit data and stabilizer data used in the paper in
one place, with no dependence on orbit computation at run time. Through
`n <= 7`, the orbit representatives and stabilizer descriptions are hard-coded
from the paper tables.

## What is included

- `src/paper_data.m`
  Hard-coded geometric and finite-field orbit lists through `n <= 7`, together
  with representative pairs and stabilizer descriptions matching the paper
  appendices.
- `src/local_matrix_models.m`
  Explicit matrix-generator models for the local families already available in
  verified closed form:
  `S_d` for `d = 2,3,4`, `S_2^2`, and `S_3 + J_{a,1}`.
- `src/stabilizer_generators.m`
  Hard-coded kernel-generator families, explicit quotient lifts on the
  associated `2`-spaces, and full stabilizer generators in `GL(V)` through
  `n <= 7`.
- `examples/print_geometric_data.m`
  Prints the geometric representatives and stabilizer descriptions.
- `examples/print_finite_data.m`
  Prints the finite-field representatives and stabilizer descriptions.
- `examples/print_finite_full_generators.m`
  Prints the actual stabilizer generators in `GL(V)` for a chosen sample field
  and dimension.
- `examples/check_counts_and_local_models.m`
  Checks the orbit counts and the currently implemented local matrix models on
  small fields.
- `examples/check_finite_stabilizer_specs.m`
  Verifies the hard-coded full finite-field stabilizer generators through
  `n <= 7`: full order checks on `q = 2,3` and quotient-lift image checks on
  `q = 5,7`.
- `examples/generate_geometric_orbit_tables_tex.m`
  Regenerates the appendix table of geometric orbit representatives.
- `examples/generate_geometric_stabilizer_tables_tex.m`
  Regenerates the appendix table of algebraic stabilizers over `k`.
- `examples/generate_finite_field_stabilizer_tables_tex.m`
  Regenerates the appendix table of finite-field stabilizers.
- `paper_tables/`
  Bundled copies of the generated appendix tables from the paper, together with
  the small paper snippets used by the audit script to verify that the repo
  outputs still match the manuscript.
- `tools/audit_paper_tables.py`
  Runs the repo-side table generators and checks that the bundled paper tables
  match them exactly.
- `tools/check_orbit_count_polynomials.m`
  Checks symbolically, for `n = 4,5,6,7`, that the sum of the orbit-stabilizer
  contributions from the hard-coded finite-field tables equals the Gaussian
  polynomial for the number of `2`-spaces in `\Lambda^2 V`.

## Scope

This package is intended to mirror the paper data exactly.

- The orbit representatives are generated directly from the block lists used in
  the paper.
- The stabilizer descriptions are hard-coded from the paper tables rather than
  computed by TameGenus or a subspace-transporter calculation.
- The full finite-field stabilizer-generator layer in `GL(V)` is explicit
  through `n <= 7`; nothing is computed on the fly from a transporter problem.

## Quick start

From the repository root:

```text
magma -b < examples/print_geometric_data.m
magma -b < examples/print_finite_data.m
magma -b < examples/print_finite_full_generators.m
magma -b < examples/check_counts_and_local_models.m
magma -b < examples/check_finite_stabilizer_specs.m
magma -b < tools/check_orbit_count_polynomials.m
python3 tools/audit_paper_tables.py
```

## Relation to the paper

The names, counts, and stabilizer strings returned by `src/paper_data.m` are
meant to match:

- `paper_tables/geometric_orbit_tables_generated.tex`
- `paper_tables/geometric_stabilizer_tables_generated.tex`
- `paper_tables/finite_field_stabilizer_tables_generated.tex`

This makes the repository a good place to keep the paper tables and the Magma
data layer in sync while also exposing actual full stabilizer generators in
`GL(V)` orbit-by-orbit through `n <= 7`.

## Accessing the stabilizers

For `F := GF(q)` and `4 <= n <= 7`,

```text
specs := FiniteStabilizerSpecsNLe7(F, n);
```

returns one record per finite-field orbit. Each record contains:

- `pair`
  the representative pair of alternating matrices
- `wedge1`, `wedge2`
  the wedge-coordinate representatives used in the paper tables
- `kernel`, `quotient`
  the paper-facing stabilizer strings
- `kernel_generators`
  generators for `K_L`
- `quotient_lifts`
  explicit lifts of quotient generators from `GL_2(q)` back into `GL(V)`
- `full_generators`
  generators for the full stabilizer subgroup `Stab_{GL(V)}(L)`
