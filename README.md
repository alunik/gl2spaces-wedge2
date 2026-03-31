# `gl2spaces-wedge2`

This repository is the Magma companion package for the geometry-first paper on
the action of `GL(V)` on `Gr(2, \wedge^2 V)`.

Through `n <= 7`, it provides:
- hard-coded geometric and finite-field orbit representatives;
- hard-coded stabilizer descriptions matching the paper;
- full hard-coded stabilizer generators in `GL(V)`.

Nothing here computes stabilizers by transporter search at runtime.

## Repository layout

- `src/`
  Core data and generator code.
- `examples/`
  The three user-facing entry points.
- `dev/`
  Verification, table-regeneration, and audit scripts for maintaining the
  package.
- `paper_tables/`
  Bundled TeX tables from the paper, used only by the maintenance audit.

## If You Just Want To Use It

From the repository root:

```text
magma -b < examples/geometric_orbits.m
magma -b < examples/finite_orbits.m
magma -b < examples/finite_stabilizers.m
```

These scripts are the intended public entry points.

- `examples/geometric_orbits.m`
  prints the geometric orbit representatives and stabilizer descriptions.
- `examples/finite_orbits.m`
  prints the finite-field orbit representatives and stabilizer descriptions.
- `examples/finite_stabilizers.m`
  prints the actual stabilizer generators in `GL(V)` for a chosen field and
  dimension.

Each script has the field and dimension choices at the top, so the easiest way
to use the package is simply to edit those parameters and rerun.

## Core API

The main entry points in `src/` are:

- `GeometricOrbitDataNLe7(F, n)` in `src/paper_data.m`
- `FiniteOrbitDataNLe7(F, n)` in `src/paper_data.m`
- `FiniteStabilizerSpecsNLe7(F, n)` in `src/stabilizer_generators.m`

For `F := GF(q)` and `4 <= n <= 7`,

```text
specs := FiniteStabilizerSpecsNLe7(F, n);
```

returns one record per finite-field orbit. Each record contains:

- `pair`
  the representative pair of alternating matrices;
- `wedge1`, `wedge2`
  the wedge-coordinate representatives used in the paper tables;
- `kernel`, `quotient`
  the paper-facing stabilizer strings;
- `kernel_generators`
  generators for `K_L`;
- `quotient_lifts`
  explicit lifts of quotient generators from `GL_2(q)` back into `GL(V)`;
- `full_generators`
  generators for the full stabilizer subgroup `Stab_{GL(V)}(L)`.

## Relation to the Paper

The names, counts, and stabilizer strings returned by `src/paper_data.m` are
meant to match:

- `paper_tables/geometric_orbit_tables_generated.tex`
- `paper_tables/geometric_stabilizer_tables_generated.tex`
- `paper_tables/finite_field_stabilizer_tables_generated.tex`

So the repository can be used both as a companion for the paper and as a
standalone finite-range data package.

## For Maintainers

The scripts in `dev/` are not needed for normal use. They are there to keep the
package synchronized with the paper:

- `dev/verify_counts_and_local_models.m`
  checks the orbit counts and the currently implemented local matrix models.
- `dev/verify_finite_stabilizer_specs.m`
  checks the full hard-coded finite-field stabilizer generators.
- `dev/check_orbit_count_polynomials.m`
  checks symbolically that the orbit-stabilizer sum over the hard-coded tables
  equals the Grassmann polynomial for the number of `2`-spaces in `\Lambda^2 V`
  for `n = 4,5,6,7`.
- `dev/generate_*_tables_tex.m`
  regenerate the TeX appendix tables.
- `dev/audit_paper_tables.py`
  checks that the repo-generated tables still match the bundled paper tables.

Typical maintenance commands:

```text
magma -b < dev/verify_counts_and_local_models.m
magma -b < dev/verify_finite_stabilizer_specs.m
magma -b < dev/check_orbit_count_polynomials.m
python3 dev/audit_paper_tables.py
```
