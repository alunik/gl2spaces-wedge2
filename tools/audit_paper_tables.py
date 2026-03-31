#!/usr/bin/env python3
from __future__ import annotations

import subprocess
from pathlib import Path
import shlex

REPO = Path(__file__).resolve().parents[1]
PAPER = REPO / "paper_tables"


def run_magma(script: Path, output: Path) -> None:
    cmd = f"magma -b < {shlex.quote(str(script))}"
    with output.open("w", encoding="utf-8") as fh:
        subprocess.run(
            cmd,
            cwd=REPO,
            stdout=fh,
            stderr=subprocess.PIPE,
            text=True,
            check=True,
            shell=True,
            executable="/bin/bash",
        )


def assert_equal(expected: Path, actual: Path, label: str) -> None:
    def normalize(text: str) -> str:
        kept = []
        for line in text.splitlines():
            if line.startswith("Loading "):
                continue
            if line.startswith("Magma V2."):
                continue
            if line.startswith("Type ? for help."):
                continue
            if line.startswith("Total time:"):
                continue
            kept.append(line)
        return "\n".join(kept).strip() + "\n"

    exp = normalize(expected.read_text(encoding="utf-8"))
    act = normalize(actual.read_text(encoding="utf-8"))
    if exp != act:
        raise SystemExit(f"{label} mismatch:\n  expected: {expected}\n  actual:   {actual}")
    print(f"[ok] {label}")


def assert_contains_lines(tex_file: Path, lines: list[str], label: str) -> None:
    text = tex_file.read_text(encoding="utf-8")
    missing = [line for line in lines if line not in text]
    if missing:
        msg = "\n".join(missing)
        raise SystemExit(f"{label} is missing expected lines:\n{msg}")
    print(f"[ok] {label}")


def main() -> None:
    outdir = REPO / "build_audit"
    outdir.mkdir(exist_ok=True)

    generated = [
        (
            REPO / "examples" / "generate_geometric_orbit_tables_tex.m",
            outdir / "geometric_orbit_tables_generated.tex",
            PAPER / "geometric_orbit_tables_generated.tex",
            "geometric orbit tables",
        ),
        (
            REPO / "examples" / "generate_geometric_stabilizer_tables_tex.m",
            outdir / "geometric_stabilizer_tables_generated.tex",
            PAPER / "geometric_stabilizer_tables_generated.tex",
            "geometric stabilizer tables",
        ),
        (
            REPO / "examples" / "generate_finite_field_stabilizer_tables_tex.m",
            outdir / "finite_field_stabilizer_tables_generated.tex",
            PAPER / "finite_field_stabilizer_tables_generated.tex",
            "finite-field stabilizer tables",
        ),
    ]

    for script, tmp, target, label in generated:
        run_magma(script, tmp)
        assert_equal(target, tmp, label)

    assert_contains_lines(
        PAPER / "local_models_appendix.tex",
        [
            r"\(S_d\) &",
            r"\(\mathbf G_a^{2d-2}\rtimes \mathbf G_m\) &",
            r"\(J_{a,2}\oplus J_{b,1}\) &",
            r"\(\SL_2(k[t]/(t^2))\times \SL_2(k)\) &",
            r"\(S_2\oplus L_0\) &",
            r"\(U_L\rtimes (k^\times\times K_{L_0})\) &",
            r"\(S_1^r\oplus L_0\) &",
            r"\(\Hom(R,W)\rtimes (\GL(W)\times \Stab_{\GL(R)}(L_0))\) &",
        ],
        "local model table (algebraic side)",
    )

    assert_contains_lines(
        PAPER / "local_models_appendix.tex",
        [
            r"pure singular &",
            r"\(\GL_2(q)\) &",
            r"one rational support point &",
            r"\(B(q)\) &",
            r"weighted rational pair &",
            r"\(T_s(q)\) &",
            r"one cubic orbit &",
            r"\(\widetilde{C}_3(q)\) &",
        ],
        "local model table (finite-field side)",
    )

    assert_contains_lines(
        PAPER / "lang_steinberg.tex",
        [
            r"pure singular &",
            r"\(\GL_2(q)\) &",
            r"quadratic pair &",
            r"\(N_{ns}(q)\) &",
            r"one cubic orbit &",
            r"\(\widetilde{C}_3(q)\) &",
        ],
        "finite-field family summary table",
    )

    print("[ok] all audited paper tables match the repo")


if __name__ == "__main__":
    main()
