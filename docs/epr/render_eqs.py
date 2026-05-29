#!/usr/bin/env python3
"""Render LaTeX equations to SVG (white card) for embedding in markdown.

Warp's markdown viewer has no math engine, so we pre-render equations to images.
Usage: define equations in the EQS dict of a doc, or import render() and call it.

    from render_eqs import render
    render("01_eq07", r"\varphi_m^2 = p_m\,\frac{\hbar\omega_m}{2E_J}")

Outputs docs/epr/img/<name>.svg (paths embedded, white background, black text).
"""
import os
import subprocess
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, "img")

TEX = r"""\documentclass[border=6pt]{standalone}
\usepackage{amsmath,amssymb,xcolor}
\begin{document}
\colorbox{white}{\color{black}$\displaystyle %s$}
\end{document}
"""


def render(name: str, latex: str) -> str:
    """Render one display equation to img/<name>.svg. Returns the svg path."""
    os.makedirs(IMG, exist_ok=True)
    with tempfile.TemporaryDirectory() as d:
        tex = os.path.join(d, "eq.tex")
        with open(tex, "w") as f:
            f.write(TEX % latex)
        subprocess.run(["latex", "-interaction=nonstopmode", "eq.tex"],
                       cwd=d, check=True, capture_output=True)
        out = os.path.join(IMG, f"{name}.svg")
        subprocess.run(["dvisvgm", "--no-fonts", "--exact",
                        "-o", out, os.path.join(d, "eq.dvi")],
                       cwd=d, check=True, capture_output=True)
    return out


def render_all(eqs: dict):
    for name, latex in eqs.items():
        render(name, latex)
        print(f"  img/{name}.svg")


if __name__ == "__main__":
    import json
    import sys
    # Read {name: latex} JSON from argv[1] (file) or stdin.
    src = open(sys.argv[1]).read() if len(sys.argv) > 1 else sys.stdin.read()
    render_all(json.loads(src))
