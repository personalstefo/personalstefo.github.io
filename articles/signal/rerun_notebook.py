"""Run the synthesis notebook's own cells in order and capture the figures.

Stefan's DSP notebook (`synthesis_draft1/synthesis.ipynb`) is the source. Its cells are
executed verbatim; the only change is the plotting backend, so `plt.show()` becomes a
`savefig` into `figures/`.

This exists so the article's figures are his output from his code, not a redrawing.
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")

NB = r"D:\_incoming\ext\DTU\DSP\synthesis_draft1\synthesis.ipynb"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures")
os.makedirs(OUT, exist_ok=True)

import matplotlib.pyplot as plt

_shown = [0]
_real_show = plt.show


def _show(*a, **k):
    """Replace plt.show with a save, numbered in call order."""
    _shown[0] += 1
    fig = plt.gcf()
    fig.savefig(os.path.join(OUT, f"cell{_shown[0]:02d}.png"), dpi=140, bbox_inches="tight")
    plt.close(fig)


def main():
    cells = json.load(open(NB, encoding="utf-8"))["cells"]
    code = ["".join(c["source"]) for c in cells if c["cell_type"] == "code"]
    assert len(code) == 8, f"expected 8 code cells, found {len(code)}"

    os.chdir(os.path.dirname(OUT))  # so figures/ and synth.wav land next to the article
    plt.show = _show

    env = {"__name__": "__main__"}
    for i, src in enumerate(code):
        try:
            exec(compile(src, f"<cell {i}>", "exec"), env)
        except Exception as exc:
            print(f"cell {i} failed: {type(exc).__name__}: {exc}")
            raise

    print(f"{_shown[0]} figures written to {OUT}")
    assert _shown[0] >= 5, "not all the notebook's plots were captured"
    assert os.path.exists(os.path.join(os.path.dirname(OUT), "synth.wav")), "no audio written"


if __name__ == "__main__":
    main()
