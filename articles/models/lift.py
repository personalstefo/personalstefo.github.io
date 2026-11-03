"""Lift the assignment code out of a 02451 notebook, verbatim, so nothing is retyped."""
import json
import os
import sys

import numpy as np

BASE = r"D:\_incoming\ext\DTU\IntroToML\02451exercises"


def assignment_source(path):
    """Everything from the '# Assignment N' markdown cell onwards, as source + code."""
    nb = json.load(open(path, encoding="utf-8"))
    cells = nb["cells"]
    start = next(i for i, c in enumerate(cells)
                 if c["cell_type"] == "markdown"
                 and "# Assignment" in "".join(c["source"]))
    return cells[start:]


def code_only(cells):
    return "\n".join("".join(c["source"]) for c in cells if c["cell_type"] == "code")


def outputs_of(cells):
    for c in cells:
        if c["cell_type"] != "code":
            continue
        text = ""
        for o in c.get("outputs", []):
            if o.get("output_type") == "stream":
                text += "".join(o.get("text", []))
            elif "data" in o:
                text += "".join(o["data"].get("text/plain", []))
        if text.strip():
            yield "".join(c["source"]), text


def extract_functions(path, names):
    """Pull the source of the named top-level functions out of a notebook, verbatim."""
    cells = assignment_source(path)
    code = code_only(cells)
    out = {}
    for name in names:
        head = f"def {name}("
        i = code.find(head)
        assert i >= 0, f"{name} not found in {os.path.basename(path)}"
        # walk forward until a line starts at column zero and is not a continuation
        j = i
        lines = code[i:].splitlines(keepends=True)
        body = [lines[0]]
        for line in lines[1:]:
            if line.strip() and not line[0].isspace():
                break
            body.append(line)
        src = "".join(body).rstrip()
        assert "NotImplementedError" not in src, f"{name} is not implemented"
        out[name] = src
    return out


MINE = ("create_polynomial_features", "fit_pm", "predict_pm", "kfold_np")


if __name__ == "__main__":
    rel = sys.argv[1] if len(sys.argv) > 1 else r"02451exercise9\exercise9.ipynb"
    cells = assignment_source(os.path.join(BASE, rel))
    code = code_only(cells)
    print(f"{len(cells)} cells from the assignment onwards, {len(code)} chars of code")
    for name in ("train_error_rate", "test_error_rate", "opt_lambda", "min_error",
                 "coefficient_norm", "lambdas"):
        print(f"  {name}: {'defined' if name in code else 'MISSING'}")
    assert "LogisticRegression" in code, "assignment code did not lift"
    assert "np.logspace" in code, "lambda sweep missing"

    # the four functions the second article reuses, lifted verbatim
    fns = extract_functions(os.path.join(BASE, r"02451exercise07\exercise7.ipynb"), MINE)
    for name, src in fns.items():
        print(f"\n--- {name} ({len(src.splitlines())} lines)")
    print(f"\n{len(fns)} functions lifted")
    assert len(fns) == 4
