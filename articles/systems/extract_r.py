"""Extract every code chunk out of the 22150 notebooks and write it out unmodified.

The notebooks are the source. This only reads them: it strips the chunk header, finds the
matching closing fence, and writes the body byte for byte. If a chunk is empty it is skipped.
Nothing is reformatted, renamed or re-indented.
"""
import json
import os
import re

import matplotlib
matplotlib.use("Agg")

NB_DIR = r"C:\Users\User\dev\stefole-site\nbs"
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "chunks")

# order matters: it is the order the notebooks were written in, which is the course order
NOTEBOOKS = [
    ("week1-igraph.qmd", "01_setup_and_poldelta.R"),
    ("week2-louvian.qmd", "02_network_topology.R"),
    ("Week-3-gsea.qmd", "03_go_overrepresentation.R"),
    ("Week-4-virtual-pulldown.qmd", "04_virtual_pulldown.R"),
    ("Week-5.qmd", "05_deseq2.R"),
    ("Week-6-DEXSeq.qmd", "06_dexseq.R"),
]


def chunks(text):
    """Every ```{r ...} block in the document, in order, as (header, body)."""
    out = []
    for m in re.finditer(r"^```\{r([^}]*)\}\s*$", text, re.M):
        start = m.end() + 1
        end = text.find("\n```", start)
        body = text[start:end]
        out.append((m.group(1).strip(), body.rstrip("\n")))
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    total = 0
    for src, dst in NOTEBOOKS:
        text = open(os.path.join(NB_DIR, src), encoding="utf-8").read()
        blocks = chunks(text)
        assert blocks, f"no chunks found in {src}"
        body = "\n\n".join(b for _, b in blocks if b.strip())
        path = os.path.join(OUT_DIR, dst)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(body + "\n")
        n_lines = body.count("\n") + 1
        total += n_lines
        print(f"{src:32} -> {dst:28} {len(blocks):2} chunks, {n_lines:4} lines")

    print(f"\n{total} lines of R across {len(NOTEBOOKS)} files")
    # the last notebook is The Last One; keep the assert cheap and real
    final = open(os.path.join(OUT_DIR, "06_dexseq.R"), encoding="utf-8").read()
    assert "DEXSeqDataSet" in final, "week 6 chunk extraction lost the DEXSeq object"
    assert "fgseaMultilevel" in final, "week 6 chunk extraction lost the GSEA call"
    print("ok")


if __name__ == "__main__":
    main()
