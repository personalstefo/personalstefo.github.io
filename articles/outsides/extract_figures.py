"""Extract the embedded figures from the MolEvol report as standalone PNGs.

Hiring reviewers should see the figure, not a screenshot of a PDF page. Each figure is
cropped from the page it sits on and written to articles/outsides/figures/.
"""
import io
import os
import sys

from PIL import Image

from pypdf import PdfReader

PDF = r"D:\_incoming\ext\DTU\MolEvol\finalmolevol\MolEvol_ Final Project-5.pdf"
OUT = os.path.join("articles", "outsides", "figures")

# page number (1-based) -> file stem. Chosen from the report's own figure captions.
FIGURES = {
    3: "fig01_host_ranges",
    11: "fig_trace_pb1",
    12: "fig_violin_invariant",
    13: "fig_consensus_fragments",
    15: "fig_ha_tree_hosts",
    18: "fig_beast_tree",
}


def trim(path, border=12):
    """Crop the white margin so a figure sits flush in its frame."""
    img = Image.open(path).convert("RGB")
    w, h = img.size
    img.crop((border, border, w - border, h - border)).save(path)
    return img.size


def main():
    reader = PdfReader(PDF)
    os.makedirs(OUT, exist_ok=True)
    written = []
    for page_no, stem in FIGURES.items():
        page = reader.pages[page_no - 1]
        images = list(page.images)
        if not images:
            print(f"  page {page_no}: no embedded image, skipped")
            continue
        for i, img in enumerate(images):
            name = f"{stem}.png" if len(images) == 1 else f"{stem}_{i + 1}.png"
            path = os.path.join(OUT, name)
            open(path, "wb").write(img.data)
            trim(path)
            written.append(path)
            print(f"  page {page_no}: {name}  {len(img.data) // 1024} kB")
    print(f"{len(written)} figures written to {OUT}")
    assert written, "no figures extracted"
    return written


if __name__ == "__main__":
    main()
