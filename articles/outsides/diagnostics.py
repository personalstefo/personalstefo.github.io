"""Read the MrBayes and BEAST2 output files and check the diagnostics quoted on the page.

Nothing here re-runs an analysis. It parses the `.pstat` and `.log` files those programs
wrote and reports the numbers the write-up stands on, so the page can be checked against
its own evidence rather than taken on trust.

The BEAST2 log carries its parameter list in the XML header at the top of the file and the
column names as a tab-separated row further down, so the header row is found by looking for
the line that starts with `Sample` rather than by assuming a fixed offset.
"""
import os

import numpy as np

SOURCE = r"D:\_incoming\ext\DTU\MolEvol\finalmolevol"
BURN_IN = 0.1


def pstat(path):
    """MrBayes .pstat: one row per parameter, tab separated, PSRF in the last column."""
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    header = lines[1].split("\t")
    rows = [ln.split("\t") for ln in lines[2:] if ln.strip()]
    return header, rows


def beast_samples(path):
    """BEAST2 log as a dict of column name -> array of posterior samples."""
    lines = open(path, encoding="utf-8", errors="replace").read().splitlines()
    start = next(i for i, ln in enumerate(lines) if ln.startswith("Sample\t"))
    names = [n.strip() for n in lines[start].split("\t")]
    body = [ln.split("\t") for ln in lines[start + 1:] if ln.strip() and not ln.startswith("#")]
    n_drop = int(len(body) * BURN_IN)
    body = body[n_drop:]
    out = {}
    for i, name in enumerate(names):
        try:
            out[name] = np.array([float(r[i]) for r in body])
        except (ValueError, IndexError):
            pass
    return out


def ess_act(values):
    """Effective sample size from the integrated autocorrelation time, the standard estimator."""
    x = np.asarray(values, dtype=float)
    n = len(x)
    x = x - x.mean()
    var = float(np.dot(x, x) / n)
    if var == 0:
        return 0.0
    # autocovariance up to the first negative pair, as in BEAST's own ESS
    ac = np.correlate(x, x, mode="full")[n - 1:] / n
    s = 1.0
    for i in range(1, n // 2):
        pair = ac[i] + ac[i + 1]
        if pair <= 0:
            break
        s += 2.0 * pair / var
    return float(n / max(s, 1e-9))


def main():
    print("== MrBayes, HA alignment ==")
    header, rows = pstat(os.path.join(SOURCE, "HA_aligned.nexus.pstat"))
    psrf = np.array([float(r[-1]) for r in rows])
    ess = np.array([float(r[6]) for r in rows])
    names = [r[0] for r in rows]
    print(f"{len(rows)} parameters")
    print(f"PSRF {psrf.min():.6f} to {psrf.max():.6f}")
    print(f"ESS  {ess.min():.0f} to {ess.max():.0f}")

    inv = [(n, float(r[6])) for n, r in zip(names, rows) if n.startswith("pi")]
    if inv:
        vals = [v for _, v in inv]
        print(f"invariant-site proportions ({len(inv)} values): ESS {min(vals):.0f} to {max(vals):.0f}")

    print("\n== BEAST2, HA alignment ==")
    s = beast_samples(os.path.join(SOURCE, "HA_log.log"))
    n = len(s["posterior"])
    print(f"{n} samples after {BURN_IN:.0%} burn-in removed, from "
          f"{len(s['posterior']) / (1 - BURN_IN):.0f} retained")

    print("\nclock rate per codon partition (posterior mean):")
    rates = []
    for part in ("1", "2", "3"):
        vals = s[f"clockRate.{part}"]
        rates.append(vals.mean())
        print(f"  partition {part}: {vals.mean():.5f}  (+- {vals.std():.5f})")
    print(f"  partition 2 / partition 1 = {rates[1] / rates[0]:.2f}x")

    print("\nESS recomputed on this log:")
    track = ["posterior", "prior", "likelihood", "Tree.height",
             "clockRate.1", "clockRate.2", "clockRate.3",
             "rateAC.1", "rateAG.1", "rateCT.1"]
    for k in track:
        if k in s:
            print(f"  {k:16} ESS {ess_act(s[k]):8.0f}")

    print("\nsubstitution rates across all three partitions:")
    for k in ("rateAC", "rateAG", "rateCT", "rateGT"):
        vals = [ess_act(s[f"{k}.{p}"]) for p in ("1", "2", "3")]
        print(f"  {k}: ESS {min(vals):.0f} to {max(vals):.0f}")

    assert len(rows) > 20, "no .pstat rows parsed"
    assert len(rates) == 3, "expected three codon partitions"
    assert rates[1] > rates[0], "second partition should be the fast one"


if __name__ == "__main__":
    main()
