"""Two-layer cross-validation on the Islay regression data, using functions lifted verbatim.

The feature builder, the least-squares fit, the predictor and the fold splitter all come
from `lift.py`, which reads them out of Stefan's assignment notebook unchanged. The tree is
from scikit-learn because there is no reason to write one.
"""
import numpy as np
from sklearn.tree import DecisionTreeRegressor

from lift import BASE, extract_functions
import os

FNS = extract_functions(os.path.join(BASE, r"02451exercise07\exercise7.ipynb"),
                        ("create_polynomial_features", "fit_pm", "predict_pm", "kfold_np"))
env = {"np": np}
for name, src in FNS.items():
    exec(compile(src, f"<lifted {name}>", "exec"), env)
create_polynomial_features = env["create_polynomial_features"]
fit_pm = env["fit_pm"]
predict_pm = env["predict_pm"]
kfold_np = env["kfold_np"]

K_OUTER = 5
K_INNER = 5
DEGREES = [1, 2, 3, 5, 7]
DEPTHS = [4, 8, 12, 16, 20, 24, 30]
SEED = 7


def load():
    """The same frame the assignment regresses: alcohol predicted from the other columns."""
    import pandas as pd
    df = pd.read_csv(os.path.join(BASE, r"02451exercise07\data\wine2.csv"))
    X = df.drop(columns=["Color", "Alcohol"]).to_numpy()
    y = df["Alcohol"].to_numpy().reshape(-1, 1)
    return X, y


def cv_error(X_train, y_train, X_test, y_test, model_fn):
    return float(np.mean(np.abs(y_test - model_fn(X_train, y_train, X_test))))


def main():
    X, y = load()
    print(f"{X.shape[0]} samples, {X.shape[1]} input columns")

    pm_errors, tree_errors, chosen_degree, chosen_depth = [], [], [], []
    outer = kfold_np(len(X), n_splits=K_OUTER, random_state=SEED)

    for fold, split in enumerate(outer):
        X_tr, y_tr = X[split["train_idx"]], y[split["train_idx"]]
        X_te, y_te = X[split["test_idx"]], y[split["test_idx"]]

        inner = kfold_np(len(X_tr), n_splits=K_INNER, random_state=SEED + 1)

        # inner loop: pick the degree by mean absolute error on held-out inner folds
        best_deg, best_deg_err = None, np.inf
        for deg in DEGREES:
            errs = []
            for isp in inner:
                a, b = X_tr[isp["train_idx"]], y_tr[isp["train_idx"]]
                c, d = X_tr[isp["test_idx"]], y_tr[isp["test_idx"]]
                w = fit_pm(a, b, deg)
                errs.append(np.mean(np.abs(d - predict_pm(w, c, deg))))
            if np.mean(errs) < best_deg_err:
                best_deg, best_deg_err = deg, np.mean(errs)

        best_dep, best_dep_err = None, np.inf
        for dep in DEPTHS:
            errs = []
            for isp in inner:
                a, b = X_tr[isp["train_idx"]], y_tr[isp["train_idx"]]
                c, d = X_tr[isp["test_idx"]], y_tr[isp["test_idx"]]
                t = DecisionTreeRegressor(max_depth=dep, random_state=SEED).fit(a, b.ravel())
                errs.append(np.mean(np.abs(d.ravel() - t.predict(c))))
            if np.mean(errs) < best_dep_err:
                best_dep, best_dep_err = dep, np.mean(errs)

        w = fit_pm(X_tr, y_tr, best_deg)
        pm_err = float(np.mean(np.abs(y_te - predict_pm(w, X_te, best_deg))))
        tree = DecisionTreeRegressor(max_depth=best_dep, random_state=SEED).fit(X_tr, y_tr.ravel())
        tree_err = float(np.mean(np.abs(y_te.ravel() - tree.predict(X_te))))

        pm_errors.append(pm_err)
        tree_errors.append(tree_err)
        chosen_degree.append(best_deg)
        chosen_depth.append(best_dep)
        print(f"fold {fold + 1}: degree {best_deg}, depth {best_dep}, "
              f"PM {pm_err:.4f}, tree {tree_err:.4f}")

    pm_errors, tree_errors = np.array(pm_errors), np.array(tree_errors)
    diff = pm_errors - tree_errors

    # correlated t-test over the paired folds: the folds share data, so the
    # differences are paired and the variance is of the differences, not of each model
    n = len(diff)
    mean_d = float(diff.mean())
    se = float(diff.std(ddof=1) / np.sqrt(n))
    from scipy import stats
    t_stat = mean_d / se if se > 0 else 0.0
    p = float(2 * stats.t.sf(abs(t_stat), df=n - 1))
    crit = float(stats.t.ppf(0.975, df=n - 1))

    print(f"\nmean PM error   {pm_errors.mean():.4f}")
    print(f"mean tree error {tree_errors.mean():.4f}")
    print(f"difference      {mean_d:.4f}  95% CI [{mean_d - crit * se:.4f}, {mean_d + crit * se:.4f}]")
    print(f"t = {t_stat:.3f}, p = {p:.3f} over {n} paired folds")
    print(f"degrees chosen  {chosen_degree}")
    print(f"depths chosen   {chosen_depth}")

    # plot the paired differences, because a five-fold CI is being read as a finding
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    idx = np.arange(1, n + 1)
    ax.axhline(0, color="#9c948a", linewidth=0.8)
    ax.axhline(mean_d, color="#c1440e", linewidth=1.2, linestyle="--",
               label=f"mean {mean_d:.4f}")
    ax.fill_between([0.5, n + 0.5], mean_d - crit * se, mean_d + crit * se,
                    color="#c1440e", alpha=0.12, label="95% CI")
    ax.plot(idx, diff, "o", color="#1a1a1a", ms=7)
    ax.set(xlabel="outer fold", ylabel="PM error − tree error",
           title="paired difference per fold")
    ax.set_xticks(idx)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig("figures/fold_diff.png", dpi=150)
    print("figure written: figures/fold_diff.png")

    assert len(chosen_degree) == K_OUTER
    assert pm_errors.mean() > 0 and tree_errors.mean() > 0
    assert len(diff) == K_OUTER


if __name__ == "__main__":
    main()
