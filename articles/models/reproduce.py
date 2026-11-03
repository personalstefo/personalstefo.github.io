"""Re-run the assignment 9 regularisation sweep on the same data, with the notebook's own code.

The sweep, the standardisation and the model come from `lift.py`, lifted verbatim from the
notebook. Only the plotting loop is rewritten, because the original is a two-panel figure
with a hard y-limit set to [-0.1, 4] that clips the curve.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression

SEED = 123
np.random.seed(SEED)

M = 11


def run():
    df = pd.read_csv("data/wine.csv")
    X = df.drop(columns=["Color"])
    y = df["Color"]

    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        X, y, test_size=0.92, stratify=y, random_state=1234)

    lambdas = np.logspace(-8, 2, 50)

    mu_sd = np.array([X_train.mean(axis=0), X_train.std(axis=0)])
    X_train = (X_train - mu_sd[0]) / mu_sd[1]
    X_test = (X_test - mu_sd[0]) / mu_sd[1]

    train_error = np.zeros(len(lambdas))
    test_error = np.zeros(len(lambdas))
    coef_norm = np.zeros(len(lambdas))

    for i, lam in enumerate(lambdas):
        model = LogisticRegression(penalty="l2", C=1 / lam, solver="liblinear",
                                   max_iter=20000, tol=1e-7, random_state=SEED)
        model.fit(X_train, y_train)
        train_error[i] = np.sum(model.predict(X_train) != y_train) / len(y_train)
        test_error[i] = np.sum(model.predict(X_test) != y_test) / len(y_test)
        coef_norm[i] = np.linalg.norm(model.coef_)

    return lambdas, train_error, test_error, coef_norm


def main():
    lambdas, train_error, test_error, coef_norm = run()
    best = int(np.argmin(test_error))

    print(f"{len(lambdas)} values of lambda, {lambdas[0]:.1e} to {lambdas[-1]:.1e}")
    print(f"train error {train_error.min():.4f} to {train_error.max():.4f}")
    print(f"test  error {test_error.min():.4f} at lambda {lambdas[best]:.3g}")
    print(f"coefficient norm {coef_norm.max():.3f} -> {coef_norm.min():.3f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].semilogx(lambdas, train_error * 100, ".-", color="#3d5a80", label="train")
    ax[0].semilogx(lambdas, test_error * 100, ".-", color="#c1440e", label="test")
    ax[0].plot(lambdas[best], test_error[best] * 100, "*", ms=13, color="#e9c46a",
               markeredgecolor="k")
    ax[0].set(xlabel="regularisation strength", ylabel="error rate (%)",
              title="classification error")
    ax[0].legend(frameon=False)

    ax[1].semilogx(lambdas, coef_norm, color="k")
    ax[1].set(xlabel="regularisation strength", ylabel=r"$\|w\|_2$",
              title="coefficient norm")
    fig.tight_layout()
    fig.savefig("figures/regularisation.svg")
    fig.savefig("figures/regularisation.png", dpi=150)

    assert test_error[best] < 0.3, "error rate is too high, the split or labels are wrong"
    assert coef_norm[0] > coef_norm[-1], "regularisation is not shrinking the weights"
    print("ok")


if __name__ == "__main__":
    main()
