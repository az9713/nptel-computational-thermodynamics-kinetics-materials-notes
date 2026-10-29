"""Exact coin probabilities and original thermodynamic entropy ledger."""


def main():
    import os
    import sys
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent
    os.environ["MPLCONFIGDIR"] = str(ROOT / "work/matplotlib-cache")
    sys.path.insert(0, str(ROOT / "code"))
    from ch08_entropy_and_counts import coin_probability, R
    from math import sqrt, log
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor="#101827")
    for ax in axes:
        ax.set_facecolor("#101827")
        ax.grid(alpha=0.15)
    for n, color in [(20, "#7dd3fc"), (100, "#f5a97f"), (1000, "#a6da95")]:
        # Divide mass at a count by bin width1/n only for comparison as a
        # density-like scale.
        k = list(range(n + 1))
        axes[0].plot(
            [i / n for i in k],
            [n * coin_probability(n, i) for i in k],
            color=color,
            label=f"N={n}",
        )
    axes[0].set(
        xlabel="Fraction of heads H/N",
        ylabel="N times probability of exact count",
        xlim=(0.2, 0.8),
        title="Exact binomial masses; shrinking fractional width",
    )
    axes[0].legend()
    ds = R * log(2)
    labels = ["Reversible", "One-step load", "Free expansion"]
    reservoir = [-ds, -R / 2, 0]
    total = [0, ds - R / 2, ds]
    x = list(range(3))
    axes[1].bar(
        [i - 0.22 for i in x],
        [ds] * 3,
        width=0.22,
        label="Gas entropy change",
        color="#7dd3fc",
    )
    axes[1].bar(
        x,
        reservoir,
        width=0.22,
        label="Reservoir entropy change",
        color="#f5a97f",
    )
    axes[1].bar(
        [i + 0.22 for i in x],
        total,
        width=0.22,
        label="Composite generation",
        color="#a6da95",
    )
    axes[1].set(
        xticks=x,
        xticklabels=labels,
        ylabel="Entropy change / (J/K)",
        title="Same ideal-gas endpoints; different surroundings",
    )
    axes[1].axhline(0, color="white", lw=0.7)
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(ROOT / "assets/entropy-counts.svg")
    fig.savefig(ROOT / "work/reader/ch08-entropy.png", dpi=140)
    plt.close(fig)
    print(
        (
            (
                "Wrote original exact coin-concentration curves and entropy-t"
                "ransfer ledgers."
            )
        )
    )


if __name__ == "__main__":
    main()
