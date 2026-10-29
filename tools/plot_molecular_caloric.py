"""Original classical speed distribution and Einstein-model heat capacity."""


def main():
    import os
    import sys
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent
    os.environ["MPLCONFIGDIR"] = str(ROOT / "work/matplotlib-cache")
    sys.path.insert(0, str(ROOT / "code"))
    from ch06_molecular_caloric import (
        KB,
        NA,
        scaled_speed_density,
        einstein_ratio,
    )
    from math import sqrt
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor="#101827")
    for ax in axes:
        ax.set_facecolor("#101827")
        ax.grid(alpha=0.15)
    v = [4 * i for i in range(501)]
    for t, color in [(300, "#7dd3fc"), (600, "#f5a97f")]:
        v0 = sqrt(2 * KB * t / (0.020 / NA))
        axes[0].plot(
            v,
            [scaled_speed_density(q / v0) / v0 for q in v],
            label=f"{t} K",
            color=color,
        )
    axes[0].set(
        xlabel="Speed / (m/s)",
        ylabel="Probability density / (s/m)",
        title="Classical gas: stipulated molar mass 20 g/mol",
    )
    axes[0].legend()
    x = [0.01 * i for i in range(1, 501)]
    axes[1].plot(x, [einstein_ratio(t) for t in x], color="#a6da95")
    axes[1].axhline(
        1, color="white", ls="--", lw=1, label="Classical 3R limit"
    )
    axes[1].set(
        xlabel="T / Einstein temperature",
        ylabel="Molar Cv / (3R)",
        title="Independent equal-frequency harmonic oscillators",
    )
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(ROOT / "assets/molecular-caloric.svg")
    fig.savefig(ROOT / "work/reader/ch06-molecular.png", dpi=140)
    plt.close(fig)
    print(
        "Wrote original Maxwell-speed and Einstein heat-capacity model figure."
    )


if __name__ == "__main__":
    main()
