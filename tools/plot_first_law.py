"""Original work paths and endpoint comparison for the declared ideal gas."""


def main():
    import os
    import sys
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent
    os.environ["MPLCONFIGDIR"] = str(ROOT / "work/matplotlib-cache")
    sys.path.insert(0, str(ROOT / "code"))
    from ch07_first_law_paths import GAMMA, T1, CV, CP, R
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), facecolor="#101827")
    for ax in axes:
        ax.set_facecolor("#101827")
        ax.grid(alpha=0.15)
        ax.set_xlabel("V / initial volume")
        ax.set_ylabel("Pressure / initial pressure")
    x = [1 + i / 200 for i in range(201)]
    axes[0].plot(
        x, [1 / v for v in x], color="#7dd3fc", label="Reversible isothermal p"
    )
    axes[0].fill_between(
        x,
        [0.5] * len(x),
        [1 / v for v in x],
        color="#7dd3fc",
        alpha=0.18,
        label="Additional reversible work / (p1 V1)",
    )
    axes[0].plot(
        [1, 2],
        [0.5, 0.5],
        color="#f5a97f",
        ls="--",
        label="One-step external load",
    )
    axes[0].scatter([1, 2], [1, 0.5], color="white", zorder=5)
    axes[0].set(title="Same complete endpoints: T = 300 K")
    axes[0].legend(fontsize=8, loc="upper right")
    y = [1 + i * 0.7 / 200 for i in range(201)]
    vr = 2 ** (1 / GAMMA)
    axes[1].plot(
        y,
        [v ** (-GAMMA) for v in y],
        color="#a6da95",
        label="Reversible adiabatic equilibrium path",
    )
    axes[1].plot(
        [1, 1.6],
        [0.5, 0.5],
        color="#f5a97f",
        ls="--",
        label="One-step external load",
    )
    axes[1].scatter([1, vr, 1.6], [1, 0.5, 0.5], color="white", zorder=5)
    axes[1].annotate(
        "Reversible: 227.36 K",
        (vr, 0.5),
        xytext=(1.1, 0.75),
        arrowprops={"arrowstyle": "->", "color": "white"},
        fontsize=9,
    )
    axes[1].annotate(
        "One step: 240 K",
        (1.6, 0.5),
        xytext=(1.47, 0.36),
        arrowprops={"arrowstyle": "->", "color": "white"},
        fontsize=9,
    )
    axes[1].set(
        title="Same final pressure; different final states", ylim=(0.30, 1.05)
    )
    axes[1].legend(fontsize=8, loc="upper right")
    fig.suptitle(
        (
            (
                "Dashed lines describe imposed external pressure, not a gas e"
                "quilibrium trajectory"
            )
        )
    )
    fig.tight_layout()
    fig.savefig(ROOT / "assets/first-law-paths.svg")
    fig.savefig(ROOT / "work/reader/ch07-paths.png", dpi=140)
    plt.close(fig)
    print("Wrote original isothermal work and adiabatic endpoint comparison.")


if __name__ == "__main__":
    main()
