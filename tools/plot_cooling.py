"""Original specified-model cooling curves, not experimental traces."""


def main():
    import os
    import sys
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent
    os.environ["MPLCONFIGDIR"] = str(ROOT / "work/matplotlib-cache")
    sys.path.insert(0, str(ROOT / "code"))
    from ch05_cooling_curves import pure, mushy, sensor_step
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    t = list(range(0, 3501, 5))
    p = [pure(i) for i in t]
    q = [mushy(i) for i in t]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), facecolor="#101827")
    for ax in axes:
        ax.set_facecolor("#101827")
        ax.grid(alpha=0.15)
        ax.set_xlabel("Elapsed time / s")
    axes[0].plot(
        t, [a[0] for a in p], label="Pure: equilibrium arrest", color="#7dd3fc"
    )
    axes[0].plot(
        t,
        [a[0] for a in q],
        label="Prescribed freezing range",
        color="#f5a97f",
    )
    axes[0].axhline(1000, color="#7dd3fc", lw=0.6, ls=":")
    axes[0].set_ylabel("True temperature / K")
    axes[0].legend(fontsize=9)
    axes[0].annotate(
        "Latent heat removed at constant T",
        (2000, 1000),
        xytext=(900, 1150),
        arrowprops={"arrowstyle": "->", "color": "white"},
        fontsize=9,
    )
    axes[1].plot(
        t, [a[1] for a in p], color="#7dd3fc", label="Pure solid mass fraction"
    )
    axes[1].plot(
        t,
        [a[1] for a in q],
        color="#f5a97f",
        label="Range-model solid mass fraction",
    )
    axes[1].set_ylabel("Solid mass fraction")
    axes[1].set_ylim(-0.03, 1.08)
    axes[1].legend(fontsize=9)
    fig.suptitle(
        "One kilogram; 100 W heat extraction; declared caloric models"
    )
    fig.tight_layout()
    fig.savefig(ROOT / "assets/cooling-models.svg")
    fig.savefig(ROOT / "work/reader/ch05-cooling.png", dpi=140)
    plt.close(fig)
    print("Wrote original cooling-curve and phase-amount figure.")


if __name__ == "__main__":
    main()
