"""Original ternary balance geometry; the phase endpoints are stipulated."""


def main():
    import os
    import sys
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent
    os.environ["MPLCONFIGDIR"] = str(ROOT / "work/matplotlib-cache")
    sys.path.insert(0, str(ROOT / "code"))
    from ch04_ternary_balances import cartesian, weighted, F
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.style.use("dark_background")
    outer = [cartesian((1, 0, 0)), cartesian((0, 1, 0)), cartesian((0, 0, 1))]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5), facecolor="#101827")

    def base(ax, title):
        ax.set_facecolor("#101827")
        p = outer + [outer[0]]
        ax.plot(*zip(*p), color="#bbc7dc")
        for xy, name in zip(outer, ["A", "B", "C"]):
            ax.annotate(
                name,
                xy,
                xytext=(0, 8 if name == "C" else -18),
                textcoords="offset points",
                ha="center",
                fontsize=13,
            )
        for k in [0.2, 0.4, 0.6, 0.8]:
            for i in range(3):
                a = [0.0, 0.0, 0.0]
                b = [0.0, 0.0, 0.0]
                a[i] = b[i] = k
                a[(i + 1) % 3] = 1 - k
                b[(i + 2) % 3] = 1 - k
                ax.plot(
                    *zip(cartesian(a), cartesian(b)),
                    color="#6b7b92",
                    lw=0.5,
                    alpha=0.35
                )
        ax.set_aspect("equal")
        ax.set(xlim=(-0.06, 1.06), ylim=(-0.07, 0.97), title=title)
        ax.axis("off")

    p = [
        (F(4, 5), F(3, 20), F(1, 20)),
        (F(1, 10), F(4, 5), F(1, 10)),
        (F(3, 20), F(3, 20), F(7, 10)),
    ]
    f = [F(1, 5), F(3, 10), F(1, 2)]
    q = weighted(p, f)
    xy = list(map(cartesian, p))
    z = cartesian(q)
    base(axes[0], "Three fixed phase compositions; changing amounts")
    for i, color in enumerate(["#7dd3fc", "#f5a97f", "#a6da95"]):
        a, b = xy[(i + 1) % 3], xy[(i + 2) % 3]
        axes[0].fill(*zip(z, a, b), color=color, alpha=0.18)
        axes[0].plot(*zip(z, xy[i]), color="#bbc7dc", lw=1)
        axes[0].scatter(*xy[i], color=color, s=45)
        axes[0].annotate(
            ["α", "β", "γ"][i],
            xy[i],
            xytext=[(-12, -5), (8, -8), (6, 8)][i],
            textcoords="offset points",
        )
    axes[0].plot(*zip(*(xy + [xy[0]])), color="white", lw=1.5)
    axes[0].scatter(*z, color="white", s=45)
    axes[0].annotate(
        "Q",
        z,
        xytext=(-18, 3),
        textcoords="offset points",
        fontsize=10,
    )
    axes[0].text(
        0.12, -0.055, "Amounts: α 0.2; β 0.3; γ 0.5", fontsize=9
    )
    base(axes[1], "A vertical section can miss both tie endpoints")
    ends = [(F(7, 10), F(1, 10), F(1, 5)), (F(1, 10), F(3, 5), F(3, 10))]
    q = weighted(ends, [F(2, 5), F(3, 5)])
    section = [(F(3, 5), F(2, 5), 0), (0, F(2, 5), F(3, 5))]
    axes[1].plot(
        *zip(*map(cartesian, section)),
        color="#c4b5fd",
        lw=2,
        ls="--",
        label="Section: x_B = 0.4"
    )
    axes[1].plot(
        *zip(*map(cartesian, ends)),
        color="#f5a97f",
        lw=2,
        label="Actual tie line"
    )
    for item, name in zip(ends, ["α", "β"]):
        pos = cartesian(item)
        axes[1].scatter(*pos, color="#f5a97f", s=45)
        axes[1].annotate(name, pos, xytext=(0, 9), textcoords="offset points")
    pos = cartesian(q)
    axes[1].scatter(*pos, color="white", s=45)
    axes[1].annotate("Q", pos, xytext=(5, -16), textcoords="offset points")
    axes[1].legend(loc="upper right", fontsize=9, framealpha=0.2)
    fig.suptitle(
        "Composition geometry on one consistent atomic-amount basis",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(ROOT / "assets/ternary-balances.svg")
    fig.savefig(ROOT / "work/reader/ch04-ternary.png", dpi=140)
    plt.close(fig)
    print("Wrote original ternary tie-triangle and off-section schematic.")


if __name__ == "__main__":
    main()
