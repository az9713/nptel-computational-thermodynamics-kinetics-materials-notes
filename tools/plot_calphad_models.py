"""Original synthetic ternary-interaction counterexample; no measured data."""


def main():
    from pathlib import Path
    import os

    ROOT = Path(__file__).resolve().parent.parent
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib-cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as tri
    import numpy as np

    pts = [(i / 40, j / 40) for i in range(41) for j in range(41 - i)]
    a = np.array([p[0] for p in pts])
    b = np.array([p[1] for p in pts])
    c = 1 - a - b
    X = b + 0.5 * c
    Y = np.sqrt(3) * c / 2
    mesh = tri.Triangulation(X, Y)
    base = a * b * 1000 + a * c * 2000 + b * c * (-3000)
    correction = -100000 * a * b * c
    plt.rcParams.update(
        {
            "text.color": "#e7edf5",
            "axes.labelcolor": "#e7edf5",
            "font.size": 11,
            "svg.fonttype": "none",
        }
    )
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.8), facecolor="#111827")
    for ax, values, title in zip(
        axs,
        [base, base + correction, correction],
        [
            "Binary-based continuation",
            "Same binary edges, different interior",
            "Added ternary interaction",
        ],
    ):
        ax.set_facecolor("#111827")
        m = ax.tricontourf(
            mesh,
            values,
            levels=np.linspace(-4000, 600, 24),
            cmap="viridis",
            extend="both",
        )
        ax.plot(
            [0, 1, 0.5, 0], [0, 0, np.sqrt(3) / 2, 0], color="#e7edf5", lw=0.8
        )
        for xx, yy, label in [
            (0, -0.04, "A"),
            (1, -0.04, "B"),
            (0.5, 0.91, "C"),
        ]:
            ax.text(xx, yy, label, ha="center")
        ax.set_aspect("equal")
        ax.axis("off")
        ax.set_title(title, pad=18, color="#e7edf5")
    fig.subplots_adjust(
        left=0.03, right=0.87, bottom=0.17, top=0.82, wspace=0.18
    )
    cb = fig.colorbar(m, cax=fig.add_axes([0.9, 0.22, 0.018, 0.50]))
    cb.set_label("Excess Gibbs energy (J/mol)")
    cb.ax.tick_params(colors="#e7edf5")
    fig.text(
        0.46,
        0.035,
        (
            (
                "Original teaching models: pair coefficients 1000, 2000, -300"
                "0; added LABC = -100000 J/mol"
            )
        ),
        ha="center",
        color="#e7edf5",
    )
    fig.savefig(
        ROOT / "assets/calphad-continuation.svg", facecolor=fig.get_facecolor()
    )
    fig.savefig(
        ROOT / "work/author/ch13-figure.png",
        facecolor=fig.get_facecolor(),
        dpi=130,
    )


if __name__ == "__main__":
    main()
