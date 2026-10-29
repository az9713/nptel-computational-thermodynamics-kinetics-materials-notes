def main():
    from pathlib import Path
    import os

    ROOT = Path(__file__).resolve().parent.parent
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib-cache"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    R = 8.3145
    D = 4000.0
    T = 1080.0
    delta = 10 * (1000 - T)
    u = np.exp((D + delta) / (R * T))
    v = np.exp((delta - D) / (R * T))
    b = (u - 1) / (u - v)
    a = v * b
    xx = np.linspace(0.18, 0.63, 600)
    q = (1 - xx) * np.log(1 - xx) + xx * np.log(xx)
    muA = R * T * np.log(1 - a)
    m = D + R * T * np.log(a / (1 - a))
    z = 0.4
    cross = D * z + R * T * ((1 - z) * np.log(1 - z) + z * np.log(z))
    tie = muA + m * z
    plt.rcParams.update(
        {
            "text.color": "#e7edf5",
            "axes.labelcolor": "#e7edf5",
            "xtick.color": "#c3d0df",
            "ytick.color": "#c3d0df",
            "axes.edgecolor": "#8394a8",
            "font.size": 11,
            "svg.fonttype": "none",
            "svg.hashsalt": "thermo-ch15",
        }
    )
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.1), facecolor="#111827")
    for ax in axs:
        ax.set_facecolor("#111827")
        ax.grid(alpha=0.15, color="#94a3b8")
    ax = axs[0]
    ax.plot(xx, D * xx + R * T * q, color="#67e8f9", label="Homogeneous P")
    ax.plot(
        xx,
        D * (1 - xx) + delta + R * T * q,
        color="#fbbf24",
        label="Homogeneous Q",
    )
    ax.plot(xx, muA + m * xx, color="#f472b6", label="Common supporting line")
    ax.scatter([z, z], [cross, tie], color="#e7edf5", s=24)
    ax.annotate(
        "",
        xy=(z, cross),
        xytext=(z, tie),
        arrowprops={"arrowstyle": "<->", "color": "#e7edf5"},
    )
    ax.text(
        0.415, (cross + tie) / 2, "212.05 J/mol", color="#e7edf5", fontsize=9
    )
    ax.set(
        xlim=(0.18, 0.63),
        ylim=(-5600, -3000),
        xlabel="x(NI), teaching variable",
        ylabel="Gibbs energy (J/mol)",
        title="Equal homogeneous energies at 1080 K",
    )
    ax.legend(facecolor="#182335", labelcolor="#e7edf5", fontsize=9)
    temps = np.linspace(650, 1350, 700)
    dd = 10 * (1000 - temps)
    uu = np.exp((D + dd) / (R * temps))
    vv = np.exp((dd - D) / (R * temps))
    bb = (uu - 1) / (uu - vv)
    aa = vv * bb
    ax = axs[1]
    ax.fill_betweenx(
        temps, aa, bb, color="#365676", alpha=0.6, label="Two-phase interval"
    )
    ax.plot(aa, temps, color="#67e8f9", label="P contact")
    ax.plot(bb, temps, color="#fbbf24", label="Q contact")
    ax.plot(
        (1400 - temps) / 800,
        temps,
        color="#f472b6",
        ls="--",
        label="T0: equal composition",
    )
    ax.scatter([z], [T], color="#e7edf5", s=24)
    ax.set(
        xlim=(0.02, 0.98),
        ylim=(650, 1350),
        xlabel="Atomic mole fraction x(NI)",
        ylabel="Temperature (K)",
        title="T0 differs from coexistence boundaries",
    )
    ax.legend(
        facecolor="#182335", labelcolor="#e7edf5", fontsize=9, loc="lower left"
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.98))
    fig.text(
        0.5,
        0.012,
        (
            (
                "Original synthetic two-branch model; analytic curves audited"
                " against executed pycalphad at the marked T0 point"
            )
        ),
        ha="center",
        color="#c3d0df",
        fontsize=10,
    )
    fig.savefig(
        ROOT / "assets/tzero-paths.svg",
        facecolor=fig.get_facecolor(),
        metadata={"Date": None},
    )
    fig.savefig(
        ROOT / "work/author/ch15-figure.png",
        facecolor=fig.get_facecolor(),
        dpi=130,
    )


if __name__ == "__main__":
    main()
