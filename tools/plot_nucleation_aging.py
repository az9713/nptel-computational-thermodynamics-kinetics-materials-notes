def main():
    from pathlib import Path
    import os
    import runpy
    import io
    import contextlib

    ROOT = Path(__file__).resolve().parent.parent
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib-cache"))
    with contextlib.redirect_stdout(io.StringIO()):
        v = runpy.run_path(str(ROOT / "code/ch21_nucleation_aging.py"))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    plt.rcParams.update(
        {
            "text.color": "#e7edf5",
            "axes.labelcolor": "#e7edf5",
            "xtick.color": "#c3d0df",
            "ytick.color": "#c3d0df",
            "axes.edgecolor": "#8394a8",
            "font.size": 11,
            "svg.fonttype": "none",
            "svg.hashsalt": "thermo-ch21",
        }
    )
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2), facecolor="#111827")
    for ax in axs:
        ax.set_facecolor("#111827")
        ax.grid(alpha=0.15, color="#94a3b8")
    ts = np.linspace(350, 900, 500)
    for fraction, color in [
        (0.01, "#67e8f9"),
        (0.5, "#fbbf24"),
        (0.99, "#c4b5fd"),
    ]:
        axs[0].semilogx(
            [v["t_fraction"](t, fraction) for t in ts],
            ts,
            color=color,
            label=f"{fraction:.0%} normalized extent",
        )
    axs[0].set(
        xlabel="Isothermal time (s)",
        ylabel="Temperature (K)",
        title="Stipulated Poisson nucleation / linear growth",
    )
    axs[0].legend(facecolor="#182335", labelcolor="#e7edf5", fontsize=9)
    rr = np.geomspace(1, 64, 500)
    cut = 100 * np.sqrt(rr / 8)
    bypass = 100 * 8 / rr
    axs[1].semilogx(rr, cut, "--", color="#fbbf24", label="Cutting stress")
    axs[1].semilogx(rr, bypass, ":", color="#c4b5fd", label="Bypass stress")
    axs[1].semilogx(
        rr,
        np.minimum(cut, bypass),
        color="#67e8f9",
        lw=2,
        label="Easier alternative",
    )
    axs[1].set(
        xlabel="Radius (nm)",
        ylabel="Obstacle shear stress (MPa)",
        title="Fixed fraction; original illustrative strength law",
        ylim=(0, 300),
    )
    axs[1].legend(facecolor="#182335", labelcolor="#e7edf5", fontsize=9)
    fig.tight_layout()
    fig.savefig(ROOT / "assets/nucleation-aging.svg", metadata={"Date": None})
    if __name__ == "__main__":
        print("Wrote assets/nucleation-aging.svg")


if __name__ == "__main__":
    main()
