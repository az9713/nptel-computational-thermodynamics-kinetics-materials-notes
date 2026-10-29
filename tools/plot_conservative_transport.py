def main():
    from pathlib import Path
    import os
    import runpy
    import io
    import contextlib

    ROOT = Path(__file__).resolve().parent.parent
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib-cache"))
    with contextlib.redirect_stdout(io.StringIO()):
        v = runpy.run_path(str(ROOT / "code/ch20_conservative_transport.py"))
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
            "svg.hashsalt": "thermo-ch20",
        }
    )
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.2), facecolor="#111827")
    for ax in axs:
        ax.set_facecolor("#111827")
        ax.grid(alpha=0.15, color="#94a3b8")
    nn = np.array([40, 80, 160])
    ee = np.array(v["errors"])
    axs[0].loglog(
        nn,
        ee,
        "o-",
        color="#67e8f9",
        label="Actual maximum cell-average error",
    )
    axs[0].loglog(
        nn,
        ee[0] * (40 / nn) ** 2,
        "--",
        color="#fbbf24",
        label="N⁻² reference",
    )
    axs[0].set(
        xlabel="Number of cells N",
        ylabel="Maximum error in mole fraction",
        title="Combined space/time refinement at Dt/L²=0.1",
    )
    axs[0].legend(facecolor="#182335", labelcolor="#e7edf5", fontsize=9)
    t = np.geomspace(0.01, 1e5, 400)
    aa, ell = v["aa"], v["ell"]
    s = 2 * aa * t / (np.sqrt(ell**2 + 2 * aa * t) + ell)
    axs[1].loglog(t, s * 1e6, color="#67e8f9", label="Mixed resistance")
    axs[1].loglog(
        t,
        aa * t / ell * 1e6,
        "--",
        color="#fbbf24",
        label="Early interface limit",
    )
    axs[1].loglog(
        t,
        np.sqrt(2 * aa * t) * 1e6,
        ":",
        color="#c4b5fd",
        label="Late diffusion limit",
    )
    axs[1].set(
        xlabel="Time (s)",
        ylabel="Interface position (µm)",
        title="One specified crossover law; synthetic parameters",
    )
    axs[1].legend(facecolor="#182335", labelcolor="#e7edf5", fontsize=9)
    fig.tight_layout()
    fig.savefig(
        ROOT / "assets/conservative-transport.svg", metadata={"Date": None}
    )
    if __name__ == "__main__":
        print("Wrote assets/conservative-transport.svg")


if __name__ == "__main__":
    main()
