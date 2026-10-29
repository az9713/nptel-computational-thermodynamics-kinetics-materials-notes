def main():
    from pathlib import Path
    import os

    ROOT = Path(__file__).resolve().parent.parent
    os.environ.setdefault("MPLCONFIGDIR", str(ROOT / "work/matplotlib-cache"))
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
            "svg.hashsalt": "thermo-ch17",
        }
    )
    fig, axs = plt.subplots(1, 2, figsize=(12, 5.1), facecolor="#111827")
    for ax in axs:
        ax.set_facecolor("#111827")
        ax.grid(alpha=0.15, color="#94a3b8")
    aa = np.linspace(0, 16000, 500)
    ax = axs[0]
    ax.axvline(9000, color="#67e8f9", lw=2, label="Calorimetry fixes A")
    ax.plot(
        aa,
        (6000 - aa) / 1000,
        color="#fbbf24",
        lw=2,
        label="Activity fixes A+1000 B",
    )
    ax.scatter([9000], [-3], color="#f472b6", s=50, zorder=5)
    ax.set(
        xlim=(0, 16000),
        ylim=(-11, 8),
        xlabel="A (J/mol atoms)",
        ylabel="B (J/(mol atoms K))",
        title="Complementary synthetic observations",
    )
    ax.legend(facecolor="#182335", labelcolor="#e7edf5", fontsize=9)
    nn = np.arange(1, 65)
    correct = np.sqrt(4 / nn + 9)
    wrong = np.sqrt(13 / nn)
    ax = axs[1]
    ax.plot(nn, correct, color="#67e8f9", label="Shared calibration retained")
    ax.plot(
        nn,
        wrong,
        color="#fbbf24",
        label="Incorrect independent-total assumption",
    )
    ax.axhline(
        3, color="#f472b6", ls="--", label="Common standard-uncertainty floor"
    )
    ax.scatter([4], [np.sqrt(10)], color="#e7edf5")
    ax.set(
        xlabel="Number of repeated readings",
        ylabel="Mean standard uncertainty (arbitrary units)",
        title="A shared contribution does not average away",
        xlim=(1, 64),
        ylim=(0, 4),
    )
    ax.legend(facecolor="#182335", labelcolor="#e7edf5", fontsize=9)
    fig.tight_layout()
    fig.savefig(
        ROOT / "assets/assessment-identifiability.svg", metadata={"Date": None}
    )
    if __name__ == "__main__":
        print("Wrote assets/assessment-identifiability.svg")


if __name__ == "__main__":
    main()
