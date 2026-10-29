"""Reproduce the explicitly synthetic two-dimensional phase map."""


def main():
    from pathlib import Path

    root = Path(__file__).resolve().parent.parent
    sets = [{(2, 2), (2, 3), (3, 2), (3, 3)}, {(1, 1), (1, 4), (4, 1), (4, 4)}]
    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 680 310"',
        ' role="img" aria-labelledby="micro-title micro-desc">',
        '<title id="micro-title">Equal phase area, unequal interface</title>',
        '<desc id="micro-desc">Both constructed six-by-six maps have four',
        " B cells. One cluster has eight internal A-B edges; four isolated",
        " cells have sixteen. This is not a micrograph.</desc>",
        '<rect width="680" height="310" fill="#11151d"/>',
    ]
    for offset, points, label in zip([65, 385], sets, ["Cluster", "Islands"]):
        parts.append(
            f'<text x="{offset}" y="27" fill="#dce4ef" '
            f'font-size="19">{label}: B area = 4/36</text>'
        )
        for x in range(6):
            for y in range(6):
                color = "#f2bf78" if (x, y) in points else "#506882"
                parts.append(
                    f'<rect x="{offset+32*x}" y="{45+32*y}" '
                    f'width="32" height="32" fill="{color}" '
                    'stroke="#11151d" stroke-width="1"/>'
                )
        edges = 8 if offset == 65 else 16
        parts.append(
            f'<text x="{offset}" y="267" fill="#dce4ef" '
            f'font-size="17">A-B interface = {edges} cell edges</text>'
        )
    parts.append(
        '<text x="65" y="299" fill="#dce4ef" font-size="16">'
        "Blue: phase A. Gold: phase B. Grid-cell edge length = 1 μm."
        "</text></svg>"
    )
    (root / "assets/microstructure-map.svg").write_text(
        "\n".join(parts), encoding="utf8"
    )


if __name__ == "__main__":
    main()
