"""Draw an original occupancy schematic; no external drawing dependency."""


def main():
    from pathlib import Path

    ROOT = Path(__file__).resolve().parent.parent
    parts = [
        (
            (
                '<svg xmlns="http://www.w3.org/2000/svg" width="900" height="'
                '450" viewBox="0 0 900 450" role="img" aria-labelledby="title'
                ' desc">'
            )
        ),
        (
            '<title id="title">Same composition, different sublattice order</'
            "title>"
        ),
        (
            (
                '<desc id="desc">Each constructed pattern has eight sites, fo'
                "ur A and four B. Left: two A on each sublattice, eta zero. R"
                "ight: all four A on alpha and all four B on beta, eta one. P"
                "ositions depict labelled site groups, not a spatial projecti"
                "on of a crystal.</desc>"
            )
        ),
        '<rect width="900" height="450" rx="14" fill="#101827"/>',
    ]
    for base, alpha, beta, label in [
        (40, "AABB", "ABAB", "No sublattice bias: η = 0"),
        (480, "AAAA", "BBBB", "Perfect sublattice order: η = 1"),
    ]:
        parts.append(
            (
                f'<text x="'
                f"{base}"
                f'" y="46" fill="#e6edf3" font-size="21" font-family="sans-ser'
                f'if">'
                f"{label}"
                f"</text>"
            )
        )
        for y, atoms, name in [(125, alpha, "α"), (230, beta, "β")]:
            parts.append(
                (
                    f'<text x="'
                    f"{base}"
                    f'" y="'
                    f"{y - 35}"
                    f'" fill="#bbc7dc" font-size="18" font-family="sans-serif"'
                    f">"
                    f"{name}"
                    f" sites</text>"
                )
            )
            for i, atom in enumerate(atoms):
                color = "#7dd3fc" if atom == "A" else "#f5a97f"
                cx = base + 32 + i * 85
                parts.append(
                    (
                        f'<circle cx="'
                        f"{cx}"
                        f'" cy="'
                        f"{y}"
                        f'" r="26" fill="'
                        f"{color}"
                        f'"/><text x="'
                        f"{cx}"
                        f'" y="'
                        f"{y + 7}"
                        f'" text-anchor="middle" fill="#101827" font-size="23"'
                        f' font-family="sans-serif">'
                        f"{atom}"
                        f"</text>"
                    )
                )
        parts.append(
            (
                f'<text x="'
                f"{base}"
                f'" y="310" fill="#e6edf3" font-size="21" font-family="sans-se'
                f'rif">4 A + 4 B: atomic A fraction = 1/2</text>'
            )
        )
    parts.extend(
        [
            (
                (
                    '<text x="40" y="370" fill="#bbc7dc" font-size="19" font-'
                    'family="sans-serif">Composition fixes the total count; o'
                    "rder distinguishes the two site groups.</text>"
                )
            ),
            (
                (
                    '<text x="40" y="410" fill="#bbc7dc" font-size="17" font-'
                    'family="sans-serif">Constructed labelled-site groups; ne'
                    "ither a micrograph nor a crystal projection.</text>"
                )
            ),
            "</svg>",
        ]
    )
    (ROOT / "assets/sublattice-order.svg").write_text(
        "".join(parts), encoding="utf8"
    )
    print("Wrote original eight-site occupancy schematic.")


if __name__ == "__main__":
    main()
