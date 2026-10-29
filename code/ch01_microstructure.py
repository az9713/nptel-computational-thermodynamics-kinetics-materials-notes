"""Original phase-budget and section-geometry checks; standard library only."""
from fractions import Fraction

water, sugar, capacity = 100, 250, 220
solid = max(0, sugar - capacity)
liquid = water + min(sugar, capacity)
total = water + sugar
w_liquid = Fraction(min(sugar, capacity), liquid)
w_overall = Fraction(sugar, total)
f_solid = Fraction(solid, total)
f_liquid = Fraction(liquid, total)
assert f_liquid + f_solid == 1
assert f_liquid * w_liquid + f_solid == w_overall
assert (w_overall - w_liquid) / (1 - w_liquid) == f_solid
cluster = {(2, 2), (2, 3), (3, 2), (3, 3)}
islands = {(1, 1), (1, 4), (4, 1), (4, 4)}


def boundary_edges(points, size=6):
    count = 0
    for x in range(size):
        for y in range(size):
            for dx, dy in [(1, 0), (0, 1)]:
                if x + dx < size and y + dy < size:
                    count += ((x, y) in points) != (
                        (x + dx, y + dy) in points)
    return count


assert len(cluster) == len(islands) == 4
assert boundary_edges(cluster) == 8
assert boundary_edges(islands) == 16
print(f'FCC atoms per conventional cell: {8/8 + 6/2:g}')
print(f'BCC atoms per conventional cell: {8/8 + 1:g}')
print('H atomic fraction in H2O:', Fraction(2, 3))
print('H mass fraction for stipulated masses H=1, O=16:', Fraction(2, 18))
print(f'overall sugar mass fraction: {float(w_overall):.8f}')
print(f'liquid sugar mass fraction: {float(w_liquid):.8f}')
print(f'solid sugar mass g: {solid}')
print(f'solid phase mass fraction: {float(f_solid):.8f}')
print('section B area fractions:', Fraction(len(cluster), 36),
      Fraction(len(islands), 36))
print('section interface lengths in cell-edge units:',
      boundary_edges(cluster), boundary_edges(islands))
