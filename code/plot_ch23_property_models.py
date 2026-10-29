"""Generate the chapter's analytic toy-model figure with standard Python.
No material data or proprietary property model is evaluated.
"""
from pathlib import Path
from math import exp, log, sqrt

ROOT = Path(__file__).resolve().parents[1]
R, TEMP = 8.31446261815324, 1000.
WIDTH, HEIGHT = 1100, 510
svg = [f'<svg xmlns="http://www.w3.org/2000/svg" '
       f'viewBox="0 0 {WIDTH} {HEIGHT}" role="img" '
       'aria-labelledby="title desc">',
       '<title id="title">Synthetic property assessment and obstacle crossover</title>',
       '<desc id="desc">Log-viscosity models agree at the midpoint but differ '
       'elsewhere. A separate cutting and bypass toy model has a finite peak.</desc>',
       '<rect width="1100" height="510" fill="#0c1220"/>',
       '<style>text{font:15px system-ui,sans-serif;fill:#e6edf6}'
       '.axis{stroke:#9aa9bd;stroke-width:1.2}.grid{stroke:#26344a;'
       'stroke-width:1}.curve{fill:none;stroke-width:2.6}</style>']


def text(x, y, value, size=15, anchor='start', fill='#e6edf6'):
    svg.append(f'<text x="{x:.2f}" y="{y:.2f}" text-anchor="{anchor}" '
               f'style="font-size:{size}px;fill:{fill}">{value}</text>')


def line(x1, y1, x2, y2, color='#9aa9bd', dashed=False):
    dash = ' stroke-dasharray="5 4"' if dashed else ''
    svg.append(f'<line x1="{x1:.3f}" y1="{y1:.3f}" x2="{x2:.3f}" '
               f'y2="{y2:.3f}" stroke="{color}"{dash}/>')


def curve(points, color, dashed=False):
    coords = ' '.join(f'{x:.3f},{y:.3f}' for x, y in points)
    dash = ' stroke-dasharray="6 4"' if dashed else ''
    svg.append(f'<polyline class="curve" points="{coords}" '
               f'stroke="{color}"{dash}/>')


text(70, 32, 'A. Independent log-viscosity assessment', 19)
text(620, 32, 'B. Competing obstacle mechanisms: toy increments', 19)
left, right, top, bottom = 75., 500., 151., 425.
px = lambda x: left+(right-left)*x
py = lambda y: bottom-(bottom-top)*(y-1.8)/(4.3-1.8)
for x in (0., .25, .5, .75, 1.):
    line(px(x), top, px(x), bottom, '#26344a')
    text(px(x), bottom+24, f'{x:g}', anchor='middle')
for y in (2., 2.5, 3., 3.5, 4.):
    line(left, py(y), right, py(y), '#26344a')
    text(left-10, py(y)+5, f'{y:g}', anchor='end')
line(left, top, left, bottom)
line(left, bottom, right, bottom)
text(left, top-13, 'Dynamic viscosity (mPa s)')
text((left+right)/2, bottom+57, 'B atomic fraction x', anchor='middle')
xs = [i/200 for i in range(201)]
for l0, l1, color, dash in [(0., 0., '#b9c2cf', True),
                           (6000., -2000., '#62d5e8', False),
                           (6000., 2000., '#f8b66b', False)]:
    values = [1000*exp((1-x)*log(.004)+x*log(.002)
                       + x*(1-x)*(l0+l1*(1-2*x))/(R*TEMP)) for x in xs]
    assert abs(values[0]-4)<1e-12 and abs(values[-1]-2)<1e-12
    curve([(px(x), py(y)) for x, y in zip(xs, values)], color, dash)
text(75, 65, 'Dashed: no excess interaction', 14, fill='#b9c2cf')
text(75, 87, 'Cyan: L0=6000, L1=−2000 J/mol', 14, fill='#62d5e8')
text(75, 109, 'Orange: same L0, L1=+2000 J/mol', 14, fill='#f8b66b')

left, right = 625., 1040.
px = lambda r: left+(right-left)*log(r/.5)/log(100)
py = lambda s: bottom-(bottom-top)*s/300
for r in (.5, 1., 2., 5., 10., 20., 50.):
    line(px(r), top, px(r), bottom, '#26344a')
    text(px(r), bottom+24, f'{r:g}', anchor='middle')
for s in (0, 50, 100, 150, 200, 250, 300):
    line(left, py(s), right, py(s), '#26344a')
    text(left-10, py(s)+5, str(s), anchor='end')
line(left, top, left, bottom)
line(left, bottom, right, bottom)
text(left, top-13, 'Strength increment (MPa)')
text((left+right)/2, bottom+57, 'Radius r (nm), logarithmic axis', anchor='middle')
radii = [.5*exp(log(100)*i/500) for i in range(501)]
curve([(px(r), py(40*sqrt(r))) for r in radii], '#62d5e8')
curve([(px(r), py(800/r)) for r in radii if 800/r <= 300], '#f8b66b')
curve([(px(r), py(min(40*sqrt(r), 800/r))) for r in radii], '#f276a4')
critical = 20**(2/3)
peak = 40*sqrt(critical)
line(px(critical), py(peak), px(critical), bottom, '#f276a4', True)
svg.append(f'<circle cx="{px(critical):.3f}" cy="{py(peak):.3f}" '
           'r="4" fill="#f276a4"/>')
line(px(critical)+8, py(peak)-8, 910, 215, '#9aa9bd')
text(920, 190, 'r* = 7.368 nm', 13)
text(920, 210, 'peak = 108.577 MPa', 13)
text(625, 65, 'Cyan: cutting = 40 sqrt(r / 1 nm)', 14, fill='#62d5e8')
text(625, 87, 'Orange: bypass = 800 / (r / 1 nm)', 14, fill='#f8b66b')
text(625, 109, 'Pink: lower threshold; no measured alloy fit', 14, fill='#f276a4')
text(550, 505, 'Original synthetic models; independent panels, not Thermo-Calc output',
     13, 'middle', '#b6c3d7')
svg.append('</svg>')
(ROOT/'assets').mkdir(exist_ok=True)
(ROOT/'assets/ch23_property_models.svg').write_bytes('\n'.join(svg).encode())
print(f'Synthetic crossover: {critical:.9f} nm, {peak:.6f} MPa')
