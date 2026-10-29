"""Draw an original occupancy schematic; no external drawing dependency."""
from pathlib import Path
ROOT=Path(__file__).resolve().parent.parent
parts=['<svg xmlns="http://www.w3.org/2000/svg" width="900" height="450" viewBox="0 0 900 450" role="img" aria-labelledby="title desc">', '<title id="title">Same composition, different sublattice order</title>', '<desc id="desc">Each constructed pattern has eight sites, four A and four B. Left: two A on each sublattice, eta zero. Right: all four A on alpha and all four B on beta, eta one. Positions depict labelled site groups, not a spatial projection of a crystal.</desc>', '<rect width="900" height="450" rx="14" fill="#101827"/>']
for base,alpha,beta,label in [(40,'AABB','ABAB','No sublattice bias: eta = 0'),(480,'AAAA','BBBB','Perfect sublattice order: eta = 1')]:
 parts.append(f'<text x="{base}" y="46" fill="#e6edf3" font-size="21" font-family="sans-serif">{label}</text>')
 for y,atoms,name in [(125,alpha,'alpha'),(230,beta,'beta')]:
  parts.append(f'<text x="{base}" y="{y-35}" fill="#bbc7dc" font-size="18" font-family="sans-serif">{name} sites</text>')
  for i,atom in enumerate(atoms):
   color='#7dd3fc' if atom=='A' else '#f5a97f';cx=base+32+i*85
   parts.append(f'<circle cx="{cx}" cy="{y}" r="26" fill="{color}"/><text x="{cx}" y="{y+7}" text-anchor="middle" fill="#101827" font-size="23" font-family="sans-serif">{atom}</text>')
 parts.append(f'<text x="{base}" y="310" fill="#e6edf3" font-size="21" font-family="sans-serif">4 A + 4 B: atomic x_A = 1/2</text>')
parts.extend(['<text x="40" y="370" fill="#bbc7dc" font-size="19" font-family="sans-serif">Composition fixes the total count; order distinguishes the two site groups.</text>','<text x="40" y="410" fill="#bbc7dc" font-size="17" font-family="sans-serif">Constructed labelled-site groups; neither a micrograph nor a crystal projection.</text>','</svg>'])
(ROOT/'assets/sublattice-order.svg').write_text(''.join(parts),encoding='utf8')
print('Wrote original eight-site occupancy schematic.')
