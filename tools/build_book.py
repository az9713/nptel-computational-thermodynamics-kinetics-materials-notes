"""Build public navigation from public manifests; no private-source dependency."""
from pathlib import Path
from html import escape,unescape
import json,re
ROOT=Path(__file__).resolve().parent.parent
TITLE='Computational Thermodynamics and Kinetics of Materials'
manifest=json.loads((ROOT/'source-manifest.json').read_bytes().decode('utf8'))
chapters=sorted((json.loads(p.read_bytes().decode('utf8')) for p in ROOT.glob('chapter-*.meta.json')),key=lambda c:c['chapter'])
def write_page(name,title,body):
    nav='<a href="index.html">Chapters</a> <a href="source-map.html">Source map</a> <a href="references.html">References</a> <a href="numerical-lab.html">Numerical lab</a>'
    page=('<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><meta name="color-scheme" content="dark light"><title>'
      +escape(title)+' · '+TITLE+'</title><link rel="stylesheet" href="assets/book.css"><script src="assets/math-config.js"></script><script defer src="assets/mathjax-tex-svg.js"></script><script defer src="assets/book.js"></script></head><body>'
      '<a class="skip" href="#main">Skip to content</a><header class="site-header"><a class="brand" href="index.html">'+TITLE+'</a><nav class="course-nav">'+nav+'</nav></header><div class="layout"><aside class="toc" id="chapter-toc"><h2>On this page</h2></aside><main id="main"><div class="reading-tools"><button id="theme-toggle">Light theme</button><button id="text-toggle">Larger text</button><button id="solutions-toggle">Expand solutions</button><button id="print-page">Print with solutions</button></div><article id="chapter-content"><h1>'+escape(title)+'</h1>'+body+'</article></main></div><footer class="site-footer">Independent educational edition. Prof. K. Guruvidyathri, University of Hyderabad; NPTEL coordinating institution IIT Madras. Attribution does not imply endorsement.</footer></body></html>')
    (ROOT/name).write_bytes(page.encode('utf8'))
accepted=sum(c.get('independent_review',{}).get('accepted') is True for c in chapters)
items=''.join('<li><a href="'+escape(c['file'],quote=True)+'">'+escape(c['title'])+'</a> — '+escape(c['status'])+'</li>' for c in chapters)
write_page('index.html',TITLE,
 '<p class="lead">An independent study edition connecting microstructure and phase diagrams to thermodynamic models, reproducible equilibrium calculations, diffusion, precipitation and material design.</p>'
 '<p>The full playlist contains sixty substantive lessons and a separate introduction. This edition has '+str(len(chapters))+' authored chapters; '+str(accepted)+' currently carry two accepted independent scientific reviews at their recorded raw body digests. Public validation records describe the release checks actually completed.</p>'
 '<p>Read entirely offline: mathematics, theme and text controls, navigation, complete expandable solutions and downloadable numerical examples are bundled locally. Code checks test stipulated mathematical models; they do not establish experimental validity.</p>'
 '<h2 id="chapters">Chapters</h2><ol>'+items+'</ol><p><a href="book-overview.html">Study guide</a> · <a href="WORKFLOW.md">Production workflow</a> · <a href="README.md">Reproduction instructions</a></p>')
mapping={r['id']:[c for c in chapters if any(q['video_id']==r['id'] for q in c.get('recordings',[]))] for r in manifest['recordings']}
rows=''
for r in manifest['recordings']:
    targets=mapping[r['id']]
    chlinks='; '.join('<a href="'+escape(c['file'],quote=True)+'">Chapter '+str(c['chapter'])+'</a>' for c in targets) or 'Orientation; no substantive lesson omitted'
    rows+='<tr><td>'+str(r['position'])+'</td><td>'+escape(str(r.get('official_lesson') or 'Introduction'))+'</td><td><a href="'+escape(r['url'],quote=True)+'">'+escape(r['title'])+'</a></td><td>'+chlinks+'</td></tr>'
write_page('source-map.html','Recording source map',
 '<p>All 61 playlist identities retain their original order, including the first four entries before the supplied URL’s starting index. The introduction is orientation; all 60 official lessons are mapped below. The instructor is Prof. K. Guruvidyathri of the University of Hyderabad; IIT Madras coordinates NPTEL distribution.</p>'
 '<p>Full local machine ASRs were read sequentially by each chapter author. Official download metadata did not provide usable captions and YouTube caption requests were rate-limited; local small.en recognition was used. Machine recognition and chunk-boundary errors remain possible. No complete audiovisual or board-equation audit is claimed. Timestamp links identify topic spans; original derivations and teaching examples are distinguished in each chapter.</p>'
 '<table><thead><tr><th>Position</th><th>Lesson</th><th>Recording</th><th>Chapter disposition</th></tr></thead><tbody>'+rows+'</tbody></table>')
links={}
programs=[]
for c in chapters:
    text=(ROOT/c['file']).read_bytes().decode('utf8')
    for url,label in re.findall(r'<a[^>]+href="(https?://[^"]+)"[^>]*>(.*?)</a>',text,flags=re.S):
        if 'youtube.com/' not in url and 'youtu.be/' not in url:
            links.setdefault(unescape(url),re.sub('<[^>]+>','',label))
    for program in re.findall(r'<pre[^>]+data-run="([^"]+)"',text):
        programs.append((c,program))
write_page('references.html','Primary sources and access scope',
 '<p>These are the primary research, official documentation and university teaching links used in the chapters. Each local citation states the supported claim and, where appropriate, the actual pages or abstract inspected. A listed paper is not a claim that its entire text, supplementary material, cited experiments or numerical implementation was audited.</p>'
 '<p>Source access limits, unavailable documents and original synthetic constructions are disclosed in the exposition. Recording links are collected separately in the <a href="source-map.html">source map</a>. No downloaded third-party paper or raw lecture transcript is redistributed.</p><ul>'
 +''.join('<li><a href="'+escape(u,quote=True)+'">'+escape(t)+'</a></li>' for u,t in links.items())+'</ul>')
write_page('numerical-lab.html','Numerical laboratory',
 '<p>Every chapter includes an executable program, actual displayed output and complete worked solutions. The programs test original stipulated models, conservation laws and mathematical derivations; the commercial course examples and proprietary databases are not reproduced.</p>'
 '<h2 id="programs">Chapter programs</h2><ol>'
 +''.join('<li><a href="'+escape(p,quote=True)+'">'+escape(p)+'</a> — <a href="'+escape(c['file'],quote=True)+'">'+escape(c['title'])+'</a></li>' for c,p in programs)
 +'</ol><h2 id="runtime">Reproducible runtime</h2><p>Install the pinned public <a href="requirements.txt">requirements</a> for the complete laboratory. Chapters 14, 15 and 22 execute actual pycalphad 0.11.1 calculations against their supplied synthetic TDB files; each chapter states the engine scope and limitations. Other examples use the Python standard library or the pinned numerical stack.</p>'
 '<h2 id="extension">Regular-solution extension</h2><p><a href="numerical_lab.py">The additional symmetric regular-solution model</a> checks a supporting common tangent, phase-mole balances and spinodal endpoints. It is a teaching model rather than an assessed alloy database.</p><img src="assets/regular-solution.svg" alt="Computed common tangent for a synthetic symmetric regular solution">')
groups=[('Physical metallurgy and phase diagrams',1,5),('Thermodynamic foundations',6,9),('Mixing and global equilibrium',10,12),('CALPHAD models, workflows and assessment',13,18),('Diffusion and transformation kinetics',19,21),('Screening, properties and additional reservoirs',22,24)]
guide=''
for title,lo,hi in groups:
    group=[c for c in chapters if lo<=c['chapter']<=hi]
    guide+='<h2 id="group-'+str(lo)+'">'+title+'</h2><p>'+', '.join('<a href="'+escape(c['file'],quote=True)+'">Chapter '+str(c['chapter'])+'</a>' for c in group)+'.</p>'
write_page('book-overview.html','Study guide',
 '<p>Begin with phase identity, amount bases and conservation. Build the thermodynamic potentials and constrained derivatives before using equilibrium engines or kinetic models. The later chapters repeatedly distinguish a positive mobility from a thermodynamic curvature, a stable phase from a kinetic pathway and an experimentally assessed model from a synthetic consistency check.</p>'
 +guide+'<h2 id="practice">How to use this edition</h2><p>Read each derivation with its assumptions, attempt the problems before expanding solutions and then execute the exact downloaded program. Follow topic timestamps when comparing with recordings. Treat source-specific ambiguities and stated access limits as part of the evidence. Scientific review and numerical checks are different forms of verification; neither alone certifies a real material prediction.</p>')
print('Built five navigation pages for',len(chapters),'authored chapters;',accepted,'accepted in two scientific rounds.')
