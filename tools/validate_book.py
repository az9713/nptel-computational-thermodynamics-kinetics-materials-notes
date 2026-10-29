"""Validate local references, frozen content digests and executable examples."""
from pathlib import Path
from html.parser import HTMLParser
from html import unescape
import re
from urllib.parse import urlsplit, unquote
import hashlib
import json
import subprocess
import sys

ROOT=Path(__file__).resolve().parent.parent


class Page(HTMLParser):
    def __init__(self):
        super().__init__();self.links=[];self.ids=[]
    def handle_starttag(self,tag,attrs):
        attrs=dict(attrs)
        if 'id' in attrs:self.ids.append(attrs['id'])
        for key in ('href','src'):
            if key in attrs:self.links.append(attrs[key])


errors=[];checked=0;listings={}
for chapter in ROOT.glob('chapter-*.html'):
    text=chapter.read_bytes().decode('utf8')
    # Rendering can silently accept escape damage as ordinary letters.
    if re.search(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', text):
        errors.append('Forbiddencontrol:'+chapter.name)
    prose=re.sub(r'<pre\b[^>]*>.*?</pre>', '', text, flags=re.S)
    for expression in re.findall(r'\$\$.*?\$\$|\$(?!\$).*?\$', prose, flags=re.S):
        if '\t' in expression or '\r' in expression:
            errors.append('Mathcontrol:'+chapter.name)
        if re.search(r'(?<!\\)\b(?:frac|qquad|theta|langle|rangle|mathrm|operatorname)\b',expression):
            errors.append('BareTeXcommand:'+chapter.name+':'+repr(expression[:90]))
    pattern=r'<pre data-run="([^"]+)"><code[^>]*>(.*?)</code></pre>\s*<pre class="output">(.*?)</pre>'
    for name,code,output in re.findall(pattern,text,re.S):
        listings[name]=(unescape(code).strip().splitlines(),unescape(output).strip().splitlines())
for path in ROOT.glob('*.html'):
    page=Page();page.feed(path.read_text(encoding='utf8'))
    if len(page.ids)!=len(set(page.ids)):errors.append(str(path.name)+':duplicateid')
    for link in page.links:
        u=urlsplit(link)
        if u.scheme or u.netloc or not u.path:continue
        dest=(ROOT/unquote(u.path)).resolve()
        if not dest.is_relative_to(ROOT.resolve()) or not dest.exists():
            errors.append(path.name+':missingorunsafelink:'+link)
        checked+=1
for sidecar in ROOT.glob('chapter-*.meta.json'):
    meta=json.loads(sidecar.read_text(encoding='utf8'))
    body=(ROOT/meta['file']).read_bytes().decode('utf8').split(
        '<!-- BOOK-CONTENT-START -->')[1].split('<!-- BOOK-CONTENT-END -->')[0]
    digest=hashlib.sha256(body.encode()).hexdigest()
    if digest!=meta['authored_body_sha256']:errors.append('Digest:'+sidecar.name)
    for span in meta['sections']:
        duration=next(r['duration_seconds'] for r in meta['recordings'] if r['video_id']==span['video_id'])
        if not 0<=span['start']<span['end']<=duration:errors.append('Timestamp:'+str(span))
for script in ROOT.glob('code/ch*.py'):
    p=subprocess.run([sys.executable,str(script)],cwd=ROOT,capture_output=True,
                     text=True,encoding='utf8',timeout=180)
    if p.returncode:errors.append('Codefailure:'+script.name+':'+p.stderr)
    rel=script.relative_to(ROOT).as_posix()
    if rel not in listings:errors.append('Missinglisting:'+rel)
    else:
        source,expected=listings[rel]
        if source!=script.read_text(encoding='utf8').strip().splitlines():
            errors.append('Listingsourcemismatch:'+rel)
        if expected!=p.stdout.strip().splitlines():
            errors.append('Displayedoutputmismatch:'+rel)
assert not errors,errors
print(json.dumps({'local_references_checked':checked,'errors':errors,
                  'chapter_digests':'pass','timestamps':'pass','code_exits':'pass','displayed_source_and_output':'pass'}))
