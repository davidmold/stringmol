#!/usr/bin/env python3
"""Build a self-contained viewer page for a Stringmol run directory.
Usage: build.py RUN_DIR [OUT.html] [EVERY]   (default OUT = RUN_DIR/viewer.html, EVERY = 500 steps)
Each frame: one byte per grid cell. 0 = empty; 1..3 = class of a non-top species
(replicase / binding-only fragment / other); 4+ = one of the top-N lineages."""
import sys, os, re, json, base64, collections
run = sys.argv[1].rstrip('/')
out = sys.argv[2] if len(sys.argv) > 2 else os.path.join(run, 'viewer.html')
every = int(sys.argv[3]) if len(sys.argv) > 3 else 500
NTOP = 12

def cls(s):
    if '$' not in s: return 2                      # binding region only, no code
    pre = s[:s.find('$')]
    if '=?>' in s and '%}' in s and len(pre) >= 15: return 1   # replicase
    return 3                                       # other (code-only, broken, run-ons)

def slim_frames():
    # slim output: species.tsv (id -> sequence) + grid.snap (zlib frames of uint32 species ids)
    import struct, zlib
    seqs = {}
    for line in open(os.path.join(run, 'species.tsv')):
        a = line.rstrip('\n').split('\t')
        if len(a) >= 5: seqs[int(a[0])] = a[4]
    out = []; gx = gy = None
    with open(os.path.join(run, 'grid.snap'), 'rb') as f:
        while True:
            h = f.read(16)
            if len(h) < 16: break
            st, gx, gy, n = struct.unpack('<4I', h); z = f.read(n)
            if len(z) < n: break                     # partial frame while the run is writing
            if st % every: continue
            ids = struct.unpack('<%dI' % (gx * gy), zlib.decompress(z))
            out.append((st, [(i % gx, i // gx, seqs.get(v & 0x7fffffff, '')) for i, v in enumerate(ids) if v]))
    return out, gx, gy

steps = sorted(int(m.group(1)) for f in os.listdir(run) if (m := re.match(r'out1_(\d+)\.conf$', f)))
steps = [s for s in steps if s % every == 0]
frames = []
gx = gy = None
abund = collections.Counter(); first = {}
slim = os.path.exists(os.path.join(run, 'grid.snap'))
if slim:
    loaded, gx, gy = slim_frames(); steps = []
else:
    loaded = None
for s in (steps if not slim else []):
    mols = []; seq = None
    for line in open(os.path.join(run, f'out1_{s:05d}.conf' if s < 100000 else f'out1_{s}.conf')):
        if line.startswith('GRIDX'): gx = int(line.split()[1])
        elif line.startswith('GRIDY'): gy = int(line.split()[1])
        elif line.startswith('AGENT'): seq = line.split()[1]
        elif line.startswith('GRIDPOS') and seq is not None:
            _, x, y = line.split(); mols.append((int(x), int(y), seq)); seq = None
        # bound pairs are written as reaction blocks; the '#####' line carries the
        # active molecule's species sequence (its REACTION_ACTIVE string is mid-copy)
        elif line.startswith('###############'): seq = line.split()[2]
        elif line.startswith('REACTION_ACTIVE') or line.startswith('REACTION_PASSIVE'):
            p = line.split(); gi = p.index('grid:'); x, y = int(p[gi + 1]), int(p[gi + 2])
            q = seq if p[0] == 'REACTION_ACTIVE' else p[2]
            mols.append((x, y, q)); seq = None
    frames.append((s, mols))
    if not mols: break          # extinct: stop at the first empty frame
    for _, _, q in mols:
        abund[q] += 1; first.setdefault(q, s)
if slim:
    for s, mols in loaded:
        frames.append((s, mols))
        for _, _, q in mols:
            abund[q] += 1; first.setdefault(q, s)
        if not mols: break
top = [q for q, _ in abund.most_common(NTOP)]
topidx = {q: i + 4 for i, q in enumerate(top)}
enc, series = [], []
for s, mols in frames:
    g = bytearray(gx * gy)
    c = collections.Counter(); tc = collections.Counter(); lens = []
    for x, y, q in mols:
        k = cls(q); c[k] += 1
        if k == 1: lens.append(len(q))
        v = topidx.get(q, k); g[y * gx + x] = v
        if v >= 4: tc[v - 4] += 1
    enc.append(base64.b64encode(bytes(g)).decode())
    series.append({'t': s, 'rep': c[1], 'frag': c[2], 'oth': c[3],
                   'replen': round(sum(lens) / len(lens), 1) if lens else None,
                   'top': [tc[i] for i in range(len(top))]})
data = json.dumps({'name': os.path.basename(os.path.abspath(run)), 'gx': gx, 'gy': gy, 'every': every,
           'species': [{'seq': q, 'cls': cls(q), 'len': len(q), 'first': first[q]} for q in top],
           'series': series, 'frames': enc}, separators=(',', ':'))
tpl = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'viewer.html')).read()
open(out, 'w').write(tpl.replace('/*DATA*/', data.replace('</', '<\\/')))
print(len(frames), 'frames', gx, 'x', gy, '->', out, os.path.getsize(out) // 1024, 'KB')
