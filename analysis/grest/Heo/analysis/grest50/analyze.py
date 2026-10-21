"""50 ns gREST inventory, exchange QC, and physical-state conversion only.

All outputs are confined to this directory. Raw files are opened read-only.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import hashlib
import json
import os
import re
import struct
import subprocess
import sys
import numpy as np

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[2]
SNAP = ROOT / 'Fugaku_Heo_20260928'
TEMPS = [298.15, 331.0, 365.0, 400.0]
NSTEPS = 12500000
PERIOD = 2500
DT = 0.004
AKMA_PS = 0.0488882129


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')


def save_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        while chunk := f.read(8 * 1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


class DCD:
    def __init__(self, path):
        self.path = Path(path)
        self.f = self.path.open('rb')
        magic = self.f.read(4)
        self.endian = '<' if magic == struct.pack('<i', 84) else '>'
        require(struct.unpack(self.endian + 'i', magic)[0] == 84, f'DCD header invalid: {path}')
        self.f.seek(0)
        h = self.record()
        require(h[:4] == b'CORD', f'Not coordinate DCD: {path}')
        values = struct.unpack(self.endian + '20i', h[4:])
        self.nheader, self.start, self.stride, self.last = values[:4]
        require(values[8] == 0 and values[11] == 0, f'Unsupported fixed atoms/4D DCD: {path}')
        self.delta = struct.unpack(self.endian + 'f', h[40:44])[0]
        self.box = bool(values[10])
        self.title = self.record()
        self.natoms = struct.unpack(self.endian + 'i', self.record())[0]
        self.offsets = []
        size = self.path.stat().st_size
        while self.f.tell() < size:
            self.offsets.append(self.f.tell())
            if self.box:
                self.record(expected=48, skip=True)
            for _ in range(3):
                self.record(expected=4 * self.natoms, skip=True)
        require(self.f.tell() == size, f'DCD byte/frame mismatch: {path}')
        require(len(self.offsets) == self.nheader, f'DCD frame/header mismatch: {path}')

    def record(self, expected=None, skip=False):
        marker = self.f.read(4)
        require(len(marker) == 4, f'Truncated record: {self.path}')
        n = struct.unpack(self.endian + 'i', marker)[0]
        require(n >= 0 and (expected is None or n == expected), f'DCD atom/record mismatch: {self.path}')
        if skip:
            self.f.seek(n, 1)
            data = None
        else:
            data = self.f.read(n)
            require(len(data) == n, f'Truncated payload: {self.path}')
        require(self.f.read(4) == marker, f'DCD record marker mismatch: {self.path}')
        return data

    def frame(self, index):
        self.f.seek(self.offsets[index])
        box = self.record() if self.box else None
        xyz = [np.frombuffer(self.record(), dtype=self.endian + 'f4') for _ in range(3)]
        return box, xyz

    def info(self):
        return dict(path=str(self.path), atom_count=self.natoms, frame_count=len(self.offsets),
                    header_start_step=self.start, header_stride=self.stride, header_last_step=self.last,
                    header_delta_akma=self.delta, header_timestep_ps=self.delta * AKMA_PS,
                    header_frame_interval_ps=self.delta * AKMA_PS * self.stride,
                    unit_cell=self.box, byte_size=self.path.stat().st_size)


def pdb_atoms(path):
    return [l for l in path.read_text().splitlines() if l.startswith(('ATOM  ', 'HETATM'))]


def pdb_key(line):
    return (line[72:76].strip(), int(line[22:26]), line[17:21].strip(), line[12:16].strip())


def psf_atoms(path):
    lines = path.read_text().splitlines()
    i = next(i for i, l in enumerate(lines) if '!NATOM' in l)
    return [l.split() for l in lines[i + 1:i + 1 + int(lines[i].split()[0])]]


def paths(system):
    raw = SNAP / ('HMR_' + system)
    suffix = '_A' if system == 'ANC' else ''
    return raw, SNAP / 'input' / f'ionized_fmc{suffix}.pdb', SNAP / 'input' / f'ionized_fmc{suffix}.psf'


def prepare():
    inventory, summaries, dcd_rows = [], {}, []
    for system in ['WT', 'ANC']:
        raw, pdb, psf = paths(system)
        qdir = OUT / 'qc' / system
        job = '51762878' if system == 'WT' else '51762851'
        stdout = raw / f'output.{job}/0/1/stdout.1.0'
        files = [(kind, r, raw / f'grest_50ns_4r_{system}{r}.{kind}')
                 for kind in ['dcd', 'rem', 'log', 'rst'] for r in range(1, 5)]
        files += [('pdb', 0, pdb), ('psf', 0, psf), ('stdout', 0, stdout),
                  ('input', 0, raw / f'INP_gREST_pro_{system}')]
        for kind, r, p in files:
            require(p.is_file() and p.stat().st_size > 0, f'Missing/empty raw file: {p}')
            st = p.stat()
            inventory.append(dict(system=system, kind=kind, replica=r, path=str(p),
                                  bytes=st.st_size, mtime_ns=st.st_mtime_ns, sha256=digest(p)))
        print(system, 'inventory and SHA256 complete', flush=True)
        atoms, topology = pdb_atoms(pdb), psf_atoms(psf)
        require(len(atoms) == len(topology), f'{system} PDB/PSF atom count mismatch')
        require([pdb_key(l) for l in atoms] == [(a[1], int(a[2]), a[3], a[4]) for a in topology],
                f'{system} PDB/PSF atom ordering mismatch')
        selected = [i for i, l in enumerate(atoms) if pdb_key(l)[0] == 'PROA']
        save_csv(qdir / 'atom_mapping.csv', [dict(reduced_index=i + 1, raw_index=j + 1,
                 segment=pdb_key(atoms[j])[0], residue=pdb_key(atoms[j])[1],
                 residue_name=pdb_key(atoms[j])[2], atom_name=pdb_key(atoms[j])[3])
                 for i, j in enumerate(selected)])
        for r in range(1, 5):
            d = DCD(raw / f'grest_50ns_4r_{system}{r}.dcd')
            require(d.natoms == len(atoms) and len(d.offsets) == 5000,
                    f'{system} replica {r} atom/frame mismatch')
            require((d.start, d.stride, d.last) == (PERIOD, PERIOD, NSTEPS), 'DCD step mismatch')
            require(abs(d.delta * AKMA_PS - DT) < 1e-8, 'DCD timestep mismatch')
            dcd_rows.append(dict(system=system, replica=r, **d.info()))
            d.f.close()
        events = defaultdict(list)
        rep_accept, final_logs = [], []
        for r in range(1, 5):
            step, pending, count, accepted = 0, False, 0, 0
            for line in (raw / f'grest_50ns_4r_{system}{r}.log').read_text().splitlines():
                fields = line.split()
                if len(fields) > 2 and fields[0] == 'INFO:' and fields[1].isdigit():
                    step, time_ps = int(fields[1]), float(fields[2])
                if pending:
                    m = re.match(r'\s*(\d+)\s+(\d+)\s+(\d+).*\s([TF])\s*$', line)
                    require(m is not None, 'Unparseable exchange record')
                    dim, a, b, result = m.groups()
                    a, b = int(a), int(b)
                    events[(step, min(a, b), max(a, b))].append((r, a, b, result == 'T'))
                    count += 1
                    accepted += result == 'T'
                    pending = False
                if 'REST Exchange>' in line:
                    pending = True
            require(step == NSTEPS and not pending, f'{system} log incomplete: {r}')
            final_logs.append(dict(replica=r, final_step=step, final_time_ps=time_ps))
            rep_accept.append(dict(replica=r, accepted=accepted, attempts=count, percent=100*accepted/count))
        pair_counts = defaultdict(Counter)
        by_step = defaultdict(list)
        event_rows = []
        for (step, a, b), entries in sorted(events.items()):
            require(len(entries) == 2 and entries[0][3] == entries[1][3]
                    and entries[0][1:3] == entries[1][2:0:-1], 'Exchange duplicate inconsistency')
            accept = entries[0][3]
            pair_counts[(a, b)]['attempts'] += 1
            pair_counts[(a, b)]['accepted'] += accept
            by_step[step].append(entries)
            event_rows.append(dict(step=step, time_ns=step*DT/1000, state_low=a, state_high=b,
                                   replica_a=entries[0][0], replica_b=entries[1][0], accepted=int(accept)))
        require(len(events) == 7500 and all(v['attempts'] == 2500 for v in pair_counts.values()),
                'Exchange trial count mismatch')
        pairs = [dict(state_low=a, state_high=b, temperature_low_K=TEMPS[a-1], temperature_high_K=TEMPS[b-1],
                      **v, percent=100*v['accepted']/v['attempts']) for (a,b),v in sorted(pair_counts.items())]
        rems = []
        for r in range(1, 5):
            rr = [tuple(map(int, l.split())) for l in (raw / f'grest_50ns_4r_{system}{r}.rem').read_text().splitlines() if l.strip()]
            require([s for s, state in rr] == list(range(0, NSTEPS, PERIOD)), 'REM frame/step mismatch')
            rems.append(rr)
        trace = [[rems[r][i][1] for r in range(4)] for i in range(5000)]
        for i in range(1, 5001):
            previous = trace[i-1]
            following = previous.copy()
            for entries in by_step[i*PERIOD]:
                for r, a, b, accepted in entries:
                    require(previous[r-1] == a, 'Log/REM state mismatch')
                    if accepted:
                        following[r-1] = b
            require(sorted(following) == [1,2,3,4], 'State permutation mismatch')
            if i < 5000:
                require(following == trace[i], 'Log/REM transition mismatch')
            else:
                trace.append(following)
        occupancy = [dict(replica=r+1, state=s, temperature_K=TEMPS[s-1],
                         intervals=sum(row[r] == s for row in trace[:-1]),
                         time_ns=sum(row[r] == s for row in trace[:-1])*0.01,
                         percent=sum(row[r] == s for row in trace[:-1])/50)
                     for r in range(4) for s in range(1,5)]
        trips, rt_summary, censored = [], [], []
        for r in range(4):
            anchor, start, opposite = None, None, None
            for i, row in enumerate(trace):
                state = row[r]
                if anchor is None and state in (1,4):
                    anchor, start = state, i
                elif anchor is not None and opposite is None and state == 5-anchor:
                    opposite = i
                elif opposite is not None and state == anchor:
                    trips.append(dict(replica=r+1, direction=f'{anchor}-{5-anchor}-{anchor}',
                                      start_ns=start*0.01, opposite_ns=opposite*0.01, end_ns=i*0.01,
                                      duration_ns=(i-start)*0.01))
                    start, opposite = i, None
            durations = [t['duration_ns'] for t in trips if t['replica'] == r+1]
            rt_summary.append(dict(replica=r+1, count=len(durations), min_ns=min(durations) if durations else None,
                                   median_ns=float(np.median(durations)) if durations else None,
                                   mean_ns=float(np.mean(durations)) if durations else None,
                                   max_ns=max(durations) if durations else None))
            censored.append(dict(replica=r+1, anchor=anchor, start_ns=start*0.01 if start is not None else None,
                                 opposite_ns=opposite*0.01 if opposite is not None else None, end_ns=50.0))
        save_csv(qdir/'exchange_events.csv', event_rows)
        save_csv(qdir/'exchange_acceptance.csv', pairs)
        save_csv(qdir/'replica_acceptance.csv', rep_accept)
        save_csv(qdir/'parameter_state_occupancy.csv', occupancy)
        save_csv(qdir/'parameter_state_trace.csv', [dict(step=i*PERIOD, time_ns=i*0.01,
                  replica_1=row[0], replica_2=row[1], replica_3=row[2], replica_4=row[3],
                  source='REM' if i < 5000 else 'log_final_exchange') for i,row in enumerate(trace)])
        save_csv(qdir/'round_trips.csv', trips, ['replica','direction','start_ns','opposite_ns','end_ns','duration_ns'])
        save_csv(qdir/'round_trip_summary.csv', rt_summary)
        save_csv(qdir/'round_trip_censored.csv', censored)
        durations = [t['duration_ns'] for t in trips]
        hist, bins = np.histogram(durations, bins=np.arange(0, 52, 2))
        save_csv(qdir/'round_trip_distribution.csv', [dict(left_ns=float(a), right_ns=float(b), count=int(n))
                  for a,b,n in zip(bins[:-1],bins[1:],hist)])
        summaries[system] = dict(raw_atoms=len(atoms), reduced_atoms=len(selected), frames=5000,
                  frame_interval_ps=10.0, duration_ns=50.0, first_frame_ns=0.01, last_frame_ns=50.0,
                  analysis_jalpha='119-141' if system == 'ANC' else '120-142',
                  core_fitting_region='1-118' if system == 'ANC' else '1-119',
                  production_solute='120-142', fitting_applied=False,
                  exchange=pairs, overall_acceptance_percent=100*sum(v['accepted'] for v in pair_counts.values())/7500,
                  round_trip_summary=rt_summary, round_trip_total=len(trips),
                  round_trip_quantiles_ns={str(q):float(np.quantile(durations,q)) for q in [0,.25,.5,.75,1]} if durations else {},
                  final_logs=final_logs, log_rem_mismatches=0)
        work = OUT/'work'/system
        work.mkdir(parents=True, exist_ok=True)
        (OUT/'trajectories'/system).mkdir(parents=True, exist_ok=True)
        relative = lambda p: os.path.relpath(p, work)
        text = f'''[INPUT]
psffile = {relative(psf)}
reffile = {relative(pdb)}
dcdfile = {relative(raw / ('grest_50ns_4r_' + system + '{}.dcd'))}
remfile = {relative(raw / ('grest_50ns_4r_' + system + '{}.rem'))}

[OUTPUT]
pdbfile = ../../trajectories/{system}/physical_state.pdb
trjfile = ../../trajectories/{system}/physical_state_{{}}.dcd

[SELECTION]
group1 = sid:PROA

[FITTING]
fitting_method = NO

[OPTION]
check_only      = NO
convert_type    = PARAMETER
convert_ids     = 1
num_replicas    = 4
nsteps          = 12500000
exchange_period = 2500
crdout_period   = 2500
trjout_period   = 2500
trjout_format   = DCD
trjout_type     = COOR+BOX
trjout_atom     = 1
centering      = NO
pbc_correct    = NO
'''
        (work/'remd_convert.inp').write_text(text)
        print(system, 'raw QC passed', json.dumps(summaries[system]), flush=True)
    save_csv(OUT/'inventory.csv', inventory)
    save_csv(OUT/'qc'/'raw_dcd_headers.csv', dcd_rows)
    save_json(OUT/'qc'/'summary.json', summaries)


def convert(system):
    work = OUT/'work'/system
    executable = ROOT/'tools/genesis-2.1.6/bin/remd_convert'
    env = dict(os.environ, OMP_NUM_THREADS='1', TMPDIR=str(OUT/'tmp'))
    (OUT/'tmp').mkdir(exist_ok=True)
    command = [str(executable), 'remd_convert.inp']
    save_json(work/'command.json', dict(command=command, cwd=str(work), OMP_NUM_THREADS='1',
                                      executable_sha256=digest(executable)))
    with (work/'remd_convert.log').open('w') as f:
        result = subprocess.run(command, cwd=work, env=env, stdout=f, stderr=subprocess.STDOUT)
    require(result.returncode == 0, f'GENESIS command failed ({result.returncode}); stop: {work}/remd_convert.log')
    log = (work/'remd_convert.log').read_text()
    require('[STEP4] Deallocate memory' in log and not re.search(r'(?im)^\s*(?:Error>|FATAL|Fortran runtime error)', log),
            f'GENESIS reported failure; stop: {work}/remd_convert.log')
    print(system, 'GENESIS completed', flush=True)


def verify(system):
    raw, pdb, psf = paths(system)
    target = OUT/'trajectories'/system
    original = pdb_atoms(pdb)
    selected = [i for i,l in enumerate(original) if pdb_key(l)[0] == 'PROA']
    reduced = pdb_atoms(target/'physical_state.pdb')
    require([pdb_key(l) for l in reduced] == [pdb_key(original[i]) for i in selected],
            'Reduced reference atom count/order mismatch; stop')
    output = DCD(target/'physical_state_1.dcd')
    require(output.natoms == len(selected) and len(output.offsets) == 5000, 'Reduced DCD atom/frame mismatch; stop')
    inputs = [DCD(raw/f'grest_50ns_4r_{system}{r}.dcd') for r in range(1,5)]
    rems = [[tuple(map(int,l.split())) for l in (raw/f'grest_50ns_4r_{system}{r}.rem').read_text().splitlines() if l.strip()] for r in range(1,5)]
    mapping = []
    for i in range(5000):
        replica = next(r for r in range(4) if rems[r][i][1] == 1)
        source_box, source_xyz = inputs[replica].frame(i)
        out_box, out_xyz = output.frame(i)
        require(all(np.array_equal(a[selected], b) for a,b in zip(source_xyz,out_xyz)),
                f'Physical-state coordinates/atom order mismatch at frame {i+1}; stop')
        require(np.array_equal(np.frombuffer(source_box, dtype=inputs[replica].endian+'f8'),
                               np.frombuffer(out_box,dtype=output.endian+'f8')), 'Box mismatch; stop')
        mapping.append(dict(frame=i+1, step=(i+1)*PERIOD, time_ps=(i+1)*10.0, time_ns=(i+1)*0.01,
                            source_replica=replica+1, rem_step=i*PERIOD, parameter_state=1, temperature_K=298.15))
    save_csv(OUT/'qc'/system/'physical_frame_mapping.csv', mapping)
    report = dict(**output.info(), physical_temperature_K=298.15, reference_pdb=str(target/'physical_state.pdb'),
                  expected_frame_interval_ps=10.0, simulation_duration_ns=50.0,
                  first_frame_time_ns=0.01, last_frame_time_ns=50.0, first_to_last_span_ns=49.99,
                  all_frames_coordinate_and_atom_order_exact_match=True, all_boxes_match=True,
                  reduced_pdb_psf_order_match=True, fitting='NO', pbc_correct='NO',
                  pdb_sha256=digest(target/'physical_state.pdb'), dcd_sha256=digest(target/'physical_state_1.dcd'))
    save_json(OUT/'qc'/system/'physical_trajectory_qc.json', report)
    for d in inputs + [output]:
        d.f.close()
    print(system, 'physical QC passed', json.dumps(report), flush=True)


def check_raw_unchanged():
    with (OUT/'inventory.csv').open() as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        p = Path(row['path'])
        require(p.stat().st_size == int(row['bytes']) and p.stat().st_mtime_ns == int(row['mtime_ns'])
                and digest(p) == row['sha256'], f'Raw file changed: {p}')
    save_json(OUT/'qc'/'raw_immutability.json', dict(files_checked=len(rows), sha256_unchanged=True,
                                                    sizes_and_mtimes_unchanged=True))
    print('All raw input SHA256/size/mtime unchanged', flush=True)


if __name__ == '__main__':
    try:
        if sys.argv[1] == 'prepare':
            prepare()
        elif sys.argv[1] == 'convert':
            convert(sys.argv[2])
            verify(sys.argv[2])
        elif sys.argv[1] == 'immutability':
            check_raw_unchanged()
        else:
            raise ValueError('Unknown action')
    except Exception as exc:
        save_json(OUT/'STOP_ERROR.json', dict(error=str(exc), action=sys.argv[1:]))
        raise
