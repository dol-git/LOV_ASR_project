"""Global medoid plus deterministic farthest-point sampling; final PDBs only."""
from pathlib import Path
import csv,json,hashlib,warnings
import numpy as np
import MDAnalysis as m
from scipy.spatial.distance import pdist,squareform
warnings.filterwarnings('ignore')
R=Path(__file__).resolve().parent;PRE=R.parent/'Ja_core_visual_20261019'
def read(p):
 assert 'grest' not in str(p).lower()
 with Path(p).open() as f:return list(csv.DictReader(f))
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def table(n,rr):
 with (R/n).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
old=read(PRE/'selected_frames.csv');mp=read(PRE/'residue_correspondence.csv');tr={(x['system'],int(x['frame_0'])):x for x in read(PRE/'all_frame_transforms.csv')};hashes={x['path']:x['sha256'] for x in read(R.parent/'Ja_literature_review_20261019/source_manifest.csv')}
protected={str(p):sha(p) for p in [PRE/'selected_frames.csv',PRE/'generate.py',PRE/'residue_correspondence.csv',PRE/'all_frame_transforms.csv']+[Path(x['pdb']) for x in old]}
(R/'protected_originals.json').write_text(json.dumps(protected,indent=2))
prior8=R.parent/'Ja_core_visual_8_20261019'
protected.update({str(p):sha(p) for p in prior8.rglob('*') if p.is_file()})
(R/'protected_originals.json').write_text(json.dumps(protected,indent=2))
manifest=[];qc=[];selection=[];sources=[];pairwise=[];raw_qc=[]
settings=dict(frames_each=8,frame_pool='HMR production run01, all 1000 stored frames/system, 0.01–10.00 ns',selection='global medoid minimizing sum Ja CA RMSD, fixed; then 7 farthest-point additions maximizing distance to nearest selected frame',distance='sqrt(sum squared common-core-aligned Ja CA differences /23), no Ja refit',tie='earliest original frame index',comparison_optimization='same farthest-point algorithm independently within each system; structural diversity emphasized; no between-system contrast objective',alignment='reuse existing common-core transforms to centered WT frame0; WT11–110 / ANC10–109 CA100; Ja WT120–142 / ANC119–141 CA23',coordinates='actual raw source PROA, existing rigid rotation/translation only; no interpolation or averaging; all stored frames checked for finite coordinates, bond continuity and aligned Ja cache agreement before selection; no PBC repair applied',PDB_output='final 8 per system only; no clustering and no intermediate PDB',meaning='Frames selected to emphasize structural diversity; do not reflect observation frequencies or equilibrium populations')
(R/'selection_definition.json').write_text(json.dumps(settings,indent=2,ensure_ascii=False))
for system in ['WT','ANC']:
 F=np.load(PRE/f'{system}_common_core_Ja_CA.npy');assert F.shape==(1000,23,3)
 source=next(x for x in old if x['system']==system);top=source['source_topology'];dcd=source['source_dcd'];assert 'grest' not in (top+dcd).lower()
 for role,p in [('topology',top),('production',dcd)]:
  h=sha(p);assert hashes[p]==h;sources.append(dict(system=system,run='run01',role=role,path=p,sha256=h,prior_match=True))
 u=m.Universe(top,dcd);ag=u.select_atoms('segid PROA');assert len(u.trajectory)==1000 and abs(u.trajectory.dt-10)<1e-4
 ja=u.select_atoms('segid PROA and name CA and resid '+('120:142' if system=='WT' else '119:141'));jidx=ja.indices
 # Check every raw stored frame before frame selection. Fail rather than silently exclude.
 assert np.isfinite(F).all()
 bonds=ag.bonds.indices
 for ts in u.trajectory:
  f=ts.frame;raw=ag.positions.astype(float)
  finite=bool(np.isfinite(raw).all());lengths=np.linalg.norm(u.atoms[bonds[:,0]].positions-u.atoms[bonds[:,1]].positions,axis=1)
  t=tr[system,f];c=np.array(json.loads(t['center']));rot=np.array(json.loads(t['rotation']))
  jerr=float(np.max(abs((ja.positions-c)@rot-F[f])))
  maxbond=float(lengths.max());minbond=float(lengths.min())
  passed=finite and minbond>.5 and maxbond<3 and jerr<1e-5
  raw_qc.append(dict(system=system,frame_0=f,frame_1=f+1,time_ns=(f+1)*.01,finite=finite,min_covalent_bond_A=minbond,max_covalent_bond_A=maxbond,Ja_cache_error_A=jerr,status='PASS' if passed else 'FAIL'))
  assert passed,('Raw coordinate/PBC/cache anomaly; stop without excluding frame',system,f)
 D=squareform(pdist(F.reshape(1000,-1))/np.sqrt(23));central=int(np.argmin(D.sum(1)));chosen=[central];near=D[:,central].copy();objectives=[float(near.sum())];marginals=[None]
 for k in range(7):
  candidate=near.copy();candidate[chosen]=-np.inf;pick=int(np.argmax(candidate));marginals.append(float(near[pick]));chosen.append(pick);near=np.minimum(near,D[:,pick]);objectives.append(float(near.sum()))
 assert len(set(chosen))==8
 for a,fa in enumerate(chosen):
  for b,fb in enumerate(chosen):
   pairwise.append(dict(system=system,rank_i=a+1,rank_j=b+1,frame1_i=fa+1,frame1_j=fb+1,time_ns_i=(fa+1)*.01,time_ns_j=(fb+1)*.01,Ja_CA_RMSD_A=float(D[fa,fb])))
 for rank,f in enumerate(chosen):
  ts=u.trajectory[f];raw=ag.positions.astype(float).copy();bonds=ag.bonds.indices;maxbond=float(np.linalg.norm(u.atoms[bonds[:,0]].positions-u.atoms[bonds[:,1]].positions,axis=1).max());assert maxbond<3
  t=tr[system,f];c=np.array(json.loads(t['center']));rot=np.array(json.loads(t['rotation']));assert abs(np.linalg.det(rot)-1)<1e-10
  jerr=float(np.max(abs((u.atoms[jidx].positions-c)@rot-F[f])));assert jerr<1e-5
  xyz=(raw-c)@rot;ag.positions=xyz;ts.dimensions=None
  out=R/'final_pdb'/f'{system}_{rank+1:02d}_f{f+1:04d}_t{(f+1)*.01:.2f}ns.pdb'
  with m.Writer(str(out),n_atoms=len(ag),bonds='conect',reindex=True) as w:w.write(ag)
  ck=m.Universe(str(out));err=float(np.max(abs(ck.atoms.positions-xyz)));inv=float(np.max(abs(ck.atoms.positions@rot.T+c-raw)))
  assert err<.00052 and inv<.0009 and np.array_equal(ck.atoms.names,ag.names) and np.array_equal(ck.atoms.resnames,ag.resnames) and np.array_equal(ck.atoms.resids,ag.resids)
  assert len(ck.select_atoms('resname FMC'))==61 and np.isfinite(ck.atoms.positions).all()
  assert set(ck.select_atoms('name CA').chainIDs)=={'X'}
  manifest.append(dict(system=system,run='run01',rank=rank+1,role='central_medoid' if rank==0 else 'additional_exemplar',source_topology=top,source_dcd=dcd,source_log=source['source_log'] if 'source_log' in source else 'see original source_selection.md',frame_0=f,frame_1=f+1,production_step=(f+1)*2500,production_time_ns=(f+1)*.01,method=settings['selection'],tie_break=settings['tie'],nearest_distance_sum_A=objectives[rank],nearest_previous_selected_RMSD_A='' if rank==0 else marginals[rank],alignment_source=str(PRE/'all_frame_transforms.csv'),core_RMSD_A=t['core_RMSD_A'],pdb=str(out),chain='X',opacity=1 if rank==0 else .60))
  qc.append(dict(system=system,pdb=str(out),atoms=len(ck.atoms),FMC_atoms=61,max_PDB_rounding_error_A=err,inverse_rigid_error_A=inv,Ja_cache_error_A=jerr,raw_max_bond_A=maxbond,status='PASS'))
 for f in range(1000):selection.append(dict(system=system,frame_0=f,frame_1=f+1,time_ns=(f+1)*.01,mean_RMSD_to_other_frames_A=float(D[f].sum()/999),nearest_selected_RMSD_A=float(near[f]),selected_rank=chosen.index(f)+1 if f in chosen else 0))
 print(system,'selected frame1',[f+1 for f in chosen],'times',[(f+1)*.01 for f in chosen],flush=True)
table('raw_all_frame_QC.csv',raw_qc);table('selected_pairwise_distances.csv',pairwise);table('selected_frames.csv',manifest);table('validation.csv',qc);table('selection_scores.csv',selection);table('sources.csv',sources)
(R/'chimera_load_list.txt').write_text('\n'.join(x['pdb'] for x in manifest)+'\n')
assert len(list((R/'final_pdb').glob('*.pdb')))==16
for p,h in protected.items():assert sha(p)==h
(R/'verification.json').write_text(json.dumps(dict(status='PASS',PDB_count=16,per_system=8,original_10_PDB_preserved=True,raw_frames_checked=len(raw_qc),all_raw_frames_passed=True,original_files_hash_unchanged=len(protected),coordinate_policy='actual frame rigid transform only',intermediate_PDB_count=0),indent=2))
print('PASS: exactly 16 final PDBs; previous outputs preserved; all 2000 raw frames checked; medoid plus FPS')
