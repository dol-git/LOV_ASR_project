"""Bounded HMR reanalysis. Never opens gREST paths and never writes PDB."""
from pathlib import Path
import sys,os,csv,json,gzip,hashlib,warnings
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R/'software'))
import mdtraj as md
warnings.filterwarnings('ignore')
C=R.parent/'Ja_core_visual_20261019';T=R.parent/'Ja_sheet_tilt_20261019';OLD=R.parent/'Heo/analysis/reanalysis_20261006'
def safe(p):
 p=Path(p);assert 'grest' not in str(p).lower();return p

def read(p):
 with safe(p).open() as f:return list(csv.DictReader(f))
def table(name,rows):
 if not rows:return
 with (R/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def sha(p):
 h=hashlib.sha256()
 with safe(p).open('rb') as f:
  for c in iter(lambda:f.read(1048576),b''):h.update(c)
 return h.hexdigest()
mp=read(C/'residue_correspondence.csv');sources=read(C/'selected_frames.csv');trs=read(C/'all_frame_transforms.csv');ident=read(OLD/'00_provenance/atom_identity/atom_identity.csv')
protected=[C/'residue_correspondence.csv',C/'all_frame_transforms.csv',C/'selected_frames.csv',T/'secondary_structure.csv',T/'all_frame_angles.csv',OLD/'03_homology_mapping/mapping_validation.md',R.parent/'Ja_frames_20261019/source_selection.md']+list((C/'final_pdb').glob('*.pdb'))+list((T/'final_pdb').glob('*.pdb'))
checks={str(p):sha(p) for p in protected};(R/'protected_before.json').write_text(json.dumps(checks,indent=2))
# No metric-based tuning; fixed before new contact calculations.
settings=dict(systems=['WT','ANC'],method='HMR NPT 298.15 K; 4 fs; one production run/system',frames=1000,window_ns=[.01,10.],halves_ns=[[.01,5.],[5.01,10.]],contact_core_pairs=100,Ja_pairs=23,contact_thresholds_A=[4.,4.5,5.],primary_threshold_A=4.5,primary_atoms='validated protein heavy atoms, including FMC cysteine portion only',atom_selection_sensitivity='backbone N CA C O only at same distance thresholds; geometric control, not equivalent side-chain interaction test',union='any frame with protein-heavy minimum distance <5 A in either system; all 2300 pairs computed',exclusions='water, ions, FMC flavin portion, nonmapped residues; no frame exclusion',pair_sequence_separation_exclusion='abs(WT core-Ja resid)<4; none of current 2300 pairs',stable_contact_rule='primary full |WT-Anc fraction| >=0.20; all 4/4.5/5 A full differences same sign with |difference|>=0.10; both halves at4.5 A same sign and |difference|>=0.10',distance_display_rule='top3 union pairs by symmetric absolute full4.5 A difference; ties coreWT then JaWT; display regardless robustness',candidate_rule='at most5 changed direct-MD corresponding sites outside WT119-124; touch a stable contact; select each site max stable |contact difference|; reciprocal WT<->Anc substitutions; exploratory only; original ASR provenance unresolved',DSSP='new MDTraj 8-state, H only alpha; old C3 and alpha4 motif reused as diagnostic crosschecks, cannot replace alpha-specific call',RMSF='same existing common-core transforms; unweighted CA fluctuation about own full/half mean; no half refit',angles='not primary; no angle search',new_PDB=False,mdtraj=md.__version__,MDAnalysis=mda.__version__,numpy=np.__version__)
(R/'definitions.json').write_text(json.dumps(settings,indent=2))
# Discovery is name-only and prunes excluded directories before visiting them.
inventory=[]
for root in [R.parent/'Heo',R.parent/'Fugaku_Heo_20260928']:
 for base,dirs,files in os.walk(root,followlinks=False):
  dirs[:]=[d for d in dirs if 'grest' not in d.lower() and d not in ['analysis','tools','regression_test'] and not Path(base,d).is_symlink()]
  for n in files:
   if 'grest' not in n.lower() and n.endswith('.dcd'):inventory.append(dict(path=str(Path(base,n)),method='HMR' if 'hmr' in str(Path(base,n)).lower() else 'unclassified',role='equilibration' if '_eq' in n else 'production',independent_run='run01' if '_prod' in n else 'not independent'))
table('trajectory_inventory.csv',inventory)
windows=[('full',0,1000),('first_half',0,500),('second_half',500,1000)];data={};provenance=[];atomrows=[];sserows=[];ssefrac=[];rmsfrows=[];placement=[];exclusions=[]
coremaps=[x for x in mp if x['region']=='core'];jamaps=[x for x in mp if x['region']=='Ja'];assert len(coremaps)==100 and len(jamaps)==23
pairrows=[]
for j,jm in enumerate(jamaps):
 for c,cm in enumerate(coremaps):pairrows.append(dict(pair_id=f"c{cm['WT_resid']}_j{jm['WT_resid']}",core_index=c,Ja_index=j,WT_core=int(cm['WT_resid']),ANC_core=int(cm['ANC_resid']),WT_core_aa=cm['WT_aa'],ANC_core_aa=cm['ANC_aa'],WT_Ja=int(jm['WT_resid']),ANC_Ja=int(jm['ANC_resid']),WT_Ja_aa=jm['WT_aa'],ANC_Ja_aa=jm['ANC_aa'],mapping_scope='direct MD pairwise correspondence; original ASR mapping unresolved',boundary_unresolved=119<=int(jm['WT_resid'])<=124))
heavycounts={}
for system,sname in [('WT','Extant'),('ANC','Ancestral')]:
 source=next(x for x in sources if x['system']==system);top=safe(source['source_topology']);dcd=safe(source['source_dcd']);u=mda.Universe(str(top),str(dcd));ag=u.select_atoms('segid PROA');assert len(u.trajectory)==1000 and abs(u.trajectory.dt-10)<1e-4
 verified={x['path']:x['sha256'] for x in read(T/'source_verification.csv')}
 for role,p in [('topology',top),('production',dcd)]:
  h=sha(p);assert h==verified[str(p)];provenance.append(dict(system=system,method='HMR',run='run01',role=role,path=str(p),sha256=h,reused_prior_hash=True))
 ii=[x for x in ident if x['system']==sname];assert len(ii)==len(u.atoms)
 for x,a in zip(ii,u.atoms):assert int(x['index'])==a.index+1 and x['atom_name']==a.name and x['resname']==a.resname
 heavy={int(x['index'])-1 for x in ii if x['segment']=='PROA' and x['is_hydrogen']=='False' and x['portion'] in ['protein','protein_cysteine']}
 residues={int(r.resid):r for r in ag.residues};ca={int(a.resid):a.index for a in ag.select_atoms('name CA')}
 cr=[int(x[system+'_resid']) for x in coremaps];jr=[int(x[system+'_resid']) for x in jamaps];allr=cr+jr
 for r in residues:
  if r not in allr:exclusions.append(dict(system=system,segid='PROA',resid=r,resname=residues[r].resname,reason='outside reused core100 and Ja23 correspondence; not forced into new contact/RMSF map'))
 groups={}
 for selection in ['protein_heavy','backbone_heavy']:
  coregroups=[np.array([a.index for a in residues[i].atoms if a.index in heavy and (selection=='protein_heavy' or a.name in ['N','CA','C','O'])]) for i in cr]
  jagroups=[np.array([a.index for a in residues[i].atoms if a.index in heavy and (selection=='protein_heavy' or a.name in ['N','CA','C','O'])]) for i in jr]
  assert all(len(x)>0 for x in coregroups+jagroups)
  ci=np.concatenate(coregroups);ji=np.concatenate(jagroups);co=np.r_[0,np.cumsum([len(x) for x in coregroups])[:-1]];jo=np.r_[0,np.cumsum([len(x) for x in jagroups])[:-1]];groups[selection]=(ci,ji,co,jo)
  heavycounts[system,selection]=([len(x) for x in coregroups],[len(x) for x in jagroups])
  for role,rr,gg in [('core',cr,coregroups),('Ja',jr,jagroups)]:
   for resid,g in zip(rr,gg):
    for idx in g:
     a=u.atoms[idx];atomrows.append(dict(system=system,selection=selection,role=role,psf_index_0=int(idx),segid=a.segid,resid=resid,resname=a.resname,name=a.name))
 # Explicit analysis-only backbone topology preserves the FMC cysteine chain position.
 bb=np.array([[res.atoms.select_atoms('name '+n)[0].index for n in ['N','CA','C','O OT1 O1']] for res in ag.residues]);mt=md.Topology();ch=mt.add_chain();at=[]
 for res in ag.residues:
  rn={'FMC':'CYS','HSD':'HIS','HSE':'HIS','HSP':'HIS'}.get(res.resname,res.resname);mr=mt.add_residue(rn,ch,resSeq=int(res.resid));aa=[]
  for name,element in [('N',md.element.nitrogen),('CA',md.element.carbon),('C',md.element.carbon),('O',md.element.oxygen)]:aa.append(mt.add_atom(name,element,mr))
  mt.add_bond(aa[0],aa[1]);mt.add_bond(aa[1],aa[2]);mt.add_bond(aa[2],aa[3])
  if at:mt.add_bond(at[-1][2],aa[0])
  at.append(aa)
 coords=np.empty((1000,len(allr),3));bbcoords=np.empty((1000,len(ag.residues)*4,3),dtype=np.float32);dist={s:np.empty((1000,23,100)) for s in groups};tr=[x for x in trs if x['system']==system];assert len(tr)==1000
 for ts in u.trajectory:
  f=ts.frame;row=tr[f];assert int(row['frame_0'])==f
  raw=u.atoms.positions.astype(float);cc=np.array(json.loads(row['center']));rot=np.array(json.loads(row['rotation']));assert abs(np.linalg.det(rot)-1)<1e-9
  coords[f]=(raw[[ca[i] for i in allr]]-cc)@rot;bbcoords[f]=raw[bb.ravel()]/10
  for selection,(ci,ji,co,jo) in groups.items():
   dd=distance_array(raw[ji],raw[ci]);dist[selection][f]=np.minimum.reduceat(np.minimum.reduceat(dd,co,axis=1),jo,axis=0)
 np.savez_compressed(R/f'{system}_frame_arrays.npz',aligned_CA=coords,backbone_nm=bbcoords,protein_heavy_distances_A=dist['protein_heavy'],backbone_heavy_distances_A=dist['backbone_heavy'])
 cache=np.load(C/f'{system}_common_core_Ja_CA.npy');cache_err=float(np.max(abs(coords[:,100:]-cache)));assert cache_err<1e-5
 sse=md.compute_dssp(md.Trajectory(bbcoords,mt),simplified=False);assert sse.shape==(1000,len(ag.residues)) and not np.isin(sse,['NA']).any();np.save(R/f'{system}_DSSP8.npy',sse)
 jss=sse[:,np.array(jr)-1];data[system]=dict(coords=coords,dist=dist,sse=jss,cr=cr,jr=jr,cache_error=cache_err)
 for f in range(1000):
  for k,j in enumerate(jr):sserows.append(dict(system=system,run='run01',frame_0=f,frame_1=f+1,time_ns=(f+1)*.01,WT_resid=int(jamaps[k]['WT_resid']),ANC_resid=int(jamaps[k]['ANC_resid']),own_resid=j,DSSP8='C' if jss[f,k]==' ' else jss[f,k],alpha_H=bool(jss[f,k]=='H'),boundary_unresolved=119<=int(jamaps[k]['WT_resid'])<=124))
 for window,start,end in windows:
  for k,x in enumerate(jamaps):
   states=jss[start:end,k];ssefrac.append(dict(system=system,window=window,WT_resid=int(x['WT_resid']),ANC_resid=int(x['ANC_resid']),WT_aa=x['WT_aa'],ANC_aa=x['ANC_aa'],alpha_H_fraction=float(np.mean(states=='H')),any_HGI_fraction=float(np.mean(np.isin(states,['H','G','I']))),coil_turn_fraction=float(np.mean(np.isin(states,[' ','T','S']))),frames=end-start))
  q=coords[start:end];rf=np.sqrt(np.mean(np.sum((q-q.mean(0))**2,axis=2),axis=0))
  for k,x in enumerate(coremaps+jamaps):rmsfrows.append(dict(system=system,window=window,region=x['region'],WT_resid=int(x['WT_resid']),ANC_resid=int(x['ANC_resid']),WT_aa=x['WT_aa'],ANC_aa=x['ANC_aa'],RMSF_A=float(rf[k]),frames=end-start))
 ref=coords[0,100:];ref0=ref-ref.mean(0)
 oldangles={int(x['frame_0']):x for x in read(T/'all_frame_angles.csv') if x['system']==system}
 for f,q in enumerate(coords[:,100:]):
  qc=q-q.mean(0);L,_,V=np.linalg.svd(qc.T@ref0);rot=L@np.diag([1,1,np.linalg.det(L@V)])@V
  placement.append(dict(system=system,frame_0=f,time_ns=(f+1)*.01,alpha_H_fraction=float(np.mean(jss[f]=='H')),core_aligned_Ja_RMSD_to_own_frame0_A=float(np.sqrt(np.mean(np.sum((q-ref)**2,axis=1)))),Ja_self_fit_RMSD_to_own_frame0_A=float(np.sqrt(np.mean(np.sum((qc@rot-ref0)**2,axis=1)))),Ja_centroid_displacement_to_own_frame0_A=float(np.linalg.norm(q.mean(0)-ref.mean(0))),Ja_half_bend_deg=float(oldangles[f]['Ja_half_bend_deg'])))
 print(system,'full DSSP computed; prior aligned Ja cache error',cache_err,flush=True)
# All statistics use the same pair grid and symmetric union.
union=np.any(data['WT']['dist']['protein_heavy']<5,axis=0)|np.any(data['ANC']['dist']['protein_heavy']<5,axis=0)
comparisons=[];freqrows=[]
for p in pairrows:
 j,c=p['Ja_index'],p['core_index'];p['union_any_frame_lt5_A']=bool(union[j,c])
 for system in ['WT','ANC']:
  p[system+'_core_heavy_atoms']=heavycounts[system,'protein_heavy'][0][c];p[system+'_Ja_heavy_atoms']=heavycounts[system,'protein_heavy'][1][j]
 for selection in ['protein_heavy','backbone_heavy']:
  for cutoff in [4.,4.5,5.]:
   for window,start,end in windows:
    a=data['WT']['dist'][selection][start:end,j,c];b=data['ANC']['dist'][selection][start:end,j,c];pa=float(np.mean(a<cutoff));pb=float(np.mean(b<cutoff))
    row=dict(**p,selection=selection,threshold_A=cutoff,window=window,WT_contact_fraction=pa,ANC_contact_fraction=pb,WT_minus_ANC_fraction=pa-pb,WT_median_distance_A=float(np.median(a)),ANC_median_distance_A=float(np.median(b)))
    if p['union_any_frame_lt5_A']:freqrows.append(row)
 # Same symmetric criterion for both directions; no one-system-only ranking.
 diffs={}
 for cutoff in [4.,4.5,5.]:
  for window,start,end in windows:
   diffs[cutoff,window]=float(np.mean(data['WT']['dist']['protein_heavy'][start:end,j,c]<cutoff)-np.mean(data['ANC']['dist']['protein_heavy'][start:end,j,c]<cutoff))
 delta=diffs[4.5,'full'];sign=np.sign(delta)
 stable=abs(delta)>=.2 and all(np.sign(diffs[t,'full'])==sign and abs(diffs[t,'full'])>=.1 for t in [4.,4.5,5.]) and all(np.sign(diffs[4.5,w])==sign and abs(diffs[4.5,w])>=.1 for w in ['first_half','second_half'])
 comparisons.append(dict(**p,delta_full=delta,delta_first=diffs[4.5,'first_half'],delta_second=diffs[4.5,'second_half'],delta_4A=diffs[4.,'full'],delta_5A=diffs[5.,'full'],stable_within_observed_run=bool(stable),candidate_eligible=bool(stable and not p['boundary_unresolved']),backbone_delta_full=float(np.mean(data['WT']['dist']['backbone_heavy'][:,j,c]<4.5)-np.mean(data['ANC']['dist']['backbone_heavy'][:,j,c]<4.5))))
table('all_pair_correspondence.csv',pairrows);table('contact_fractions.csv',freqrows);table('contact_robustness.csv',comparisons);table('contact_atom_selections.csv',atomrows);table('excluded_residues.csv',exclusions);table('source_manifest.csv',provenance);table('Ja_DSSP8_timeseries.csv',sserows);table('Ja_helix_fractions.csv',ssefrac);table('common_core_aligned_RMSF.csv',rmsfrows);table('Ja_placement_vs_deformation.csv',placement)
# All observed union pairs, not just the displayed top3, retain complete distance traces.
with gzip.open(R/'union_pair_distances.csv.gz','wt') as f:
 w=csv.writer(f);w.writerow(['system','frame_0','time_ns','pair_id','protein_heavy_min_A','backbone_min_A'])
 for system in ['WT','ANC']:
  for p in pairrows:
   if not p['union_any_frame_lt5_A']:continue
   j,c=p['Ja_index'],p['core_index']
   for fr in range(1000):w.writerow([system,fr,(fr+1)*.01,p['pair_id'],data[system]['dist']['protein_heavy'][fr,j,c],data[system]['dist']['backbone_heavy'][fr,j,c]])
ranked=sorted([x for x in comparisons if x['union_any_frame_lt5_A']],key=lambda x:(-abs(x['delta_full']),x['WT_core'],x['WT_Ja']));table('distance_display_pairs.csv',ranked[:3])
# Save observations, not mechanistic causal claims.
summary=dict(run_scope='HMR WT1/ANC1 only; no separate conventional production found',union_pairs=int(union.sum()),all_pairs=2300,stable_pairs=sum(x['stable_within_observed_run'] for x in comparisons),candidate_eligible_pairs=sum(x['candidate_eligible'] for x in comparisons),DSSP_alpha_fraction={s:{w:float(np.mean(data[s]['sse'][a:b]=='H')) for w,a,b in windows} for s in data},Ja_mean_RMSF_A={s:{w:float(np.mean([x['RMSF_A'] for x in rmsfrows if x['system']==s and x['window']==w and x['region']=='Ja'])) for w,_,_ in windows} for s in data},aligned_Ja_cache_error_A={s:data[s]['cache_error'] for s in data},display_pairs=ranked[:3],new_PDB=0)
(R/'summary.json').write_text(json.dumps(summary,indent=2))
for p,h in checks.items():assert sha(Path(p))==h
(R/'validation.json').write_text(json.dumps(dict(status='PASS',protected_files_unchanged=len(checks),PDB_created=0,all_frames=2000,paired_residues=123,contacts_all_pairs=2300,MDtraj_version=md.__version__,old_aligned_Ja_coordinates_reused_and_matched=summary['aligned_Ja_cache_error_A']),indent=2))
print(json.dumps(summary,indent=2),flush=True)
