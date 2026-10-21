from pathlib import Path
import csv,json,gzip
import numpy as np
import MDAnalysis as m
from MDAnalysis.lib.distances import distance_array
R=Path(__file__).resolve().parent
read=lambda p:list(csv.DictReader(open(p)))
mp=read(R.parent/'Ja_core_visual_20261019/residue_correspondence.csv');src=read(R.parent/'Ja_core_visual_20261019/selected_frames.csv');ids=read(R.parent/'Heo/analysis/reanalysis_20261006/00_provenance/atom_identity/atom_identity.csv')
core=[x for x in mp if x['region']=='core'];ja=[x for x in mp if x['region']=='Ja'];UU={};groups={};rows=[]
for system in ['WT','ANC']:
 s=next(x for x in src if x['system']==system);assert 'grest' not in s['source_dcd'].lower();UU[system]=m.Universe(s['source_topology'],s['source_dcd'])
for role,maps in [('core',core),('Ja',ja)]:
 groups[role]={s:[] for s in UU}
 for row in maps:
  di={}
  for system,old in [('WT','Extant'),('ANC','Ancestral')]:
   r=int(row[system+'_resid']);aa={x['atom_name']:x for x in ids if x['system']==old and x['segment']=='PROA' and int(x['residue'])==r and x['is_hydrogen']=='False' and x['portion'] in ['protein','protein_cysteine']};di[system]=aa
  common=sorted(n for n in set(di['WT'])&set(di['ANC']) if (di['WT'][n]['element'],di['WT'][n]['atom_type'])==(di['ANC'][n]['element'],di['ANC'][n]['atom_type']))
  assert len(common)>=4
  for system in UU:
   group=np.array([int(di[system][n]['index'])-1 for n in common]);groups[role][system].append(group)
   for n,index in zip(common,group):rows.append(dict(system=system,role=role,resid=int(row[system+'_resid']),WT_resid=int(row['WT_resid']),ANC_resid=int(row['ANC_resid']),name=n,CHARMM_type=di[system][n]['atom_type'],psf_index_0=int(index)))
with (R/'matched_atom_control_definition.json').open('w') as f:json.dump(dict(rule='intersection by atom name AND element AND CHARMM type, separately at each direct MD corresponding residue; no atom-mass criterion',interpretation='geometry control retaining common chemistry/types; mutation-specific distal functional groups are removed, so failure to persist does not prove contact is fictitious or side-chain energetic effects absent',thresholds_A=[4,4.5,5],windows=['full','first_half','second_half'],frame_exclusions=0),f,indent=2)
with (R/'matched_heavy_atom_selections.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
data={}
for system,u in UU.items():
 cg=groups['core'][system];jg=groups['Ja'][system];ci=np.concatenate(cg);ji=np.concatenate(jg);co=np.r_[0,np.cumsum([len(x) for x in cg])[:-1]];jo=np.r_[0,np.cumsum([len(x) for x in jg])[:-1]];a=np.empty((1000,23,100))
 for ts in u.trajectory:a[ts.frame]=np.minimum.reduceat(np.minimum.reduceat(distance_array(u.atoms[ji].positions,u.atoms[ci].positions),co,axis=1),jo,axis=0)
 data[system]=a
np.savez_compressed(R/'matched_heavy_distances.npz',**data)
pairs=read(R/'contact_robustness.csv');out=[];new=[]
for p in pairs:
 j,c=int(p['Ja_index']),int(p['core_index']);delta={}
 for cutoff in [4.,4.5,5.]:
  for window,start,end in [('full',0,1000),('first_half',0,500),('second_half',500,1000)]:
   aw=float(np.mean(data['WT'][start:end,j,c]<cutoff));ab=float(np.mean(data['ANC'][start:end,j,c]<cutoff));delta[cutoff,window]=aw-ab
   if p['union_any_frame_lt5_A']=='True':out.append(dict(pair_id=p['pair_id'],selection='matched_heavy_name_element_type',threshold_A=cutoff,window=window,WT_contact_fraction=aw,ANC_contact_fraction=ab,WT_minus_ANC_fraction=aw-ab,WT_median_distance_A=float(np.median(data['WT'][start:end,j,c])),ANC_median_distance_A=float(np.median(data['ANC'][start:end,j,c]))))
 sign=np.sign(float(p['delta_full']));p['matched_atom_delta_full']=delta[4.5,'full'];p['matched_atom_delta_first']=delta[4.5,'first_half'];p['matched_atom_delta_second']=delta[4.5,'second_half'];p['matched_atom_delta_4A']=delta[4.,'full'];p['matched_atom_delta_5A']=delta[5.,'full'];p['matched_atom_same_sign_all_thresholds']=all(np.sign(delta[t,'full'])==sign and abs(delta[t,'full'])>=.1 for t in [4.,4.5,5.]);p['matched_atom_same_sign_both_halves']=all(np.sign(delta[4.5,w])==sign and abs(delta[4.5,w])>=.1 for w in ['first_half','second_half']);new.append(p)
for fn,rr in [('matched_atom_contact_fractions.csv',out),('contact_robustness_with_atom_control.csv',new)]:
 with (R/fn).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
ss=[x for x in new if x['stable_within_observed_run']=='True'];print('Primary stable pairs:',len(ss));print('Common-atom control same direction all thresholds and both halves:',sum(x['matched_atom_same_sign_all_thresholds'] and x['matched_atom_same_sign_both_halves'] for x in ss))
for p in ss:print(p['pair_id'],p['delta_full'],p['matched_atom_delta_full'],p['matched_atom_same_sign_all_thresholds'],p['matched_atom_same_sign_both_halves'])
