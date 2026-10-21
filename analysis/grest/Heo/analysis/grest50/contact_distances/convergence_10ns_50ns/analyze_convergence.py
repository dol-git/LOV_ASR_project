from pathlib import Path
import csv,json,hashlib,struct,os,itertools
import numpy as np
OUT=Path(__file__).resolve().parent;CONTACT=OUT.parent;BASE=CONTACT.parent;ROOT=BASE.parents[2];PILOT=ROOT/'Heo/analysis/grest10_pilot'
PAIRS=[(4,139,5,140),(14,140,15,141),(96,141,97,142),(107,137,108,138),(90,129,91,130),(90,126,91,127),(74,126,75,127)]
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def save(name,data):
 with (OUT/name).open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def js(name,data):
 with (OUT/name).open('x') as f:json.dump(data,f,indent=2);f.write('\n')
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(8*1024*1024):h.update(b)
 return h.hexdigest()
class DCD:
 def __init__(self,path):
  self.f=path.open('rb');self.endian='<' if struct.unpack('<i',self.f.read(4))[0]==84 else '>';self.f.seek(0)
  h=self.rec();assert h[:4]==b'CORD';self.ctrl=struct.unpack(self.endian+'20i',h[4:]);assert self.ctrl[8]==self.ctrl[11]==0 and self.ctrl[10]
  self.rec();self.n=struct.unpack(self.endian+'i',self.rec())[0];self.start=self.f.tell();self.stride=56+3*(8+4*self.n)
  assert path.stat().st_size==self.start+self.stride*self.ctrl[0]
 def rec(self):
  tag=self.f.read(4);n=struct.unpack(self.endian+'i',tag)[0];b=self.f.read(n);assert len(b)==n and self.f.read(4)==tag;return b
 def frame(self,i):
  self.f.seek(self.start+i*self.stride);cell=np.frombuffer(self.rec(),dtype=self.endian+'f8');assert len(cell)==6 and np.allclose(cell[[1,3,4]],0)
  xyz=np.stack([np.frombuffer(self.rec(),dtype=self.endian+'f4') for _ in range(3)],axis=1).astype(float);assert np.isfinite(xyz).all();return cell[[0,2,5]],xyz

def stat(a):
 q1,q3=np.quantile(a,[.25,.75],method='linear')
 return dict(n=len(a),median_A=float(np.median(a)),mean_A=float(a.mean()),SD_A=float(a.std(ddof=1)),Q1_A=float(q1),Q3_A=float(q3),IQR_A=float(q3-q1))
def ecdf_d(a,b):
 a=np.sort(a);b=np.sort(b);grid=np.union1d(a,b)
 return float(np.max(abs(np.searchsorted(a,grid,side='right')/len(a)-np.searchsorted(b,grid,side='right')/len(b))))
protected=[p for folder in [CONTACT,BASE/'rmsd',BASE/'wt_structure_audit',PILOT] for p in folder.rglob('*') if p.is_file() and OUT not in p.parents]
for sys in ['WT','ANC']:protected.extend([ROOT/'Heo/input'/('ionized_fmc.psf' if sys=='WT' else 'ionized_fmc_A.psf')])
before={str(p):sha(p) for p in protected}
qc={};pilot_data={};mapping_rows=[];selection_rows=[]
for system in ['WT','ANC']:
 pdb=PILOT/'results/remd_convert'/system/'physical_state.pdb';pd=DCD(pdb.with_name('physical_state_1.dcd'));assert pd.ctrl[0]==1000
 lines=[l for l in pdb.read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))];keys=[(l[72:76].strip(),int(l[22:26]),l[17:21].strip(),l[12:16].strip()) for l in lines]
 psf=ROOT/'Heo/input'/('ionized_fmc.psf' if system=='WT' else 'ionized_fmc_A.psf');text=psf.read_text().splitlines();k=next(i for i,l in enumerate(text) if '!NATOM' in l);n=int(text[k].split()[0]);atoms=[l.split() for l in text[k+1:k+1+n]];sel=[a for a in atoms if a[1]=='PROA'];assert keys==[(a[1],int(a[2]),a[3],a[4]) for a in sel]
 raw_idx=np.array([int(a[0])-1 for a in sel]);raw_to_red={int(a[0]):i for i,a in enumerate(sel)}
 k=next(i for i,l in enumerate(text) if '!NBOND' in l);nb=int(text[k].split()[0]);nums=[]
 for l in text[k+1:]:
  if len(nums)>=2*nb:break
  nums.extend(map(int,l.split()))
 bonds=np.array([(raw_to_red[a],raw_to_red[b]) for a,b in zip(nums[::2],nums[1::2]) if a in raw_to_red and b in raw_to_red])
 rawbase=ROOT/'Heo'/f'HMR_{system}'/'grest/prod_10ns';raw=[DCD(rawbase/f'grest_10ns_4r_{system}{r}.dcd') for r in range(1,5)]
 assert all(d.ctrl[0]==1000 and d.ctrl[1]==2500 and d.ctrl[2]==2500 and d.n==n for d in raw)
 rem=[np.loadtxt(rawbase/f'grest_10ns_4r_{system}{r}.rem',dtype=int) for r in range(1,5)]
 assert all(a.shape==(1000,2) and np.array_equal(a[:,0],np.arange(1000)*2500) for a in rem)
 state=np.stack([a[:,1] for a in rem],axis=1);assert np.all(np.sort(state,axis=1)==np.arange(1,5));sources=np.argmax(state==1,axis=1)
 sels=[]
 for j,(a,b,w,x) in enumerate(PAIRS):
  ra,rb=(w,x) if system=='WT' else (a,b);indices=[]
  for residue in [ra,rb]:
   ii=[i for i,key in enumerate(keys) if key[1]==residue and key[3].lstrip('0123456789')[0]!='H'];assert ii;indices.append(ii)
   for i in ii:selection_rows.append(dict(system=system,pair_id=j+1,residue=residue,atom_index=i+1,residue_name=keys[i][2],atom_name=keys[i][3]))
  sels.append(indices)
 # Verify exact homologous selection identity against the established 50-ns analysis.
 existing=read(CONTACT/'heavy_atom_selections.csv')
 for j in range(7):
  old={(int(r['residue']),r['residue_name'],r['atom_name']) for r in existing if r['system']==system and int(r['pair_id'])==j+1}
  new={(r['residue'],r['residue_name'],r['atom_name']) for r in selection_rows if r['system']==system and r['pair_id']==j+1};assert old==new
 values=np.empty((1000,7));error=0.;pbc_difference=0.;bmin=np.inf;bmax=0.;bondmic=0.;qcrows=[]
 for i in range(1000):
  box,px=pd.frame(i);rawbox,full=raw[sources[i]].frame(i);rx=full[raw_idx];assert px.shape==rx.shape and np.allclose(box,rawbox)
  delta=rx[bonds[:,0]]-rx[bonds[:,1]];length=np.linalg.norm(delta,axis=1);mic=np.linalg.norm(delta-rawbox*np.rint(delta/rawbox),axis=1)
  bmin=min(bmin,float(length.min()));bmax=max(bmax,float(length.max()));bd=float(abs(length-mic).max());bondmic=max(bondmic,bd)
  localerror=0.;localpbc=0.
  for j,(left,right) in enumerate(sels):
   d=rx[left,None,:]-rx[None,right,:];direct=float(np.linalg.norm(d,axis=2).min());minimum=float(np.linalg.norm(d-rawbox*np.rint(d/rawbox),axis=2).min())
   p=float(np.linalg.norm(px[left,None,:]-px[None,right,:],axis=2).min());localerror=max(localerror,abs(p-direct));localpbc=max(localpbc,abs(direct-minimum));values[i,j]=p
  error=max(error,localerror);pbc_difference=max(pbc_difference,localpbc)
  mapping_rows.append(dict(system=system,frame=i+1,step=(i+1)*2500,time_ns=(i+1)*.01,source_replica=int(sources[i]+1),parameter_state=1,rem_step=i*2500))
  qcrows.append(dict(frame=i+1,max_pilot_vs_raw_distance_difference_A=localerror,max_raw_direct_MIC_difference_A=localpbc,bond_min_A=float(length.min()),bond_max_A=float(length.max()),max_bond_direct_MIC_difference_A=bd))
 assert error<1e-4 and pbc_difference<1e-6 and bondmic<1e-6 and bmin>.5 and bmax<2.5
 save(f'{system}_pilot_pbc_and_raw_validation.csv',qcrows);pilot_data[system]=values
 qc[system]=dict(frames=1000,frame_interval_ps=10,existing_pilot_rigid_body_fit=True,additional_fitting=False,max_pilot_vs_source_raw_distance_difference_A=error,max_raw_direct_MIC_difference_A=pbc_difference,bond_min_A=bmin,bond_max_A=bmax,max_bond_direct_MIC_difference_A=bondmic,selection_identical_to_50ns=True)
 pd.f.close()
 for d in raw:d.f.close()
 print('PILOT_QC',system,json.dumps(qc[system]),flush=True)
js('pilot_validation.json',qc);save('pilot_frame_mapping.csv',mapping_rows);save('pilot_heavy_atom_selections.csv',selection_rows)
allstats=[];changes=[];ecdfrows=[];ecdfcomparisons=[];summaries=[];occ=[];datasets={}
for system in ['WT','ANC']:
 existing=read(CONTACT/f'{system}_minimum_distances.csv');full=np.empty((5000,7));seen=set()
 for r in existing:
  i,j=int(r['frame'])-1,int(r['pair_id'])-1;assert (i,j) not in seen;seen.add((i,j));full[i,j]=float(r['minimum_heavy_atom_distance_A']);assert abs(float(r['time_ns'])-(i+1)*.01)<1e-12
 assert len(seen)==35000
 pdata=pilot_data[system];dataset={'pilot_0_10':pdata,'full_0_50':full,**{f'block_{b*10}_{(b+1)*10}':full[b*1000:(b+1)*1000] for b in range(5)}};datasets[system]=dataset
 save(f'{system}_pilot_minimum_distances.csv',[dict(pair_id=j+1,frame=i+1,time_ns=(i+1)*.01,minimum_heavy_atom_distance_A=pdata[i,j]) for j in range(7) for i in range(1000)])
 for j,(a,b,w,x) in enumerate(PAIRS):
  ra,rb=(w,x) if system=='WT' else (a,b);common=dict(system=system,pair_id=j+1,residue_a=ra,residue_b=rb);local={}
  for name,data in dataset.items():
   arr=data[:,j];local[name]=stat(arr);allstats.append(dict(**common,dataset=name,**local[name]));occ.append(dict(**common,dataset=name,n=len(arr),cutoff_A=4.5,occupancy_percent=float(100*np.mean(arr<=4.5))))
   for rank,value in enumerate(np.sort(arr),1):ecdfrows.append(dict(**common,dataset=name,distance_A=value,empirical_CDF=rank/len(arr)))
  names=[f'block_{b*10}_{(b+1)*10}' for b in range(5)];med=[local[n]['median_A'] for n in names];iqr=[local[n]['IQR_A'] for n in names]
  for b,name in enumerate(names):changes.append(dict(**common,block=name,median_A=med[b],IQR_A=iqr[b],median_change_from_previous_A=med[b]-med[b-1] if b else '',IQR_change_from_previous_A=iqr[b]-iqr[b-1] if b else '',median_change_from_first_A=med[b]-med[0],IQR_change_from_first_A=iqr[b]-iqr[0]))
  comparisons=[('pilot_0_10','full_0_50')]+list(itertools.combinations(names,2))+[('pilot_0_10',n) for n in names]+[(n,'full_0_50') for n in names]
  for left,right in comparisons:ecdfcomparisons.append(dict(**common,dataset_a=left,dataset_b=right,ECDF_max_absolute_difference=ecdf_d(dataset[left][:,j],dataset[right][:,j])))
  ds=[ecdf_d(dataset[left][:,j],dataset[right][:,j]) for left,right in itertools.combinations(names,2)]
  summaries.append(dict(**common,pilot_median_A=local['pilot_0_10']['median_A'],full50_median_A=local['full_0_50']['median_A'],full50_minus_pilot_median_A=local['full_0_50']['median_A']-local['pilot_0_10']['median_A'],block_median_min_A=min(med),block_median_max_A=max(med),block_median_span_A=max(med)-min(med),block_IQR_min_A=min(iqr),block_IQR_max_A=max(iqr),block_IQR_span_A=max(iqr)-min(iqr),max_pairwise_block_ECDF_difference=max(ds),pilot_vs_full50_ECDF_difference=ecdf_d(pdata[:,j],full[:,j]),pilot_vs_50_firstblock_ECDF_difference=ecdf_d(pdata[:,j],full[:1000,j]),pilot_vs_50_firstblock_max_frame_distance_difference_A=float(abs(pdata[:,j]-full[:1000,j]).max())))
save('distribution_statistics.csv',allstats);save('block_changes.csv',changes);save('ECDF_values.csv',ecdfrows);save('ECDF_comparisons.csv',ecdfcomparisons);save('convergence_summary.csv',summaries);save('auxiliary_4_5_occupancy.csv',occ)
os.environ['MPLCONFIGDIR']=str(OUT/'matplotlib_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
(OUT/'figures').mkdir(exist_ok=True)
colors=['#2166ac','#67a9cf','#8c510a','#dfc27d','#b2182b']
for j,(a,b,w,x) in enumerate(PAIRS):
 fig,axes=plt.subplots(2,2,figsize=(12,8))
 for col,system in enumerate(['WT','ANC']):
  data=datasets[system]
  for name,color,lw in [('pilot_0_10','black',1.5),('full_0_50','#88419d',1.5)]+[(f'block_{i*10}_{(i+1)*10}',colors[i],.9) for i in range(5)]:
   arr=np.sort(data[name][:,j]);axes[0,col].step(arr,np.arange(1,len(arr)+1)/len(arr),where='post',label=name,color=color,lw=lw)
  axes[0,col].set_title(system);axes[0,col].set_xlabel('Minimum heavy-atom distance (A)');axes[0,col].set_ylabel('ECDF');axes[0,col].legend(fontsize=7)
  med=[];lo=[];hi=[]
  for i in range(5):
   st=stat(data[f'block_{i*10}_{(i+1)*10}'][:,j]);med.append(st['median_A']);lo.append(st['median_A']-st['Q1_A']);hi.append(st['Q3_A']-st['median_A'])
  axes[1,col].errorbar(np.arange(1,6),med,yerr=[lo,hi],fmt='o-',capsize=4,label='block median [Q1,Q3]')
  axes[1,col].axhline(np.median(data['pilot_0_10'][:,j]),color='black',ls='--',label='pilot median');axes[1,col].axhline(np.median(data['full_0_50'][:,j]),color='#88419d',ls=':',label='50-ns median');axes[1,col].set_xticks(np.arange(1,6),['0-10','10-20','20-30','30-40','40-50']);axes[1,col].set_xlabel('50-ns trajectory block (ns)');axes[1,col].set_ylabel('Distance (A)');axes[1,col].legend(fontsize=7)
 fig.suptitle(f'ANC {a}-{b} / WT {w}-{x}');fig.tight_layout();fig.savefig(OUT/'figures'/f'pair_{j+1}_ECDF_block_comparison.png',dpi=180);plt.close(fig)
assert all(sha(Path(p))==h for p,h in before.items());js('protected_sha256.json',{'unchanged':True,'files':before})
js('metadata.json',dict(pairs=PAIRS,pilot_frames=1000,full50_frames=5000,block_frames=1000,primary='continuous distance distributions',existing_50ns_distances_reused=True,block_convention='(start,end] ns',SD_ddof=1,quantiles='linear',IQR='Q3-Q1',change_sign='later minus earlier; full50 minus pilot',ECDF_difference='sup_x |F_a(x)-F_b(x)|, descriptive only, no p-values or independence assumption',new_alignment=False,additional_RMSD=False,convergence_claim='finite-window stability assessment only; no proof of convergence',pilot_source_mapping='REM state before exchange at step (frame-1)*2500 identifies source for saved frame at frame*2500; verified against all 1000 converted frames through pair distances'))
print('SUMMARIES',json.dumps(summaries),flush=True)
