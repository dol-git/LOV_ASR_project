"""Population-only analysis of frozen global cluster labels; no clustering."""
from pathlib import Path
import csv,json,hashlib,os
import numpy as np
from scipy.spatial.distance import jensenshannon
O=Path(__file__).resolve().parent; C=O.parent; B=C.parent
os.environ['MPLCONFIGDIR']=str(O/'matplotlib_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def save(n,rows):
 with (O/n).open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def js(n,x):
 with (O/n).open('x') as f:json.dump(x,f,indent=2,ensure_ascii=False)
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for x in iter(lambda:f.read(8388608),b''):h.update(x)
 return h.hexdigest()
def entropy(p):
 nz=p>0;h=float(-np.sum(p[nz]*np.log(p[nz])));return h,float(np.exp(h))
def metrics(p,q):
 m=(p+q)/2
 def kl(a,b):
  k=a>0;return np.sum(a[k]*np.log2(a[k]/b[k]))
 jsd=float((kl(p,m)+kl(q,m))/2);tv=float(.5*np.abs(p-q).sum())
 assert abs(jsd-jensenshannon(p,q,base=2)**2)<1e-12
 assert -1e-12<=jsd<=1+1e-12 and 0<=tv<=1+1e-12
 return jsd,tv
assert metrics(np.array([1.,0]),np.array([0.,1]))==(1.,1.)
assert metrics(np.array([.5,.5]),np.array([.5,.5]))==(0.,0.)
protected=[p for p in B.rglob('*') if p.is_file() and O not in p.parents]
before={str(p):sha(p) for p in protected};js('protected_sha256_before.json',before)
summary=[];population=[];comparisons=[];toprows=[];allvectors={};valid={};frameblocks=[]
for cut in [1.,1.25]:
 for s in ['WT','ANC']:
  label=np.load(C/f'{s}_{cut:.2f}A_membership.npy',allow_pickle=False);globalrows=read(C/f'{s}_{cut:.2f}A_sweep_population.csv');maps=read(B/'qc'/s/'physical_frame_mapping.csv')
  assert label.shape==(5000,) and np.issubdtype(label.dtype,np.integer)
  k=len(globalrows);assert np.array_equal(np.unique(label),np.arange(1,k+1))
  assert len(maps)==5000 and all(int(r['frame'])==i+1 for i,r in enumerate(maps))
  times=np.array([float(r['time_ns']) for r in maps]);assert np.allclose(times,np.arange(1,5001)*.01,rtol=0,atol=1e-12)
  assert all(int(r['parameter_state'])==1 and float(r['temperature_K'])==298.15 for r in maps)
  count=np.bincount(label,minlength=k+1)[1:];assert np.array_equal(count,[int(r['population_frames']) for r in globalrows]);p=count/5000
  if cut==1.25:
   existing=read(C/'cutoff_1.25A'/f'{s}_cluster_membership_per_frame.csv');assert np.array_equal(label,[int(r['cluster_id']) for r in existing])
  vectors=[p];counts=[count];windows=['full_0_50'];nframes=[5000]
  for b in range(5):
   mask=(times>b*10)&(times<=(b+1)*10);assert mask.sum()==1000
   ct=np.bincount(label[mask],minlength=k+1)[1:];vectors.append(ct/1000);counts.append(ct);windows.append(f'({b*10},{(b+1)*10}]');nframes.append(1000)
   frameblocks.extend(dict(system=s,cutoff_A=cut,frame=int(i+1),time_ns=float(times[i]),block_start_ns=b*10,block_end_ns=(b+1)*10,global_cluster_id=int(label[i])) for i in np.flatnonzero(mask))
  assert np.array_equal(np.sum(counts[1:],axis=0),count) and np.allclose(np.mean(vectors[1:],axis=0),p,atol=1e-15)
  for wi,(win,v,n,ct) in enumerate(zip(windows,vectors,nframes,counts)):
   ranking=np.lexsort((np.arange(1,k+1),-v));h,eff=entropy(v);jsd,tv=metrics(v,p);imax=int(np.argmax(abs(v-p)))
   row=dict(system=s,cutoff_A=cut,window=win,frame_count=n,global_cluster_count=k,occupied_cluster_count=int((ct>0).sum()),largest_cluster_id=int(ranking[0]+1),largest_population_percent=float(v[ranking[0]]*100),block_ranked_top3_percent=float(v[ranking[:3]].sum()*100),block_ranked_top5_percent=float(v[ranking[:5]].sum()*100),fixed_global_top3_percent=float(v[:3].sum()*100),fixed_global_top5_percent=float(v[:5].sum()*100),JSD_vs_full_bits=jsd,TV_vs_full=tv,H_population_nats=h,effective_cluster_number=eff,max_abs_population_difference_pp=float(abs(v-p)[imax]*100),max_difference_cluster_id=imax+1)
   summary.append(row)
   for rank,idx in enumerate(ranking[:5],1):toprows.append(dict(system=s,cutoff_A=cut,window=win,rank=rank,global_cluster_id=int(idx+1),population_frames=int(ct[idx]),population_percent=float(v[idx]*100),full_population_percent=float(p[idx]*100),difference_vs_full_pp=float((v[idx]-p[idx])*100)))
   for idx in range(k):population.append(dict(system=s,cutoff_A=cut,window=win,global_cluster_id=idx+1,population_frames=int(ct[idx]),frame_count=n,population_percent=float(v[idx]*100),full_population_percent=float(p[idx]*100),difference_vs_full_pp=float((v[idx]-p[idx])*100)))
  matrix_js=np.zeros((6,6));matrix_tv=np.zeros((6,6))
  for i in range(6):
   for j in range(i+1,6):
    jsd,tv=metrics(vectors[i],vectors[j]);matrix_js[i,j]=matrix_js[j,i]=jsd;matrix_tv[i,j]=matrix_tv[j,i]=tv
    comparisons.append(dict(system=s,cutoff_A=cut,window_i=windows[i],window_j=windows[j],Jensen_Shannon_divergence_bits=jsd,total_variation_distance=tv))
  allvectors[(cut,s)]=(windows,vectors,matrix_js,matrix_tv)
  valid[f'{s}_{cut:.2f}A']=dict(frames=5000,block_frames=[1000]*5,global_counts_match_existing=True,block_counts_sum_to_full=True,cluster_count=k)
  print(s,cut,'H/eff full',entropy(p),'block JS',[round(metrics(v,p)[0],5) for v in vectors[1:]],flush=True)
save('population_summary.csv',summary);save('all_cluster_block_populations_and_differences.csv',population);save('top5_cluster_populations.csv',toprows);save('population_vector_pairwise_distances.csv',comparisons);save('frame_block_membership.csv',frameblocks)
hetero=[]
for cut in [1.,1.25]:
 for win in ['full_0_50','(0,10]','(10,20]','(20,30]','(30,40]','(40,50]']:
  w=next(r for r in summary if r['cutoff_A']==cut and r['system']=='WT' and r['window']==win);a=next(r for r in summary if r['cutoff_A']==cut and r['system']=='ANC' and r['window']==win)
  hetero.append(dict(cutoff_A=cut,window=win,WT_H_nats=w['H_population_nats'],ANC_H_nats=a['H_population_nats'],WT_minus_ANC_H_nats=w['H_population_nats']-a['H_population_nats'],WT_effective_cluster_number=w['effective_cluster_number'],ANC_effective_cluster_number=a['effective_cluster_number'],WT_minus_ANC_effective_number=w['effective_cluster_number']-a['effective_cluster_number'],WT_occupied_clusters=w['occupied_cluster_count'],ANC_occupied_clusters=a['occupied_cluster_count']))
save('WT_ANC_heterogeneity_comparison.csv',hetero)
plt.rcParams.update({'font.size':9})
for cut in [1.,1.25]:
 fig,axs=plt.subplots(2,2,figsize=(12,10))
 for row,s in enumerate(['WT','ANC']):
  windows,vectors,mj,mt=allvectors[(cut,s)]
  for col,(matrix,title,lim) in enumerate([(mj,'Jensen–Shannon divergence (bits)',1),(mt,'Total variation distance',1)]):
   ax=axs[row,col];im=ax.imshow(matrix,vmin=0,vmax=max(.05,float(matrix.max())),cmap='viridis');ax.set_xticks(range(6),['Full','0–10','10–20','20–30','30–40','40–50'],rotation=45);ax.set_yticks(range(6),['Full','0–10','10–20','20–30','30–40','40–50']);ax.set_title(f'{s}, {cut:.2f} Å: {title}')
   for i in range(6):
    for j in range(6):ax.text(j,i,f'{matrix[i,j]:.3f}',ha='center',va='center',fontsize=8,color='white' if matrix[i,j]<matrix.max()*.5 else 'black')
   fig.colorbar(im,ax=ax,shrink=.8)
 fig.tight_layout();fig.savefig(O/f'{cut:.2f}A_population_distance_matrices.png',dpi=180);plt.close(fig)
 fig,axs=plt.subplots(2,1,figsize=(10,7),sharey=True)
 for ax,s in zip(axs,['WT','ANC']):
  windows,vec,_,_=allvectors[(cut,s)];arr=np.array(vec);vals=np.column_stack([arr[:,:5],1-arr[:,:5].sum(1)]);bottom=np.zeros(6)
  for j in range(6):ax.bar(np.arange(6),vals[:,j]*100,bottom=bottom,label=f'Global cluster {j+1}' if j<5 else 'Other clusters');bottom+=vals[:,j]*100
  ax.set_xticks(range(6),['Full','(0,10]','(10,20]','(20,30]','(30,40]','(40,50]']);ax.set(ylabel='Population (%)',title=f'{s}, cutoff {cut:.2f} Å',ylim=(0,100));ax.legend(ncol=3,fontsize=8,loc='upper center',bbox_to_anchor=(.5,1.01))
 fig.tight_layout();fig.savefig(O/f'{cut:.2f}A_fixed_cluster_populations.png',dpi=180);plt.close(fig)
fig,axs=plt.subplots(1,2,figsize=(11,4),sharey=True)
for ax,cut in zip(axs,[1.,1.25]):
 for s in ['WT','ANC']:
  rr=[r for r in summary if r['system']==s and r['cutoff_A']==cut];ax.plot(range(6),[r['effective_cluster_number'] for r in rr],marker='o',label=s)
 ax.set_xticks(range(6),['Full','0–10','10–20','20–30','30–40','40–50'],rotation=30);ax.set(title=f'Cutoff {cut:.2f} Å',ylabel='Effective cluster number exp(H)');ax.legend()
fig.tight_layout();fig.savefig(O/'effective_cluster_number.png',dpi=180);plt.close(fig)
assert all(sha(Path(p))==h for p,h in before.items());valid['protected_file_count']=len(before);valid['all_previous_grest50_files_SHA256_unchanged']=True;js('validation.json',valid)
js('metadata.json',dict(inputs='Existing *_1.00A_membership.npy and *_1.25A_membership.npy only; global population CSV and physical_frame_mapping.csv used for validation',new_clustering=False,blocks='(0,10],(10,20],(20,30],(30,40],(40,50] ns; 1000 frames each',cluster_support='all labels in full trajectory, including zero-population labels in each block',JSD='0.5 KL_base2(p||m)+0.5 KL_base2(q||m), m=(p+q)/2; divergence, not square-root distance; range 0..1 bits',TV='0.5 sum(abs(p-q)); range 0..1',difference='block percentage minus full-trajectory percentage, percentage points',top_populations='both within-window ranked top3/top5 and fixed global cluster IDs1..3/1..5; no relabeling',entropy='H=-sum(p ln p), 0 ln0=0; nats; exp(H) effective cluster number',interpretation='clustering-based population heterogeneity index only; not thermodynamic/conformational entropy; no comparison with experimental ΔS‡',cross_system_JSD=False,cross_system_reason='WT and ANC labels were defined independently; their population vectors do not share structural categories',full_vs_block_dependence='full distribution includes the block itself; descriptive distances, no independence or significance claim',kinetics=False))
report=['# Block population stability with fixed cluster labels','','Only existing 1.0 and 1.25 Å memberships were used. Clustering and trajectory fitting were not repeated. Each block is (start,end] ns with 1,000 frames; the full dataset has 5,000 frames. The same global cluster labels are used in every block, including labels with zero population.','','Jensen–Shannon divergence (JSD) uses log₂ and is reported in bits, not as a square-root distance. Total variation (TV) is ½Σ|p−q|. Both range from 0 to 1, with zero for identical distributions. Full-versus-block comparisons are descriptive because the full distribution includes the block. WT/ANC labels were defined separately, so cross-system JSD/TV was not calculated.','','## Leading clusters and differences from the full distribution','','Entries show cluster ID: population% (percentage-point difference from the full distribution). IDs remain fixed; only their population ranks are reordered within each row.']
for cut in [1.,1.25]:
 for s in ['WT','ANC']:
  report+=['',f'### {s}, {cut:.2f} Å','','| Block ns | Top five clusters: population% (Δpp) | Ranked top3 % | Ranked top5 % | Fixed ID1–3 % | Fixed ID1–5 % | JSD vs full (bits) | TV vs full |','|---|---|---:|---:|---:|---:|---:|---:|']
  for r in [r for r in summary if r['system']==s and r['cutoff_A']==cut and r['window']!='full_0_50']:
   tr=[v for v in toprows if v['system']==s and v['cutoff_A']==cut and v['window']==r['window']]
   text='; '.join(f"C{v['global_cluster_id']}: {v['population_percent']:.2f} ({v['difference_vs_full_pp']:+.2f})" for v in tr)
   report.append(f"| {r['window']} | {text} | {r['block_ranked_top3_percent']:.2f} | {r['block_ranked_top5_percent']:.2f} | {r['fixed_global_top3_percent']:.2f} | {r['fixed_global_top5_percent']:.2f} | {r['JSD_vs_full_bits']:.6f} | {r['TV_vs_full']:.6f} |")
  pp=[r for r in comparisons if r['system']==s and r['cutoff_A']==cut and r['window_i']!='full_0_50'];mjs=max(pp,key=lambda r:r['Jensen_Shannon_divergence_bits']);mtv=max(pp,key=lambda r:r['total_variation_distance'])
  report+=['',f"Pairwise block JSD range {min(r['Jensen_Shannon_divergence_bits'] for r in pp):.6f}–{mjs['Jensen_Shannon_divergence_bits']:.6f} bits; maximum pair {mjs['window_i']} vs {mjs['window_j']}. TV range {min(r['total_variation_distance'] for r in pp):.6f}–{mtv['total_variation_distance']:.6f}; maximum pair {mtv['window_i']} vs {mtv['window_j']}."]
report+=['','## Population heterogeneity index','','H=−Σp ln p (nats), effective cluster number=exp(H). **This Shannon entropy is a clustering-based heterogeneity index, not thermodynamic or conformational entropy. Do not compare it directly with experimental ΔS‡.**','','| Cutoff Å | Window ns | WT H | ANC H | WT−ANC H | WT exp(H) | ANC exp(H) |','|---|---|---:|---:|---:|---:|---:|']
for r in hetero:report.append(f"| {r['cutoff_A']:.2f} | {r['window']} | {r['WT_H_nats']:.6f} | {r['ANC_H_nats']:.6f} | {r['WT_minus_ANC_H_nats']:+.6f} | {r['WT_effective_cluster_number']:.6f} | {r['ANC_effective_cluster_number']:.6f} |")
report+=['','## Block dependence of WT/ANC differences','']
for cut in [1.,1.25]:
 rr=[r for r in hetero if r['cutoff_A']==cut and r['window']!='full_0_50'];pos=[r['window'] for r in rr if r['WT_minus_ANC_H_nats']>0];neg=[r['window'] for r in rr if r['WT_minus_ANC_H_nats']<0]
 report.append(f"- {cut:.2f} Å: WT H > ANC H in blocks {len(pos)}/5 ({', '.join(pos)}); WT H < ANC H in blocks {len(neg)}/5 ({', '.join(neg)}). WT−ANC H range across blocks {min(r['WT_minus_ANC_H_nats'] for r in rr):+.6f}–{max(r['WT_minus_ANC_H_nats'] for r in rr):+.6f} nats.")
report+=['','These results quantify population differences using fixed labels. Finite-block differences do not establish statistical significance, convergence or a mechanism. Parameter-sorted sample order was not interpreted as kinetics; transition rates and lifetimes were not calculated.','','## Validation and outputs','',f'SHA-256 values for {len(before)} existing grest50 files were unchanged. Validation checked 1,000 frames per block, 5,000 in total, agreement with existing full-dataset cluster counts, and equality of mean block populations with full populations. JSD was checked against squared SciPy Jensen–Shannon distance.','', '- population_summary.csv: full and block summaries, top3/top5, JSD/TV and H/exp(H).','- all_cluster_block_populations_and_differences.csv: populations for all fixed labels and Δpp relative to the full distribution.','- top5_cluster_populations.csv: five leading IDs and populations per window.','- population_vector_pairwise_distances.csv: JSD/TV for full versus block and all ten block pairs.','- WT_ANC_heterogeneity_comparison.csv: H/exp(H) comparisons at matching cutoffs and blocks.','- frame_block_membership.csv: existing labels, block assignments and frame provenance.','- PNG: fixed-label populations, JSD/TV matrices and effective cluster numbers.']
with (O/'REPORT.ko.md').open('x') as f:f.write('\n'.join(report)+'\n')
print('COMPLETE',json.dumps(hetero),flush=True)
