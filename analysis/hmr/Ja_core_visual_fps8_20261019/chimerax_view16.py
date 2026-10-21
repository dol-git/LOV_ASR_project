"""ChimeraX 1.11.1 viewer. Open via File > Open; commands: lovview WT full / lovview Anc ja.
Reads only the adjacent manifest/config/final_pdb. No coordinate edits, alignment, or files saved.
"""
from pathlib import Path,PurePosixPath
import csv,json,hashlib
import numpy as np
from chimerax.core.commands import run,register,CmdDesc,EnumOf
from chimerax.core.models import Model
from chimerax.geometry import Place
from chimerax.graphics import OrthographicCamera
from chimerax.std_commands.cartoon import cartoon
from chimerax.label.label2d import label_create,label_delete

ROOT=Path(__file__).resolve().parent
CFG=json.loads((ROOT/'chimerax_display_config.json').read_text())
with (ROOT/'selected_frames.csv').open() as handle:
    rows=list(csv.DictReader(handle))
expected=CFG['frames_per_system']
for system in ['WT','ANC']:
    rr=sorted((x for x in rows if x['system']==system),key=lambda x:int(x['rank']))
    if len(rr)!=expected or [int(x['rank']) for x in rr]!=list(range(1,expected+1)):
        raise ValueError('Manifest does not contain the expected ordered frames')
    if sum(x['role']=='central_medoid' for x in rr)!=1 or rr[0]['role']!='central_medoid':
        raise ValueError('Expected rank 1 to be the sole existing global medoid')
old=getattr(session,'lov_visualization',None)
if old:
    for lab in old['labels']:
        if lab is not None:label_delete(session,[lab])
    opened=session.models.list()
    session.models.close([g for g in old['groups'] if g in opened])
state=dict(groups=[],systems={},labels=[],mode='full',config=CFG,models=[])
for system in ['WT','ANC']:
    group=Model('LOV '+system,session)
    rr=sorted((x for x in rows if x['system']==system),key=lambda x:int(x['rank']))
    info=CFG['systems'][system]
    for row in rr:
        filename=PurePosixPath(row['pdb']).name
        path=ROOT/'final_pdb'/filename
        if 'grest' in str(path).lower():raise ValueError('Excluded path')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=CFG['PDB_sha256'][filename]:
            raise ValueError('PDB content differs from verified configuration: '+filename)
        models,status=session.open_command.open_data(str(path))
        if len(models)!=1:raise ValueError('Expected one PDB model: '+filename)
        mol=models[0]
        coords_before=mol.atoms.coords.copy()
        bonds_before=len(mol.bonds)
        position_before=mol.scene_position.matrix.copy()
        res=mol.residues
        chain=np.asarray(res.chain_ids)==info['chain']
        ja=chain & (res.numbers>=info['Ja'][0]) & (res.numbers<=info['Ja'][1])
        core=chain & (res.numbers>=info['core'][0]) & (res.numbers<=info['core'][1])
        if int(ja.sum())!=23 or set(res[ja].numbers)!=set(range(info['Ja'][0],info['Ja'][1]+1)):
            raise ValueError('Actual PDB chain/residue selection mismatch: '+filename)
        group.add([mol])
        session.models.add([group]) if group.id is None else None
        # Native cartoon smoothing=0: original atom positions remain the control points.
        cartoon(session,atoms=mol.atoms,smooth=0,suppress_backbone_display=True,spine=False)
        mol.atoms.displays=False
        mol.bonds.displays=False
        res.ribbon_displays=False
        res.ribbon_colors=tuple(CFG['core_RGBA'])
        is_medoid=row['role']=='central_medoid'
        if is_medoid:res[core].ribbon_displays=True
        res[ja].ribbon_displays=True
        alpha=255 if is_medoid else int(2.56*(100-CFG['additional_Ja_transparency_percent']))
        res[ja].ribbon_colors=tuple(info['Ja_RGB'])+(alpha,)
        # Never create/link pseudobonds; hide all reader-provided pseudobond drawings.
        mol.ribbon_tether_scale=0
        for pbg in mol.pbg_map.values():pbg.display=False
        if not np.array_equal(coords_before,mol.atoms.coords) or len(mol.bonds)!=bonds_before:
            raise RuntimeError('Display operation changed atomic coordinates or bonds')
        if not np.array_equal(position_before,mol.scene_position.matrix):
            raise RuntimeError('Unexpected model transform')
        state['models'].append(mol)
    group.display=False
    state['groups'].append(group);state['systems'][system]=group
session.lov_visualization=state

def show_view(session,system,mode=None):
    state=session.lov_visualization;cfg=state['config']
    key=system.upper()
    if key not in state['systems']:raise ValueError('Use WT or Anc')
    if mode is None:mode=state['mode']
    if mode not in ['full','ja']:raise ValueError('Use full or ja')
    state['mode']=mode
    for s,g in state['systems'].items():g.display=(s==key)
    for lab in state['labels']:
        if lab is not None:label_delete(session,[lab])
    state['labels'].clear()
    if session.ui.is_gui:run(session,'windowsize %d %d'%tuple(cfg['viewport_pixels']))
    preset=cfg['camera_presets'][mode]
    view=session.main_view
    view.camera=OrthographicCamera(field_width=preset['field_width_A'])
    view.camera.position=Place(axes=cfg['camera_axes'],origin=preset['camera_origin'])
    view.center_of_rotation=np.array(preset['center'])
    view.background_color=(1,1,1,1)
    view.clip_planes.clear()
    if session.ui.is_gui:
        title='%s | medoid + %d farthest-point frames | %s'%(key,cfg['frames_per_system']-1,'Full structure' if mode=='full' else 'Jalpha close-up')
        footer='Frames selected to emphasize structural diversity; do not reflect observation frequencies or equilibrium populations.'
        for name,text,y,size in [('lov_title',title,.94,22),('lov_footer',footer,.035,12)]:
            state['labels'].append(label_create(session,name,text=text,color=(30,30,30,255),xpos=.03,ypos=y,size=size))
    session.logger.info('LOV %s, %s: fixed shared camera center %s, orthographic width %.4f A; no fitting.'%(key,mode,preset['center'],preset['field_width_A']))

if not old:
    register('lovview',CmdDesc(required=[('system',EnumOf(('WT','Anc','ANC')))],optional=[('mode',EnumOf(('full','ja')))],synopsis='Switch WT/Anc with identical full or Jalpha camera'),show_view,logger=session.logger)
show_view(session,'WT','full')
session.logger.info('Commands: lovview WT full; lovview Anc full; lovview WT ja; lovview Anc ja. Without a mode, retain the current preset.')
