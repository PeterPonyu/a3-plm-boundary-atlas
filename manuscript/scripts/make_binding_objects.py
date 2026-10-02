#!/usr/bin/env python3
"""A3 Fig. 5: actual binding object and locked measurements, no mutant 3D model."""
from pathlib import Path
import json,hashlib,os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch
from PIL import Image
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
ROOT=Path(__file__).resolve().parents[1]
OBJ=ROOT/'objects'
DATA=ROOT.parent/'analysis/xfit_calibration_20260918/inputs/scored_3b_wt.csv'

def main():
    # A substituted font changes both printed size and geometry. Refuse it.
    for weight in ('normal', 'bold'):
        face = font_manager.findfont(
            font_manager.FontProperties(family='Arial', weight=weight),
            fallback_to_default=False)
        resolved = font_manager.get_font(face)
        style = resolved.style_name.lower()
        if (resolved.family_name != 'Arial' or
                ('bold' in style) != (weight == 'bold') or
                'italic' in style or 'oblique' in style):
            raise RuntimeError('A3 binding figure requires regular and bold Arial.')
    format_only = os.environ.get('A3_FORMAT_ONLY') == '1'
    figdir = Path(os.environ.get('A3_FIGURE_OUTPUT', str(ROOT/'figures')))
    if not figdir.is_dir():
        raise RuntimeError('A3 binding figure requires an existing output directory.')
    d=pd.read_csv(DATA)
    case=d[d.mutkey.isin(['1A22_A_B|A|I58A','1A22_A_B|A|I58L'])].sort_values('mutkey')
    assert len(case)==2 and set(case.map_method)=={'direct(pos=58)'} and set(case.n_rep)=={1}
    if not format_only:
        case.to_csv(OBJ/'binding_examples.csv',index=False)
    agg=pd.read_csv(ROOT/'figures/plotdata/F4_binding.csv')
    plt.rcParams.update({'font.family':'Arial','font.size':8,'axes.labelsize':9,'xtick.labelsize':8,'ytick.labelsize':8,'pdf.fonttype':42,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(6.5,5.5),facecolor='white')
    for x,cx,y,t,title in [(.025,.265,.985,'A','Deposited binding complex'),(.54,.765,.985,'B','Side-chain substitutions'),(.025,.285,.500,'C','Binding rank agreement'),(.54,.805,.500,'D','Same-site 3B scores')]:
        fig.text(x,y,t,weight='bold',fontsize=11,va='top')
        fig.text(cx,y,title,fontsize=9,ha='center',va='top')
    fig.text(.275,.947,'PDB 1A22 · G120R background',ha='center',fontsize=8)
    a=fig.add_axes([.015,.585,.32,.345]);a.imshow(Image.open(OBJ/'complex.png'));a.set_axis_off()
    a=fig.add_axes([.33,.610,.18,.26]);a.imshow(Image.open(OBJ/'site.png'));a.set_axis_off()
    fig.text(.425,.585,'Ile58 local context',fontsize=8,ha='center')
    fig.text(.255,.550,'Blue: hormone  |  Grey: receptor  |  Orange: Ile58',ha='center',fontsize=8)
    fig.text(.77,.938,'Side-chain connectivity (schematic)',ha='center',fontsize=8)
    for smi,x,y,w,h,lab in [('*[C@@H](C)CC',.63,.797,.26,.12,'Ile58'),('*C',.55,.625,.18,.12,'Ala'),('*CC(C)C',.795,.625,.18,.12,'Leu')]:
        m=Chem.MolFromSmiles(smi);im=Draw.MolToImage(m,size=(600,230));aa=fig.add_axes([x,y,w,h]);aa.imshow(im);aa.set_axis_off();fig.text(x+w/2,y-.008,lab,ha='center',fontsize=8)
    for end in [(.645,.756),(.875,.756)]:fig.add_artist(FancyArrowPatch((.76,.791),end,arrowstyle='->',mutation_scale=10,color='#777777',transform=fig.transFigure,lw=.9))
    for row,x in zip(case.itertuples(),[.645,.885]):
        fig.text(x,.570,f'ΔΔG +{row.ddg:.2f} kcal/mol',ha='center',fontsize=8,color='#A04D1B')
    fig.text(.77,.530,'Positive ΔΔG: weaker binding at 298 K',ha='center',fontsize=8)
    # Original 650M aggregate cells: error bars and denominators retained.
    ax=fig.add_axes([.225,.095,.325,.305]);labels=['All','Non-antibody','Antibody–antigen','Interface','Non-interface']
    for i,r in enumerate(agg.itertuples()):
        c='#0072B2' if i not in [2,4] else '#D55E00'
        ax.errorbar(r.spearman,4-i,xerr=[[r.spearman-r.ci_lo],[r.ci_hi-r.spearman]],fmt='o',ms=4,color=c,capsize=2,lw=1.1)
    ax.set_yticks(range(5));ax.set_yticklabels([f'{s}  (n={int(r.n):,})' for s,r in zip(labels,agg.itertuples())][::-1],fontsize=8)
    ax.set_xlim(-.12,.35);ax.set_ylim(-.6,4.6);ax.set_xticks([-.1,0,.1,.2,.3]);ax.axvline(0,color='#aaaaaa',lw=.7,ls='--');ax.set_xlabel('Spearman(LLR, −ΔΔG)',fontsize=8);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    fig.text(.285,.440,'650M aggregates: 95% bootstrap CI',ha='center',fontsize=8)
    ax=fig.add_axes([.715,.095,.255,.305]);ys=[1,0]
    for r,y in zip(case.itertuples(),ys):
        ax.scatter(r.llr_3b_wt,y,c='#0072B2',s=28);ax.text(r.llr_3b_wt,y+.17,f'{r.llr_3b_wt:.2f}',ha='center',fontsize=8)
    ax.set_ylim(-.5,1.5);ax.set_xlim(-13,0);ax.set_yticks(ys);ax.set_yticklabels(['I58A','I58L']);ax.set_xticks([-12,-8,-4,0]);ax.axvline(0,color='#888888',lw=.7);ax.set_xlabel('ESM-2 3B WT-marginal LLR',fontsize=8);ax.spines['left'].set_visible(False);ax.tick_params(axis='y',length=0)
    fig.text(.81,.440,'Existing illustrative cases',ha='center',fontsize=8)
    for ext in ['pdf','png']:fig.savefig(figdir/f'F7_binding_boundary.{ext}',dpi=450)
    plt.close(fig)
    # The illustrative rows and structural provenance are locked inputs, not
    # products of a typography render. Do not rewrite their sidecars here.
    if not format_only:
        record={'selection':'Existing two directly mapped substitutions at hormone residue 58; same background, temperature and one record each. Illustrative, not an independent validation cohort.','structure_source':'https://files.rcsb.org/download/1A22.cif','structure_identity':'Human growth hormone G120R (chain A), receptor (chain B). Auth residue A58 is ILE; residue A120 is ARG. No hypothetical mutant coordinates drawn.','source_table':str(DATA.relative_to(ROOT.parent)),'source_sha256':hashlib.sha256(DATA.read_bytes()).hexdigest(),'structure_sha256':hashlib.sha256((OBJ/'1a22.cif').read_bytes()).hexdigest(),'rows':case.to_dict('records'),'scope':'Physical binding changes use existing ddg values. 3B sensitivity scores must not be substituted for the 650M aggregate.'}
        (OBJ/'provenance.json').write_text(json.dumps(record,indent=2,allow_nan=True)+'\n')
if __name__=='__main__':main()
