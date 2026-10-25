#!/usr/bin/env python3
"""Broad, adaptive Gaussian amplitude survey; all offsets are raw ratios."""
import argparse,csv,hashlib,json,math,os,re,subprocess,threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR','/tmp/tdk-window-mpl')
import numpy as np
HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
DATA=ROOT/'data';OUT=ROOT/'results';CACHE_DIR=ROOT/'.cache'
for path in [DATA,OUT,CACHE_DIR]:path.mkdir(exist_ok=True)
REFERENCE=json.loads((ROOT/'critical_amplitude/reference.json').read_text())
AC=REFERENCE['critical_amplitude'];AC_ERROR=REFERENCE['reported_uncertainty']
SOURCE=hashlib.sha256((HERE/'solver.cpp').read_bytes()).hexdigest()
JOURNAL=DATA/'runs.jsonl';CACHE={};LOCK=threading.Lock()
if JOURNAL.exists():
    for line in JOURNAL.read_text().splitlines():
        try:
            row=json.loads(line)
            if row['source']==SOURCE:CACHE[row['key']]=row
        except (ValueError,KeyError):pass

def write_csv(path,rows):
    if not rows:return
    with path.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

def build():
    executable=CACHE_DIR/'solver';stamp=CACHE_DIR/'source.sha256'
    if executable.exists() and stamp.exists() and stamp.read_text().strip()==SOURCE:return
    subprocess.run(['g++','-O3','-std=c++17','-DSPATIAL_ORDER=4',str(HERE/'solver.cpp'),'-o',str(executable)],check=True)
    stamp.write_text(SOURCE+'\n')

def simulate(amplitude,n=999,cfl=.4,horizon=100,radius=25,start=20,trajectory=None):
    settings=dict(amplitude=float(amplitude),n=int(n),cfl=float(cfl),horizon=float(horizon),radius=float(radius),start=float(start),level=48)
    key=json.dumps(settings,sort_keys=True)
    if trajectory is None and key in CACHE:return CACHE[key]
    args=[str(CACHE_DIR/'solver'),repr(float(amplitude)),str(n),str(cfl),str(horizon),'48',str(radius),str(start),'7','.02',str(trajectory) if trajectory else '-']
    p=subprocess.run(args,check=True,text=True,capture_output=True)
    event,time,steps,sign,peak,final,actual,t12,t24,zeros=p.stdout.split()
    row=dict(key=key,source=SOURCE,**settings,relative_offset=(amplitude-AC)/AC,
             collapsed=int(event),end_time=float(time),steps=int(steps),event_sign=int(sign),
             zero_crossings=int(zeros),t12=float(t12),t24=float(t24),peak_center=float(peak),
             final_center=float(final),actual_radius=float(actual))
    if row['collapsed'] and row['event_sign']!=(-1 if row['zero_crossings']%2 else 1):
        raise RuntimeError('Collapse sign inconsistent with zero crossings')
    with LOCK:
        CACHE[key]=row
        with JOURNAL.open('a') as f:f.write(json.dumps(row)+'\n')
    return row

def state(row):
    return f'{row["event_sign"]:+d}/{row["zero_crossings"]}' if row['collapsed'] else 'not_crossed'

def primary(positive=False):
    rows=[r for r in CACHE.values() if r['n']==999 and r['cfl']==.4 and r['horizon']==100 and r['radius']==25 and r['start']==20]
    if positive:rows=[r for r in rows if r['relative_offset']>=1e-8*(1-1e-7)]
    return sorted(rows,key=lambda r:r['amplitude'])

def batch(amplitudes,tag,workers=4,**settings):
    amplitudes=sorted(set(float(a) for a in amplitudes))
    def work(a):return simulate(a,**settings)
    rows=[]
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for row in pool.map(work,amplitudes):
            rows.append(row)
            if len(rows)%100==0:print(f'{tag}: {len(rows)}/{len(amplitudes)}',flush=True)
    print(f'{tag}: completed {len(rows)} configurations',flush=True)
    return rows

def primary_signature():
    payload='\n'.join(f'{r["amplitude"]:.17g}:{state(r)}' for r in primary())
    return hashlib.sha256(payload.encode()).hexdigest()

def groups(rows):
    output=[]
    for row in rows:
        if not output or state(output[-1][-1])!=state(row):output.append([])
        output[-1].append(row)
    return output

def survey():
    positive=np.unique(np.r_[np.geomspace(1e-8,3,2001),np.linspace(.01,3,601)])
    negative=-np.geomspace(1e-8,.1,121)
    amplitudes=np.r_[np.linspace(.2,.9*AC,81),AC*(1+negative),AC*(1+positive)]
    batch(amplitudes,'broad survey')

def refine():
    # Probe the interior of each initially detected non-crossing gap, as well as edges.
    for round_index in range(2):
        rows=primary(True);gs=groups(rows);extra=[]
        for i,g in enumerate(gs):
            if not g[0]['collapsed'] and i>0 and i<len(gs)-1:
                lo=gs[i-1][-1]['relative_offset'];hi=gs[i+1][0]['relative_offset']
                extra.extend(AC*(1+np.geomspace(lo,hi,97 if round_index==0 else 49)))
        if extra:batch(extra,f'gap interior pass {round_index+1}')
    refine_edges()

def refine_edges():
    # Classify quarter points around every change; intermediate branches are retained.
    for round_index in range(8):
        rows=primary(True);extra=[]
        for a,b in zip(rows,rows[1:]):
            if state(a)!=state(b) and math.log10(b['relative_offset']/a['relative_offset'])>1e-6:
                extra.extend(AC*(1+np.geomspace(a['relative_offset'],b['relative_offset'],5)[1:-1]))
        if not extra:break
        batch(extra,f'adaptive edge pass {round_index+1}')
    gs=groups(primary(True))
    print('Resolved primary groups:',[(state(g[0]),len(g)) for g in gs],flush=True)

def representative(g):
    mid=math.sqrt(g[0]['relative_offset']*g[-1]['relative_offset'])
    return min(g,key=lambda r:abs(math.log(r['relative_offset']/mid)))

def validate():
    # Newly discovered narrow windows also receive interior samples before validation.
    extra=[]
    for g in groups(primary(True)):
        if not g[0]['collapsed'] and len(g)<65 and g[-1]['relative_offset']>g[0]['relative_offset']:
            extra.extend(AC*(1+np.geomspace(g[0]['relative_offset'],g[-1]['relative_offset'],65)))
    if extra:batch(extra,'new-window interior refinement')
    refine_edges()
    reps=[representative(g) for g in groups(primary(True))]
    def work(item):
        i,base=item
        middle=simulate(base['amplitude'],n=1997)
        fine=simulate(base['amplitude'],n=3993) if base['collapsed'] else None
        large=simulate(base['amplitude'],n=1997,horizon=200,radius=75,start=60) if not base['collapsed'] else None
        return dict(group=i+1,amplitude=base['amplitude'],relative_offset=base['relative_offset'],
                    primary_state=state(base),n1997_state=state(middle),n3993_state=state(fine) if fine else '',
                    extended_state=state(large) if large else '',
                    n999_t12=base['t12'],n1997_t12=middle['t12'],n3993_t12=fine['t12'] if fine else '',
                    interior_confirmed=int(state(base)==state(middle) and (not fine or state(base)==state(fine)) and (not large or state(large)=='not_crossed')))
    with ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(work,enumerate(reps)))
    write_csv(OUT/'resolution_checks.csv',rows)
    # Extra near-critical extension is kept separate from the primary window count.
    batch(AC*(1+np.geomspace(1e-10,1e-8,81)),'near-critical extension',n=1997)
    trajectories={}
    for i,base in enumerate(reps):
        path=CACHE_DIR/f'trajectory_{i}.csv'
        simulate(base['amplitude'],n=1997,trajectory=path)
        array=np.genfromtxt(path,delimiter=',',names=True)
        trajectories[f'group_{i+1:02d}']=np.column_stack([array[n] for n in array.dtype.names])
        path.unlink()
    np.savez_compressed(DATA/'representative_trajectories.npz',**trajectories)
    (DATA/'validation_manifest.json').write_text(json.dumps(dict(primary_signature=primary_signature(),source=SOURCE),indent=2)+'\n')
    print('Resolution checks:',rows,flush=True)

def analyze():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import MaxNLocator
    from matplotlib.colors import Normalize
    plt.style.use(HERE/'tdk.mplstyle')
    rows=primary();positive=primary(True);gs=groups(positive)
    checks={}
    validation_current=False
    if (DATA/'validation_manifest.json').exists():
        validation_current=json.loads((DATA/'validation_manifest.json').read_text()).get('primary_signature')==primary_signature()
    if validation_current and (OUT/'resolution_checks.csv').exists():
        with (OUT/'resolution_checks.csv').open() as f:checks={int(x['group']):x for x in csv.DictReader(f)}
    clean=[{k:v for k,v in r.items() if k not in ['key','source']} for r in rows]
    write_csv(OUT/'amplitude_samples.csv',clean)
    windows=[];wid=0;gid=0
    for i,g in enumerate(gs):
        if g[0]['collapsed']:wid+=1
        else:gid+=1
        before=gs[i-1][-1] if i else None;after=gs[i+1][0] if i+1<len(gs) else None
        rep=representative(g);check=checks.get(i+1,{})
        windows.append(dict(group=i+1,window=f'B{wid:02d}' if g[0]['collapsed'] else f'G{gid:02d}',
                            kind='collapse_band' if g[0]['collapsed'] else 'noncrossing_window',
                            state=state(g[0]),event_sign=g[0]['event_sign'],zero_crossings=g[0]['zero_crossings'] if g[0]['collapsed'] else '',
                            sample_count=len(g),delta_first=g[0]['relative_offset'],delta_last=g[-1]['relative_offset'],
                            left_outer_delta=before['relative_offset'] if before else '',right_outer_delta=after['relative_offset'] if after else '',
                            amplitude_first=g[0]['amplitude'],amplitude_last=g[-1]['amplitude'],
                            representative_amplitude=rep['amplitude'],representative_delta=rep['relative_offset'],
                            resolution_confirmed=check.get('interior_confirmed','not_checked')))
    write_csv(OUT/'windows.csv',windows)
    collapse_windows=[w for w in windows if w['kind']=='collapse_band']
    gaps=[w for w in windows if w['kind']=='noncrossing_window']
    fig,axes=plt.subplots(2,1,figsize=(12,7),layout='constrained',sharex=True)
    cmap=plt.get_cmap('tab10')
    for i,(g,w) in enumerate(zip(gs,windows)):
        delta=[r['relative_offset'] for r in g]
        if g[0]['collapsed']:
            color=cmap(int(w['window'][1:])-1);time=[r['t12'] for r in g]
            axes[0].scatter(delta,time,s=9,color=color,marker='o' if g[0]['event_sign']>0 else 'v',linewidths=0)
            axes[1].scatter(delta,[r['zero_crossings'] for r in g],s=8,color=color,linewidths=0)
            mid=math.sqrt(delta[0]*delta[-1]);axes[0].annotate(w['window'],(mid,np.median(time)),xytext=(0,-16 if w['window']=='B01' else 9),textcoords='offset points',ha='center',fontsize=8)
        else:
            axes[0].scatter(delta,[.95]*len(delta),transform=axes[0].get_xaxis_transform(),marker='x',s=14,color='#8064a2')
            for ax in axes:ax.axvspan(delta[0],delta[-1],color='#8064a2',alpha=.10)
    axes[0].legend(handles=[Line2D([],[],marker='o',ls='',color='#405169',label='Positive-field collapse'),Line2D([],[],marker='v',ls='',color='#405169',label='Negative-field collapse')],loc='upper right',fontsize=8)
    axes[0].text(.015,.985,r'$\times$: 48 not reached by $t=100$ (upper band)',transform=axes[0].transAxes,va='top',fontsize=8)
    axes[0].set(ylabel=r'Time to $|\phi(0,t)|=12$, $t_{12}$',title='(a) Cutoff time: collapse bands and non-crossing windows')
    axes[1].set(xscale='log',xlabel=r'Raw amplitude excess, $\delta=(A-A_c)/A_c$',ylabel='Central-field zero crossings',title='(b) Oscillation count; shading marks non-crossing windows')
    axes[1].set_yticks(sorted({r['zero_crossings'] for r in positive if r['collapsed']}))
    for ax in axes:ax.grid(which='minor',alpha=.12)
    fig.suptitle(f'Broad Gaussian amplitude search: {gid} non-crossing windows, {wid} collapse bands\nFourth-order solver | N=999 | classification to field 48 or t=100',fontsize=12,fontweight='bold')
    fig.savefig(OUT/'window_map.png',dpi=300);plt.close(fig)
    fig,ax=plt.subplots(figsize=(11,4.5),layout='constrained')
    for sign,marker,color,label in [(1,'o','#1f77b4','Positive-field collapse'),(-1,'v','#e2642a','Negative-field collapse'),(0,'x','#8064a2','48 not reached by t=100')]:
        selected=[r for r in rows if r['event_sign']==sign]
        ax.scatter([r['amplitude'] for r in selected],[r['t12'] if r['collapsed'] else 100 for r in selected],s=10,marker=marker,color=color,label=label)
    ax.axvline(AC,color='#243449',ls=':',lw=1,label=r'$A_c$ reference')
    ax.set(xlabel='Initial Gaussian amplitude, A',ylabel='Cutoff time / observation horizon',title='Full amplitude range: subcritical data and strongly supercritical collapse')
    ax.legend(fontsize=8);fig.savefig(OUT/'amplitude_overview.png',dpi=300);plt.close(fig)
    if gaps:
        count=len(gaps);cols=min(3,count);nr=math.ceil(count/cols)
        fig,axes=plt.subplots(nr,cols,figsize=(4.4*cols,3*nr),layout='constrained',squeeze=False)
        for ax,w in zip(axes.flat,gaps):
            lo=float(w['left_outer_delta']);hi=float(w['right_outer_delta']);span=max(hi-lo,hi*1e-6);lo=max(1e-8,lo-.4*span);hi+=.4*span
            chosen=[r for r in positive if lo<=r['relative_offset']<=hi]
            for sign,marker,color in [(1,'o','#1f77b4'),(-1,'v','#e2642a'),(0,'x','#8064a2')]:
                sub=[r for r in chosen if r['event_sign']==sign]
                if sign:ax.scatter([r['relative_offset'] for r in sub],[r['t12'] for r in sub],s=9,marker=marker,color=color)
                else:ax.scatter([r['relative_offset'] for r in sub],[.94]*len(sub),transform=ax.get_xaxis_transform(),s=12,marker=marker,color=color)
            ax.set(xlabel=r'Raw excess, $\delta$',ylabel=r'$t_{12}$',title=f'{w["window"]}: raw delta near {w["representative_delta"]:.4g}')
            ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0),useOffset=False)
            ax.xaxis.set_major_locator(MaxNLocator(4))
        for ax in list(axes.flat)[count:]:ax.set_visible(False)
        fig.suptitle('Refined gaps: crosses mark no crossing to 48 by t=100',fontsize=12,fontweight='bold')
        fig.savefig(OUT/'gap_details.png',dpi=300);plt.close(fig)
    trajectory_file=DATA/'representative_trajectories.npz'
    if validation_current and trajectory_file.exists():
        arrays=np.load(trajectory_file);count=len(windows);cols=3;nr=math.ceil(count/cols)
        fig,axes=plt.subplots(nr,cols,figsize=(13,2.5*nr),layout='constrained',squeeze=False)
        for i,(ax,w) in enumerate(zip(axes.flat,windows)):
            d=arrays[f'group_{i+1:02d}'];ax.plot(d[:,0],d[:,1],lw=.9,color='#1f77b4')
            ax.set(title=f'{w["window"]} | raw delta={w["representative_delta"]:.4g}',xlabel='Time, t',ylabel=r'$\phi(0,t)$')
            if w['kind']=='noncrossing_window':ax.set_yscale('symlog',linthresh=.01)
            else:ax.set_ylim(-12,12)
        for ax in list(axes.flat)[count:]:ax.set_visible(False)
        fig.suptitle('Representative evolutions on the finer N=1997 grid\nCollapse panels show the field range -12 to +12; gaps show the complete t=100 run',fontsize=12,fontweight='bold')
        fig.savefig(OUT/'representative_trajectories.png',dpi=300);plt.close(fig)
    confirmed=sum(w['resolution_confirmed']=='1' for w in collapse_windows)
    near_rows=[r for r in CACHE.values() if r['n']==1997 and r['cfl']==.4 and r['horizon']==100 and r['radius']==25 and 0<r['relative_offset']<=1.000001e-8]
    settings=dict(near_critical_extension_runs=len(near_rows),near_critical_extension_states=sorted({state(r) for r in near_rows}),reference=AC,reference_error=AC_ERROR,amplitude_min=min(r['amplitude'] for r in rows),amplitude_max=max(r['amplitude'] for r in rows),delta_min_positive=positive[0]['relative_offset'],delta_max=positive[-1]['relative_offset'],n_primary=999,spatial_order=4,cfl=.4,horizon=100,classification_level=48,plotted_time_level=12,primary_runs=len(rows),positive_runs=len(positive),transient_12_crossings_without_48=sum(not r['collapsed'] and r['t12']>=0 for r in rows),sampled_collapse_bands=wid,noncrossing_windows=len(gaps),confirmed_collapse_band_interiors=confirmed,confirmed_noncrossing_window_interiors=sum(w['resolution_confirmed']=='1' for w in gaps),total_cached_configurations=len(CACHE),source_sha256=SOURCE)
    (OUT/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    report(settings,windows,checks)
    guide=ROOT/'README.md'
    if guide.exists():
        block=f'<!-- survey-summary -->\n**Current result: {len(gaps)} non-crossing windows and {wid} collapse bands**, from {len(rows)} primary amplitudes. All counts are finite-horizon numerical observations. Validated window interiors: {settings["confirmed_noncrossing_window_interiors"]}/{len(gaps)}.\n<!-- /survey-summary -->'
        body=guide.read_text()
        if '<!-- survey-summary -->' in body:body=re.sub(r'<!-- survey-summary -->.*?<!-- /survey-summary -->',lambda m:block,body,flags=re.S)
        else:body=body.replace('# Gaussian collapse windows\n','# Gaussian collapse windows\n\n'+block+'\n',1)
        guide.write_text(body)
    print(json.dumps(settings,indent=2),flush=True)

def report(settings,windows,checks):
    collapse=[w for w in windows if w['kind']=='collapse_band']
    text='# Broad amplitude search: resolved windows\n\n'
    text+=f'**{settings["noncrossing_windows"]} finite-horizon non-crossing windows and {len(collapse)} collapse bands were sampled. {settings["confirmed_noncrossing_window_interiors"]} window interiors pass the finer-grid and larger-domain checks; {settings["confirmed_collapse_band_interiors"]} collapse-band interiors agree across N=999,1997,3993.**\n\n'
    text+=f'The search covers A={settings["amplitude_min"]:.6g} to {settings["amplitude_max"]:.6g}, including subcritical starts and raw positive offsets {settings["delta_min_positive"]:.3g} to {settings["delta_max"]:.3g}. There are {settings["primary_runs"]} primary configurations, with adaptive sampling around state changes and two extra passes through non-crossing gaps.\n\n'
    text+='A non-crossing window (G labels) is a sampled interval where field 48 is not reached by t=100, bounded by collapse samples on both sides. A collapse band (B labels) has a fixed final field sign and number of central-field zero crossings. Both are ordered by amplitude. A sign change or oscillation-count change is reported as a branch change, separately from a non-crossing window. This finite scan is not a proof of the exact number of continuum windows; smaller unresolved windows may exist.\n\n'
    text+='The reference is Ac=1.2026097332023025 with estimated uncertainty 5e-12. All offset values are raw `(A-Ac)/Ac`, never amplitude percentages. The primary method is fourth order at N=999, CFL=0.4, Rmax=25, sponge starts at 20, R0=2, zero initial velocity. Collapse classification uses field 48; the plot retains t12 to match the original cutoff-time figure. Neither t12 nor t48 is the singularity time tc.\n\n'
    text+=f'Among primary non-crossing runs, {settings["transient_12_crossings_without_48"]} reached field 12 without reaching 48. A zero count means the observed classification also agrees with the earlier field-12 criterion.\n\n'
    text+='## Collapse bands\n\n| Band | Sampled raw delta interval | Final sign | Zero crossings | Finer-grid interior confirmed |\n| --- | --- | ---: | ---: | --- |\n'
    for w in collapse:text+=f'| {w["window"]} | [{w["delta_first"]:.9g}, {w["delta_last"]:.9g}] | {w["event_sign"]:+d} | {w["zero_crossings"]} | {w["resolution_confirmed"]} |\n'
    text+='\n## Non-crossing windows\n\n| Window | Sampled raw delta interval | Validation confirmed |\n| --- | --- | --- |\n'
    for w in windows:
        if w['kind']=='noncrossing_window':text+=f'| {w["window"]} | [{w["delta_first"]:.9g}, {w["delta_last"]:.9g}] | {w["resolution_confirmed"]} |\n'
    text+='\nThe endpoints above are the first/last sampled amplitudes of each group. `windows.csv` also lists the neighboring outside samples that bracket each boundary. They are N=999 sampled boundaries, not extrapolated continuum edge locations. Resolution checks validate representative interiors, not every boundary or the exact width of every window.\n\n'
    text+=f'\nThe extra near-critical extension contains {settings["near_critical_extension_runs"]} runs, with observed states {settings["near_critical_extension_states"]}. State +1/6 means positive-field collapse after six central-field zero crossings. No additional window was detected in this extension if only that state is present.\n\n'
    text+='## Validation and limitations\n\nEach collapse-window representative is rerun at the same amplitude on N=1997 and N=3993. Each gap representative is checked on N=1997 and in a larger domain (radius at least 75, sponge starts at 60) through t=200. `resolution_checks.csv` records exact agreement or disagreement; a discrepancy is retained rather than hidden. Extra near-critical runs at N=1997 cover raw delta=1e-10..1e-8 separately; they do not inflate the main window count. A no-crossing result means only no field-48 event by the observation horizon.\n\n'
    text+='The tested range is broader than the old endpoint 1.2 Ac: it extends to 4 Ac. Sampling is logarithmic near criticality, supplemented by a linear broad scan. Adaptive subdivision probes changes in sign/oscillation count and fills known gaps; a finite search cannot rule out arbitrarily narrow islands between equal-outcome samples.\n\n'
    text+='## Files\n\n- `window_map.png`: main cutoff-time/log(raw offset) figure and oscillation counts.\n- `amplitude_overview.png`: the full linear amplitude range.\n- `gap_details.png`: individual non-crossing intervals.\n- `representative_trajectories.png`: finer-grid examples.\n- `amplitude_samples.csv`: all primary amplitudes and outcomes.\n- `windows.csv`: sampled window extents and exterior boundary brackets.\n- `resolution_checks.csv`: finer-grid and larger-domain checks.\n- `../data/runs.jsonl`: complete cached computations, including validation runs.\n- `../data/representative_trajectories.npz`: compressed representative data; columns time, phi_center, pi_center, r50.\n'
    (OUT/'REPORT.md').write_text(text)

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--stage',choices=['all','survey','refine','edges','validate','plots'],default='all');args=parser.parse_args();build()
    if args.stage in ['all','survey']:survey()
    if args.stage in ['all','refine']:refine()
    if args.stage=='edges':refine_edges()
    if args.stage in ['all','validate']:validate()
    analyze()
if __name__=='__main__':main()
