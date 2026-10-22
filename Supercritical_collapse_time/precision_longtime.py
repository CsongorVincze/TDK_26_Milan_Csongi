#!/usr/bin/env python3
"""Local threshold precision audit, close-offset scan and t=100 inspection."""
import csv,json,subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LogNorm
import study as s
OUT=s.OUT/'precision_longtime'
OUT.mkdir(exist_ok=True)

def refine(n,cfl,lo,hi,tol):
    trials=[]
    def f(a):
        r=s.simulate(a,n=n,cfl=cfl);trials.append(r)
        return r['collapsed']
    assert not f(lo) and f(hi), (n,cfl,lo,hi)
    while hi-lo>tol:
        mid=(lo+hi)/2
        if mid==lo or mid==hi:break
        if f(mid):hi=mid
        else:lo=mid
    s.write_csv(OUT/f'search_N{n}_cfl{cfl}.csv',trials)
    result=dict(n=n,cfl=cfl,lower_noncollapsing=lo,upper_collapsing=hi,
                midpoint=(lo+hi)/2,bracket_width=hi-lo,half_width=(hi-lo)/2,
                relative_half_width=(hi-lo)/(hi+lo),horizon=100)
    print('precision',result,flush=True)
    return result

def figures(old, reruns, near, precision, representatives):
    plt.style.use(s.HERE/'tdk.mplstyle')
    ac=precision[0]['midpoint']
    # Full trajectories for every original run, without truncating at t=20.
    fig,axes=plt.subplots(2,1,figsize=(11,7),layout='constrained')
    norm=LogNorm(1e-9,.2);cmap=plt.get_cmap('viridis')
    for r in old:
        d=np.genfromtxt(s.OUT/r['trajectory'],delimiter=',',names=True)
        color='#243449' if r['relative_offset']==0 else cmap(norm(r['relative_offset']))
        axes[0].plot(d['time'],d['phi_center'],color=color,lw=.55,alpha=.5)
    for r in reruns:
        d=np.genfromtxt(OUT/r['trajectory'],delimiter=',',names=True)
        axes[1].plot(d['time'],d['phi_center'],lw=.6,alpha=.5,
                     color='#243449' if r['relative_offset']==0 else cmap(norm(r['relative_offset'])))
    axes[0].set(title='(a) All 502 original amplitudes: full evolution',xlim=(0,100),ylim=(-13,13),ylabel=r'Central field, $\phi(0,t)$')
    axes[1].set(title='(b) Runs surviving past t = 20: late-time detail',xlim=(20,100),
                yscale='symlog',ylabel=r'$\phi(0,t)$ (symmetric log)',xlabel=r'Time, $t$ (dimensionless)')
    axes[1].set_yscale('symlog',linthresh=.01)
    fig.suptitle('Evolution through t = 100; trajectories stop only at the field cutoff',fontsize=12,fontweight='bold')
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes,pad=.02,fraction=.025)
    cb.set_label(r'Raw offset, $(A-A_c)/A_c$ (original reference)')
    fig.savefig(OUT/'full_evolution_t100.png',dpi=300);plt.close(fig)
    # Reuse existing sampled amplitudes after transforming to the refined reference.
    near=list(near)
    existing={r['amplitude'] for r in near}
    for row in old:
        delta=(row['amplitude']-ac)/ac
        if 0<delta<=1e-3 and row['amplitude'] not in existing:
            near.append(dict(row,relative_offset=delta))
    near.sort(key=lambda r:r['relative_offset'])
    # New close-offset sweep and uncertainty scale.
    fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained')
    for sign,marker,color in [(-1,'v','#e2642a'),(1,'o','#1f77b4')]:
        group=[r for r in near if r['relative_offset']>0 and r['collapsed'] and r['event_sign']==sign]
        ax.scatter([r['relative_offset'] for r in group],[r['event_time'] for r in group],s=15,marker=marker,color=color,label=f'Field sign {sign:+d}')
    missing=[r for r in near if r['relative_offset']>0 and not r['collapsed']]
    if missing:ax.scatter([r['relative_offset'] for r in missing],[.95]*len(missing),transform=ax.get_xaxis_transform(),marker='x',color='#8064a2',label='No crossing by t = 100 (upper band)')
    oldc=json.loads((s.OUT/'critical_summary.json').read_text())
    olderr=oldc['bracket_width']/2/oldc['midpoint']
    ax.axvline(olderr,color='#8064a2',ls='--',lw=1,label='Previous relative bracket half-width')
    shifts=[abs(p['midpoint']/ac-1) for p in precision[1:] if p['n']==500]
    for value,label in [(precision[0]['relative_half_width'],'Refined bracket half-width'),(max(shifts),'Largest tested time-step shift')]:
        ax.axvline(value,ls=':',lw=1,label=label)
    ax.set(xscale='log',xlabel=r'Raw offset, $\delta=(A-A_c)/A_c$',ylabel=r'Cutoff time, $t_{12}$',title='Near-critical scan with the refined N = 500 reference')
    ax.legend(fontsize=7,loc='lower left');ax.grid(which='minor',alpha=.15)
    fig.savefig(OUT/'nearcritical_precision.png',dpi=300);plt.close(fig)
    # Spatial evolution shows whether a central transient travels away.
    fig,axes=plt.subplots(3,2,figsize=(12,9),layout='constrained',sharex=True,sharey=True)
    radius=25*np.sinh(5*np.linspace(0,1,500))/np.sinh(5);use=radius<=8
    for ax,rep in zip(axes.flat,representatives):
        field=np.genfromtxt(OUT/rep['fields'],delimiter=',',skip_header=1)
        im=ax.pcolormesh(field[:,0],radius[use],field[:,1:][:,use].T,cmap='RdBu_r',
                         norm=SymLogNorm(linthresh=.01,vmin=-12,vmax=12),shading='auto',rasterized=True)
        ax.set(title=rep['label'],xlim=(0,100),ylim=(0,8))
    for ax in axes[-1]:ax.set_xlabel('Time, t')
    for ax in axes[:,0]:ax.set_ylabel('Radius, r')
    cb=fig.colorbar(im,ax=axes,pad=.02,fraction=.025,extend='both');cb.set_label(r'Field, $\phi(r,t)$ (symmetric log)')
    fig.suptitle('Spatial evolution to t = 100 (blank after cutoff)',fontsize=12,fontweight='bold')
    for ax in axes.flat:ax.grid(False)
    fig.savefig(OUT/'spatial_evolution_t100.png',dpi=300);plt.close(fig)
    fig,axes=plt.subplots(3,1,figsize=(10,7),layout='constrained',sharex=True,sharey=True)
    for ax,rep in zip(axes,[representatives[i] for i in [0,2,3]]):
        d=np.genfromtxt(OUT/rep['trajectory'],delimiter=',',names=True)
        ax.plot(d['time'],d['phi_center'],lw=1,color='#1f77b4')
        ax.set_yscale('symlog',linthresh=.01)
        ax.set(title=rep['label'],ylabel=r'$\phi(0,t)$',xlim=(0,100))
        ax.axvline(20,color='#8064a2',ls=':',lw=.8)
    axes[-1].set_xlabel('Time, t (dimensionless)')
    fig.suptitle('Three surviving trajectories: symmetric-log field scale',fontsize=12,fontweight='bold')
    fig.savefig(OUT/'surviving_oscillons_t100.png',dpi=300);plt.close(fig)

def main():
    subprocess.run(['g++','-O3','-std=c++17',str(s.HERE/'solver.cpp'),'-o',str(s.HERE/'solver')],check=True)
    c=json.loads((s.OUT/'critical_summary.json').read_text())
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures=[pool.submit(refine,500,cfl,c['lower_noncollapsing'],c['upper_collapsing'],1e-14) for cfl in [.4,.2,.1]]
        precision=[f.result() for f in futures]
    # Same broad seed, reported as a bisection-selected local transition, not global.
    precision.append(refine(1000,.4,.1,1.8,1e-9))
    s.write_csv(OUT/'precision_summary.csv',precision)
    ac=precision[0]['midpoint']
    old=s.load_saved_rows()
    survivors=[r for r in old if r['end_time']>20]
    (OUT/'trajectories').mkdir(exist_ok=True)
    def rerun(item):
        i,r=item;path=f'trajectories/survivor_{i:03d}.csv'
        new=s.simulate(r['amplitude'],horizon=100,path=OUT/path)
        new.update(relative_offset=r['relative_offset'],trajectory=path,source_trajectory=r['trajectory'])
        assert new['collapsed']==r['collapsed']
        if new['collapsed']:assert abs(new['event_time']-r['event_time'])<1e-9
        return new
    with ThreadPoolExecutor(max_workers=6) as pool:reruns=list(pool.map(rerun,enumerate(survivors)))
    s.write_csv(OUT/'survivors_t100.csv',reruns)
    print('survivors rerun',len(reruns),flush=True)
    offsets=np.unique(np.r_[0,np.geomspace(1e-13,1e-3,201)])
    def close_run(item):
        i,delta=item;path=f'trajectories/near_{i:03d}.csv';a=float(ac*(1+delta))
        new=s.simulate(a,horizon=100,path=OUT/path)
        d=np.genfromtxt(OUT/path,delimiter=',',names=True)
        new.update(requested_offset=float(delta),relative_offset=(a-ac)/ac,trajectory=path,
                   event_sign=int(np.sign(d['phi_center'][-1])) if new['collapsed'] else 0)
        return new
    with ThreadPoolExecutor(max_workers=6) as pool:near=list(pool.map(close_run,enumerate(offsets)))
    s.write_csv(OUT/'nearcritical_sweep.csv',near)
    print('nearcritical done',len(near),flush=True)
    noncross=[r for r in old if not r['collapsed'] and r['relative_offset']>0]
    bands=[[r for r in noncross if r['relative_offset']<1e-4],[r for r in noncross if r['relative_offset']>=1e-4]]
    choices=[(c['midpoint'],'Original threshold estimate'),(c['upper_collapsing'],'Original upper edge')]
    for i,band in enumerate(bands,1):
        choices.append((band[len(band)//2]['amplitude'],f'Non-crossing band {i}: representative'))
    choices.extend([(ac*(1+1e-10),'Refined reference + raw offset 1e-10'),(ac*(1+1e-6),'Refined reference + raw offset 1e-6')])
    def field_run(item):
        i,(a,label)=item;fields=f'fields_{i}.csv';path=f'trajectories/spatial_{i}.csv'
        proc=subprocess.run([str(s.HERE/'solver'),repr(a),'500','.4','100','.02',str(OUT/path),str(OUT/fields)],capture_output=True,text=True,check=True)
        event,t,steps=proc.stdout.split()
        return dict(amplitude=a,relative_offset=(a-ac)/ac,label=label,fields=fields,trajectory=path,collapsed=int(event),end_time=float(t))
    with ThreadPoolExecutor(max_workers=6) as pool:representatives=list(pool.map(field_run,enumerate(choices)))
    s.write_csv(OUT/'spatial_runs.csv',representatives)
    np.savetxt(OUT/'spatial_radius.csv',25*np.sinh(5*np.linspace(0,1,500))/np.sinh(5),delimiter=',',header='radius',comments='')
    figures(old,reruns,near,precision,representatives)
    report=dict(refined_reference=ac,precision=precision,old_reference=c['midpoint'],
                survivor_runs=len(reruns),survivors_reaching_horizon=sum(not r['collapsed'] for r in reruns),
                nearcritical_runs=len(near),nearcritical_noncrossing=sum(not r['collapsed'] for r in near),
                no_pdf=True,horizon=100)
    (OUT/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
if __name__=='__main__':
    import sys
    if sys.argv[1:]==['--plot-only']:
        def read(name):
            with (OUT/name).open() as f:rows=list(csv.DictReader(f))
            for r in rows:
                for k,v in r.items():
                    if v=='':r[k]=None
                    else:
                        try:r[k]=float(v)
                        except ValueError:pass
            return rows
        figures(s.load_saved_rows(),read('survivors_t100.csv'),read('nearcritical_sweep.csv'),read('precision_summary.csv'),read('spatial_runs.csv'))
    else:main()
