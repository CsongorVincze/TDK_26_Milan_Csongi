#!/usr/bin/env python3
"""Reproducible Module 4 threshold search and dense 0..20% sweep."""
import csv
import json
import os
from pathlib import Path
import subprocess
from concurrent.futures import ThreadPoolExecutor
import sys
import importlib.util
os.environ.setdefault('MPLCONFIGDIR', '/tmp/supercritical-mpl')
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
HERE = Path(__file__).resolve().parent
OUT = HERE / 'results'
OUT.mkdir(exist_ok=True)

def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def simulate(a, n=500, cfl=.4, horizon=100, path=None):
    result = subprocess.run([str(HERE/'solver'), repr(float(a)), str(n), str(cfl),
                             str(horizon), '.02', str(path) if path else '-'],
                            check=True, capture_output=True, text=True)
    event, time, steps = result.stdout.split()
    return dict(amplitude=float(a), n=n, cfl=cfl, horizon=horizon,
                collapsed=int(event), event_time=float(time) if int(event) else None,
                end_time=float(time), steps=int(steps))

def search(n=500, cfl=.4):
    trials=[]
    def classify(a):
        row=simulate(a,n,cfl)
        trials.append(row)
        print(f'search N={n} CFL={cfl}: A={a:.12f}, event={row["collapsed"]}, end={row["end_time"]:.6f}', flush=True)
        return row['collapsed']
    lo,hi=.1,1.8
    assert not classify(lo) and classify(hi)
    while hi-lo>1e-10:
        mid=(lo+hi)/2
        if classify(mid): hi=mid
        else: lo=mid
    write_csv(OUT/f'critical_search_N{n}_cfl{cfl}.csv',trials)
    return dict(n=n,cfl=cfl,lower_noncollapsing=lo,upper_collapsing=hi,
                midpoint=(lo+hi)/2, bracket_width=hi-lo)

def load_saved_rows():
    with (OUT/'sweep_summary.csv').open() as handle:
        rows=list(csv.DictReader(handle))
    for row in rows:
        for key in ['amplitude','cfl','horizon','end_time','relative_offset','offset_percent']:
            row[key]=float(row[key])
        for key in ['n','collapsed','steps','event_sign']:
            row[key]=int(row[key])
        row['event_time']=float(row['event_time']) if row['event_time'] else None
    return rows


def export_figure(fig, name):
    fig.savefig(OUT/f'{name}.png',dpi=300)
    plt.close(fig)


def plot(rows, critical):
    """Apply TDK_plot_skill to saved results; no numerical data are changed."""
    from matplotlib.ticker import LogFormatterMathtext
    plt.style.use(HERE/'tdk.mplstyle')
    acrit=critical['midpoint']
    positives=[r for r in rows if r['relative_offset']>0]
    norm=LogNorm(positives[0]['relative_offset'],.2)
    cmap=plt.get_cmap('viridis')
    max_event=max(r['event_time'] for r in rows if r['collapsed'])
    missing=[r for r in positives if not r['collapsed']]
    settings=r'3D Gaussian initial data  |  $N=500$, $R_0=2$  |  $|\phi(0,t)|=12$ cutoff'
    fig,axes=plt.subplots(1,2,figsize=(14,5.4),layout='constrained')
    fig.suptitle(f'Supercritical amplitude sweep\n'+settings,fontsize=12,fontweight='bold')
    for row in positives:
        color=cmap(norm(row['relative_offset']))
        data=np.genfromtxt(OUT/row['trajectory'],delimiter=',',names=True)
        axes[0].plot(data['time'],data['phi_center'],color=color,lw=.55,alpha=.65,rasterized=True)
        if row['collapsed']:
            axes[1].scatter(row['relative_offset'],row['event_time'],color=color,s=11,linewidths=0,zorder=3)
        else:
            axes[1].scatter(row['relative_offset'],.95,color=color,s=25,marker='x',
                            transform=axes[1].get_xaxis_transform(),zorder=4)
    data=np.genfromtxt(OUT/rows[0]['trajectory'],delimiter=',',names=True)
    axes[0].plot(data['time'],data['phi_center'],color='#243449',lw=1.6,
                 label=r'$A=A_{\mathrm{crit}}$ estimate',zorder=5)
    for level in [-12,12]:
        axes[0].axhline(level,color='#e2642a',ls=':',lw=1.1)
    axes[0].set(xlim=(0,1.05*max_event),ylim=(-13.2,13.2),
                xlabel=r'Time, $t$ (dimensionless)',ylabel=r'Central field, $\phi(0,t)$',
                title='(a) Central-field trajectories')
    axes[0].legend(loc='lower left',fontsize=8)
    axes[1].set(xscale='log',ylim=(0,1.19*max_event),
                xlabel=r'Amplitude excess, $(A-A_{\mathrm{crit}})/A_{\mathrm{crit}}$ (dimensionless)',
                ylabel=r'Cutoff time, $t_{12}$ (dimensionless)',
                title='(b) Time to the collapse cutoff')
    axes[1].xaxis.set_major_formatter(LogFormatterMathtext())
    axes[1].grid(which='minor',alpha=.18,lw=.4)
    axes[1].text(.025,.985,r'$\times$: no crossing by $t=100$ (upper band)',
                 transform=axes[1].transAxes,va='top',fontsize=8)
    axes[1].text(.97,.74,f'$A_{{\\mathrm{{crit}}}}={acrit:.12f}$\nLocal numerical threshold',
                 transform=axes[1].transAxes,ha='right',va='bottom',fontsize=8,
                 bbox=dict(boxstyle='round,pad=.45',facecolor='white',edgecolor='#d6dfea',lw=.7))
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=cmap),ax=axes[0],pad=.025,fraction=.045)
    cb.set_label('Amplitude excess (dimensionless)',fontsize=9)
    cb.ax.tick_params(labelsize=8)
    cb.outline.set_linewidth(.5)
    export_figure(fig,'supercritical_sweep')

    fig,ax=plt.subplots(figsize=(10,4.9),layout='constrained')
    for sign,marker,color,label in [(1,'o','#1f77b4',r'Positive: $\phi(0,t_{12})=+12$'),
                                    (-1,'v','#e2642a',r'Negative: $\phi(0,t_{12})=-12$')]:
        chosen=[r for r in positives if r['collapsed'] and r['event_sign']==sign]
        ax.scatter([r['relative_offset'] for r in chosen],[r['event_time'] for r in chosen],
                   s=15,marker=marker,color=color,label=label,linewidths=0,zorder=3)
    if missing:
        ax.scatter([r['relative_offset'] for r in missing],[.95]*len(missing),
                   marker='x',s=28,color='#8064a2',transform=ax.get_xaxis_transform(),zorder=4)
    ax.set(xscale='log',ylim=(0,1.2*max_event),xlabel='Amplitude excess (dimensionless)',
           ylabel=r'Cutoff time, $t_{12}$ (dimensionless)',title='Collapse-cutoff time and final field sign')
    ax.text(.02,.985,r'$\times$: no crossing by $t=100$ (upper band)',transform=ax.transAxes,va='top',fontsize=8)
    ax.grid(which='minor',alpha=.18,lw=.4)
    ax.legend(loc='upper right',bbox_to_anchor=(.99,.88),fontsize=8)
    export_figure(fig,'collapse_time_detail')
    plot_noncrossing(rows)


def plot_noncrossing(rows):
    """Separate non-crossing families so their full-horizon behavior is legible."""
    plt.style.use(HERE/'tdk.mplstyle')
    missing=[r for r in rows if not r['collapsed']]
    groups=[([r for r in missing if r['relative_offset']==0], '(a) Threshold estimate', '#243449'),
            ([r for r in missing if 0<r['relative_offset']<.0001], '(b) First non-crossing band', '#1f77b4'),
            ([r for r in missing if r['relative_offset']>=.0001], '(c) Second non-crossing band', '#8064a2')]
    fig,axes=plt.subplots(1,3,figsize=(14,4.1),layout='constrained',sharex=True,sharey=True)
    fig.suptitle(r'Non-crossing trajectories through $t=100$'+'\nNo run shown reaches '+r'$|\phi(0,t)|=12$',fontsize=12,fontweight='bold')
    for ax,(group,title,color) in zip(axes,groups):
        for row in group:
            data=np.genfromtxt(OUT/row['trajectory'],delimiter=',',names=True)
            ax.plot(data['time'],data['phi_center'],color=color,lw=1 if len(group)==1 else .65,
                    alpha=1 if len(group)==1 else .4)
        ax.set(title=title,xlabel=r'Time, $t$ (dimensionless)',xlim=(0,100),ylim=(-13,13))
        for level in [-12,12]:ax.axhline(level,color='#e2642a',ls=':',lw=.8)
        if group:
            offsets=[r['relative_offset'] for r in group]
            label=(r'$A=A_{\mathrm{crit}}$ estimate' if len(group)==1 else
                   f'{len(group)} runs | excess {min(offsets):.6g}–{max(offsets):.6g}')
            ax.text(.04,.96,label,transform=ax.transAxes,va='top',fontsize=8,
                    bbox=dict(boxstyle='round,pad=.4',facecolor='white',edgecolor='#d6dfea',lw=.6))
    axes[0].set_ylabel(r'Central field, $\phi(0,t)$')
    export_figure(fig,'noncrossing_trajectories')

def main():
    subprocess.run(['g++','-O3','-std=c++17',str(HERE/'solver.cpp'),'-o',str(HERE/'solver')],check=True)
    # Independent reference implementation check, without output-dependent step clipping.
    spec=importlib.util.spec_from_file_location('module4_reference',HERE.parent/'module4/milan_alapkod_module4.py')
    mod=importlib.util.module_from_spec(spec);sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
    validation=[]
    for a in [1.203812310596742,1.3]:
        c=mod.Config();g=mod.make_grid(c)
        event,time,steps,_=mod.simulate(a,c,g,100,save=False)
        fast=simulate(a)
        error=abs(time-fast['event_time'])
        assert event==bool(fast['collapsed']) and error<1e-7
        validation.append(dict(amplitude=a,reference_time=time,cpp_time=fast['event_time'],absolute_difference=error))
    write_csv(OUT/'solver_validation.csv',validation)
    critical=search()
    (OUT/'critical_summary.json').write_text(json.dumps(critical,indent=2)+'\n')
    acrit=critical['midpoint']
    # 181 logarithmic near-critical points, 161 linear broad points, plus endpoints.
    offsets=np.unique(np.r_[0.,np.geomspace(1e-9,.01,181),np.linspace(.01,.2,161)])
    def work(item):
        i,offset=item
        path=Path('trajectories')/f'run_{i:04d}.csv'
        row=simulate(acrit*(1+offset),path=OUT/path)
        data=np.genfromtxt(OUT/path,delimiter=',',names=True)
        row.update(relative_offset=float(offset),offset_percent=float(100*offset),
                   event_sign=int(np.sign(data['phi_center'][-1])) if row['collapsed'] else 0,
                   trajectory=str(path))
        return row
    (OUT/'trajectories').mkdir(exist_ok=True)
    rows=[]
    with ThreadPoolExecutor(max_workers=6) as pool:
        for row in pool.map(work,enumerate(offsets)):
            rows.append(row)
            if len(rows)%20==0:print(f'sweep {len(rows)}/{len(offsets)}',flush=True)
    # Resolve narrow non-crossing features, preserving every completed run.
    extra=[]
    for row in rows:
        if row['relative_offset']>0 and not row['collapsed']:
            extra.extend(np.linspace(.9*row['relative_offset'],1.1*row['relative_offset'],81))
    existing={r['relative_offset'] for r in rows}
    extra=sorted(set(extra)-existing)
    with ThreadPoolExecutor(max_workers=6) as pool:
        rows.extend(pool.map(work,enumerate(extra,start=len(rows))))
    rows.sort(key=lambda r:r['relative_offset'])
    write_csv(OUT/'sweep_summary.csv',rows)
    plot(rows,critical)
    for label in ['lower_noncollapsing','upper_collapsing']:
        simulate(critical[label],horizon=200,path=OUT/f'critical_{label}.csv')
    # Check the finite-horizon bracket and time-step sensitivity at its two edges.
    checks=[]
    for cfl in [.4,.2]:
        for horizon in [100,200]:
            for a in [critical['lower_noncollapsing'],critical['upper_collapsing']]:
                checks.append(simulate(a,cfl=cfl,horizon=horizon))
    write_csv(OUT/'critical_sensitivity.csv',checks)
    settings=dict(model='phi_tt=radial_laplacian(phi)-phi+phi^3-gamma*phi_t',n=500,
                  r_max=25,stretch=5,r0=2,initial_profile='A exp(-r^2/(2 R0^2))',initial_velocity=0,
                  cfl=.4,sponge_fraction=.8,sponge_strength=7,sponge_power=3,
                  critical_level=12,critical_horizon=100,critical_tolerance=1e-10,
                  output_interval=.02,output_method='linear interpolation; no step modification',
                  sweep_count=len(rows),maximum_relative_offset=.2,
                  interpretation='Local finite-horizon numerical transition, not a continuum critical amplitude')
    (OUT/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
    print(json.dumps(dict(critical=critical,runs=len(rows),noncrossing=sum(not r['collapsed'] for r in rows)),indent=2),flush=True)
if __name__=='__main__':
    if sys.argv[1:]==['--plot-only']:
        plot(load_saved_rows(),json.loads((OUT/'critical_summary.json').read_text()))
    else:
        main()
