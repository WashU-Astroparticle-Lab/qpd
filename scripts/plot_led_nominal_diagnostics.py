"""Plot controlled nominal efficiency/purity sweeps and fixed 0.5s example."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from qpd.mlebench.led_generalization_score import timing_settings
from qpd.mlebench.led_reporting import strict_event_count, c1_truth_states, preservation_counts

def estimate(rows, method, denominator):
    numerator = np.array([r['strict_counts'][method]['native'] for r in rows])
    den = np.array([r['metrics'][method][denominator] for r in rows])
    value = 100*numerator.sum()/den.sum() if den.sum() else np.nan
    ix = np.random.default_rng(71116001).integers(0, len(rows), (1000, len(rows)))
    ns, ds = numerator[ix].sum(axis=1), den[ix].sum(axis=1)
    bootstrap = 100*ns[ds>0]/ds[ds>0]
    interval = np.percentile(bootstrap, [2.5,97.5]) if len(bootstrap) else [np.nan,np.nan]
    return value, interval, int(numerator.sum()), int(den.sum())


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    a=p.parse_args()
    plan=json.loads((a.data/'plan.json').read_text())
    assert json.loads((a.data/'COMPLETE.json').read_text())['complete']
    records=json.loads((a.data/'records.json').read_text())
    if plan.get('split') != 'display':
        p.error('Use evaluation results from the fixed display split')
    methods = [(m['key'], m['label'], m['color'], m['marker']) for m in plan['plot_methods']]
    fs = plan['nominal']['sample_rate_hz']
    nominal_snr = plan['nominal']['snr_max']
    nominal_mean = plan['nominal']['mean_tunnels']
    tol = timing_settings(fs)[0]
    # Compute both strict diagnostics directly from archived event lists.
    # The primary C1 assignment supplies the late/precise truth states.
    paired = {}
    for row in records:
        truth = row['truth']
        row['strict_counts'] = {}
        states = {}
        for key, _, _, _ in methods:
            pred = row['predictions'][key]
            row['strict_counts'][key] = dict(native=strict_event_count(truth,pred,tol),
                                              fixed50=strict_event_count(truth,pred,50e-6))
            states[key] = c1_truth_states(truth,pred,fs)
        condition = f"mean={row['mean_tunnels']:g},snr={row['snr_max']:g}"
        point = paired.setdefault(condition,{})
        for key in states:
            totals = point.setdefault(key, dict(records=0, truth=0, predictions=0,
                seconds=0., c1_fp=0, strict=0, strict50=0, c1_states_vs_hmmsf=np.zeros((3,3),dtype=int)))
            totals['records'] += 1
            totals['truth'] += len(truth)
            totals['predictions'] += len(row['predictions'][key])
            totals['c1_fp'] += row['metrics'][key]['false_positives']
            totals['seconds'] += row['metrics'][key]['duration_s']
            totals['strict'] += row['strict_counts'][key]['native']
            totals['strict50'] += row['strict_counts'][key]['fixed50']
            totals['c1_states_vs_hmmsf'] += preservation_counts(states['hmmsf'],states[key])
    for point in paired.values():
        for totals in point.values():
            totals['c1_states_vs_hmmsf'] = totals['c1_states_vs_hmmsf'].tolist()
            totals['c1_fp_per_s'] = totals['c1_fp']/totals['seconds']
    a.out.mkdir(parents=True,exist_ok=True)
    (a.out/'paired_nominal.json').write_text(json.dumps(dict(
        state_order=['unmatched','late','precise'], rows='native HMM-SF',columns='compared method',
        strict_matching='maximum one-to-one count within inclusive tolerance',
        tolerance_s=tol,points=paired),indent=2)+'\n')
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    condition='10 kHz · 1 s records · LED ON/OFF 50 µs / 30 ms · background 20 Hz · resets 100/min'
    curves={}
    for metric,den in [('efficiency','n_truth'),('purity','n_pred')]:
        fig,axs=plt.subplots(1,2,figsize=(13,5.8),sharey=True)
        fig.subplots_adjust(top=.74,bottom=.26,left=.07,right=.98,wspace=.12)
        curves[metric]={}
        for ax,axis,grid,fixed,label in [
            (axs[0],'mean_tunnels',plan['mean_grid'],('snr_max',nominal_snr),'Mean LED flips per pulse'),
            (axs[1],'snr_max',plan['snr_grid'],('mean_tunnels',nominal_mean),'Maximum SNR')]:
            curves[metric][axis]={}
            for method,name,color,marker in methods:
                values=[];intervals=[];points=[]
                for value in grid:
                    selected=[r for r in records if r[axis]==value and r[fixed[0]]==fixed[1]]
                    assert len(selected)==plan['replicates'], (axis,value,len(selected))
                    y,ci,n,d=estimate(selected,method,den)
                    values.append(y);intervals.append(ci)
                    points.append(dict(x=value,percent=None if not np.isfinite(y) else y,
                        ci95=[None if not np.isfinite(v) else float(v) for v in ci],precise_matches=n,denominator=d))
                intervals=np.array(intervals)
                ax.plot(grid,values,color=color,marker=marker,label=name,lw=1.8,ls='--' if method.endswith('_previous') else '-')
                ax.fill_between(grid,intervals[:,0],intervals[:,1],color=color,alpha=.1)
                curves[metric][axis][method]=points
            ax.set_xscale('symlog',linthresh=.1) if axis=='mean_tunnels' else ax.set_xscale('log')
            ax.set_xticks(grid,[f'{x:g}' for x in grid]);ax.minorticks_off()
            ax.set_xlabel(label);ax.set_ylim(0,100);ax.grid(alpha=.2)
            ax.set_title(f'Max SNR fixed at {nominal_snr}' if axis=='mean_tunnels' else f'Mean LED flips fixed at {nominal_mean:g}',fontsize=12)
        axs[0].set_ylabel(metric.capitalize()+' (%)')
        fig.suptitle(f'Detection {metric} — otherwise nominal conditions',fontsize=18,y=.98)
        fig.text(.5,.905,condition,ha='center',fontsize=10)
        handles,labels=axs[0].get_legend_handles_labels()
        fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.87),ncol=len(methods),frameon=False)
        fig.text(.07,.16,'Strict one-to-one matches within ±200 µs (2 samples). All parity flips included.',fontsize=10)
        note = ' High purity may reflect very few detections; zero-prediction points are undefined.' if metric=='purity' else ''
        fig.text(.07,.115,f"{plan['replicates']} records/point; pooled counts; shading = 95% record-bootstrap interval."+note,fontsize=8.5)
        decoder_note=plan['plot_note']
        fig.text(.07,.07,decoder_note,fontsize=9)
        fig.text(.07,.025,'EMG μ/σ/τ = 0.2/0.4/0.8 ms; fixed onset 10 ms; identity IQ augmentation. Development diagnostic, not reserve evaluation.',fontsize=9)
        for ext in ['png','pdf']:fig.savefig(a.out/f'{metric}_nominal.{ext}',dpi=180)
        plt.close(fig)
    (a.out/'curves.json').write_text(json.dumps(curves,indent=2,allow_nan=False)+'\n')

    with np.load(a.data/'nominal_example.npz') as f:
        z=f['I'].astype(float)+1j*f['Q'].astype(float)
        keys=['truth']+[m[0] for m in methods]
        events={k:f[k].copy() for k in keys}
        pulses=f['pulse_times'].copy()
    fs=plan['nominal']['sample_rate_hz']
    noise=np.median(abs(np.diff(z)))/(2*np.sqrt(np.log(2)))
    z=(z-z.mean())/noise
    fig=plt.figure(figsize=(15,17),layout='constrained')
    gs=fig.add_gridspec(15,1,height_ratios=[1,1,1.35]*5,hspace=.13)
    raster_names=['Truth']+[m[1] for m in methods]
    raster_keys=['truth']+[m[0] for m in methods]
    raster_colors=['#333333']+[m[2] for m in methods]
    top_row=len(raster_keys)-1
    example_counts={key:int(np.count_nonzero((v>=0)&(v<.5))) for key,v in events.items()}
    for block in range(5):
        lo,hi=block*.1,(block+1)*.1
        start,stop=round(lo*fs),round(hi*fs)
        t=np.arange(start,stop)/fs*1000
        axes=[fig.add_subplot(gs[3*block+j]) for j in range(3)]
        for channel,ax in enumerate(axes[:2]):
            signal=z.real if channel==0 else z.imag
            ax.plot(t,signal[start:stop],color='#a6aeb9',lw=.5,label='Observed IQ')
            smooth=np.convolve(signal,np.ones(3)/3,mode='same')
            ax.plot(t,smooth[start:stop],color='#355f8d',lw=.75,label='3-sample mean (display only)')
            for event in events['truth'][(events['truth']>=lo)&(events['truth']<hi)]:
                ax.axvline(event*1000,color='#333333',ls='--',lw=.7,alpha=.65)
            ax.set_ylabel(('I' if channel==0 else 'Q')+' / σ̂')
            ax.tick_params(labelbottom=False)
        axes[0].set_title(f'{lo*1000:.0f}–{hi*1000:.0f} ms',loc='left',fontsize=12)
        if block==0:axes[0].legend(loc='upper right',fontsize=8,ncol=2)
        for y,(key,color) in enumerate(zip(raster_keys,raster_colors)):
            v=events[key];visible=v[(v>=lo)&(v<hi)]
            axes[2].vlines(visible*1000,top_row-y-.27,top_row-y+.27,color=color,lw=1.5)
            if len(visible)==0:axes[2].text((lo+hi)*500,top_row-y,'no flips',color=color,ha='center',va='center',fontsize=8)
        axes[2].set_yticks(list(range(top_row,-1,-1)),raster_names,fontsize=9);axes[2].set_ylim(-.55,top_row+.55)
        axes[2].set_xlabel('Time from record start (ms)')
        for ax in axes:
            for pulse in pulses[(pulses>=lo)&(pulses<hi)]:
                ax.axvline(pulse*1000,color='#bd8a16',lw=.8,alpha=.7)
                ax.axvspan(pulse*1000,(pulse+plan['nominal']['on_time_s'])*1000,color='#e2b140',alpha=.25)
            ax.set_xlim(lo*1000,hi*1000);ax.grid(axis='x',alpha=.15)
    fig.suptitle('Nominal 0.5-second I/Q waveform and flip reconstruction\n'
        f'Max SNR = {nominal_snr} · mean LED flips/pulse = {nominal_mean:g} · 10 kHz · EMG μ/σ/τ = 0.2/0.4/0.8 ms\n'
        'Black dashed lines: true flips · gold lines: LED onset · same waveform for every method',fontsize=14)
    fig.supxlabel('Fixed first realization (seed 106111000); 0–0.5 s of a 1 s record; no example selection after prediction.\n'
        +'Window counts: '+', '.join(f'{name} {example_counts[key]}' for key,name in zip(raster_keys,raster_names))+
        ('\n'+decoder_note),fontsize=10)
    for ext in ['png','pdf']:fig.savefig(a.out/f'waveform_nominal_0p5s.{ext}',dpi=170)
    plt.close(fig)
    print(json.dumps(dict(figures=str(a.out),example_counts=example_counts),indent=2))


if __name__=='__main__':main()
