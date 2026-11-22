"""Independent assignment oracle and aggregation regression checks for C1."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from qpd.mlebench.led_generalization_score import score_record, timing_settings, grade_generalization


def oracle(t,p,fs):
    eps,cut = timing_settings(fs)
    n,m = len(t),len(p)
    if not n or not m:
        return n+3*m
    d = np.abs(np.asarray(t)[:,None]-np.asarray(p)[None,:])
    savings = np.where(d<cut,4-4*(np.maximum(d-eps,0)/(cut-eps))**2,0)
    i,j = linear_sum_assignment(savings,maximize=True)
    return n+3*m-savings[i,j].sum()


def main():
    rng=np.random.default_rng(801)
    for fs in (5000,10000,40000,500000):
        for _ in range(150):
            # Quantized draws deliberately include multiplicities and ties.
            t=np.sort(rng.integers(0,80,rng.integers(0,9)))*.00005
            p=np.sort(rng.integers(0,80,rng.integers(0,9)))*.00005
            r=score_record(t,p,fs)
            assert np.isclose(r['total_cost'],oracle(t,p,fs),atol=1e-10),(t,p,r)
            assert r['n_matches']+r['false_negatives']==len(t)
            assert r['n_matches']+r['false_positives']==len(p)
    for fs in (5000,10000,40000,500000):
        eps,cut=timing_settings(fs)
        assert score_record([0],[eps],fs)['total_cost']==0
        assert np.isclose(score_record([0],[cut],fs)['total_cost'],4)
        assert score_record([0],[cut],fs)['n_matches']==0
        assert score_record([0],[np.nextafter(cut,np.inf)],fs)['total_cost']==4
    assert score_record([.001]*4,[.001]*3,10000)['total_cost']==1
    assert score_record([],[],10000)['total_cost']==0
    assert score_record([],[.001,.002],10000)['total_cost']==6
    for bad in ([np.nan],[np.inf],[.1,.01]):
        try: score_record(bad,[],10000)
        except ValueError: pass
        else: raise AssertionError('Invalid times accepted')
    rows=[]; truths=[]; preds=[]
    for component in ('broad','nominal','stress'):
        for zero in (True,False):
            rid=f'{component}{zero}'
            rows.append(dict(record_id=rid,component=component,zero_led=zero,duration_s=1.,sample_rate_hz=10000.))
            truths.append(dict(record_id=rid,t=.1))
            preds.append(dict(record_id=rid,t=.1))
    empty=dict(record_id='empty',component='broad',zero_led=True,duration_s=1.,sample_rate_hz=10000.)
    rows.append(empty); preds.append(dict(record_id='empty',t=.2))
    manifest=pd.DataFrame(rows); answers=pd.DataFrame(truths); submission=pd.DataFrame(preds)
    result=grade_generalization(submission,answers,manifest)
    assert result['primary_valid'] and np.isclose(result['error'],.6*.1*3)
    assert result['empty_truth']['false_positives_per_second']==1
    copies=[]
    for table in (submission,answers,manifest):
        duplicate=table.copy(); duplicate.record_id=duplicate.record_id+'copy'
        copies.append(pd.concat([table,duplicate],ignore_index=True))
    assert np.isclose(grade_generalization(*copies)['error'],result['error'])
    partial=grade_generalization(submission[submission.record_id.str.startswith('broad')],answers[answers.record_id.str.startswith('broad')],manifest[manifest.component=='broad'])
    assert not partial['primary_valid'] and partial['error'] is None
    shuffled=grade_generalization(submission.sample(frac=1,random_state=1),answers.sample(frac=1,random_state=2),manifest)
    assert shuffled['error']==result['error']
    # Many disconnected pairs: ensures no dense global matrix allocation.
    t=np.arange(20000)*.01
    assert score_record(t,t,10000)['total_cost']==0
    print('PASS: 600 independent Hungarian-oracle cases; boundary/multiplicity/empty/aggregation/scaling checks')


if __name__=='__main__': main()
