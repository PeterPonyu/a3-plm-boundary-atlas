import importlib.util
from pathlib import Path
import numpy as np
import pytest

spec = importlib.util.spec_from_file_location("clinical", Path(__file__).with_name("analyze_clinical.py"))
c = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c)


def test_isotonic_equal_scores_weighted_and_permutation_invariant():
    x=np.array([.3,.3,.1,.8,.8]); y=np.array([0,1,1,0,1]); w=np.array([1,3,1,2,1])
    expected=c.isotonic_fit(x,y,w)(x)
    for seed in range(20):
        order=np.random.default_rng(seed).permutation(len(x))
        assert np.array_equal(c.isotonic_fit(x[order],y[order],w[order])(x),expected)
    assert expected[0] == expected[1]


def test_perfect_discrimination_percentile_interval_degenerates():
    y=np.array([0,0,1,1]); x=np.array([-2.,-1.,1.,2.])
    assert c.auc_interval(y,x,np.random.default_rng(5),100) == [1.,1.]


def test_weighted_metrics_match_expanded_copies():
    p=np.array([.2,.5,.9]); y=np.array([0,1,1]); w=np.array([3,1,2])
    a=c.metrics(p,y,w); b=c.metrics(np.repeat(p,w),np.repeat(y,w))
    assert a == pytest.approx(b)


def test_refit_uses_original_identity_weights_and_rejects_duplicate_rows():
    rows=[dict(example_id=f"v{i:02}",gene="G",primary_llr_fp32=str((i%7)-3),label_benign=str(i%2)) for i in range(40)]
    p,y,fold=c.crossfit(rows,np.array([3 if i%3==0 else 1 for i in range(40)]))
    assert set(fold)==set(range(5))
    assert all(np.all(np.isfinite(values)) for values in p.values())
    with pytest.raises(ValueError,match="One row per original"):
        c.crossfit(rows+[rows[0]])


def test_prevalence_baseline_low_ece_does_not_imply_discrimination():
    y=np.array([0,0,0,1]); p=np.repeat(.25,4)
    assert c.metrics(p,y)["ece"] == 0.
    assert c.metrics(p,y)["brier"] == .1875


def test_pav_fitted_values_match_independent_sklearn_implementation():
    from sklearn.isotonic import IsotonicRegression
    for seed in range(30):
        rng = np.random.default_rng(seed)
        x = rng.integers(0, 12, 60).astype(float)
        y = rng.integers(0, 2, 60)
        w = rng.integers(1, 4, 60)
        expected = IsotonicRegression(out_of_bounds="clip").fit(x, y, sample_weight=w)
        assert c.isotonic_fit(x,y,w)(np.unique(x)) == pytest.approx(
            np.clip(expected.predict(np.unique(x)),1e-6,1-1e-6), abs=1e-12)


def test_crossfit_fit_only_standardization_and_identity_multiplicity(monkeypatch):
    rows=[dict(example_id=f"v{i:02}",gene="G",primary_llr_fp32=str((i%7)-3),
               label_benign=str(i%2)) for i in range(40)]
    weights=np.array([0 if i==0 else 3 if i%3==0 else 1 for i in range(40)])
    fitted=[]
    original=c.platt_fit
    def record(z,y,w):
        fitted.append((z.copy(),y.copy(),w.copy()))
        return original(z,y,w)
    monkeypatch.setattr(c,"platt_fit",record)
    predictions,y,fold=c.crossfit(rows,weights)
    score=-np.array([float(r["primary_llr_fp32"]) for r in rows])
    assert fold[0] == -1
    for f,(z,yt,wt) in enumerate(fitted):
        train=(fold>=0)&(fold!=f)
        test=fold==f
        assert np.array_equal(wt,weights[train])
        assert np.array_equal(yt,y[train])
        mu=np.average(score[train],weights=weights[train])
        sd=np.sqrt(np.average((score[train]-mu)**2,weights=weights[train]))
        assert z == pytest.approx((score[train]-mu)/sd)
        assert predictions["sigmoid_z"][test] == pytest.approx(c.expit((score[test]-mu)/sd))
