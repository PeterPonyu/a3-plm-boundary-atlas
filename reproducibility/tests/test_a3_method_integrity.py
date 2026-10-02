"""Regression checks for original-identity isolation, without fitting new paper estimates."""
from pathlib import Path
import random
import sys
import pytest

CAP = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(CAP/'analysis/xfit_calibration_20260918'))
import run_xfit_calibration as xfit
import run_xfit_round3 as r3


@pytest.fixture(scope='module')
def rows():
    return xfit.load_manifest(xfit.INPUTS/'manifest_gof_650M.csv')


def test_original_unique_row_partition_is_unchanged(rows):
    assert len(rows) == len({r['example_id'] for r in rows}) == 4677
    assert xfit.assign_identity_folds(rows, 5, xfit.SEED) == xfit.assign_folds(
        [r['y'] for r in rows], 5, xfit.SEED)


def test_all_200_original_draws_expose_old_leakage_and_new_folds_isolate(rows):
    rng = random.Random(r3.NESTED_SEED)
    for b in range(r3.N_NESTED):
        sample = [rows[rng.randrange(len(rows))] for _ in rows]
        old = xfit.assign_folds([r['y'] for r in sample], 5, r3.NESTED_SEED+b)
        with pytest.raises(ValueError, match='crosses'):
            xfit.assert_identity_separation(sample, old)
        fixed = xfit.assign_identity_folds(sample, 5, r3.NESTED_SEED+b)
        xfit.assert_identity_separation(sample, fixed)
        for fold in range(5):
            train = {r['example_id'] for r, f in zip(sample, fixed) if f != fold}
            test = {r['example_id'] for r, f in zip(sample, fixed) if f == fold}
            assert train and test and train.isdisjoint(test)


@pytest.mark.parametrize('change', [{'s': 999.0}, {'y': 3}, {'gene': 'conflicting-gene'}])
def test_conflicting_identity_is_rejected(rows, change):
    with pytest.raises(ValueError, match='conflicting'):
        xfit.assign_identity_folds(rows+[dict(rows[0], **change)], 5, 1)


def test_missing_identity_and_too_few_units_rejected(rows):
    with pytest.raises(ValueError, match='example_id'):
        xfit.assign_identity_folds([dict(rows[0], example_id='')]+rows[1:], 5, 1)
    with pytest.raises(ValueError, match='fewer'):
        xfit.assign_identity_folds(rows[:2]*10, 5, 1)


def test_nested_fit_refuses_silent_single_class_fallback():
    sample = [dict(example_id=str(i), gene='toy', s=float(i), y=1) for i in range(10)]
    with pytest.raises(ValueError, match='one class'):
        r3.nested_crossfit(sample, 1)


def test_corrected_nested_output_cannot_replace_historical_ledger():
    assert r3.OUT.name == 'round3_identity_grouped_candidate'
    assert r3.OUT != r3.HERE/'outputs/round3'
