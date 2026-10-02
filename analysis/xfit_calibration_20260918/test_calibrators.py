"""Isotonic tie handling on clinical scores must be row-order invariant."""
import random

import pytest

from calibrators import fit_isotonic


def test_two_equal_scores_share_weighted_mean_in_either_order():
    for labels in ([0, 1], [1, 0]):
        predictor = fit_isotonic([.3, .3], labels)
        assert predictor([.3, .3, .2, .4]) == [.5, .5, .5, .5]


def test_ties_cannot_be_split_by_pav():
    scores = [.1, .3, .3, .3, .6, .9]
    labels = [0, 1, 0, 1, 0, 1]
    expected = fit_isotonic(scores, labels)(scores)
    rng = random.Random(9)
    for _ in range(100):
        order = list(range(len(scores)))
        rng.shuffle(order)
        actual = fit_isotonic([scores[i] for i in order], [labels[i] for i in order])(scores)
        assert actual == expected
    assert expected[1] == expected[2] == expected[3]


@pytest.mark.parametrize("scores,labels", [([.1, .2], [1]), ([], []), ([float("nan")], [1]), ([.1], [2])])
def test_bad_training_input_is_rejected(scores, labels):
    with pytest.raises(ValueError):
        fit_isotonic(scores, labels)
