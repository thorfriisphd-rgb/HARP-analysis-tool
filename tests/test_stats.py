import numpy as np

from harp.stats import _plus_one_p, centered_phase_signature, panel_leave_one_out_test


def test_plus_one_p_minimum():
    null = np.array([0.0, 0.1, 0.2, 0.3])
    assert _plus_one_p(1.0, null) == 1 / 5


def test_centered_signature_is_centered_and_unit_length():
    signal = np.arange(1.0, 8.0)
    phase_index = np.arange(7, dtype=float)
    sig = centered_phase_signature(signal, phase_index)
    assert np.isclose(sig.mean(), 0.0, atol=1e-12)
    assert np.isclose(np.linalg.norm(sig), 1.0, atol=1e-12)


def test_panel_test_is_seed_deterministic():
    rng = np.random.default_rng(4)
    x = rng.normal(size=(5, 7))
    x = x - x.mean(axis=1, keepdims=True)
    x = x / np.linalg.norm(x, axis=1, keepdims=True)
    a, scores_a, null_a = panel_leave_one_out_test(x, n_permutations=99, seed=17)
    b, scores_b, null_b = panel_leave_one_out_test(x, n_permutations=99, seed=17)
    assert a == b
    assert np.array_equal(scores_a, scores_b)
    assert np.array_equal(null_a, null_b)
