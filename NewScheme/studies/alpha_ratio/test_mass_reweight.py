import unittest
import numpy as np
from pt_reweight_core import Estimate
from test_pt_reweight import fixture
from mass_reweight_dy import mass_grid, fit_mass, expand_fit, mass_predictions, calibration_predictions


class TestMassReweight(unittest.TestCase):
    def setUp(self):
        self.space = fixture()
        self.groups = [np.array([0]), np.array([1])]

    def fits(self, regions):
        return {ch: fit_mass(self.space, 1, ch, regions, self.groups) for ch in (0, 1)}

    def test_mass_calibration_values_and_jacobian(self):
        result = calibration_predictions(self.space, 1, [0, 1], [0, 2], self.groups, self.fits([0, 2]))
        np.testing.assert_allclose(result['dy_reweighted'].value, result['observed'].value)
        np.testing.assert_allclose(result['dy_reweighted'].jac, result['observed'].jac, atol=1e-12)

    def test_alpha_cancels_nonconstant_mass_factors(self):
        fits = self.fits([0, 2])
        for ch, fit in fits.items():
            fit['factor'] = Estimate(np.array([.7+ch*.1, 1.6+ch*.2]), fit['factor'].jac*1.4)
        result = mass_predictions(self.space, 1, [0, 1], [1], [0, 2], self.groups, fits)
        np.testing.assert_allclose(result['alpha_reweighted'].value, result['alpha_nominal'].value)
        np.testing.assert_allclose(result['alpha_reweighted'].jac, result['alpha_nominal'].jac, atol=1e-12)

    def test_direct_equals_alpha_with_matched_calibration(self):
        for target, source in [([1], [0, 2]), ([2], [0])]:
            result = mass_predictions(self.space, 1, [0, 1], target, source, self.groups, self.fits(source))
            np.testing.assert_allclose(result['dy_reweighted'].value, result['alpha_nominal'].value)
            np.testing.assert_allclose(result['dy_reweighted'].jac, result['alpha_nominal'].jac, atol=1e-12)

    def test_channel_isolation_and_no_sr_fit(self):
        for ch, fit in self.fits([0]).items():
            for proc in (0, 1):
                self.assertTrue(np.all(fit['factor'].jac[:, self.space.indices(proc, [1], [1-ch], [0, 1, 2])] == 0))
            self.assertTrue(np.all(fit['factor'].jac[:, self.space.indices(0, [1], [ch], [1, 2])] == 0))

    def test_joint_finite_difference(self):
        def evaluate():
            return mass_predictions(self.space, 1, [0, 1], [1], [0, 2], self.groups, self.fits([0, 2]))
        central = evaluate()
        for index in np.random.default_rng(4).choice(self.space.flat.size, 24, replace=False):
            old, eps = self.space.flat[index], 1e-3
            self.space.flat[index] = old+eps
            up = evaluate()
            self.space.flat[index] = old-eps
            down = evaluate()
            self.space.flat[index] = old
            for name in ('dy_reweighted', 'alpha_reweighted', 'dy_normonly'):
                np.testing.assert_allclose((up[name].value-down[name].value)/(2*eps), central[name].jac[:, index], rtol=2e-6, atol=2e-8)

    def test_grid_does_not_use_sr_observation(self):
        before = mass_grid(self.space, 1, self.groups)
        self.space.flat[self.space.indices(0, [1], [0, 1], [1])] = 1e9
        after = mass_grid(self.space, 1, self.groups)
        self.assertEqual([g.tolist() for g in before], [g.tolist() for g in after])

    def test_alpha_overlap_rejected(self):
        with self.assertRaises(ValueError):
            mass_predictions(self.space, 1, [0], [2], [0, 2], self.groups, self.fits([0, 2]))

    def test_coarse_correction_on_original_alpha_grid(self):
        self.space.flat[self.space.indices(0, [1], [0], [0, 2], [0])] *= 1.5
        coarse_groups = [np.array([0, 1])]
        coarse = {ch: fit_mass(self.space, 1, ch, [0, 2], coarse_groups) for ch in (0, 1)}
        expanded = {ch: expand_fit(fit, coarse_groups, self.groups) for ch, fit in coarse.items()}
        fine = mass_predictions(self.space, 1, [0], [1], [0, 2], self.groups, expanded)
        np.testing.assert_allclose(fine['alpha_reweighted'].value, fine['alpha_nominal'].value)
        np.testing.assert_allclose(fine['alpha_reweighted'].jac, fine['alpha_nominal'].jac, atol=1e-12)
        self.assertFalse(np.allclose(fine['dy_reweighted'].value, fine['alpha_nominal'].value))
        summed = mass_predictions(self.space, 1, [0], [1], [0, 2], coarse_groups, coarse)
        np.testing.assert_allclose(summed['dy_reweighted'].value, summed['alpha_nominal'].value)
        calibration = calibration_predictions(self.space, 1, [0, 1], [0, 2], self.groups, expanded)
        np.testing.assert_allclose(calibration['dy_reweighted'].total().value, calibration['observed'].total().value)
        np.testing.assert_allclose(calibration['dy_reweighted'].total().jac, calibration['observed'].total().jac, atol=1e-12)


if __name__ == '__main__':
    unittest.main()
