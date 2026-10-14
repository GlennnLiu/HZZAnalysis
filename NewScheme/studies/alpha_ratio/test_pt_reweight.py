import unittest
import numpy as np
from pt_reweight_core import Space, fit_pt, predictions, Estimate


def fixture():
    arrays={'mass_edges__resolved':np.array([300,500,np.inf]),
            'pt_edges':np.array([0,100,200,np.inf])}
    for tag in ('btag','untag','vbftag'):
        shape=(2,3,2,3)
        dy=np.ones(shape)*50
        dy[:,1,:,0]=80;dy[:,2,:,2]=60
        for proc in ('DY','Data','TTplusWW','VZ'):
            if proc=='DY': value=dy;var=dy*1.3
            elif proc=='Data':
                value=dy*np.array([.8,1.1,1.4])+5;var=value
            else: value=np.ones(shape)*2.5;var=value*.4
            arrays['resolved__'+tag+'__'+proc]=np.stack([value,var],axis=-1)
    return Space(arrays,'resolved')


class TestPTReweight(unittest.TestCase):
    def test_adaptive_ratio_and_cr_closure(self):
        space=fixture();fit=fit_pt(space,[0,1,2],[0,2])
        weighted=space.weighted([0,1,2],[0,1],[0,2],
                    [np.arange(2)],fit['groups'],fit['factor'])
        observed=space.linear(0,[0,1,2],[0,1],[0,2])
        np.testing.assert_allclose(weighted.value,observed.value)
        np.testing.assert_allclose(weighted.jac,observed.jac,atol=1e-12)

    def test_joint_gradients_finite_difference(self):
        space=fixture();groups=[np.array([i]) for i in range(3)]
        def evaluate():
            z=space.linear(0,[0,1,2],[0,1],[0,2],ptgroups=groups)
            d=space.linear(1,[0,1,2],[0,1],[0,2],ptgroups=groups)
            fit={'groups':groups,'factor':z/d,'observation':z,'dy':d}
            return predictions(space,1,[0,1],2,[np.array([0]),np.array([1])],fit)
        result=evaluate()
        indices=np.random.default_rng(1).choice(space.flat.size,25,replace=False)
        for index in indices:
            old=space.flat[index];epsilon=1e-3
            space.flat[index]=old+epsilon;up=evaluate()
            space.flat[index]=old-epsilon;down=evaluate()
            space.flat[index]=old
            for key in ('dy_reweighted','alpha_reweighted','dy_normonly'):
                numeric=(up[key].value-down[key].value)/(2*epsilon)
                np.testing.assert_allclose(result[key].jac[:,index],numeric,atol=2e-8,rtol=2e-6)

    def test_constant_sf_cancels_alpha(self):
        space=fixture();fit=fit_pt(space,[1],[0,2])
        fit['factor']=Estimate(np.ones(len(fit['groups']))*1.7,
                               np.zeros_like(fit['factor'].jac))
        r=predictions(space,1,[0,1],1,[np.array([0]),np.array([1])],fit)
        np.testing.assert_allclose(r['alpha_reweighted'].value,r['alpha_nominal'].value)
        np.testing.assert_allclose(r['alpha_reweighted'].jac,r['alpha_nominal'].jac,atol=1e-12)

    def test_no_sr_data_in_fit_or_prediction(self):
        space=fixture();fit=fit_pt(space,[1],[0,2])
        r=predictions(space,1,[0,1],1,[np.array([0]),np.array([1])],fit)
        indices=space.indices(0,[0,1,2],[0,1],[1])
        self.assertTrue(np.all(fit['factor'].jac[:,indices]==0))
        self.assertTrue(np.all(r['alpha_reweighted'].jac[:,indices]==0))

    def test_shared_fit_channel_covariance(self):
        space=fixture();fit=fit_pt(space,[1],[0,2])
        groups=[np.array([0]),np.array([1])]
        a=predictions(space,1,[0],1,groups,fit)['dy_reweighted']
        b=predictions(space,1,[1],1,groups,fit)['dy_reweighted']
        joint=(a+b).covariance(space.variance)
        independent=a.covariance(space.variance)+b.covariance(space.variance)
        self.assertGreater(float(joint.sum()),float(independent.sum()))

    def test_inclusive_source_gradient(self):
        space=fixture();fit=fit_pt(space,[1],[0,2]);groups=fit['groups']
        def evaluate():
            z=space.linear(0,[1],[0,1],[0,2],ptgroups=groups)
            d=space.linear(1,[1],[0,1],[0,2],ptgroups=groups)
            current={'groups':groups,'factor':z/d,'observation':z,'dy':d}
            return predictions(space,1,[0,1],1,[np.array([0]),np.array([1])],current,cr_regions=(0,2))['alpha_reweighted']
        central=evaluate()
        for index in space.indices(0,[1],[0],[2],[0],[0]):
            value=space.flat[index];epsilon=1e-3
            space.flat[index]=value+epsilon;up=evaluate().value
            space.flat[index]=value-epsilon;down=evaluate().value
            space.flat[index]=value
            np.testing.assert_allclose((up-down)/(2*epsilon),central.jac[:,index],rtol=2e-6,atol=2e-8)

    def test_overlapping_alpha_source_is_rejected(self):
        space=fixture();fit=fit_pt(space,[1],[0,2])
        with self.assertRaises(ValueError):
            predictions(space,1,[0,1],2,[np.array([0]),np.array([1])],fit,cr_regions=(0,2))


if __name__=='__main__':unittest.main()
