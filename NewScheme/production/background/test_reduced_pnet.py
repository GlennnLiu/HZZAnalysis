import unittest
import awkward as ak
import numpy as np
import make_reduced_trees as p


class Weights(unittest.TestCase):
    def setUp(self):
        self.e=ak.Array(dict(EventWeight=[2.,-1.],ptmerged=[250.,650.],
            particleNetZvsQCD=[.92,.97],particleNetZbbvslight=[.2,.95],
            isbjet=[True,False],iscjet=[False,False],islightjet=[False,True],w_deepjet=[.9,1.1]))
        self.c=dict(lumi={'2018':1.},samples_inf={'DY_test':['MC',1.],'TTTo2L2Nu':['MC',1.]})
        self.s={tag:{flav:{pur:{pt:(.8 if tag=='ZvsQCD' else 1.2)
            for pt in ['200To400','400To600','600ToInf']} for pur in ['LP','MP','HP']}
            for flav in ['bhadron','chadron','lighthadron']} for tag in ['ZvsQCD','ZvsLight']}

    def test_dy_full_and_no_pnet(self):
        w=p.mc_weights(self.e,self.c,'2018','DY_test','DY','merged','untag',1000,self.s,True)
        np.testing.assert_allclose(w['weight_no_pnet'],[2.,-1.])
        np.testing.assert_allclose(w['weight'],np.array([2.,-1.])*.8*1.2)

    def test_vbf_only_qcd(self):
        w=p.mc_weights(self.e,self.c,'2018','DY_test','DY','merged','vbftag',1000,self.s,True)
        np.testing.assert_allclose(w['pnet_sf_zvslight'],[1.,1.])
        np.testing.assert_allclose(w['weight'],[1.6,-.8])

    def test_resolved_keeps_deepjet(self):
        w=p.mc_weights(self.e,self.c,'2018','TTTo2L2Nu','TTplusWW','resolved','untag',1000,self.s,True)
        np.testing.assert_allclose(w['weight_no_pnet'],[1.8,-1.1])
        np.testing.assert_array_equal(w['weight'],w['weight_no_pnet'])

    def test_uncovered_raises(self):
        e=ak.with_field(self.e,ak.Array([150.,650.]),'ptmerged')
        with self.assertRaises(ValueError): p.particle_net_sf('2018',e,'ZvsQCD',self.s)

    def test_ambiguous_flavour_raises(self):
        e=ak.with_field(self.e,ak.Array([True,True]),'islightjet')
        with self.assertRaises(ValueError): p.particle_net_sf('2018',e,'ZvsQCD',self.s)


if __name__=='__main__': unittest.main()
