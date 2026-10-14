import numpy as np
import os 
import correctionlib
import yaml
import uproot
import hist
import awkward as ak

class DeepJet:
    def __init__(self,year,array,config = None,WP = 'M'):
        self.year = year
        self.array = array
        self.WP = WP
        self.file_sf_path = f"/cms/user/liugeliang/CMS/HZZAnalysis/deepjetSFsfiles/{year}/btagging.json.gz"
        self.file_sf = correctionlib.CorrectionSet.from_file(self.file_sf_path)

        self.bcjet_sf_name = 'deepJet_comb'
        self.lightjet_sf_name = 'deepJet_incl'

        year_btagWP_map = {
           '2016preAPV': 0.2598,
           '2016postAPV': 0.2489,
           '2017': 0.3040,
           '2018': 0.2783
        }
        self.btagWP = year_btagWP_map.get(self.year, 0.0)



        if config == None:
            with open(f"/cms/user/liugeliang/CMS/HZZAnalysis/cards/config_UL{self.year}.yml") as f:
                self.config = yaml.safe_load(f)
        else:
            self.config = config

    def get_deepjet_sf(self,type):
        r"""
            Get the deepjet SF for each jet in the array.
            SF is separated by jet flavor, pt, abs(eta), and WP.

            input:
            type: related to the uncertainty type, such as "central", "up", "down".

            output:
            sf_jet1: the deepjet SF for jet 1 in the array.
            sf_jet2: the deepjet SF for jet 2 in the array.
        """
        # Get the SF for each jet in the array for jet 1 and jet 2.
        #jet 1
        #bcjet
        sf_jet1_bcjet = self.file_sf[self.bcjet_sf_name].evaluate(
            f"{type}",
            f"{self.WP}",
            self.array.jet_1_hadronflavor[self.array.jet_1_hadronflavor != 0],
            abs(self.array.jet_1_eta[self.array.jet_1_hadronflavor != 0]),
            self.array.jet_1_pt[self.array.jet_1_hadronflavor != 0],
            )
        #lightjet
        sf_jet1_lightjet = self.file_sf[self.lightjet_sf_name].evaluate(
            f"{type}",
            f"{self.WP}",
            self.array.jet_1_hadronflavor[self.array.jet_1_hadronflavor == 0],
            abs(self.array.jet_1_eta[self.array.jet_1_hadronflavor == 0]),
            self.array.jet_1_pt[self.array.jet_1_hadronflavor == 0],
            )
        #combine
        sf_jet1 = np.append(sf_jet1_bcjet,sf_jet1_lightjet)

        #jet 2
        #bcjet
        sf_jet2_bcjet = self.file_sf[self.bcjet_sf_name].evaluate(
            f"{type}",
            f"{self.WP}",
            self.array.jet_2_hadronflavor[self.array.jet_2_hadronflavor != 0],
            abs(self.array.jet_2_eta[self.array.jet_2_hadronflavor != 0]),
            self.array.jet_2_pt[self.array.jet_2_hadronflavor != 0],
            )
        #lightjet
        sf_jet2_lightjet = self.file_sf[self.lightjet_sf_name].evaluate(
            f"{type}",
            f"{self.WP}",
            self.array.jet_2_hadronflavor[self.array.jet_2_hadronflavor == 0],
            abs(self.array.jet_2_eta[self.array.jet_2_hadronflavor == 0]),
            self.array.jet_2_pt[self.array.jet_2_hadronflavor == 0],
            )
        #combine
        sf_jet2 = np.append(sf_jet2_bcjet,sf_jet2_lightjet)

        return sf_jet1,sf_jet2
    
    #make cut
    def make_cut(self,selection,arr):
        r'''
        input: selection, arr as ak.Aarry()
        output: cuted arr with selection
        '''
        cut  = ak.numexpr.evaluate(selection,arr)
        return arr[cut]
    
    def get_untag_deepjet_weight(self,sample):

        # Determine process based on sample
        def determine_process(sample):
            if sample in ['TTTo2L2Nu', 'tZq']:
                return 'top'
            elif sample in ['ZZTo2Q2L', 'WZTo2Q2L','WWTo2L2Nu']:
                return 'zv'
            return None
        
        process = determine_process(sample)
        if not process:
            raise ValueError("Sample not recognized")
        
        #open the eff file
        f_eff = uproot.open(f"/cms/user/liugeliang/CMS/HZZAnalysis/deepjetSFsfiles/{self.year}/btag_eff_{process}.root")
        #get 2D eff hist
        heff_b = f_eff['eff_b'].to_hist()
        heff_c = f_eff['eff_c'].to_hist()
        heff_udsg = f_eff['eff_udsg'].to_hist()

        #add four columns to the array. one is jet_1_sfs and the other is jet_2_sfs and the last one is w_btag, and eff_btag
        sfs1,sfs2 = self.get_deepjet_sf("central")
        self.array['jet_1_sfs'] = sfs1
        self.array['jet_2_sfs'] = sfs2
        #self.array['w_btag'] = ak.ones_like(self.array.jet_1_pt)
        self.array['eff_btag'] = ak.ones_like(self.array.jet_1_pt)
    

        #start to produce weight
        #case 1: jet 1 pass WP, jet 2 fail WP, and jet 2 is bjet
        selection = f"(jet_1_deepbtag >= {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_2_hadronflavor == 5)"
        arr1 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr1.jet_2_pt,abs(arr1.jet_2_eta))])
        w_btag = (arr1.jet_1_sfs*(1-arr1.jet_2_sfs*eff_fail))/(1-eff_fail)
        #add a new column to the array
        arr1['w_btag'] = w_btag

        #case 2: jet 1 pass WP, jet 2 fail WP, and jet 2 is cjet
        selection = f"(jet_1_deepbtag >= {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_2_hadronflavor == 4)"
        arr2 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr2.jet_2_pt,abs(arr2.jet_2_eta))])
        w_btag = (arr2.jet_1_sfs*(1-arr2.jet_2_sfs*eff_fail))/(1-eff_fail)
        #add a new column to the array
        arr2['w_btag'] = w_btag

        #case 3: jet 1 pass WP, jet 2 fail WP, and jet 2 is lightjet
        selection = f"(jet_1_deepbtag >= {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_2_hadronflavor == 0)"
        arr3 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr3.jet_2_pt,abs(arr3.jet_2_eta))])
        w_btag = (arr3.jet_1_sfs*(1-arr3.jet_2_sfs*eff_fail))/(1-eff_fail)
        #add a new column to the array
        arr3['w_btag'] = w_btag

        #case 4: jet 1 fail WP, jet 2 pass WP, and jet 1 is bjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag >= {self.btagWP}) & (jet_1_hadronflavor == 5)"
        arr4 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr4.jet_1_pt,abs(arr4.jet_1_eta))])
        w_btag = (arr4.jet_2_sfs*(1-arr4.jet_1_sfs*eff_fail))/(1-eff_fail)
        #add a new column to the array
        arr4['w_btag'] = w_btag

        #case 5: jet 1 fail WP, jet 2 pass WP, and jet 1 is cjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag >= {self.btagWP}) & (jet_1_hadronflavor == 4)"
        arr5 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr5.jet_1_pt,abs(arr5.jet_1_eta))])
        w_btag = (arr5.jet_2_sfs*(1-arr5.jet_1_sfs*eff_fail))/(1-eff_fail)
        #add a new column to the array
        arr5['w_btag'] = w_btag

        #case 6: jet 1 fail WP, jet 2 pass WP, and jet 1 is lightjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag >= {self.btagWP}) & (jet_1_hadronflavor == 0)"
        arr6 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr6.jet_1_pt,abs(arr6.jet_1_eta))])
        w_btag = (arr6.jet_2_sfs*(1-arr6.jet_1_sfs*eff_fail))/(1-eff_fail)
        #add a new column to the array
        arr6['w_btag'] = w_btag

        #case 7: both jet 1 and jet 2 fail WP, and jet 1 is bjet and jet 2 is bjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 5) & (jet_2_hadronflavor == 5)"
        arr7 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr7.jet_1_pt,abs(arr7.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr7.jet_2_pt,abs(arr7.jet_2_eta))])
        w_btag = (1-arr7.jet_1_sfs*eff_fail1)*(1-arr7.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr7['w_btag'] = w_btag

        #case 8: both jet 1 and jet 2 fail WP, and jet 1 is bjet and jet 2 is cjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 5) & (jet_2_hadronflavor == 4)"
        arr8 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr8.jet_1_pt,abs(arr8.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr8.jet_2_pt,abs(arr8.jet_2_eta))])
        w_btag = (1-arr8.jet_1_sfs*eff_fail1)*(1-arr8.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr8['w_btag'] = w_btag

        #case 9: both jet 1 and jet 2 fail WP, and jet 1 is bjet and jet 2 is lightjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 5) & (jet_2_hadronflavor == 0)"
        arr9 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr9.jet_1_pt,abs(arr9.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr9.jet_2_pt,abs(arr9.jet_2_eta))])
        w_btag = (1-arr9.jet_1_sfs*eff_fail1)*(1-arr9.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr9['w_btag'] = w_btag

        #case 10: both jet 1 and jet 2 fail WP, and jet 1 is cjet and jet 2 is bjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 4) & (jet_2_hadronflavor == 5)"
        arr10 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr10.jet_1_pt,abs(arr10.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr10.jet_2_pt,abs(arr10.jet_2_eta))])
        w_btag = (1-arr10.jet_1_sfs*eff_fail1)*(1-arr10.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr10['w_btag'] = w_btag

        #case 11: both jet 1 and jet 2 fail WP, and jet 1 is cjet and jet 2 is cjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 4) & (jet_2_hadronflavor == 4)"
        arr11 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr11.jet_1_pt,abs(arr11.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr11.jet_2_pt,abs(arr11.jet_2_eta))])
        w_btag = (1-arr11.jet_1_sfs*eff_fail1)*(1-arr11.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr11['w_btag'] = w_btag

        #case 12: both jet 1 and jet 2 fail WP, and jet 1 is cjet and jet 2 is lightjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 4) & (jet_2_hadronflavor == 0)"
        arr12 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr12.jet_1_pt,abs(arr12.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr12.jet_2_pt,abs(arr12.jet_2_eta))])
        w_btag = (1-arr12.jet_1_sfs*eff_fail1)*(1-arr12.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr12['w_btag'] = w_btag

        #case 13: both jet 1 and jet 2 fail WP, and jet 1 is lightjet and jet 2 is bjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 0) & (jet_2_hadronflavor == 5)"
        arr13 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr13.jet_1_pt,abs(arr13.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_b[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr13.jet_2_pt,abs(arr13.jet_2_eta))])
        w_btag = (1-arr13.jet_1_sfs*eff_fail1)*(1-arr13.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr13['w_btag'] = w_btag

        #case 14: both jet 1 and jet 2 fail WP, and jet 1 is lightjet and jet 2 is cjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 0) & (jet_2_hadronflavor == 4)"
        arr14 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr14.jet_1_pt,abs(arr14.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_c[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr14.jet_2_pt,abs(arr14.jet_2_eta))])
        w_btag = (1-arr14.jet_1_sfs*eff_fail1)*(1-arr14.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr14['w_btag'] = w_btag

        #case 15: both jet 1 and jet 2 fail WP, and jet 1 is lightjet and jet 2 is lightjet
        selection = f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 0) & (jet_2_hadronflavor == 0)"
        arr15 = self.make_cut(selection,self.array)
        #get the eff and weight
        eff_fail1 = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr15.jet_1_pt,abs(arr15.jet_1_eta))])
        eff_fail2 = 1-ak.Array([heff_udsg[hist.loc(pt),hist.loc(eta)] for pt,eta in zip(arr15.jet_2_pt,abs(arr15.jet_2_eta))])
        w_btag = (1-arr15.jet_1_sfs*eff_fail1)*(1-arr15.jet_2_sfs*eff_fail2)/(1-eff_fail1)/(1-eff_fail2)
        #add a new column to the array
        arr15['w_btag'] = w_btag

        #combine all the arrays
        arr = ak.concatenate([arr1,arr2,arr3,arr4,arr5,arr6,arr7,arr8,arr9,arr10,arr11,arr12,arr13,arr14,arr15])

        return arr
    
    def get_untag_deepjet_array(self, sample):

        # Determine process based on sample
        def determine_process(sample):
            if sample in ['TTTo2L2Nu', 'tZq']:
                return 'top'
            elif sample in ['ZZTo2Q2L', 'WZTo2Q2L','WWTo2L2Nu']:
                return 'zv'
            return None

        # Open the efficiency file and get histograms
        def get_eff_histograms(year, process):
            with uproot.open(f"/cms/user/liugeliang/CMS/HZZAnalysis/deepjetSFsfiles/{year}/btag_eff_{process}.root") as f_eff:
                heff_b = f_eff['eff_b'].to_hist()
                heff_c = f_eff['eff_c'].to_hist()
                heff_udsg = f_eff['eff_udsg'].to_hist()
            return heff_b, heff_c, heff_udsg

        # Calculate weights based on the selection and histogram data
        def calculate_weights_pass_fail(array, selection, histogram, pass_sfs, fail_sfs,failpt, faileta):
            arr = self.make_cut(selection,array)
            pt, eta = arr[failpt], abs(arr[faileta])
            eff_fail = 1 - ak.Array([histogram[hist.loc(thispt), hist.loc(thiseta)] for thispt, thiseta in zip(pt, eta)])
            sfs_pass = arr[pass_sfs]
            sfs_fail = arr[fail_sfs]
            w_btag = (sfs_pass*(1 - sfs_fail*eff_fail))/(1 - eff_fail)
            arr['w_btag'] = w_btag
            return arr

        def calculate_weights_both_fail(array,selection,h1,h2):
            arr = self.make_cut(selection,array)
            pt1, eta1 = arr.jet_1_pt, abs(arr.jet_1_eta)
            pt2, eta2 = arr.jet_2_pt, abs(arr.jet_2_eta)
            eff_fail1 = 1-ak.Array([h1[hist.loc(thispt),hist.loc(thiseta)] for thispt,thiseta in zip(pt1,eta1)])
            eff_fail2 = 1-ak.Array([h2[hist.loc(thispt),hist.loc(thiseta)] for thispt,thiseta in zip(pt2,eta2)])
            sfs1 = arr.jet_1_sfs
            sfs2 = arr.jet_2_sfs
            w_btag = ((1 - sfs1 * eff_fail1) * (1 - sfs2 * eff_fail2)) / ((1 - eff_fail1) * (1 - eff_fail2))
            arr['w_btag'] = w_btag
            return arr

        process = determine_process(sample)
        if not process:
            raise ValueError("Sample not recognized")

        heff_b, heff_c, heff_udsg = get_eff_histograms(self.year, process)

        sfs1, sfs2 = self.get_deepjet_sf("central")
        self.array['jet_1_sfs'] = sfs1
        self.array['jet_2_sfs'] = sfs2
        self.array['eff_btag'] = ak.ones_like(self.array.jet_1_pt)

        # Define selections for different scenarios
        scenarios_pass_fail = [
            {"selection": f"(jet_1_deepbtag >= {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_2_hadronflavor == 5)", "histogram": heff_b, "pass_sfs":"jet_1_sfs", "fail_sfs":"jet_2_sfs","failpt":"jet_2_pt", "faileta":"jet_2_eta"},
            {"selection": f"(jet_1_deepbtag >= {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_2_hadronflavor == 4)", "histogram": heff_c, "pass_sfs":"jet_1_sfs", "fail_sfs":"jet_2_sfs","failpt":"jet_2_pt", "faileta":"jet_2_eta"},
            {"selection": f"(jet_1_deepbtag >= {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_2_hadronflavor == 0)", "histogram": heff_udsg, "pass_sfs":"jet_1_sfs", "fail_sfs":"jet_2_sfs","failpt":"jet_2_pt", "faileta":"jet_2_eta"},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag >= {self.btagWP}) & (jet_1_hadronflavor == 5)", "histogram": heff_b, "pass_sfs":"jet_2_sfs", "fail_sfs":"jet_1_sfs","failpt":"jet_1_pt", "faileta":"jet_1_eta"},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag >= {self.btagWP}) & (jet_1_hadronflavor == 4)", "histogram": heff_c, "pass_sfs":"jet_2_sfs", "fail_sfs":"jet_1_sfs","failpt":"jet_1_pt", "faileta":"jet_1_eta"},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag >= {self.btagWP}) & (jet_1_hadronflavor == 0)", "histogram": heff_udsg, "pass_sfs":"jet_2_sfs", "fail_sfs":"jet_1_sfs","failpt":"jet_1_pt", "faileta":"jet_1_eta"},            
            # Add other scenarios here following the same structure
        ]

        scenarios_both_fail = [
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 5) & (jet_2_hadronflavor == 5)", "hist1": heff_b, "hist2": heff_b},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 5) & (jet_2_hadronflavor == 4)", "hist1": heff_b, "hist2": heff_c},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 5) & (jet_2_hadronflavor == 0)", "hist1": heff_b, "hist2": heff_udsg},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 4) & (jet_2_hadronflavor == 5)", "hist1": heff_c, "hist2": heff_b},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 4) & (jet_2_hadronflavor == 4)", "hist1": heff_c, "hist2": heff_c},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 4) & (jet_2_hadronflavor == 0)", "hist1": heff_c, "hist2": heff_udsg},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 0) & (jet_2_hadronflavor == 5)", "hist1": heff_udsg, "hist2": heff_b},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 0) & (jet_2_hadronflavor == 4)", "hist1": heff_udsg, "hist2": heff_c},
            {"selection": f"(jet_1_deepbtag < {self.btagWP}) & (jet_2_deepbtag < {self.btagWP}) & (jet_1_hadronflavor == 0) & (jet_2_hadronflavor == 0)", "hist1": heff_udsg, "hist2": heff_udsg}
        ]

        # Process all scenarios
        results_pass_fail = [calculate_weights_pass_fail(self.array, scenario["selection"], scenario["histogram"],scenario['pass_sfs'],scenario['fail_sfs'],scenario['failpt'],scenario['faileta']) for scenario in scenarios_pass_fail]
        results_both_fail = [calculate_weights_both_fail(self.array, scenario["selection"],scenario["hist1"], scenario["hist2"]) for scenario in scenarios_both_fail]
        results = results_pass_fail + results_both_fail

        # Combine all the arrays
        final_array = ak.concatenate(results)

        return final_array


                               







        