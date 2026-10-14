import ROOT
import sys
sys.path.append("/cms/user/guojl/ME_test/CMSSW_10_6_26/src/HZZAnalysis/lib")
from array import array
from logset import *

######################################## BtagEff class ########################################
class BtagEff:
    #initialize the class
    def __init__(self,year,tree,name):
        self.sample_name = name
        self.tree = tree
        self.year = year
        self.save_path = f'/cms/user/guojl/ME_test/CMSSW_10_6_26/src/HZZAnalysis/deepjetSFsfiles/{self.year}'

        #self.ptbins_bc = array('d',[30,50,70,100,140,200,300,600,1000])
        #self.ptbins_udsg = array('d',[30,50,70,100,140,200,300,600,1000])
        self.ptbins_bc = array('d',[30,50,70,100,140,200,300,1000])
        self.ptbins_udsg = array('d',[30,50,70,100,140,200,300,1000])
        #self.ptbins_udsg = array('d',[20,1000])
        #self.etabins = array('d',[-2.5,-1.5,0.0,1.5,2.5])
        self.etabins = array('d',[0.0,1.5,2.5])
        self.bins_bc = (len(self.ptbins_bc)-1,self.ptbins_bc,len(self.etabins)-1,self.etabins)
        self.bins_udsg = (len(self.ptbins_udsg)-1,self.ptbins_udsg,len(self.etabins)-1,self.etabins)

        self.hist = {}
        self.hist['pass_b'] = ROOT.TH2F('pass_b','pass_b',*self.bins_bc)
        self.hist['total_b'] = ROOT.TH2F('total_b','total_b',*self.bins_bc)
        self.hist['pass_c'] = ROOT.TH2F('pass_c','pass_c',*self.bins_bc)
        self.hist['total_c'] = ROOT.TH2F('total_c','total_c',*self.bins_bc)
        self.hist['pass_udsg'] = ROOT.TH2F('pass_udsg','pass_udsg',*self.bins_udsg)
        self.hist['total_udsg'] = ROOT.TH2F('total_udsg','total_udsg',*self.bins_udsg)
        #create the eff hist, the eff hist is the ratio of pass/total. will be filled after the event loop
        self.hist['eff_b'] = None
        self.hist['eff_c'] = None
        self.hist['eff_udsg'] = None

        #setup the btagging working point according to the year
        year_btagWP_map = {
           '2016preAPV': 0.2598,
           '2016postAPV': 0.2489,
           '2017': 0.3040,
           '2018': 0.2783
        }
        self.btagWP = year_btagWP_map.get(self.year, 0.0)

        #based on 

    def loop(self):
        r'''
        loop over the tree to fill the histograms
        '''
        for ievent,event in enumerate(self.tree):
            if ievent % 10000 == 0:
                logger.info(f'Processing event {ievent}')
            
            #make sure the event pass the selection expect the btagging requirement
            #btag:     '(Met<=150) & ((mass2jet>70) & (mass2jet<105)) & (foundZ1LCandidate==True) & (foundZ2JCandidate==True) & ((KD_jjVBF<0.5) & (KD_jVBF<0.7))  & (jet_1_deepbtag>0.2598) & (jet_2_deepbtag>0.2598)'
            #if not (event.Met<=150 and (event.mass2jet>70 and event.mass2jet<105) and event.foundZ1LCandidate==True and event.foundZ2JCandidate==True and (event.KD_jjVBF<0.5 and event.KD_jVBF<0.7)):
            if not (event.Met<=150 and event.foundZ1LCandidate==True and event.foundZ2JCandidate==True and (event.KD_jjVBF<0.5 and event.KD_jVBF<0.7)):
                continue

            #check if the event pass btag requirement then compute btag eff based on the jet flavor information
            #jet 1
            if event.jet_1_hadronflavor == 5: #b jet
                self.hist['total_b'].Fill(event.jet_1_pt,abs(event.jet_1_eta))
                if event.jet_1_deepbtag>self.btagWP:
                    self.hist['pass_b'].Fill(event.jet_1_pt,abs(event.jet_1_eta))
            elif event.jet_1_hadronflavor == 4: #c jet
                self.hist['total_c'].Fill(event.jet_1_pt,abs(event.jet_1_eta))
                if event.jet_1_deepbtag>self.btagWP:
                    self.hist['pass_c'].Fill(event.jet_1_pt,abs(event.jet_1_eta))
            else: #light jet
                self.hist['total_udsg'].Fill(event.jet_1_pt,abs(event.jet_1_eta))
                if event.jet_1_deepbtag>self.btagWP:
                    self.hist['pass_udsg'].Fill(event.jet_1_pt,abs(event.jet_1_eta))
            #jet 2
            if event.jet_2_hadronflavor == 5: #b jet
                self.hist['total_b'].Fill(event.jet_2_pt,abs(event.jet_2_eta))
                if event.jet_2_deepbtag>self.btagWP:
                    self.hist['pass_b'].Fill(event.jet_2_pt,abs(event.jet_2_eta))
            elif event.jet_2_hadronflavor == 4: #c jet
                self.hist['total_c'].Fill(event.jet_2_pt,abs(event.jet_2_eta))
                if event.jet_2_deepbtag>self.btagWP:
                    self.hist['pass_c'].Fill(event.jet_2_pt,abs(event.jet_2_eta))
            else: #light jet
                self.hist['total_udsg'].Fill(event.jet_2_pt,abs(event.jet_2_eta))
                if event.jet_2_deepbtag>self.btagWP:
                    self.hist['pass_udsg'].Fill(event.jet_2_pt,abs(event.jet_2_eta))
    
    def get_eff(self):
        r'''
        compute the btag eff and return the eff and error
        '''
        self.hist['eff_b'] = self.hist['pass_b'].Clone('eff_b')
        self.hist['eff_b'].Divide(self.hist['total_b'])
        self.hist['eff_c'] = self.hist['pass_c'].Clone('eff_c')
        self.hist['eff_c'].Divide(self.hist['total_c'])
        self.hist['eff_udsg'] = self.hist['pass_udsg'].Clone('eff_udsg')
        self.hist['eff_udsg'].Divide(self.hist['total_udsg'])
        
        #self.hist['pass_b'].Divide(self.hist['total_b'])
        #self.hist['pass_udsg'].Divide(self.hist['total_udsg'])

    def plot_eff(self):
        r'''
        plot the btag eff
        '''

        #plot for b jet
        c = ROOT.TCanvas('c','c',800,600)
        self.hist['eff_b'].SetStats(0) #remove the stat box
        self.hist['eff_b'].SetTitle('b jet btag eff')
        self.hist['eff_b'].GetXaxis().SetTitle('p_{T} [GeV]')
        self.hist['eff_b'].GetYaxis().SetTitle('#eta')
        self.hist['eff_b'].Draw('colz')
        c.SaveAs(f'{self.save_path}/bjet_btag_eff_{self.sample_name}.png')

        #plot for c jet
        c = ROOT.TCanvas('c','c',800,600)
        self.hist['eff_c'].SetStats(0)
        self.hist['eff_c'].SetTitle('c jet btag eff')
        self.hist['eff_c'].GetXaxis().SetTitle('p_{T} [GeV]')
        self.hist['eff_c'].GetYaxis().SetTitle('#eta')
        self.hist['eff_c'].Draw('colz')
        c.SaveAs(f'{self.save_path}/cjet_btag_eff_{self.sample_name}.png')

        #plot for light jet
        c = ROOT.TCanvas('c','c',800,600)
        self.hist['eff_udsg'].SetStats(0) #remove the stat box
        self.hist['eff_udsg'].SetTitle('light jet btag eff')
        self.hist['eff_udsg'].GetXaxis().SetTitle('p_{T} [GeV]')
        self.hist['eff_udsg'].GetYaxis().SetTitle('#eta')
        self.hist['eff_udsg'].Draw('colz')
        c.SaveAs(f'{self.save_path}/lightjet_btag_eff_{self.sample_name}.png')

        #plot total events for b c and light 
        for name in ['total_b','total_c','total_udsg']:
            c = ROOT.TCanvas('c','c',800,600)
            self.hist[name].SetStats(0)
            self.hist[name].SetTitle(f'{name}')
            self.hist[name].GetXaxis().SetTitle('p_{T} [GeV]')
            self.hist[name].GetYaxis().SetTitle('#eta')
            self.hist[name].Draw('colz')
            self.hist[name].GetZaxis().SetLabelSize(0.02)
            c.SaveAs(f'{self.save_path}/{name}_{self.sample_name}.png')
            
    
    def save_eff_histo(self):
        r'''
        save the eff histos
        '''
        file = ROOT.TFile(f'{self.save_path}/btag_eff_{self.sample_name}.root','recreate')
        self.hist['pass_b'].Write()
        self.hist['total_b'].Write()

        self.hist['pass_udsg'].Write()
        self.hist['total_udsg'].Write()

        self.hist['pass_c'].Write()
        self.hist['total_c'].Write()

        self.hist['eff_c'].Write()
        self.hist['eff_b'].Write()
        self.hist['eff_udsg'].Write()
        file.Close()

    def run(self):
        r'''
        run the btag eff
        '''
        self.loop()
        self.get_eff()
        self.plot_eff()
        self.save_eff_histo()


###########################################################################################
######################################## SFs class ########################################
###########################################################################################
#SFs class is subclass of BtagEff class        
class SFs(BtagEff):
    #initialize the class
    def __init__(self,year,file_path,name,WP = 'M'):
        #call the __init__ function of the parent class
        super().__init__(year,file_path,name)
        self.WP = WP
        self.type = 'central'
        self.bcjet_sf_name = 'deepJet_comb'
        self.lightjet_sf_name = 'deepJet_incl'

        self.file_sf_path = f"/cms/user/guojl/ME_test/CMSSW_10_6_26/src/HZZAnalysis/deepjetSFsfiles/{year}/btagging.json.gz"
        self.file_sf = correctionlib.CorrectionSet.from_file(self.file_sf_path)

        process = self.determine_process(name)
        self.file_eff_path = f"/cms/user/guojl/ME_test/CMSSW_10_6_26/src/HZZAnalysis/deepjetSFsfiles/{year}/btag_eff_{process}.root"
        #self.file_eff_path = f"/cms/user/guojl/ME_test/CMSSW_10_6_26/src/HZZAnalysis/deepjetSFsfiles/{year}/btag_eff_{name}.root"
        self.file_eff = ROOT.TFile(self.file_eff_path)

        #define bins for SFs. range: [0,1], bins: 50
        self.bins_sfs = array('d',[0,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.42,0.44,0.46,0.48,0.5,0.52,0.54,0.56,0.58,0.6,0.62,0.64,0.66,0.68,0.7,0.72,0.74,0.76,0.78,0.8,0.82,0.84,0.86,0.88,0.9,0.92,0.94,0.96,0.98,1])
        self.bins_2d_sfs = (len(self.ptbins_bc)-1,self.ptbins_bc,len(self.bins_sfs)-1,self.bins_sfs)

        self.hist = {}
        #extract the eff hist from the eff file
        self.hist['eff_b'] = self.file_eff.Get('eff_b')
        self.hist['eff_udsg'] = self.file_eff.Get('eff_udsg')
        self.hist['eff_c'] = self.file_eff.Get('eff_c')
        #define SFs hist
        self.hist['SFs_jet1'] = ROOT.TH2F('SFs_jet1','SFs_jet1',*self.bins_2d_sfs)
        self.hist['SFs_jet1_udsg'] = ROOT.TH2F('SFs_jet1_udsg','SFs_jet1_udsg',*self.bins_2d_sfs)
        self.hist['SFs_jet2'] = ROOT.TH2F('SFs_jet2','SFs_jet2',*self.bins_2d_sfs)
        self.hist['SFs_jet2_udsg'] = ROOT.TH2F('SFs_jet2_udsg','SFs_jet2_udsg',*self.bins_2d_sfs)

        self.hist['SFs_weight'] = ROOT.TH1D('SFs_weight','SFs_weight',100,0.5,1.5)

        self.if_test_event = False
        self.test_event_num = 50000

        self.updatefile_path = file_path

    def determine_process(self,sample):
        if sample in ['TTTo2L2Nu', 'tZq','ZZTo2Q2L', 'WZTo2Q2L','ggh3000','zv']:
            return sample
        elif sample in ['WWTo2L2Nu']:
            return 'TTTo2L2Nu'
        return None


    def loop(self):
        r'''
        loop over the tree to fill the histograms
        The scale factors are applied as a multiplicative event weight, 
        considering all jets where the b-tagging information is used in the analysis
        update the file by adding weight into the tree as a new branch
        '''
        #open the file and get the tree
        file = ROOT.TFile(self.updatefile_path,'update')
        tree = file.Get('passedEvents')

        #define the new branch
        b_weight = array('f',[0])
        b_weight_up = array('f',[0])
        b_weight_dn = array('f',[0])
        #if the branch is already exist, update the value for very entry, else ceate a new branch (branch name: w_deepjet)
        if tree.GetBranch('w_deepjet'):
            #branch = tree.GetBranch('w_deepjet')
            #branch.SetAddress(b_weight)
            tree.SetBranchAddress('w_deepjet',b_weight)
            tree.SetBranchAddress('w_deepjet_up',b_weight_up)
            tree.SetBranchAddress('w_deepjet_dn',b_weight_dn)
            updatebranch = True
            logger.info(f"update the branch")
        else:
            branch = tree.Branch('w_deepjet',b_weight,'w_deepjet/F')
            branch_up = tree.Branch('w_deepjet_up',b_weight_up,'w_deepjet_up/F')
            branch_dn = tree.Branch('w_deepjet_dn',b_weight_dn,'w_deepjet_dn/F')
            updatebranch = False
            logger.info(f"create a new branch")



        #branch = tree.Branch('w_deepjet',b_weight,'w_deepjet/F')


        for ievent,event in enumerate(tree):
            if ievent % 10000 == 0:
                logger.info(f'Processing event {ievent}')
            #test with 5000 events if needed
            if self.if_test_event:
                if ievent == self.test_event_num:
                    break
            #initialize the weight
            b_weight[0] = -999

            #make sure the event pass the selection expect the btagging requirement
            #if not (event.Met<=150 and (event.mass2jet>70 and event.mass2jet<105) and event.foundZ1LCandidate==True and event.foundZ2JCandidate==True and (event.KD_jjVBF<0.5 and event.KD_jVBF<0.7)):
            if not (event.Met<=150 and event.foundZ1LCandidate==True and event.foundZ2JCandidate==True and (event.KD_jjVBF<0.5 and event.KD_jVBF<0.7)):
                if not updatebranch:
                    branch.Fill()
                    branch_up.Fill()
                    branch_dn.Fill()
                continue
            
            #with each btag event, get SFs based on jet flavor information
            #check with b tagged cut 
            #if event.jet_1_deepbtag>self.btagWP and event.jet_2_deepbtag>self.btagWP:
            #jet 1
            if event.jet_1_hadronflavor == 5 or event.jet_1_hadronflavor == 4: #bcjet
                sf_jet1 = self.file_sf[self.bcjet_sf_name].evaluate('central',self.WP,event.jet_1_hadronflavor,abs(event.jet_1_eta),event.jet_1_pt)
                sf_jet1_up = self.file_sf[self.bcjet_sf_name].evaluate('up',self.WP,event.jet_1_hadronflavor,abs(event.jet_1_eta),event.jet_1_pt)
                sf_jet1_dn = self.file_sf[self.bcjet_sf_name].evaluate('down',self.WP,event.jet_1_hadronflavor,abs(event.jet_1_eta),event.jet_1_pt)

                logger.debug(f"This bc case jet1: sf_jet1: {sf_jet1}, sf_jet1_up: {sf_jet1_up}, sf_jet1_dn: {sf_jet1_dn}")
            else:
                sf_jet1 = self.file_sf[self.lightjet_sf_name].evaluate('central',self.WP,event.jet_1_hadronflavor,abs(event.jet_1_eta),event.jet_1_pt)
                sf_jet1_up = self.file_sf[self.lightjet_sf_name].evaluate('up',self.WP,event.jet_1_hadronflavor,abs(event.jet_1_eta),event.jet_1_pt)
                sf_jet1_dn = self.file_sf[self.lightjet_sf_name].evaluate('down',self.WP,event.jet_1_hadronflavor,abs(event.jet_1_eta),event.jet_1_pt)

                logger.debug(f"This light case jet1: sf_jet1: {sf_jet1}, sf_jet1_up: {sf_jet1_up}, sf_jet1_dn: {sf_jet1_dn}")

            self.hist['SFs_jet1'].Fill(event.jet_1_pt,sf_jet1)

            #jet 2
            if event.jet_2_hadronflavor == 5 or event.jet_2_hadronflavor == 4: #bcjet
                sf_jet2 = self.file_sf[self.bcjet_sf_name].evaluate('central',self.WP,event.jet_2_hadronflavor,abs(event.jet_2_eta),event.jet_2_pt)
                sf_jet2_up = self.file_sf[self.bcjet_sf_name].evaluate('up',self.WP,event.jet_2_hadronflavor,abs(event.jet_2_eta),event.jet_2_pt)
                sf_jet2_dn = self.file_sf[self.bcjet_sf_name].evaluate('down',self.WP,event.jet_2_hadronflavor,abs(event.jet_2_eta),event.jet_2_pt)

                logger.debug(f"This bc case jet2: sf_jet2: {sf_jet2}, sf_jet2_up: {sf_jet2_up}, sf_jet2_dn: {sf_jet2_dn}")
            else:
                sf_jet2 = self.file_sf[self.lightjet_sf_name].evaluate('central',self.WP,event.jet_2_hadronflavor,abs(event.jet_2_eta),event.jet_2_pt)
                sf_jet2_up = self.file_sf[self.lightjet_sf_name].evaluate('up',self.WP,event.jet_2_hadronflavor,abs(event.jet_2_eta),event.jet_2_pt)
                sf_jet2_dn = self.file_sf[self.lightjet_sf_name].evaluate('down',self.WP,event.jet_2_hadronflavor,abs(event.jet_2_eta),event.jet_2_pt)

                logger.debug(f"This light case jet2: sf_jet2: {sf_jet2}, sf_jet2_up: {sf_jet2_up}, sf_jet2_dn: {sf_jet2_dn}")
                
            self.hist['SFs_jet2'].Fill(event.jet_2_pt,sf_jet2)

            #with each event, copmute the weight. based on jet pass or fail the btagging requirement
            #case 1: both jet pass the btagging requirement
            if event.jet_1_deepbtag>=self.btagWP and event.jet_2_deepbtag>=self.btagWP:
                logger.debug(f"Both jet pass the btagging requirement")
                weight = sf_jet1*sf_jet2
                weight_up = sf_jet1_up*sf_jet2_up
                weight_dn = sf_jet1_dn*sf_jet2_dn
                logger.debug(f"weight: {weight}")
            #case 2: jet 1 pass the btagging requirement, jet 2 fail
            elif event.jet_1_deepbtag>=self.btagWP and event.jet_2_deepbtag<self.btagWP:
                logger.debug(f"Jet 1 pass the btagging requirement, jet 2 fail")
                if event.jet_2_hadronflavor == 5: #bjet
                    #hist['eff_*'] is the eff of passing the btagging requirement, so the eff_fail is 1-eff
                    #eff_fail = 1-self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    eff_fail = self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    logger.debug(f"eff_fail: {eff_fail} with jet2 is bjet")
                elif event.jet_2_hadronflavor == 4: #cjet
                    #eff_fail = 1-self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    eff_fail = self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    logger.debug(f"eff_fail: {eff_fail} with jet2 is cjet")
                else: #light jet
                    #eff_fail = 1-self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    eff_fail = self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    logger.debug(f"eff_fail: {eff_fail} with jet2 is light jet")

                weight = (sf_jet1*(1-sf_jet2*eff_fail))/(1-eff_fail)
                weight_up = (sf_jet1_up*(1-sf_jet2_up*eff_fail))/(1-eff_fail)
                weight_dn = (sf_jet1_dn*(1-sf_jet2_dn*eff_fail))/(1-eff_fail)
                logger.debug(f"weight: {weight}")

            #case 3: jet 1 fail the btagging requirement, jet 2 pass
            elif event.jet_1_deepbtag<self.btagWP and event.jet_2_deepbtag>=self.btagWP:
                logger.debug(f"Jet 1 fail the btagging requirement, jet 2 pass")
                if event.jet_1_hadronflavor == 5: #bcjet
                    #eff_fail = 1-self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    eff_fail = self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    logger.debug(f"eff_fail: {eff_fail} with jet1 is bcjet")
                elif event.jet_1_hadronflavor == 4: #cjet
                    #eff_fail = 1-self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    eff_fail = self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    logger.debug(f"eff_fail: {eff_fail} with jet1 is cjet")
                else: #light jet
                    #eff_fail = 1-self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    eff_fail = self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    logger.debug(f"eff_fail: {eff_fail} with jet1 is light jet")

                weight = (sf_jet2*(1-sf_jet1*eff_fail))/(1-eff_fail)
                weight_up = (sf_jet2_up*(1-sf_jet1_up*eff_fail))/(1-eff_fail)
                weight_dn = (sf_jet2_dn*(1-sf_jet1_dn*eff_fail))/(1-eff_fail)
                logger.debug(f"weight: {weight}")

            #case 4: both jet fail
            else:
                logger.debug(f"Both jet fail the btagging requirement")
                if event.jet_1_hadronflavor == 5:
                    #eff_fail_jet1 = 1-self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    eff_fail_jet1 = self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    logger.debug(f"eff_fail_jet1: {eff_fail_jet1} with jet1 is bjet")
                elif event.jet_1_hadronflavor == 4:
                    #eff_fail_jet1 = 1-self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    eff_fail_jet1 = self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    logger.debug(f"eff_fail_jet1: {eff_fail_jet1} with jet1 is cjet")
                else:
                    #eff_fail_jet1 = 1-self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    eff_fail_jet1 = self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_1_pt,abs(event.jet_1_eta)))
                    logger.debug(f"eff_fail_jet1: {eff_fail_jet1} with jet1 is light jet")

                if event.jet_2_hadronflavor == 5 or event.jet_2_hadronflavor == 4:
                    #eff_fail_jet2 = 1-self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    eff_fail_jet2 = self.hist['eff_b'].GetBinContent(self.hist['eff_b'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    logger.debug(f"eff_fail_jet2: {eff_fail_jet2} with jet2 is bcjet")
                elif event.jet_2_hadronflavor == 4:
                    #eff_fail_jet2 = 1-self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    eff_fail_jet2 = self.hist['eff_c'].GetBinContent(self.hist['eff_c'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    logger.debug(f"eff_fail_jet2: {eff_fail_jet2} with jet2 is cjet")
                else:
                    #eff_fail_jet2 = 1-self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    eff_fail_jet2 = self.hist['eff_udsg'].GetBinContent(self.hist['eff_udsg'].FindBin(event.jet_2_pt,abs(event.jet_2_eta)))
                    logger.debug(f"eff_fail_jet2: {eff_fail_jet2} with jet2 is light jet")

                weight = (1-sf_jet1*eff_fail_jet1)*(1-sf_jet2*eff_fail_jet2)/(1-eff_fail_jet1)/(1-eff_fail_jet2)
                weight_up = (1-sf_jet1_up*eff_fail_jet1)*(1-sf_jet2_up*eff_fail_jet2)/(1-eff_fail_jet1)/(1-eff_fail_jet2)
                weight_dn = (1-sf_jet1_dn*eff_fail_jet1)*(1-sf_jet2_dn*eff_fail_jet2)/(1-eff_fail_jet1)/(1-eff_fail_jet2)
                logger.debug(f"weight: {weight}")

            #fill the weigt hist
            self.hist['SFs_weight'].Fill(weight)
            #change value of b_weight to weight
            b_weight[0] = weight
            b_weight_up[0] = weight_up
            b_weight_dn[0] = weight_dn
            #fill the new branch
            if not updatebranch:
                branch.Fill()
                branch_up.Fill()
                branch_dn.Fill()

        #write the tree
        tree.Write("",ROOT.TObject.kOverwrite) # save only the new version of the tree
        file.Close()
            

    def plot_sfs(self):
        r'''
        plot the SFs, and SFs_weight
        '''
        #plot for SFs jet
        for name in ['SFs_jet1','SFs_jet2']:
            c = ROOT.TCanvas('c','c',800,600)
            self.hist[name].SetStats(0)
            self.hist[name].SetTitle(f'{name}')
            self.hist[name].GetXaxis().SetTitle('p_{T} [GeV]')
            self.hist[name].GetYaxis().SetTitle('SFs')
            self.hist[name].Draw('colz')
            c.SaveAs(f'{self.save_path}/{name}_{self.sample_name}.png')
        #plot for SFs weight
        c = ROOT.TCanvas('c','c',800,600)
        self.hist['SFs_weight'].SetStats(0)
        self.hist['SFs_weight'].SetTitle('SFs_weight')
        self.hist['SFs_weight'].GetXaxis().SetTitle('weight')
        self.hist['SFs_weight'].Draw()
        c.SaveAs(f'{self.save_path}/SFs_weight_{self.sample_name}.png')

        pass
        
    def run(self):
        r'''
        run the SFs
        '''
        self.loop()
        self.plot_sfs()

#define this scrpit as main function
if __name__ == '__main__':
    #define some inpit arugments. 1.year 2.sample name 3.log level
    import argparse
    parser = argparse.ArgumentParser(description='btageff. Usage: python3 btageff.py -y 2018 -n ggh -l info --option eff')
    parser.add_argument('-y','--year',dest='year',default='2016',help='year to run or run for all three year. Options: 2016, 2016preAPV, 2016postAPV,2017, 2018, all')
    parser.add_argument('-n','--name',dest='name',default='all',help='sample name that you want to run. only DY_pt50To100,DY_pt100To250,DY_pt400To650,DY_pt650ToInf,TTJets,TTTo2L2Nu,WW_TuneCP5,WZTo2Q2L,ZZTo2Q2L,Signal are available. if you want to run all you can skip this argument')
    parser.add_argument('-l','--log',dest='log',default='debug',help='log level. Options: debug, info, warning, error, critical.')
    parser.add_argument('-o','--option',dest = 'option',default='eff',help='run eff or sfs')
    args = parser.parse_args()

    #set log level
    if args.log == 'debug':
        logger.setLevel(logging.DEBUG)
        logger.debug("set log level to debug")
    elif args.log == 'info':
        logger.setLevel(logging.INFO)
        print("set log level to info")
    elif args.log == 'warning':
        logger.setLevel(logging.WARNING)
    elif args.log == 'error':
        logger.setLevel(logging.ERROR)
    elif args.log == 'critical':
        logger.setLevel(logging.CRITICAL)
    else:
        logger.setLevel(logging.ERROR)

    samples_dirs = {'DY_pt50To100':'MC/DYJetsToLL_LHEFilterPtZ-50To100_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/DY_pt50To100.root',
                    'DY_pt100To250':'MC/DYJetsToLL_LHEFilterPtZ-100To250_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/DY_pt100To250.root',
                    'DY_pt250To400':'MC/DYJetsToLL_LHEFilterPtZ-250To400_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/DY_pt250To400.root',
                    'DY_pt400To650':'MC/DYJetsToLL_LHEFilterPtZ-400To650_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/DY_pt400To650.root',
                    'DY_pt650ToInf':'MC/DYJetsToLL_LHEFilterPtZ-650ToInf_MatchEWPDG20_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/DY_pt650ToInf.root',
                    #'TTJets':'MC/TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8',
                    'TTTo2L2Nu':'MC/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/skimed/TTTo2L2Nu.root',
                    #'WW_TuneCP5':'MC/WW_TuneCP5_13TeV-pythia8',
                    'WWTo2L2Nu':'MC/WWTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/skimed/WWTo2L2Nu.root',
                    'WZTo2Q2L':'MC/WZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/WZTo2Q2L.root',
                    'ZZTo2Q2L':'MC/ZZTo2Q2L_mllmin4p0_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/ZZTo2Q2L.root',
                    'ZZTo2L2Nu':'MC/ZZTo2L2Nu_TuneCP5_13TeV_powheg_pythia8/skimed/ZZTo2L2Nu.root',
                    'tZq':"MC/tZq_ll_4f_ckm_NLO_TuneCP5_13TeV-amcatnlo-pythia8/skimed/tZq.root",
                    'top':"MC/top.root",
                    'zv':"MC/zv.root",
                    'WZTo1L1Nu2Q':'MC/WZTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/WZTo1L1Nu2Q.root',
                    'WWTo1L1Nu2Q':'MC/WWTo1L1Nu2Q_4f_TuneCP5_13TeV-amcatnloFXFX-pythia8/skimed/WWTo1L1Nu2Q.root',
                    #'ZH_HToBB':'MC/ZH_HToBB_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8',
                    #'ZH_HToCC':'MC/ZH_HToCC_ZToLL_M-125_TuneCP5_13TeV-powheg-pythia8',
                    #'ZH_HToGG':'MC/ZH_HToGG_ZToAll_M125_TuneCP5_13TeV-powheg-pythia8',
                    #'ggh':'MC/ggh/skimed/GluGluHToZZTo2L2Q_M1000_TuneCP5_13TeV_powheg2_JHUGenV7011_pythia8__v16_L1v1-v1_0.root',
                    'ggh':'MC/ggh/skimed/ggh.root',
                    #'ggh125':'MC/ggh125',
                    #'ggh300':'MC/ggh300',
                    #'ggh350':'MC/ggh350',
                    #'ggh400':'MC/ggh400',
                    #'ggh450':'MC/ggh450',
                    #'ggh500':'MC/ggh500',
                    #'ggh550':'MC/ggh550',
                    #'ggh600':'MC/ggh600',
                    #'ggh700':'MC/ggh700',
                    #'ggh800':'MC/ggh800',
                    #'ggh900':'MC/ggh900',
                    #'ggh1000':'MC/ggh1000',
                    #'ggh1500':'MC/ggh1500',
                    'ggh3000':'MC/ggh3000/skimed/ggh3000.root',
                    #'ggh2500':'MC/ggh2500',
                    'vbf':'MC/vbf',
                    #'test':'MC/test',
                    'Data':'Data',
                    }

    #fine the root file path acorrind to the input arguments,only for background MC samples
    path_sample = '/cms/user/guojl/Sample/2L2Q/UL_Legacy/'+args.year+'/'+samples_dirs[args.name]

    logger.info(f'Processing sample {args.name} in year {args.year}')


    #run eff or SFs
    if args.option == 'eff':
        #from root file to get the tree named "passedEvents"
        file = ROOT.TFile(path_sample)
        tree = file.Get('passedEvents')

        #create a objet of class BtagEff
        btag = BtagEff(args.year,tree,args.name)

        #run the btag eff
        btag.run()
    elif args.option == 'sfs':
        #import required packages
        import correctionlib
        import numpy as np
        #create a object of class SFs
        sfs = SFs(args.year,path_sample,args.name)
    
        #run the SFs
        sfs.run()

    else:
        logger.error('Invalid option. Please choose eff or sfs')
        sys.exit(1)