import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
bphEfficiency = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/*"),
    verbose        = cms.untracked.uint32(0), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'efficiency vs #phi^{\mu}; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'efficiency vs #eta^{#mu}; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'efficiency vs p_{T}^{#mu}; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
        "effic_phPhi       'efficiency vs #phi^{#gamma}; ph phi [rad]; efficiency' phPhi_numerator       phPhi_denominator",
        "effic_phEta       'efficiency vs #eta^{#gamma}; ph eta [rad]; efficiency' phEta_numerator       phEta_denominator",
        "effic_phPt       'efficiency vs p_{T}^{#gamma}; ph pt [GeV]; efficiency' phPt_numerator       phPt_denominator",
        "effic_trPhi       'efficiency vs #phi^{tk}; tr phi [rad]; efficiency' trPhi_numerator       trPhi_denominator",
        "effic_trEta       'efficiency vs #eta^{tk}; tr eta [rad]; efficiency' trEta_numerator       trEta_denominator",
        "effic_trPt       'efficiency vs p_{T}^{tk}; tr pt [GeV]; efficiency' trPt_numerator       trPt_denominator",
        "effic_mu1Phi       'efficiency vs #phi^{#mu_1}; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'efficiency vs eta^{#mu_1}; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'efficiency vs p_{T}^{#mu_1}; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'efficiency vs #phi^{#mu_2}; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'efficiency vs #eta^{#mu_2}; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'efficiency vs p_{T}^{#mu_2}; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_mu3Phi       'efficiency vs #phi^{#mu_3}; mu3 phi [rad]; efficiency' mu3Phi_numerator       mu3Phi_denominator",
        "effic_mu3Eta       'efficiency vs #eta^{#mu_3}; mu3 eta [rad]; efficiency' mu3Eta_numerator       mu3Eta_denominator",
        "effic_mu3Pt       'efficiency vs p_{T}^{#mu_3}; mu3 pt [GeV]; efficiency' mu3Pt_numerator       mu3Pt_denominator",
        "effic_DiMuPhi       'efficiency vs #phi^{dimuon}; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'efficiency vs #eta{dimuon}; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'efficiency vs p_{T}^{dimuon}; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuPVcos       'efficiency vs cosPV; DiMu cosPV ; efficiency' DiMuPVcos_numerator       DiMuPVcos_denominator",
        "effic_DiMuProb       'efficiency vs prob; DiMu prob ; efficiency' DiMuProb_numerator       DiMuProb_denominator",
        "effic_DiMuDS       'efficiency vs DS; DiMu DS; efficiency' DiMuDS_numerator       DiMuDS_denominator",
        "effic_DiMuDCA       'efficiency vs Dimuon DCA; DiMu DCA [cm]; efficiency' DiMuDCA_numerator       DiMuDCA_denominator",
        "effic_DiMuMass       'efficiency vs Dimuon Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_BMass       'efficiency vs B Mass; B Mass[GeV]; efficiency' BMass_numerator       BMass_denominator",
        "effic_DiMudR       'efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",
        "effic_tr1Phi       'efficiency vs #phi^{tk1}; tr1 phi [rad]; efficiency' tr1Phi_numerator       tr1Phi_denominator",
        "effic_tr1Eta       'efficiency vs #eta^{tk1}; tr1 eta [rad]; efficiency' tr1Eta_numerator       tr1Eta_denominator",
        "effic_tr1Pt       'efficiency vs p_{T}^{tk1}; tr1 pt [GeV]; efficiency' tr1Pt_numerator       tr1Pt_denominator",
        "effic_tr2Phi       'efficiency vs #phi^{tk2}; tr2 phi [rad]; efficiency' tr2Phi_numerator       tr2Phi_denominator",
        "effic_tr2Eta       'efficiency vs #eta^{tk2}; tr2 eta [rad]; efficiency' tr2Eta_numerator       tr2Eta_denominator",
        "effic_tr2Pt       'efficiency vs p_{T}^{tk2}; tr2 pt [GeV]; efficiency' tr2Pt_numerator       tr2Pt_denominator",
    ),
#    efficiencyProfile = cms.untracked.vstring(
#        "effic_met_vs_LS 'MET efficiency vs LS; LS; PF MET efficiency' metVsLS_numerator metVsLS_denominator"
#    ),
  
)

bphClient = cms.Sequence(
    bphEfficiency
)

##

##
