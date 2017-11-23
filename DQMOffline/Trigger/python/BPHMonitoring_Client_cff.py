import FWCore.ParameterSet.Config as cms
from DQMServices.Core.DQMEDHarvester import DQMEDHarvester
bphEfficiency1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_Jpsi_L1_NO_OS_denTrack2/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency1_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu25_Jpsi_noCorr/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency1_2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_Upsilon_L1_NO_OS_denTrack2/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency1_3 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_Upsilon_L1_NO_OS_denTrack7"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency1_4 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_Jpsi_L1_NO_OS_denTrack7"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)



bphEfficiency1_5 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_Jpsi_L1_NO_OS_denTrack3p5/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency1_6 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_Upsilon_L1_NO_OS_denTrack3p5/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency1_7 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_L1_L3TnP_Jpsi/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency1_8 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/1_DiMu0_L1_L3TnP_Upsilon/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)


bphEfficiency2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/2_DiMu0_Jpsi_L1_OS/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",

    ),
)


bphEfficiency2_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/2_DiMu0_Jpsi_L1_HLT_OS/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",

    ),
)

bphEfficiency2_2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/2_DiMu0_Jpsi_L1_HLT_OS1/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",

    ),
)
#bphEfficiency2_3 = DQMEDHarvester("DQMGenericClient",
#    subDirs        = cms.untracked.vstring("HLT/BPH/2_DimuX_HLT_OS_Vtx/"),
#    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
#    resolution     = cms.vstring(),
#    efficiency     = cms.vstring(
#        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
#        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
#        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
#        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
#        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
#        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
#        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
#        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
#        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
#
#    ),
#)


bphEfficiency3 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/3_DiMu0_Upsilon_L1_er/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",

    ),
)

bphEfficiency3_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/3_DiMu0_Lowmass_L1_er/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",

    ),
)

bphEfficiency4 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu0_Lowmass_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu0_Lowmass_L1_dR_low/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)


bphEfficiency4_2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_3_Bs_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_3 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_3_Jpsi_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_4 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu14_Phi_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)
bphEfficiency4_5 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu20_Jpsi_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_6 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu10_PsiPrime_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)
bphEfficiency4_7 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_LowMassNonResonantTrk_Displaced_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)
bphEfficiency4_8 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_LowMassNonResonantTrk_Displaced_L1_dR_low/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_9 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_JpsiTrk_Displaced_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_10 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_JpsiTrk_Displaced_L1_dR_low/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_11 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_PsiPrimeTrk_Displaced_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_12 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DMu4_PsiPrimeTrk_Displaced_L1_dR_low/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_13 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu25_Jpsi_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_14 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu18_PsiPrime_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)
bphEfficiency4_15 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu12_Upsilon_L1_dR/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_16 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu25_Jpsi_L1_dR_low/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency4_17 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/4_DiMu18_Jpsi_L1_dR_low/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)



bphEfficiency5 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/5_DiMu0_Upsilon_L1_masscut/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency5_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/5_DiMu20_Upsilon_L1_masscut1/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)
bphEfficiency5_2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/5_DiMu12_Upsilon_L1_masscut2/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency5_3 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/5_DiMu0_Lowmass_L1_masscut3/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)
bphEfficiency5_4 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/5_TripleMu2_Upsilon_L1_masscut4"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency5_5 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/5_DoubleMu3_Trk_L1_masscut5"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)

bphEfficiency5_6 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/5_DoubleMu3_Trk_L1_masscut6"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",

    ),
)











bphEfficiency6 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/6_DiMu0_Lowmass_L1_tripleMu1/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_mu3Phi       'mu3 efficiency vs phi; mu3 phi [rad]; efficiency' mu3Phi_numerator       mu3Phi_denominator",
        "effic_mu3Eta       'mu3 efficiency vs eta; mu3 eta [rad]; efficiency' mu3Eta_numerator       mu3Eta_denominator",
        "effic_mu3Pt       'mu3 efficiency vs pt; mu3 pt [GeV]; efficiency' mu3Pt_numerator       mu3Pt_denominator",

    ),
)


bphEfficiency6_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/6_DiMu0_Lowmass_L1_tripleMu2/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_mu3Phi       'mu3 efficiency vs phi; mu3 phi [rad]; efficiency' mu3Phi_numerator       mu3Phi_denominator",
        "effic_mu3Eta       'mu3 efficiency vs eta; mu3 eta [rad]; efficiency' mu3Eta_numerator       mu3Eta_denominator",
        "effic_mu3Pt       'mu3 efficiency vs pt; mu3 pt [GeV]; efficiency' mu3Pt_numerator       mu3Pt_denominator",

    ),
)

bphEfficiency6_2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/6_DiMu0_Lowmass_L1_tripleMu3/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
        "effic_mu3Phi       'mu3 efficiency vs phi; mu3 phi [rad]; efficiency' mu3Phi_numerator       mu3Phi_denominator",
        "effic_mu3Eta       'mu3 efficiency vs eta; mu3 eta [rad]; efficiency' mu3Eta_numerator       mu3Eta_denominator",
        "effic_mu3Pt       'mu3 efficiency vs pt; mu3 pt [GeV]; efficiency' mu3Pt_numerator       mu3Pt_denominator",

    ),
)


#bphEfficiency7 = DQMEDHarvester("DQMGenericClient",
#    subDirs        = cms.untracked.vstring("HLT/BPH/7_DiMu0_Lowmass_L1_photon1/"),
#    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
#    resolution     = cms.vstring(),
#    efficiency     = cms.vstring(
#        "effic_phPhi       'ph efficiency vs phi; ph phi [rad]; efficiency' phPhi_numerator       phPhi_denominator",
#        "effic_phEta       'ph efficiency vs eta; ph eta [rad]; efficiency' phEta_numerator       phEta_denominator",
#        "effic_phPt       'ph efficiency vs pt; ph pt [GeV]; efficiency' phPt_numerator       phPt_denominator",
#    ),
#)

#bphEfficiency7 = DQMEDHarvester("DQMGenericClient",
#    subDirs        = cms.untracked.vstring("HLT/BPH/7*"),
#    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
#    resolution     = cms.vstring(),
#    efficiency     = cms.vstring(
#        "effic_phPhi       'ph efficiency vs phi; ph phi [rad]; efficiency' phPhi_numerator       phPhi_denominator",
#        "effic_phEta       'ph efficiency vs eta; ph eta [rad]; efficiency' phEta_numerator       phEta_denominator",
#        "effic_phPt       'ph efficiency vs pt; ph pt [GeV]; efficiency' phPt_numerator       phPt_denominator",
#    ),
#)


bphEfficiency8 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/8_DiMu0_L1_looseVtx_Jpsi/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuPVcos       'DiMu efficiency vs cosPV; DiMu cosPV ; efficiency' DiMuPVcos_numerator       DiMuPVcos_denominator",
        "effic_DiMuProb       'DiMu efficiency vs prob; DiMu prob ; efficiency' DiMuProb_numerator       DiMuProb_denominator",
        "effic_DiMuDS       'DiMu efficiency vs DS; DiMu DS; efficiency' DiMuDS_numerator       DiMuDS_denominator",
        "effic_DiMuDCA       'DiMu efficiency vs DCA; DiMu DCA [cm]; efficiency' DiMuDCA_numerator       DiMuDCA_denominator",

    ),
)

bphEfficiency8_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/8_DiMu0_L1_looseVtx_Upsilon/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuPVcos       'DiMu efficiency vs cosPV; DiMu cosPV ; efficiency' DiMuPVcos_numerator       DiMuPVcos_denominator",
        "effic_DiMuProb       'DiMu efficiency vs prob; DiMu prob ; efficiency' DiMuProb_numerator       DiMuProb_denominator",
        "effic_DiMuDS       'DiMu efficiency vs DS; DiMu DS; efficiency' DiMuDS_numerator       DiMuDS_denominator",
        "effic_DiMuDCA       'DiMu efficiency vs DCA; DiMu DCA [cm]; efficiency' DiMuDCA_numerator       DiMuDCA_denominator",

    ),
)

bphEfficiency8_2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/8_DiMu0_L1_tightVtx_Jpsi/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
        "effic_DiMuPVcos       'DiMu efficiency vs cosPV; DiMu cosPV ; efficiency' DiMuPVcos_numerator       DiMuPVcos_denominator",
        "effic_DiMuProb       'DiMu efficiency vs prob; DiMu prob ; efficiency' DiMuProb_numerator       DiMuProb_denominator",
        "effic_DiMuDS       'DiMu efficiency vs DS; DiMu DS; efficiency' DiMuDS_numerator       DiMuDS_denominator",
        "effic_DiMuDCA       'DiMu efficiency vs DCA; DiMu DCA [cm]; efficiency' DiMuDCA_numerator       DiMuDCA_denominator",

    ),
)



bphEfficiency9 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/9_DiMu0_L1_addTrack_Jpsi/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_trPhi       'tr efficiency vs phi; tr phi [rad]; efficiency' trPhi_numerator       trPhi_denominator",
        "effic_trEta       'tr efficiency vs eta; tr eta [rad]; efficiency' trEta_numerator       trEta_denominator",
        "effic_trPt       'tr efficiency vs pt; tr pt [GeV]; efficiency' trPt_numerator       trPt_denominator",

    ),
)


bphEfficiency10 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/10_DiMu0_L1_addTrackMu_Phi/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency10_1 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/10_DiMu0_L1_addTrackMu_Onia/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)
bphEfficiency10_2 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/10_DiMu0_L1_addTrackMu_Phi1/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)
bphEfficiency10_3 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/10_DiMu0_L1_addTrackMu_Onia1/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_muPhi       'mu efficiency vs phi; mu phi [rad]; efficiency' muPhi_numerator       muPhi_denominator",
        "effic_muEta       'mu efficiency vs eta; mu eta [rad]; efficiency' muEta_numerator       muEta_denominator",
        "effic_muPt       'mu efficiency vs pt; mu pt [GeV]; efficiency' muPt_numerator       muPt_denominator",
    ),
)

bphEfficiency11 = DQMEDHarvester("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/BPH/11_DiMu0_L1_addTrackTrack_Jpsi/"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(),
    efficiency     = cms.vstring(
        "effic_tr1Phi       'tr1 efficiency vs phi; tr1 phi [rad]; efficiency' tr1Phi_numerator       tr1Phi_denominator",
        "effic_tr1Eta       'tr1 efficiency vs eta; tr1 eta [rad]; efficiency' tr1Eta_numerator       tr1Eta_denominator",
        "effic_tr1Pt       'tr1 efficiency vs pt; tr1 pt [GeV]; efficiency' tr1Pt_numerator       tr1Pt_denominator",
        "effic_tr2Phi       'tr2 efficiency vs phi; tr2 phi [rad]; efficiency' tr2Phi_numerator       tr2Phi_denominator",
        "effic_tr2Eta       'tr2 efficiency vs eta; tr2 eta [rad]; efficiency' tr2Eta_numerator       tr2Eta_denominator",
        "effic_tr2Pt       'tr2 efficiency vs pt; tr2 pt [GeV]; efficiency' tr2Pt_numerator       tr2Pt_denominator",

    ),
)

#        "effic_phPhi       'ph efficiency vs phi; ph phi [rad]; efficiency' phPhi_numerator       phPhi_denominator",
#        "effic_phEta       'ph efficiency vs eta; ph eta [rad]; efficiency' phEta_numerator       phEta_denominator",
#        "effic_phPt       'ph efficiency vs pt; ph pt [GeV]; efficiency' phPt_numerator       phPt_denominator",
#        "effic_trPhi       'tr efficiency vs phi; tr phi [rad]; efficiency' trPhi_numerator       trPhi_denominator",
#        "effic_trEta       'tr efficiency vs eta; tr eta [rad]; efficiency' trEta_numerator       trEta_denominator",
#        "effic_trPt       'tr efficiency vs pt; tr pt [GeV]; efficiency' trPt_numerator       trPt_denominator",
#        "effic_mu1Phi       'mu1 efficiency vs phi; mu1 phi [rad]; efficiency' mu1Phi_numerator       mu1Phi_denominator",
#        "effic_mu1Eta       'mu1 efficiency vs eta; mu1 eta [rad]; efficiency' mu1Eta_numerator       mu1Eta_denominator",
#        "effic_mu1Pt       'mu1 efficiency vs pt; mu1 pt [GeV]; efficiency' mu1Pt_numerator       mu1Pt_denominator",
#        "effic_mu2Phi       'mu2 efficiency vs phi; mu2 phi [rad]; efficiency' mu2Phi_numerator       mu2Phi_denominator",
#        "effic_mu2Eta       'mu2 efficiency vs eta; mu2 eta [rad]; efficiency' mu2Eta_numerator       mu2Eta_denominator",
#        "effic_mu2Pt       'mu2 efficiency vs pt; mu2 pt [GeV]; efficiency' mu2Pt_numerator       mu2Pt_denominator",
#        "effic_mu3Phi       'mu3 efficiency vs phi; mu3 phi [rad]; efficiency' mu3Phi_numerator       mu3Phi_denominator",
#        "effic_mu3Eta       'mu3 efficiency vs eta; mu3 eta [rad]; efficiency' mu3Eta_numerator       mu3Eta_denominator",
#        "effic_mu3Pt       'mu3 efficiency vs pt; mu3 pt [GeV]; efficiency' mu3Pt_numerator       mu3Pt_denominator",
#        "effic_DiMuPhi       'DiMu efficiency vs phi; DiMu phi [rad]; efficiency' DiMuPhi_numerator       DiMuPhi_denominator",
#        "effic_DiMuEta       'DiMu efficiency vs eta; DiMu eta [rad]; efficiency' DiMuEta_numerator       DiMuEta_denominator",
#        "effic_DiMuPt       'DiMu efficiency vs pt; DiMu pt [GeV]; efficiency' DiMuPt_numerator       DiMuPt_denominator",
#        "effic_DiMuPVcos       'DiMu efficiency vs cosPV; DiMu cosPV ; efficiency' DiMuPVcos_numerator       DiMuPVcos_denominator",
#        "effic_DiMuProb       'DiMu efficiency vs prob; DiMu prob ; efficiency' DiMuProb_numerator       DiMuProb_denominator",
#        "effic_DiMuDS       'DiMu efficiency vs DS; DiMu DS; efficiency' DiMuDS_numerator       DiMuDS_denominator",
#        "effic_DiMuDCA       'DiMu efficiency vs DCA; DiMu DCA [cm]; efficiency' DiMuDCA_numerator       DiMuDCA_denominator",
#        "effic_DiMuMass       'DiMu efficiency vs Mass; DiMu Mass[GeV]; efficiency' DiMuMass_numerator       DiMuMass_denominator",
#        "effic_DiMudR       'DiMu efficiency vs dR; DiMu dR; efficiency' DiMudR_numerator       DiMudR_denominator",
#        "effic_tr_d0       'tr efficiency vs d0; tr d0 [cm]; efficiency' tr_d0_numerator       tr_d0_denominator",
#        "effic_tr_z0       'tr efficiency vs z0; tr z0 [cm]; efficiency' tr_z0_numerator       tr_z0_denominator",


#    ),
#    efficiencyProfile = cms.untracked.vstring(
#        "effic_met_vs_LS 'MET efficiency vs LS; LS; PF MET efficiency' metVsLS_numerator metVsLS_denominator"
#    ),
  
#)

bphClient = cms.Sequence(
    bphEfficiency1
    + bphEfficiency2
    + bphEfficiency2_1
    + bphEfficiency2_2
#    + bphEfficiency2_3
    + bphEfficiency3
    + bphEfficiency3_1
    + bphEfficiency4
    + bphEfficiency4_1
    + bphEfficiency4_2
    + bphEfficiency4_3
    + bphEfficiency4_4
    + bphEfficiency4_5
    + bphEfficiency4_7
    + bphEfficiency4_8
    + bphEfficiency4_9
    + bphEfficiency4_10
    + bphEfficiency4_11
    + bphEfficiency4_12
    + bphEfficiency4_13
    + bphEfficiency4_15
    + bphEfficiency4_16
    + bphEfficiency4_17
    + bphEfficiency5
    + bphEfficiency5_1
    + bphEfficiency5_2
    + bphEfficiency5_3
    + bphEfficiency5_4
    + bphEfficiency5_5
    + bphEfficiency5_6
    + bphEfficiency6
    + bphEfficiency6_1
    + bphEfficiency6_2
##    + bphEfficiency7
    + bphEfficiency8
    + bphEfficiency8_1
    + bphEfficiency8_2
    + bphEfficiency9
    + bphEfficiency10
    + bphEfficiency10_1
    + bphEfficiency10_2
    + bphEfficiency10_3
    + bphEfficiency11
)

##

##
