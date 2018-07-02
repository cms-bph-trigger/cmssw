#include "DQMOffline/Trigger/plugins/BPHMonitor.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"


// -----------------------------
//  constructors and destructor
// -----------------------------

BPHMonitor::BPHMonitor( const edm::ParameterSet& iConfig ) : 
  folderName_ ( iConfig.getParameter<std::string>("FolderName") )
  , muoToken_ ( mayConsume<reco::MuonCollection> (iConfig.getParameter<edm::InputTag>("muons") ) )  
  , bsToken_ ( mayConsume<reco::BeamSpot> (iConfig.getParameter<edm::InputTag>("beamSpot")))
  , trToken_ ( mayConsume<reco::TrackCollection> (iConfig.getParameter<edm::InputTag>("tracks")))
  , phToken_ ( mayConsume<reco::PhotonCollection> (iConfig.getParameter<edm::InputTag>("photons")))
  , vtxToken_( mayConsume<reco::VertexCollection>( iConfig.getParameter<edm::InputTag>("offlinePVs") ) )
  , phi_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("phiPSet") ) )
  , pt_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("ptPSet") ) )
  , dMu_pt_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("dMu_ptPSet") ) )
  , eta_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("etaPSet") ) )
  , d0_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("d0PSet") ) )
  , z0_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("z0PSet") ) )
  , dR_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("dRPSet") ) )
  , mass_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("massPSet") ) )
  , Bmass_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("BmassPSet") ) )
  , dca_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("dcaPSet") ) )
  , ds_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("dsPSet") ) )
  , cos_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("cosPSet") ) )
  , prob_binning_ ( getHistoPSet (iConfig.getParameter<edm::ParameterSet>("histoPSet").getParameter<edm::ParameterSet> ("probPSet") ) )
  , verbosity_( iConfig.getParameter<unsigned int>("verbosityLevel") )
  , num_genTriggerEventFlag_(new GenericTriggerEventFlag(iConfig.getParameter<edm::ParameterSet>("numGenericTriggerEventPSet"),consumesCollector(), *this))
  , den_genTriggerEventFlag_(new GenericTriggerEventFlag(iConfig.getParameter<edm::ParameterSet>("denGenericTriggerEventPSet"),consumesCollector(), *this))
  , hltPrescale_ (new HLTPrescaleProvider(iConfig, consumesCollector(), *this))
  , muPt ( iConfig.getParameter<std::string>("muPt") )
  , muEta ( iConfig.getParameter<std::string>("muEta") )
  , tkPt ( iConfig.getParameter<std::string>("tkPt") )
  , tkEta ( iConfig.getParameter<std::string>("tkEta") )
  , dMuEta ( iConfig.getParameter<std::string>("dMuEta") )
  , dMuPt ( iConfig.getParameter<std::string>("dMuPt") )
  , dMuMass ( iConfig.getParameter<std::string>("dMuMass") )
  , muQual ( iConfig.getParameter<std::string>("muQual") )
  , muoSelection_ ( iConfig.getParameter<std::string>("muoSelection") )
  , muoSelection_ref ( iConfig.getParameter<std::string>("muoSelection_ref") )
  , muoSelection_tag ( iConfig.getParameter<std::string>("muoSelection_tag") )
  , muoSelection_probe ( iConfig.getParameter<std::string>("muoSelection_probe") )
  , nmuons_ ( iConfig.getParameter<int>("nmuons" ) )
  , tnp_ ( iConfig.getParameter<bool>("tnp" ) )
  , L3_ ( iConfig.getParameter<int>("L3" ) )
  , ptCut_ ( iConfig.getParameter<double>("ptCut" ) )
  , displaced_ ( iConfig.getParameter<int>("displaced" ) )
  , trOrMu_ ( iConfig.getParameter<int>("trOrMu" ) )
  , Jpsi_ ( iConfig.getParameter<int>("Jpsi" ) )
  , Upsilon_ ( iConfig.getParameter<int>("Upsilon" ) ) // if ==1 path with Upsilon constraint
  , enum_ ( iConfig.getParameter<int>("enum" ) )
  , seagull_ ( iConfig.getParameter<int>("seagull" ) )
  , maxmass_ ( iConfig.getParameter<double>("maxmass" ) )
  , minmass_ ( iConfig.getParameter<double>("minmass" ) )
  , maxmassJpsi ( iConfig.getParameter<double>("maxmassJpsi" ) )
  , minmassJpsi ( iConfig.getParameter<double>("minmassJpsi" ) )
  , maxmassUpsilon ( iConfig.getParameter<double>("maxmassUpsilon" ) )
  , minmassUpsilon ( iConfig.getParameter<double>("minmassUpsilon" ) )
  , maxmassTkTk ( iConfig.getParameter<double>("maxmassTkTk" ) )
  , minmassTkTk ( iConfig.getParameter<double>("minmassTkTk" ) )
  , maxmassJpsiTk ( iConfig.getParameter<double>("maxmassJpsiTk" ) )
  , minmassJpsiTk ( iConfig.getParameter<double>("minmassJpsiTk" ) )
  , kaon_mass ( iConfig.getParameter<double>("kaon_mass" ) )
  , mu_mass ( iConfig.getParameter<double>("mu_mass" ) )
  , min_dR ( iConfig.getParameter<double>("min_dR" ) )
  , max_dR ( iConfig.getParameter<double>("max_dR" ) )
  , minprob ( iConfig.getParameter<double>("minprob" ) )
  , mincos ( iConfig.getParameter<double>("mincos" ) )
  , minDS ( iConfig.getParameter<double>("minDS" ) )
  , hltInputTag_1 ( iConfig.getParameter<edm::InputTag>("hltTriggerSummaryAOD"))
  , hltInputTag_ (mayConsume<trigger::TriggerEvent>( iConfig.getParameter<edm::InputTag>("hltTriggerSummaryAOD")))
  , hltpaths_num ( iConfig.getParameter<edm::ParameterSet>("numGenericTriggerEventPSet").getParameter<std::vector<std::string>>("hltPaths"))
  , hltpaths_den ( iConfig.getParameter<edm::ParameterSet>("denGenericTriggerEventPSet").getParameter<std::vector<std::string>>("hltPaths"))
  , trSelection_ ( iConfig.getParameter<std::string>("muoSelection") )
  , trSelection_ref ( iConfig.getParameter<std::string>("trSelection_ref") )
  , DMSelection_ref ( iConfig.getParameter<std::string>("DMSelection_ref") )
{

  muPhi_.numerator = nullptr;
  muPhi_.denominator = nullptr;
  muEta_.numerator = nullptr;
  muEta_.denominator = nullptr;
  muPt_.numerator = nullptr;
  muPt_.denominator = nullptr;
  mud0_.numerator = nullptr;
  mud0_.denominator = nullptr;
  muz0_.numerator = nullptr;
  muz0_.denominator = nullptr;

  mu1Phi_.numerator = nullptr;
  mu1Phi_.denominator = nullptr;
  mu1Eta_.numerator = nullptr;
  mu1Eta_.denominator = nullptr;
  mu1Pt_.numerator = nullptr;
  mu1Pt_.denominator = nullptr;

  mu2Phi_.numerator = nullptr;
  mu2Phi_.denominator = nullptr;
  mu2Eta_.numerator = nullptr;
  mu2Eta_.denominator = nullptr;
  mu2Pt_.numerator = nullptr;
  mu2Pt_.denominator = nullptr;

  mu3Phi_.numerator = nullptr;
  mu3Phi_.denominator = nullptr;
  mu3Eta_.numerator = nullptr;
  mu3Eta_.denominator = nullptr;
  mu3Pt_.numerator = nullptr;
  mu3Pt_.denominator = nullptr;

  phPhi_.numerator = nullptr;
  phPhi_.denominator = nullptr;
  phEta_.numerator = nullptr;
  phEta_.denominator = nullptr;
  phPt_.numerator = nullptr;
  phPt_.denominator = nullptr;


  DiMuPhi_.numerator = nullptr;
  DiMuPhi_.denominator = nullptr;
  DiMuEta_.numerator = nullptr;
  DiMuEta_.denominator = nullptr;
  DiMuPt_.numerator = nullptr;
  DiMuPt_.denominator = nullptr;
  DiMuPVcos_.numerator = nullptr;
  DiMuPVcos_.denominator = nullptr;
  DiMuProb_.numerator = nullptr;
  DiMuProb_.denominator = nullptr;
  DiMuDS_.numerator = nullptr;
  DiMuDS_.denominator = nullptr;
  DiMuDCA_.numerator = nullptr;
  DiMuDCA_.denominator = nullptr;
  DiMuMass_.numerator = nullptr;
  DiMuMass_.denominator = nullptr;
  BMass_.numerator = nullptr;
  BMass_.denominator = nullptr;
  DiMudR_.numerator = nullptr;
  DiMudR_.denominator = nullptr;


}

BPHMonitor::~BPHMonitor()
{
  if (num_genTriggerEventFlag_) delete num_genTriggerEventFlag_;
  if (den_genTriggerEventFlag_) delete den_genTriggerEventFlag_;
  delete hltPrescale_;
}

MEbinning BPHMonitor::getHistoPSet(edm::ParameterSet pset)
{
  // Due to the setup of the fillDescription only one of the
  // two cases is possible at this point.
  if (pset.existsAs<std::vector<double>>("edges")) {
    return MEbinning{pset.getParameter<std::vector<double>>("edges")};
  }

  return MEbinning {
    pset.getParameter<int32_t>("nbins"),
      pset.getParameter<double>("xmin"),
      pset.getParameter<double>("xmax"),
      };
}

// MEbinning BPHMonitor::getHistoLSPSet(edm::ParameterSet pset)
// {
//   return MEbinning{
//     pset.getParameter<int32_t>("nbins"),
//       0.,
//       double(pset.getParameter<int32_t>("nbins"))
//       };
// }

void BPHMonitor::setMETitle(METME& me, std::string titleX, std::string titleY)
{
  me.numerator->setAxisTitle(titleX,1);
  me.numerator->setAxisTitle(titleY,2);
  me.denominator->setAxisTitle(titleX,1);
  me.denominator->setAxisTitle(titleY,2);

}

void BPHMonitor::bookME(DQMStore::IBooker &ibooker, METME& me, std::string& histname, std::string& histtitle, int& nbins, double& min, double& max)
{
  me.numerator   = ibooker.book1D(histname+"_numerator",   histtitle+" (numerator)",   nbins, min, max);
  me.denominator = ibooker.book1D(histname+"_denominator", histtitle+" (denominator)", nbins, min, max);
}
void BPHMonitor::bookME(DQMStore::IBooker &ibooker, METME& me, std::string& histname, std::string& histtitle, std::vector<double> binning)
{
  int nbins = binning.size()-1;
  std::vector<float> fbinning(binning.begin(),binning.end());
  float* arr = &fbinning[0];
  me.numerator   = ibooker.book1D(histname+"_numerator",   histtitle+" (numerator)",   nbins, arr);
  me.denominator = ibooker.book1D(histname+"_denominator", histtitle+" (denominator)", nbins, arr);
}
void BPHMonitor::bookME(DQMStore::IBooker &ibooker, METME& me, std::string& histname, std::string& histtitle, int& nbinsX, double& xmin, double& xmax, double& ymin, double& ymax)
{
  me.numerator   = ibooker.bookProfile(histname+"_numerator",   histtitle+" (numerator)",   nbinsX, xmin, xmax, ymin, ymax);
  me.denominator = ibooker.bookProfile(histname+"_denominator", histtitle+" (denominator)", nbinsX, xmin, xmax, ymin, ymax);
}
void BPHMonitor::bookME(DQMStore::IBooker &ibooker, METME& me, std::string& histname, std::string& histtitle, int& nbinsX, double& xmin, double& xmax, int& nbinsY, double& ymin, double& ymax)
{
  me.numerator   = ibooker.book2D(histname+"_numerator",   histtitle+" (numerator)",   nbinsX, xmin, xmax, nbinsY, ymin, ymax);
  me.denominator = ibooker.book2D(histname+"_denominator", histtitle+" (denominator)", nbinsX, xmin, xmax, nbinsY, ymin, ymax);
}
void BPHMonitor::bookME(DQMStore::IBooker &ibooker, METME& me, std::string& histname, std::string& histtitle, std::vector<double> binningX, std::vector<double> binningY)
{
  int nbinsX = binningX.size()-1;
  std::vector<float> fbinningX(binningX.begin(),binningX.end());
  float* arrX = &fbinningX[0];
  int nbinsY = binningY.size()-1;
  std::vector<float> fbinningY(binningY.begin(),binningY.end());
  float* arrY = &fbinningY[0];

  me.numerator   = ibooker.book2D(histname+"_numerator",   histtitle+" (numerator)",   nbinsX, arrX, nbinsY, arrY);
  me.denominator = ibooker.book2D(histname+"_denominator", histtitle+" (denominator)", nbinsX, arrX, nbinsY, arrY);
}

void BPHMonitor::bookME(DQMStore::IBooker &ibooker, METME &me, std::string &histname, std::string &histtitle, /*const*/ MEbinning& binning)
{
  // If the vector in the binning is filled use the bins defined there
  // otherwise use a linear binning between min and max
  if (binning.edges.empty()) {
    this->bookME(ibooker, me, histname, histtitle, binning.nbins, binning.xmin, binning.xmax);
  } else {
    this->bookME(ibooker, me, histname, histtitle, binning.edges);
  }
}


void BPHMonitor::bookHistograms(DQMStore::IBooker     & ibooker,
				edm::Run const        & iRun,
				edm::EventSetup const & iSetup) 
{  
  std::string histname, histtitle, istnp, trMuPh;
  bool Ph_ = false; if (enum_ == 7) Ph_ = true;
  if (tnp_) istnp = "Tag_and_Probe/"; else istnp = "";
  std::string currentFolder = folderName_ + istnp;
  ibooker.setCurrentFolder(currentFolder);
  if (trOrMu_) trMuPh = "tr"; else if (Ph_) trMuPh = "ph"; else trMuPh = "mu";

  if (enum_ == 7 || enum_ == 1 || enum_ == 9 || enum_ == 10) {  
    histname = trMuPh+"Pt"; histtitle = trMuPh+"_P_{t}";
    bookME(ibooker,muPt_,histname,histtitle, pt_binning_);
    setMETitle(muPt_,trMuPh+"_Pt[GeV]","events / 1 GeV");

    histname = trMuPh+"Phi"; histtitle = trMuPh+"Phi";
    bookME(ibooker,muPhi_,histname,histtitle, phi_binning_);
    setMETitle(muPhi_,trMuPh+"_#phi","events / 0.1 rad");

    histname = trMuPh+"Eta"; histtitle = trMuPh+"_Eta";
    bookME(ibooker,muEta_,histname,histtitle, eta_binning_);
    setMETitle(muEta_,trMuPh+"_#eta","events / 0.2");
    
    if (enum_ ==9)
    {
      histname = "BMass"; histtitle = "BMass";
      bookME(ibooker,BMass_,histname,histtitle, Bmass_binning_);
      setMETitle(BMass_,"B_#mass","events /");

    }
  }
  else {
    if (enum_ !=8)
    {
    histname = trMuPh+"1Pt"; histtitle = trMuPh+"1_P_{t}";
    bookME(ibooker,mu1Pt_,histname,histtitle, pt_binning_);
    setMETitle(mu1Pt_,trMuPh+"_Pt[GeV]","events / 1 GeV");

    histname = trMuPh+"1Phi"; histtitle = trMuPh+"1Phi";
    bookME(ibooker,mu1Phi_,histname,histtitle, phi_binning_);
    setMETitle(mu1Phi_,trMuPh+"_#phi","events / 0.1 rad");
  
    histname = trMuPh+"1Eta"; histtitle = trMuPh+"1_Eta";
    bookME(ibooker,mu1Eta_,histname,histtitle, eta_binning_);
    setMETitle(mu1Eta_,trMuPh+"_#eta","events / 0.2");

    histname = trMuPh+"2Pt"; histtitle = trMuPh+"2_P_{t}";
    bookME(ibooker,mu2Pt_,histname,histtitle, pt_binning_);
    setMETitle(mu2Pt_,trMuPh+"_Pt[GeV]","events / 1 GeV");

    histname = trMuPh+"2Phi"; histtitle = trMuPh+"2Phi";
    bookME(ibooker,mu2Phi_,histname,histtitle, phi_binning_);
    setMETitle(mu2Phi_,trMuPh+"_#phi","events / 0.1 rad");

    histname = trMuPh+"2Eta"; histtitle = trMuPh+"2_Eta";
    bookME(ibooker,mu2Eta_,histname,histtitle, eta_binning_);
    setMETitle(mu2Eta_,trMuPh+"_#eta","events / 0.2");
    if (enum_ ==11)
    {
      histname = "BMass"; histtitle = "BMass";
      bookME(ibooker,BMass_,histname,histtitle, Bmass_binning_);
      setMETitle(BMass_,"B_#mass","events /");

    }

    }
    if (enum_ == 6) {
      histname = trMuPh+"3Eta"; histtitle = trMuPh+"3Eta";
      bookME(ibooker,mu3Eta_,histname,histtitle, eta_binning_);
      setMETitle(mu3Eta_,trMuPh+"3#eta","events / 0.2");

      histname = trMuPh+"3Pt"; histtitle = trMuPh+"3_P_{t}";
      bookME(ibooker,mu3Pt_,histname,histtitle, pt_binning_);
      setMETitle(mu3Pt_,trMuPh+"3_Pt[GeV]","events / 1 GeV");

      histname = trMuPh+"3Phi"; histtitle = trMuPh+"3Phi";
      bookME(ibooker,mu3Phi_,histname,histtitle, phi_binning_);
      setMETitle(mu3Phi_,trMuPh+"3_#phi","events / 0.1 rad");

    }
    else if (enum_ == 2 || enum_ == 4 || enum_ == 5 || enum_ == 8) {
      histname = "DiMuEta"; histtitle = "DiMuEta";
      bookME(ibooker,DiMuEta_,histname,histtitle, eta_binning_);
      setMETitle(DiMuEta_,"DiMu#eta","events / 0.2");

      histname = "DiMuPt"; histtitle = "DiMu_P_{t}";
      bookME(ibooker,DiMuPt_,histname,histtitle, dMu_pt_binning_);
      setMETitle(DiMuPt_,"DiMu_Pt[GeV]","events / 1 GeV");

      histname = "DiMuPhi"; histtitle = "DiMuPhi";
      bookME(ibooker,DiMuPhi_,histname,histtitle, phi_binning_);
      setMETitle(DiMuPhi_,"DiMu_#phi","events / 0.1 rad");

      if (enum_ == 4 || enum_ == 5) {
	histname = "DiMudR"; histtitle = "DiMudR";
	bookME(ibooker,DiMudR_,histname,histtitle, dR_binning_);
	setMETitle(DiMudR_,"DiMu_#dR","events /");

	if (enum_ == 4) {
	  histname = "DiMuMass"; histtitle = "DiMuMass";
	  bookME(ibooker,DiMuMass_,histname,histtitle, mass_binning_);
	  setMETitle(DiMuMass_,"DiMu_#mass","events /");

	}
      } else if (enum_ == 8) {
	histname = "DiMuProb"; histtitle = "DiMuProb";
	bookME(ibooker,DiMuProb_,histname,histtitle, prob_binning_);
	setMETitle(DiMuProb_,"DiMu_#prob","events /");

	histname = "DiMuPVcos"; histtitle = "DiMuPVcos";
	bookME(ibooker,DiMuPVcos_,histname,histtitle, cos_binning_);
	setMETitle(DiMuPVcos_,"DiMu_#cosPV","events /");

	histname = "DiMuDS"; histtitle = "DiMuDS";
	bookME(ibooker,DiMuDS_,histname,histtitle, ds_binning_);
	setMETitle(DiMuDS_,"DiMu_#ds","events /");

	histname = "DiMuDCA"; histtitle = "DiMuDCA";
	bookME(ibooker,DiMuDCA_,histname,histtitle, dca_binning_);
	setMETitle(DiMuDCA_,"DiMu_#dca","events /");

      }
    } // if (enum_ == 2 || enum_ == 4 || enum_ == 5 || enum_ == 8)
  }

  // Initialize the GenericTriggerEventFlag
  if ( num_genTriggerEventFlag_ && num_genTriggerEventFlag_->on() ) num_genTriggerEventFlag_->initRun( iRun, iSetup );
  if ( den_genTriggerEventFlag_ && den_genTriggerEventFlag_->on() ) den_genTriggerEventFlag_->initRun( iRun, iSetup );
  bool changed = true;

  if(!hltPrescale_->init(iRun,iSetup,"HLT",changed) )
  {
    // std::cout<<"prescale module init failed"<<std::endl;
    // return;
  }
  // else std::cout<<"prescale module init succeded"<<std::endl;

  hltConfig_ = hltPrescale_->hltConfigProvider();
}

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "TLorentzVector.h"
void BPHMonitor::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup)  {

  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken( bsToken_,  beamSpot);

  edm::Handle<reco::MuonCollection> muoHandle;
  iEvent.getByToken( muoToken_, muoHandle );


  edm::Handle<reco::TrackCollection> trHandle;
  iEvent.getByToken( trToken_, trHandle );

  edm::Handle<reco::PhotonCollection> phHandle;
  iEvent.getByToken( phToken_, phHandle );


  edm::Handle<edm::TriggerResults> handleTriggerTrigRes; 

  edm::ESHandle<MagneticField> bFieldHandle;

 

  int PrescaleWeight =1;
  const std::string & hltpath = getTriggerName(hltpaths_den[0]);
  const std::string & hltpath1 = getTriggerName(hltpaths_num[0]);



  if (den_genTriggerEventFlag_->on() &&  den_genTriggerEventFlag_->accept( iEvent, iSetup) && num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) ) 
  {
    PrescaleWeight = Prescale(hltpath1, hltpath, iEvent, iSetup, hltPrescale_);

  }
  if (tnp_>0) {//TnP method 
    if (den_genTriggerEventFlag_->on() && ! den_genTriggerEventFlag_->accept( iEvent, iSetup) ) return;
    iEvent.getByToken( hltInputTag_, handleTriggerEvent);
    if (handleTriggerEvent->sizeFilters()== 0) return;
    std::vector<reco::Muon> tagMuons;
    for ( auto const & m : *muoHandle ) { // applying tag selection 
      if ( !matchToTrigger(hltpath,m))continue;
      if ( muoSelection_ref( m ) ) tagMuons.push_back(m);
    }
    for (int i = 0; i<int(tagMuons.size());i++) {
      for ( auto const & m : *muoHandle ) { 
        if ( !matchToTrigger(hltpath,m))continue;
        if ((tagMuons[i].pt() == m.pt())) continue; //not the same  
        if ((tagMuons[i].p4()+m.p4()).M() >minmass_&& (tagMuons[i].p4()+m.p4()).M() <maxmass_) { //near to J/psi mass
          muPhi_.denominator->Fill(m.phi());
          muEta_.denominator->Fill(m.eta());
          muPt_.denominator ->Fill(m.pt());
          if (muoSelection_( m ) && num_genTriggerEventFlag_->on() && num_genTriggerEventFlag_->accept( iEvent, iSetup)) {
            muPhi_.numerator->Fill(m.phi(),PrescaleWeight);
            muEta_.numerator->Fill(m.eta(),PrescaleWeight);
            muPt_.numerator ->Fill(m.pt(),PrescaleWeight);
          }
        }      
      }
 
    }
 

  }  
  else { // reference method
    if (den_genTriggerEventFlag_->on() && ! den_genTriggerEventFlag_->accept( iEvent, iSetup) ) return;
    iEvent.getByToken( hltInputTag_, handleTriggerEvent);
    if (handleTriggerEvent->sizeFilters()== 0) return;
    for (auto const & m : *muoHandle ) {
      if ( !matchToTrigger(hltpath,m))continue;
      if (!muQual(m)) continue; 
      for (auto const & m1 : *muoHandle ) {
      if (!muQual(m1)) continue; 
        if (!(m1.pt() > m.pt())) continue; //to get rid of double counting
	iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);
	const reco::BeamSpot& vertexBeamSpot = *beamSpot;
	std::vector<reco::TransientTrack> j_tks;
	reco::TransientTrack mu1TT(m.track(), &(*bFieldHandle));
	reco::TransientTrack mu2TT(m1.track(), &(*bFieldHandle));
	j_tks.push_back(mu1TT);
	j_tks.push_back(mu2TT);
	KalmanVertexFitter jkvf;
	TransientVertex jtv = jkvf.vertex(j_tks);
	if (!jtv.isValid()) continue;
	reco::Vertex jpsivertex = jtv;
	float dimuonCL = 0;
	if ( (jpsivertex.chi2() >= 0) && (jpsivertex.ndof() > 0) ) // I think these values are "unphysical"(no one will need to change them ever)so the can be fixed
	  dimuonCL = TMath::Prob(jpsivertex.chi2(), jpsivertex.ndof() );
	math::XYZVector jpperp(m.px() + m1.px() ,
			       m.py() + m1.py() ,
			       0.);
	GlobalPoint jVertex = jtv.position();
	GlobalError jerr = jtv.positionError();
	GlobalPoint displacementFromBeamspotJpsi( -1*((vertexBeamSpot.x0() - jVertex.x()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()),
						  -1*((vertexBeamSpot.y0() - jVertex.y()) + (jVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()),
						  0);
        reco::Vertex::Point vperpj(displacementFromBeamspotJpsi.x(), displacementFromBeamspotJpsi.y(), 0.);
        float jpsi_cos = vperpj.Dot(jpperp) / (vperpj.R()*jpperp.R());
        TrajectoryStateClosestToPoint mu1TS = mu1TT.impactPointTSCP();
        TrajectoryStateClosestToPoint mu2TS = mu2TT.impactPointTSCP();
        ClosestApproachInRPhi cApp;
        if (mu1TS.isValid() && mu2TS.isValid()) 
        {
          if (!cApp.calculate(mu1TS.theState(), mu2TS.theState()))continue;
        } 
        else continue;

        double DiMuMass = (m1.p4()+m.p4()).M();
        bool detach=true;
        switch(enum_) { // enum_ = 1...9, represents different sets of variables for different paths, we want to have different hists for different paths

        case 1: tnp_=true; // already filled hists for tnp method

        case 2:
          if ((Jpsi_) && (!Upsilon_)) {
            if (DiMuMass> maxmassJpsi || DiMuMass< minmassJpsi) continue;
          }
          if ((!Jpsi_) && (Upsilon_)) {
            if (DiMuMass> maxmassUpsilon || DiMuMass< minmassUpsilon) continue;
          }
          if (dimuonCL < minprob)continue;
          if ( muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) ) 
          {
            mu1Pt_.denominator ->Fill(m.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Pt_.numerator ->Fill(m.pt(),PrescaleWeight);
            }
          }
          if ( muPt(m) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu2Pt_.denominator ->Fill(m1.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Pt_.numerator ->Fill(m1.pt(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu1Eta_.denominator->Fill(m.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Eta_.numerator ->Fill(m.eta(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && muEta(m) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu2Eta_.denominator->Fill(m1.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Eta_.numerator ->Fill(m1.eta(),PrescaleWeight);
            }

          } 
          if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && dMuEta(m1.p4()+m.p4()) )
          {
            DiMuPt_.denominator ->Fill((m1.p4()+m.p4()).Pt() );
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuPt_.numerator ->Fill((m1.p4()+m.p4()).Pt(),PrescaleWeight);
            }
          } 
          if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4())  )
          {
            DiMuEta_.denominator ->Fill((m1.p4()+m.p4()).Eta() );
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuEta_.numerator ->Fill((m1.p4()+m.p4()).Eta(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu1Phi_.denominator->Fill(m.phi());
            mu2Phi_.denominator->Fill(m1.phi());
            DiMuPhi_.denominator ->Fill((m1.p4()+m.p4()).Phi());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) 
              {
                mu1Phi_.numerator->Fill(m.phi(),PrescaleWeight);
                mu2Phi_.numerator->Fill(m1.phi(),PrescaleWeight);
                DiMuPhi_.numerator ->Fill((m1.p4()+m.p4()).Phi(),PrescaleWeight);
              }
            }
          }
          
          break;

        case 3:
          if ((Jpsi_) && (!Upsilon_)){
            if (DiMuMass> maxmassJpsi || DiMuMass< minmassJpsi) continue;
          }
          if ((!Jpsi_) && (Upsilon_)){
            if (DiMuMass> maxmassUpsilon || DiMuMass< minmassUpsilon) continue;
          }
          if (dimuonCL<minprob) continue;
          if ( muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu1Pt_.denominator ->Fill(m.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Pt_.numerator ->Fill(m.pt(),PrescaleWeight);
            }
          }
          if ( muPt(m) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu2Pt_.denominator ->Fill(m1.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Pt_.numerator ->Fill(m1.pt(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu1Eta_.denominator->Fill(m.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Eta_.numerator ->Fill(m.eta(),PrescaleWeight);
            }
          }
          if ( muPt(m) && muPt(m1) && muEta(m) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu2Eta_.denominator->Fill(m1.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Eta_.numerator ->Fill(m1.eta(),PrescaleWeight);
            }
          }
          if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu1Phi_.denominator->Fill(m.phi());
            mu2Phi_.denominator->Fill(m1.phi());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m))
              {
                mu1Phi_.numerator->Fill(m.phi(),PrescaleWeight);
                mu2Phi_.numerator->Fill(m1.phi(),PrescaleWeight);
              }

            }
          }
          
          break; 

        case 4:
          if (dimuonCL<minprob) continue;
          if ( muPt(m) && muPt(m1) && (m1.pt()>ptCut_) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            DiMuMass_.denominator ->Fill(DiMuMass);
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup))
      	    {
      	      if (seagull_ && ((m.charge()* deltaPhi(m.phi(), m1.phi())) > 0.) ) continue;
      	      if( !matchToTrigger(hltpath1,m1))continue;          
      	      if( !matchToTrigger(hltpath1,m))continue;          
      	      DiMuMass_.numerator ->Fill(DiMuMass, PrescaleWeight);      
      	    }
          }
          if ((Jpsi_) && (!Upsilon_)){
            if (DiMuMass> maxmassJpsi || DiMuMass< minmassJpsi) continue;
          }
          if ((!Jpsi_) && (Upsilon_)){
            if (DiMuMass> maxmassUpsilon || DiMuMass< minmassUpsilon) continue;
          }

          if ( muPt(m1) && (m1.pt()>ptCut_) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) && dMuMass(m1.p4()+m.p4())) 
          {
            mu1Pt_.denominator ->Fill(m.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
    	        if (seagull_ && ((m.charge()* deltaPhi(m.phi(), m1.phi())) > 0.) ) continue;
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Pt_.numerator ->Fill(m.pt(),PrescaleWeight);
            }
          }
          if ( muPt(m) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) && dMuMass(m1.p4()+m.p4()))
          {
            mu2Pt_.denominator ->Fill(m1.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Pt_.numerator ->Fill(m1.pt(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && (m1.pt()>ptCut_) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) && dMuMass(m1.p4()+m.p4()))
          {
            mu1Eta_.denominator->Fill(m.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Eta_.numerator ->Fill(m.eta(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && (m1.pt()>ptCut_) && muEta(m) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) && dMuMass(m1.p4()+m.p4()))
          {
            mu2Eta_.denominator->Fill(m1.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Eta_.numerator ->Fill(m1.eta(),PrescaleWeight);
            }

          } 
          if ( muPt(m) && muPt(m1) && (m1.pt()>ptCut_) && muEta(m) && muEta(m1) && dMuEta(m1.p4()+m.p4()) && dMuMass(m1.p4()+m.p4()))
          {
            DiMuPt_.denominator ->Fill((m1.p4()+m.p4()).Pt() );
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuPt_.numerator ->Fill((m1.p4()+m.p4()).Pt(),PrescaleWeight);
            }
          } 
          if ( muPt(m) && muPt(m1) && (m1.pt()>ptCut_) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuMass(m1.p4()+m.p4()) )
          {
            DiMuEta_.denominator ->Fill((m1.p4()+m.p4()).Eta() );
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuEta_.numerator ->Fill((m1.p4()+m.p4()).Eta(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && (m1.pt()>ptCut_) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) && dMuMass(m1.p4()+m.p4()) )
          {
            mu1Phi_.denominator->Fill(m.phi());
            mu2Phi_.denominator->Fill(m1.phi());
            DiMuPhi_.denominator ->Fill((m1.p4()+m.p4()).Phi());
            DiMudR_.denominator ->Fill(reco::deltaR(m,m1));
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) 
              {
                mu1Phi_.numerator->Fill(m.phi(),PrescaleWeight);
                mu2Phi_.numerator->Fill(m1.phi(),PrescaleWeight);
                DiMuPhi_.numerator ->Fill((m1.p4()+m.p4()).Phi(),PrescaleWeight);
    	          DiMudR_.numerator ->Fill(reco::deltaR(m,m1),PrescaleWeight);
              }

            }


          }
          
          break;

        case 5:
          if (dimuonCL<minprob) continue;
          if ((Jpsi_) && (!Upsilon_)){
            if (DiMuMass> maxmassJpsi || DiMuMass< minmassJpsi) continue;
          }
  
          if ((!Jpsi_) && (Upsilon_)){
            if (DiMuMass> maxmassUpsilon || DiMuMass< minmassUpsilon) continue;
          }
          
          if ( muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) ) 
          {
            mu1Pt_.denominator ->Fill(m.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
    	        if (seagull_ && ((m.charge()* deltaPhi(m.phi(), m1.phi())) > 0.) ) continue;
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Pt_.numerator ->Fill(m.pt(),PrescaleWeight);
            }
          }
          if ( muPt(m) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu2Pt_.denominator ->Fill(m1.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Pt_.numerator ->Fill(m1.pt(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu1Eta_.denominator->Fill(m.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu1Eta_.numerator ->Fill(m.eta(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && muEta(m) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu2Eta_.denominator->Fill(m1.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) mu2Eta_.numerator ->Fill(m1.eta(),PrescaleWeight);
            }

          } 
          if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && dMuEta(m1.p4()+m.p4()) )
          {
            DiMuPt_.denominator ->Fill((m1.p4()+m.p4()).Pt() );
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuPt_.numerator ->Fill((m1.p4()+m.p4()).Pt(),PrescaleWeight);
            }
          } 
          if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4())  )
          {
            DiMuEta_.denominator ->Fill((m1.p4()+m.p4()).Eta() );
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuEta_.numerator ->Fill((m1.p4()+m.p4()).Eta(),PrescaleWeight);
            }

          }
          if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
          {
            mu1Phi_.denominator->Fill(m.phi());
            mu2Phi_.denominator->Fill(m1.phi());
            DiMuPhi_.denominator ->Fill((m1.p4()+m.p4()).Phi());
            DiMudR_.denominator ->Fill(reco::deltaR(m,m1));
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) 
              {
                mu1Phi_.numerator->Fill(m.phi(),PrescaleWeight);
                mu2Phi_.numerator->Fill(m1.phi(),PrescaleWeight);
                DiMuPhi_.numerator ->Fill((m1.p4()+m.p4()).Phi(),PrescaleWeight);
    	          DiMudR_.numerator ->Fill(reco::deltaR(m,m1),PrescaleWeight);
              }

            }


          }
          break;

        case 6: 
          if (dimuonCL<minprob) continue;
          if ((Jpsi_) && (!Upsilon_)){
            if (DiMuMass> maxmassJpsi || DiMuMass< minmassJpsi) continue;
          }
          if ((!Jpsi_) && (Upsilon_)){
            if (DiMuMass> maxmassUpsilon || DiMuMass< minmassUpsilon) continue;
          }
          for (auto const & m2 : *muoHandle) 
          {//triple muon paths
            if( !matchToTrigger(hltpath,m2))continue;
            if (m2.pt() == m.pt()) continue;

            if ( muPt(m1) && muEta(m) && muEta(m1) && muPt(m2) && muEta(m2)  && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
            {
              mu1Pt_.denominator ->Fill(m.pt());
              if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
              {
                if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,m2)) mu1Pt_.numerator ->Fill(m.pt(),PrescaleWeight);
              }
            }
            if ( muPt(m) && muEta(m) && muEta(m1) && muPt(m2) && muEta(m2) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
            {
              mu2Pt_.denominator ->Fill(m1.pt());
              if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
              {
                if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,m2)) mu2Pt_.numerator ->Fill(m1.pt(),PrescaleWeight);
              }
  
            }
            if ( muPt(m) && muPt(m1) && muEta(m1) && muPt(m2) && muEta(m2) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
            {
              mu1Eta_.denominator->Fill(m.eta());
              if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
              {
                if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,m2)) mu1Eta_.numerator ->Fill(m.eta(),PrescaleWeight);
              }
            }
            if ( muPt(m) && muPt(m1) && muEta(m) && muPt(m2) && muEta(m2) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
            {
              mu2Eta_.denominator->Fill(m1.eta());
              if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
              {
                if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,m2)) mu2Eta_.numerator ->Fill(m1.eta(),PrescaleWeight);
              }
            }
   
            if ( muPt(m) && muEta(m) && muPt(m1) && muEta(m1) && muEta(m2) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
            {
              mu3Pt_.denominator ->Fill(m2.pt());
              if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
              {
                if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,m2)) mu3Pt_.numerator ->Fill(m2.pt(),PrescaleWeight);
              }
  
            }
   
            if ( muPt(m) &&  muEta(m1) && muPt(m1) && muEta(m1) && muPt(m2) &&  dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
            {
              mu3Eta_.denominator->Fill(m2.eta());
              if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
              {
                if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,m2)) mu3Eta_.numerator ->Fill(m2.eta(),PrescaleWeight);
              }
            }
  
            if ( muPt(m) && muPt(m1) && muEta(m) && muEta(m1) && muPt(m2) && muEta(m2) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
            {
              mu1Phi_.denominator->Fill(m.phi());
              mu2Phi_.denominator->Fill(m1.phi());
              mu3Phi_.denominator->Fill(m2.phi());
              if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
              {
                if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,m2))
                {
                  mu1Phi_.numerator->Fill(m.phi(),PrescaleWeight);
                  mu2Phi_.numerator->Fill(m1.phi(),PrescaleWeight);
                  mu3Phi_.numerator->Fill(m2.phi(),PrescaleWeight);
                }
  
              }
            }
            
          } 
          break;    
  
        case 7:// the hists for photon monitoring will be filled on 515 line
          if(!(muoSelection_ref( m1 ) &&  muoSelection_ref( m )))continue;
          if(phHandle->size()>0)
    	    { 
    	      for (auto const & p : *phHandle)
        		{
        		  if( !matchToTrigger(hltpath,p))continue;
        		  phPhi_.denominator->Fill(p.phi());
        		  phEta_.denominator->Fill(p.eta());
        		  phPt_.denominator ->Fill(p.pt());
        		  if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
        		  {
        		      if( !matchToTrigger(hltpath1,p))continue;
        		      if( !matchToTrigger(hltpath1,m))continue;
        		      if( !matchToTrigger(hltpath1,m1))continue;
        		      phPhi_.numerator->Fill(p.phi(),PrescaleWeight);
        		      phEta_.numerator->Fill(p.eta(),PrescaleWeight);
        		      phPt_.numerator ->Fill(p.pt(),PrescaleWeight);
                  
        		  }
              
        		}
    	    } 
          break;
          
        case 8://vtx monitoring, filling probability, DS, DCA, cos of pointing angle to the PV, eta, pT of dimuon
              if(!(muoSelection_ref( m1 ) &&  muoSelection_ref( m )))continue;

              if ((Jpsi_) && (!Upsilon_)){
                if (DiMuMass> maxmassJpsi || DiMuMass< minmassJpsi) continue;
              }
      
           	  if ((!Jpsi_) && (Upsilon_)) {
           	    if (DiMuMass> maxmassUpsilon || DiMuMass< minmassUpsilon) continue;
           	  }
      
              if (displaced_)
              {
                if (displacementFromBeamspotJpsi.perp()/sqrt(jerr.rerr(displacementFromBeamspotJpsi))>minDS)detach = true;
                else detach=false;
      
              }
              else detach = true;

              if ( jpsi_cos>mincos && dimuonCL>minprob && (detach) && dMuEta(m1.p4()+m.p4()) )
              {
                DiMuPt_.denominator ->Fill((m1.p4()+m.p4()).Pt() );
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuPt_.numerator ->Fill((m1.p4()+m.p4()).Pt(),PrescaleWeight);
                }
              } 

              if ( jpsi_cos>mincos && dimuonCL>minprob && (detach) && dMuPt(m1.p4()+m.p4())  )
              {
                DiMuEta_.denominator ->Fill((m1.p4()+m.p4()).Eta() );
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) DiMuEta_.numerator ->Fill((m1.p4()+m.p4()).Eta(),PrescaleWeight);
                }
    
              }
 
              if ( jpsi_cos>mincos && (detach) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
              { 
                DiMuProb_.denominator ->Fill( dimuonCL);
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if( !matchToTrigger(hltpath1,m1))continue;
                  if( !matchToTrigger(hltpath1,m))continue;
                  DiMuProb_.numerator ->Fill( dimuonCL,PrescaleWeight);
                }
              }
              if (  jpsi_cos>mincos && dimuonCL>minprob && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
              {
                DiMuDS_.denominator ->Fill( displacementFromBeamspotJpsi.perp()/sqrt(jerr.rerr(displacementFromBeamspotJpsi)));
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if( !matchToTrigger(hltpath1,m1))continue;
                  if( !matchToTrigger(hltpath1,m))continue;
                  DiMuDS_.denominator ->Fill( displacementFromBeamspotJpsi.perp()/sqrt(jerr.rerr(displacementFromBeamspotJpsi)),PrescaleWeight);
 
                }
              }
              if ( dimuonCL>minprob && (detach) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
              {
                DiMuPVcos_.denominator ->Fill(jpsi_cos );
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m))
                  {
                    DiMuPVcos_.numerator ->Fill(jpsi_cos ,PrescaleWeight);
                  }
                }
              }

 
              if ( jpsi_cos>mincos && dimuonCL>minprob && (detach) && dMuPt(m1.p4()+m.p4()) && dMuEta(m1.p4()+m.p4()) )
              {
                DiMuPhi_.denominator ->Fill((m1.p4()+m.p4()).Phi());
                DiMuDCA_.denominator ->Fill( cApp.distance());
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m)) 
                  {
                    DiMuPhi_.numerator ->Fill((m1.p4()+m.p4()).Phi(),PrescaleWeight);
            	      DiMuDCA_.numerator ->Fill( cApp.distance(),PrescaleWeight);
                  } 
                }
              }
              
              break;

      	case 9:
          if(!(muoSelection_ref( m1 ) &&  muoSelection_ref( m )))continue;
          if (dimuonCL<minprob) continue;
          if (jpsi_cos<mincos) continue;
          if ((displacementFromBeamspotJpsi.perp()/sqrt(jerr.rerr(displacementFromBeamspotJpsi)))<minDS) continue;
      	  if (trHandle.isValid())
    	    { 
            for (auto const & t : *trHandle) 
            {
          		if(!trSelection_ref(t)) continue;
          		const reco::Track& itrk1       = t ;                                                
          		if((reco::deltaR(t,m1) <= min_dR)) continue; // checking overlapping
          		if((reco::deltaR(t,m) <= min_dR)) continue;
          		if (! itrk1.quality(reco::TrackBase::highPurity))     continue;
          		reco::Particle::LorentzVector pB, p1, p2, p3;
          		double trackMass2 = kaon_mass * kaon_mass;
          		double MuMass2 = mu_mass * mu_mass;//0.1056583745 *0.1056583745;
          		double e1   = sqrt(m.momentum().Mag2()  + MuMass2          );
          		double e2   = sqrt(m1.momentum().Mag2()  + MuMass2          );
          		double e3   = sqrt(itrk1.momentum().Mag2() + trackMass2  );
          		p1   = reco::Particle::LorentzVector(m.px() , m.py() , m.pz() , e1  );
          		p2   = reco::Particle::LorentzVector(m1.px() , m1.py() , m1.pz() , e2  );
          		p3   = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3  );
          		pB   = p1 + p2 + p3;
          		reco::TransientTrack trTT(itrk1, &(*bFieldHandle));
          		std::vector<reco::TransientTrack> t_tks;
          		t_tks.push_back(mu1TT);
          		t_tks.push_back(mu2TT);
          		t_tks.push_back(trTT);
          		KalmanVertexFitter kvf;
          		TransientVertex tv  = kvf.vertex(t_tks);
          		reco::Vertex vertex = tv;
          		if (!tv.isValid()) continue;
          		float JpsiTkCL = 0;
          		if ((vertex.chi2()>=0.0) && (vertex.ndof()>0) )   
          		  JpsiTkCL = TMath::Prob(vertex.chi2(), vertex.ndof() );
          		math::XYZVector pperp(m.px() + m1.px() + itrk1.px(),
          				      m.py() + m1.py() + itrk1.py(),
          				      0.);
          		GlobalPoint secondaryVertex = tv.position();
          		GlobalError err             = tv.positionError();
          		GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0() - secondaryVertex.x()) + 
          							  (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), 
          						      -1*((vertexBeamSpot.y0() - secondaryVertex.y()) + 
          							  (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()), 
          						      0);
          		reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
          		float jpsiKcos = vperp.Dot(pperp)/(vperp.R()*pperp.R());
          		if (JpsiTkCL<minprob) continue;
          		if (jpsiKcos<mincos) continue;
          		if ((displacementFromBeamspot.perp()/sqrt(err.rerr(displacementFromBeamspot)))<minDS) continue;
              if (tkPt(t) && tkEta(t))
              {
            		BMass_.denominator ->Fill(pB.mass());
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                    if( !matchToTrigger(hltpath1,m1))continue;
                    if( !matchToTrigger(hltpath1,m))continue;
                    if( !matchToTrigger(hltpath1,t))continue;
                    BMass_.numerator ->Fill(pB.mass(),PrescaleWeight);
                }
              }
          		if( pB.mass()> maxmassJpsiTk || pB.mass()< minmassJpsiTk) continue;
              if (tkEta(t))
              {
                muPt_.denominator ->Fill(t.pt());
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t)) muPt_.numerator ->Fill(t.pt(),PrescaleWeight);
                }
              }
              if (tkPt(t))
              {
                muEta_.denominator->Fill(t.eta());
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t)) muEta_.numerator ->Fill(t.eta(),PrescaleWeight);
                }
              }
              if (tkPt(t) && tkEta(t))
              {
                muPhi_.denominator->Fill(t.phi());
                if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
                {
                  if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t)) muPhi_.numerator ->Fill(t.phi(),PrescaleWeight);
                }
              }

  	        }
  	      }
    	  break;

    	case 10:
      if(!( muoSelection_ref( m )))continue;
      if(!( muoSelection_ref( m1 )))continue;
  	  if (trHandle.isValid())
      {
  	    for (auto const & t : *trHandle)
        {
  	      if(!trSelection_ref(t)) continue;
  	      const reco::Track& itrk1       = t ;                                                
  	      if((reco::deltaR(t,m1) <= min_dR)) continue;//checking overlaping
  	      if((reco::deltaR(t,m) <= min_dR)) continue;
  	      if (! itrk1.quality(reco::TrackBase::highPurity))     continue;
  	      reco::Particle::LorentzVector pB, p2, p3;
  	      double trackMass2 = kaon_mass * kaon_mass;
  	      double MuMass2 = mu_mass * mu_mass;//0.1056583745 *0.1056583745;
  	      double e2   = sqrt(m1.momentum().Mag2()  + MuMass2          );
  	      double e3   = sqrt(itrk1.momentum().Mag2() + trackMass2  );
  	      p2   = reco::Particle::LorentzVector(m1.px() , m1.py() , m1.pz() , e2  );
  	      p3   = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3  );
  	      pB   = p2 + p3;
  	      if( pB.mass()> maxmassJpsiTk || pB.mass()< minmassJpsiTk) continue;
  	      reco::TransientTrack trTT(itrk1, &(*bFieldHandle));
  	      std::vector<reco::TransientTrack> t_tks;
  	      t_tks.push_back(mu2TT);
  	      t_tks.push_back(trTT);
  	      KalmanVertexFitter kvf;
  	      TransientVertex tv  = kvf.vertex(t_tks);
  	      reco::Vertex vertex = tv;
  	      if (!tv.isValid()) continue;
  	      float JpsiTkCL = 0;
  	      if ((vertex.chi2()>=0.0) && (vertex.ndof()>0) )   
      		JpsiTkCL = TMath::Prob(vertex.chi2(), vertex.ndof() );
  	      math::XYZVector pperp(m1.px() + itrk1.px(),
  				    m1.py() + itrk1.py(),
  				    0.);
  	      GlobalPoint secondaryVertex = tv.position();
  	      GlobalError err             = tv.positionError();
  	      GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0() - secondaryVertex.x()) + 
  							(secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()), 
  						    -1*((vertexBeamSpot.y0() - secondaryVertex.y()) + 
  							(secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()), 
  						    0);
  	      reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
  	      if (JpsiTkCL<minprob) continue;
          if (tkEta(t))
          {
            muPt_.denominator ->Fill(t.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t)) muPt_.numerator ->Fill(t.pt(),PrescaleWeight);
            }
          }
          if (tkPt(t))
          {
            muEta_.denominator->Fill(t.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t)) muEta_.numerator ->Fill(t.eta(),PrescaleWeight);
            }
          }
          if (tkPt(t) && tkEta(t))
          {
            muPhi_.denominator->Fill(t.phi());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t)) muPhi_.numerator ->Fill(t.phi(),PrescaleWeight);
            }
          }
	    }
	  }
  	break;

    case 11:
  	  if (dimuonCL < minprob) continue;
  	  if (fabs(jpsi_cos) < mincos) continue;
  	  if ((displacementFromBeamspotJpsi.perp() / sqrt(jerr.rerr(displacementFromBeamspotJpsi))) < minDS) continue;
  	  if (trHandle.isValid()) {
 	    for (auto const & t : *trHandle) 
      {
 	      if (!trSelection_ref(t)) continue;
 	      if ((reco::deltaR(t,m) <= min_dR)) continue; // checking overlapping
 	      if ((reco::deltaR(t,m1) <= min_dR)) continue;  // checking overlapping
 	      for (auto const & t1 : *trHandle) 
        {
      		if (!(t.pt()>t1.pt())) continue; 
      		if (!trSelection_ref(t1)) continue;
      		if ((reco::deltaR(t1,m) <= min_dR)) continue;  // checking overlapping
      		if ((reco::deltaR(t1,m1) <= min_dR)) continue; // checking overlapping
      		if ((reco::deltaR(t,t1) <= min_dR)) continue;  // checking overlapping
      		const reco::Track& itrk1 = t ;
      		const reco::Track& itrk2 = t1 ;
      		if (! itrk1.quality(reco::TrackBase::highPurity)) continue;
      		if (! itrk2.quality(reco::TrackBase::highPurity)) continue;
      		reco::Particle::LorentzVector pB, pTkTk, p1, p2, p3, p4;
      		double trackMass2 = kaon_mass * kaon_mass;
      		double MuMass2 = mu_mass * mu_mass; // 0.1056583745 *0.1056583745;
      		double e1 = sqrt(m.momentum().Mag2()  + MuMass2 );
      		double e2 = sqrt(m1.momentum().Mag2()  + MuMass2 );
      		double e3 = sqrt(itrk1.momentum().Mag2() + trackMass2  );
      		double e4 = sqrt(itrk2.momentum().Mag2() + trackMass2  );
      		p1 = reco::Particle::LorentzVector(m.px() , m.py() , m.pz() , e1  );
      		p2 = reco::Particle::LorentzVector(m1.px() , m1.py() , m1.pz() , e2  );
      		p3 = reco::Particle::LorentzVector(itrk1.px(), itrk1.py(), itrk1.pz(), e3  );
      		p4 = reco::Particle::LorentzVector(itrk2.px(), itrk2.py(), itrk2.pz(), e4  );
      		pTkTk = p3 + p4;
      		if (pTkTk.mass() > maxmassTkTk || pTkTk.mass() < minmassTkTk) continue;
      		pB = p1 + p2 + p3 + p4;
      		reco::TransientTrack mu1TT(m.track(), &(*bFieldHandle));
      		reco::TransientTrack mu2TT(m1.track(), &(*bFieldHandle));
      		reco::TransientTrack trTT(itrk1, &(*bFieldHandle));
      		reco::TransientTrack tr1TT(itrk2, &(*bFieldHandle));
      		std::vector<reco::TransientTrack> t_tks;
      		t_tks.push_back(mu1TT);
      		t_tks.push_back(mu2TT);
      		t_tks.push_back(trTT);
      		t_tks.push_back(tr1TT);
      		KalmanVertexFitter kvf;
      		TransientVertex tv  = kvf.vertex(t_tks); // this will compare the tracks
      		reco::Vertex vertex = tv;
      		if (!tv.isValid()) continue;
      		float JpsiTkCL = 0;
      		if ((vertex.chi2() >= 0.0) && (vertex.ndof() > 0) )
      		  JpsiTkCL = TMath::Prob(vertex.chi2(), vertex.ndof() );
      		math::XYZVector pperp(m.px() + m1.px() + itrk1.px() + itrk2.px(),
      				      m.py() + m1.py() + itrk1.py() + itrk2.py(),
      				      0.);
      		GlobalPoint secondaryVertex = tv.position();
      		GlobalError err             = tv.positionError();
      		GlobalPoint displacementFromBeamspot( -1*((vertexBeamSpot.x0() - secondaryVertex.x()) +
      							  (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dxdz()),
      						      -1*((vertexBeamSpot.y0() - secondaryVertex.y()) +
      							  (secondaryVertex.z() - vertexBeamSpot.z0()) * vertexBeamSpot.dydz()),
      						      0);
      		reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
      		float jpsiKcos = vperp.Dot(pperp) / (vperp.R()*pperp.R());
      		if (JpsiTkCL < minprob) continue;
      		if (fabs(jpsiKcos) < mincos) continue;
      		if ((displacementFromBeamspot.perp() / sqrt(err.rerr(displacementFromBeamspot))) < minDS) continue;

          if (tkPt(t) && tkEta(t) && tkPt(t1) && tkEta(t1))
          {      
            BMass_.denominator ->Fill(pB.mass());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if( !matchToTrigger(hltpath1,m1))continue;
              if( !matchToTrigger(hltpath1,m))continue;
              if( !matchToTrigger(hltpath1,t))continue;
              if( !matchToTrigger(hltpath1,t1))continue;
              BMass_.numerator ->Fill(pB.mass(),PrescaleWeight);
            }
          }

      		if ( pB.mass() > maxmassJpsiTk || pB.mass()< minmassJpsiTk) continue;
          if (tkEta(t) && tkPt(t1) && tkEta(t1))
          {
            mu1Pt_.denominator ->Fill(t.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t) && matchToTrigger(hltpath1,t1)) mu1Pt_.numerator ->Fill(t.pt(),PrescaleWeight);
            }
          }
          if (tkPt(t) && tkPt(t1) && tkEta(t1))
          {
            mu1Eta_.denominator->Fill(t.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t)&& matchToTrigger(hltpath1,t1)) mu1Eta_.numerator ->Fill(t.eta(),PrescaleWeight);
            }
          }

          if (tkEta(t1) && tkPt(t) && tkEta(t))
          {
            mu2Pt_.denominator ->Fill(t1.pt());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t) && matchToTrigger(hltpath1,t1)) mu2Pt_.numerator ->Fill(t1.pt(),PrescaleWeight);
            }
          }
          if (tkPt(t1) && tkPt(t) && tkEta(t))
          {
            mu2Eta_.denominator->Fill(t1.eta());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t) && matchToTrigger(hltpath1,t1)) mu2Eta_.numerator ->Fill(t1.eta(),PrescaleWeight);
            }
          }


          if (tkPt(t) && tkEta(t) && tkPt(t1) && tkEta(t1)) 
          {
            mu1Phi_.denominator->Fill(t.phi());
            mu2Phi_.denominator->Fill(t1.phi());
            if (num_genTriggerEventFlag_->on() &&  num_genTriggerEventFlag_->accept( iEvent, iSetup) )
            {
              if(matchToTrigger(hltpath1,m1) && matchToTrigger(hltpath1,m) && matchToTrigger(hltpath1,t) && matchToTrigger(hltpath1,t1)) 
              {
                mu1Phi_.numerator ->Fill(t.phi(),PrescaleWeight);
                mu2Phi_.numerator ->Fill(t1.phi(),PrescaleWeight);
              }
            }
          }
 	      } // for (auto const & t1 : *trHandle)
 	    } // for (auto const & t : *trHandle)
   	  } // if (trHandle.isValid())
  	  break;
	    } 
      }
    }
  } 
}

void BPHMonitor::fillHistoPSetDescription(edm::ParameterSetDescription & pset)
{
  pset.addNode((edm::ParameterDescription<int>("nbins", true) and
		edm::ParameterDescription<double>("xmin", true) and
		edm::ParameterDescription<double>("xmax", true)) xor
	       edm::ParameterDescription<std::vector<double>>("edges", true));
}

void BPHMonitor::fillHistoLSPSetDescription(edm::ParameterSetDescription & pset)
{
  pset.add<int> ( "nbins", 2500);
}

void BPHMonitor::fillDescriptions(edm::ConfigurationDescriptions & descriptions)
{
  edm::ParameterSetDescription desc;
  desc.add<std::string>  ( "FolderName", "HLT/BPH/" );
  desc.add<edm::InputTag>( "tracks",  edm::InputTag("generalTracks") );
  desc.add<edm::InputTag>( "photons",  edm::InputTag("photons") );
  desc.add<edm::InputTag>( "offlinePVs", edm::InputTag("offlinePrimaryVertices") );
  desc.add<edm::InputTag>( "beamSpot",edm::InputTag("offlineBeamSpot") );
  desc.add<edm::InputTag>( "muons", edm::InputTag("muons") );
  desc.add<edm::InputTag>( "hltTriggerSummaryAOD", edm::InputTag("hltTriggerSummaryAOD","","HLT") );
  desc.add<std::string>("muoSelection", "");
  desc.add<std::string>("muPt", "");
  desc.add<std::string>("muEta", "");
  desc.add<std::string>("tkPt", "");
  desc.add<std::string>("tkEta", "");
  desc.add<std::string>("dMuEta", "");
  desc.add<std::string>("dMuPt", "");
  desc.add<std::string>("dMuMass", "");
  desc.add<std::string>("muoSelection_ref", "");
  desc.add<std::string>("muQual", "isPFMuon & isGlobalMuon  & innerTrack.hitPattern.trackerLayersWithMeasurement>5 & innerTrack.hitPattern.numberOfValidPixelHits> 0");
  desc.add<std::string>("muoSelection_tag",  ""); // tight selection for tag muon
  desc.add<std::string>("muoSelection_probe", "");
  desc.add<std::string>("trSelection_ref", "");
  desc.add<std::string>("DMSelection_ref", "Pt>4 & abs(Eta)");

  desc.add<int>("nmuons", 1);
  desc.add<bool>( "tnp", false );
  desc.add<int>( "L3", 0 );
  desc.add<double>( "ptCut", 0 );
  desc.add<int>( "displaced", 0 );
  desc.add<int>( "trOrMu", 0 ); // if =0, track param monitoring
  desc.add<int>( "Jpsi", 0 );
  desc.add<int>( "Upsilon", 0 );
  desc.add<int>( "enum", 1 ); // 1...9, 9 sets of variables to be filled, depends on the hlt path
  desc.add<int>( "seagull", 1 ); 
  desc.add<double>( "maxmass", 3.596 );
  desc.add<double>( "minmass", 2.596 );
  desc.add<double>( "maxmassJpsi", 3.2 );
  desc.add<double>( "minmassJpsi", 3. );
  desc.add<double>( "maxmassUpsilon", 10.0 );
  desc.add<double>( "minmassUpsilon", 8.8 );
  desc.add<double>( "maxmassTkTk", 10);
  desc.add<double>( "minmassTkTk", 0);
  desc.add<double>( "maxmassJpsiTk", 5.46 );
  desc.add<double>( "minmassJpsiTk", 5.1 );
  desc.add<double>( "kaon_mass", 0.493677 );
  desc.add<double>( "mu_mass", 0.1056583745);
  desc.add<double>( "min_dR", 0.001);
  desc.add<double>( "max_dR", 1.4);
  desc.add<double>( "minprob", 0.005 );
  desc.add<double>( "mincos", 0.95 );
  desc.add<double>( "minDS", 3. );

  edm::ParameterSetDescription genericTriggerEventPSet;
  genericTriggerEventPSet.add<bool>("andOr");
  genericTriggerEventPSet.add<edm::InputTag>("dcsInputTag", edm::InputTag("scalersRawToDigi") );
  genericTriggerEventPSet.add<edm::InputTag>("hltInputTag", edm::InputTag("TriggerResults::HLT") );
  genericTriggerEventPSet.add<std::vector<int> >("dcsPartitions",{});
  genericTriggerEventPSet.add<bool>("andOrDcs", false);
  genericTriggerEventPSet.add<bool>("errorReplyDcs", true);
  genericTriggerEventPSet.add<std::string>("dbLabel","");
  genericTriggerEventPSet.add<bool>("andOrHlt", true);
  genericTriggerEventPSet.add<bool>("andOrL1", true);
  genericTriggerEventPSet.add<std::vector<std::string> >("hltPaths",{});
  genericTriggerEventPSet.add<std::vector<std::string> >("l1Algorithms",{});
  genericTriggerEventPSet.add<std::string>("hltDBKey","");
  genericTriggerEventPSet.add<bool>("errorReplyHlt",false);
  genericTriggerEventPSet.add<bool>("errorReplyL1",true);
  genericTriggerEventPSet.add<bool>("l1BeforeMask",true);
  genericTriggerEventPSet.add<unsigned int>("verbosityLevel",0);
  desc.add<edm::ParameterSetDescription>("numGenericTriggerEventPSet", genericTriggerEventPSet);
  desc.add<edm::ParameterSetDescription>("denGenericTriggerEventPSet", genericTriggerEventPSet);

//  edm::ParameterSetDescription PrescaleTriggerEventPSet;
//  PrescaleTriggerEventPSet.add<unsigned int>("prescaleWeightVerbosityLevel",0);
//  PrescaleTriggerEventPSet.add<edm::InputTag>("prescaleWeightTriggerResults",edm::InputTag("TriggerResults::HLT"));
//  PrescaleTriggerEventPSet.add<edm::InputTag>("prescaleWeightL1GtTriggerMenuLite",edm::InputTag("l1GtTriggerMenuLite"));
//  PrescaleTriggerEventPSet.add<std::vector<std::string>>("prescaleWeightHltPaths",{});
//  desc.add<edm::ParameterSetDescription>("PrescaleTriggerEventPSet", PrescaleTriggerEventPSet);

  edm::ParameterSetDescription histoPSet;
  edm::ParameterSetDescription phiPSet;
  edm::ParameterSetDescription etaPSet;
  edm::ParameterSetDescription ptPSet;
  edm::ParameterSetDescription dMu_ptPSet;
  edm::ParameterSetDescription d0PSet;
  edm::ParameterSetDescription z0PSet;
  edm::ParameterSetDescription dRPSet;
  edm::ParameterSetDescription massPSet;
  edm::ParameterSetDescription BmassPSet;
  edm::ParameterSetDescription dcaPSet;
  edm::ParameterSetDescription dsPSet;
  edm::ParameterSetDescription cosPSet;
  edm::ParameterSetDescription probPSet;
  edm::ParameterSetDescription TCoPSet;
  edm::ParameterSetDescription PUPSet;
  fillHistoPSetDescription(phiPSet);
  fillHistoPSetDescription(ptPSet);
  fillHistoPSetDescription(dMu_ptPSet);
  fillHistoPSetDescription(etaPSet);
  fillHistoPSetDescription(z0PSet);
  fillHistoPSetDescription(d0PSet);
  fillHistoPSetDescription(dRPSet);
  fillHistoPSetDescription(massPSet);
  fillHistoPSetDescription(BmassPSet);
  fillHistoPSetDescription(dcaPSet);
  fillHistoPSetDescription(dsPSet);
  fillHistoPSetDescription(cosPSet);
  fillHistoPSetDescription(probPSet);
  histoPSet.add<edm::ParameterSetDescription>("d0PSet", d0PSet);
  histoPSet.add<edm::ParameterSetDescription>("etaPSet", etaPSet);
  histoPSet.add<edm::ParameterSetDescription>("phiPSet", phiPSet);
  histoPSet.add<edm::ParameterSetDescription>("ptPSet", ptPSet);
  histoPSet.add<edm::ParameterSetDescription>("dMu_ptPSet", dMu_ptPSet);
  histoPSet.add<edm::ParameterSetDescription>("z0PSet", z0PSet);
  histoPSet.add<edm::ParameterSetDescription>("dRPSet", dRPSet);
  histoPSet.add<edm::ParameterSetDescription>("massPSet", massPSet);
  histoPSet.add<edm::ParameterSetDescription>("BmassPSet", BmassPSet);
  histoPSet.add<edm::ParameterSetDescription>("dcaPSet", dcaPSet);
  histoPSet.add<edm::ParameterSetDescription>("dsPSet", dsPSet);
  histoPSet.add<edm::ParameterSetDescription>("cosPSet", cosPSet);
  histoPSet.add<edm::ParameterSetDescription>("probPSet", probPSet);
  desc.add<edm::ParameterSetDescription>("histoPSet",histoPSet);
  desc.add<unsigned int>("verbosityLevel",0);
  descriptions.add("bphMonitoring", desc);
}

std::string BPHMonitor::getTriggerName(std::string partialName) {

  const std::string trigger_name_tmp = partialName.substr(0,partialName.find("v*"));
  const unsigned int Ntriggers(hltConfig_.size());
  std::string trigger_name = "";
  for (unsigned int i=0;i<Ntriggers;i++) {
    trigger_name = hltConfig_.triggerName(i);
    if ( trigger_name.find(trigger_name_tmp) != std::string::npos ) break;
  }
  // if ( trigger_name=="" ) {
  //   std::cout << "WARNING: Could not find the trigger name." << std::endl;
  // }

  return trigger_name;
}

template <typename T>
bool BPHMonitor::matchToTrigger(const std::string  &theTriggerName , T t){
  //verbosity levels: >42 to print details // >1337 to print the objects of all filters in the trigger path

  bool matched = false;
  //validity check(s)
  if ( !hltConfig_.inited() ) 
    {std::cout << "WARNING: hltConfig uninitialised." << std::endl;return false;}
  
  if ( verbosity_>42 ) 
    {
      std::cout << std::endl << "Performing trigger matching." << std::endl;  
      std::cout << "trying to match: pt/eta/phi: " << t.pt() << "/" << t.eta() << "/" << t.phi() << std::endl;
    }

  //Find the precise trigger name
  std::string trigger_name = getTriggerName(theTriggerName);
  const unsigned int trigger_index = hltConfig_.triggerIndex(trigger_name);
  if (verbosity_>42 ) 
    {std::cout << "trigger name: " << trigger_name << " index: " << trigger_index << " pre-scale: " << hltConfig_.prescaleValue(1,trigger_name) << std::endl;}

  //loop over all the modules for this trigger
  //by default use the last one
  unsigned int Nmodules = hltConfig_.size(trigger_index);
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(trigger_index));
  if (verbosity_>42 )
    {std::cout << "#filters: " << handleTriggerEvent->sizeFilters() << std::endl;}
  unsigned int fIdx=0;
  unsigned int module_label_index=0;
  for (unsigned int i=0;i<Nmodules;i++)
    {
      const unsigned int tmp_fIdx = handleTriggerEvent->filterIndex(edm::InputTag(moduleLabels[i],"",hltInputTag_1.process()));
      if ( tmp_fIdx< handleTriggerEvent->sizeFilters() ) //index of not used filters are set to sizeFilters()
  {
    fIdx = tmp_fIdx;
    module_label_index = i;
//    if ( filterName==moduleLabels[i] ) {break;}
    if ( verbosity_>1337 )
      {
        //print out the objects for all the filters
        std::cout << "module label: " << moduleLabels[module_label_index] << " (type: " << hltConfig_.moduleType(moduleLabels[module_label_index]) << ")" << " index: " << fIdx <<std::endl;
        const trigger::Keys& KEYS(handleTriggerEvent->filterKeys(fIdx));
        const trigger::size_type nK(KEYS.size());
        const trigger::TriggerObjectCollection& TOC(handleTriggerEvent->getObjects());
        for (trigger::size_type i=0; i!=nK; ++i) 
          {
            const trigger::TriggerObject& TO(TOC[KEYS[i]]); 
            std::cout << "   index: " << i << " key: " << KEYS[i] << ": id: "
            << TO.id() << " pt: " << TO.pt() << " eta: " << TO.eta() << " phi: " << TO.phi() << " mass: " << TO.mass()
            << std::endl;
          }
      }//verbosity
  }//good index
    }

  //loop over all the objects in the filter of choice
  if (verbosity_>42 && verbosity_<=1337 )
    {std::cout << "module label: " << moduleLabels[module_label_index] << " (type: " << hltConfig_.moduleType(moduleLabels[module_label_index]) << ")" << " index: " << fIdx <<std::endl;}
  const trigger::Keys& KEYS(handleTriggerEvent->filterKeys(fIdx));
  const trigger::size_type nK(KEYS.size());
  const trigger::TriggerObjectCollection& TOC(handleTriggerEvent->getObjects());
  for (trigger::size_type i=0; i!=nK; ++i) 
    {
      const trigger::TriggerObject& TO(TOC[KEYS[i]]);
      if (verbosity_>42 && verbosity_<=1337 )
  { 
    std::cout << "   index: " << i << " key: " << KEYS[i] << ": id: "
        << TO.id() << " pt: " << TO.pt() << " eta: " << TO.eta() << " phi: " << TO.phi() << " mass: " << TO.mass()
        << std::endl;
  }
      //perform matching: deltaR and pt check
      if( (reco::deltaR(t.eta(), t.phi(),TO.eta(),TO.phi()) <= 0.2) && (TMath::Abs(t.pt()-TO.pt()) < 0.12) )
  {
    if ( verbosity_>42 ) {std::cout << "matched!" << std::endl;}
    matched = true;
  }
      
    }
  return matched;
}

int BPHMonitor::Prescale(const std::string  hltpath1, const std::string  hltpath, edm::Event const& iEvent, edm::EventSetup const& iSetup,  HLTPrescaleProvider* hltPrescale_)
{
  int PrescaleHLT_num = 1;
  int PrescaleHLT_den = 1;
  double Prescale_num = 1;
  int TotalPrescale =1;
  double  L1P=1, HLTP=1;
  bool flag=true;
  std::vector<bool> theSame_den;    
  std::vector<bool> theSame_num;    
//retrieving HLT prescale
  PrescaleHLT_den = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).second;
  PrescaleHLT_num = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).second;
  HLTP =PrescaleHLT_num/std::__gcd(PrescaleHLT_num, PrescaleHLT_den);
  std::cout<<"PrescaleHLT_num"<<PrescaleHLT_num<<"; PrescaleHLT_den = "<<PrescaleHLT_den<<std::endl;
//
//retrieving L1 prescale
//Checking if we have the same l1 seeds in den and num taking into account that they can be written in different order in num and den and some of them can be also switched off

//check if for each den l1 there is the same l1 seed in num 
  if ( (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).first.size() > 0 )
  {  
    for (size_t iSeed=0; iSeed < (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).first.size(); ++iSeed)
    {
      std::string l1_den = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).first.at(iSeed).first;
      int l1_denp = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).first.at(iSeed).second;
      std::cout<<"denominator l1 = "<<l1_denp<<"; at i = "<<iSeed<<l1_den<<std::endl;
      if (l1_denp<1) continue;
      flag = false;
      for (size_t iSeed1=0; iSeed1 < (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.size(); ++iSeed1)
      {
        std::string l1_num = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.at(iSeed1).first;
        int l1_nump= (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.at(iSeed1).second;
        if (!l1_num.compare(l1_den) && (l1_nump>=1))//the same seed  
        {
          flag = true;
          break;
        }
      }
      theSame_den.push_back(flag);
    }
  }
//check if for each num l1 there is the same l1 seed in den 
  if ( (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.size() > 0 )
  {
    for (size_t iSeed=0; iSeed < (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.size(); ++iSeed)
    {
      std::string l1_num = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.at(iSeed).first;
      int l1_nump = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.at(iSeed).second;
      std::cout<<"num l1 = "<<l1_nump<<"; at i = "<<iSeed<<l1_num<<std::endl;
      if (l1_nump<1) continue;
      flag = false;
      for (size_t iSeed1=0; iSeed1 < (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).first.size(); ++iSeed1)
      {
        std::string l1_den = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).first.at(iSeed1).first;
        int l1_denp= (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath)).first.at(iSeed1).second;
        if (!l1_den.compare(l1_num) && (l1_denp>=1))//the same seed
        {
          flag = true;
          break;
        }
      }
      theSame_num.push_back(flag);
    }
  } 
  flag = true;    

  if (theSame_num.size() == theSame_den.size())
  {
    for(size_t i=0; i<theSame_num.size() ; ++i)
    {
      if ((!theSame_num.at(i)) || (!theSame_den.at(i)))
      {
        flag = false;break;     

      }

    }
  }

  if (flag && (theSame_num.size() == theSame_den.size())) 
  {
    L1P = 1; //den and num have the same set of l1 seeds
  }
  else //
  {
///Numerator
    if ( (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.size() > 0 )   
    {
      Prescale_num =1;
      std::cout<<"number of num seeds"<<(hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.size()<<std::endl;
      for (size_t iSeed=0; iSeed < (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.size(); ++iSeed) 
      {
        
        int l1 = (hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.at(iSeed).second;
        if (l1<1) continue; 
        std::cout<<"l1 = "<<l1<<"; at i = "<<iSeed<<(hltPrescale_->prescaleValuesInDetail(iEvent, iSetup, hltpath1)).first.at(iSeed).first<<std::endl;
        if (l1==1){
          Prescale_num =1;
          break;
        }
        else Prescale_num *= 1 - (1.0/(l1));
      }
      if (Prescale_num!=1 )Prescale_num = 1.0 / (1 - Prescale_num);
    }
    L1P = Prescale_num;
  }    
  TotalPrescale = (int)(L1P*HLTP + 0.5);
  std::cout<<"L1P="<<L1P<<"; HLTP= "<<HLTP<<"total = "<<TotalPrescale<<std::endl;
  return TotalPrescale;
}

// Define this as a plug-in


#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(BPHMonitor);
