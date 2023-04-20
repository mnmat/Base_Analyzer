// -*- C++ -*-
//
// Package:    PlayGround/Base_Analyzer
// Class:      Base_Analyzer
//
/**\class Base_Analyzer Base_Analyzer.cc PlayGround/Base_Analyzer/plugins/Base_Analyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Mark Nicholas Matthewman
//         Created:  Fri, 14 Apr 2023 18:15:48 GMT
//
//

// system include files
#include <memory>
#include <numeric>
#include <sstream>
#include <any>
#include <iomanip>
#include <cmath>
#include <typeinfo>
#include <iostream>
#include <fstream>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TProfile.h"

#include "SimDataFormats/CaloAnalysis/interface/CaloParticle.h"
#include "SimDataFormats/CaloAnalysis/interface/CaloParticleFwd.h"

#include "DataFormats/HGCalReco/interface/Trackster.h"
#include "DataFormats/HGCalReco/interface/KFHit.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/ForwardDetId/interface/HGCalDetId.h"

#include "SimDataFormats/CaloAnalysis/interface/SimClusterFwd.h"
#include "SimDataFormats/CaloAnalysis/interface/SimCluster.h"

#include "DataFormats/HGCRecHit/interface/HGCRecHitCollections.h"

#include "DataFormats/ParticleFlowReco/interface/HGCalMultiCluster.h"

#include "RecoLocalCalo/HGCalRecAlgos/interface/RecHitTools.h"

#include "Geometry/CaloTopology/interface/HGCalTopology.h"
#include "Geometry/HGCalGeometry/interface/HGCalGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/Math/interface/Vector3D.h"


//ROOT includes
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <TFile.h>
#include <TROOT.h>
#include <TVector.h>
#include "TBranch.h"
#include <string>
#include <vector>
#include "TSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "THStack.h"
#include <algorithm>
#include "TString.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TPaveStats.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class Base_Analyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Base_Analyzer(const edm::ParameterSet&);
  ~Base_Analyzer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;
  virtual void fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap, 
    const HGCRecHitCollection& rechitsEE, 
    const HGCRecHitCollection& rechitsFH,
    const HGCRecHitCollection& rechitsBH) const;

  hgcal::RecHitTools recHitTools_;


  // ----------member data ---------------------------
  edm::EDGetTokenT<std::vector<CaloParticle>> caloParticlesToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsEEToken_; 
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsFHToken_;
  edm::EDGetTokenT<HGCRecHitCollection> hgcalRecHitsBHToken_;

  edm::EDGetTokenT<std::vector<KFHit>> KFHitsToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeomToken_;

  // TTree *tree = new TTree("tree","tree");

  std::vector<std::string> detectors, objects, positions, hittypes;

  // variables

  int eventnr =0;
  std::string eta_;
  std::string energy_;
  std::string outdir_;
  std::shared_ptr<hgcal::RecHitTools> recHitTools;


  // KF

  /*
  std::vector<float> vec_x;
  std::vector<float> vec_y;
  std::vector<float> vec_z;
  std::vector<float> vec_e;
  std::vector<float> vec_cov_xx;
  std::vector<float> vec_cov_yy;
  std::vector<float> vec_cov_yy;
  std::vector<int> vec_detid;
  std::vector<int> kf_charge;
  std::vector<int> kf_layer;
  std::vector<std::string> vec_dtype;
  std::vector<int> vec_evt;
  */

  edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

void Base_Analyzer::fillHitMap(std::map<DetId, const HGCRecHit*>& hitMap,
                                const HGCRecHitCollection& rechitsEE,
                                const HGCRecHitCollection& rechitsFH,
                                const HGCRecHitCollection& rechitsBH) const {
  hitMap.clear();
  for (const auto& hit : rechitsEE) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsFH) {
    hitMap.emplace(hit.detid(), &hit);
  }

  for (const auto& hit : rechitsBH) {
    hitMap.emplace(hit.detid(), &hit);
  }
} // end of EfficiencyStudies::fillHitMap

//
// constructors and destructor
//
Base_Analyzer::Base_Analyzer(const edm::ParameterSet& iConfig)
    : 
      caloParticlesToken_(consumes<std::vector<CaloParticle> >(iConfig.getParameter<edm::InputTag>("caloParticles"))), 
      hgcalRecHitsEEToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsEE"))),
      hgcalRecHitsFHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsFH"))),
      hgcalRecHitsBHToken_(consumes<HGCRecHitCollection>(iConfig.getParameter<edm::InputTag>("hgcalRecHitsBH"))),
      KFHitsToken_(consumes<std::vector<KFHit>>(iConfig.getParameter<edm::InputTag>("KFHits"))),
      caloGeomToken_(esConsumes<CaloGeometry, CaloGeometryRecord>()),
      eta_(iConfig.getParameter<std::string>("eta")),
      energy_(iConfig.getParameter<std::string>("energy")){

  detectors = {"", "Si", "Si 120", "Si 200", "Si 300", "Sc"};
  hittypes = {"Simhits","Rechits","KF"};
  objects = {"Simhits", "Rechits"};

  /*

  usesResource("TFileService");
  edm::Service<TFileService> file;

  tree = file->make<TTree>("tree","tree");

  */

  /*

  // KF

  tree->Branch("x", &vec_x);
  tree->Branch("y", &vec_y);
  tree->Branch("z", &vec_z);
  tree->Branch("e", &vec_e);
  tree->Branch("layer", &kf_layer);
  tree->Branch("detid", &vec_detid);
  tree->Branch("dtype", &vec_dtype);
  tree->Branch("cov_xx", &vec_cov_xx);
  tree->Branch("cov_xy", &vec_cov_yy);
  tree->Branch("cov_yy", &vec_cov_yy);
  tree->Branch("evt", &vec_evt);

  */

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

Base_Analyzer::~Base_Analyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void Base_Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  edm::Handle<std::vector<CaloParticle>> CaloParticles;
  iEvent.getByToken(caloParticlesToken_, CaloParticles);
  const CaloParticleCollection& cps = *CaloParticles;

  edm::Handle<HGCRecHitCollection> recHitHandleEE;  
  iEvent.getByToken(hgcalRecHitsEEToken_, recHitHandleEE);

  edm::Handle<HGCRecHitCollection> recHitHandleFH;
  iEvent.getByToken(hgcalRecHitsFHToken_, recHitHandleFH);

  edm::Handle<HGCRecHitCollection> recHitHandleBH;
  iEvent.getByToken(hgcalRecHitsBHToken_, recHitHandleBH);

  std::map<DetId, const HGCRecHit*> hitMap;
  fillHitMap(hitMap, *recHitHandleEE, *recHitHandleFH, *recHitHandleBH);

  edm::Handle<std::vector<KFHit>> KFHitsHandle;
  iEvent.getByToken(KFHitsToken_, KFHitsHandle);
  const std::vector<KFHit> hits = *KFHitsHandle;

  // init vars
  std::map<float, GlobalPoint> map_gps;
  // std::map<float, float> map_xx_kf, map_xy_kf, map_yy_kf;

  const CaloGeometry &geom = iSetup.getData(caloGeomToken_);
  recHitTools_.setGeometry(geom);

    // KF

  float kf_energy = 0;
  int kf_hits = 0;

  for(int i = 0;i<int(hits.size());i++){
    std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(hits[i].detid);
    float e = 0;
    if (itcheck != hitMap.end()){
      e = hitMap[hits[i].detid]->energy();
    }
    unsigned int layer_ = recHitTools_.getLayerWithOffset(hits[i].detid);


    /*
    std::string detector;
    std::string thickness;
    std::string tmp;

    if(recHitTools_.isSilicon(hits[i].detid)){
      detector = "Si";
      thickness = std::to_string(int(recHitTools_.getSiThickness(hits[i].detid))); 
      tmp = detector+" "+thickness;
    }
    else{
      detector = "Sc";
      thickness = "None";
      tmp = "Sc";
    } 
    */

    kf_energy+=e;
    map_gps[hits[i].center.z()]=hits[i].center;

    /*
    vec_x.push_back(hits[i].center.x());
    vec_y.push_back(hits[i].center.y());
    vec_z.push_back(hits[i].center.z());
    vec_e.push_back(e);
    vec_detid.push_back(hits[i].detid);
    vec_layer.push_back(layer_);
    vec_dtype.push_back(tmp);
    vec_evt.push_back(eventnr);

    */
  }

  // Loop over Caloparticles 

  // std::vector<DetId> tmprechits_; tmprechits_.clear();
  int comp_hits = 0;
  int rhit_e = 0;
  int rhits = 0;

  for (const auto& it_cp : cps) {
    // do something with track parameters, e.g, plot the charge.
    const CaloParticle& cp = ((it_cp)); 
    const SimClusterRefVector& simclusters = cp.simClusters();

    for (const auto& it_simc : simclusters){
      const SimCluster& simc = (*(it_simc));
      const auto& sc_hae = simc.hits_and_energies();

      for (const auto& it_sc_hae : sc_hae){

        DetId detid_ = (it_sc_hae.first);
        std::map<DetId,const HGCRecHit *>::const_iterator itcheck = hitMap.find(detid_);
        unsigned int layer_ = recHitTools_.getLayerWithOffset(detid_);

        std::string detector;
        std::string thickness;
        std::string tmp;

        if(recHitTools_.isSilicon(detid_)){
          detector = "Si";
          thickness = std::to_string(int(recHitTools_.getSiThickness(detid_))); 
          tmp = detector+" "+thickness;
        }
        else{
          detector = "Sc";
          thickness = "None";
          tmp = "Sc";
        } 


        if (itcheck != hitMap.end()){
          rhits+=1;
          rhit_e+=hitMap[detid_]->energy();

          auto gp = map_gps[recHitTools_.getPosition(detid_).z()];
          DetId closest_detid;
          if (detector == "Sc") closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCell(gp);
          else closest_detid = static_cast<const HGCalGeometry*>(recHitTools_.getSubdetectorGeometry(detid_))->getClosestCellHex(gp, true);
          if(detid_==closest_detid){
            comp_hits+=1;
          }
        }
      }
    }
  }



  std::cout << energy_ << "," << comp_hits << "," << rhits << "," << kf_energy << "," <<rhit_e << std::endl;
  eventnr=eventnr+1;

  /*

  tree->Fill();

  // KF

  vec_x.clear();
  vec_y.clear();
  vec_z.clear();
  vec_e.clear();
  vec_detid.clear();
  kf_layer.clear();
  vec_dtype.clear();
  vec_cov_xx.clear();
  vec_cov_yy.clear();
  vec_cov_yy.clear();
  vec_evt.clear();


  //clear_arrays();

  */

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void Base_Analyzer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void Base_Analyzer::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Base_Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Base_Analyzer);
