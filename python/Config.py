import FWCore.ParameterSet.Config as cms
import argparse

from Configuration.Eras.Era_Phase2C17I13M9_cff import Phase2C17I13M9
process = cms.Process('PROD',Phase2C17I13M9)
process.load('Configuration.Geometry.GeometryExtended2026D99Reco_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

# Specify Eta

import sys

eta = sys.argv[2]
energy = sys.argv[3]
#eta = "16"
#energy = "10"
mb = "mb_ngun"
nevents = "500"
cap = "zpos"

fname = '/eos/home-m/mmatthew/Data/KF/CMSSW_13_1_0_pre1/Ntuplizer/' + cap+'/n'+nevents+'/Eta_'+eta+'/singlemuon_flatEGun_hgcalCenter/step3/'

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring('file:'+fname + 'step3_singlemuon_e'+energy+'GeV_eta'+eta+'_'+cap+'_events'+nevents+'_nopu.root'))


process.demo = cms.EDAnalyzer('Base_Analyzer',
   caloParticles = cms.InputTag("mix", "MergedCaloTruth"),
   hgcalRecHitsEE = cms.InputTag("HGCalRecHit", "HGCEERecHits"),
   hgcalRecHitsFH = cms.InputTag("HGCalRecHit", "HGCHEFRecHits"),
   hgcalRecHitsBH = cms.InputTag("HGCalRecHit", "HGCHEBRecHits"),
   KFHits = cms.InputTag("ticlTrackstersKF","KFHits","RECO"),
   eta = cms.string(eta),
   energy = cms.string(energy))
process.p = cms.Path(process.demo)

