import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                #'file:/afs/cern.ch/cms/Tutorials/workbook_twiki2021/MinBias_pythia8_14TeV_100events.root'
                                #'file:/afs/cern.ch/work/s/sanjeev/private/2022/a2/CMSSW_10_6_28/src/BPH-RunIISummer20UL18MiniAODv2-00122_110.root'
                                #'root://se01.indiacms.res.in//cms/store/data/Run2018D/ParkingBPH3/MINIAOD/UL2018_MiniAODv2-v1/2520000/12B1CCAB-39FF-5D43-B7A6-77214E8877EC.root'
                                #'/store/data/Run2018D/ParkingBPH3/MINIAOD/UL2018_MiniAODv2-v1/2520000/12B1CCAB-39FF-5D43-B7A6-77214E8877EC.root',
				'/store/data/Run2018D/ParkingBPH2/MINIAOD/UL2018_MiniAODv2-v1/2560006/ED5D2CE5-1D5D-6445-AACB-EC0F6B6FA89A.root',
                        )
)

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v35' ) #auto:run2_data')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.demo = cms.EDAnalyzer('BsToEMu',
                              #OutputFileName = cms.string("output.root"),
                              #tracks    = cms.untracked.InputTag('generalTracks'),
                              muons     = cms.InputTag("slimmedMuons"),
                              beamSpot = cms.InputTag("offlineBeamSpot"),
                              electrons = cms.InputTag("slimmedElectrons"),
                              pfCands   = cms.InputTag("packedPFCandidates"),
                              primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              ltrk = cms.InputTag('lostTracks'),
                              trg_res = cms.InputTag('TriggerResults', '', 'HLT'),
                              #prg_ps = cms.InputTag('patTrigger'),
                              tobjs = cms.InputTag('slimmedPatTrigger'),

                              #genpr = cms.InputTag('prunedGenParticles'),
                              #genpk = cms.InputTag('packedGenParticles'),
                              #genInfo = cms.InputTag('generator'),
                              isMC = cms.bool(False),
                              #isMC = cms.bool(True),
                              OnlyGen = cms.bool(False),
                              #TriggerNames = cms.vstring('HLT_Mu7_IP4*',),
                              #hltResults = cms.InputTag( "TriggerResults", "", "HLT" ),

                          )

#process.triggerSelection = cms.EDFilter("TriggerResultsFilter", triggerConditions = cms.vstring(
#v2.0 'HLT_Mu9_IP6_part', 'HLT_Mu8p5_IP3p5', 'HLT_Mu10p5_IP3p5', 'HLT_Mu8_IP3',
#v2.2  'HLT_Mu12_IP6', 'HLT_Mu9_IP5', 'HLT_Mu7_IP4*', 'HLT_Mu9_IP4', 'HLT_Mu8_IP5', 'HLT_Mu8_IP6', 
#v3.5  'HLT_Mu9_IP3', 'HLT_Mu9_IP0'
#),  hltResults = cms.InputTag( "TriggerResults", "", "HLT" ), l1tResults = cms.InputTag( "" ), throw = cms.bool(False) )

'''
v2.0
  HLT_Mu9_IP6_part : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9Q
  HLT_Mu8p5_IP3p5 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8p5Q
  HLT_Mu10p5_IP3p5 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered10p5Q
  HLT_Mu8_IP3 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8Q

v2.2
  HLT_Mu12_IP6 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered12Q
  HLT_Mu9_IP5 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP5Q
  HLT_Mu7_IP4 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered7IP4Q
  HLT_Mu9_IP4 :hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP4Q
  HLT_Mu8_IP5 :hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP5Q
  HLT_Mu8_IP6 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered8IP6Q

v3.5
  HLT_Mu9_IP3 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP3Q
  HLT_Mu9_IP0 : hltL3fL1sMu22OrParkL1f0L2f10QL3Filtered9IP0Q
'''

process.load('Configuration.StandardSequences.Services_cff')
#process.TFileService = cms.Service("TFileService", fileName = cms.string('Bs2EMu_Feb23_mc.root'))
process.TFileService = cms.Service("TFileService", fileName = cms.string('Bs2EMu_Feb23_data.root'))

process.p = cms.Path(process.demo)
