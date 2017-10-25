# ######################################################################
# t3ui01
# /mnt/t3nfs01/data01/shome/wiederkehr_s/DQM_930pre3/src/DQMOffline/Trigger/test/jobs/HLT_300816
# file list contains 23155 events
# mkPyFiles -t ../hltBPHdqm-Run300816-XXXX.py -f ../catalogues/Run300816 -n 1
# ./hltBPHdqm-Run300816-Run300816-0000.py with 23155 events
# ######################################################################
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: hltBPHdqm -s DQM:bphMonitorHLT --conditions=92X_dataRun2_Prompt_v8 --geometry DB:Extended --eventcontent DQM --datatier DQMIO --data --era Run2_2017 --filetype=EDM --filein file:/afs/cern.ch/work/a/aboletti/public/BPH_DQM/CMSSW_930p1_test/src/DQMOffline/Trigger/test/Run2017C-Charmonium-AOD-PRv1.root -n -1 --fileout DQM_onlyBPHHLT.root --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

# process = cms.Process('DQM',eras.Run2_2017)
process = cms.Process('DQM')


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
# process.load("Configuration.StandardSequences.Reconstruction_cff") #bmm
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
# process.load("Configuration.Geometry.GeometryDB_cff") #bmm
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
# process.load("Configuration.StandardSequences.MagneticField_38T_cff") #bmm
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# process.MessageLogger.categories.append('HLTrigReport')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
 "PoolSource",
  fileNames = cms.untracked.vstring(
         "/store/data/Run2017C/Charmonium/AOD/PromptReco-v3/000/301/046/00000/18D47622-0E83-E711-AF10-02163E014535.root"
  )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('hltBPHdqm nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition 
rootFileName = "hltBPHdqm-Run301046-Run301046-0000.root"
process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string(rootFileName),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v9', '')

# Path and EndPath definitions
process.dqmoffline_step = cms.EndPath(process.bphMonitorHLT)
process.dqmofflineOnPAT_step = cms.EndPath(process.bphMonitorHLT)
process.DQMoutput_step = cms.EndPath(process.DQMoutput)

#HLT Summary
# process.load( "HLTrigger.HLTanalyzers.hlTrigReport_cfi" )
# process.hlTrigReport.HLTriggerResults   = cms.InputTag("TriggerResults", "", "HLT")
# process.hlTrigReport.ReferencePath      = cms.untracked.string( "HLTriggerFinalPath" )
# process.hlTrigReport.ReferenceRate      = cms.untracked.double( 100.0 )

# Schedule definition
process.schedule = cms.Schedule(process.dqmoffline_step,process.dqmofflineOnPAT_step,process.DQMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
