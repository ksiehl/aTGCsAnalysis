#! /bin/bash

WANTHISTS=--WantToWriteHists
WITHSYSTEMATICS=--withSystematics

echo ./draw --CR ttbar --channel ele --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_ele.root hists_signal-ttbar_ele.root
echo mv hists_ele.root hists-ttbar_ele.root

echo ./draw --CR WJets --channel ele --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_ele.root hists_signal-WJets_ele.root
echo mv hists_ele.root hists-WJets_ele.root

echo ./draw --CR TTbarEnrichedInclusive --channel ele --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_ele.root hists_signal-ttbar-incl_ele.root
echo mv hists_ele.root hists-ttbar-incl_ele.root

echo ./draw --CR TTbarEnrichedBTagVeto --channel ele --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_ele.root hists_signal-ttbar-btagveto_ele.root
echo mv hists_ele.root hists-ttbar-btagveto_ele.root

echo ./draw --CR signal --channel ele --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_ele.root hists_signal-signalBOTH_ele.root
echo mv hists_ele.root hists-signalBOTH_ele.root

echo ./draw --CR signalWW --channel ele --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_ele.root hists_signal-signalWW_ele.root
echo mv hists_ele.root hists-signalWW_ele.root

echo ./draw --CR signalWZ --channel ele --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_ele.root hists_signal-signalWZ_ele.root
echo mv hists_ele.root hists-signalWZ_ele.root

#############################################

echo ./draw --CR ttbar --channel mu --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_mu.root hists_signal-ttbar_mu.root
echo mv hists_mu.root hists-ttbar_mu.root

echo ./draw --CR WJets --channel mu --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_mu.root hists_signal-WJets_mu.root
echo mv hists_mu.root hists-WJets_mu.root

echo ./draw --CR TTbarEnrichedInclusive --channel mu --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_mu.root hists_signal-ttbar-incl_mu.root
echo mv hists_mu.root hists-ttbar-incl_mu.root

echo ./draw --CR TTbarEnrichedBTagVeto --channel mu --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_mu.root hists_signal-ttbar-btagveto_mu.root
echo mv hists_mu.root hists-ttbar-btagveto_mu.root

echo ./draw --CR signal --channel mu --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_mu.root hists_signal-signalBOTH_mu.root
echo mv hists_mu.root hists-signalBOTH_mu.root

echo ./draw --CR signalWW --channel mu --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_mu.root hists_signal-signalWW_mu.root
echo mv hists_mu.root hists-signalWW_mu.root

echo ./draw --CR signalWZ --channel mu --input /afs/cern.ch/work/k/ksiehl/public/ansar_project/ntuple_output_storage/ --withSignal --withMC --withData $WANTHISTS $WITHSYSTEMATICS
echo mv hists_signal_mu.root hists_signal-signalWZ_mu.root
echo mv hists_mu.root hists-signalWZ_mu.root
