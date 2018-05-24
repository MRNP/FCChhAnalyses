

# build the reproducer of the heppy analysis with TDataFrame
# todo: cleaner directory structure, cmake, etc

all: tdataframe_analysis

tdataframe_analysis:  heppyanalysis.cxx
	g++ -g -o heppyanalysis heppyanalysis.cxx `root-config --cflags --glibs` -lROOTVecOps -lROOTDataFrame  -I/afs/cern.ch/work/v/vavolkl/public/heppy_analysis_samples/fcc-edm/install/include/ -I/afs/cern.ch/work/v/vavolkl/public/heppy_analysis_samples/podio/install/include/ 


# run the analysis

run:
	./heppyanalysis /afs/cern.ch/work/v/vavolkl/public/heppy_analysis_samples/heppy/FCChhAnalyses/FCChh/tth_4l/local/mgp8_pp_tth01j_5f_hllll/events_000001000.root

# todo: comparison with expected output  (/afs/cern.ch/work/v/vavolkl/public/heppy_analysis_samples/localtest/mgp8_pp_tth01j_5f_hllll_1/heppy.FCChhAnalyses.FCChh.tth_4l.TreeProducer.TreeProducer_1/tree.root)
