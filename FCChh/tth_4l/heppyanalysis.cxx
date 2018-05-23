
#include <ROOT/TDataFrame.hxx>
#include "TLorentzVector.h"
#include <TSystem.h>
#include <ROOT/TVec.hxx>

// FCC event datamodel includes
//#include "datamodel/MCParticleCollection.h"
//#include "datamodel/MCParticleData.h"
#include "datamodel/ParticleData.h"
#include "datamodel/TaggedParticleData.h"

/*
 * Running from lxplus:
 *   source /cvmfs/sft.cern.ch/lcg/views/dev3/latest/x86_64-slc6-gcc62-opt/setup.sh
 *   g++ -g -o analysis heppyAnalysis.cxx `root-config --cflags --glibs` -lROOTVecOps -lROOTDataFrame
 */




// Test to reproduce Heppy analysis
int main(int argc, char* argv[]){


   // get input file from command line
   std::string fname;
   if (argc > 1) {
    fname = argv[1]; 
   } else { // no command line arguments given, use some default
     fname = "FCCDelphesOutput.root";
   }
   TFile f1(fname.c_str());

   std::cout << "Read file " << fname << std::endl;

   // fcc edm libraries
   gSystem->Load("libdatamodel.so");
   //f1.MakeProject("dictsDelphes", "*", "RECREATE+");
   //gSystem->Load("dictsDelphes/dictsDelphes.so");

   using ints   = ROOT::Experimental::VecOps::TVec<int>;
   //using fourvectors   = ROOT::Experimental::VecOps::TVec<TLorentzVector>;
   using fourvectors   = std::vector<TLorentzVector>;
   using floats = ROOT::Experimental::VecOps::TVec<float>;
   using u_ints = ROOT::Experimental::VecOps::TVec<unsigned int>;
   //using particles = ROOT::Experimental::VecOps::TVec<fcc::ParticleData>;
   
   std::cout << "Creating TDataFrame ..." << std::endl;
   ROOT::Experimental::TDataFrame df("events", fname);
   //auto d_0_30 = df.Range(0, 10);
   //
   
   auto getInvMass = [] (floats& pxs, floats& pys, floats& pzs, floats& pes) {
     return sqrt(pes*pes - pxs*pxs - pys*pys - pzs * pzs);
     };

   auto getPts = [](floats &pxs, floats &pys){
       auto all_pts = sqrt(pxs * pxs + pys * pys);
       //std::cout << pxs[0] << "\t" << pxs.size() << std::endl;
       return all_pts;
   };

   auto getVecSize = [](ints& pxs) -> Int_t{
       return (Int_t)pxs.size();
   };


   auto selectLeptons = [](std::vector<fcc::ParticleData> in, std::vector<fcc::TaggedParticleData> iso) {
    std::vector<fcc::ParticleData> result;
    result.reserve(in.size());
    for (int i = 0; i < in.size(); ++i) {
      auto & p = in[i];
      if (std::sqrt(std::pow(p.core.p4.px,2) + std::pow(p.core.p4.py,2)) > 20) {
        if (iso[i].tag  < 0.4) {
          result.emplace_back(p);

        }
      }
    }
    return result;
   };

   auto getFourVectorWithCut = [](floats& pxs, floats& pys, floats& pzs, floats& pes, floats& iso) {
      fourvectors fv;
       for (int i = 0; i < pxs.size(); ++i) {
         if ((sqrt(pxs[i] * pxs[i] + pys[i] * pys[i]) > 20) && (iso[i] < 0.4)) {

          fv.emplace_back();
          fv.back().SetXYZM(pxs[i], pys[i], pzs[i], pes[i]);
         }
       }
       return fv;
   };
   auto getChargewithCut = [](ints& charge, floats& pxs, floats& pys, floats& iso) {
     ints chargewithcut;
     for (int i = 0; i < charge.size(); i++) {
       if ((sqrt(pxs[i]*pxs[i] + pys[i]*pys[i]) > 20.) && (iso[i] < 0.4)) {
         chargewithcut.emplace_back(charge[i]);
       }

     }

    return chargewithcut;
   };

   auto getPtFromFourVector = [](fourvectors & fv) ->floats {
     floats pts;
      for (int i = 0; i < fv.size(); ++i) {
        pts.emplace_back(fv[i].Pt());
      }
      return pts;

   };
   auto getMassFromFourVector = [](fourvectors & fv) ->floats {
     floats pts;
      for (int i = 0; i < fv.size(); ++i) {
        pts.emplace_back(fv[i].M());
      }
      return pts;

   };



   auto getPts2 = [](std::vector<fcc::ParticleData> in){
     std::vector<float> result;
       for (int i = 0; i < in.size(); ++i) {
         result.push_back(sqrt(in[i].core.p4.px * in[i].core.p4.px + in[i].core.p4.py * in[i].core.p4.py));
       }
       return result;
   };

   auto id = [](floats &pxs){
       return pxs;
   };

   auto id3 = [](std::vector<float> pxs){
       return pxs;
   };

   auto id2 = [](std::vector<fcc::ParticleData> x) {
     return x;
   };

   auto mergeElectronsAndMuons = [](std::vector<fcc::ParticleData> x, std::vector<fcc::ParticleData> y) {
     std::vector<fcc::ParticleData> result;
     result.reserve(x.size() + y.size());
     result.insert( result.end(), x.begin(), x.end() );
     result.insert( result.end(), y.begin(), y.end() );
     return result;

   };

    //auto findZLegs = [](ints& echarge, ints& mcharge) {
    //std::vector<bool> sel(n);
    //
    //}



   auto buildResonances = [](fourvectors& es, fourvectors & ms, ints& echarge, ints& mcharge) {
     fourvectors higgses;
     //// eeee
     //
     if (es.size() >  3) {
       for (int i = 0; i < es.size() - 3; ++i) {
         for(int j=i+1; j < es.size() - 2; ++j) {
           for(int k=j+1; k < es.size() - 1; ++k) {
             for(int l=k+1; l < es.size(); ++l) {
               if(echarge[i] + echarge[j] + echarge[k] + echarge[l] == 0) {
                 higgses.push_back(TLorentzVector(es[i] + es[j] + es[k] + es[l]));
               }
             }
           }
         }
       }
     }
    //// mumumumu
     if (ms.size() >  3) {
       for (int i = 0; i < ms.size() - 3; ++i) {
         for(int j=i+1; j < ms.size() - 2; ++j) {
           for(int k=j+1; k < ms.size() - 1; ++k) {
             for(int l=k+1; l < ms.size(); ++l) {
               if(mcharge[i] + mcharge[j] + mcharge[k] + mcharge[l] == 0) {
                 higgses.push_back(TLorentzVector(ms[i] + ms[j] + ms[k] + ms[l]));
               }
             }
           }
         }
       }
     }
    //// eemumu
     if (ms.size() >  1 && es.size() > 1) {
       for (int i = 0; i < es.size() - 1; ++i) {
         for(int j=i+1; j < es.size(); ++j) {
           for(int k=0; k < ms.size() - 1; ++k) {
             for(int l=k+1; l < ms.size(); ++l) {
               if((echarge[i] + echarge[j] == 0) && (mcharge[k] + mcharge[l] == 0)) { 
               higgses.push_back(TLorentzVector(es[i] + es[j] + ms[k] + ms[l]));
               }
             }
           }
         }
       }
     }
     auto  higgsresonancesort = [] (TLorentzVector i ,TLorentzVector j) { return (abs( 125. -i.M())<abs(125.-j.M())); };
     std::sort(higgses.begin(), higgses.end(), higgsresonancesort);
     fourvectors finalhiggses;
     if (higgses.size() > 0) {
       finalhiggses.push_back(TLorentzVector(higgses[0]));
     }
     return finalhiggses;

   };

   std::cout << "Apply simple selectors ..." << std::endl;
   auto selectors =  df
                      .Define("selected_electrons", selectLeptons, {"electrons", "electronITags"})
                      .Define("selected_muons", selectLeptons, {"muons", "muonITags"})
                      .Define("selected_leptons", mergeElectronsAndMuons, {"selected_electrons", "selected_muons"})


                    ;
  auto nentries = selectors.Count();
  std::cout << "Count events: " <<  *nentries << std::endl;
  std::cout << "Writing snapshot to disk ..." << std::endl;
  selectors.Snapshot("events", "tree.root",
    { 
      
      "selected_leptons_pt",
      "selected_electrons",
      "selected_muons",
      "selected_leptons",
      }
    );

   return 0;
}
