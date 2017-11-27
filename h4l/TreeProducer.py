from heppy.framework.analyzer import Analyzer
from heppy.statistics.tree import Tree
from heppy.analyzers.ntuple import *
from heppy.particles.tlv.resonance import Resonance2 as Resonance

from ROOT import TFile

class TreeProducer(Analyzer):

    def beginLoop(self, setup):
        super(TreeProducer, self).beginLoop(setup)
        self.rootfile = TFile('/'.join([self.dirName,
                                        'tree.root']),
                              'recreate')
        self.tree = Tree( 'events', '')
        self.tree.var('weight', float)
        #self.tree.var('met', float)

        bookParticle(self.tree, 'higgs')
        bookParticle(self.tree, 'zed1')
        bookParticle(self.tree, 'zed2')
        bookParticle(self.tree, 'l1')
        bookParticle(self.tree, 'l2')
        bookParticle(self.tree, 'l3')
        bookParticle(self.tree, 'l4')
        bookParticle(self.tree, 'Jets1')
        bookParticle(self.tree, 'bJets1')
        bookMet(self.tree, 'met')

        self.tree.var('nljets', float)
        self.tree.var('nbjets', float)
        self.tree.var('njets', float)



    def process(self, event):
        self.tree.reset()
        zeds = getattr(event, self.cfg_ana.zeds)
        zeds.sort(key=lambda x: abs(x.m()-91.))
        higgses = getattr(event, self.cfg_ana.higgses)
        selected_bs = getattr(event, self.cfg_ana.selected_bs)
        selected_lights = getattr(event, self.cfg_ana.selected_lights)
        leptons = []


        if len(higgses) > 0:  
            self.tree.fill('weight' , event.weight )

            # Reco Higgs
            higgs = higgses[0]
            fillParticle(self.tree, 'higgs', higgs)
            fillParticle(self.tree, 'zed1', zeds[0])
            fillParticle(self.tree, 'zed2', zeds[1])
            fillMet(self.tree, 'met', event.met)

            leptons.append(zeds[0].legs[0])
            leptons.append(zeds[0].legs[1])
            leptons.append(zeds[1].legs[0])
            leptons.append(zeds[1].legs[1])

            leptons.sort(key=lambda x: x.pt(), reverse=True)

            fillLepton(self.tree, 'l1', leptons[0])
            fillLepton(self.tree, 'l2', leptons[1])
            fillLepton(self.tree, 'l3', leptons[2])
            fillLepton(self.tree, 'l4', leptons[3])

	    self.tree.fill('nbjets' , len(event.selected_bs) )
            self.tree.fill('nljets' , len(event.selected_lights) )
            self.tree.fill('njets' , len(event.selected_lights) + len(event.selected_bs))

            if len(event.selected_lights) > 0:
              fillParticle(self.tree, 'Jets1', event.selected_lights[0])
            if len(event.selected_bs) > 0:
              fillParticle(self.tree, 'bJets1', event.selected_bs[0])



            self.tree.tree.Fill()
        
    def write(self, setup):
        self.rootfile.Write()
        self.rootfile.Close()

