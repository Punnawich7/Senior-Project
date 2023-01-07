import sys 

import ROOT
from ROOT import TFile, TTree, TChain, gROOT, addressof

gROOT.ProcessLine(
  "struct mssmhbb_str {\
    Int_t    njets;\
    Int_t    nbjets;\
    Float_t  mj1;\
    Float_t  ptj1;\
    Float_t  etaj1;\
    Float_t  phij1;\
    Int_t    btagj1;\
    Int_t    flavorj1;\
    Float_t  mj2;\
    Float_t  ptj2;\
    Float_t  etaj2;\
    Float_t  phij2;\
    Int_t    btagj2;\
    Int_t    flavorj2;\
    Float_t  mj3;\
    Float_t  ptj3;\
    Float_t  etaj3;\
    Float_t  phij3;\
    Int_t    btagj3;\
    Int_t    flavorj3;\
    Float_t  mj4;\
    Float_t  ptj4;\
    Float_t  etaj4;\
    Float_t  phij4;\
    Int_t    btagj4;\
    Int_t    flavorj4;\
    Float_t  pt12;\
    Float_t  eta12;\
    Float_t  phi12;\
    Float_t  m12;\
    Float_t  dR12;\
    Float_t  dEta12;\
    Float_t  dPhi12;\
    Float_t  pt13;\
    Float_t  eta13;\
    Float_t  phi13;\
    Float_t  m13;\
    Float_t  dR13;\
    Float_t  dEta13;\
    Float_t  dPhi13;\
    Float_t  pt14;\
    Float_t  eta14;\
    Float_t  phi14;\
    Float_t  m14;\
    Float_t  dR14;\
    Float_t  dEta14;\
    Float_t  dPhi14;\
    Float_t  pt23;\
    Float_t  eta23;\
    Float_t  phi23;\
    Float_t  m23;\
    Float_t  dR23;\
    Float_t  dEta23;\
    Float_t  dPhi23;\
    Float_t  pt24;\
    Float_t  eta24;\
    Float_t  phi24;\
    Float_t  m24;\
    Float_t  dR24;\
    Float_t  dEta24;\
    Float_t  dPhi24;\
    Float_t  pt34;\
    Float_t  eta34;\
    Float_t  phi34;\
    Float_t  m34;\
    Float_t  dR34;\
    Float_t  dEta34;\
    Float_t  dPhi34;\
    Float_t  pt1234;\
    Float_t  eta1234;\
    Float_t  phi1234;\
    Float_t  m1234;\
    Float_t  ht;\
  };"); 

if len(sys.argv) != 3:
    print('Usage: %s input_file output_file' % sys.argv[0])
    sys.exit(1)

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass

ROOT.gSystem.Load("libDelphes")

inputFile = sys.argv[1]
outputFile = sys.argv[2]
print('Reading from', inputFile, 'and writing to', outputFile)

chain = ROOT.TChain('Delphes')
chain.Add(inputFile)

treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

branchJet = treeReader.UseBranch('Jet')
branchHT  = treeReader.UseBranch('ScalarHT')

outfile = ROOT.TFile.Open(outputFile, 'RECREATE')
outfile.cd()

mssmhbb = ROOT.TTree('mssmhbb', 'mssmhbb')
tree_out = ROOT.mssmhbb_str()

mssmhbb.Branch('njets', addressof(tree_out, 'njets' ), 'njets/I')
mssmhbb.Branch('nbjets', addressof(tree_out, 'nbjets'), 'nbjets/I')
# kinematic Jet1
mssmhbb.Branch('mj1', addressof(tree_out, 'mj1'), 'mj1/F')
mssmhbb.Branch('ptj1', addressof(tree_out, 'ptj1'), 'ptj1/F')
mssmhbb.Branch('etaj1', addressof(tree_out, 'etaj1'), 'etaj1/F')
mssmhbb.Branch('phij1', addressof(tree_out, 'phij1'), 'phij1/F')
mssmhbb.Branch('btagj1', addressof(tree_out, 'btagj1'), 'btagj1/F')
mssmhbb.Branch('flavorj1', addressof(tree_out, 'flavorj1'), 'flavorj1/I')
# kinematic Jet2
mssmhbb.Branch('mj2', addressof(tree_out, 'mj2'), 'mj2/F')
mssmhbb.Branch('ptj2', addressof(tree_out, 'ptj2'), 'ptj2/F')
mssmhbb.Branch('etaj2', addressof(tree_out, 'etaj2'), 'etaj2/F')
mssmhbb.Branch('phij2', addressof(tree_out, 'phij2'), 'phij2/F')
mssmhbb.Branch('btagj2', addressof(tree_out, 'btagj2'), 'btagj2/F')
mssmhbb.Branch('flavorj2', addressof(tree_out, 'flavorj2'), 'flavorj2/I')
# kinematic Jet3
mssmhbb.Branch('mj3', addressof(tree_out, 'mj3'), 'mj3/F')
mssmhbb.Branch('ptj3', addressof(tree_out, 'ptj3'), 'ptj3/F')
mssmhbb.Branch('etaj3', addressof(tree_out, 'etaj3'), 'etaj3/F')
mssmhbb.Branch('phij3', addressof(tree_out, 'phij3'), 'phij3/F')
mssmhbb.Branch('btagj3', addressof(tree_out, 'btagj3'), 'btagj3/F')
mssmhbb.Branch('flavorj3', addressof(tree_out, 'flavorj3'), 'flavorj3/I')
# kinematic Jet4
mssmhbb.Branch('mj4', addressof(tree_out, 'mj4'), 'mj4/F')
mssmhbb.Branch('ptj4', addressof(tree_out, 'ptj4'), 'ptj4/F')
mssmhbb.Branch('etaj4', addressof(tree_out, 'etaj4'), 'etaj4/F')
mssmhbb.Branch('phij4', addressof(tree_out, 'phij4'), 'phij4/F')
mssmhbb.Branch('btagj4', addressof(tree_out, 'btagj4'), 'btagj4/F')
mssmhbb.Branch('flavorj4', addressof(tree_out, 'flavorj4'), 'flavorj4/I')
# 4 vector Jet(1,2)
mssmhbb.Branch('pt12', addressof(tree_out, 'pt12'), 'pt12/F')
mssmhbb.Branch('eta12', addressof(tree_out, 'eta12'), 'eta12/F')
mssmhbb.Branch('phi12', addressof(tree_out, 'phi12'), 'phi12/F')
mssmhbb.Branch('m12', addressof(tree_out, 'm12'), 'm12/F')
mssmhbb.Branch('dR12', addressof(tree_out, 'dR12'), 'dR12/F')
mssmhbb.Branch('dEta12', addressof(tree_out, 'dEta12'), 'dEta12/F')
mssmhbb.Branch('dPhi12', addressof(tree_out, 'dPhi12'), 'dPhi12/F')
# 4 vector Jet(1,3)
mssmhbb.Branch('pt13', addressof(tree_out, 'pt13'), 'pt13/F')
mssmhbb.Branch('eta13', addressof(tree_out, 'eta13'), 'eta13/F')
mssmhbb.Branch('phi13', addressof(tree_out, 'phi13'), 'phi13/F')
mssmhbb.Branch('m13', addressof(tree_out, 'm13'), 'm13/F')
mssmhbb.Branch('dR13', addressof(tree_out, 'dR13'), 'dR13/F')
mssmhbb.Branch('dEta13', addressof(tree_out, 'dEta13'), 'dEta13/F')
mssmhbb.Branch('dPhi13', addressof(tree_out, 'dPhi13'), 'dPhi13/F')
# 4 vector Jet(1,4)
mssmhbb.Branch('pt14', addressof(tree_out, 'pt14'), 'pt14/F')
mssmhbb.Branch('eta14', addressof(tree_out, 'eta14'), 'eta14/F')
mssmhbb.Branch('phi14', addressof(tree_out, 'phi14'), 'phi14/F')
mssmhbb.Branch('m14', addressof(tree_out, 'm14'), 'm14/F')
mssmhbb.Branch('dR14', addressof(tree_out, 'dR14'), 'dR14/F')
mssmhbb.Branch('dEta14', addressof(tree_out, 'dEta14'), 'dEta14/F')
mssmhbb.Branch('dPhi14', addressof(tree_out, 'dPhi14'), 'dPhi14/F')
# 4 vector Jet(2,3)
mssmhbb.Branch('pt23', addressof(tree_out, 'pt23'), 'pt23/F')
mssmhbb.Branch('eta23', addressof(tree_out, 'eta23'), 'eta23/F')
mssmhbb.Branch('phi23', addressof(tree_out, 'phi23'), 'phi23/F')
mssmhbb.Branch('m23', addressof(tree_out, 'm23'), 'm23/F')
mssmhbb.Branch('dR23', addressof(tree_out, 'dR23'), 'dR23/F')
mssmhbb.Branch('dEta23', addressof(tree_out, 'dEta23'), 'dEta23/F')
mssmhbb.Branch('dPhi23', addressof(tree_out, 'dPhi23'), 'dPhi23/F')
# 4 vector Jet(2,4)
mssmhbb.Branch('pt24', addressof(tree_out, 'pt24'), 'pt24/F')
mssmhbb.Branch('eta24', addressof(tree_out, 'eta24'), 'eta24/F')
mssmhbb.Branch('phi24', addressof(tree_out, 'phi24'), 'phi24/F')
mssmhbb.Branch('m24', addressof(tree_out, 'm24'), 'm24/F')
mssmhbb.Branch('dR24', addressof(tree_out, 'dR24'), 'dR24/F')
mssmhbb.Branch('dEta24', addressof(tree_out, 'dEta24'), 'dEta24/F')
mssmhbb.Branch('dPhi24', addressof(tree_out, 'dPhi24'), 'dPhi24/F')
# 4 vector Jet(3,4)
mssmhbb.Branch('pt34', addressof(tree_out, 'pt34'), 'pt34/F')
mssmhbb.Branch('eta34', addressof(tree_out, 'eta34'), 'eta34/F')
mssmhbb.Branch('phi34', addressof(tree_out, 'phi34'), 'phi34/F')
mssmhbb.Branch('m34', addressof(tree_out, 'm34'), 'm34/F')
mssmhbb.Branch('dR34', addressof(tree_out, 'dR34'), 'dR34/F')
mssmhbb.Branch('dEta34', addressof(tree_out, 'dEta34'), 'dEta34/F')
mssmhbb.Branch('dPhi34', addressof(tree_out, 'dPhi34'), 'dPhi34/F')
# 4 vector Jet(1, 2, 3, 4)
mssmhbb.Branch('pt1234', addressof(tree_out, 'pt1234'), 'pt1234/F')
mssmhbb.Branch('eta1234', addressof(tree_out, 'eta1234'), 'eta1234/F')
mssmhbb.Branch('phi1234', addressof(tree_out, 'phi1234'), 'phi1234/F')
mssmhbb.Branch('m1234', addressof(tree_out, 'm1234'), 'm1234/F')
# ScalarHT
mssmhbb.Branch('ht', addressof(tree_out, 'ht'), 'ht/F')

for entry in range(numberOfEntries):
    treeReader.ReadEntry(entry)

    tree_out.njets = int(branchJet.GetEntries())
    tree_out.nbjets = 0

    for i in range(branchJet.GetEntries()):
        jet = branchJet.At(i)
        if jet.BTag == 1: 
            tree_out.nbjets += 1

    if branchJet.GetEntries() > 4: 
        jet1 = branchJet.At(0)
        jet2 = branchJet.At(1)
        jet3 = branchJet.At(2)
        jet4 = branchJet.At(3)

        pjet1 = ROOT.TLorentzVector()
        pjet2 = ROOT.TLorentzVector()
        pjet3 = ROOT.TLorentzVector()
        pjet4 = ROOT.TLorentzVector()

        pjet1.SetPtEtaPhiM(jet1.PT, jet1.Eta, jet1.Phi, jet1.Mass)
        pjet2.SetPtEtaPhiM(jet2.PT, jet2.Eta, jet2.Phi, jet2.Mass)
        pjet3.SetPtEtaPhiM(jet3.PT, jet3.Eta, jet3.Phi, jet3.Mass)
        pjet4.SetPtEtaPhiM(jet4.PT, jet4.Eta, jet4.Phi, jet4.Mass)

        dijet12 = pjet1 + pjet2
        dijet13 = pjet1 + pjet3
        dijet14 = pjet1 + pjet4
        dijet23 = pjet2 + pjet3
        dijet24 = pjet2 + pjet4
        dijet34 = pjet3 + pjet4
        jet1234 = pjet1 + pjet2 + pjet3 + pjet4

        tree_out.mj1      = float(jet1.Mass)
        tree_out.ptj1     = float(jet1.PT)
        tree_out.etaj1    = float(jet1.Eta)
        tree_out.phij1    = float(jet1.Phi)
        tree_out.btagj1   = float(jet1.Btag)
        tree_out.flavorj1 = int(jet1.Flavor)

        tree_out.mj2      = float(jet2.Mass)
        tree_out.ptj2     = float(jet2.PT)
        tree_out.etaj2    = float(jet2.Eta)
        tree_out.phij2    = float(jet2.Phi)
        tree_out.btagj2   = float(jet2.Btag)
        tree_out.flavorj2 = int(jet2.Flavor)

        tree_out.mj3      = float(jet3.Mass)
        tree_out.ptj3     = float(jet3.PT)
        tree_out.etaj3    = float(jet3.Eta)
        tree_out.phij3    = float(jet3.Phi)
        tree_out.btagj3   = float(jet3.Btag)
        tree_out.flavorj3 = int(jet3.Flavor)

        tree_out.mj4      = float(jet4.Mass)
        tree_out.ptj4     = float(jet4.PT)
        tree_out.etaj4    = float(jet4.Eta)
        tree_out.phij4    = float(jet4.Phi)
        tree_out.btagj4   = float(jet4.Btag)
        tree_out.flavorj4 = int(jet4.Flavor)

        tree_out.pt12     = float(dijet12.Pt())
        tree_out.eta12    = float(dijet12.Eta())
        tree_out.phi12    = float(dijet12.Phi())
        tree_out.m12      = float(abs(dijet12.M()))
        tree_out.dR12     = float(pjet1.DeltaR(pjet2))
        tree_out.dEta12   = float(jet1.Eta - jet2.Eta)
        tree_out.dPhi12   = float(pjet1.DeltaPhi(pjet2))

        tree_out.pt13     = float(dijet13.Pt())
        tree_out.eta13    = float(dijet13.Eta())
        tree_out.phi13    = float(dijet13.Phi())
        tree_out.m13      = float(abs(dijet13.M()))
        tree_out.dR13     = float(pjet1.DeltaR(pjet3))
        tree_out.dEta13   = float(jet1.Eta - jet3.Eta)
        tree_out.dPhi13   = float(pjet1.DeltaPhi(pjet3))

        tree_out.pt14     = float(dijet14.Pt())
        tree_out.eta14    = float(dijet14.Eta())
        tree_out.phi14    = float(dijet14.Phi())
        tree_out.m14      = float(abs(dijet14.M()))
        tree_out.dR14     = float(pjet1.DeltaR(pjet4))
        tree_out.dEta14   = float(jet1.Eta - jet4.Eta)
        tree_out.dPhi14   = float(pjet1.DeltaPhi(pjet4))

        tree_out.pt23     = float(dijet23.Pt())
        tree_out.eta23    = float(dijet23.Eta())
        tree_out.phi23    = float(dijet23.Phi())
        tree_out.m23      = float(abs(dijet23.M()))
        tree_out.dR23     = float(pjet2.DeltaR(pjet3))
        tree_out.dEta23   = float(jet2.Eta - jet3.Eta)
        tree_out.dPhi23   = float(pjet2.DeltaPhi(pjet3))

        tree_out.pt24     = float(dijet24.Pt())
        tree_out.eta24    = float(dijet24.Eta())
        tree_out.phi24    = float(dijet24.Phi())
        tree_out.m24      = float(abs(dijet24.M()))
        tree_out.dR24     = float(pjet2.DeltaR(pjet4))
        tree_out.dEta24   = float(jet2.Eta - jet4.Eta)
        tree_out.dPhi24   = float(pjet2.DeltaPhi(pjet4))

        tree_out.pt34     = float(dijet34.Pt())
        tree_out.eta34    = float(dijet34.Eta())
        tree_out.phi34    = float(dijet34.Phi())
        tree_out.m34      = float(abs(dijet34.M()))
        tree_out.dR34     = float(pjet3.DeltaR(pjet4))
        tree_out.dEta34   = float(jet3.Eta - jet4.Eta)
        tree_out.dPhi34   = float(pjet3.DeltaPhi(pjet4))

        tree_out.pt1234  = float(jet1234.Pt())
        tree_out.eta1234 = float(jet1234.Eta())
        tree_out.phi1234 = float(jet1234.Phi())
        tree_out.m1234   = abs(float(jet1234.M()))

    for i in range(branchHT.GetEntries()):
        scalar = branchHT.At(i)
        tree_out.ht = float(scalar.HT())

    mssmhbb.Fill()

mssmhbb.Write()
outfile.Close()

