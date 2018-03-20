#!/usr/bin/python

'''
Script to ... 
'''

from ROOT import * 
gROOT.SetBatch(1)
gROOT.SetBatch(1)
gStyle.SetGridStyle(3) 
gStyle.SetPadLeftMargin(0.15) # increase space for left margin
gStyle.SetPadBottomMargin(0.15) # increase space for left margin
gStyle.SetGridColor(kGray)
gStyle.SetPadTickX(1) # add tics on top x
gStyle.SetPadTickY(1) # add tics on right y

from Colours import Colours
from functions import appendSlash, getReverseCumulativeHisto, prepareLegend

resultsPath = '/atlas/data4/userdata/wfawcett/delphes/results/'
resultsPath = '/Users/Will/Documents/fcc/delphes/results/'
nEvents = 1000.0

OUTPUT_DIR = 'plots/'

def main(verbose):

    # Define canvas 
    can = TCanvas('can', 'can', 500, 500)
    can.SetLogy()
    can.SetGrid()

    # Open sample files 
    ifile_minbias = TFile.Open(resultsPath+'hist_MinBias_mu200_s80_n1000.root')
    ifile_ttbar = TFile.Open(resultsPath+'hist_ttbar_mu200_s80_n1000.root') 

    colours = Colours()
    
    for nJet in range(1,8):

        nominalName = 'nominalJet{0}Pt'.format(nJet) # jets reconstructed using all tracks with pT > 1GeV
        suppressedName = 'associatedJet{0}Pt'.format(nJet)

        # get jet pT histos from minbias
        mb_nominal = ifile_minbias.Get(nominalName)
        mb_PUsupp  = ifile_minbias.Get(suppressedName)

        # get jet pT histos for ttbar
        tt_nominal = ifile_ttbar.Get(nominalName)
        tt_PUsupp  = ifile_ttbar.Get(suppressedName)

        tt_nominal.GetXaxis().SetRangeUser(0,100)
        tt_nominal.Draw()
        tt_PUsupp.SetLineColor(kRed)
        tt_PUsupp.Draw('same')
        #can.SaveAs('test2.pdf')

        # get cumulative histograms for minbias 
        cmb_nominal = getReverseCumulativeHisto(mb_nominal)
        cmb_PUsupp  = getReverseCumulativeHisto(mb_PUsupp)

        # get cumulative histograms for ttbar
        ctt_nominal = getReverseCumulativeHisto(tt_nominal)
        ctt_PUsupp  = getReverseCumulativeHisto(tt_PUsupp)

        # Style
        xaxis = cmb_nominal.GetXaxis()
        yaxis = cmb_nominal.GetYaxis()
        yaxis.SetTitle('Relative Rate')
        setStyle(cmb_nominal, colours.red)
        setStyle(cmb_PUsupp, colours.blue)

        setStyle(ctt_nominal, colours.red)
        setStyle(ctt_PUsupp, colours.blue)


        #######################
        # Draw for minimum bias 
        #######################
        xaxis.SetRangeUser(0, 80)
        cmb_nominal.Draw()
        cmb_PUsupp.Draw('same')
        
        # Add a legend
        predefined = [0.6, 0.7, 0.9, 0.9]
        leg = prepareLegend('topRight', predefined)
        leg.AddEntry(cmb_nominal, 'All tracks (minbias)', 'lp')
        leg.AddEntry(cmb_PUsupp, 'Tracks from PB (minbias)', 'lp')
        leg.Draw()
        can.SaveAs(OUTPUT_DIR+'triggerRate{0}jets_minbias.pdf'.format(nJet))

        #######################
        # Draw for ttbar        
        #######################
        ctt_nominal.Draw()
        ctt_PUsupp.Draw('same')
        leg.Clear()
        leg.AddEntry(ctt_nominal, 'All tracks (ttbar)', 'lp')
        leg.AddEntry(ctt_PUsupp, 'Tracks from PB (ttbar)', 'lp')
        leg.Draw()
        can.SaveAs(OUTPUT_DIR+'triggerRate{0}jets_ttbar.pdf'.format(nJet))


def setStyle(hist, c):
    hist.SetLineColor(c)
    hist.SetMarkerColor(c)
    hist.Scale(1.0/nEvents)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", help="Turn on verbose messages", action="store_true", default=False)
    args = parser.parse_args()
    verbose = args.verbose

    main(verbose)


