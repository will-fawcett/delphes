
from ROOT import * 
gROOT.SetBatch(1)
gROOT.SetBatch(1)
gStyle.SetGridStyle(3) 
gStyle.SetPadLeftMargin(0.15) # increase space for left margin
gStyle.SetPadBottomMargin(0.15) # increase space for left margin
gStyle.SetGridColor(kGray)
gStyle.SetPadTickX(1) # add tics on top x
gStyle.SetPadTickY(1) # add tics on right y

from extractTrackParameters import prepareLegend
from Colours import Colours
from functions import appendSlash, getReverseCumulativeHisto  


def main(inputFile, outputDir, text):

    outputDir = appendSlash(outputDir)

    can = TCanvas('can', 'can', 500, 500)
    ifile = TFile.Open(inputFile)

    colours = Colours()
    
    for nJet in range(1,8):

        can.SetLogy(0)
    
        # extract jet pT histogram
        hName = 'nominalJet{0}Pt'.format(nJet) # jets reconstructed using all tracks with pT > 1GeV
        print hName
        h = ifile.Get(hName)
        h.SetLineColor(colours.blue)
        c = getReverseCumulativeHisto(h)
        c.SetMarkerColor(colours.blue)
        c.SetLineColor(colours.blue)

        # get associated jet pT histogram
        h2 = ifile.Get('associatedJet{0}Pt'.format(nJet))
        h2.SetLineColor(colours.red)
        c2 = getReverseCumulativeHisto(h2)
        c2.SetMarkerColor(colours.red)
        c2.SetLineColor(colours.red)

        # Prepare legend
        leg = prepareLegend('topRight')
        leg.SetHeader(text)
        leg.AddEntry(c, 'All tracks', 'lp')
        leg.AddEntry(c2, 'Tracks from PB', 'lp')

        #################
        # draw jet pT histo     
        #################
        h.GetXaxis().SetRangeUser(0, 200)
        h.Draw()
        h2.Draw('same')
        leg.Draw()
        #saveName = outputDir+'jet{0}Pt{1}.pdf'.format(nJet, text)
        saveName = outputDir+hName+'_{0}.pdf'.format(text)
        can.SaveAs(saveName)
        
        #################
        # rate plots
        #################
        can.SetLogy(1)
        c.GetXaxis().SetRangeUser(0, 200)
        c.GetYaxis().SetTitle('Relative Rate')

        #c.DrawNormalized()
        c.Draw()
        c2.Draw('same')
        leg.Draw()
        saveName = outputDir+'jet{0}RateCumulativePt{1}.pdf'.format(nJet, text)
        can.SaveAs(saveName)


    

if __name__ == "__main__":

    from optparse import OptionParser
    import sys

    parser = OptionParser()
    parser.add_option("-i", "--inputFile", action="store", type="string", help="Input ROOT file")
    parser.add_option("-o", "--outputDir", action="store", type="string",  default='./', help="Output directory for plots")
    parser.add_option("-t", "--text", action="store", type="string",  default="", help="text to add to plot names")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    options, args = parser.parse_args()
    option_dict = dict( (k, v) for k, v in vars(options).iteritems() if v is not None)
    main(**option_dict)
