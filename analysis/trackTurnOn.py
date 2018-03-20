
from extractTrackParameters import prepareLegend
from Colours import Colours


from ROOT import * 
gROOT.SetBatch(1)
gStyle.SetGridStyle(3) 
gStyle.SetPadLeftMargin(0.15) # increase space for left margin
gStyle.SetPadBottomMargin(0.15) # increase space for left margin
gStyle.SetGridColor(kGray)
gStyle.SetPadTickX(1) # add tics on top x
gStyle.SetPadTickY(1) # add tics on right y

OUTPUT_DIR = 'plots/'
REBIN = 2

def main():

    ifile = TFile.Open('/Users/Will/Documents/fcc/trackerSW/delphes/output_ttbar_mu1000.root')
    colourDef = Colours()

    truthTrackPt = ifile.Get('truthTrack100')
    truthTrackPt.Rebin(REBIN)
    #truthTrackPt = TH1D('tracks', '', 100, 0, 100)

    '''
    for bin in range(truthTrackPt_1000.GetNbinsX()):
        if bin > 100: continue
        truthTrackPt.SetBinContent(bin, truthTrackPt_1000.GetBinContent(bin))
    truthTrackPt_1000.GetXaxis().SetRangeUser(0,200)
    truthTrackPt_1000.Draw()
    truthTrackPt.SetLineColor(kGreen)
    truthTrackPt.Draw('same')
    can.SaveAs('test.pdf')
    '''
    can = TCanvas('can', 'can', 500, 500)

    line = TF1('line', '1', 0, 100)
    line.SetLineColor(kGray)

    tGraphs = {}

    leg = prepareLegend('bottomRight', [0.7, 0.15, 0.9, 0.35])

    for i in range(0, 6):
        ptCut = (i+1)*5
        hName = 'truthTrackPt{0}'.format(ptCut)
        print hName
        ptAfterCut = ifile.Get(hName)
        ptAfterCut.SetLineColor(kRed)
        ptAfterCut.Rebin(REBIN)

        can.SetLogy()
        truthTrackPt.Draw()
        ptAfterCut.Draw('same')
        can.SaveAs(OUTPUT_DIR+'tracksPt{0}.pdf'.format(ptCut))

        # to make turn on to TGraphAsymmErrors(numerator, denominator)
        ratio = TGraphAsymmErrors(ptAfterCut, truthTrackPt)

        can.SetLogy(0)
        ratio.Draw('AP')
        line.Draw('same')
        xaxis = ratio.GetXaxis()
        xaxis.SetRangeUser(0, ptCut*3)
        xaxis.SetTitle('Truth track p_{T} [GeV]')
        yaxis = ratio.GetYaxis()
        yaxis.SetTitle('Efficiency')

        can.SaveAs(OUTPUT_DIR+'turnOnPt{0}.pdf'.format(ptCut))

        tGraphs[ptCut] = ratio

    # now draw series of TGraphs
    ptCuts = [5, 10, 15, 20]
    colours = [colourDef.blue, colourDef.red, colourDef.orange, colourDef.purple]
    for i, cut in enumerate(ptCuts):
        gr = tGraphs[cut]
        gr.SetLineColor(colours[i])
        gr.SetMarkerColor(colours[i])
        leg.AddEntry(gr, 'p_{T} > '+str(cut)+' GeV')
        if i==0:
            gr.Draw('APl')
            gr.SetMinimum(0)
            gr.GetXaxis().SetRangeUser(0, 45)
            line.Draw('same')
            gr.Draw('Psame')
        else:
            gr.Draw('Plsame')
    leg.Draw()
    
    can.SaveAs(OUTPUT_DIR+'trackTurnOn.pdf')

    

if __name__ == "__main__":
    main()
