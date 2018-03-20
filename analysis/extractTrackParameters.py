'''
Script to read histograms produced by trackParameters.C
And calculate the track parameters from fits to these histograms
'''
from ROOT import *
gROOT.SetBatch(1)
gStyle.SetGridStyle(3) 
gStyle.SetPadLeftMargin(0.15) # increase space for left margin
gStyle.SetGridColor(kGray)


colours = [
    865, # blue
    801, # orange
    629, # red
    418,  # green
    15,  # grey
    618, # purple
    1, # black
]

PI = 3.14159265

# Convenience function to check if an object is a TH3
def is_TH3(h):
        return isinstance(h, TH3)

def appendSlash(path):
    if path[-1] != '/':
        return path+'/'
    else:
        return path

def main(inputFile, outputDir):

    outputDir = appendSlash(outputDir)

    # open file 
    iFile = TFile.Open(inputFile)

    ptRanges = [
            [0,2],
            [9, 11],
            [50, 101]
            ]
    etaRanges = [
            [0, 0.2],
            [0.2, 0.4],
            [0.4, 0.8],
            [0.8, 1.0],
            [1.0, 1.2],
            [1.2, 1.4],
            [1.4, 1.6],
            [1.6, 1.8],
            [1.8, 2.0]
            ]

    parameters = {
            'z0Res'   : {'title' : 'z0 resolution' , 'units' : '[mm]', 'label': '#deltaz_{0}'},
            'd0Res'   : {'title' : 'd0 resolution' , 'units' : '[mm]', 'label': '#deltad_{0}'},
            'ptRes'   : {'title' : 'pt resolution' , 'units' : '[%]',    'label': '#deltap_{T}/p_{T}'},
            'CtgThetaRes': {'title' : 'cot(#theta) resolution', 'units' : '',      'label' : '#deltacot(#theta)'},
            'ptResRaw': {'title' : 'p_{T} resolution',       'units' : '[GeV]', 'label' : '#deltap_{T}'},
            'phiRes'  : {'title' : '#phi resolution',        'units' : '[deg]',      'label' : '#delta#phi'},
            }

    for par in parameters.keys():

        # Get 3D plot
        plotName = par+'_pt_eta'
        resPlot = iFile.Get(plotName)
        if not is_TH3(resPlot):
            print plotName
            print 'ERROR, {0} is not of type TH3, it is {1}'.format(resPlot.GetName(), type(resPlot))
            sys.exti()

        # Create legend 
        #leg = prepareLegend('bottomRight')
        leg = prepareLegend('topLeft')
        leg.SetHeader('Track p_{T} [GeV]')

        # Fit resolution, extract sigma  
        graphs = []
        for n, ptRange in enumerate(ptRanges):
            graphs.append(TGraphErrors())
            for i, etaRange in enumerate(etaRanges):

                xposition = (etaRange[0] + etaRange[1])/2.0
                fitResults = extractResolution(resPlot, ptRange, etaRange, par, outputDir)
                yVal = fitResults['Sigma'][0]
                yErr = fitResults['Sigma'][1]
                if par == "phiRes": # want phi resolution to be in degrees (to compare with tkLayout)
                    yVal *= 180.0/PI
                    yErr *= 180.0/PI 
                graphs[n].SetPoint(i, xposition, yVal)
                graphs[n].SetPointError(i, 0, yErr)
                graphs[n].SetMarkerColor(colours[n])
                graphs[n].SetLineColor(colours[n])
            leg.AddEntry(graphs[n], '{0} < pT < {1}'.format(ptRange[0], ptRange[1]), 'lp')

        # Draw graphs 
        can = TCanvas('can', 'can', 500, 500)
        can.SetGrid()
        can.SetLogy()
        mg = TMultiGraph()
        for g in graphs:
            mg.Add(g, 'p')
        mg.SetTitle('{0} from Delphes; #eta;{1} {2}'.format(parameters[par]['title'], parameters[par]['label'], parameters[par]['units']) )
        mg.Draw('a')
        leg.Draw()

        # Change the axis limits
        mgMin = mg.GetHistogram().GetYaxis().GetXmin()
        mgMax = mg.GetHistogram().GetYaxis().GetXmax()
        if mgMax/10 < mgMin:
            mg.GetHistogram().GetYaxis().SetRangeUser(mgMin*0.5, mgMax*5) # make space for legend
        
        can.SaveAs(outputDir+par+'_graphs.pdf')

def extractResolution(plot3D, ptRange, etaRange, parameterName, outputDir):
    # projection of a slice (resolution, pT, eta) 

    yaxis = plot3D.GetYaxis()
    zaxis = plot3D.GetZaxis()

    # find bins corresponding to pT range
    binLowPt  = yaxis.FindBin(ptRange[0])
    binHighPt = yaxis.FindBin(ptRange[1]) - 1 # want up-to ptRangeMax

    binLowEta  = zaxis.FindBin(etaRange[0])
    binHighEta = zaxis.FindBin(etaRange[1]) - 1 # want up-to etaRangeMax

    # use the SetRange function to pick out a slice in pT (SetRange uses bin numbers) 
    yaxis.SetRange(binLowPt, binHighPt)
    zaxis.SetRange(binLowEta, binHighEta)

    # project  
    res = plot3D.Project3D('x') # project onto the resolution axis

    # fit with gaussian 
    res.Fit("gaus", "MQO") # M: "improve fit", Q: "quiet" (no printing), O: "Don't draw fit or histo
    theFit = res.GetFunction("gaus")

    # Extract fit results
    fitResults = {}
    for i in range(0, 3):
        parName = theFit.GetParName(i)
        value = theFit.GetParameter(i)
        error = theFit.GetParError(i)
        fitResults[parName] = [value, error]

    sigma = fitResults['Sigma'][0]
    sigmaErr = fitResults['Sigma'][1]

    # save plot for diagnostic purposes
    name = '{0}pt{1}_{2}eta{3}'.format(ptRange[0], ptRange[1], int(etaRange[0]*10), int(etaRange[1]*10))
    res.SetTitle('{0} sigma = {1:.2f}#pm{2:.2f}'.format(name, sigma, sigmaErr))
    tempCan = TCanvas('tcan', 'can', 500, 500)
    res.GetXaxis().SetRangeUser(-sigma*4, sigma*4)

    res.Draw()
    tempCan.SaveAs(outputDir+'{0}_{1}.pdf'.format(parameterName, name))

    return fitResults




if __name__ == "__main__":
    from optparse import OptionParser
    import sys

    parser = OptionParser()
    parser.add_option("-i", "--inputFile", action="store", type="string", help="Input ROOT file")
    parser.add_option("-o", "--outputDir", action="store", type="string",  default='./', help="Output directory for plots")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    options, args = parser.parse_args()
    option_dict = dict( (k, v) for k, v in vars(options).iteritems() if v is not None)
    main(**option_dict)

