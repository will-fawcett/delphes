#!/usr/bin/python

#___________________________________________________________________________
def prepareLegend(position, predefined=None):
    import ROOT

    bottomLeft  = [0.15, 0.1, 0.35, 0.3]
    bottomRight = [0.7, 0.1, 0.9, 0.3]
    topRight    = [0.7, 0.7, 0.9, 0.9]
    topLeft     = [0.15, 0.7, 0.35, 0.9]

    if (position == "topLeft"):
        myPosition = topLeft
    if (position == "topRight"):
        myPosition = topRight
    if (position == "bottomLeft"):
        myPosition = bottomLeft
    if (position == "bottomRight"):
        myPosition = bottomRight

    if predefined:
        myPosition = predefined

    return ROOT.TLegend(myPosition[0], myPosition[1], myPosition[2], myPosition[3])

#____________________________________________________________________________
def appendSlash(path):
    if path[-1] != '/':
        return path+'/'
    else:
        return path

#____________________________________________________________________________
def getReverseCumulativeHisto(histo):
    import ROOT
    '''
    Function to return the cumulative efficiency curve from a
    differential efficiency curve
    '''
    nbins = histo.GetNbinsX()
    error = ROOT.Double()
    cumulHisto = histo.Clone()
    cumulHisto.Reset()
    for bin in range(1,nbins+1):
        integral = histo.IntegralAndError(bin,nbins,error)
        integralErr = error
        cumulHisto.SetBinContent(bin,integral)
        cumulHisto.SetBinError(bin,integralErr)
    return cumulHisto
