#!/usr/bin/python

from ROOT import * 
import json
import sys

def addTrack(aDict, cut):

    try:
        aDict[cut] += 1
    except KeyError:
        aDict[cut] = 1


def main(inputFile, outputFile, treeName):


    '''
    resultsDir = "/atlas/data4/userdata/wfawcett/delphes/results/"
    baseDir = resultsDir+"hits_phiEtaSeg_tolerance05mm_phi2GeV_curvature0005_nVertexSigma4/"
    iFile = TFile.Open(baseDir+"hits_ttbar_pu0_multiGeometry_tracks.root", "READ")
    '''
    iFile = TFile.Open(inputFile) 

    
    oFile = TFile(outputFile, "RECREATE")
    # define histograms
    pTcuts = [0, 2, 10, 20, 30, 40, 50]
    deltaKappaTrueHists = {}
    deltaKappaFakeHists = {}
    for pTcut in pTcuts:
        deltaKappaTrueHists[pTcut] = TH1D("deltaKappaTruePt{0}".format(pTcut), "", 1000, -0.03, 0.03)
        deltaKappaFakeHists[pTcut] = TH1D("deltaKappaFakePt{0}".format(pTcut), "", 1000, -0.03, 0.03)


    # get relevant tree
    tree = iFile.Get(treeName)

    numberOfHits = tree.GetEntries()
    print 'There are {0} entries'.format(numberOfHits)

    nFakesOriginal = 0.0

    # 
    kappaCuts = [0.0050, 0.0045, 0.0040, 0.0035, 0.0030, 0.0025]
    beamLineInterceptCuts = [] 

    nTracksSurviveKappaCut = {}
    nFakesSurviveKappaCut = {}

    # high pT fake rate
    nFakesPt50 = 0.0
    nFakesSurviveKappaCutPt50 = {}
    nTracksSurviveKappaCutPt50 = {}

    # Loop over all events in tree
    iEvent = 0 
    for event in tree:
        
        #if iEvent>100000: break
            
        if event.isFake:
            nFakesOriginal += 1 
            if event.hit3pT > 50:
                nFakesPt50 += 1

        deltaKappa = event.kappa_123 - event.kappa_013

        # Fill delta kappa histograms 
        for pTcut in pTcuts:
            if event.pT > pTcut:
                if event.isFake:
                    deltaKappaFakeHists[pTcut].Fill(deltaKappa)
                else:
                    deltaKappaTrueHists[pTcut].Fill(deltaKappa)

        for cut in kappaCuts:

            if abs(deltaKappa) < cut: 
                addTrack(nTracksSurviveKappaCut, cut)

                if event.hit3pT > 50:
                    addTrack(nTracksSurviveKappaCutPt50, cut)

                if event.isFake:
                    addTrack(nFakesSurviveKappaCut, cut)

                    if event.hit3pT > 50:
                        addTrack(nFakesSurviveKappaCutPt50, cut) 

        iEvent += 1

    iFile.Close()
    ################ 
    # End event loop
    ################ 

    # Write histograms
    # open output histogram
    oFile.cd()
    for pTcut in pTcuts:
        deltaKappaTrueHists[pTcut].Write()
        deltaKappaFakeHists[pTcut].Write()

    print 'Original number of fakes', nFakesOriginal
    print 'Original fake rate', nFakesOriginal/float(numberOfHits)

    resultsDict = {}

    for cut in kappaCuts:

        overallTracksRemoved = numberOfHits - nTracksSurviveKappaCut[cut]
        fakesRemoved         = nFakesOriginal - nFakesSurviveKappaCut[cut]

        fakesRemaining       = nFakesOriginal - fakesRemoved
        trueTracksRemoved    = overallTracksRemoved - fakesRemoved

        newFakeRate          = fakesRemaining / nTracksSurviveKappaCut[cut]

        efficiency           = nTracksSurviveKappaCut[cut] / float(numberOfHits)

        resultsDict[cut] = {
                'numberOfHits' : numberOfHits, 
                'nFakesOriginal' : nFakesOriginal,
                'overallTracksRemoved' : overallTracksRemoved,
                'fakesRemaining' : fakesRemaining, 
                'trueTracksRemoved' : trueTracksRemoved, 
                'newFakeRate' : newFakeRate,
                'efficiency' : efficiency
                }

        print 'Kappa cut: ', cut
        print '\tTotal number of tracks removed by cut: {0}'.format(overallTracksRemoved)
        print '\tNumber of fakes removed: {0}'.format(fakesRemoved)
        print '\tNumber of fakes remaining: {0}'.format(fakesRemaining)
        print '\tNumber of true tracks removed: {0}'.format(trueTracksRemoved)
        print '\tEfficiency: {0:.2f}, efficiency loss {1:.2f}'.format(efficiency*100, (1-efficiency)*100)
        print '\tNew fake rate: {0}'.format(newFakeRate)
        print ''

    # Save JSON file
    with open(outputFile.replace('root', 'json'), 'w') as fp:
        json.dump(resultsDict, fp, sort_keys=True, indent=4)


if __name__ == "__main__":


    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--inputFile", action="store", type="string", help="input file name")
    parser.add_option("-o", "--outputFile", action="store", type="string", help="output file name")
    parser.add_option("-t", "--treeName", action="store", type="string", help="Name of the tree to use")
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    options, args = parser.parse_args()
    option_dict = dict( (k, v) for k, v in vars(options).iteritems() if v is not None)

    print 'Input arguments:'
    for opt in option_dict.keys():
        print '\t', opt, option_dict[opt]

    main(**option_dict)
