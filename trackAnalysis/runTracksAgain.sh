PILEUP=200
resultsPath=/atlas/data4/userdata/wfawcett/delphes/results/hits_phiEtaSeg_tolerance05mm_phi2GeV_curvature0005_nVertexSigma5
outputPath=/atlas/data4/userdata/wfawcett/delphes/results/kappaOptimisation

for spacing in {10..50..10}; do
  echo python tracksAgain.py -i $resultsPath/hits_ttbar_pu${PILEUP}_multiGeometry_tracks.root -o $outputPath/hits_ttbar_pu${PILEUP}_multiGeometry_tracks${spacing}.root -t Tracks${spacing}

done
