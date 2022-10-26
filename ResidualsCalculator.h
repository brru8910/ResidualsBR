/**
  \file
  Analysis script for calculating residual corrections within TPCs. Chamber volumes
  are binned into 8 x 8 x 8 cm^3 voxels. The corrections are written to a 
  database-compatible file and applied by the residuals applicator.

  \author B. Rumberger
  \version $Id:    $
  \date 31 Aug 2020
*/

#include <map>
#include <cmath>
#include <iostream>
#include <TH3I.h>
#include <TFile.h>

#include <modutils/ResidualsTools.h>

typedef std::vector<std::vector<std::vector<float>>> ResidualsCorrections;
ResidualsCorrections fPreviousResiduals;

double fMinPotentialPointsRatio = 0.6;
double fMaxPotentialPointsRatio = 1.1;
double fMaxAverageResidual = 0.06*utl::cm;
double fMinEntries = 1;
double fMaxChi2 = 200.0;
double fMaxBX = 4.0;
double fMaxBY = 2.0;
std::map<det::TPCConst::EId,unsigned int> fMinClustersByTPCId;

/// Main function.
int main(int argc, char* argv[]);

/// Configuration file parsing function.
void ParseConfigFile(const std::string& configFile);

// Display usage.
void DisplayUsage()
{
  std::cerr << "\nUsage:\nResidualsCalculator -i rootFiles -o outputPrefix "
    "[--updateResiuals / -u] [--localTrackResiduals / -l] [--useLegacyBins]\n" << std::endl;
  exit(-1);
}
