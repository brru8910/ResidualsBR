/**
  \file
  Analysis script for calculating residual corrections within TPCs. Chamber volumes
  are binned into user-defined voxels. The corrections are written to a 
  database-compatible file and applied by the residuals applicator.

  \author B. Rumberger
  \version $Id:    $
  \date 5 Oct 2022
*/

#include"ResidualsCalculator.h"

#include <fwk/CentralConfig.h>
#include <det/Detector.h>
#include <det/TPC.h>

#include <TFile.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TPaletteAxis.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TPaveText.h>
#include <TPaletteAxis.h>

#include <fstream>

#include <boost/filesystem.hpp>

using namespace std;
using namespace modutils;
using namespace modutils::ResidualsTools;

/// Main function.
int main(int argc, char* argv[])
{
  int exitCode = 0;

  //Option for updating previously-calculated residuals,
  //to run multi-stage residual calibration.
  bool updateResiduals = false;

  //Option for calculating local track residuals only. Removes vertex
  //track cuts.
  bool localTrackResiduals = false;

  //Use legacy-style binning scheme for residuals.
  bool useLegacyResidualBins = false;

  const vector<string> argumentsVector(argv + 1, argv + argc);
  vector<string> filenamesVector;

  string configFilename = "Config.txt";
  string outputPrefix;
  for (auto it = argumentsVector.begin(), itEnd = argumentsVector.end(); 
       it != itEnd; ++it) {
    if (*it == string("-h") || *it == string("--help")) {
      DisplayUsage();
    }
    else if (*it == string("-o")) {
      if (next(it) == itEnd || next(it)->rfind("-",0) != string::npos) {
	cout << "[ERROR] No output prefix provided with argument -o!" << endl;
	DisplayUsage();
      }
      advance(it,1);
      outputPrefix = "Residuals-" + *it;
    }
    else if (*it == string("-c") || *it == string("--config")) {
      if (next(it) == itEnd || next(it)->rfind("-",0) != string::npos) {
	cout << "[ERROR] No config filename provided with argument -c!" << endl;
	DisplayUsage();
      }
      advance(it,1);
      configFilename = *it;
      cout << "[INFO] User-provided config file: " << configFilename << endl;
    }
    else if (*it == string("-u") || *it == string("--updateResiduals")) {
      updateResiduals = true;
      cout << "[INFO] Updating previously-calcuated residuals." << endl;
    }
    else if (*it == string("-l") || *it == string("--localTrackResiduals")) {
      cout << "[INFO] Calculating local track residuals. Dropping vertex track cuts." << endl;
      localTrackResiduals = true;
    }
    else if (*it == string("--useLegacyBins")) {
      cout << "[INFO] Using Legacy binning scheme." << endl;
      useLegacyResidualBins = true;
    }
    else if (*it == string("-i") || *it == string("--inputFiles")) {
      filenamesVector.assign(next(it),itEnd);
      break;
    }
    else {
      cout << "[ERROR] Invalid argument " << *it << "!" << endl;
      DisplayUsage();
    }
  }

  if (filenamesVector.size() == 0) {
    cout << "[ERROR] No input filenames provided!" << endl;
    DisplayUsage();
  }
  if (outputPrefix.size() == 0) {
    cout << "[ERROR] No output prefix provided!" << endl;
    DisplayUsage();
  }
  cout << "[INFO] Number of input files: " << filenamesVector.size()
       << ". Config file: " << configFilename 
       << ". Update previously-calculated residuals? " << updateResiduals << endl;
  
  //Name and create output folder.
  boost::filesystem::create_directory(outputPrefix);

  //Name output file.
  TString outputFilename = outputPrefix + "/averages-" + outputPrefix + ".root";
  //Open the file for filling and get the TTree.
  TFile* inputFile1 = new TFile(filenamesVector.front().c_str(),"READ");
  if (inputFile1->IsZombie()) {
    cout << "Error opening input file 1! Something went wrong." << endl;
    return 1;
  }
  //Create output file.
  TFile* outputFile = new TFile(outputFilename,"RECREATE");  

  //Parse configuration file. FIXME! Put cuts, etc in here.
  ParseConfigFile(configFilename);
  
  //Container for holding residuals and performing calculations.
  ResidualsAccumulator residualsAccumulator;

  //Get run number for updating detector.
  TTree* fRunNumberTree = (TTree*)inputFile1->Get("fRunNumber");
  int runNumber = 0;
  fRunNumberTree->SetBranchAddress("fRunNumber",&runNumber);
  if (fRunNumberTree->GetEntries() > 0)
    fRunNumberTree->GetEntry(0);

  //Get detector and event interfaces.
  fwk::CentralConfig::GetInstance("bootstrap.xml");
  det::Detector& detector  = det::Detector::GetInstance();
  const utl::TimeStamp dummyTime = utl::TimeStamp(1);
  detector.Update(dummyTime,runNumber);
  const det::TPC& tpc = detector.GetTPC();
  //Map of coordinate system pointers for caching.
  map<det::TPCConst::EId,map<int,utl::CoordinateSystemPtr> > sectorCoordinateSystems;
    
  //Loop through chambers and initialize container.
  for (auto it = tpc.ChambersBegin(), itEnd = tpc.ChambersEnd(); it != itEnd; ++it) {
    const det::TPCChamber& chamber = *it;
    for (unsigned int sectorId = 1; sectorId <= chamber.GetNSectors(); ++sectorId) {
      const det::TPCSector& sector = chamber.GetSector(sectorId);
      //Sector local coordinate system is used.
      const utl::CoordinateSystemPtr& sectorCS = sector.GetComponentCoordinateSystem();
      sectorCoordinateSystems[chamber.GetId()][sectorId] = sectorCS;
    }
  }

  //Get previously-calculated residuals, if updating.
  if (updateResiduals) {
    fPreviousResiduals = tpc.GetResCorrArrays();
    const string residualsPath = tpc.GetResCorrString();
    cout << "[INFO] Updating residuals using files found in " << residualsPath << endl;
  }

  //Construct residual correction bins within each TPC.
  if (useLegacyResidualBins)  {
    cout << "[INFO] Using Legacy binning scheme." << endl;
    ResidualsBinTool binTool;
    binTool.CreateLegacyBinningHistograms();
    binTool.SetBinningHistograms(residualsAccumulator);
  }
  else {
    cout << "[INFO] Using SHINE-generated binning scheme." << endl;
    ResidualsBinTool binTool;
    binTool.CreateBinningHistograms(tpc);
    binTool.SetBinningHistograms(residualsAccumulator);
  }
    

  //Loop over input files. Give progress percentage.
  double filesProcessed = 0;
  for (auto fileIt = filenamesVector.begin(), fileEnd = filenamesVector.end();
       fileIt != fileEnd; ++fileIt) {
    ++filesProcessed;
    const double progressPercentage = round(1000*(filesProcessed - 1)/filenamesVector.size()/10.);
    cout << "[INFO] Processing file " << filesProcessed
         << " / " << filenamesVector.size() 
	 << " (" << progressPercentage << "% complete)." << endl;

    //Open the file for filling and get the TTree.
    const string& filename = *fileIt;
    TFile* inputFile = new TFile(filename.c_str(),"READ");
    if (inputFile == NULL || inputFile->IsZombie()) {
      cout << "[ERROR] Error opening input file! Something went wrong." << endl;
      continue;
    }

    if (inputFile->GetNkeys() == 0) {
      cout << "[INFO] " << filename << " has no keys. Skipping." << endl;
      continue;
    }
      
    //Get the trees. Set branch address.
    TTree* fResidualsData = (TTree*)inputFile->Get("fResiduals");
    modutils::ResidualsTools::ResidualEntry entry;
    fResidualsData->SetBranchAddress("fResidualEntry",&entry);
  
    //For printing progress.
    double previousPercentage = 0;
    long unsigned int nTotalEntries = fResidualsData->GetEntries();
  
    //Loop through all data to calculate total average for sectors.
    for (long unsigned int i = 0; i < nTotalEntries; ++i) {
      fResidualsData->GetEntry(i);
    
      //Cuts. For vertex tracks only.
      if (!localTrackResiduals) {
	if (!entry.fIsRST)
	  continue;
	if (entry.fPotentialPointsRatio < fMinPotentialPointsRatio ||
	    entry.fPotentialPointsRatio > fMaxPotentialPointsRatio  )
	  continue;
	if (entry.fAverageResidual > fMaxAverageResidual)
	  continue;
	if (entry.fChi2 > fMaxChi2)
	  continue;
	if (entry.fBX > fMaxBX ||
	    entry.fBY > fMaxBY  )
	  continue;
      }

      //Cuts that apply to local tracks.
      if (localTrackResiduals) {
	if (entry.fAverageResidual > fMaxAverageResidual)
	  continue;
	if (entry.fNClusters < fMinClustersByTPCId[(det::TPCConst::EId)entry.fTPCId])
	  continue;
      }
    
      //Residuals are stored and applied in the local TPC coordinate systems.
      const det::TPCConst::EId tpcId = (det::TPCConst::EId)entry.fTPCId;
      const int sectorId = entry.fSectorId;
      const utl::Point position(entry.fX,
				entry.fY,
				entry.fZ);
      cout << det::TPCConst::GetName(tpcId) << " point " << position << endl;
      const utl::Point
	localPosition(position.GetCoordinates(sectorCoordinateSystems[tpcId][sectorId]).get<0>(),
		      position.GetCoordinates(sectorCoordinateSystems[tpcId][sectorId]).get<1>(),
		      entry.fPadrowId);
      const utl::Vector residual(entry.fXResidual,
				 entry.fYResidual,
				 0);

      //Add to accumulator.
      residualsAccumulator.AddEntry(tpcId,
				    entry.fSectorId,
				    entry.fPadrowId,
				    localPosition,
				    residual);

      //Let the user know what's going on.
      int progressPercentage = 100*i/fResidualsData->GetEntries();
      if (progressPercentage % 5 == 0 &&
	  progressPercentage != previousPercentage) {
	ostringstream info;
	info << "[INFO] Calculating...." << progressPercentage << "% finished.";
	INFO(info);
	previousPercentage = progressPercentage;
      }
    } //End TTree loop.
    inputFile->Close();
  } //End filename loop.
  
  ostringstream info;
  info << "[INFO] Resduals calculation complete. Writing corrections and QA files.";
  INFO(info);
  
  //TTree for storing results.
  TTree* fResultTree = new TTree("fResultTree","Analysis Results");
  double fX;
  fResultTree->Branch("fX",&fX);
  double fY;
  fResultTree->Branch("fY",&fY);
  double fZ;
  fResultTree->Branch("fZ",&fZ);
  double fLocalX;
  fResultTree->Branch("fLocalX",&fLocalX);
  double fLocalY;
  fResultTree->Branch("fLocalY",&fLocalY);
  double fLocalZ;
  fResultTree->Branch("fLocalZ",&fLocalZ);
  double fXResidual;
  fResultTree->Branch("fXResidual",&fXResidual);
  double fYResidual;
  fResultTree->Branch("fYResidual",&fYResidual);
  double fTPCId;
  fResultTree->Branch("fTPCId",&fTPCId);
  int fSectorId;
  fResultTree->Branch("fSectorId",&fSectorId);
  int fPadrowId;
  fResultTree->Branch("fPadrowId",&fPadrowId);
  unsigned int fAccumulatorIndex;
  fResultTree->Branch("fAccumulatorIndex",&fAccumulatorIndex);

  //QA histograms.
  const unsigned int nZPlanes = tpc.GetZPlaneIdRange().GetEnd();
  const double maxPlotResidual = 0.2*utl::cm;
  map<det::TPCConst::EId,map<unsigned int,TProfile2D*> > xResidualPlots;
  map<det::TPCConst::EId,map<unsigned int,TProfile2D*> > yResidualPlots;
  map<det::TPCConst::EId,TH1D*> xResiduals1D;
  map<det::TPCConst::EId,TH1D*> yResiduals1D;
  TProfile* x1DResidualPlots =
    new TProfile("x1DResiduals","X Residuals;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DResidualPlots =
    new TProfile("y1DResiduals","Y Residuals;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);

  TProfile* x1DResidualPlotsYMinus10XNegative =
    new TProfile("x1DResidualsYMinus10XNegative",
		 "X Residuals, Y < -10 cm & X < 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DResidualPlotsYMinus10XNegative =
    new TProfile("y1DResidualsYMinus10XNegative",
		 "Y Residuals, Y < -10 cm & X < 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);
  TProfile* x1DResidualPlotsYMinus10XPositive =
    new TProfile("x1DResidualsYMinus10XPositive",
		 "X Residuals, Y < -10 cm & X > 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DResidualPlotsYMinus10XPositive =
    new TProfile("y1DResidualsYMinus10XPositive",
		 "Y Residuals, Y < -10 cm & X > 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);

  TProfile* x1DResidualPlotsYPlus10XNegative =
    new TProfile("x1DResidualsYPlus10XNegative",
		 "X Residuals, Y > 10 cm & X < 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DResidualPlotsYPlus10XNegative =
    new TProfile("y1DResidualsYPlus10XNegative",
		 "Y Residuals, Y > 10 cm & X < 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);
  TProfile* x1DResidualPlotsYPlus10XPositive =
    new TProfile("x1DResidualsYPlus10XPositive",
		 "X Residuals, Y > 10 cm & X > 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DResidualPlotsYPlus10XPositive =
    new TProfile("y1DResidualsYPlus10XPositive",
		 "Y Residuals, Y > 10 cm & X > 0 cm;Padrow ID;Residual [cm]",nZPlanes,0,nZPlanes);

  TProfile* x1DResidualPlotsYCentral10XNegative =
    new TProfile("x1DResidualsYCentral10XNegative",
		 "X Residuals, -10 < Y < 10 cm & X < 0 cm;Padrow ID;Residual [cm]",
		 nZPlanes,0,nZPlanes);
  TProfile* y1DResidualPlotsYCentral10XNegative =
    new TProfile("y1DResidualsYCentral10XNegative",
		 "Y Residuals, -10 < Y > 10 cm & X < 0 cm;Padrow ID;Residual [cm]",
		 nZPlanes,0,nZPlanes);
  TProfile* x1DResidualPlotsYCentral10XPositive =
    new TProfile("x1DResidualsYCentral10XPositive",
		 "X Residuals, -10 < Y < 10 cm & X > 0 cm;Padrow ID;Residual [cm]",
		 nZPlanes,0,nZPlanes);
  TProfile* y1DResidualPlotsYCentral10XPositive =
    new TProfile("y1DResidualsYCentral10XPositive",
		 "Y Residuals, -10 < Y < 10 cm & X > 0 cm;Padrow ID;Residual [cm]",
		 nZPlanes,0,nZPlanes);

  TProfile* x1DCorrectionPlots =
    new TProfile("x1DCorrections","X Corrections;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DCorrectionPlots =
    new TProfile("y1DCorrections","Y Corrections;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);

  TProfile* x1DCorrectionPlotsYMinus10XNegative =
    new TProfile("x1DCorrectionsYMinus10XNegative",
		 "X Corrections, Y < -10 cm & X < 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DCorrectionPlotsYMinus10XNegative =
    new TProfile("y1DCorrectionsYMinus10XNegative",
		 "Y Corrections, Y < -10 cm & X < 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);
  TProfile* x1DCorrectionPlotsYMinus10XPositive =
    new TProfile("x1DCorrectionsYMinus10XPositive",
		 "X Corrections, Y < -10 cm & X > 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DCorrectionPlotsYMinus10XPositive =
    new TProfile("y1DCorrectionsYMinus10XPositive",
		 "Y Corrections, Y < -10 cm & X > 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);

  TProfile* x1DCorrectionPlotsYPlus10XNegative =
    new TProfile("x1DCorrectionsYPlus10XNegative",
		 "X Corrections, Y > 10 cm & X < 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DCorrectionPlotsYPlus10XNegative =
    new TProfile("y1DCorrectionsYPlus10XNegative",
		 "Y Corrections, Y > 10 cm & X < 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);
  TProfile* x1DCorrectionPlotsYPlus10XPositive =
    new TProfile("x1DCorrectionsYPlus10XPositive",
		 "X Corrections, Y > 10 cm & X > 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);
  TProfile* y1DCorrectionPlotsYPlus10XPositive =
    new TProfile("y1DCorrectionsYPlus10XPositive",
		 "Y Corrections, Y > 10 cm & X > 0 cm;Padrow ID;Correction [cm]",nZPlanes,0,nZPlanes);

  TProfile* x1DCorrectionPlotsYCentral10XNegative =
    new TProfile("x1DCorrectionsYCentral10XNegative",
		 "X Corrections, -10 < Y < 10 cm & X < 0 cm;Padrow ID;Correction [cm]",
		 nZPlanes,0,nZPlanes);
  TProfile* y1DCorrectionPlotsYCentral10XNegative =
    new TProfile("y1DCorrectionsYCentral10XNegative",
		 "Y Corrections, -10 < Y > 10 cm & X < 0 cm;Padrow ID;Correction [cm]",
		 nZPlanes,0,nZPlanes);
  TProfile* x1DCorrectionPlotsYCentral10XPositive =
    new TProfile("x1DCorrectionsYCentral10XPositive",
		 "X Corrections, -10 < Y < 10 cm & X > 0 cm;Padrow ID;Correction [cm]",
		 nZPlanes,0,nZPlanes);
  TProfile* y1DCorrectionPlotsYCentral10XPositive =
    new TProfile("y1DCorrectionsYCentral10XPositive",
		 "Y Corrections, -10 < Y < 10 cm & X > 0 cm;Padrow ID;Correction [cm]",
		 nZPlanes,0,nZPlanes);
  
  
  //Loop through averagers to calculate residuals.
  for (auto chamberIt = residualsAccumulator.Begin(), chamberEnd = residualsAccumulator.End();
       chamberIt != chamberEnd; ++chamberIt) {
    const det::TPCConst::EId tpcId = chamberIt->first;

    const string& tpcName = det::TPCConst::GetName(tpcId);
    xResiduals1D[tpcId] = new TH1D(Form("%sXResiduals",tpcName.data()),
				   Form("%s X-Residuals;Residual [cm];Entries",tpcName.data()),
				   500,-1*maxPlotResidual,maxPlotResidual);
    yResiduals1D[tpcId] = new TH1D(Form("%sYResiduals",tpcName.data()),
				   Form("%s Y-Residuals;Residual [cm];Entries",tpcName.data()),
				   500,-1*maxPlotResidual,maxPlotResidual);
    
    //Make sure detector has valid histograms.
    if (!residualsAccumulator.HasDetectorHistograms(tpcId))
      continue;

    //Initialize QA histograms.
    const det::TPCChamber& chamber = tpc.GetChamber(tpcId);
    for (unsigned int i = 1; i <= chamber.GetNSectors(); ++i) {
      TString titleStringX = det::TPCConst::GetName(tpcId) + Form("Sector%iXResiduals",i);
      TString titleStringY = det::TPCConst::GetName(tpcId) + Form("Sector%iYResiduals",i);

      TH3I histogram = residualsAccumulator.GetSectorHistogram(tpcId,i);

      const TArrayD* xBins = histogram.GetXaxis()->GetXbins();
      // const TArrayD* yBins = residualsAccumulator.GetSectorYBins(tpcId,i);
      const TArrayD* zBins = histogram.GetZaxis()->GetXbins();
      xResidualPlots[tpcId][i] = new TProfile2D(titleStringX,
						"",
						zBins->GetSize()-1,
						zBins->GetArray(),
						xBins->GetSize()-1,
						xBins->GetArray());
      yResidualPlots[tpcId][i] = new TProfile2D(titleStringY,
						"",
						zBins->GetSize()-1,
						zBins->GetArray(),
						xBins->GetSize()-1,
						xBins->GetArray());
    }
    
    
    //Create corrections file.
    //Naming and indexing conventions are put here for compatibility
    //with the Legacy corrections reader. They can be changed if the
    //manager that reads Residuals files is updated.
    string correctionsFilename = outputPrefix;
    if (tpcId == det::TPCConst::eVTPC1) 
      correctionsFilename += "/vt1.corr";
    else if (tpcId == det::TPCConst::eVTPC2) 
      correctionsFilename += "/vt2.corr";
    else if (tpcId == det::TPCConst::eGTPC) 
      correctionsFilename += "/gt.corr";
    else if (tpcId == det::TPCConst::eMTPCL) 
      correctionsFilename += "/mtl.corr";
    else if (tpcId == det::TPCConst::eMTPCR) 
      correctionsFilename += "/mtr.corr";
    else if (tpcId == det::TPCConst::eFTPC1) 
      correctionsFilename += "/ft1.corr";
    else if (tpcId == det::TPCConst::eFTPC2) 
      correctionsFilename += "/ft2.corr";
    else if (tpcId == det::TPCConst::eFTPC3) 
      correctionsFilename += "/ft3.corr";
    else {
      ostringstream msg;
      msg << "[ERROR] " << GetName(tpcId)
	  << " does not have a correction file name "
	"/ index defined. Not creating a file for this detector!";
      WARNING(msg);
      continue;
    }

    ofstream correctionsFileStream;
    correctionsFileStream.open(correctionsFilename+"-NotUpdated");

    ofstream updatedCorrectionsFileStream;
    updatedCorrectionsFileStream.open(correctionsFilename);

    const int maxIndex = residualsAccumulator.GetMaxCorrectionIndex(tpcId);
    
    //Loop through bins.
    for (int correctionIndex = 0; correctionIndex <= maxIndex; ++correctionIndex) {
      
      //Get previously-calculated residual.
      const double previousXCorrection = (updateResiduals) ?
	fPreviousResiduals.at(tpcId-1).at(0).at(correctionIndex) : 0;
      const double previousYCorrection = (updateResiduals) ?
	fPreviousResiduals.at(tpcId-1).at(1).at(correctionIndex) : 0;
      
      const double nEntries = residualsAccumulator.GetNumberOfEntries(tpcId,correctionIndex);
      
      double xCorrection = 0;
      if (!isnan(residualsAccumulator.GetXResidual(tpcId,correctionIndex)) &&
	  nEntries >= fMinEntries)
        xCorrection = residualsAccumulator.GetXResidual(tpcId,correctionIndex);
      double yCorrection = 0;
      if (!isnan(residualsAccumulator.GetYResidual(tpcId,correctionIndex)) &&
	  nEntries >= fMinEntries)
        yCorrection = residualsAccumulator.GetYResidual(tpcId,correctionIndex);

      correctionsFileStream << xCorrection << " "
			    << to_string(correctionIndex) << " 0 "
			    << to_string(tpcId) << "\n";
      correctionsFileStream << yCorrection << " "
			    << to_string(correctionIndex) << " 1 "
			    << to_string(tpcId) << "\n";
      
      updatedCorrectionsFileStream << xCorrection + previousXCorrection << " "
				   << to_string(correctionIndex) << " 0 "
				   << to_string(tpcId) << "\n";
      updatedCorrectionsFileStream << yCorrection + previousYCorrection << " "
				   << to_string(correctionIndex) << " 1 "
				   << to_string(tpcId) << "\n";
      
      if (xCorrection == 0 && yCorrection == 0)
	continue;
      
      if (!residualsAccumulator.HasIndex(tpcId,correctionIndex)) {
	cout << "ERROR! No correction index " << correctionIndex
	     << " for " << det::TPCConst::GetName(tpcId) << endl;
	continue;
      }
      if (residualsAccumulator.HasIndex(tpcId,correctionIndex)) {

	fSectorId = residualsAccumulator.GetSector(tpcId,correctionIndex);
	fPadrowId = residualsAccumulator.GetPadrow(tpcId,correctionIndex);

	const utl::Point& localPosition =
	  residualsAccumulator.GetBinCenterLocal(tpcId,correctionIndex);

	const det::TPCPadrow& padrow =
	  tpc.GetChamber(tpcId).GetSector(fSectorId).GetPadrow(fPadrowId);
	const double zPlaneId = padrow.GetZPlaneId();
	const double padrowLocalZ = padrow.GetLocalZ();

	const utl::Point position(localPosition.GetX(),
				  localPosition.GetY(),
				  padrowLocalZ,
				  sectorCoordinateSystems[tpcId][fSectorId]);
	fX = position.GetX();
	fY = position.GetY();
	fZ = position.GetZ();
	fLocalX = localPosition.GetX();
	fLocalY = localPosition.GetY();
	fLocalZ = localPosition.GetZ();
	fXResidual = xCorrection;
	fYResidual = yCorrection;
	fTPCId = tpcId;
	fAccumulatorIndex = correctionIndex;
	
	// if (tpcId == det::TPCConst::eFTPC1 || tpcId == det::TPCConst::eFTPC2)
	//   cout << det::TPCConst::GetName(tpcId)
	//        << " index " << correctionIndex
	//        << " Bin center: " << position
	//        << ". Local bin center: " << localPosition
	//        << ". correction: (" << xCorrection << "," << yCorrection << ")" << endl;

	
	fResultTree->Fill();
	
	xResiduals1D[(det::TPCConst::EId)fTPCId]->Fill(xCorrection);
	yResiduals1D[(det::TPCConst::EId)fTPCId]->Fill(yCorrection);
	
	xResidualPlots[(det::TPCConst::EId)fTPCId][fSectorId]->Fill(fPadrowId,
								    localPosition.GetX(),
								    xCorrection);
	yResidualPlots[(det::TPCConst::EId)fTPCId][fSectorId]->Fill(fPadrowId,
								    localPosition.GetX(),
								    yCorrection);
	
	x1DResidualPlots->Fill(zPlaneId,xCorrection);
	y1DResidualPlots->Fill(zPlaneId,yCorrection);	

	x1DCorrectionPlots->Fill(zPlaneId,xCorrection + previousXCorrection);
	y1DCorrectionPlots->Fill(zPlaneId,yCorrection + previousYCorrection);	

	if (fX < 0) {
	  if (fY < -10) {
	    x1DResidualPlotsYMinus10XNegative->Fill(zPlaneId,xCorrection);
	    y1DResidualPlotsYMinus10XNegative->Fill(zPlaneId,yCorrection);
	    x1DCorrectionPlotsYMinus10XNegative->Fill(zPlaneId,xCorrection + previousXCorrection);
	    y1DCorrectionPlotsYMinus10XNegative->Fill(zPlaneId,yCorrection + previousYCorrection);
	  }
	  if (fY >= -10 && fY <= 10) {
	    x1DResidualPlotsYCentral10XNegative->Fill(zPlaneId,xCorrection);
	    y1DResidualPlotsYCentral10XNegative->Fill(zPlaneId,yCorrection);
	    x1DCorrectionPlotsYCentral10XNegative->Fill(zPlaneId,xCorrection + previousXCorrection);
	    y1DCorrectionPlotsYCentral10XNegative->Fill(zPlaneId,yCorrection + previousYCorrection);
	  }
	  if (fY > 10) {
	    x1DResidualPlotsYPlus10XNegative->Fill(zPlaneId,xCorrection);
	    y1DResidualPlotsYPlus10XNegative->Fill(zPlaneId,yCorrection);
	    x1DCorrectionPlotsYPlus10XNegative->Fill(zPlaneId,xCorrection + previousXCorrection);
	    y1DCorrectionPlotsYPlus10XNegative->Fill(zPlaneId,yCorrection + previousYCorrection);
 	  }
	}
	if (fX > 0) {
	  if (fY < -10) {
	    x1DResidualPlotsYMinus10XPositive->Fill(zPlaneId,xCorrection);
	    y1DResidualPlotsYMinus10XPositive->Fill(zPlaneId,yCorrection);
	    x1DCorrectionPlotsYMinus10XPositive->Fill(zPlaneId,xCorrection + previousXCorrection);
	    y1DCorrectionPlotsYMinus10XPositive->Fill(zPlaneId,yCorrection + previousYCorrection);
	  }
	  if (fY >= -10 && fY <= 10) {
	    x1DResidualPlotsYCentral10XPositive->Fill(zPlaneId,xCorrection);
	    y1DResidualPlotsYCentral10XPositive->Fill(zPlaneId,yCorrection);
	    x1DCorrectionPlotsYCentral10XPositive->Fill(zPlaneId,xCorrection + previousXCorrection);
	    y1DCorrectionPlotsYCentral10XPositive->Fill(zPlaneId,yCorrection + previousYCorrection);
	  }
	  if (fY > 10) {
	    x1DResidualPlotsYPlus10XPositive->Fill(zPlaneId,xCorrection);
	    y1DResidualPlotsYPlus10XPositive->Fill(zPlaneId,yCorrection);
	    x1DCorrectionPlotsYPlus10XPositive->Fill(zPlaneId,xCorrection + previousXCorrection);
	    y1DCorrectionPlotsYPlus10XPositive->Fill(zPlaneId,yCorrection + previousYCorrection);
	  }
	}

      }
      
    }
    correctionsFileStream.close();
  }

  ostringstream finish;
  finish << "[INFO] Correction and QA file creation complete. Thanks!";
  INFO(finish);

  TString residualsPDFName = outputPrefix + "/" + outputPrefix + ".pdf";
  TString residualsOpenString = residualsPDFName + "[";
  TString residualsCloseString = residualsPDFName + "]";
  TCanvas dummy;
  dummy.SaveAs(residualsOpenString);

  gStyle->SetOptStat(0);
  
  //Create z-scale palette using dummy TH2D.
  TH2D* dummy2D = new TH2D("dummy","dummy",100,0,1,100,0,1);
  dummy2D->Fill(0.1,0.1,-1*maxPlotResidual);
  dummy2D->Fill(0.9,0.9,maxPlotResidual);
  dummy2D->GetZaxis()->SetLabelSize(0.02);
  dummy2D->Draw("COLZ");
  // dummy.SaveAs(residualsPDFName);
  dummy.Update();
  TPaletteAxis* palette =
    (TPaletteAxis*)dummy2D->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.1);
  palette->SetX2NDC(0.11);
  palette->SetY1NDC(0.7);
  palette->SetY2NDC(0.9);

  TPaveText* label1 = new TPaveText(.15,.88,.19,.92,"xResiduals");
  label1->AddText("[cm]");
  label1->SetFillColor(0);
  TPaveText* label2 = new TPaveText(.15,.68,.19,.72,"xResiduals");
  label2->AddText("[cm]");
  label2->SetFillColor(0);

  
  const double minZ = -500*utl::cm;
  const double maxZ =  910*utl::cm;
  const double minX = -400*utl::cm;
  const double maxX =  400*utl::cm;

  const double zRange = maxZ - minZ;
  const double xRange = maxX - minX;
  
  TCanvas xOverviewCanvas;
  TPaveText* xTitle = new TPaveText(.3,.95,.6,.99,"xResiduals");
  xTitle->AddText("X-Residuals");
  xTitle->Draw();
  
  for (auto it = xResidualPlots.begin(), itEnd = xResidualPlots.end(); it != itEnd; ++it) {
    map<unsigned int, TProfile2D*> sectorPlots = it->second;
    for (auto sectorIt = sectorPlots.begin(), sectorEnd = sectorPlots.end();
	 sectorIt != sectorEnd; ++sectorIt) {
      
      //Determine pad coordinates in canvas.
      const det::TPCSector& sector = tpc.GetChamber(it->first).GetSector(sectorIt->first);
      const det::TPCPadrow& firstPadrow = sector.GetPadrow(1);
      const det::TPCPadrow& lastPadrow = sector.GetPadrow(sector.GetNPadrows());
      const double minSectorX = firstPadrow.GetPadPosition(1).GetX();
      const double maxSectorX = firstPadrow.GetPadPosition(firstPadrow.GetNPads()).GetX();
      const double minSectorZ = firstPadrow.GetCenterPosition().GetZ();
      const double maxSectorZ = lastPadrow.GetCenterPosition().GetZ();
      
      TProfile2D* xPlot = sectorIt->second;
      
      xPlot->SetMinimum(-1*maxPlotResidual);
      xPlot->SetMaximum(maxPlotResidual);

      //Translate to NDC [0,1] coordinates in ROOT's canvas. Global Z is plotted along the x-axis.
      const double ndcX0 = (minSectorZ - minZ)/zRange;
      const double ndcX1 = (maxSectorZ - minZ)/zRange;
      const double ndcY0 = (minSectorX - minX)/xRange;
      const double ndcY1 = (maxSectorX - minX)/xRange;

      //Pad constructor coordinate system: TPad(x0,y0,x1,y1) (in NDC)
      xOverviewCanvas.cd();
      TPad* pad = new TPad(xPlot->GetName(),xPlot->GetName(),ndcX0,ndcY0,ndcX1,ndcY1);
      pad->Draw();
      pad->cd();
      xPlot->Draw("COLA");
    }    
  } 
  xOverviewCanvas.cd();
  palette->Draw();
  label1->Draw();
  label2->Draw();
  xOverviewCanvas.SaveAs(residualsPDFName);

  TCanvas yOverviewCanvas;
  TPaveText* yTitle = new TPaveText(.3,.95,.6,.99,"yResiduals");
  yTitle->AddText("Y-Residuals");
  yTitle->Draw();
  
  for (auto it = yResidualPlots.begin(), itEnd = yResidualPlots.end(); it != itEnd; ++it) {
    map<unsigned int, TProfile2D*> sectorPlots = it->second;
    for (auto sectorIt = sectorPlots.begin(), sectorEnd = sectorPlots.end();
  	 sectorIt != sectorEnd; ++sectorIt) {

      //Determine pad coordinates in canvas.
      const det::TPCSector& sector = tpc.GetChamber(it->first).GetSector(sectorIt->first);
      const det::TPCPadrow& firstPadrow = sector.GetPadrow(1);
      const det::TPCPadrow& lastPadrow = sector.GetPadrow(sector.GetNPadrows());
      const double minSectorX = firstPadrow.GetPadPosition(1).GetX();
      const double maxSectorX = firstPadrow.GetPadPosition(firstPadrow.GetNPads()).GetX();
      const double minSectorZ = firstPadrow.GetCenterPosition().GetZ();
      const double maxSectorZ = lastPadrow.GetCenterPosition().GetZ();
      
      TProfile2D* yPlot = sectorIt->second;
      
      yPlot->SetMinimum(-1*maxPlotResidual);
      yPlot->SetMaximum(maxPlotResidual);

      //Translate to NDC [0,1] coordinates in ROOT's canvas. Global Z is plotted along the x-axis.
      const double ndcX0 = (minSectorZ - minZ)/zRange;
      const double ndcX1 = (maxSectorZ - minZ)/zRange;
      const double ndcY0 = (minSectorX - minX)/xRange;
      const double ndcY1 = (maxSectorX - minX)/xRange;

      //Pad constructor coordinate system: TPad(x0,y0,x1,y1) (in NDC)
      yOverviewCanvas.cd();
      TPad* pad = new TPad(yPlot->GetName(),yPlot->GetName(),ndcX0,ndcY0,ndcX1,ndcY1);
      pad->Draw();
      pad->cd();
      yPlot->Draw("COLA");
    }    
  }
  
  yOverviewCanvas.cd();
  palette->Draw();
  label1->Draw();
  label2->Draw();
  yOverviewCanvas.SaveAs(residualsPDFName);

  TCanvas oneDimensionalCanvas;
  oneDimensionalCanvas.cd();
  x1DResidualPlots->SetMinimum(-1*maxPlotResidual);
  x1DResidualPlots->SetMaximum(maxPlotResidual);
  x1DResidualPlots->Draw();
  x1DResidualPlots->Write();
  oneDimensionalCanvas.Write();
  oneDimensionalCanvas.SaveAs(residualsPDFName);
  y1DResidualPlots->SetMinimum(-1*maxPlotResidual);
  y1DResidualPlots->SetMaximum(maxPlotResidual);
  y1DResidualPlots->Draw();
  y1DResidualPlots->Write();
  oneDimensionalCanvas.SaveAs(residualsPDFName);
  oneDimensionalCanvas.Write();

  TCanvas oneDimensionalSliceCanvasX;
  oneDimensionalSliceCanvasX.Divide(2,3);
  oneDimensionalSliceCanvasX.cd(1);
  x1DResidualPlotsYMinus10XNegative->SetMinimum(-1*maxPlotResidual);
  x1DResidualPlotsYMinus10XNegative->SetMaximum(maxPlotResidual);
  x1DResidualPlotsYMinus10XNegative->Draw();
  x1DResidualPlotsYMinus10XNegative->Write();
  oneDimensionalSliceCanvasX.cd(2);
  x1DResidualPlotsYMinus10XPositive->SetMinimum(-1*maxPlotResidual);
  x1DResidualPlotsYMinus10XPositive->SetMaximum(maxPlotResidual);
  x1DResidualPlotsYMinus10XPositive->Draw();
  x1DResidualPlotsYMinus10XPositive->Write();
  oneDimensionalSliceCanvasX.cd(3);
  x1DResidualPlotsYCentral10XNegative->SetMinimum(-1*maxPlotResidual);
  x1DResidualPlotsYCentral10XNegative->SetMaximum(maxPlotResidual);
  x1DResidualPlotsYCentral10XNegative->Draw();
  x1DResidualPlotsYCentral10XNegative->Write();
  oneDimensionalSliceCanvasX.cd(4);
  x1DResidualPlotsYCentral10XPositive->SetMinimum(-1*maxPlotResidual);
  x1DResidualPlotsYCentral10XPositive->SetMaximum(maxPlotResidual);
  x1DResidualPlotsYCentral10XPositive->Draw();
  x1DResidualPlotsYCentral10XPositive->Write();
  oneDimensionalSliceCanvasX.cd(5);
  x1DResidualPlotsYPlus10XNegative->SetMinimum(-1*maxPlotResidual);
  x1DResidualPlotsYPlus10XNegative->SetMaximum(maxPlotResidual);
  x1DResidualPlotsYPlus10XNegative->Draw();
  x1DResidualPlotsYPlus10XNegative->Write();
  oneDimensionalSliceCanvasX.cd(6);
  x1DResidualPlotsYPlus10XPositive->SetMinimum(-1*maxPlotResidual);
  x1DResidualPlotsYPlus10XPositive->SetMaximum(maxPlotResidual);
  x1DResidualPlotsYPlus10XPositive->Draw();
  x1DResidualPlotsYPlus10XPositive->Write();
  oneDimensionalSliceCanvasX.SaveAs(residualsPDFName);
  oneDimensionalSliceCanvasX.Write();
  
  TCanvas oneDimensionalSliceCanvasY;
  oneDimensionalSliceCanvasY.Divide(2,3);
  oneDimensionalSliceCanvasY.cd(1);
  y1DResidualPlotsYMinus10XNegative->SetMinimum(-1*maxPlotResidual);
  y1DResidualPlotsYMinus10XNegative->SetMaximum(maxPlotResidual);
  y1DResidualPlotsYMinus10XNegative->Draw();
  y1DResidualPlotsYMinus10XNegative->Write();
  oneDimensionalSliceCanvasY.cd(2);
  y1DResidualPlotsYMinus10XPositive->SetMinimum(-1*maxPlotResidual);
  y1DResidualPlotsYMinus10XPositive->SetMaximum(maxPlotResidual);
  y1DResidualPlotsYMinus10XPositive->Draw();
  y1DResidualPlotsYMinus10XPositive->Write();
  oneDimensionalSliceCanvasY.cd(3);
  y1DResidualPlotsYCentral10XNegative->SetMinimum(-1*maxPlotResidual);
  y1DResidualPlotsYCentral10XNegative->SetMaximum(maxPlotResidual);
  y1DResidualPlotsYCentral10XNegative->Draw();
  y1DResidualPlotsYCentral10XNegative->Write();
  oneDimensionalSliceCanvasY.cd(4);
  y1DResidualPlotsYCentral10XPositive->SetMinimum(-1*maxPlotResidual);
  y1DResidualPlotsYCentral10XPositive->SetMaximum(maxPlotResidual);
  y1DResidualPlotsYCentral10XPositive->Draw();
  y1DResidualPlotsYCentral10XPositive->Write();
  oneDimensionalSliceCanvasY.cd(5);
  y1DResidualPlotsYPlus10XNegative->SetMinimum(-1*maxPlotResidual);
  y1DResidualPlotsYPlus10XNegative->SetMaximum(maxPlotResidual);
  y1DResidualPlotsYPlus10XNegative->Draw();
  y1DResidualPlotsYPlus10XNegative->Write();
  oneDimensionalSliceCanvasY.cd(6);
  y1DResidualPlotsYPlus10XPositive->SetMinimum(-1*maxPlotResidual);
  y1DResidualPlotsYPlus10XPositive->SetMaximum(maxPlotResidual);
  y1DResidualPlotsYPlus10XPositive->Draw();
  y1DResidualPlotsYPlus10XPositive->Write();
  oneDimensionalSliceCanvasY.SaveAs(residualsPDFName);
  oneDimensionalSliceCanvasY.Write();

  TCanvas fullTPCResidualsCanvas;
  fullTPCResidualsCanvas.Divide(4,2);
  int padId = 1;
  for (auto it = xResiduals1D.begin(),
	 itEnd = xResiduals1D.end();
       it != itEnd; ++it) {
    if (it->first == det::TPCConst::eLMPDJU ||
	it->first == det::TPCConst::eLMPDJD ||
	it->first == det::TPCConst::eLMPDSU ||
	it->first == det::TPCConst::eLMPDSD ||
	it->first == det::TPCConst::eLMPDP   )
      continue;
    fullTPCResidualsCanvas.cd(padId);  
    it->second->Draw();
    ++padId;
  }
  fullTPCResidualsCanvas.Write();
  fullTPCResidualsCanvas.SaveAs(residualsPDFName);
  padId = 1;
  for (auto it = yResiduals1D.begin(),
	 itEnd = yResiduals1D.end();
       it != itEnd; ++it) {
    if (it->first == det::TPCConst::eLMPDJU ||
	it->first == det::TPCConst::eLMPDJD ||
	it->first == det::TPCConst::eLMPDSU ||
	it->first == det::TPCConst::eLMPDSD ||
	it->first == det::TPCConst::eLMPDP   )
      continue;
    fullTPCResidualsCanvas.cd(padId);  
    it->second->Draw();
    ++padId;
  }
  fullTPCResidualsCanvas.Write();
  fullTPCResidualsCanvas.SaveAs(residualsPDFName);
  
  TCanvas correctionCanvas;
  correctionCanvas.cd();
  x1DCorrectionPlots->SetMinimum(-1*maxPlotResidual);
  x1DCorrectionPlots->SetMaximum(maxPlotResidual);
  x1DCorrectionPlots->Draw();
  x1DCorrectionPlots->Write();
  correctionCanvas.Write();
  correctionCanvas.SaveAs(residualsPDFName);
  y1DCorrectionPlots->SetMinimum(-1*maxPlotResidual);
  y1DCorrectionPlots->SetMaximum(maxPlotResidual);
  y1DCorrectionPlots->Draw();
  y1DCorrectionPlots->Write();
  correctionCanvas.SaveAs(residualsPDFName);
  correctionCanvas.Write();

  TCanvas correctionSliceCanvasX;
  correctionSliceCanvasX.Divide(2,3);
  correctionSliceCanvasX.cd(1);
  x1DCorrectionPlotsYMinus10XNegative->SetMinimum(-1*maxPlotResidual);
  x1DCorrectionPlotsYMinus10XNegative->SetMaximum(maxPlotResidual);
  x1DCorrectionPlotsYMinus10XNegative->Draw();
  x1DCorrectionPlotsYMinus10XNegative->Write();
  correctionSliceCanvasX.cd(2);
  x1DCorrectionPlotsYMinus10XPositive->SetMinimum(-1*maxPlotResidual);
  x1DCorrectionPlotsYMinus10XPositive->SetMaximum(maxPlotResidual);
  x1DCorrectionPlotsYMinus10XPositive->Draw();
  x1DCorrectionPlotsYMinus10XPositive->Write();
  correctionSliceCanvasX.cd(3);
  x1DCorrectionPlotsYCentral10XNegative->SetMinimum(-1*maxPlotResidual);
  x1DCorrectionPlotsYCentral10XNegative->SetMaximum(maxPlotResidual);
  x1DCorrectionPlotsYCentral10XNegative->Draw();
  x1DCorrectionPlotsYCentral10XNegative->Write();
  correctionSliceCanvasX.cd(4);
  x1DCorrectionPlotsYCentral10XPositive->SetMinimum(-1*maxPlotResidual);
  x1DCorrectionPlotsYCentral10XPositive->SetMaximum(maxPlotResidual);
  x1DCorrectionPlotsYCentral10XPositive->Draw();
  x1DCorrectionPlotsYCentral10XPositive->Write();
  correctionSliceCanvasX.cd(5);
  x1DCorrectionPlotsYPlus10XNegative->SetMinimum(-1*maxPlotResidual);
  x1DCorrectionPlotsYPlus10XNegative->SetMaximum(maxPlotResidual);
  x1DCorrectionPlotsYPlus10XNegative->Draw();
  x1DCorrectionPlotsYPlus10XNegative->Write();
  correctionSliceCanvasX.cd(6);
  x1DCorrectionPlotsYPlus10XPositive->SetMinimum(-1*maxPlotResidual);
  x1DCorrectionPlotsYPlus10XPositive->SetMaximum(maxPlotResidual);
  x1DCorrectionPlotsYPlus10XPositive->Draw();
  x1DCorrectionPlotsYPlus10XPositive->Write();
  correctionSliceCanvasX.SaveAs(residualsPDFName);
  correctionSliceCanvasX.Write();
  
  TCanvas correctionSliceCanvasY;
  correctionSliceCanvasY.Divide(2,3);
  correctionSliceCanvasY.cd(1);
  y1DCorrectionPlotsYMinus10XNegative->SetMinimum(-1*maxPlotResidual);
  y1DCorrectionPlotsYMinus10XNegative->SetMaximum(maxPlotResidual);
  y1DCorrectionPlotsYMinus10XNegative->Draw();
  y1DCorrectionPlotsYMinus10XNegative->Write();
  correctionSliceCanvasY.cd(2);
  y1DCorrectionPlotsYMinus10XPositive->SetMinimum(-1*maxPlotResidual);
  y1DCorrectionPlotsYMinus10XPositive->SetMaximum(maxPlotResidual);
  y1DCorrectionPlotsYMinus10XPositive->Draw();
  y1DCorrectionPlotsYMinus10XPositive->Write();
  correctionSliceCanvasY.cd(3);
  y1DCorrectionPlotsYCentral10XNegative->SetMinimum(-1*maxPlotResidual);
  y1DCorrectionPlotsYCentral10XNegative->SetMaximum(maxPlotResidual);
  y1DCorrectionPlotsYCentral10XNegative->Draw();
  y1DCorrectionPlotsYCentral10XNegative->Write();
  correctionSliceCanvasY.cd(4);
  y1DCorrectionPlotsYCentral10XPositive->SetMinimum(-1*maxPlotResidual);
  y1DCorrectionPlotsYCentral10XPositive->SetMaximum(maxPlotResidual);
  y1DCorrectionPlotsYCentral10XPositive->Draw();
  y1DCorrectionPlotsYCentral10XPositive->Write();
  correctionSliceCanvasY.cd(5);
  y1DCorrectionPlotsYPlus10XNegative->SetMinimum(-1*maxPlotResidual);
  y1DCorrectionPlotsYPlus10XNegative->SetMaximum(maxPlotResidual);
  y1DCorrectionPlotsYPlus10XNegative->Draw();
  y1DCorrectionPlotsYPlus10XNegative->Write();
  correctionSliceCanvasY.cd(6);
  y1DCorrectionPlotsYPlus10XPositive->SetMinimum(-1*maxPlotResidual);
  y1DCorrectionPlotsYPlus10XPositive->SetMaximum(maxPlotResidual);
  y1DCorrectionPlotsYPlus10XPositive->Draw();
  y1DCorrectionPlotsYPlus10XPositive->Write();
  correctionSliceCanvasY.SaveAs(residualsPDFName);
  correctionSliceCanvasY.Write();

  
  //Close PDF.
  dummy.SaveAs(residualsCloseString.Data());
  
  outputFile->cd();
  fResultTree->Write();
  outputFile->Close();
  
  return exitCode; 
}

void ParseConfigFile(const std::string& configFile) {
  //Open file.  
  ifstream file(configFile);
  //Parse lines in file.
  std::string line;
  cout << "[INFO] Parsing config file:" << endl;
  while (std::getline(file, line))
    {
      std::istringstream lineString(line);
      //Ignore lines beginning with a "#".
      if (lineString.str().front() == '#')
        continue;
      //Parse everything else.
      string variableName = "";
      int foundPosition = 0;
      string tpcName;
      int minClusters;
      if (lineString.str().find("minPotentialPointsRatio",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinPotentialPointsRatio)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Minimum potential points ratio for Vertex Tracks: "
	     << fMinPotentialPointsRatio << endl;
      }
      else if (lineString.str().find("maxPotentialPointsRatio",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxPotentialPointsRatio)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Maximum potential points ratio for Vertex Tracks: "
	     << fMaxPotentialPointsRatio << endl;
      }
      else if (lineString.str().find("maxBX",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxBX)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Maximum impact parameter X for track consideration (Vertex Tracks): "
	     << fMaxBX << endl;
      }
      else if (lineString.str().find("maxBY",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxBY)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Maximum impact parameter Y for track consideration (Vertex Tracks): "
	     << fMaxBY << endl;
      }
      else if (lineString.str().find("minEntries",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMinEntries)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Minimum bin entries for calculation: "
	     << fMinEntries << endl;
      }
      else if (lineString.str().find("maxAverageResidual",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxAverageResidual)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Maximum average residual for track consideration: "
	     << fMaxAverageResidual << endl;
      }
      else if (lineString.str().find("maxChi2",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> fMaxChi2)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
        cout << "[INFO] Maximum fit chi2 for track consideration (Vertex Tracks): "
	     << fMaxChi2 << endl;
      }
      else if (lineString.str().find("minClustersByTPC",foundPosition) != string::npos) {
        if (!(lineString >> variableName >> tpcName >> minClusters)) {
          cout << "[ERROR] File parsing failed! Line: " << lineString.str() << endl;
        }
	const det::TPCConst::EId tpcId = det::TPCConst::GetId(tpcName);
	fMinClustersByTPCId[tpcId] = minClusters;
        cout << "[INFO] Minimum " << tpcName << " clusters for track consideration (Local Tracks): "
	     << minClusters << endl;
      }
      
    } //End parsing.
  return;  
}
