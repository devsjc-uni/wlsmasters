// Program that plots and superposes an arbitrary number of excel files or root files as graphs
// has the option to include TPB and PEN emission and absorption spectra
// attmepts to model rayleigh scattering as well
// wlsmodeller 4 attempts to take into account the fact that the sample is at an angle
// AND takes a histogram as an input for the distribution as opposed to a TF1
// AND reads in all hist objects from a root file
// Author: Sol Cotton

#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <math.h>
#include <vector> 
#include <string>
#include "TFile.h"
#include <sstream>
#include "TChain.h"
#include "TRandom.h"
#include "TSystemDirectory.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TString.h"
#include "TAxis.h"
#include "TPaveLabel.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TGraphErrors.h"
#define _USE_MATH_DEFINES

using namespace std;

float midpoint(int a, int b){
	return (a + b)/2.;
}

float toRadians(float a){
	return M_PI * a / 180;
}

void model(bool verbosity, string wavelengthShifterType, string filterManufacterer, int filterCutoff, int mfpToLengthPercent){

	bool spectra = 0;

	// open root file containing hists
	TFile *file = TFile::Open("wlsModellerHists.root", "READ");
	file->Print();

	// ---------------------------------- GET CHOSEN FILTER
	
	// choose the filter we want by placing the relevant angular data according to which filter we want to model
	if (verbosity) cout << "\n\nSelecting hists of angular filter data for " << filterCutoff << " " << filterManufacterer << "filter...\n";
	// define arrays containing the angles of incidence the filters were measured at
	vector<int> filterAngles_Edmund = {8, 25, 45, 50, 55, 60, 65, 70, 75};
	vector<int> filterAngles_THORLABS = {8, 25, 45, 55, 65, 75};
	vector<int> filterAngles_Asahi = {8, 25, 45, 50, 55, 60, 65, 70, 75};
	// select the relevant filter
	vector<int> filterAngles;
	if (filterManufacterer == "Edmund") for (int i = 0; i < filterAngles_Edmund.size(); i++) { filterAngles.push_back(filterAngles_Edmund[i]); }
	if (filterManufacterer == "THORLABS") for (int i = 0; i < filterAngles_THORLABS.size(); i++) { filterAngles.push_back(filterAngles_THORLABS[i]); }
	if (filterManufacterer == "Asahi") for (int i = 0; i < filterAngles_Asahi.size(); i++) { filterAngles.push_back(filterAngles_Asahi[i]); }


	// create an array of histograms to store the chosen filter's angular data
	TH1F *h_AngularFilterData[filterAngles.size()+1];
	for (int i = 0; i < filterAngles.size(); i++){
		// loop over all the different measured angles for the chosen filter
		// create a stringstream object for the chosen filter's specific angles
		ostringstream AngleDataHistogramNameStream;
		if (filterManufacterer == "Asahi") AngleDataHistogramNameStream << filterManufacterer << "_" << filterCutoff << "_" << filterAngles[i] << "_rp";
		else AngleDataHistogramNameStream << filterManufacterer << "_" << filterCutoff << "_" << filterAngles[i];
		// import the relevant histogram from the TFile
		h_AngularFilterData[i] = (TH1F*)file->Get(AngleDataHistogramNameStream.str().c_str());
		if (verbosity) cout << "Imported " << AngleDataHistogramNameStream.str().c_str() << endl;
	}

	// create a final angular filter data histogram for light deflected beyond 90 which just sees the sphere
	h_AngularFilterData[filterAngles.size()] = new TH1F("backscatter", "", 800, 0, 800);
	for (int j = 0; j < h_AngularFilterData[filterAngles.size()]->GetNbinsX(); j++){
		// set all values to 100, as if light just hits sphere
		h_AngularFilterData[filterAngles.size()]->SetBinContent(j, 100); 
	}

	// ---------------------------------- CARRY OUT RAYLEIGH SCATTERING
	
	if (verbosity) cout << "\n\nCarrying out rayleigh scattering..." << endl;
	// read in the effect of TPB on an ideal reflector
	TH1F *h_TPBEffectOnIdealReflector = (TH1F*)file->Get("TPB_EffectOnIdeal_R");
	if (verbosity) cout << "Imported TPB_EffectOnIdeal_R" << endl;
	// read in angular distribution generated from MC - number of bins must match sim hist and TH1F object when uploading to root file
	ostringstream AngularDistributionHistogramNameStream;
	AngularDistributionHistogramNameStream << "RaylighDistribution_" << mfpToLengthPercent << "%L";
	TH1F *h_RayleighAngularDistribution = (TH1F*)file->Get(AngularDistributionHistogramNameStream.str().c_str()); 
	if (verbosity) cout << "Imported " << AngularDistributionHistogramNameStream.str().c_str() << endl;

	// Calculate, from the angular distribution, the fraction of photons hitting the filter at each available measured angle
	if (verbosity) cout << "\nCaluclating midpoints of filter angles..." << endl;
	// create a vector for the midpoints
	vector<float> angleMidpoints;
	// push back the first boundary as zero 
	angleMidpoints.push_back(0);
	for (int i = 0; i < filterAngles.size() - 1; i++){
		// push back the next boundaries as the midpoints between the filter angles
		angleMidpoints.push_back(toRadians(midpoint(filterAngles[i], filterAngles[i+1])));
	}
	// the penultimate boundary is 90 degrees, and the final boundary is 180 degrees
	angleMidpoints.push_back(toRadians(90));
	angleMidpoints.push_back(toRadians(180));
	if (verbosity) {for (int i = 0; i < angleMidpoints.size(); i++) cout << angleMidpoints[i] << endl;}
	// define a set of floats to contain the fraction of photons binned to each filter angle
	float fractionAtAngle[filterAngles.size() + 1];
	// since rayleigh histogram has 1000 bins spanning from 0 to pi, every bin is pi/1000, so multiply midpoints by 1000/pi to get bins
	// get total integral of rayleigh distribution
	float rayleighTotalIntegral = h_RayleighAngularDistribution->Integral(0, 1000);
	if (verbosity) cout << "\nCalculating max an min angles to be binned to each filter angle..." << endl;
	for (int i = 0; i < angleMidpoints.size() - 1; i++){
		// loop over all midpoints, get fraction of distribution between each set of boundaries
		int lowerBin = int (angleMidpoints[i]*1000/M_PI);
		int upperBin = int (angleMidpoints[i+1]*1000/M_PI); 
		fractionAtAngle[i] = h_RayleighAngularDistribution->Integral(lowerBin, upperBin) / rayleighTotalIntegral;
		if (verbosity) cout << "Lower bin: " << lowerBin << setw(3) << ", Upper bin: " << upperBin << endl;
	}
	// normalise the fractions to 1
	float normalfracTot = 0;
	int numberOfAngles = sizeof(fractionAtAngle)/sizeof(fractionAtAngle[0]);
	for (int i = 0; i < numberOfAngles; i++){ normalfracTot += fractionAtAngle[i]; }
	float normalScale = 1/normalfracTot; 
	for (int i = 0; i < numberOfAngles; i++){ fractionAtAngle[i] = fractionAtAngle[i]*normalScale; }
	if (verbosity) {
		cout << "\nRayleigh distribution fractions: " << endl;
		for (int i = 0; i < numberOfAngles - 1; i++){
			cout << filterAngles[i] << "deg: " << setw(3) << fractionAtAngle[i] << endl;
		}
		cout << "backscatter: " << setw(3) << fractionAtAngle[numberOfAngles - 1] << endl;
	}
	// at each wavelength, calculate the average R/T that is scattered to over all filter angles
	TH1F *h_RayleighModel = new TH1F("h_RayleighModel", "", 800, 0, 800);
	for (int j = 0; j < 800; j++){
		// loop over all wavelengths
		float currentavg = 0;
		for (int i = 0; i < numberOfAngles; i++){
			// at each wavelength, loop over all the angles, and add to the running average their value of T/R, scaled by their involvement fraction from the MC dist
			currentavg += h_AngularFilterData[i]->GetBinContent(j)*(fractionAtAngle[i]);
		}
		// set the value of the model to be the average multiplied by the modifier that is the effect on an ideal reflector/transmitter
		h_RayleighModel->SetBinContent(j, currentavg*h_TPBEffectOnIdealReflector->GetBinContent(j)*0.01);
		//h_RayleighModel->SetBinContent(j, currentavg);
	}

	// ---------------------------------- ADD THE WLS EFFECT

	if (verbosity) cout << "\n\nAdding WLS effect..." << endl;
	// import excitation and emmission
	ostringstream ExcitationHistogramNameStream, EmissionHistogramNameStream;
	ExcitationHistogramNameStream << wavelengthShifterType << "_Excitation";
	EmissionHistogramNameStream << wavelengthShifterType << "_Emission";
	TH1F *h_EmissionSpectrum = (TH1F*)file->Get(EmissionHistogramNameStream.str().c_str());
	TH1F *h_ExcitationSpectrum = (TH1F*)file->Get(ExcitationHistogramNameStream.str().c_str());
	if (verbosity) cout << "Imported " << EmissionHistogramNameStream.str().c_str() << "\nImported " << ExcitationHistogramNameStream.str().c_str() << endl;
	// calculate the average transmission/reflection percentage of transmitted light
	// by carrying out sum over i of [filter_T/R_spectrum_bin[i] * (emission_spectrum_bin[i] / 100)] / integral of (em_spect_bin[i]/100)
	// and adding this on to the attempted recreated graph*the absorption spectra value at bin [i] 
	if (verbosity) cout << "\n\nCalculating average T/R from emission spec for WLS effect -------------------\n\n" ;
	float emissionAverage[numberOfAngles];
	// calculate average shifted-to reflection/transmission for each angle
	for (int i = 0; i < numberOfAngles; i++){
		// loop over the different angles and calculate the average shifted-to wavelength for each
		emissionAverage[i] = 0;
		for (int j = 0; j < 800; j++){
			// add to the average the product of the transmission/reflection value and the emission spectrum as a probability
			emissionAverage[i] += h_AngularFilterData[i]->GetBinContent(j)*h_EmissionSpectrum->GetBinContent(j)/100.;
		}
		// normalise to the total integral of the emission spectrum
		emissionAverage[i] = emissionAverage[i] / (h_EmissionSpectrum->Integral(0,800)/100.);
		if (verbosity) cout << " \nEmission Avg " << i << ": " << emissionAverage[i] << endl;
	}
	// emission of shifted light is isotropic, so define isotropic distribution
	TH1F *h_IsotropicAngularDistribution = new TH1F("h_IsotropicAngularDistribution", "", 1000, 0, 3.141529);
	// initialise it as 1 everywhere (isotropic)
	for (int i = 0; i < h_IsotropicAngularDistribution->GetNbinsX(); i++) h_IsotropicAngularDistribution->SetBinContent(i, 1);
	// calculate the fractions at each angle in the same manner as for rayleigh
	float isotropicTotalIntegral = h_IsotropicAngularDistribution->Integral(0, 1000);
	if (verbosity) cout << "\nCalculating max and min angles to be binned to each filter angle..." << endl;
	for (int i = 0; i < angleMidpoints.size() - 1; i++){
		// loop over all midpoints, get fraction of distribution between each set of boundaries
		int lowerBin = int (angleMidpoints[i]*1000/M_PI);
		int upperBin = int (angleMidpoints[i+1]*1000/M_PI); 
		fractionAtAngle[i] = h_IsotropicAngularDistribution->Integral(lowerBin, upperBin) / isotropicTotalIntegral;
		if (verbosity) cout << "Lower bin: " << lowerBin << setw(3) << ", Upper bin: " << upperBin << endl;
	}
	// normalise the fractions to 1
	normalfracTot = 0;
	for (int i = 0; i < numberOfAngles; i++){ normalfracTot += fractionAtAngle[i]; }
	normalScale = 1/normalfracTot; 
	for (int i = 0; i < numberOfAngles; i++){ fractionAtAngle[i] = fractionAtAngle[i]*normalScale; }
	if (verbosity) {
		cout << "\nIsotropic distribution fractions: " << endl;
		for (int i = 0; i < numberOfAngles - 1; i++){
			cout << filterAngles[i] << "deg: " << setw(3) << fractionAtAngle[i] << endl;
		}
		cout << "backscatter: " << setw(3) << fractionAtAngle[numberOfAngles - 1] << endl;
	}
	// get average shifted-to reflection over all angles
	float overallShiftedToR = 0;
	for (int i = 0; i < numberOfAngles; i++){
		overallShiftedToR += emissionAverage[i]*(fractionAtAngle[i]);
	}
	if (verbosity) cout << "\nCalculating effect from wavelength shifting..." << endl;
	TH1F *h_RayleighAndWLSModel = (TH1F*)h_RayleighModel->Clone("h_RayleighAndWLSModel");
	for (int i = 0; i < 800; i++){
		// define the effect of a wavelength shift as the overall shifted-to reflection value, scaled by the value of the excitation spectrum at each wavelength 
		float effectOfShift = (h_ExcitationSpectrum->GetBinContent(i)/100.)*overallShiftedToR;
		// and the effect of a not shift as 1 minus this probability, multiplied by the already existing model value
		float effectOfNotShift = (1-(h_ExcitationSpectrum->GetBinContent(i)/100.))*h_RayleighModel->GetBinContent(i);
		// set the model values as the sum of these effects
		h_RayleighAndWLSModel->SetBinContent(i, effectOfShift + effectOfNotShift);
	}

	// ---------------------------------- CORRECT FOR THE SPHERE TO URA DIFFERENCE

	if (verbosity) cout << "\n\nCorrecting for sphere difference (to URA)..." << endl;
	TH1F *h_RayleighAndWLSModel_SphereCorrected = (TH1F*)h_RayleighAndWLSModel->Clone("h_RayleighAndWLSModel_SphereCorrected");
	// import the relevant correction from the root file
	ostringstream SphereDiffHistogramNameStream;
	SphereDiffHistogramNameStream << filterManufacterer << "_" << filterCutoff << "_SphereCorrection";
	TH1F *h_SphereCorrection = (TH1F*)file->Get(SphereDiffHistogramNameStream.str().c_str());
	if (verbosity) cout << "Imported " << SphereDiffHistogramNameStream.str().c_str() << endl;
	// at all wavelengths, add on the sphere correction
	for (int i = 0; i < 800; i++){
		h_RayleighAndWLSModel_SphereCorrected->SetBinContent(i, h_RayleighAndWLSModel->GetBinContent(i) + h_SphereCorrection->GetBinContent(i));
	}

	// --------------------------------- PLOT ALL THE GRAPHS 

	// define style for plots
	gStyle->SetOptTitle(kFALSE);
	gStyle->SetOptStat(kFALSE);
    gStyle->SetPalette(91); 
    gStyle->SetHatchesSpacing(1.3);
    gStyle->SetErrorX(0);
    // create a new canvas for histograms to plot to
	TCanvas *c1 = new TCanvas("c1", "Title", 900, 700);
	c1->SetLeftMargin(0.15);
	c1->SetBottomMargin(0.15);
	
	// import measured spectra of TPB on filter to compare with model
	ostringstream MeasuredHistogramNameStream;
	MeasuredHistogramNameStream << wavelengthShifterType << "on_" << filterManufacterer << "_" << filterCutoff;
	TH1F *h_MeasuredWLSOnFilterSpectrum = (TH1F*)file->Get(MeasuredHistogramNameStream.str().c_str());
	// Set plotting parameters
	h_MeasuredWLSOnFilterSpectrum->GetXaxis()->SetRangeUser(380, 800);
	h_MeasuredWLSOnFilterSpectrum->GetYaxis()->SetRangeUser(0, 100);
	h_MeasuredWLSOnFilterSpectrum->GetXaxis()->SetLabelSize(0.06);
	h_MeasuredWLSOnFilterSpectrum->GetXaxis()->SetLabelOffset(0.006);
	h_MeasuredWLSOnFilterSpectrum->GetXaxis()->SetTitleSize(0.06);
	h_MeasuredWLSOnFilterSpectrum->GetYaxis()->SetLabelSize(0.06);
	h_MeasuredWLSOnFilterSpectrum->GetYaxis()->SetLabelOffset(0.006);
	h_MeasuredWLSOnFilterSpectrum->GetYaxis()->SetTitleSize(0.06);
	h_MeasuredWLSOnFilterSpectrum->GetXaxis()->SetTitleOffset(1);
	h_MeasuredWLSOnFilterSpectrum->GetYaxis()->SetTitleOffset(1);
	h_MeasuredWLSOnFilterSpectrum->GetXaxis()->SetTitle("Wavelength [nm]");
	h_MeasuredWLSOnFilterSpectrum->GetYaxis()->SetTitle("Reflection [%]");
	h_MeasuredWLSOnFilterSpectrum->SetLineWidth(2);

	// draw histograms
	h_MeasuredWLSOnFilterSpectrum->Draw("PLC L");
	h_RayleighAndWLSModel_SphereCorrected->SetLineColor(kBlack);
	h_RayleighAndWLSModel_SphereCorrected->SetLineWidth(2);
	h_RayleighAndWLSModel_SphereCorrected->Draw("L SAME");
	//h_RayleighAndWLSModel->Draw("PLC L SAME");
	if (spectra) h_ExcitationSpectrum->Draw("SAME");
	if (spectra) h_EmissionSpectrum->Draw("SAME");
	
	// draw one legend
	TLegend *legend = new TLegend(0.5,0.6,0.8,0.8);
	legend->SetFillColor(4000);
   	legend->SetLineColor(0);
	// if desired, plot emission and absorption spectra
	if (spectra) legend->AddEntry(h_EmissionSpectrum, "Emission spectrum", "f");
	if (spectra) legend->AddEntry(h_ExcitationSpectrum, "Excitation spectrum", "f");
	legend->AddEntry(h_RayleighAndWLSModel_SphereCorrected, "Model", "l");
	legend->AddEntry(h_MeasuredWLSOnFilterSpectrum, MeasuredHistogramNameStream.str().c_str(), "l");
   	legend->Draw("SAME");

   	// DRAW TITLE
   	string titlestring = "wlsmodeller4: " + wavelengthShifterType + " " + filterManufacterer + " " 
   						 + to_string(filterCutoff) + "nm, WLS, sphere corrected, mfp=" + to_string(mfpToLengthPercent) + "%L";    
   	TPaveLabel *title = new TPaveLabel(.01,.93,.99,1.0, titlestring.c_str(), "brNDC");
	title->SetBorderSize(0);
	title->SetFillColor(0);
	title->Draw("SAME");

}