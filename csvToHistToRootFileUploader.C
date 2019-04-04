// file that converts csv files from URA or ellipsometery that contain data in two column form
// wavelength (nm), %reflection (out of 100)
// use python file csvConverter.py to generate these if do not already exist
// Author: Sol Cotton 26/03/19

#include <iostream>
#include <vector>
#include <cstring>
#include <string>
#include <sstream>
#include <TObject.h>


string toUpper(string input){
	// uppercases all chars in input string
	for (auto& x : input)
		x = toupper(x);
	return input;
}

bool getStringAOrB(string a, string b){
	// returns true on validated "a" input, false on "b" input
	// uppercase everything for consistency
	string aOrB;
	while (aOrB.compare(toUpper(a)) && aOrB.compare(toUpper(b))){
		cout << "Please enter " << a << " or " << b << ": " << endl;
		cin >> aOrB;
		cin.ignore(100, '\n');
		aOrB = toUpper(aOrB);
	}
	if (aOrB == toUpper(a)) return true;
	else if (aOrB == toUpper(b)) return false;
	else exit(0);
}

void convert(bool verbosity){
	
	
	// Set these values before running! DON'T FORGET CSV
	string inputCSVfile = "RaleighDist_2p0L.txt";
	string outputHistName = "RaylighDistribution_200%L";

	/*
	// FOR NORMAL HISTS
	// create a 2d vector that will store float data from each file
	// Values[linenumber j in file][0: xvalue, 1: yvalue]
	vector< vector< float > > dataVector;
	//loop over all values
    ifstream in(inputCSVfile.c_str());
    if(!in) cout << "File not found" << endl;
    string line, sector;
    while(getline(in, line)){ 
    	//loop over all lines, make a 1d vector for line float data
        vector<float> v;
        stringstream ss(line);
         while (getline(ss, sector, ',')) { 
         	// loop over all sectors in line, deliminated by commas
         	stringstream converter(sector);
         	float floatValue = 0.000000;
         	converter >> floatValue;
            v.push_back(floatValue);  
            // add each sector as a float to the 1D vector
            if(!ss.good() ) break;
        }
        // add the 1D vector to the 2D vector
        dataVector.push_back(v); 
    }
	// Generate a title for the histograms 
	ostringstream histogramNameStream;
	histogramNameStream << outputHistName;
	
	// Create a new histogram to add the read in data to
	TH1F *h_added = new TH1F(histogramNameStream.str().c_str(), "", 800, 0, 800);

	cout << "nm HISTOGRAM UPLOADER:\nConverting CSV " << setw(20) << inputCSVfile << "\n to TH1F object " << setw(20)
		 << histogramNameStream.str().c_str() << endl;

	cout << "Does this look right? \nIf not, change the infile and filter information in the script.\n";
	if (!getStringAOrB("y", "n")) exit(0);

	for (int q = 1; q<801; q++){
		// loop over all bins in new histogram
		double localsum = 0, localnum = 0;
		for (int j = 0; j < dataVector.size(); j++){
			// loop over all values in read in array
			if (dataVector[j][0] >= q-1.5 && dataVector[j][0] < q+1.5){
				// q runs through the wavelength (bin) values in h_added
				// sum any yValues whose xValues lie within a region 3nm around q 
				localsum += dataVector[j][1];
				localnum += 1;
			}
        }
        // set the value of bin (wavelength) q in h_added to the average of the yValues in the region
        if(localnum != 0) h_added->SetBinContent(q, localsum/localnum);
		if(localnum == 0) h_added->SetBinContent(q, 0);
		
		if (verbosity) cout << "Bin: " << q << setw(5) << " Value: " << h_added->GetBinContent(q) << endl;
	}
	*/

	// FOR RAYLEIGH DISTRIBUTIONS
	///*
	TH1F *h_added = new TH1F("h_added", "", 1000, 0, 3.141529);

	// generate a title for the histograms
	ostringstream histogramNameStream;
	histogramNameStream << outputHistName;

	cout << "DISTRIBUTION UPLOADER:\nConverting TXT " << setw(20) << inputCSVfile << "\n to TH1F object " << setw(20)
		 << histogramNameStream.str().c_str() << endl;

	cout << "Does this look right? \nIf not, change the infile and filter information in the script.\n";
	if (!getStringAOrB("y", "n")) exit(0);

	ifstream angfile;
	angfile.open(inputCSVfile);
	string line;
	int counter = 0;
	angfile.seekg(0);
	while (getline(angfile, line)){
		istringstream iss(line);
		float a, b;
		if (!(iss >> a >> b)) { cout << "bleargh"; break; } // error
		h_added->SetBinContent(counter, b);
		if(verbosity) cout << "Bin: " << counter << setw(3) << " Value: " << h_added->GetBinContent(counter);
		counter++;
	}

	//*/


	// open the root file
	TFile *MyFile = new TFile("wlsModellerHists.root","UPDATE");

	// uncomment for excitation spectra
	/*h_added->SetLineWidth(0);
	h_added->SetFillStyle(3002);
	h_added->SetFillColor(kRed-7);
	h_added->SetLineColor(kRed-7);*/

	// uncomment for emission spectra
	//h_added->SetLineWidth(0);
	//h_added->SetFillStyle(3003);
	//h_added->SetFillColor(kBlue-7);
	//h_added->SetLineColor(kBlue-7);
	
	// add the histogram object, overwrite if already exists
	h_added->Write(histogramNameStream.str().c_str(), TObject::kOverwrite);
	cout << "\nSuccessfully added";
}