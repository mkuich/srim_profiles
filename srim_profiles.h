#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"

using namespace std;

// simple input files from SRIM. Energy vs range
string file_alpha = "/mnt/27817d28-ee1d-4124-ab41-54755c9ee50d/work/EliTPC/Higs_analysis/vdg_data/srim/alpha_in_carbon01+oxygen02d0.422.srim.dat";
string file_carbon = "/mnt/27817d28-ee1d-4124-ab41-54755c9ee50d/work/EliTPC/Higs_analysis/vdg_data/srim/carbon_in_carbon01+oxygen02d0.422.srim.dat";
string file_proton = "/mnt/27817d28-ee1d-4124-ab41-54755c9ee50d/work/EliTPC/Higs_analysis/vdg_data/srim/proton_in_carbon01+oxygen02d0.422.srim.dat";
string file_nitrogen = "/mnt/27817d28-ee1d-4124-ab41-54755c9ee50d/work/EliTPC/Higs_analysis/vdg_data/srim/nitrogen_in_carbon01+oxygen02d0.422.srim.dat";

// scales the particle range to mm
double check_range_unit_conversion(double &range, string &unit){
    double factor;
    if(unit.compare("um")==0)
            factor=0.001;
    if(unit.compare("mm")==0)
            factor=1;
    if(unit.compare("cm")==0)
            factor=10;
    if(unit.compare("m")==0)
            factor=1000;
        
    return range*factor;
}

// scales the particle energy to keV
double check_energy_unit_conversion(double &energy, string &unit){
    double factor;
    if(unit.compare("eV")==0)
            factor=0.001;
    if(unit.compare("keV")==0)
            factor=1;
    if(unit.compare("MeV")==0)
            factor=1000;

    return energy*factor;
}

// reads SRIM files and returns vector of vectors result.at(0) is for energy, result.at(1) is for range
vector<vector<double>> srim_file_read(string &file_name){
    vector<vector<double>> result;
    vector<double> vec_energy, vec_range;
    ifstream infile(file_name.c_str());
    if(!infile.is_open()){
        cerr<<"Error! SRIM file not opened!"<<endl;
    }
    else
        cout<<"Reading SRIM file"<<endl;
    
    string Trash;
    for(int i=0; i<=27; i++){
        getline(infile, Trash);
    }
    while(!infile.eof()){
        double energy, energy_keV, range, range_mm;
        string unit_energy, unit_range;
        infile>>energy>>unit_energy>>Trash>>Trash>>range>>unit_range>>Trash>>Trash>>Trash>>Trash;
        if(energy==0.0){
            cout<<"File read successfully"<<endl;
            break;
        }
        range_mm=check_range_unit_conversion(range, unit_range);
        energy_keV=check_energy_unit_conversion(energy, unit_energy);
        vec_energy.push_back(energy_keV);
        vec_range.push_back(range_mm);
    }
    result.push_back(vec_energy);
    result.push_back(vec_range);
    return result;
}

// fits range(energy) from SRIM to later interpolate between the simulated points
TF1 *interpolate_range(vector<vector<double>> &result, string title){
    TCanvas tmp("tmp", "", 800, 600);
    TH2F frame("frame", "", 1000, 0, result.at(0).back(), 1000, 0, result.at(1).back());
    frame.GetXaxis()->SetTitle("Energy [keV]");
    frame.GetYaxis()->SetTitle("range [mm]");
    frame.SetTitle(title.c_str());
    frame.SetStats(0);
    frame.Draw("");
    TGraph graph(result.at(0).size());
    graph.SetName(title.c_str());
    for(int i=0; i<result.at(0).size(); i++){
        graph.SetPoint(i, result.at(0).at(i), result.at(1).at(i));
    }
    TF1 *range_interpolation = new TF1(TString::Format("%s_range_interpolation", title.c_str()), "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x+[6]*x*x*x*x*x*x+[7]*x*x*x*x*x*x*x+[8]*x*x*x*x*x*x*x*x", 0, result.at(0).back());
    graph.SetMarkerStyle(kFullSquare);
    graph.Draw("same p");
    graph.Fit(TString::Format("%s_range_interpolation", title.c_str()),"emir");
//     tmp.Print(TString::Format("%s_range.png", title.c_str()));
    
    return range_interpolation;
}

// gets range of particle of a given energy
double range_from_energy(TF1 *interpolation, double energy){
    return interpolation->Eval(energy);
}

// gets energy of particle of a given range 
double energy_from_range(TF1 *interpolation, double range){
    return interpolation->GetX(range);
}

// creates dE/dx profile with bin size of default_step or slightly smaller (if the rage/default_step doesn't give an intiger)
TH1D energy_loss_profile(TF1 *interpolation, double energy, double default_step=1.0){
    double tmp_steps = range_from_energy(interpolation, energy)/default_step;
    int n_steps=1;
    if(!((tmp_steps-(int)tmp_steps)==0))
        n_steps = (int)tmp_steps+1;
    else
        n_steps=(int)tmp_steps;

    TH1D hprofile(TString::Format("hprof_energy_%d", (int)energy), "", n_steps, 0, range_from_energy(interpolation, energy));
    hprofile.Sumw2();
    int j=1;
    for(int i=hprofile.GetNbinsX(); i>0; i--){
        if(hprofile.GetBinLowEdge(i)==0){
            hprofile.SetBinContent(j, (energy_from_range(interpolation, hprofile.GetBinLowEdge(i)+hprofile.GetBinWidth(i))));
            break;
        }
        
        hprofile.SetBinContent(j, (energy_from_range(interpolation, hprofile.GetBinLowEdge(i)+hprofile.GetBinWidth(i))-energy_from_range(interpolation, hprofile.GetBinLowEdge(i))));
        j++;
    }
    
    return hprofile;
}

// reverses profile (optional)
TH1D reverse_profile(TH1D original_profile){
     TH1D hprofile = original_profile;
     hprofile.Reset();
     int j=1;
     for(int i=original_profile.GetNbinsX(); i>0; i--){
         hprofile.SetBinContent(j,original_profile.GetBinContent(i));
         j++;
    }
    hprofile.SetName(TString::Format("%s_reversed", original_profile.GetName()));
    return hprofile;
}

// normalizes profile (optional)
TH1D normalise_profile(TH1D original_profile, double norm_value=1.0){
    TH1D hprofile = original_profile;
    hprofile.Scale(norm_value/hprofile.GetMaximum());
    hprofile.SetName(TString::Format("%s_normalised", original_profile.GetName()));
    
    return hprofile;
}
