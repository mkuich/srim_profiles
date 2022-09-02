#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TGraph.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TF1.h"

#include "./srim_profiles.h"

using namespace std;

// example of usage
int srim_profiles(){
    
// read SRIM file for alpha particle
    vector<vector<double>> vec_alpha = srim_file_read(file_alpha);
    
// make interpolation of alpha energy and range between the simulated points
    TF1 *alpha_interpolation = interpolate_range(vec_alpha, "Alpha");
    
// example of getting energy and range    
    cout<<"Range if alpha energy is 5003 keV "<<range_from_energy(alpha_interpolation, 5003)<<endl;
    cout<<"Energy if alpha ragne is 82 mm "<<energy_from_range(alpha_interpolation, 82)<<endl;

// example of dE/dx profile preparation, if we know only alpha length
    TH1D alpha_profile = energy_loss_profile(alpha_interpolation, energy_from_range(alpha_interpolation, 82), 25);
    
    return 0;
}
