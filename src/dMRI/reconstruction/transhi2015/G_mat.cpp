#include "recon_transhi2015.h"

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> makeGmat(double D_inAx, double D_trapped, int b)
{

    // maxl=8, D_inAx and D_trapped are eigenvalues at a voxel
    double b1 = b * (D_inAx - D_trapped);
    double b2 = erf(sqrt(b1));
    double pi = M_PI;
    
    // Temp int e^{-b1*t^2}t^{l}, l= 0,2,4,6,8,10,12,14,16,18
    // using Wolframalpha to calculate a0,a2,a4,a6,a8 and G(:)
    double a0 = sqrt(pi/b1)*b2;
    double a2 = sqrt(pi)*b2/(2*pow(b1,(1.5)))-1.0/(b1*exp(b1));
    double a4 = -(3+2*b1)/(exp(b1)*2*pow(b1,2)) + (3*sqrt(pi)*b2)/(4*pow(b1,2.5));
    double a6 = 15*sqrt(pi)*b2/(8*pow(b1,3.5)) - (4*pow(b1,2)+10*b1+15)/(4*pow(b1,3)*exp(b1));
    double a8 = 105*sqrt(pi)*b2/(16*pow(b1,4.5)) - (2*b1*(2*b1*(2*b1+7)+35)+105)/(8*pow(b1,4)*exp(b1));
    double a10 = 945*sqrt(pi)*b2/(32*pow(b1,5.5)) - (2*b1*(2*b1*(2*b1*(2*b1+9)+63)+315)+945)/(16*pow(b1,5)*exp(b1));
    double a12 = 10395*sqrt(pi)*b2/(64*pow(b1,6.5)) - (2*b1*(2*b1*(2*b1*(4*pow(b1,2)+22*b1+99)+693)+3465)+10395)/(32*pow(b1,6)*exp(b1));
    double a14 = 135135*sqrt(pi)*b2/(128*pow(b1,7.5)) - (64*pow(b1,6) + 416*pow(b1,5) + 2288*pow(b1,4) + 10296*pow(b1,3)+36036*pow(b1,2) + 90090*b1+135135)/(64*pow(b1,7)*exp(b1));
    double a16 = 2027025*sqrt(pi)*b2/(256*pow(b1,8.5)) - (128*pow(b1,7) + 960*pow(b1,6) + 6240*pow(b1,5) + 34320*pow(b1,4) + 154440*pow(b1,3) + 540540*pow(b1,2) + 1351350*b1 + 2027025)/(128*pow(b1,8)*exp(b1));

    // cout << "a0: " << a0 << '\n';
    // cout << "a2: " << a2 << '\n';
    // cout << "a4: " << a4 << '\n';
    // cout << "a6: " << a6 << '\n';
    // cout << "a8: " << a8 << '\n';
    // cout << "a10: " << a10 << '\n';
    // cout << "a12: " << a12 << '\n';
    // cout << "a14: " << a14 << '\n';
    // cout << "a16: " << a16 << '\n';

    // G_l(D_inAx,D_trapped) := Ratio between the SH coefficients of s and those of the fod, l= 0,2,4,6,8
    std::vector<double> G(9, 0);
    G[0] = a0; // l = 0
    G[1] = 1.5*a2-0.5*a0; // l= 2
    G[2] = 35/8.0*a4-15/4.0*a2+3/8.0*a0; // l = 4
    G[3] = (231.0*a6 - 315*a4 + 105*a2 - 5*a0)/16.0; // l = 6
    G[4] = (6435.0*a8 - 12012*a6 + 6930*a4 - 1260*a2 + 35*a0)/128.0; // l = 8
    G[5] = (46189.0*a10 - 109395*a8 + 90090*a6 - 30030*a4 + 3465*a2 - 63*a0)/256.0; // l = 10
    G[6] = (676039.0*a12 - 1939938*a10 + 2078505*a8 - 1021020*a6 + 225225*a4 - 18018*a2 + 231*a0)/1024.0; // l = 12
    G[7] = (5014575.0*a14 - 16900975*a12 + 22309287*a10 - 14549535*a8 + 4849845*a6 - 765765*a4 + 45045*a2 - 429*a0)/2048.0; //  l = 14
    G[8] = (300540195.0*a16 - 1163381400.0*a14 + 1825305300.0*a12 - 1487285800.0*a10 + 669278610.0*a8 - 162954792.0*a6 + 19399380.0*a4 - 875160.0*a2 + 6435.0*a0)/32768.0; // l = 16
    
    // cout << "G0: " << G[0] << '\n';
    // cout << "G1: " << G[1] << '\n';
    // cout << "G2: " << G[2] << '\n';
    // cout << "G3: " << G[3] << '\n';
    // cout << "G4: " << G[4] << '\n';
    // cout << "G5: " << G[5] << '\n';
    // cout << "G6: " << G[6] << '\n';
    // cout << "G7: " << G[7] << '\n';
    // cout << "G8: " << G[8] << '\n';

    for(std::size_t i = 0; i < G.size(); i++) {
        G[i] = G[i]*2*pi*exp(-b*D_trapped);
        if(isnan(G[i])) NIBR::disp(MSG_ERROR,"G_matrix G_vector error");
    }

    // H1 = \dfrac{\partial G}{\partial \lambda_1}
    std::vector<double> H1(G.size(), 0);
    H1[0] = a2; // l = 0
    H1[1] = (3*a4-a2)/2; // l = 2
    H1[2] = (35*a6-30*a4+3*a2)/8; // l = 4
    H1[3] = (231*a8 - 315*a6 + 105*a4 - 5*a2)/16; //  l = 6
    H1[4] = (6435*a10 - 12012*a8 + 6930*a6 - 1260*a4 + 35*a2)/128; // l = 8
    H1[5] = (46189*a12 - 109395*a10 + 90090*a8 - 30030*a6 + 3465*a4 - 63*a2)/256.0; // l = 10
    H1[6] = (676039*a14 - 1939938*a12 + 2078505*a10 - 1021020*a8 + 225225*a6 - 18018*a4 + 231*a2)/1024; // l = 12

    for(std::size_t i = 0; i < H1.size(); i++) {
        H1[i] = -2*b*pi*exp(-b*D_trapped)*H1[i];
        if(isnan(H1[i])) NIBR::disp(MSG_ERROR,"G_matrix H1_vector error");
    }
    // H2 = \dfrac{\partial G}{\partial \lambda_2}
    std::vector<double> H2(G.size(), 0);
    for(std::size_t i = 0; i < H2.size(); i++) {
        H2[i] = -b*G[i]-H1[i];
        if(isnan(H2[i])) NIBR::disp(MSG_ERROR,"G_matrix H1_vector error");
    }

    return make_tuple(G, H1, H2);
}

std::vector<std::vector<double>> NIBR::G_mat(std::vector<float> bval, double D_inAx, double D_trapped, int maxOrder)
{
    std::vector<std::vector<double>> output(bval.size(), std::vector<double>(153, 0));

    // Selecting unique elements from the bvalues
    std::vector<float> uniqueBvalues = bval;
    sort(uniqueBvalues.begin(), uniqueBvalues.end());
    auto uniqueEnd = unique(uniqueBvalues.begin(), uniqueBvalues.end());
    uniqueBvalues.resize(distance(uniqueBvalues.begin(), uniqueEnd));
    
    for(auto i : uniqueBvalues) {
        std::vector<std::size_t> ind;
        for(std::size_t j = 0; j < bval.size(); j++) {
            if(bval[j] == i) ind.push_back(j);
        }

        disp(MSG_DEBUG,"%f %f %f", i, D_inAx, D_trapped);

        std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> G = makeGmat(D_inAx, D_trapped, i);

        std::vector<double> g = std::get<0>(G);
        g.resize(maxOrder / 2 + 1);
        
        for(auto index : ind) {
            output[index][0] = g[0];
            for(std::size_t g_ind = 1; g_ind < 6; g_ind++) output[index][g_ind] = g[1];
            for(std::size_t g_ind = 6; g_ind < 15; g_ind++) output[index][g_ind] = g[2];

            if(maxOrder > 4) {
                for(std::size_t g_ind = 15; g_ind < 28; g_ind++) output[index][g_ind] = g[3];
            }

            if(maxOrder > 6) {
                for(std::size_t g_ind = 28; g_ind < 45; g_ind++) output[index][g_ind] = g[4];
            }

            if(maxOrder > 8) {
                for(std::size_t g_ind = 45; g_ind < 66; g_ind++) output[index][g_ind] = g[5];
            }

            if(maxOrder > 10) {
                for(std::size_t g_ind = 66; g_ind < 91; g_ind++) output[index][g_ind] = g[6];
            }

            if(maxOrder > 12) {
                for(std::size_t g_ind = 91; g_ind < 120; g_ind++) output[index][g_ind] = g[7];
            }

            if(maxOrder > 14) {
                for(std::size_t g_ind = 120; g_ind < 153; g_ind++) output[index][g_ind] = g[8];
            }
        }
    }

    // Rmoving unnecessary padding
    for(std::size_t row = 0; row < output.size(); row++) {
        while(!output[row].empty() && output[row].back() == 0) output[row].pop_back();
    }

    return output;
}