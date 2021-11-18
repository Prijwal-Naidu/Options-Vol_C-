#include "Functions.h"

using namespace std;


int main() {
    double S = 7420.40;  // Spot price
    double r = 0.001;   // Risk-free rate (0.1%)
    double T = 30.0/365.0;    // One year until expiry
    double tol = 0.05;

    
    vector<double> strikePrices = {7425.00, 7650.00,
                                    7700.00, 7800.00, 7825.00};

    vector<double> marketPrices = {100.00, 25.00,
                                    15.00, 4.00, 9.00};
    
    vector<double> impVols;

    for (int i = 0; i < 6; i++ ) {
    
    double K = strikePrices[i];
    double marketPrice = marketPrices[i];
    double vol_estimate = impliedVolatility(S, K, r, T, tol, marketPrice);


    impVols.push_back(vol_estimate);
    std::cout << "The Implied Volatilty for XJO Call with Strike " << K << " is " << impVols[i] << endl; 

        }

    return 0;
}