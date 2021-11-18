#include "Functions.h"

using namespace std;
 
int main() {
    double S = 7420.40;  // Spot price
    double r = 0.001;   // Risk-free rate (0.1%)
    double T = 30.0/365.0;    // One year until expiry

    vector<double> optionsToPrice = {7450.00, 7500.00, 7550.00,
                                7600.00, 7700.00, 7750.00};

    vector<double> strikePrices = {7425.00, 7650.00,
                                    7750.00, 7800.00, 7825.00};

    vector<double> impVols = {0.123118, 0.122606, 0.109341,
                                0.134298, 0.108216};

    for (int i = 0; i < 6; i++){
    
    double strikeWanted = optionsToPrice[i];


    double volInt = volInterpolation(strikePrices, impVols, strikeWanted);

    tuple<double, double, double, double, double, double,double> callCalculations = monte_carlo_call_price(100000, S, strikeWanted, r, volInt, T, 0.1);

    cout << "Price: " << get<0>(callCalculations) << "\n";
    cout << "Vol:   " << volInt	<< "\n";
    cout << "Delta: " << get<1>(callCalculations) << "\n";
    cout << "Gamma: " << get<2>(callCalculations) << "\n";
    cout << "Vega:  " << get<3>(callCalculations) << "\n";
    cout << "Theta: " << get<4>(callCalculations) << "\n";
    cout << "Rho:   " << get<5>(callCalculations) << "\n";
    cout << "d:" << get<6>(callCalculations) << "\n";


    cout << "\n";

    }

    return 0;
}