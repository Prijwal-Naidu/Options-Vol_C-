#include "Functions.h"
#include <tuple>
using namespace std;
// (const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T, const double d, const double steps)
int main () {

  double S = 7420.40;  // Spot price
  double r = 0.001;   // Risk-free rate (0.1%)
  double T = 30.0/365.0;    // One year until expiry

  vector<double> optionsToPrice = {7450.00, 7500.00, 7550.00,
                              7600.00, 7700.00, 7750.00};

  vector<double> strikePrices = {7425.00, 7650.00,
                                  7700.00, 7800.00, 7825.00};

  vector<double> impVols = {0.123118, 0.122606, 0.109341,
                              0.134298, 0.108216};

  for (int i = 0; i < 6; i++){
    
    double strikeWanted = optionsToPrice[i];


    double volInt = volInterpolation(strikePrices, impVols, strikeWanted);

    tuple<double, double, double, double, double, double, double, double> AsianCalculations = asianOption(10000, S, strikeWanted, r, volInt, 1.0, 0.01, 252);
    cout << "Call Price :" << get<0>(AsianCalculations) << "\n";
    cout << "Call Price STD:" << get<1>(AsianCalculations) << "\n"; 

    cout << "\n";

    cout << "Call Price CV:" << get<2>(AsianCalculations) << "\n";
    cout << "Call Price CV STD:" << get<3>(AsianCalculations) << "\n";

    cout << "\n";

    cout << "Put Price:" << get<4>(AsianCalculations) << "\n";  
    cout << "Put Price Std:" << get<5>(AsianCalculations) << "\n"; 

    cout << "\n";

    cout << "Put Price CV:" << get<6>(AsianCalculations) << "\n"; 
    cout << "Put Price CV STD:" << get<7>(AsianCalculations) << "\n";

    cout << "\n";
}

 return 0;
}