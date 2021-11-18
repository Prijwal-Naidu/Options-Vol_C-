#ifndef __Functions_H
#define __Functions_H

#include <algorithm>    
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>
#include <numeric>

#include "spline.h"

using namespace std;

// A simple implementation of the Box-Muller algorithm, used to generate
// gaussian random numbers - necessary for the Monte Carlo method 

double gaussian_box_muller() {
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  // Continue generating two uniform random variables
  // until the square of their "euclidean distance" 
  // is less than unity
  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX)-1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2*log(euclid_sq)/euclid_sq);
}

double norm_pdf(const double& x) {
    return (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x);
}

// An approximation to the cumulative distribution function
// for the standard normal distribution
// Note: This is a recursive function
double norm_cdf(const double& x) {
    double k = 1.0/(1.0 + 0.2316419*x);
    double k_sum = k*(0.319381530 + k*(-0.356563782 + k*(1.781477937 + k*(-1.821255978 + 1.330274429*k))));

    if (x >= 0.0) {
        return (1.0 - (1.0/(pow(2*M_PI,0.5)))*exp(-0.5*x*x) * k_sum);
    } else {
        return 1.0 - norm_cdf(-x);
    }
}

// This calculates d_j, for j in {1,2}. This term appears in the closed
// form solution for the European call or put price
double d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
    return (log(S/K) + (r + (pow(-1,j-1))*0.5*v*v)*T)/(v*(pow(T,0.5)));
}

// Calculate the European vanilla call price based on
// underlying S, strike K, risk-free rate r, volatility of
// underlying sigma and time to maturity T
double call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return S * norm_cdf(d_j(1, S, K, r, v, T))-K*exp(-r*T) * norm_cdf(d_j(2, S, K, r, v, T));
}
// call_price = 9.6482

// Calculate the European vanilla put price based on
// underlying S, strike K, risk-free rate r, volatility of
// underlying sigma and time to maturity T
double put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
    return -S*norm_cdf(-d_j(1, S, K, r, v, T))+K*exp(-r*T) * norm_cdf(-d_j(2, S, K, r, v, T));
}

tuple<double, double> optionBS(const double& S, const double& K, const double& r, const double& v, const double& T){
    double callPrice = call_price(S, K, r, v, T);
    double putPrice = put_price(S, K, r, v, T);

    return make_tuple(callPrice, putPrice);
}

double vega(double S, double K, double T, double v, double r)
{
	const double PI = 3.141592;
	double d1 = (log(S/K) + (r + (pow(v, 2)/2)*T))/(v*sqrt(T));
	double vegaC = (1.0/sqrt(2*PI))*(S*exp(-T))*(exp(-(pow(d1, 2))))*sqrt(T);
	return vegaC;
}

// Pricing a European vanilla call option with a Monte Carlo method
tuple<double, double, double, double, double, double,double> monte_carlo_call_price(const int& num_sims, const double& S, const double& K, 
                                      const double& r, const double& v, const double& T, const double d) {
    double S_pv = S * exp(T*(r - 0.5*v*v));
    double S_cur = 0.0;
    double call_price_sum = 0.0;

    //this is for the greeks
    double S_pv_change_S_up = (S+d) * exp(T*(r - 0.5*v*v));
    double S_cur_change_S_up = 0.0;
    double call_price_S_up_sum = 0.0;
    double S_pv_change_S_down = (S-d) * exp(T*(r - 0.5*v*v));
    double S_cur_change_S_down = 0.0;
    double call_price_S_down_sum = 0.0;

    double S_pv_change_v_up = S * exp(T*(r - 0.5*(v+d)*(v+d)));
    double S_cur_change_v_up = 0.0;
    double call_price_v_up_sum = 0.0;
    double S_pv_change_v_down = S * exp(T*(r - 0.5*(v-d)*(v-d)));
    double S_cur_change_v_down = 0.0;
    double call_price_v_down_sum = 0.0;

    double S_pv_change_T_up = S * exp((T+d)*(r - 0.5*v*v));
    double S_cur_change_T_up = 0.0;
    double call_price_T_up_sum = 0.0;
    double S_pv_change_T_down = S * exp((T-d)*(r - 0.5*v*v));
    double S_cur_change_T_down = 0.0;
    double call_price_T_down_sum = 0.0;

    double S_pv_change_r_up = S * exp(T*((r+d) - 0.5*v*v));
    double S_cur_change_r_up = 0.0;
    double call_price_r_up_sum = 0.0;
    double S_pv_change_r_down = S * exp(T*((r-d) - 0.5*v*v));
    double S_cur_change_r_down = 0.0;
    double call_price_r_down_sum = 0.0;

      

    for (int i=0; i<num_sims; i++){
      double gaus_bm = gaussian_box_muller();
      S_cur = S_pv * exp(sqrt(v*v*T)*gaus_bm);
      call_price_sum += max(-K+ S_cur, 0.0);

      S_cur_change_S_up = S_pv_change_S_up * exp(sqrt(v*v*T)*gaus_bm);
      S_cur_change_S_down = S_pv_change_S_down * exp(sqrt(v*v*T)*gaus_bm);

      S_cur_change_v_up = S_pv_change_v_up * exp(sqrt((v+d)*(v+d)*T)*gaus_bm);
      S_cur_change_v_down = S_pv_change_v_down * exp(sqrt((v-d)*(v-d)*T)*gaus_bm);

      S_cur_change_T_up = S_pv_change_T_up * exp(sqrt((T+d)*v*v)*gaus_bm);
      S_cur_change_T_down = S_pv_change_T_down * exp(sqrt((T+d)*v*v)*gaus_bm);

      S_cur_change_r_up = S_pv_change_r_up * exp(sqrt(v*v*T)*gaus_bm);
      S_cur_change_r_down = S_pv_change_r_down * exp(sqrt(v*v*T)*gaus_bm);

      call_price_S_up_sum += max(-K+ S_cur_change_S_up, 0.0);
      call_price_S_down_sum += max(-K+ S_cur_change_S_down, 0.0);

      call_price_v_up_sum += max(-K+ S_cur_change_v_up, 0.0);
      call_price_v_up_sum += max(-K+ S_cur_change_v_down, 0.0);

      call_price_T_up_sum += max(-K+ S_cur_change_T_up, 0.0);
      call_price_T_down_sum += max(-K+ S_cur_change_T_down, 0.0);

      call_price_r_up_sum += max(-K+ S_cur_change_r_up, 0.0);
      call_price_r_down_sum += max(-K+ S_cur_change_r_down, 0.0);

      

    };

    double call_price = (call_price_sum /num_sims)* exp(-r*T);

    double call_delta =  (call_price_S_up_sum - call_price_S_down_sum)/(2*d) * exp(-r*T) / num_sims;
    double call_gamma = (call_price_S_up_sum + call_price_S_down_sum - 2*call_price_sum)/(d*d) * exp(-r*T) / num_sims;
    double call_vega = sqrt((call_price_v_up_sum - call_price_v_down_sum)/(2*d*100) * d*exp(-r*T) / num_sims);

    // we have to individually calculate the rho and theta prices as exp(-r*T) will be different 

    double call_price_T_up = call_price_T_up_sum / num_sims * exp(-r*(T+d));
    double call_price_T_down = call_price_T_down_sum / num_sims * exp(-r*(T-d));

    double call_price_r_up = call_price_r_up_sum / num_sims * exp(-(r+d)*T);
    double call_price_r_down = call_price_r_down_sum / num_sims * exp(-(r-d)*T);

    double call_theta = (call_price_T_up -  call_price_T_down)/(d*2*100);
    double call_rho = (call_price_r_up - call_price_r_down)/(d*2*100);

    return make_tuple(call_price, call_delta, call_gamma, call_vega, call_theta, call_rho, d);

}


double impliedVolatility(const double& S, const double& K, const double& r, const double& T, double tol, double C){
  double ivOld = 1.00; //initial estimate
  double ivNew = 0.5;
  int count = 0;
  while (abs(ivNew - ivOld) > tol && count<10) {  

      ivOld = ivNew;
      ivNew = ivNew - ((call_price(S, K, r, ivNew, T) - C)/vega(S, K, T, ivNew, r));

      count ++;
  }

  return ivNew;
}

double variance(vector<double> samples)
{
     int size = samples.size();

     double variance = 0;
     double t = samples[0];
     for (int i = 1; i < size; i++)
     {
          t += samples[i];
          double diff = ((i + 1) * samples[i]) - t;
          variance += (diff * diff) / ((i + 1.0) *i);
     }

     return variance / (size - 1);
}

tuple<double, double, double> numberVector_description(vector<double> input_vector) {
  double sum = accumulate(input_vector.begin(), input_vector.end(), 0.0);
  double mean = sum/input_vector.size();
  double std_dev = sqrt(variance(input_vector));

  return make_tuple(sum, mean, std_dev);

}

double mean_array(double arr[], int n){
   double sum = 0;
   for(int i = 0; i < n; i++)
   sum = sum + arr[i];
   return sum / n;
}

double covariance(vector<double> vec1, vector<double> vec2)
{   
    int n = vec1.size();

    double* arr1 = &vec1[0];
    double* arr2 = &vec2[0];

    double sum = 0.0;
    for(int i = 0; i < n; i++)
        sum = sum + (arr1[i] - mean_array(arr1, n)) *
                    (arr2[i] - mean_array(arr2, n));
    return sum / (n - 1);
}

tuple<double, double, double, double, double, double, double, double> asianOption(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T, const double d, const double steps) {
 
  double dt = T/steps; //time step
  double drift = exp(dt*(r-0.5*v*v)); //drift
  double vol = sqrt(v*v*dt); //vol per step

  double callPrice = 0.0; //init call price for asian option
  double putPrice = 0.0; //init put price for asian option

  double callPrice_vanilla = 0.0; //init call price for vanilla option
  double putPrice_vanilla = 0.0; //init put price for vanilla option

  double price_T = 0.0; //init the price at expiry
  
  double arithmeticMean;
  
  vector<double> spot_path; // create vector that will contain the price of the underlying at each timestep

  vector<double> callPrices_vector; // create vector of callPrices for each sim
  vector<double> putPrices_vector; //create vector of put prices for each sim

  vector<double> callPrices_vanilla_vector; // create vector of callPrices for each sim
  vector<double> putPrices_vanilla_vector; //create vector of put prices for each sim


  for (int j=0; j <num_sims; j++) { 

    //vector<double> spot_path = {S*1.0};
    spot_path.clear();
    spot_path = {S*1.0};

    for (int i=1; i<steps; i++) {

        double gaus_bm = gaussian_box_muller();
        spot_path.push_back(spot_path.back()  * drift * exp(vol*gaus_bm)); //generates the next price in the time series 
      
      };

    
    double sum = accumulate(spot_path.begin(), spot_path.end(), 0.0); //get the sum of the vector
    double arithmeticMean = sum / steps; //get mean of the vector

    price_T = spot_path.back(); // get the price at time T

    // double arithmeticMean = sum ;
  
    callPrice += max(arithmeticMean - K, 0.0)* exp(-r*T); //call price of the asian option being added
    putPrice += max(-arithmeticMean + K, 0.0)* exp(-r*T); // same as above for put

    callPrice_vanilla += max(price_T - K, 0.0)* exp(-r*T); //cal price of vanilla option (euro)
    putPrice_vanilla += max(-price_T + K, 0.0)* exp(-r*T);  // put price of vanilla option (euro)

    callPrices_vector.push_back(max(arithmeticMean - K, 0.0)* exp(-r*T)); //call option price of asian into a vector
    putPrices_vector.push_back(max(-arithmeticMean + K, 0.0)* exp(-r*T)); //put option price of asian into a vector

    callPrices_vanilla_vector.push_back(max(arithmeticMean - K, 0.0)* exp(-r*T)); //call option price of asian into a vector
    putPrices_vanilla_vector.push_back(max(-arithmeticMean + K, 0.0)* exp(-r*T)); //put option price of asian into a vector
  };

    tuple<double, double, double> callPrice_descriptors = numberVector_description(callPrices_vector); //this function provides a tuple of sum, mean, std
    double callPrice_mean = get<1>(callPrice_descriptors);
    double callPrice_std = get<2>(callPrice_descriptors);

    tuple<double, double, double> putPrice_descriptors = numberVector_description(putPrices_vector);
    double putPrice_mean = get<1>(putPrice_descriptors);
    double putPrice_std = get<2>(putPrice_descriptors);

    //same as above but for vanilla optiond

    tuple<double, double, double> callPrice_vanilla_descriptors = numberVector_description(callPrices_vanilla_vector); 
    double callPrice_vanilla_mean = get<1>(callPrice_vanilla_descriptors);
    double callPrice_vanilla_std = get<2>(callPrice_vanilla_descriptors);

    tuple<double, double, double> putPrice_vanilla_descriptors = numberVector_description(putPrices_vanilla_vector);
    double putPrice_vanilla_mean = get<1>(putPrice_vanilla_descriptors);
    double putPrice_vanilla_std = get<2>(putPrice_vanilla_descriptors);


    callPrice = callPrice/num_sims; //gives us call price of asian option
    putPrice = putPrice/num_sims; //gives us put price of asian option

    callPrice_vanilla = callPrice_vanilla / num_sims; //same as above but for vanilla
    putPrice_vanilla = putPrice_vanilla / num_sims; // same as above but for vanilla

    double callBeta = covariance(callPrices_vector, callPrices_vanilla_vector) / variance(callPrices_vector);
    double callPrice_CV =  call_price(S, K, r, v, T) + callBeta*(callPrice - callPrice_vanilla);
    double callPrice_CV_std = sqrt((variance(callPrices_vector) + variance(callPrices_vanilla_vector) + 2*callBeta * covariance(callPrices_vector, callPrices_vanilla_vector) + callBeta*callBeta*variance(callPrices_vector))/num_sims);

    double putBeta = covariance(putPrices_vector, putPrices_vanilla_vector) / variance(putPrices_vector);
    double putPrice_CV =  put_price(S, K, r, v, T) + putBeta*(putPrice - putPrice_vanilla);
    double putPrice_CV_std = sqrt((variance(putPrices_vector) + variance(putPrices_vanilla_vector) + 2*putBeta * covariance(putPrices_vector, putPrices_vanilla_vector) + putBeta*putBeta*variance(putPrices_vector))/num_sims);
  

    return make_tuple(callPrice, callPrice_std, callPrice_CV, callPrice_CV_std, putPrice, putPrice_std, putPrice_CV, putPrice_CV_std);
}

double volInterpolation(vector<double> strikePrices, vector<double> impVols, double strikeWanted) {

  tk::spline s(strikePrices, impVols);
  double x = strikeWanted, y = s(x), deriv=s.deriv(1,x);

  return y;

}

#endif


