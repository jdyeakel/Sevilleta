#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
NumericMatrix SDP_beq_func(
  NumericMatrix pk,
  double Mc
) {
  
  
  //Define state variable max's and min's
  //'State' refers to the value used in the calculation...
  //Otherwise, the value is an integer and used in locating fitness values.
  
  //2750 kcal/kg
  //~10 MJ in 2750 kcal.... so 10 MJ per kg.
  double xmax_state = Mc*10; //Body size in kilograms to MegaJoules
  int xmax = (int) floor(xmax_state);
  double xc_state = 0.5*xmax_state;
  int xc = (int) floor(xc_state);
  
  double thetamax_state = 500;
  int thetamax = (int) thetamax_state;
  
  NumericMatrix Wterm(xmax,thetamax);
  for (int i=xc;i<xmax;i++) {
    for (int j=0;j<thetamax;j++) {
      double x_state = (double) i+1;
      double theta_state = (double) j+1;
      Wterm(i,j) = 0.5*(2 - (xc_state/x_state) - (1/((0.01*theta_state)+1)));
    }
  }
  double Wterm_max = Wterm(xmax-1,thetamax-1);
  for (int i=xc;i<xmax;i++) {
    for (int j=0;j<thetamax;j++) {
      Wterm(i,j) = Wterm(i,j)/Wterm_max;
    }
  }
  
  
  
  
  List Cout(5);
  return Wterm;
}
