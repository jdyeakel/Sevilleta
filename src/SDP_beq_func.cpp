#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List SDP_beq_func(
  double Mc,
  int tmax,
  NumericMatrix pk,
  NumericVector gain,
  NumericVector cost,
  IntegerVector seed
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
    //Stomach capacity
    double xs = 0.10*xmax_state;
    
    double thetamax_state = 500;
    int thetamax = (int) thetamax_state;
    
    
    //Food qualities
    int num_res = gain.size();
    int maxk = pk.size(1);
    
    //Build terminal fitness function
    //It will be important to test result sensitivities to this function...
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
    
    //Build fitness list, istar list
    List W_theta(thetamax);
    List jstar_theta(thetamax);
    for (int i=0;i<thetamax;i++) {
      NumericMatrix W_xt(xmax,tmax);
      NumericMatrix jstar_xt(xmax,tmax-1);
      for (int j=0;j<xmax;j++) {
        W_xt(j,tmax-1) = Wterm(j,i);
      }
      W_theta(i) = W_xt;
      jstar_theta(i) = jstar_xt;
    }
    
    
    //Begin backwards equation over theta
    for (int theta=0;theta<thetamax;theta++) {
      
      double theta_state = (double) theta;
      
      NumericMatrix W_xt = W_theta(theta);
      NumericMatrix istar_xt = istar_theta(theta);
      
      //Iterate backwards over t
      for (int t=(tmax-1); t --> 0;) {
        
        //Itertate over x
        for (int x=xc;x<xmax;x++) {
          
          double x_state = (double) x + 1;
          
          NumericVector value(num_res);
          
          //Iterate over potential foods j
          for(int j=0;j<num_res;j++) {
            
            double gainj = gain(j);
            
            NumericVector pkj(maxk);
            for (int k=0;k<maxk;k++) {
              pkj(k) = pk(k,j);
            }
            
            
            //Define Y_theta
            NumericVector temp_v(2); temp_v(0)=theta_state; temp_v(1)=xs;
            int minvalue = which_min(temp_v);
            double Y_theta = temp_v(minvalue);
            
            //The case of k=0
            if (theta > 0) {
              
              ////Don't find food, don't eat cache
              double x_dfdc = x_state - a*pow(Mc,b);
              double theta_dfdc = theta_state;
              
              //Boundary conditions
              if (x_dfdc < xc_state) {x_dfdc = xc_state;}
              if (x_dfdc > xmax_state) {x_dfdc = xmax_state;}
              if (theta_dfdc < 0) {theta_dfdc == 0;}
              if (theta_dfdc > thetamax_state) {theta_dfdc = thetamax_state;}
              
              //Fitness Interpolation
              int x_dfdc_low = (int) floor(x_dfdc);
              int x_dfdc_high = x_dfdc + 1;
              double x_dfdc_h = (double) x_dfdc_high;
              int theta_dfdc_low = (int) floor(theta_dfdc);
              int theta_dfdc_high = (int) theta_dfdc_low + 1;
              double theta_dfdc_h = (double) theta_dfdc_high;
              
              //Interpolation weights
              double qx_dfdc = x_dfdc_h - x_dfdc;
              double qtheta_dfdc = theta_dfdc_h - theta_dfdc;
              //Adjust to represent indices rather than state
              x_dfdc = x_dfdc - 1;
              theta_dfdc = theta_dfdc - 1;
              
              //Define fitness
              W_theta_low = W_theta(theta_dfdc_low);
              W_theta_high= W_theta(theta_dfdc_high);
              
              W_dfdc = qtheta_dfdc_h*(qx_dfdc*W_theta_low(x_dfdc_low,t+1) + 
              (1 - qx_dfdc)*W_theta_low(x_dfdc_high,t+1)) + 
              qtheta_dfdc_l*(qx_dfdc*W_theta_high(x_dfdc_low,t+1) + 
              (1 - qx_dfdc)*W_theta_high(x_dfdc_high,t+1));
              
              
              ////Don't find food, eat cache
              double x_dfc = x_state - a*pow(Mc,b);
              double theta_dfc = theta_state - Y_theta;
              
              //Boundary conditions
              if (x_dfc < xc_state) {x_dfc = xc_state;}
              if (x_dfc > xmax_state) {x_dfc = xmax_state;}
              if (theta_dfc < 0) {theta_dfc == 0;}
              if (theta_dfc > thetamax_state) {theta_dfc = thetamax_state;}
              
              //Fitness Interpolation
              int x_dfc_low = (int) floor(x_dfc);
              int x_dfc_high = x_dfc + 1;
              double x_dfc_h = (double) x_dfc_high;
              int theta_dfc_low = (int) floor(theta_dfc);
              int theta_dfc_high = (int) theta_dfc_low + 1;
              double theta_dfc_h = (double) theta_dfc_high;
              
              double qx_dfc = x_dfc_h - x_dfc;
              double qtheta_dfc = theta_dfc_h - theta_dfc;
              
              //Adjust to represent indices rather than state
              x_dfc = x_dfc - 1;
              theta_dfc = theta_dfc - 1;
              
              W_dfc = qtheta_dfc_h*(qx_dfc*W_theta_low(x_dfc_low,t+1) + 
              (1 - qx_dfc)*W_theta_low(x_dfc_high,t+1)) + 
              qtheta_dfc_l*(qx_dfc*W_theta_high(x_dfc_low,t+1) + 
              (1 - qx_dfc)*W_theta_high(x_dfc_high,t+1));

              //Which maximizes fitness over dfdc and dfc?
              NumericVector fitness_df(2);
              temp_v(0) = W_dfdc;
              temp_v(1) = W_dfc;
              int max_dec_df = which_max(fitness_df);
            } else {
              //If there is no cache, can't each from cache
              //Automatically the decision is dfdc
              max_dec_df = 0;
            }
            
            //Iterate over nonzero values of k
            for (k=1;k<kmax;k++) {
              
              //Not plus 1 bc we are 
              double k_state = (double) k;
              
              
              //Define Y_k
              temp_v(0) = k_state*gainj; temp_v(1)=xs;
              int minvalue = which_min(temp_v);
              double Y_k = temp_v(minvalue);
              
              if (theta_state > 0) {
                
              }
              //Different foraging decisions:
              //Find food, reject, don't store, eat cache
              x_fc
              
              //Find food, reject, store, eat cache
              x_fsc
              
              //Find food, accept
              x_f
              
              //Boundary Conditions
              
              //Interpolation
              
              
              
              //Find maximum of {W(xfc), W(xfsc), W(xf)}
              NumericVector fitness_f(5);
              temp_v(0) = W_fdsdc;
              temp_v(1) = W_fdsc;
              temp_v(2) = W_fsdc;
              temp_v(3) = W_fsc;
              temp_v(4) = W_f;
              int max_dec_f = which_max(fitness_f);
              
              
            } //End k
            
            
            
            
            
            
          } //End j
          
        } //End x
        
        
        
        
        
      } //End t iterations
      
      
      
    } //End theta iterations
    
    
    
    
    //List Cout(5);
    return W_theta;
  }
  