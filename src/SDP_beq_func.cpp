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
    double Y_decay = 5;
    
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
    
    //Build fitness list, jstar list
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
      
      IntegerMatrix dec_df(xmax,tmax-1);
      
      //Iterate backwards over t
      for (int t=(tmax-1); t --> 0;) {
        
        //Itertate over x
        for (int x=xc;x<xmax;x++) {
          
          double x_state = (double) x + 1;
          
          NumericVector value(num_res);
          IntegerVector dec_df_j(num_res);
          IntegerMatrix dec_f(kmax,num_res);
          
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
            
            ////////////
            //Decisions
            ///////////
            
            ////Don't find food, don't eat cache
            ////
            double x_dfdc = x_state - a*pow(Mc,b);
            double theta_dfdc = theta_state - Y_decay;
            ////
            
            ////Don't find food, eat cache
            ////
            double x_dfc = x_state - a*pow(Mc,b);
            double theta_dfc = theta_state - Y_theta - Y_decay;
            ////
            
            //Boundary conditions
            if (x_dfdc < xc_state) {x_dfdc = xc_state;}
            if (x_dfdc > xmax_state) {x_dfdc = xmax_state;}
            if (theta_dfdc < 0) {theta_dfdc=0;}
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
            x_dfdc_low = x_dfdc_low - 1;
            x_dfdc_high = x_dfdc_high - 1;
            theta_dfdc_low = theta_dfdc_low - 1;
            theta_dfdc_low = theta_dfdc_low - 1;
            
            //Define fitness
            NumericMatrix W_theta_low = W_theta(theta_dfdc_low);
            NumericMatrix W_theta_high= W_theta(theta_dfdc_high);
            
            double W_dfdc = qtheta_dfdc*(qx_dfdc*W_theta_low(x_dfdc_low,t+1) + (1 - qx_dfdc)*W_theta_low(x_dfdc_high,t+1)) + 
            (1 - qtheta_dfdc)*(qx_dfdc*W_theta_high(x_dfdc_low,t+1) +  (1 - qx_dfdc)*W_theta_high(x_dfdc_high,t+1));
            
            
            //Boundary conditions
            if (x_dfc < xc_state) {x_dfc = xc_state;}
            if (x_dfc > xmax_state) {x_dfc = xmax_state;}
            if (theta_dfc < 0) {theta_dfc=0;}
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
            x_dfc_low = x_dfc_low - 1;
            x_dfc_high = x_dfc_high - 1;
            theta_dfc_low = theta_dfc_low - 1;
            theta_dfc_high = theta_dfc_high - 1;
            
            //Define fitness
            NumericMatrix W_theta_low = W_theta(theta_dfc_low);
            NumericMatrix W_theta_high= W_theta(theta_dfc_high);
            
            double W_dfc = qtheta_dfc*(qx_dfc*W_theta_low(x_dfc_low,t+1) + (1 - qx_dfc)*W_theta_low(x_dfc_high,t+1)) +
            (1 - qtheta_dfc)*(qx_dfc*W_theta_high(x_dfc_low,t+1) + (1 - qx_dfc)*W_theta_high(x_dfc_high,t+1));
            
            
            //Which maximizes fitness over dfdc and dfc?
            NumericVector fitness_df(2);
            fitness_df(0) = W_dfdc;
            fitness_df(1) = W_dfc;
            int max_dec_df = which_max(fitness_df);
            dec_df_j(j) = max_dec_df;
            
            //Weighted by the probability of k=0 for food type j
            double W_max_df = pkj(0)*fitness_df(max_dec_df);
            
            
            IntegerVector max_dec_f(kmax);
            NumericVector W_max_f(kmax);
            
            
            //Iterate over nonzero values of k
            for (k=1;k<kmax;k++) {
              
              //Not plus 1 bc we are 
              double k_state = (double) k;
              
              
              //Define Y_k
              temp_v(0) = k_state*gainj; temp_v(1)=xs;
              int minvalue = which_min(temp_v);
              double Y_k = temp_v(minvalue);
              double Y_remainder = k_state*gainj - Y_k;
              if (Y_remainder < 0) {Y_remainder = 0;}
              
              //Define Y_remain
              NumericVector temp_remain(2); temp_remain(0)=Y_remainder; temp_remain(1)=xs;
              int minvalue_remain = which_min(temp_remain);
              double Y_remain = temp_remain(minvalue);
              
              ////////////
              //Decisions
              ///////////
              
              ////Find food, reject, don't store, don't eat cache
              //// FDSDC
              double x_fdsdc = x_state - a*pow(Mc,b);
              double theta_fdsdc = theta_state - Y_decay;
              ////
              
              ////Find food, reject, don't store, eat cache
              //// FDSC
              double x_fdsc = x_state - a*pow(Mc,b) + Y_theta;
              double theta_fdsc = theta_state - Y_theta - Y_decay;
              //// 
              
              ////Find food, store, don't each cache
              //// FSDC
              double x_fsdc = x_state - a*pow(Mc,b);
              double theta_fsdc = theta_state + Y_k - Y_decay;
              //// 
              
              ////Find food, store, eat cache
              //// FSC
              double x_fsc = x_state - a*pow(Mc,b) + Y_theta;
              double theta_fsc = theta_state + Y_k - Y_theta - Y_decay;
              ////
              
              ////Find food, accept, store remainder
              //// FS
              double x_fs = x_state - a*pow(Mc,b) + Y_k;
              double theta_fs = theta_state + Y_remain - Y_decay;
              ////
              
              ////find food, accept, don't store remainder
              //// FDS
              double x_fds = x_state - a*pow(Mc,b) + Y_k;
              double theta_fds = theta_state - Y_decay;
              ////
              
              //Boundary conditions
              if (x_fdsdc < xc_state) {x_fdsdc = xc_state;}
              if (x_fdsdc > xmax_state) {x_fdsdc = xmax_state;}
              if (theta_fdsdc < 0) {theta_fdsdc=0;}
              if (theta_fdsdc > thetamax_state) {theta_fdsdc = thetamax_state;}
              
              //Fitness Interpolation
              int x_fdsdc_low = (int) floor(x_fdsdc);
              int x_fdsdc_high = x_fdsdc + 1;
              double x_fdsdc_h = (double) x_fdsdc_high;
              int theta_fdsdc_low = (int) floor(theta_fdsdc);
              int theta_fdsdc_high = (int) theta_fdsdc_low + 1;
              double theta_fdsdc_h = (double) theta_fdsdc_high;
              
              double qx_fdsdc = x_fdsdc_h - x_fdsdc;
              double qtheta_fdsdc = theta_fdsdc_h - theta_fdsdc;
              
              //Adjust to represent indices rather than state
              x_fdsdc_low = x_fdsdc_low - 1;
              x_fdsdc_high = x_fdsdc_high - 1;
              theta_fdsdc_low = theta_fdsdc_low - 1;
              theta_fdsdc_high = theta_fdsdc_high - 1;
              
              //Define fitness
              NumericMatrix W_theta_low = W_theta(theta_fdsdc_low);
              NumericMatrix W_theta_high= W_theta(theta_fdsdc_high);
              
              double W_fdsdc = qtheta_fdsdc*(qx_fdsdc*W_theta_low(x_fdsdc_low,t+1) + (1 - qx_fdsdc)*W_theta_low(x_fdsdc_high,t+1)) +
              (1 - qtheta_fdsdc)*(qx_fdsdc*W_theta_high(x_fdsdc_low,t+1) + (1 - qx_fdsdc)*W_theta_high(x_fdsdc_high,t+1));
              
              
              //Boundary conditions
              if (x_fdsc < xc_state) {x_fdsc = xc_state;}
              if (x_fdsc > xmax_state) {x_fdsc = xmax_state;}
              if (theta_fdsc < 0) {theta_fdsc=0;}
              if (theta_fdsc > thetamax_state) {theta_fdsc = thetamax_state;}
              
              //Fitness Interpolation
              int x_fdsc_low = (int) floor(x_fdsc);
              int x_fdsc_high = x_fdsc + 1;
              double x_fdsc_h = (double) x_fdsc_high;
              int theta_fdsc_low = (int) floor(theta_fdsc);
              int theta_fdsc_high = (int) theta_fdsc_low + 1;
              double theta_fdsc_h = (double) theta_fdsc_high;
              
              double qx_fdsc = x_fdsc_h - x_fdsc;
              double qtheta_fdsc = theta_fdsc_h - theta_fdsc;
              
              //Adjust to represent indices rather than state
              x_fdsc_low = x_fdsc_low - 1;
              x_fdsc_high = x_fdsc_high - 1;
              theta_fdsc_low = theta_fdsc_low - 1;
              theta_fdsc_high = theta_fdsc_high - 1;
              
              //Define fitness
              NumericMatrix W_theta_low = W_theta(theta_fdsc_low);
              NumericMatrix W_theta_high= W_theta(theta_fdsc_high);
              
              double W_fdsc = qtheta_fdsc*(qx_fdsc*W_theta_low(x_fdsc_low,t+1) + (1 - qx_fdsc)*W_theta_low(x_fdsc_high,t+1)) +
              (1 - qtheta_fdsc)*(qx_fdsc*W_theta_high(x_fdsc_low,t+1) + (1 - qx_fdsc)*W_theta_high(x_fdsc_high,t+1));
              
              
              //Boundary conditions
              if (x_fsdc < xc_state) {x_fsdc = xc_state;}
              if (x_fsdc > xmax_state) {x_fsdc = xmax_state;}
              if (theta_fsdc < 0) {theta_fsdc=0;}
              if (theta_fsdc > thetamax_state) {theta_fsdc = thetamax_state;}
              
              //Fitness Interpolation
              int x_fsdc_low = (int) floor(x_fsdc);
              int x_fsdc_high = x_fsdc + 1;
              double x_fsdc_h = (double) x_fsdc_high;
              int theta_fsdc_low = (int) floor(theta_fsdc);
              int theta_fsdc_high = (int) theta_fsdc_low + 1;
              double theta_fsdc_h = (double) theta_fsdc_high;
              
              double qx_fsdc = x_fsdc_h - x_fsdc;
              double qtheta_fsdc = theta_fsdc_h - theta_fsdc;
              
              //Adjust to represent indices rather than state
              x_fsdc_low = x_fsdc_low - 1;
              x_fsdc_high = x_fsdc_high - 1;
              theta_fsdc_low = theta_fsdc_low - 1;
              theta_fsdc_high = theta_fsdc_high - 1;
              
              //Define fitness
              NumericMatrix W_theta_low = W_theta(theta_fsdc_low);
              NumericMatrix W_theta_high= W_theta(theta_fsdc_high);
              
              double W_fsdc = qtheta_fsdc*(qx_fsdc*W_theta_low(x_fsdc_low,t+1) + (1 - qx_fsdc)*W_theta_low(x_fsdc_high,t+1)) +
              (1 - qtheta_fsdc)*(qx_fsdc*W_theta_high(x_fsdc_low,t+1) + (1 - qx_fsdc)*W_theta_high(x_fsdc_high,t+1));
              
              
              //Boundary conditions
              if (x_fsc < xc_state) {x_fsc = xc_state;}
              if (x_fsc > xmax_state) {x_fsc = xmax_state;}
              if (theta_fsc < 0) {theta_fsc=0;}
              if (theta_fsc > thetamax_state) {theta_fsc = thetamax_state;}
              
              //Fitness Interpolation
              int x_fsc_low = (int) floor(x_fsc);
              int x_fsc_high = x_fsc + 1;
              double x_fsc_h = (double) x_fsc_high;
              int theta_fsc_low = (int) floor(theta_fsc);
              int theta_fsc_high = (int) theta_fsc_low + 1;
              double theta_fsc_h = (double) theta_fsc_high;
              
              double qx_fsc = x_fsc_h - x_fsc;
              double qtheta_fsc = theta_fsc_h - theta_fsc;
              
              //Adjust to represent indices rather than state
              x_fsc_low = x_fsc_low - 1;
              x_fsc_high = x_fsc_high - 1;
              theta_fsc_low = theta_fsc_low - 1;
              theta_fsc_high = theta_fsc_high - 1;
              
              //Define fitness
              NumericMatrix W_theta_low = W_theta(theta_fsc_low);
              NumericMatrix W_theta_high= W_theta(theta_fsc_high);
              
              double W_fsc = qtheta_fsc*(qx_fsc*W_theta_low(x_fsc_low,t+1) + (1 - qx_fsc)*W_theta_low(x_fsc_high,t+1)) +
              (1-qtheta_fsc)*(qx_fsc*W_theta_high(x_fsc_low,t+1) + (1 - qx_fsc)*W_theta_high(x_fsc_high,t+1));
              
              
              //Boundary conditions
              if (x_fs < xc_state) {x_fs = xc_state;}
              if (x_fs > xmax_state) {x_fs = xmax_state;}
              if (theta_fs < 0) {theta_fs=0;}
              if (theta_fs > thetamax_state) {theta_fs = thetamax_state;}
              
              //Fitness Interpolation
              int x_fs_low = (int) floor(x_fs);
              int x_fs_high = x_fs + 1;
              double x_fs_h = (double) x_fs_high;
              int theta_fs_low = (int) floor(theta_fs);
              int theta_fs_high = (int) theta_fs_low + 1;
              double theta_fs_h = (double) theta_fs_high;
              
              double qx_fs = x_fs_h - x_fs;
              double qtheta_fs = theta_fs_h - theta_fs;
              
              //Adjust to represent indices rather than state
              x_fs_low = x_fs_low - 1;
              x_fs_high = x_fs_high - 1;
              theta_fs_low = theta_fs_low - 1;
              theta_fs_high = theta_fs_high - 1;
              
              //Define fitness
              NumericMatrix W_theta_low = W_theta(theta_fs_low);
              NumericMatrix W_theta_high= W_theta(theta_fs_high);
              
              double W_fs = qtheta_fs*(qx_fs*W_theta_low(x_fs_low,t+1) + (1 - qx_fs)*W_theta_low(x_fs_high,t+1)) +
              (1 - qtheta_fs)*(qx_fs*W_theta_high(x_fs_low,t+1) + (1 - qx_fs)*W_theta_high(x_fs_high,t+1));
              
              
              //Boundary conditions
              if (x_fds < xc_state) {x_fds = xc_state;}
              if (x_fds > xmax_state) {x_fds = xmax_state;}
              if (theta_fds < 0) {theta_fds=0;}
              if (theta_fds > thetamax_state) {theta_fds = thetamax_state;}
              
              //Fitness Interpolation
              int x_fds_low = (int) floor(x_fds);
              int x_fds_high = x_fds + 1;
              double x_fds_h = (double) x_fds_high;
              int theta_fds_low = (int) floor(theta_fds);
              int theta_fds_high = (int) theta_fds_low + 1;
              double theta_fds_h = (double) theta_fds_high;
              
              double qx_fds = x_fds_h - x_fds;
              double qtheta_fds = theta_fds_h - theta_fds;
              
              //Adjust to represent indices rather than state
              x_fds_low = x_fds_low - 1;
              x_fds_high = x_fds_high - 1;
              theta_fds_low = theta_fds_low - 1;
              theta_fds_high = theta_fds_high - 1;
              
              //Define fitness
              NumericMatrix W_theta_low = W_theta(theta_fds_low);
              NumericMatrix W_theta_high= W_theta(theta_fds_high);
              
              double W_fds = qtheta_fds*(qx_fds*W_theta_low(x_fds_low,t+1) + (1 - qx_fds)*W_theta_low(x_fds_high,t+1)) +
              (1 - qtheta_fds)*(qx_fds*W_theta_high(x_fds_low,t+1) + (1 - qx_fds)*W_theta_high(x_fds_high,t+1));
              
              
              //Find maximum of {W(xfc), W(xfsc), W(xf)}
              NumericVector fitness_f(6);
              fitness_f(0) = W_fdsdc;
              fitness_f(1) = W_fdsc;
              fitness_f(2) = W_fsdc;
              fitness_f(3) = W_fsc;
              fitness_f(4) = W_fs;
              fitness_f(5) = W_fds;
              
              //Find the maximum
              max_dec_f(k) = which_max(fitness_f);
              dec_f(k,j) = max_dec_f(k);
              
              //Fitness value weighted by the probability of k
              W_max_f(k) = pkj(k)*fitness_f(max_dec_f(k));
              
              
            } //End k
            
            
            
            
            //Weight fitness by pr(k)
            double W_k0 = W_max_df;
            double W_k = sum(W_max_f);
            
            value(j) = W_k0 + W_k;
            
          } //End j
          
          //Find maximum value over j
          int maxvalue = which_max(value);
          //Store jstar
          jstar_theta(theta)(x,t) = maxvalue;
          
          //The maximizing decision for the maximizing food if k=0
          dec_df(x,t) = dec_df_j(maxvalue);
          //MORE COMPLICATED FOR K>0... think about this...
          
          //Update fitness matrix
          W_theta(theta)(x,t) = value(maxvalue);
          
        } //End x
        
        
      } //End t iterations
      
      
    } //End theta iterations
    
    
    //List Cout(5);
    return W_theta;
  }
  