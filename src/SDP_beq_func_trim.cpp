#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar

// [[Rcpp::export]]
List SDP_beq_func_trim(
  double Mc,
  double a,
  double b,
  int thetamax,
  int tmax,
  NumericMatrix pk,
  NumericVector gain) {


    //IntegerVector seed
    //Maybe build this in later on.
    //seed: you can only store if the food is a seed...


    //Define state variable max's and min's
    //'State' refers to the value used in the calculation...
    //Otherwise, the value is an integer and used in locating fitness values.

    //2750 kcal/kg
    //~10 MJ in 2750 kcal.... so 10 MJ per kg.
    double xmax_state = Mc*10; //Body size in kilograms to MegaJoules
    int xmax = (int) floor(xmax_state);
    double xc_state = (double) floor(0.5*xmax_state); //The critical value
    int xc = (int) xc_state - 1; //The index for the critical value

    //Stomach/cheek capacity
    double xs = 0.10*xmax_state;


    //If thetamax is 50, thetamax_state is 49...
    //The vector is 50 units long... but ranges from theta=0 to theta=49.
    //So in theta[50], your chache is just 49 units.
    double thetamax_state = thetamax - 1;

    double Y_decay = 2;

    //Food qualities
    int num_res = gain.size();
    int kmax = pk.nrow();

    Rcout << "xmax " << xmax << std::endl;
    Rcout << "xc " << xc << std::endl;
    Rcout << "xs " << xs << std::endl;

    //Build terminal fitness function
    //It will be important to test result sensitivities to this function...
    NumericMatrix Wterm(xmax,thetamax);
    for (int i=xc+1;i<xmax;i++) {
      for (int j=0;j<thetamax;j++) {
        double x_state = (double) i+1;
        double theta_state = (double) j+1;
        //Wterm(i,j) = 0.5*(2 - (xc_state/x_state) - (1/((0.01*theta_state)+1)));
        //Appears to work better than the above terminal time function...
        Wterm(i,j) = 0.5*(2 - (xc_state/x_state) - (1/(theta_state+1)));
      }
    }
    double Wterm_max = Wterm(xmax-1,thetamax-1);
    for (int i=xc;i<xmax;i++) {
      for (int j=0;j<thetamax;j++) {
        Wterm(i,j) = Wterm(i,j)/Wterm_max;
      }
    }

    Rcout << "Fitness at the critical value: " << Wterm(xc_state,0) << std::endl;

    //Build fitness list, jstar list
    List W_time(tmax);
    List jstar_time(tmax-1);

    W_time(tmax-1) = Wterm;

    for (int t=0;t<(tmax-1);t++) {
      NumericMatrix W_xt(xmax,thetamax);
      NumericMatrix jstar_xt(xmax,thetamax);
      W_time(t) = W_xt;
      jstar_time(t) = jstar_xt;
    }

    //Saving decision rules
    List dec_time(tmax-1);


    //Iterate backwards over t
    for (int t=(tmax-1); t --> 0;) {

      //Values for the future time... starting with term time...
      NumericMatrix W_fut = W_time(t+1);

      //X,Theta matrix for time t... these will initially be empty
      NumericMatrix W_xt = W_time(t);
      NumericMatrix jstar_xt = jstar_time(t);



      //Saving decision rules
      List dec_theta(thetamax);

      //Over theta
      for (int theta=0;theta<thetamax;theta++) {

        //define theta
        //do not add +1 because it starts at zero
        double theta_state = (double) theta;

        //Saving decision values
        List dec_x(xmax);

        //Itertate over x
        for (int x=xc+1;x<xmax;x++) {

          //Add one so that we range through the minimal x_state (xc + 1)
          //to the maximum x_state (xmax).
          //Then it is assumed that xc results in mortality (0 fitness)
          double x_state = (double) x + 1;

          NumericVector value(num_res);

          //Saving decision values
          IntegerMatrix dec(kmax,num_res);

          //Iterate over potential foods j
          for(int j=0;j<num_res;j++) {

            double gainj = gain(j);

            NumericVector pkj(kmax);
            for (int k=0;k<kmax;k++) {
              pkj(k) = pk(k,j);
            }


            //Define Y_theta
            NumericVector temp_v(2); temp_v(0)=theta_state; temp_v(1)=xs;
            int minvalue = which_min(temp_v);
            double Y_theta = temp_v(minvalue); //What you can gain from cache

            //The case of k=0

            ////////////
            //Decisions
            ///////////


            ////Don't find food, don't eat cache
            //// DFDC
            // double x_dfdc = x_state - a*pow(Mc,b);
            // double theta_dfdc = theta_state - Y_decay;
            ////

            ////Don't find food, eat cache
            //// DFC
            double x_dfc = x_state - a*pow(Mc,b) + Y_theta;
            double theta_dfc = theta_state - Y_theta - Y_decay;
            ////
            //
            // //Boundary conditions
            // if (x_dfdc < xc_state) {x_dfdc = xc_state;}
            // if (x_dfdc >= xmax_state) {x_dfdc = xmax_state-0.001;}
            // if (theta_dfdc < 0) {theta_dfdc=0;}
            // if (theta_dfdc >= thetamax_state) {theta_dfdc = thetamax_state-0.001;}
            //
            // //Fitness Interpolation
            // int x_dfdc_low = (int) floor(x_dfdc);
            // int x_dfdc_high = x_dfdc + 1;
            // double x_dfdc_h = (double) x_dfdc_high;
            // int theta_dfdc_low = (int) floor(theta_dfdc);
            // int theta_dfdc_high = (int) theta_dfdc_low + 1;
            // double theta_dfdc_h = (double) theta_dfdc_high;
            //
            // //Interpolation weights
            // double qx_dfdc = x_dfdc_h - x_dfdc;
            // double qtheta_dfdc = theta_dfdc_h - theta_dfdc;
            // //Adjust to represent indices rather than state
            // x_dfdc_low = x_dfdc_low - 1;
            // x_dfdc_high = x_dfdc_high - 1;
            // // theta_dfdc_low = theta_dfdc_low - 1;
            // // theta_dfdc_high = theta_dfdc_high - 1;
            //
            //
            // // Rcout << "theta_dfdc_high = " << theta_dfdc_high << std::endl;
            //
            // //Define fitness
            // double W_dfdc;
            // W_dfdc = qtheta_dfdc*(qx_dfdc*W_fut(x_dfdc_low,theta_dfdc_low) + (1 - qx_dfdc)*W_fut(x_dfdc_high,theta_dfdc_low)) +
            // (1 - qtheta_dfdc)*(qx_dfdc*W_fut(x_dfdc_low,theta_dfdc_high) + (1 - qx_dfdc)*W_fut(x_dfdc_high,theta_dfdc_high));

            // NumericMatrix W_theta_dfdc_low = W_theta(theta_dfdc_low);
            // NumericMatrix W_theta_dfdc_high = W_theta(theta_dfdc_high);
            // W_dfdc = qtheta_dfdc*(qx_dfdc*W_theta_dfdc_low(x_dfdc_low,t+1) + (1 - qx_dfdc)*W_theta_dfdc_low(x_dfdc_high,t+1)) +
            // (1 - qtheta_dfdc)*(qx_dfdc*W_theta_dfdc_high(x_dfdc_low,t+1) +  (1 - qx_dfdc)*W_theta_dfdc_high(x_dfdc_high,t+1));

            //Boundary conditions
            if (x_dfc < xc_state) {x_dfc = xc_state;}
            if (x_dfc >= xmax_state) {x_dfc = xmax_state-0.001;}
            if (theta_dfc < 0) {theta_dfc=0;}
            if (theta_dfc >= thetamax_state) {theta_dfc = thetamax_state-0.001;}

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

            //Define fitness
            //The if statement is to deal with the case: theta_state = theta_max_state, meaning the theta_high will be out of bounds
            double W_dfc;
            W_dfc = qtheta_dfc*(qx_dfc*W_fut(x_dfc_low,theta_dfc_low) + (1 - qx_dfc)*W_fut(x_dfc_high,theta_dfc_low)) +
            (1 - qtheta_dfc)*(qx_dfc*W_fut(x_dfc_low,theta_dfc_high) + (1 - qx_dfc)*W_fut(x_dfc_high,theta_dfc_high));

            // NumericMatrix W_theta_dfc_low = W_theta(theta_dfc_low);
            // NumericMatrix W_theta_dfc_high = W_theta(theta_dfc_high);
            // W_dfc = qtheta_dfc*(qx_dfc*W_theta_dfc_low(x_dfc_low,t+1) + (1 - qx_dfc)*W_theta_dfc_low(x_dfc_high,t+1)) +
            // (1 - qtheta_dfc)*(qx_dfc*W_theta_dfc_high(x_dfc_low,t+1) +  (1 - qx_dfc)*W_theta_dfc_high(x_dfc_high,t+1));


            //Which maximizes fitness over dfdc and dfc?
            //NumericVector fitness_df(1);
            //fitness_df(0) = W_dfdc;
            //fitness_df(1) = W_dfc;
            double fitness_df = W_dfc;

            // if (fitness_df(0) == fitness_df(1)) {
            //   Rcout << "The df values are the same!" << theta << std::endl;
            // }


            //int max_dec_df = which_max(fitness_df);
            //dec(0,j) = max_dec_df + 1;
            dec(0,j) = 1;
            //Weighted by the probability of k=0 for food type j
            //double W_max_df = pkj(0)*fitness_df(max_dec_df);
            double W_max_df = pkj(0)*fitness_df;


            IntegerVector max_dec_f(kmax);
            NumericVector W_max_f(kmax);



            //Iterate over nonzero values of k
            for (int k=1;k<kmax;k++) {

              //Not plus 1 bc we are
              double k_state = (double) k;


              //Define Y_k
              temp_v(0) = k_state*gainj; temp_v(1)=xs;
              int minvalue = which_min(temp_v);
              double Y_k = temp_v(minvalue);
              // double Y_remainder = (k_state*gainj) - Y_k;
              // if (Y_remainder < 0) {Y_remainder = 0;}

              // //Define Y_remain
              // NumericVector temp_remain(2); temp_remain(0)=Y_remainder; temp_remain(1)=xs;
              // int minvalue_remain = which_min(temp_remain);
              // double Y_remain = temp_remain(minvalue_remain);

              ////////////
              //Decisions
              ///////////

              ////Find food, store, don't each cache
              //// FSDC
              double x_fsdc = x_state - a*pow(Mc,b);
              double theta_fsdc = theta_state + Y_k - Y_decay;
              ////

              ////Find food, accept, store remainder
              //// FS
              double x_fs = x_state - a*pow(Mc,b) + Y_k;
              double theta_fs = theta_state - Y_decay; //+ Y_remain
              ////


              //Boundary conditions
              if (x_fsdc < xc_state) {x_fsdc = xc_state;}
              if (x_fsdc >= xmax_state) {x_fsdc = xmax_state-0.001;}
              if (theta_fsdc < 0) {theta_fsdc=0;}
              if (theta_fsdc >= thetamax_state) {theta_fsdc = thetamax_state-0.001;}

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

              //Define fitness
              //The if statement is to deal with the case: theta_state = theta_max_state, meaning the theta_high will be out of bounds
              double W_fsdc;
              W_fsdc = qtheta_fsdc*(qx_fsdc*W_fut(x_fsdc_low,theta_fsdc_low) + (1 - qx_fsdc)*W_fut(x_fsdc_high,theta_fsdc_low)) +
              (1 - qtheta_fsdc)*(qx_fsdc*W_fut(x_fsdc_low,theta_fsdc_high) + (1 - qx_fsdc)*W_fut(x_fsdc_high,theta_fsdc_high));

              // NumericMatrix W_theta_fsdc_low = W_theta(theta_fsdc_low);
              // NumericMatrix W_theta_fsdc_high = W_theta(theta_fsdc_high);
              // W_fsdc = qtheta_fsdc*(qx_fsdc*W_theta_fsdc_low(x_fsdc_low,t+1) + (1 - qx_fsdc)*W_theta_fsdc_low(x_fsdc_high,t+1)) +
              // (1 - qtheta_fsdc)*(qx_fsdc*W_theta_fsdc_high(x_fsdc_low,t+1) +  (1 - qx_fsdc)*W_theta_fsdc_high(x_fsdc_high,t+1));


              //Boundary conditions
              if (x_fs < xc_state) {x_fs = xc_state;}
              if (x_fs >= xmax_state) {x_fs = xmax_state-0.001;}
              if (theta_fs < 0) {theta_fs=0;}
              if (theta_fs >= thetamax_state) {theta_fs = thetamax_state-0.001;}

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
              // theta_fs_low = theta_fs_low - 1;
              // theta_fs_high = theta_fs_high - 1;

              //Define fitness
              //The if statement is to deal with the case: theta_state = theta_max_state, meaning the theta_high will be out of bounds
              double W_fs;
              W_fs = qtheta_fs*(qx_fs*W_fut(x_fs_low,theta_fs_low) + (1 - qx_fs)*W_fut(x_fs_high,theta_fs_low)) +
              (1 - qtheta_fs)*(qx_fs*W_fut(x_fs_low,theta_fs_high) + (1 - qx_fs)*W_fut(x_fs_high,theta_fs_high));

              // NumericMatrix W_theta_fs_low = W_theta(theta_fs_low);
              // NumericMatrix W_theta_fs_high = W_theta(theta_fs_high);
              // W_fs = qtheta_fs*(qx_fs*W_theta_fs_low(x_fs_low,t+1) + (1 - qx_fs)*W_theta_fs_low(x_fs_high,t+1)) +
              // (1 - qtheta_fs)*(qx_fs*W_theta_fs_high(x_fs_low,t+1) +  (1 - qx_fs)*W_theta_fs_high(x_fs_high,t+1));

              //Alternate version... without FSDC and FDS
              NumericVector fitness_f(2);
              fitness_f(0) = W_fsdc;  //Find, store
              fitness_f(1) = W_fs; //Find, eat


              //When costs are greater than 2, you get similarity errors.
              //Check for similarities
              if (fitness_f(0) == fitness_f(1)) {
                Rcout << "The f values are the same! x = " << x << std::endl;
                Rcout << "The f values are the same! x_fsdc = " << Y_k << std::endl;
              }


              //Find the maximum
              max_dec_f(k) = which_max(fitness_f);

              //Saving this for later
              dec(k,j) = max_dec_f(k) + 2;

              //Fitness value weighted by the probability of k
              W_max_f(k) = pkj(k)*fitness_f(max_dec_f(k));


            } //End k


            double W_k0 = W_max_df;
            double W_k = sum(W_max_f);

            value(j) = W_k0 + W_k;

          } //End j



          //Find maximum value over j
          int maxvalue = which_max(value);

          //Update fitness matrix
          W_xt(x,theta) = value(maxvalue);

          //Store jstar
          jstar_xt(x,theta) = maxvalue + 1; //+1 so that the output reflects resources 1:5

          //The maximizing decision for the maximizing food if k=0
          dec_x(x) = dec;


        } //End x

        //Only updating decision matrix because dec_x is a List
        //within which is embedded dec(k,j)
        dec_theta(theta) = dec_x;

        Rcout << "Theta_state = " << theta_state << std::endl;

      } //End theta iterations

      //Update fitness matrix
      W_time(t) = W_xt;
      //Update the jstar matrix
      jstar_time(t) = jstar_xt;
      //Update the decision matrix
      dec_time(t) = dec_theta;

      Rcout << "Time= " << t << std::endl;

    } //End t iterations


    List cout(3);

    cout(0) = W_time;
    cout(1) = jstar_time;
    cout(2) = dec_time;
    return cout;
  }
