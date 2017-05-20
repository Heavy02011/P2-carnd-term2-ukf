#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5; //30; //##############################################################

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5; //30; //##########################################################

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;
/*
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
*/
  
  // predict maesurement
    //radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  //radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.1;
 
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  
  //set example state
  //VectorXd x = VectorXd(n_x);
  x_ <<   5.7441,
         1.3800,
         2.2049,
         0.5015,
         0.3528;

  //set example covariance matrix
  //MatrixXd P = MatrixXd(n_x, n_x);
/* 
  P_ <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
          -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
           0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
          -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
          -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
*/    
    // initialize state covarianve matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
         
  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 0.2;

  //define spreading parameter
  lambda_ = 3 - n_aug_;
  
  //create augmented mean vector
  x_aug_ = VectorXd(7);

  //create augmented state covariance
  P_aug_ = MatrixXd(7, 7);

  //create sigma point matrix
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  //create matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  
  is_initialized_ = false;

  previous_timestamp_ = 0;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;
  
  //add measurement noise covariance matrix
  S_ = MatrixXd(3, 3);
  
  //mean predicted measurement
  z_pred_ = VectorXd(3);
  
  //create matrix for sigma points in measurement space
  Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  //measurement covariance matrix - laser
  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << std_laspx_*std_laspx_,0,
    0,std_laspy_*std_laspy_;
/*
  R_laser_ << 0.0225, 0,
        0, 0.0225;
*/
  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
                0, 0,std_radrd_*std_radrd_;
/*
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;
  R_radar_ = MatrixXd(3, 3);
*/ 
  
  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  
  // time elapsed between the current and previous measurements
  dt = 0.001;
  
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
   if (!is_initialized_) {
    //cout << "rbx: FusionEKF initializing..." << endl;
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
/*      
    // initialize state covarianve matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;
*/
    // initialize weights
    double weight_0_ = lambda_/(lambda_+n_aug_);
    weights_(0) = weight_0_;
    double weight_ = 0.5/(n_aug_+lambda_);
    for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
      weights_(i) = weight_;      
    }
    
    //set the state with the initial location and zero velocity
    //The state vector x for the CTRV model contains x=[px,p​y,v,psi,psidot]
      
    // first measurement
    cout << "EKF: " << endl;
    //x_ = VectorXd(5);
    x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double x= rho * cos(phi);
      double y= rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);
      double v = sqrt(vx*vx+vy*vy);
      //x_ << x, y, vx, vy,0;
      x_ << x, y, v, 0 ,0;
            

      // ####### check ########

    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {

      //Initialize state.
      double px = meas_package.raw_measurements_[0]; 
      double py = meas_package.raw_measurements_[1];
      double v = 0;
      double psi = 0;
      double psidot = 0;
      // setup state vector
      x_ << px,py,v,psi,psidot;

    }
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = meas_package.timestamp_; // correction # 1
    cout << "initialization: x_" << x_ << endl; //rbxok
    cout << "initialization: P_" << P_ << endl; //diff
    cout << "initialization: weights_" << weights_ << endl; //rbxok
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float noise_ax = 9;
  float noise_ay = 9;

  //compute the time elapsed between the current and previous measurements
  dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
/*  
  // possible improvement: if ( dt > 0.001 )
  float dt_2 = dt * dt; 
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // setup matrix F
  F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;
  //Modify the F matrix so that the time is integrated
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  //set the process covariance matrix Q
  Q_ = MatrixXd(4, 4);
  Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0, 
              0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
              dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
              0, dt_3/2*noise_ay, 0, dt_2*noise_ay;  
  
  Predict();
*/
  
  cout << "ukfProcessMeasurement: prediction" << endl;
  //dt = 0.001;
  Prediction(dt);
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

    if (use_radar_) {
      // Radar updates
      cout << "Updating Radar" << endl;
      //H_ = tools.CalculateJacobian(x_);
      //R_ = R_radar_;
      //cout << "meas_package.raw_measurements_" << meas_package.raw_measurements_ << endl;
      //UpdateEKF(meas_package.raw_measurements_);
      UpdateRadar(meas_package);

    }


  } else {

    if (use_laser_) {
      // Laser updates
      //measurement update
      cout << "Updating Laser" << endl;
      //H_ = H_laser_;
      //R_ = R_laser_;
      //Update(meas_package.raw_measurements_);
    }

  }
 
  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  // 1. generate (augmented) sigma points x(k|k)
  //MatrixXd Xsig = MatrixXd(11, 5);
  //UKF::GenerateSigmaPoints(&Xsig);
  //std::cout << "Xsig = " << std::endl << Xsig << std::endl;
  
  MatrixXd Xsig_aug = MatrixXd(15, 7);
  UKF::AugmentedSigmaPoints(&Xsig_aug);
  //std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
  std::cout << "UKF::Prediction: Xsig_aug_ = " << std::endl << Xsig_aug_ << std::endl; //rbxok
  
  // 2. predict Sigma Points of next time step x(k+1|k)
  MatrixXd Xsig_pred = MatrixXd(15, 5);
  UKF::SigmaPointPrediction(&Xsig_pred);
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
  std::cout << "UKF::Prediction: Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl; // rbxok
  cout << "UKF::Prediction: dt = " << std::endl << dt << std::endl; // rbxok
  
  // 3. predict mean and covariance of next time step
  VectorXd x_pred = VectorXd(5);
  MatrixXd P_pred = MatrixXd(5, 5);
  UKF::PredictMeanAndCovariance(&x_pred, &P_pred);
  //##std::cout << "x_pred = " << std::endl << x_pred << std::endl;
  //##std::cout << "P_pred = " << std::endl << P_pred << std::endl;
  cout << "UKF::Prediction: x_ = " << std::endl << x_ << std::endl;
  cout << "UKF::Prediction: P_ = " << std::endl << P_ << std::endl;
  //cout << "predict:<< std_a_, std_yawdd_   " << std_a_ << "   " << std_yawdd_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  // predict new z
  VectorXd z_out = VectorXd(3);
  MatrixXd S_out = MatrixXd(3, 3);
  //PredictRadarMeasurement(&z_out, &S_out);
  PredictRadarMeasurement(&z_out, &S_out);
/*
  std::cout << "UpdateRadar: z_out: " << std::endl << z_out << std::endl;
  std::cout << "UpdateRadar: S_out: " << std::endl << S_out << std::endl;
  std::cout << "UpdateRadar: z_pred: " << std::endl << z_pred_ << std::endl;
  std::cout << "UpdateRadar: S: " << std::endl << S_ << std::endl;
*/  
  //std::cout << "UpdateRadar: Zsig_: " << std::endl << Zsig_ << std::endl;
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_ = 3;

  //define spreading parameter
  double lambda_ = 3 - n_aug_;

  //set vector for weights
  //VectorXd weights_ = VectorXd(2*n_aug_+1);
  double weight_0_ = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0_;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight_ = 0.5/(n_aug_+lambda_);
    weights_(i) = weight_;
  }
/*
  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_ <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create example vector for predicted state mean
  VectorXd x_ = VectorXd(n_x_);
  x_ <<
     5.93637,
     1.49035,
     2.20528,
    0.536853,
    0.353577;

  //create example matrix for predicted state covariance
  MatrixXd P_ = MatrixXd(n_x_,n_x_);
  P_ <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
 -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
 -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;

  //create example matrix with sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig_ <<
      6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
     0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
      2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;

  //create example vector for mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_ <<
      6.12155,
     0.245993,
      2.10313;

  //create example matrix for predicted measurement covariance
  MatrixXd S_ = MatrixXd(n_z_,n_z_);
  S_ <<
      0.0946171, -0.000139448,   0.00407016,
   -0.000139448,  0.000617548, -0.000770652,
     0.00407016, -0.000770652,    0.0180917;
*/
  //create example vector for incoming radar measurement
  VectorXd z_ = VectorXd(n_z_);
/*
  z_ <<
      5.9214,
      0.2187,
      2.0062;
*/  
  z_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);

/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //calculate cross correlation matrix
  Tc_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;
    //angle normalization
    while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    Tc_ = Tc_ + weights_(i) * x_diff_ * z_diff_.transpose();
  }

  //Kalman gain K;
  MatrixXd K_ = Tc_ * S_.inverse();

  //residual
  VectorXd z_diff_ = z_ - z_pred_;

  //angle normalization
  while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
  while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K_ * z_diff_;
  P_ = P_ - K_*S_*K_.transpose();

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  std::cout << "UpdateRadar: Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "UpdateRadar: Updated state covariance P: " << std::endl << P_ << std::endl;

  //write result
  //*x_out = x;
  //*P_out = P;

  /* expected result x:
     x =
 5.92276
 1.41823
 2.15593
0.489274
0.321338
    */

  /* expected result P:
     P =
  0.00361579 -0.000357881   0.00208316 -0.000937196  -0.00071727
-0.000357881   0.00539867   0.00156846   0.00455342   0.00358885
  0.00208316   0.00156846   0.00410651   0.00160333   0.00171811
-0.000937196   0.00455342   0.00160333   0.00652634   0.00669436
 -0.00071719   0.00358884   0.00171811   0.00669426   0.00881797
 
 meas_package.raw_measurements_ 1.01489
0.554329
 4.89281
 
UpdateRadar: Updated state x: 
 5.92276
 1.41823
 2.15593
0.489274
0.321338

UpdateRadar: Updated state covariance P: 
  0.00361579 -0.000357881   0.00208316 -0.000937196  -0.00071727
-0.000357881   0.00539867   0.00156846   0.00455342   0.00358885
  0.00208316   0.00156846   0.00410651   0.00160333   0.00171811
-0.000937196   0.00455342   0.00160333   0.00652634   0.00669436
-0.00071719   0.00358884   0.00171811   0.00669426   0.00881797
 
 */  

}

//#############################################################################
//### added functions 
//#############################################################################

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  //set state dimension
  n_x_ = 5;

  //define spreading parameter
  double lambda_ = 3 - n_x_;

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)      = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }

  //write result
  *Xsig_out = Xsig;

/* expected result:
   Xsig =
    5.7441  5.85768   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441
      1.38  1.34566  1.52806     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38
    2.2049  2.28414  2.24557  2.29582   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049
    0.5015  0.44339 0.631886 0.516923 0.595227   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015
    0.3528 0.299973 0.462123 0.376339  0.48417 0.418721 0.405627 0.243477 0.329261  0.22143 0.286879
*/
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
/*
  //set state dimension
  int n_x_ = 5;

  //set augmented dimension
  int n_aug_ = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_ = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_ = 0.2;

  //define spreading parameter
  double lambda_ = 3 - n_aug_;
  
  //create augmented mean vector
  VectorXd x_aug_ = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug_ = MatrixXd(7, 7);

  //create sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
*/
  //cout << "AugmentedSigmaPoints: lambda_ " << 3 - n_aug_ << endl;
  
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  //cout << "std_a_, std_yawdd_   " << std_a_ << "   " << std_yawdd_ << endl;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;

  //cout << "P_aug_ " << P_aug_ << endl;
  
  //create square root matrix
  MatrixXd L_ = P_aug_.llt().matrixL();

  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug_;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)        = x_aug_ + sqrt(lambda_+n_aug_) * L_.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug_ - sqrt(lambda_+n_aug_) * L_.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug_;

/* expected result:
   Xsig_aug =
  5.7441  5.85768   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441
    1.38  1.34566  1.52806     1.38     1.38     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38     1.38     1.38
  2.2049  2.28414  2.24557  2.29582   2.2049   2.2049   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049   2.2049   2.2049
  0.5015  0.44339 0.631886 0.516923 0.595227   0.5015   0.5015   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015   0.5015   0.5015
  0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528   0.3528 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528   0.3528
       0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641        0
       0        0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641

*/

}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out3) {
/*
  //set state dimension
  int n_x_ = 5;

  //set augmented dimension
  int n_aug_ = 7;

  //create example sigma point matrix
  MatrixXd Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
     Xsig_aug_ <<
    5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
      1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
    2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
    0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
    0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
         0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
         0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;
*/
  //create matrix with predicted sigma points as columns
  //Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double delta_t = dt; //0.1; //time diff in sec
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out3 = Xsig_pred_;

}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
/*
  //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
*/
/*
  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);
*/

/*******************************************************************************
 * Student part begin
 ******************************************************************************/
/*
  // set weights
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_);
    weights(i) = weight;
  }
*/
  //predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_ = x_+ weights_(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff_ * x_diff_.transpose() ;
  }


/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  //std::cout << "Predicted state" << std::endl;
  //std::cout << x << std::endl;
  //std::cout << "Predicted covariance matrix" << std::endl;
  //std::cout << P << std::endl;

  //write result
  *x_out = x_;
  *P_out = P_;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
/*
  //set state dimension
  n_x_ = 5;

  //set augmented dimension
  n_aug_ = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;
*/
  //define spreading parameter
  lambda_ = 3 - n_aug_;

  //set vector for weights
  VectorXd weights_ = VectorXd(2*n_aug_+1);
   double weight_0_ = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0_;
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight_ = 0.5/(n_aug_+lambda_);
    weights_(i) = weight_;
  }
/*
  //radar measurement noise standard deviation radius in m
  double std_radr_ = 0.3;

  //radar measurement noise standard deviation angle in rad
  double std_radphi_ = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  double std_radrd_ = 0.1;
*/
  

  //create example matrix with predicted sigma points
/*
  MatrixXd Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_ <<
         5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
           1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
          2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
         0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
          0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
*/
  //create matrix for sigma points in measurement space
  //MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);
//##rbx###
  
/*******************************************************************************
 * Student part begin
 ******************************************************************************/

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig_(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig_(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig_(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  //VectorXd z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  //MatrixXd S_ = MatrixXd(n_z_,n_z_);
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff_ = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

    S_ = S_ + weights_(i) * z_diff_ * z_diff_.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R_ = MatrixXd(n_z_,n_z_);
  R_ <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S_ = S_ + R_;

  
/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  //std::cout << "z_pred: " << std::endl << z_pred_ << std::endl;
  //std::cout << "S: " << std::endl << S_ << std::endl;

  //write result
  *z_out = z_pred_;
  *S_out = S_;
  
/*
z_pred: 
 6.12155
0.245993
 2.10313
S: 
   0.0946171 -0.000139448   0.00407016
-0.000139448  0.000617548 -0.000770652
0.00407016 -0.000770652    0.0180917

*/

}
