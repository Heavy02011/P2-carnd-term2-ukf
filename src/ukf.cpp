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
    //cout << "EKF: " << endl;
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

  cout << endl;
  cout << "########## Prediction ##########" << endl;
  cout << endl;
      
  float noise_ax = 9;
  float noise_ay = 9;

  //compute the time elapsed between the current and previous measurements
  dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = meas_package.timestamp_;
  
  // perform prediction step
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
      cout << endl;
      cout << "########## Updating Radar ##########" << endl;
      cout << endl;
      //cout << "meas_package.raw_measurements_" << meas_package.raw_measurements_ << endl;
      //UpdateEKF(meas_package.raw_measurements_);
      UpdateRadar(meas_package);

    }


  } else {

    if (use_laser_) {
      // Laser updates
      //measurement update
      cout << endl;
      cout << "########## Updating Lidar ##########" << endl;
      cout << endl;
      UpdateLidar(meas_package);
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
  //std::cout << "UKF::Prediction: Xsig_aug_ = " << std::endl << Xsig_aug_ << std::endl; //rbxok
  
  // 2. predict Sigma Points of next time step x(k+1|k)
  MatrixXd Xsig_pred = MatrixXd(15, 5);
  UKF::SigmaPointPrediction(&Xsig_pred);
  //std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
  //std::cout << "UKF::Prediction: Xsig_pred_ = " << std::endl << Xsig_pred_ << std::endl; // rbxok
  //cout << "UKF::Prediction: dt = " << std::endl << dt << std::endl; // rbxok
  
  // 3. predict mean and covariance of next time step
  VectorXd x_pred = VectorXd(5);
  MatrixXd P_pred = MatrixXd(5, 5);
  UKF::PredictMeanAndCovariance(&x_pred, &P_pred);
  //##std::cout << "x_pred = " << std::endl << x_pred << std::endl;
  //##std::cout << "P_pred = " << std::endl << P_pred << std::endl;
  //cout << "UKF::Prediction: x_ = " << std::endl << x_ << std::endl;
  //cout << "UKF::Prediction: P_ = " << std::endl << P_ << std::endl;
  //cout << "predict:<< std_a_, std_yawdd_   " << std_a_ << "   " << std_yawdd_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO: angles correct???

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // only 2 measurements, i.e. px & py
  //n_z_ = 2;
  
  // hint: https://discussions.udacity.com/t/lidar-prediction-and-update/243853/2
  //cout << endl;
  //cout << "UKF::UpdateLidar: Xsig_pred_ = " << Xsig_pred_ << endl;
  
  MatrixXd Zsig = Xsig_pred_.topRows(2);  
  //cout << "UKF::UpdateLidar: Zsig = " << Zsig << endl;

  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 2;
  
  //create example vector for incoming radar measurement
  VectorXd z_ = VectorXd(n_z_);
  z_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1);
  
  // Mean predicted measurement
  z_pred_  = Zsig * weights_;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);
   
  //calculate cross correlation matrix
  Tc_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    
    //residual
    VectorXd z_diff_ = Zsig.col(i) - z_pred_;  // Zsig_ --> Zsig from above !!
    
    //angle normalization
    //while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    //while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;
    
    // state difference
    VectorXd x_diff_ = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    //while (x_diff_(3)> M_PI) x_diff_(3)-=2.*M_PI;
    //while (x_diff_(3)<-M_PI) x_diff_(3)+=2.*M_PI;

    Tc_ = Tc_ + weights_(i) * x_diff_ * z_diff_.transpose();
  }
  cout << "Tc=" << Tc_ << endl;
  
  //measurement covariance matrix S (taken from radar measurement prediction)
  S_ = MatrixXd(n_z_,n_z_);
  S_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff_ = Zsig.col(i) - z_pred_;  // Zsig_ --> Zsig for lidar

    //angle normalization
    //while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
    //while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

    S_ = S_ + weights_(i) * z_diff_ * z_diff_.transpose();
  } 
  // add measurement noise for lidar
  S_ = S_ + R_lidar_;
  cout << "S_=" << S_ << endl;  
  
  //Kalman gain K;
  MatrixXd K_ = Tc_ * S_.inverse();
  
  //residual
  VectorXd z_diff_ = z_ - z_pred_;

  //angle normalization
  //while (z_diff_(1)> M_PI) z_diff_(1)-=2.*M_PI;
  //while (z_diff_(1)<-M_PI) z_diff_(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K_ * z_diff_;
  P_ = P_ - K_*S_*K_.transpose();
  
  //cout << "########## UKF::UpdateLidar +++ END +++ ##########" << endl;
    
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  
  // predict new z
  VectorXd z_out = VectorXd(3);
  MatrixXd S_out = MatrixXd(3, 3);
  PredictRadarMeasurement(&z_out, &S_out);
  
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
  n_z_ = 3;

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
  
  //create example vector for incoming radar measurement
  VectorXd z_ = VectorXd(n_z_); 
  z_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), meas_package.raw_measurements_(2);
  
  //create matrix for cross correlation Tc
  MatrixXd Tc_ = MatrixXd(n_x_, n_z_);

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

  //print result
  //std::cout << "UpdateRadar: Updated state x: " << std::endl << x_ << std::endl;
  //std::cout << "UpdateRadar: Updated state covariance P: " << std::endl << P_ << std::endl;

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

}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
    
  //create augmented mean state
  x_aug_.head(5) = x_;
  x_aug_(5) = 0;
  x_aug_(6) = 0;

  //create augmented covariance matrix
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(5,5) = P_;
  P_aug_(5,5) = std_a_*std_a_;
  P_aug_(6,6) = std_yawdd_*std_yawdd_;
  
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

}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out3) {
  
  //create matrix with predicted sigma points as columns
  //Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  double delta_t = dt; //0.1; //time diff in sec

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

  //write result
  *Xsig_out3 = Xsig_pred_;

}

void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
  
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

  //write result
  *x_out = x_;
  *P_out = P_;
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  n_z_ = 3;

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
  z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  //measurement covariance matrix S
  //MatrixXd S_ = MatrixXd(n_z_,n_z_);
  S_ = MatrixXd(n_z_,n_z_);
  
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

  //write result
  *z_out = z_pred_;
  *S_out = S_;

}
