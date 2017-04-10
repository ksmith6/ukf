#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector [px, py, v, yaw, yawRate]
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
  Hint: one or more values initialized above might be wildly off...
  */
  
  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  VectorXd weights_; //  = VectorXd(2 * n_aug_ + 1);

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  ///* the current NIS for radar
  NIS_radar_ = 0.0;

  ///* the current NIS for laser
  NIS_laser_ = 0.0;

  DebugMode_ = true;

  // Added variables

  // Container for UKF sigma points.
  MatrixXd Xsig_; 

  MatrixXd Xsig_pred_;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (is_initialized_) {
    // Yes, the filter is initialized, so process the measurement using
    // the typical predict-correct process.

    // Compute the elapsed time from the previous measurement/state.
    double dt = meas_package.timestamp_ - previous_timestamp_; 

    // First, update the state and covariance to the current timestamp of the received measurement.
    Prediction(dt);

    // Second, perform the measurement update based on the measurement type.
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
      // Radar measurment update mode.
      UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
      // Laser Measurement Update mode.
      UpdateLidar(meas_package);
    } else {
      // Error-handling code for unknown measurement types.
      // Perform no update.
      cout << "Received unknown measurement type" << endl;
    }
  } else {
    // Initialization logic.
    Initialization(meas_package);
  }

  

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
  bool DebugThis = true;

  if (DebugMode_) {
    cout << " ------- PREDICTION ------- " << endl;
  }

  /**
  1) Generate sigma points
  2) Predict sigma points
  3) Predict mean and covariance
  */

  // ======== 1) Generate Sigma Points ==============
  // Generate matrix to hold sigma points.
  if (DebugMode_ && DebugThis) { cout << "  Inititalizing Xsig_ " << endl; }
  Xsig_ = MatrixXd(n_aug_, 2*n_aug_+1);
  if (DebugMode_ && DebugThis) { cout << "    Xsig_.size()  : (" << Xsig_.rows() << ", " << Xsig_.cols() << ")" << endl; }

  
  // Generate augmented sigma points
  AugmentedSigmaPoints();

  // ======== 2) Predict Sigma Points ===============
  // MatrixXd Xsig_pred = MatrixXd(n_aug, )
  // MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  //SigmaPointPrediction(&Xsig_pred);
  SigmaPointPrediction(delta_t);

  // ======== 3) Predict mean and covariance ========
  PredictMeanAndCovariance();



}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  // Create matrix to hold sigma points
  MatrixXd Xsig = MatrixXd(n_x_, 2*n_x_ + 1);

  // Compute square root of covariance
  MatrixXd A = P_.llt().matrixL();

  // Store the mean in the first column.
  Xsig.col(0) = x_;
  
  // Compute the spreading parameter
  float sig = sqrt(lambda_ + n_x_);

  for (int i=0; i<n_x_; i++) {
      Xsig.col(i+1)      = x_ + sig * A.col(i);
      Xsig.col(i+1+n_x_) = x_ - sig * A.col(i);
  }

  // Write output matrix to memory.
  *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints() {
  if (DebugMode_) {
    cout << "  a) Generating Augmented Sigma Points" << endl;
  }

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
 
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_ * std_a_;
  P_aug(6,6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  float sig = sqrt(lambda_ + n_aug_);
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sig * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sig * L.col(i);
  }

  //write result
  Xsig_ = Xsig_aug;
}

/**
  * This predicts the sigma points forward/backward some delta_t. 
  * This algorithm has a test mode, and testing indicates that it is performing nominally.
  */
void UKF::SigmaPointPrediction(double delta_t) {
  // MatrixXd* Xsig_out, 
  // MatrixXd* Xsig_pred, 
  bool DebugThis = false;

  bool TestMode = false;
  if (TestMode) {
    delta_t = 0.1;
    Xsig_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
     Xsig_ <<
    5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
      1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
    2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
    0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
    0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
         0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
         0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;
  }



  if (DebugMode_) {
    cout << "  b) Predicting Sigma Points Forward" << endl;
  }

  // re-initialize matrix with predicted sigma points as columns
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  //predict sigma points
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_(0,i);
    double p_y = Xsig_(1,i);
    double v = Xsig_(2,i);
    double yaw = Xsig_(3,i);
    double yawd = Xsig_(4,i);
    double nu_a = Xsig_(5,i);
    double nu_yawdd = Xsig_(6,i);

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
  if (DebugMode_ && DebugThis) { cout << "    Xsig_.size()  : (" << Xsig_.rows() << ", " << Xsig_.cols() << ")" << endl; }
  if (DebugMode_ && DebugThis) { cout << "    Xsig_ : " << endl << Xsig_ << endl; }

  if (TestMode) {
    cout << "Xsig_pred_" <<  Xsig_pred_ << endl;
  }

}

/**
  * Computes the mean and covariance of the predicted sigma points.
  * Probably has a bug in it.
  */
void UKF::PredictMeanAndCovariance() {

  bool DebugThis = true;

  if (DebugMode_) {
    cout << "  c) Predicting mean and covariance" << endl;
  }

  //create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x_);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x_, n_x_);

  if (DebugMode_) { cout << "    1) Populating weights" << endl; }
  
  // set weights
  double weight_0 = lambda_/(lambda_ + n_aug_);
  double weight = 0.5/(n_aug_+lambda_);

  if (DebugMode_ && DebugThis) { cout << "       weight_0 = " << weight_0 << endl; }
  if (DebugMode_ && DebugThis) { cout << "       weight = " << weight << endl; }
  if (DebugMode_ && DebugThis) { cout << "       weights_.size() = " << weights_.size() << endl; }  


  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_ + 1; i++) {
    weights_(i) = weight;
  }
  

  if (DebugMode_ && DebugThis) { cout << "       Filled weights_ with weight " << endl; }
  if (DebugMode_ && DebugThis) { cout << "       weights_: " << endl << weights_ << endl; }

  //predicted state mean
  x.fill(0.0);
  if (DebugMode_ && DebugThis) { cout << "    Xsig_.size()  : (" << Xsig_.rows() << ", " << Xsig_.cols() << ")" << endl; }

  if (DebugMode_ && DebugThis) { cout << Xsig_.col(0) << endl; }

  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    double w = weights_(i);
    x = x + weights_(i) * Xsig_pred_.col(i);
  }
  if (DebugMode_ && DebugThis) { cout << "      Mean State: " << endl << x << endl;}

  //predicted state covariance matrix
  P.fill(0.0);
  if (DebugMode_ && DebugThis) { cout << "       Filled covariance P with zeros " << endl; }
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x;
    
    if (DebugMode_) { cout << "x_diff = " << x_diff << endl; }

    //angle normalization
    x_diff(3) = angleNormalization(x_diff(3));

    P = P + weights_(i) * x_diff * x_diff.transpose() ;
  }
  if (DebugMode_ && DebugThis) { cout << "      Covariance: " << endl << P << endl; }
  if (DebugMode_ && DebugThis) { cout << "       Writing output x_ and P_ " << endl; }
  // write result
  x_ = x;
  P_ = P;
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
  if (DebugMode_) {
    cout << " ------- LASER MEAS UPDATE ------- " << endl;
  }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  if (DebugMode_) {
    cout << " ------- RADAR MEAS UPDATE ------- " << endl;
  }

  int n_z = 3;
  VectorXd z_out = VectorXd(n_z);
  MatrixXd S_out = MatrixXd(n_z, n_z);
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred.fill(0.0); // TODO - Change code to import Xsig_pred from Prediction() method
  PredictRadarMeasurement(&z_out, &S_out, &Xsig_pred);

}

/**
 * Generates predicted radar measurements. 
 */
void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Xsig_pred) {

  if (DebugMode_) { cout << "  a) Predict Radar Measurement" << endl; }

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug_;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug_+1);
  
  // Initialize weights (non-looping)
  double weight = 0.5/(n_aug_+lambda);
  weights.fill(weight);
  weights(0) = lambda/(lambda+n_aug_);


  //radar measurement noise standard deviation radius in m
  //double std_radr = 0.3;

  //radar measurement noise standard deviation angle in rad
  //double std_radphi = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  //double std_radrd = 0.1;

  //create example matrix with predicted sigma points
  //MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  //Xsig_pred <<
  //       5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  //         1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  //        2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  //       0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  //        0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);


  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = (*Xsig_pred)(0,i);
    double p_y = (*Xsig_pred)(1,i);
    double v  = (*Xsig_pred)(2,i);
    double yaw = (*Xsig_pred)(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    // TODO - Fix this
    //while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    //while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // John Chen angle normalization logic - taken from Slack forum.
    if (z_diff(1)> M_PI) z_diff(1) = remainder(z_diff(1), (2.*M_PI)) - M_PI;
    if (z_diff(1)<-M_PI) z_diff(1) = remainder(z_diff(1), (2.*M_PI)) + M_PI;

    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R <<    std_radr_*std_radr_, 0, 0,
          0, std_radphi_*std_radphi_, 0,
          0, 0,std_radrd_*std_radrd_;
  S = S + R;

  
/*******************************************************************************
 * Student part end
 ******************************************************************************/

  //print result
  //std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  //std::cout << "S: " << std::endl << S << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
}

/**
  * Wraps the input angle into the [-pi pi] domain without resorting to loops.
  * @param {double} angle
  */
double UKF::angleNormalization(double angle) {
  if (angle > M_PI) {
    return remainder(angle, (2.*M_PI)) - M_PI;
  } else if (angle < -M_PI) {
    return remainder(angle, (2.*M_PI)) + M_PI;
  } else {
    return angle;
  }
}


/**
  * Initializes the filter state and covariance with some measurement.
  * @param {MeasurementPackage} meas_package
  */
void UKF::Initialization(MeasurementPackage meas_package) {
  if (DebugMode_) {
    cout << "Initializing filter with ";
  }
  // First, set the timestamp.
  previous_timestamp_ = meas_package.timestamp_;

  // Next, initialize the state and covariance according to the measurement type.
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    // Radar measurment update mode.
    InitFilterRadar(meas_package);
    if (DebugMode_) { 
      cout << "Radar measurement" << endl;
    }
    is_initialized_ = true;
  } else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    // Laser Measurement Update mode.
    InitFilterLaser(meas_package);
    if (DebugMode_) { 
      cout << "Laser measurement" << endl;
    }
    is_initialized_ = true;
  } else {
    // Error-handling code for unknown measurement types or all measurements disabled.
    // Perform no initialization.
    cout << "Received unknown measurement type" << endl;
    is_initialized_ = false;
  }
}

/**
  * Initializes the filter state and covariance with a radar measurement.
  * @param {MeasurementPackage} meas_package
  */
void UKF::InitFilterRadar(MeasurementPackage meas_package) {
  // TODO - Initialize the state and covariance based on a radar measurement

  /*
  Convert radar from polar to cartesian coordinates and initialize state.
  */
  float range = meas_package.raw_measurements_[0];
  float bearing = meas_package.raw_measurements_[1];
  float rangeRate = meas_package.raw_measurements_[2]; 

  // TODO - Convert polar measuremnts to Cartesian coordinates
  float px = range * cos(bearing);
  float py = range * sin(bearing);
  float vx = rangeRate * cos(bearing);
  float vy = rangeRate * sin(bearing);

  float yaw = atan2(vy,vx);
  float yawRate = 0.0;  // Assumes initial yaw rate is zero.

  // Initialize the state vector.
  x_ << px, py, sqrt(vx*vx + vy*vy), yaw, yawRate;

  // Initialize the covariance matrix with all zeros.
  P_.setZero(5,5);

  // Populate only the diagonal elements.
  P_(0,0) = std_laspx_ * std_laspx_;
  P_(1,1) = std_laspy_ * std_laspy_;
  P_(2,2) = 0; // TODO
  P_(3,3) = 0; // TODO
  P_(4,4) = 0; // TODO
}

/**
  * Initializes the filter state and covariance with a laser measurement.
  * @param {MeasurementPackage} meas_package
  */
void UKF::InitFilterLaser(MeasurementPackage meas_package) {
  /**
  Initialize states.
  NOTE: Initial velocity states initialized with the assumption that the object is stationary.
  */
  float px = meas_package.raw_measurements_[0];
  float py = meas_package.raw_measurements_[1];
  float vMag = 0.0;   // Assumes object is stationary
  float yaw = 0.0; // No yaw information directly observable from single laser measurement
  float yawRate = 0.0; // No yaw rate information directly observable from single laser measurement
  // Initialize the states with the first laser measurement.
  x_ << px, py, vMag, yaw, yawRate;

  // Initialize the covariance matrix with all zeros.
  P_.setZero(5,5);

  // Populate diagonals with initial values.
  P_(0,0) = std_laspx_ * std_laspx_;
  P_(1,1) = std_laspy_ * std_laspy_;
  P_(2,2) = 0; // TODO
  P_(3,3) = 0; // TODO
  P_(4,4) = 0; // TODO


}
