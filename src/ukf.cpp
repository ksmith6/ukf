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
  VectorXd weights_; // TODO

  ///* State dimension
  n_x_ = 5;

  ///* Augmented state dimension
  n_aug_ = n_x_ + 2;

  ///* Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  ///* the current NIS for radar
  NIS_radar_ = 0.0;

  ///* the current NIS for laser
  NIS_laser_ = 0.0;

  DebugMode_ = true;
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
  int n_aug = n_x_ + 2;
  MatrixXd Xsig = MatrixXd(n_aug, 2*n_x_+1);
  
  // Generate augmented sigma points
  AugmentedSigmaPoints(&Xsig);

  // ======== 2) Predict Sigma Points ===============
  // MatrixXd Xsig_pred = MatrixXd(n_aug, )
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug + 1);
  //SigmaPointPrediction(&Xsig_pred);
  SigmaPointPrediction(&Xsig, &Xsig_pred, delta_t);

  // ======== 3) Predict mean and covariance ========
  //PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred);



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

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //set augmented dimension
  int n_aug = n_x_ + 2;

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
 
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
  float sig = sqrt(lambda_ + n_aug);
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sig * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sig * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, MatrixXd* Xsig_pred, double delta_t) {

  //set augmented dimension
  int n_aug = n_x_ + 2;

  //create matrix with predicted sigma points as columns
  // MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug + 1);

  //double delta_t = 0.1; //time diff in sec

  //predict sigma points
  for (int i = 0; i< 2*n_aug+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_out(0,i); // <------- This chunk of code (230-236) is preventing compilation.
    double p_y = Xsig_out(1,i);
    double v = Xsig_out(2,i);
    double yaw = Xsig_out(3,i);
    double yawd = Xsig_out(4,i);
    double nu_a = Xsig_out(5,i);
    double nu_yawdd = Xsig_out(6,i);

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
    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  //write result
  //*Xsig_out = Xsig_pred;
  //return Xsig_pred;
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
