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

  // initial state vector
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

  DebugMode_ = false;
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

    double dt = 0;  // TODO

    // First, update the state and covariance to the current timestamp of the received measurement.
    Prediction(dt);



    // Second, perform the measurement update based on the measurement type.
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar measurment update mode.
      UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      // Laser Measurement Update mode.
      UpdateLidar(meas_package);
    } else {
      cout << "Received unknown measurement type" << endl;
    }
  } else {
    // Initialization logic.

    // TODO

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
