#include "FusionUKF.h"
#include "utility.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// The fusion EKF engine.  
// This object maintains the state of the target, and can update it given a nwew measurement
FusionUKF::FusionUKF() {
  // Has the state been initialized yet?
  is_initialized_ = false;

  // Keep track of teh time of the most recent measurement update
  previous_timestamp_ = 0;

  // Fixed measurement noise covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  // Fixed measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  // Fixed linear measurement matrix for laser
  H_laser_ = MatrixXd::Identity(2, 4);

  // Set the process noise uncertainty
  s2ax_ = 9.0;
  s2ay_ = 9.0;

}


// Given a new measurement (either radar or lidar), process the measurement and update the state
void FusionUKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  // Get the timestamp of the measurememnt
  const long long timestamp = measurement_pack.timestamp_;
  
  // If this is the first measurement, we have to initialize the state
  if (!is_initialized_) {

    x_ = VectorXd::Zero(4);
    P_ = MatrixXd::Identity(4,4);
    P_(2,2) = 100.0;
    P_(3,3) = 100.0;  // we assume that we are more uncertain about velocity

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Initialize using a radar measurement
      const double r = measurement_pack.raw_measurements_[0];
      const double phi = measurement_pack.raw_measurements_[1];
      const double rdot = measurement_pack.raw_measurements_[2];
      x_[0] = r*cos(phi);
      x_[1] = r*sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize using a laser measuremement
      x_.segment(0,2) = measurement_pack.raw_measurements_;
    }
    else {
      std::cout << "UNRECOGNIZED measuremement type : " << measurement_pack.sensor_type_ << std::endl;
    }

    // That's all folks.  Record the time and return.
    is_initialized_ = true;
    previous_timestamp_ = timestamp;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  // How much time elapsed since the last update in seconds
  const double dt_s = double(timestamp - previous_timestamp_) / 1.0e6;
  // Make the linear prediction matrix
  MatrixXd F = MatrixXd::Identity(4,4);
  F(0,2) = dt_s;
  F(1,3) = dt_s; 
  // And the process noise
  MatrixXd Q = utility::makeProcessNoiseCovariance(dt_s,s2ax_,s2ay_);

  // Call the predict function to push x_ and P_ forward to the current time
  unscented_kalman_filter::Predict(x_,P_,F,Q);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // unpack the measureement data
  const VectorXd z = measurement_pack.raw_measurements_;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // process a radar measurement (have to compute the Jacobian for nonlinear measurement function)
    const MatrixXd H_radar = utility::CalculateJacobian(x_); 
    unscented_kalman_filter::Update(x_,P_,H_radar,R_radar_,utility::zminushx_radar(z,x_));
  } 
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) { 
    // process a laser measurement (measurement function is linear)
    unscented_kalman_filter::Update(x_,P_,H_laser_,R_laser_,utility::zminushx_laser(z,x_));
  }

  // udpate the clock
  previous_timestamp_ = timestamp;

}
