#include "FusionUKF.h"
#include "unscented_kalman_filter.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

VectorXd radar_measurement_function(const VectorXd& x) {
    const double p_x = x[0];
    const double p_y = x[1];
    const double v  = x[2];
    const double yaw = x[3];

    const double v1 = cos(yaw)*v;
    const double v2 = sin(yaw)*v;

    VectorXd z = VectorXd::Zero(3);
    z[0] = sqrt(p_x*p_x + p_y*p_y);                        //r
    z[1] = atan2(p_y,p_x);                                 //phi
    z[2] = (p_x*v1 + p_y*v2 ) / z[0];   //r_dot

    return z;

}

VectorXd laser_measurement_function(const VectorXd& x) {
  return x.segment(0,2);
}


// The fusion EKF engine.  
// This object maintains the state of the target, and can update it given a nwew measurement
FusionUKF::FusionUKF() {
  // Has the state been initialized yet?
  is_initialized_ = false;

  // Keep track of the time of the most recent measurement update
  previous_timestamp_ = 0;

  // process noise 
  process_noise_std_ = VectorXd::Zero(2);
  process_noise_std_[0] = 1.5;
  process_noise_std_[1] = 0.6;

  radar_meas_noise_std_ = VectorXd::Zero(3);
  radar_meas_noise_std_[0] = 0.3;
  radar_meas_noise_std_[1] = 0.03;
  radar_meas_noise_std_[2] = 0.3;

  laser_meas_noise_std_ = VectorXd::Zero(2);
  laser_meas_noise_std_[0] = 0.15;
  laser_meas_noise_std_[1] = 0.15;

}


// Given a new measurement (either radar or lidar), process the measurement and update the state
void FusionUKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

  // Get the timestamp of the measurememnt
  const long long timestamp = measurement_pack.timestamp_;
  
  // If this is the first measurement, we have to initialize the state
  if (!is_initialized_) {

    x_ = VectorXd::Zero(5);
    P_ = MatrixXd::Identity(5,5);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Initialize using a radar measurement
      const double r = measurement_pack.raw_measurements_[0];
      const double phi = measurement_pack.raw_measurements_[1];
      x_[0] = r*cos(phi);
      x_[1] = r*sin(phi);

      MatrixXd JacobianH = MatrixXd::Zero(2,2);
      JacobianH(0,0) = cos(phi);
      JacobianH(0,1) = r*sin(phi);
      JacobianH(1,0) = sin(phi);
      JacobianH(1,1) = -r*cos(phi);

      MatrixXd R = MatrixXd::Zero(2,2);
      R(0,0) =  radar_meas_noise_std_[0]*radar_meas_noise_std_[0];
      R(1,1) = radar_meas_noise_std_[1]*radar_meas_noise_std_[1];

      // map the measuremement covariance to the state space
      P_.block(0,0,2,2) = JacobianH * R * JacobianH.transpose();

      P_(2,2) = 10.0;
      P_(3,3) = 5.0;
      P_(4,4) = 2.0;  // we assume that we are more uncertain about velocity
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize using a laser measuremement
      x_.segment(0,2) = measurement_pack.raw_measurements_;

      P_(0,0) = laser_meas_noise_std_[0] * laser_meas_noise_std_[0];
      P_(1,1) = laser_meas_noise_std_[1] * laser_meas_noise_std_[1];
      P_(2,2) = 10.0;
      P_(3,3) = 5.0;
      P_(4,4) = 2.0;  // we assume that we are more uncertain about velocity 
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

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  // unpack the measureement data
  const VectorXd z = measurement_pack.raw_measurements_;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // process a radar measurement 
    unscented_kalman_filter::UpdateState(dt_s,process_noise_std_, radar_meas_noise_std_, 
      z, radar_measurement_function, 
      x_, P_);

  } 
  else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) { 
    // process a laser measurement
    unscented_kalman_filter::UpdateState(dt_s,process_noise_std_, laser_meas_noise_std_, 
      z, laser_measurement_function, 
      x_, P_);
  }

  // udpate the clock
  previous_timestamp_ = timestamp;

}

  double FusionUKF::get_x() const {
    return x_[0];
  }
  double FusionUKF::get_y() const {
    return x_[1];
  }
  double FusionUKF::get_vx() const {
    return x_[2] * cos(x_[3]);
  }
  double FusionUKF::get_vy() const {
    return x_[2] * sin(x_[3]);
  }
