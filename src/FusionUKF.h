#pragma once

#include "FusionKF.h"
#include "measurement_package.h"
#include "Eigen/Dense"
#include "unscented_kalman_filter.h"

class FusionUKF : public FusionKF {
public:

  FusionUKF();

  // update the state with the given measurement
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  double get_x();
  double get_y();
  double get_vx();
  double get_vy();

private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

  // Measurement noise
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  
  // Laser measurement function
  Eigen::MatrixXd H_laser_;

  Eigen::VectorXd process_noise_std_;
  Eigen::VectorXd radar_meas_noise_std_;
  Eigen::VectorXd laser_meas_noise_std_;
};

