#pragma once

#include "FusionKF.h"
#include "measurement_package.h"
#include "Eigen/Dense"
#include "unscented_kalman_filter.h"

// The Unscented Kalman Filter object that stores and updates the state

class FusionUKF : public FusionKF {
public:

  FusionUKF();

  // update the state with the given measurement
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  double get_x() const;
  double get_y() const;
  double get_vx() const;
  double get_vy() const;

private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;
  
  Eigen::VectorXd process_noise_std_;

  Eigen::VectorXd radar_meas_noise_std_;
  Eigen::VectorXd laser_meas_noise_std_;
};

