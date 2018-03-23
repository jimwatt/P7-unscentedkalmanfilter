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

  // Process noise
  double s2ax_;
  double s2ay_;
};

