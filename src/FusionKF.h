#pragma once

#include "measurement_package.h"
#include "Eigen/Dense"

class FusionKF {
public:

  // update the state with the given measurement
  virtual void ProcessMeasurement(const MeasurementPackage &measurement_pack) = 0;

  // state (mean and covariance)
  Eigen::VectorXd x_;
  Eigen::MatrixXd P_;

};

