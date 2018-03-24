#pragma once

#include "measurement_package.h"
#include "Eigen/Dense"

// Fusion KF defines the abstract class and interface for all kalman filters.

class FusionKF {
public:

  // update the internal state with the given measurement
  virtual void ProcessMeasurement(const MeasurementPackage &measurement_pack) = 0;

  // retrieve the current state
  virtual double get_x() const = 0;
  virtual double get_y() const = 0;
  virtual double get_vx() const = 0;
  virtual double get_vy() const = 0;

  // state (mean and covariance)
  Eigen::VectorXd x_;
  Eigen::MatrixXd P_;

};

