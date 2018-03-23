#include "utility.h"

#include <iostream>

namespace utility {

using Eigen::MatrixXd;
using Eigen::VectorXd;

// compute the Jacobian fo the Radar measurement function
MatrixXd CalculateJacobian(const VectorXd& x_state) {

  MatrixXd Hj = MatrixXd::Zero(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  //pre-compute a set of terms to avoid repeated calculation
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);

  //check division by zero
  if(fabs(c1) < 0.0001){
    std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
    return Hj;
  }

  //compute the Jacobian matrix
  Hj(0,0) =  px/c2;
  Hj(0,1) =  py/c2;
  Hj(1,0) = -py/c1;
  Hj(1,1) =  px/c1;
  Hj(2,0) =  py*(vx*py - vy*px)/c3;
  Hj(2,1) =  px*(px*vy - py*vx)/c3;
  Hj(2,2) =  px/c2;
  Hj(2,3) =  py/c2;

  return Hj;
}


// Function to wrap an angle from -pi to pi
double ModPosNegPi(const double xraw) {
    const double x = xraw + M_PI;
    const double y = 2.0*M_PI;
    const double m = x - y * floor(x/y);
    if (m>=y) 
        return 0;
    if (m<0) {
        if (y+m == y)
            return 0;
        else
            return y + m; 
    }
    return m - M_PI;
}

// Create the process noise covariance for elapsed time and acceleration noise
MatrixXd makeProcessNoiseCovariance(const double dt, const double s2ax, const double s2ay) {
  MatrixXd Q = MatrixXd::Zero(4,4);
  const double dt2 = dt*dt;
  const double dt3 = dt2*dt;
  const double dt4 = dt3*dt;
  Q(0,0) = 0.25*dt4*s2ax;
  Q(1,1) = 0.25*dt4*s2ay;
  Q(2,2) = dt2*s2ax;
  Q(3,3) = dt2*s2ay;
  Q(2,0) = Q(0,2) = 0.5*dt3*s2ax;
  Q(3,1) = Q(1,3) = 0.5*dt3*s2ay;
  return Q;
}

// Measurement residul for radar
VectorXd zminushx_radar(const VectorXd& z, const VectorXd& x) {
  const double r = x.segment(0,2).norm();
  const double phi = atan2(x[1],x[0]);
  const double rdot = (x[0]*x[2] + x[1]*x[3]) / (r + 1e-8);
  VectorXd zmhx = VectorXd::Zero(3);
  zmhx[0] = z[0] - r;
  zmhx[1] = ModPosNegPi(z[1] - phi);
  zmhx[2] = z[2] - rdot;
  return zmhx;
}

// Measurement residual for laser
VectorXd zminushx_laser(const VectorXd& z, const VectorXd& x) {
  return z - x.segment(0,2);
}

}    // namespace utility