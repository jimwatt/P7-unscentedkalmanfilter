#include "Eigen/Dense"

// A place for pure functions to live undisturbed from the silent evil of side-effects and global variables

namespace utility {

// compute the Jacobian fo the Radar measurement function
Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);

// Function to wrap an angle from -pi to pi
double ModPosNegPi(const double xraw) ;

// Create the process noise covariance for elapsed time and acceleration noise
Eigen::MatrixXd makeProcessNoiseCovariance(const double dt, const double s2ax, const double s2ay);

// Measurement residul for radar
Eigen::VectorXd zminushx_radar(const Eigen::VectorXd& z, const Eigen::VectorXd& x);

// Measurement residual for laser
Eigen::VectorXd zminushx_laser(const Eigen::VectorXd& z, const Eigen::VectorXd& x);


}   //namespace utility