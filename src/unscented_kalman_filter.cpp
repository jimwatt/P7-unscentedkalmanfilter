#include "unscented_kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

namespace unscented_kalman_filter {


	void PredictMeanAndCovariance(const MatrixXd& Xsig_pred, const VectorXd& weights, const double lmbda, VectorXd& x, MatrixXd& P) {

  //set augmented dimension
		const int n_x = x.size();
		const int n_aug = (weights.size() - 1)/2;

  //predicted state mean
	x.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
  	x = x + weights[i] * Xsig_pred.col(i);
  }

  //predicted state covariance matrix
  P.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //iterate over sigma points
  	VectorXd x_diff = Xsig_pred.col(i) - x;
  	while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
  	while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
  	P = P + weights[i] * x_diff * x_diff.transpose() ;
  }

}


MatrixXd SigmaPointPrediction(const int n_x, const double delta_t, const MatrixXd& Xsig_aug) {

  //set augmented dimension
	const int n_sigmapts = Xsig_aug.cols();

  //create matrix with predicted sigma points as columns
	MatrixXd Xsig_pred = MatrixXd(n_x, n_sigmapts);

  //predict sigma points
	for (int i = 0; i< n_sigmapts; i++) {
    //extract values for better readability
		const double p_x = Xsig_aug(0,i);
		const double p_y = Xsig_aug(1,i);
		const double v = Xsig_aug(2,i);
		const double yaw = Xsig_aug(3,i);
		const double yawd = Xsig_aug(4,i);
		const double nu_a = Xsig_aug(5,i);
		const double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
		double px_p, py_p;

    //avoid division by zero
		if (fabs(yawd) > 0.001) {
			px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
			py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
		}
		else {
			px_p = p_x + v*delta_t*cos(yaw);
			py_p = p_y + v*delta_t*sin(yaw);
		}

    //write predicted sigma point into correct column
		Xsig_pred(0,i) = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
		Xsig_pred(1,i) = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
		Xsig_pred(2,i) = v + nu_a*delta_t;
		Xsig_pred(3,i) = yaw + yawd*delta_t + 0.5*nu_yawdd*delta_t*delta_t;
		Xsig_pred(4,i) = yawd + nu_yawdd*delta_t;
	}

	return Xsig_pred;

}


MatrixXd AugmentedSigmaPoints(const VectorXd& process_noise_std, const VectorXd& x, const MatrixXd& P, const double lambda) {

  //set state dimension
	const int n_x = x.size();
	const int n_noise = process_noise_std.size();
	const int n_aug = n_x + n_noise;

  //create augmented mean vector
	VectorXd x_aug = VectorXd::Zero(n_aug);
	x_aug.segment(0,n_x) = x;

  //create augmented state covariance
	MatrixXd P_aug = MatrixXd::Zero(n_aug, n_aug);
	P_aug.topLeftCorner(n_x,n_x) = P;
	for(int nn=0;nn<n_noise;++nn) {
		P_aug(n_x+nn,n_x+nn) = process_noise_std[nn] * process_noise_std[nn];
	}
	const MatrixXd L = P_aug.llt().matrixL();

  //create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd::Zero(n_aug, 2 * n_aug + 1);
	Xsig_aug.col(0)  = x_aug;
	for (int i = 0; i<n_aug; i++) {
		Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
		Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
	}

	return Xsig_aug;

}



void PredictMeasurement(const VectorXd& meas_noise_std, const MatrixXd& Xsig_pred, const VectorXd& weights, const double lambda,
	VectorXd measurement_function(const VectorXd& x),
	VectorXd& z_out, MatrixXd& S_out, MatrixXd& Zsig) {

  //set state dimension
	const int n_x = Xsig_pred.rows();

  //set augmented dimensionconst VectorXd& x
	const int n_aug = (Xsig_pred.cols()-1)/2;

  //set measurement dimension
	const int n_z = meas_noise_std.size();

	Zsig = MatrixXd::Zero(n_z, 2 * n_aug + 1);

	  //transform sigma points into measurement space
	for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
	  	Zsig.col(i) = measurement_function(Xsig_pred.col(i));
	}

	  //mean predicted measurement
	VectorXd z_pred = VectorXd::Zero(n_z);
	for (int i=0; i < 2*n_aug+1; i++) {
		z_pred = z_pred + weights(i) * Zsig.col(i);
	}

	  //innovation covariance matrix S
	MatrixXd S = MatrixXd::Zero(n_z,n_z);
	for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
		VectorXd z_diff = Zsig.col(i) - z_pred;
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
		S = S + weights(i) * z_diff * z_diff.transpose();
	}

	//add measurement noise covariance matrix
	MatrixXd R = MatrixXd::Zero(n_z,n_z);
	for(int ii=0;ii<n_z;++ii) {
		R(ii,ii) = meas_noise_std[ii]*meas_noise_std[ii];
	}
	S = S + R;

	//write result
	z_out = z_pred;
	S_out = S;


}


void ApplyKalmanGain(const MatrixXd& Xsig_pred, const VectorXd& z, const MatrixXd& Zsig, const VectorXd& z_pred, const MatrixXd& S, const VectorXd& weights, VectorXd& x, MatrixXd& P) {

	const int n_x = x.size();
	const int n_z = Zsig.rows();
	const int n_aug = (Zsig.cols()-1)/2;

	MatrixXd Tc = MatrixXd::Zero(n_x, n_z);

	for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
	//residual
		VectorXd z_diff = Zsig.col(i) - z_pred;
		while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
		while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	// state difference
		VectorXd x_diff = Xsig_pred.col(i) - x;
		while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
		while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
		Tc = Tc + weights(i) * x_diff * z_diff.transpose();
	}

	//Kalman gain K;
	const MatrixXd K = Tc * S.inverse();
	VectorXd z_diff = z - z_pred;
	while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
	while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

	//update state mean and covariance matrix
	x = x + K * z_diff;
	P = P - K*S*K.transpose();

}

VectorXd generateWeights(const int n, const double lambda) {
	VectorXd weights = VectorXd::Zero(2*n+1);
	weights(0) = lambda/(lambda+n);
	  for (int i=1; i<2*n+1; i++) {  //2n+1 weights
	  	weights(i) = 0.5/(n+lambda);
	  }
	  return weights;
}


void UpdateState(const double delta_t, 
	const VectorXd& process_noise_std, 
	const VectorXd& meas_noise_std, 
	const VectorXd& z, 
	VectorXd measurement_function(const VectorXd&),
	VectorXd& x, 
	MatrixXd& P) {

	//set state dimension
	const int n_x = x.size();
	const int n_noise = process_noise_std.size();
	const int n_aug = n_x + n_noise;
	const int n_z = z.size();

	//define spreading parameter
	const double lambda = 3 - n_aug;

	// Define Weights
	const VectorXd weights = generateWeights(n_aug,lambda);

	  // Generate sigma points
	  const MatrixXd Xsig_aug = AugmentedSigmaPoints(process_noise_std, x, P, lambda);

	  // Predict the sigma points
	  const MatrixXd Xsig_pred = SigmaPointPrediction(n_x, delta_t, Xsig_aug);

	  // Predict the mean and covariance
	  PredictMeanAndCovariance(Xsig_pred, weights, lambda, x, P);

	  // Predict the Measurement
	  MatrixXd Zsig = MatrixXd::Zero(n_z, 2 * n_aug + 1);
	  VectorXd z_pred = VectorXd::Zero(n_z);
	  MatrixXd S = MatrixXd::Zero(n_z,n_z);
	  PredictMeasurement(meas_noise_std, Xsig_pred, weights, lambda, measurement_function, z_pred, S, Zsig);

	  // Update the State
	  ApplyKalmanGain(Xsig_pred,z,Zsig,z_pred,S,weights,x,P);



}




}   // namespace unscented_kalman_filter


