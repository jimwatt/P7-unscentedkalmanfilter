#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

// Compute the RMSE error between estimate and truth
VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

	VectorXd rmse = VectorXd::Zero(4);

	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
		cout << "Invalid estimation or ground_truth data" << endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){
		VectorXd residual = estimations[i] - ground_truth[i];

		while (residual[3]> M_PI) residual[3]-=2.*M_PI;
  		while (residual[3]<-M_PI) residual[3]+=2.*M_PI;

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	return rmse;

}


