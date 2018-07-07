#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
   	VectorXd rmse(4);
	VectorXd xest(4), xtrue(4), xdiff(4), xdiff2(4);
	rmse << 0,0,0,0;

    int n = estimations.size();
	if(estimations.size() == 0 || estimations.size() != ground_truth.size())
    {
        std::cout << "Invalid estimation or ground truth" << std::endl;
	    return rmse;
    }

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
        // ... your code here
		xest = estimations[i];
		xtrue = ground_truth[i];
		xdiff = (xest - xtrue);
		xdiff2 = xdiff.array()*xdiff.array();
		//rmse = rmse +  xdiff.array() * xdiff.array();

		rmse += xdiff2;
	}

	//calculate the mean
	rmse = rmse *1.0/(float)n;

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}
