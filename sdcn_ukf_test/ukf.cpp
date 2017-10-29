#include "stdafx.h"
#include <iostream>
#include "ukf.h"

UKF::UKF() {
	//TODO Auto-generated constructor stub
	Init();
}

UKF::~UKF() {
	//TODO Auto-generated destructor stub
}

void UKF::Init() {

}


/*******************************************************************************
* Programming assignment functions:
*******************************************************************************/

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

	//set state dimension
	int n_x = 5;

	//set augmented dimension
	int n_aug = 7;

	//Process noise standard deviation longitudinal acceleration in m/s^2
	double std_a = 0.2;

	//Process noise standard deviation yaw acceleration in rad/s^2
	double std_yawdd = 0.2;

	//define spreading parameter
	double lambda = 3 - n_aug;

	//set example state
	VectorXd x = VectorXd(n_x);
	x << 5.7441,
		1.3800,
		2.2049,
		0.5015,
		0.3528;

	//create example covariance matrix
	MatrixXd P = MatrixXd(n_x, n_x);
	P << 0.0043, -0.0013, 0.0030, -0.0022, -0.0020,
		-0.0013, 0.0077, 0.0011, 0.0071, 0.0060,
		0.0030, 0.0011, 0.0054, 0.0007, 0.0008,
		-0.0022, 0.0071, 0.0007, 0.0098, 0.0100,
		-0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

	//create augmented mean vector
	VectorXd x_aug = VectorXd(n_aug);

	//create augmented state covariance
	MatrixXd P_aug = MatrixXd(n_aug, n_aug);
	P_aug.setZero();

	//create sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);

	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//create augmented mean state
	x_aug.setZero();
	x_aug.head(5) = x;

	//create augmented covariance matrix
	P_aug.topLeftCorner(n_x, n_x) = P;
	MatrixXd Q = MatrixXd(2, 2);
	Q.setZero();
	const float std_a_2 = std_a * std_a;
	const float std_yawdd_2 = std_yawdd * std_yawdd;
	Q(0, 0) = std_a_2;
	Q(1, 1) = std_yawdd_2;
	P_aug.bottomRightCorner(2, 2) = Q;
	//std::cout << "P_aug:" << std::endl << P_aug << std::endl;

	//create square root matrix
	MatrixXd A = P_aug.llt().matrixL();
	//std::cout << "A:" << std::endl << A << std::endl;
	//create augmented sigma points
	Xsig_aug.col(0) << x_aug;

	const float sqrt_lamb_plus_n_aug = sqrt(lambda + n_aug);

	for (int i = 1; i < n_aug + 1; ++i) {
		Xsig_aug.col(i) << x_aug + sqrt_lamb_plus_n_aug * A.col(i - 1);
	}

	const int upper_limit = n_aug + 1;
	for (int i = upper_limit; i < 2 * n_aug + 1; ++i) {
		Xsig_aug.col(i) << x_aug - sqrt_lamb_plus_n_aug * A.col(i - upper_limit);
	}

	/*******************************************************************************
	* Student part end
	******************************************************************************/

	//print result
	//std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
	//getchar();

	//write result
	*Xsig_out = Xsig_aug;

	/* expected result:
	Xsig_aug =
	5.7441  5.85768   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441  5.63052   5.7441   5.7441   5.7441   5.7441   5.7441   5.7441
	1.38  1.34566  1.52806     1.38     1.38     1.38     1.38     1.38  1.41434  1.23194     1.38     1.38     1.38     1.38     1.38
	2.2049  2.28414  2.24557  2.29582   2.2049   2.2049   2.2049   2.2049  2.12566  2.16423  2.11398   2.2049   2.2049   2.2049   2.2049
	0.5015  0.44339 0.631886 0.516923 0.595227   0.5015   0.5015   0.5015  0.55961 0.371114 0.486077 0.407773   0.5015   0.5015   0.5015
	0.3528 0.299973 0.462123 0.376339  0.48417 0.418721   0.3528   0.3528 0.405627 0.243477 0.329261  0.22143 0.286879   0.3528   0.3528
	0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641        0
	0        0        0        0        0        0        0  0.34641        0        0        0        0        0        0 -0.34641
	*/

}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {

	//set state dimension
	int n_x = 5;

	//set augmented dimension
	int n_aug = 7;

	//create example sigma point matrix
	MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
	Xsig_aug <<
		5.7441, 5.85768, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.63052, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441,
		1.38, 1.34566, 1.52806, 1.38, 1.38, 1.38, 1.38, 1.38, 1.41434, 1.23194, 1.38, 1.38, 1.38, 1.38, 1.38,
		2.2049, 2.28414, 2.24557, 2.29582, 2.2049, 2.2049, 2.2049, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.2049, 2.2049,
		0.5015, 0.44339, 0.631886, 0.516923, 0.595227, 0.5015, 0.5015, 0.5015, 0.55961, 0.371114, 0.486077, 0.407773, 0.5015, 0.5015, 0.5015,
		0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.3528, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.3528,
		0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641, 0,
		0, 0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641;

	//create matrix with predicted sigma points as columns
	MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

	double delta_t = 0.1; //time diff in sec
	/*******************************************************************************
	* Student part begin
	******************************************************************************/

	//predict sigma points
	//avoid division by zero
	//write predicted sigma points into right column

	//predict sigma points
	const int upper_limit = 2 * n_aug + 1;
	for (int i = 0; i < upper_limit; i++)	{
		//extract values for better readability
		double p_x = Xsig_aug(0, i);
		double p_y = Xsig_aug(1, i);
		double v = Xsig_aug(2, i);
		double yaw = Xsig_aug(3, i);
		double yawd = Xsig_aug(4, i);
		double nu_a = Xsig_aug(5, i);
		double nu_yawdd = Xsig_aug(6, i);

		//predicted state values
		double p_x_p, p_y_p;
		double yawd_delta_t = yawd * delta_t;			 
		double sin_yaw = sin(yaw);
		double cos_yaw = cos(yaw);

		//avoid division by zero
		if (fabs(yawd) > 0.001) {
			p_x_p = p_x + v / yawd * (sin(yaw + yawd_delta_t) - sin_yaw);
			p_y_p = p_y + v / yawd * (cos_yaw - cos(yaw + yawd_delta_t));
		}
		else {
			p_x_p = p_x + v * delta_t * cos_yaw;
			p_y_p = p_y + v * delta_t * sin_yaw;
		}

		double v_p = v;
		double yaw_p = yaw + yawd_delta_t;
		double yawd_p = yawd;

		//add noise
		double half_nu_a_delt_2 = 0.5 * nu_a * delta_t * delta_t;
		p_x_p = p_x_p + half_nu_a_delt_2 * cos_yaw;
		p_y_p = p_y_p + half_nu_a_delt_2 * sin_yaw;
		v_p = v_p + nu_a * delta_t;

		yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
		yawd_p = yawd_p + nu_yawdd*delta_t;

		//write predicted sigma point into right column
		Xsig_pred(0, i) = p_x_p;
		Xsig_pred(1, i) = p_y_p;
		Xsig_pred(2, i) = v_p;
		Xsig_pred(3, i) = yaw_p;
		Xsig_pred(4, i) = yawd_p;
	}

	/*******************************************************************************
	* Student part end
	******************************************************************************/

	//print result
	std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
	getwchar();

	//write result
	*Xsig_out = Xsig_pred;

}
