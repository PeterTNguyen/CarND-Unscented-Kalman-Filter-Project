#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  is_initialized_ = false;

  // initial state vector
  x_ = VectorXd(5);
  x_aug_ = VectorXd(7);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_aug_ = MatrixXd(7, 7);

  // init sigma pt matrices and lengths
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3-n_x_;
  lambda_aug_ = 3-n_aug_;
  n_sig_ = 2*n_x_ + 1;
  n_sig_aug_ = 2*n_aug_ + 1;
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);

  // Weigths calculation
  weights_ = VectorXd(n_sig_aug_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2*n_aug_).setConstant(1/(lambda_ + n_aug_)/2);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:
  //test

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_)
  {
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
      float range = meas_package.raw_measurements_[0];
      float velocity = meas_package.raw_measurements_[2];
      float bearing = meas_package.raw_measurements_[1];
      float tmp_x = range * cos(bearing);
      float tmp_y = range * sin(bearing);
      float tmp_vx = velocity * cos(bearing);
      float tmp_vy = velocity * sin(bearing);
      x_ << tmp_x, tmp_y, velocity, 0, 0;// TODO init yadd? radar ang =/= psi


      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }
    if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {

      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1], 0, 0, 0;

      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }
  }

  // Process measurements after initialization
  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  //calculate square root of P

  //Augmented state vector
  x_aug_.head(n_x_) = x_;
  x_aug_(n_x_) = 0;
  x_aug_(n_x_+1) = 0;

  //Augmented covariance matrix
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(n_x_,n_x_) = std_a_ * std_a_;
  P_aug_(n_x_+1,n_x_+1) = std_yawdd_ * std_yawdd_;

  //Generate Sigma Pts
  MatrixXd A = P_aug_.llt().matrixL();
  A *= sqrt(lambda_aug_ - n_aug_);
  Xsig_aug_ = x_aug_.replicate(1,n_sig_aug_);
  Xsig_aug_.block(0,1,n_aug_,n_aug_) += A;
  Xsig_aug_.block(0,1+n_aug_,n_aug_,n_aug_) -= A;

  //Prediction Step
  VectorXd pred1 = VectorXd(n_x_);
  pred1(2) = 0;
  pred1(4) = 0;
  VectorXd pred2 = VectorXd(n_x_);
  Xsig_pred_ = Xsig_aug_.topLeftCorner(n_x_, n_sig_aug_);
  double v, psi, psi_dot, nu_a, nu_psi;
  double delta_t2 = delta_t * delta_t;
  for(int i = 0; i < n_sig_aug_; i++)
  {
    v = Xsig_aug_(2, i);
    psi = Xsig_aug_(3, i);
    psi_dot = Xsig_aug_(4, i);
    nu_a = Xsig_aug_(5, i);
    nu_psi = Xsig_aug_(6, i);
    //avoid division by zero
    if(psi_dot == 0.0)
    {
      std::cout << "Error, can't divide by zero" << std::endl;
      return;
    }
    //first order
    pred1(0) = v/psi_dot*(sin(psi + psi_dot*delta_t ) - sin(psi));
    pred1(1) = v/psi_dot*(-cos(psi + psi_dot*delta_t ) + cos(psi));
    pred1(3) = psi_dot*delta_t;

    //second order
    pred2(0) = delta_t2*cos(psi)*nu_a/2;
    pred2(1) = delta_t2*sin(psi)*nu_a/2;
    pred2(2) = delta_t*nu_a;
    pred2(3) = delta_t2*nu_psi/2;
    pred2(4) = delta_t*nu_psi;

    //update predicted sigma points
    Xsig_pred_.col(i) += pred1 + pred2;
  }

  //Calculate Predicted state and state covariance
  //predict state mean
  x_ = Xsig_pred_*weights_;
  //predict state covariance matrix
  MatrixXd meanRep = x_.replicate(1,2*n_aug_+1);
  MatrixXd sigma_diff = Xsig_pred_ - meanRep;
  P_ = sigma_diff*weights_.asDiagonal()*sigma_diff.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
