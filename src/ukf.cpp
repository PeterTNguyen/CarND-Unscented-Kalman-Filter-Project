#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <fstream>

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

  // init sigma pt matrices and lengths
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3-n_x_;
  lambda_aug_ = 3-n_aug_;
  n_sig_ = 2*n_x_ + 1;
  n_sig_aug_ = 2*n_aug_ + 1;
  n_zr_ = 3;
  n_zl_ = 2;
  Xsig_pred_ = MatrixXd(n_x_, n_sig_aug_);

  // initial state vector
  x_ = VectorXd(5);
  x_aug_ = VectorXd(7).setZero();

  // initial covariance matrix
  P_ = MatrixXd(5, 5).setZero();
  P_aug_ = MatrixXd(7, 7).setZero();

  // Weigths calculation
  weights_ = VectorXd(n_sig_aug_);
  weights_(0) = lambda_aug_ / (lambda_aug_ + n_aug_);
  weights_.tail(2*n_aug_).setConstant(1/(lambda_aug_ + n_aug_)*0.5);
  cout << weights_<< endl;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/6;

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

  // Radar Measurement covariance
  R_lidar_      = MatrixXd(n_zl_,n_zl_).setZero();
  R_lidar_(0,0) = std_laspx_*std_laspx_;
  R_lidar_(1,1) = std_laspy_*std_laspy_;

  // Radar Measurement covariance
  R_radar_      = MatrixXd(n_zr_,n_zr_).setZero();
  R_radar_(0,0) = std_radr_*std_radr_;
  R_radar_(1,1) = std_radphi_*std_radphi_;
  R_radar_(2,2) = std_radrd_*std_radrd_;

  //Kalman filter
  H_lidar_ = MatrixXd(2,5);
  H_lidar_ << 1, 0, 0, 0, 0,
              0, 1, 0, 0, 0;

  //Q matrix
  P_aug_(n_x_,n_x_) = std_a_ * std_a_;
  P_aug_(n_x_+1,n_x_+1) = std_yawdd_ * std_yawdd_;

  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
}

UKF::~UKF(){}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_)
  {
    P_ << 2, 0, 0, 0, 0,
          0, 2, 0, 0, 0,
          0, 0, 2, 0, 0,
          0, 0, 0, 2, 0,
          0, 0, 0, 0, 2;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
    {
      float range = meas_package.raw_measurements_[0];
      float velocity = meas_package.raw_measurements_[2];
      float bearing = meas_package.raw_measurements_[1];
      float tmp_x = range * cos(bearing);
      float tmp_y = range * sin(bearing);
      x_ << tmp_x, tmp_y, velocity, 0, 0;

      cout << "RADAR: "<< tmp_x << ", " << tmp_y << ", "<< velocity<< endl;

      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
    {

      x_ << meas_package.raw_measurements_[0],
            meas_package.raw_measurements_[1], 0, 0, 0;
      cout << "LIDAR: " << meas_package.raw_measurements_[0] << ", " << 
            meas_package.raw_measurements_[1] << endl;

      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
    }
  }

  // Process measurements after initialization
  double dt = (meas_package.timestamp_ - time_us_)/1000000.0;
  Prediction(dt);
  time_us_ = meas_package.timestamp_;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_)
  {
    UpdateRadar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)
  {
    UpdateLidar(meas_package);
  }
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

  //Augmented covariance matrix
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;

  //Generate Sigma Pts
  MatrixXd A = P_aug_.llt().matrixL();
  A = sqrt(lambda_aug_ + n_aug_)*A;
  Xsig_aug_ = x_aug_.replicate(1,n_sig_aug_);
  Xsig_aug_.block(0,1,n_aug_,n_aug_) += A;
  Xsig_aug_.block(0,1+n_aug_,n_aug_,n_aug_) -= A;

  //Prediction Step
  VectorXd pred1 = VectorXd(n_x_);
  pred1(2) = 0;
  pred1(4) = 0;
  VectorXd pred2 = VectorXd(n_x_);
  Xsig_pred_ = Xsig_aug_.topLeftCorner(n_x_, n_sig_aug_);
  double v, yaw, yaw_dot, nu_a, nu_yaw;
  double delta_t2 = delta_t * delta_t;
  for(int i = 0; i < n_sig_aug_; i++)
  {
    v = Xsig_aug_(2, i);
    yaw = Xsig_aug_(3, i);
    yaw_dot = Xsig_aug_(4, i);
    nu_a = Xsig_aug_(5, i);
    nu_yaw = Xsig_aug_(6, i);
    //first order
    if(fabs(yaw_dot) < 0.001)
    {
      pred1(0) = v*cos(yaw)*delta_t;
      pred1(1) = v*sin(yaw)*delta_t;
    }
    else
    {
      pred1(0) = v/yaw_dot*(sin(yaw + yaw_dot*delta_t ) - sin(yaw));
      pred1(1) = v/yaw_dot*(cos(yaw) - cos(yaw + yaw_dot*delta_t ));
    }
    pred1(3) = yaw_dot*delta_t;

    //second order
    pred2(0) = delta_t2*cos(yaw)*nu_a*0.5;
    pred2(1) = delta_t2*sin(yaw)*nu_a*0.5;
    pred2(2) = delta_t*nu_a;
    pred2(3) = delta_t2*nu_yaw*0.5;
    pred2(4) = delta_t*nu_yaw;

    //update predicted sigma points
    Xsig_pred_.col(i) += pred1 + pred2;
  }

  //Calculate Predicted state and state covariance
  //predict state mean
  x_pred_ = Xsig_pred_*weights_;
  //predict state covariance matrix
  MatrixXd meanRep = x_pred_.replicate(1,n_sig_aug_);
  MatrixXd sigma_diff = Xsig_pred_ - meanRep;
  for(int i = 0; i < n_sig_aug_; i++)
  {
    UnwrapAngle(sigma_diff(3,i));
  }
  P_pred_ = sigma_diff*weights_.asDiagonal()*sigma_diff.transpose();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  //Kalman filter
  VectorXd z = meas_package.raw_measurements_;
  VectorXd z_pred = H_lidar_ * x_pred_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_lidar_.transpose();
  MatrixXd S = H_lidar_ * P_pred_ * Ht + R_lidar_;
  MatrixXd PHt = P_pred_ * Ht;
  MatrixXd K = PHt * S.inverse();

  //new estimate
  x_ = x_pred_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_lidar_) * P_pred_;

  //Calculate NIS
  MatrixXd nis_lidar = (z - z_pred).transpose()*S.inverse()*(z - z_pred);
  nis_lidar_ = nis_lidar(0,0);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  //transform sigma points into measurement space
  //Radar sigma pts
  MatrixXd Zsig = MatrixXd(n_zr_, n_sig_aug_).setZero();
  float px, py, v, yaw;
  for(int i = 0; i < n_sig_aug_; i++)
  {
    px = Xsig_pred_(0,i);
    py = Xsig_pred_(1,i);
    v = Xsig_pred_(2,i);
    yaw = Xsig_pred_(3,i);
    if((px*px + py*py) == 0.0)
    {
      cout << "ERROR: divide by zeros" << endl;
      return;
    }
    //Radar transformation
    Zsig(0,i) = sqrt(px*px + py*py);
    Zsig(1,i) = atan2(py, px);
    Zsig(2,i) = (v*px*cos(yaw) + v*py*sin(yaw))/Zsig(0,i);

  }
  //calculate mean predicted measurement
  VectorXd z_pred = Zsig * weights_;

  //calculate innovation covariance matrix S
  MatrixXd z_pred_mat = z_pred.replicate(1,n_sig_aug_);
  MatrixXd z_sigma_diff = Zsig - z_pred_mat;
  for(int i = 0; i < n_sig_aug_; i++)
  {
    UnwrapAngle(z_sigma_diff(1,i));
  }

  //Radar measurement covariance
  MatrixXd S = z_sigma_diff*weights_.asDiagonal()*z_sigma_diff.transpose() + R_radar_;

  //Cros correlation matrix
  MatrixXd sigma_state_diff = Xsig_pred_ - x_pred_.replicate(1,n_sig_aug_);
  MatrixXd sigma_meas_diff = (Zsig - z_pred.replicate(1,n_sig_aug_));
  for(int i = 0; i < n_sig_aug_; i++)
  {
    UnwrapAngle(sigma_state_diff(3,i));
    UnwrapAngle(sigma_meas_diff(1,i));
  }
  MatrixXd Tc = sigma_state_diff*weights_.asDiagonal()*sigma_meas_diff.transpose();

  //Kalman Gain
  MatrixXd K = Tc * S.inverse();

  //Get Measurement
  VectorXd z = meas_package.raw_measurements_;

  //Update to x state mean and P state covariance
  VectorXd z_diff = z-z_pred;
  UnwrapAngle(z_diff(1));
  x_ = x_pred_ + K*z_diff;
  P_ = P_pred_ - K*S*K.transpose();

  //Calculate NIS
  MatrixXd nis_radar = z_diff.transpose()*S.inverse()*z_diff;
  nis_radar_ = nis_radar(0,0);
}


void UKF::UnwrapAngle(double& ang)
{
  while (ang > M_PI) ang -= 2.* M_PI;
  while (ang < -M_PI) ang += 2.* M_PI;
}
