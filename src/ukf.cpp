#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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
  
  /**
   * End DO NOT MODIFY section for measurement noise values 
   */
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  // If this is the first measurement, initialize state vector (x_) and
  // covariance matrix (P_) then toggle initialized flag and return
  if (!is_initialized_) 
  {
    // Is measurement RADAR or LIDAR?
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      // Initialize state vector (x_)
      x_.fill(0);
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);

      // Initialize covariance matrix (P_)
      P_  <<  std_laspx_*std_laspx_,  0,                      0,                    0,                        0,
              0,                      std_laspy_*std_laspy_,  0,                    0,                        0,
              0,                      0,                      std_radr_*std_radr_,  0,                        0,
              0,                      0,                      0,                    std_radphi_*std_radphi_,  0,
              0,                      0,                      0,                    0,                        std_radrd_*std_radrd_;

      // Set initialized flag to true
      is_initialized_ = true;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Initialize state vector (x_)
      x_.fill(0);
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      x_(2) = meas_package.raw_measurements_(2);

      // Initialize covariance matrix (P_)
      P_  <<  1,  0,  0,  0,  0,
              0,  1,  0,  0,  0,
              0,  0,  1,  0,  0,
              0,  0,  0,  1,  0,
              0,  0,  0,  0,  1;

      // Set initialized flag to true
      is_initialized_ = true;
    }
    else
    {
      cout << "ERROR: Unknown Sensor Type." << endl;
      return;
    }
  }
  else  // Not First measurement
  {
    // Is measurement RADAR or LIDAR?
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
    }
    else
    {
      cout << "ERROR: Unknown Sensor Type." << endl;
      return;
    }

  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // Initialize an augmented state vector, covariance matrix, 
  // and Augmented Sigma Point Matrix
  VectorXd x_aug = VectorXd(7);
  MatrixXd P_aug = MatrixXd(7, 7);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Fill in augmented state vector (x_aug)
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  // Fill in augmented covariance matrix (P_aug)
  P_aug.fill(0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;

  // Calculate the square root of the Augmented Covariance Matrix
  MatrixXd PaugSqrt = P_aug.llt().matrixL();

  // Calculate and fill in the Augmented Sigma Point Matrix
  // i.e. Generate Sigma Points
  Xsig_aug.col(0) = x_aug;
  for (int i = 1; i < (n_aug_+1); i++)
  {
    Xsig_aug.col(i)         = x_aug + sqrt(lambda+n_aug_) * PaugSqrt.col(i-1);
    Xsig_aug.col(i+n_aug_)  = x_aug - sqrt(lambda+n_aug_) * PaugSqrt.col(i-1);
  }

  // Calculate (1/2)t^2 for multiple use
  float delT2 = 0.5 * delta_t * delta_t;

  // For each column in the Augmented Sigma Point Matrix, calculate the predicted sigma points
  // i.e. Predict Sigma Points
  for (int i = 0; i < (2 * n_aug_ + 1); i++)
  {
    float vk = Xsig_aug(2,i);
    float psi = Xsig_aug(3,i);
    float psiDot = Xsig_aug(4,i);
    float nu_a = Xsig_aug(5,i);
    float nu_psiDot = Xsig_aug(6,i);

    // Prevent division by 0  (NOTE:  THere should by some minimum value threshold here instead)
    if (psiDot != 0)
    {
      Xsig_pred_(0,i) = Xsig_aug(0,i) + (vk/psiDot)*(sin(psi+psiDot*delta_t)-sin(psi)) + delT2*cos(psi)*nu_a; 
      Xsig_pred_(1,i) = Xsig_aug(1,i) + (vk/psiDot)*(-cos(psi+psiDot*delta_t)+cos(psi)) + delT2*sin(psi)*nu_a;
      Xsig_pred_(2,i) = Xsig_aug(2,i) + delta_t*nu_a;
      Xsig_pred_(3,i) = Xsig_aug(3,i) + psiDot*delta_t+delT2*nu_psiDdot;
      Xsig_pred_(4,i) = Xsig_aug(4,i) + delta_t*nu_psiDdot;
    }
    else  // psiDot == 0
    {
      Xsig_pred_(0,i) = Xsig_aug(0,i) + vk*cos(psi)*delta_t + delT2*cos(psi)*nu_a;
      Xsig_pred_(1,i) = Xsig_aug(1,i) + vk*sin(psi)*delta_t + delT2*sin(psi)*nu_a;
      Xsig_pred_(2,i) = Xsig_aug(2,i) + delta_t*nu_a;
      Xsig_pred_(3,i) = Xsig_aug(3,i) + psiDot*delta_t + delT2*nu_psiDdot;
      Xsig_pred_(4,i) = Xsig_aug(4,i) + delta_t*nu_psiDdot;
    }
  }

  // Calculate weights
  weights_ = VectorXd(2*n_aug_+1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < (2*n_aug_+1); i++)
  {
    weights_(i) = 1/(2*(lambda_+n_aug_));
  }

  // Calculate predictied state vector
  // i.e. Predict mean
  x_.fill(0);
  for (int i = 0; i < (2*n_aug_+1); i++)
  {
    x_ += weights_(i)*Xsig_pred_.col(i);    
  }

  // Calculate predicted covariance matrix
  // i.e. Predict Covariance
  P_.fill(0);
  for (int i = 0; i < (2*n_aug_+1); i++)
  {
    P_ += weights_(i)*(Xsig_pred_.col(i)-x_)*(Xsig_pred_.col(i)-x_).transpose();    
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
   MatrixXd Tc = MatrixXd(n_x_, 3);

   // Normalize angle
   while (zDiff > M_PI) zDiff -= 2 * M_PI;
   while (zDiff < M_PI) zDiff += 2 * M_PI;

   // Update state and mean covariance matrices
   x_ = x_ + K * (z - z_pred);
   P_ = P_ - K * S * K.transpose();
}