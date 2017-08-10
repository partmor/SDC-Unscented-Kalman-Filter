#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;


UKF::UKF() {

  // first measurement still not processed
  is_initialized_ = false;

  // initial time, in us
  time_us_ = 0;

  /**
   * Sensor-related attributes
   */

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

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

  // Laser measurement dimensions
  n_z_laser_ = 2;

  // Radar measurement dimensions
  n_z_radar_ = 3;

  /**
   * Process-related attributes
   */

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  /**
   * Sigma points and state-related attributes
   */

  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // Sigma point spreading parameter for augmented state
  lambda_aug_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_aug_/(lambda_aug_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug_+lambda_aug_);
    weights_(i) = weight;
  }

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

/**
*  Initialization
*/

  if (!is_initialized_) {
    // first measurement

    float px, py, v, yaw, yaw_dot;
    x_ = VectorXd(n_x_);
    P_ = MatrixXd(n_x_,n_x_);

    if (use_radar_ && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {

      //recover measurement parameters
      float rho = meas_package.raw_measurements_(0);
      float phi = meas_package.raw_measurements_(1);
      float rho_dot = meas_package.raw_measurements_(2);

      // positions px,py can be recovered rigorously
      px = rho * cos(phi);
      py = rho * sin(phi);

      // rho_dot is the projection of the object's velocity v on the rho direction.
      // the v vector could be any one yielding a projection equal to rho_dot,
      // i.e, there are infinite vectors that can project onto rho_dot
      v = rho_dot;
      yaw = 0;
      yaw_dot = 0;

      P_.fill(0.0);
      P_(0,0) = std_radr_;
      P_(1,1) = std_radr_;
      P_(2,2) = rho_dot;
      P_(3,3) = M_PI;
      P_(4,4) = 0.1;

    }
    else if (use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {

      px = meas_package.raw_measurements_(0);
      py = meas_package.raw_measurements_(1);
      v = 5;
      yaw = 0;
      yaw_dot = 0;

      P_.fill(0.0);
      P_(0,0) = std_laspx_;
      P_(1,1) = std_laspy_;
      P_(2,2) = 5;
      P_(3,3) = M_PI;
      P_(4,4) = 0.1;

    }
    else{
      // it could happen that first measurement comes from radar, but user_radar_ is set to false.
      // in this case we should not mark as initialized or set timestamp until we receive the first laser measurement
      return;
    }

    x_ << px, py, v, yaw, yaw_dot;

    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    // done initialization, no need to predict or update
    return;
  }

  /**
  *  Prediction
  */

  // elapsed time between current and previous measurement (in seconds)
  float dt = (meas_package.timestamp_ - time_us_) / 1e6;
  time_us_ = meas_package.timestamp_;

  // perform prediction
  Prediction(dt);

  /*
  * Radar updates
  */
  if (use_radar_ && (meas_package.sensor_type_ == MeasurementPackage::RADAR)) {
    UpdateRadar(meas_package);
  }
  /*
  * Laser updates
  */
  else if (use_laser_ && (meas_package.sensor_type_ == MeasurementPackage::LASER)) {
    UpdateLidar(meas_package);
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // 1. generate sigma points (for augmented state)
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  AugmentedSigmaPoints(&Xsig_aug);

  // 2. predict sigma points: apply process function to the augmented sigma points
  SigmaPointPrediction(Xsig_aug, delta_t);

  // 3. predict mean and covariance: calculate a priori state, x(k+1|k) and
  //    P(k+1|k) given the predicted sigma points
  PredictMeanAndCovariance();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);
  //create vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_laser_);
  //create matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z_laser_, n_z_laser_);

  // 1. predict measurement: apply measurement model transformation to predicted
  //    state k+1|k (given by the predicted sigma points)
  PredictLidarMeasurement(&Zsig, &z_pred, &S);

  // 2. update state: calculate posterior state k+1|k+1
  UpdateState(meas_package.raw_measurements_, z_pred, S, Zsig);

  // Finally, compute Normalized Innovation Squared
  NIS_laser_ = CalculateNIS(meas_package.raw_measurements_, z_pred, S);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  //create vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  //create matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z_radar_, n_z_radar_);

  // 1. predict measurement: apply measurement model transformation to predicted
  //    state k+1|k (given by the predicted sigma points)
  PredictRadarMeasurement(&Zsig, &z_pred, &S);

  // 2. update state: calculate posterior state k+1|k+1
  UpdateState(meas_package.raw_measurements_, z_pred, S, Zsig);

  // Finally, compute Normalized Innovation Squared
  NIS_radar_ = CalculateNIS(meas_package.raw_measurements_, z_pred, S);
}

void UKF::GenerateSigmaPoints(MatrixXd* Xsig_out) {

  //create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);

  //calculate square root of P
  MatrixXd A = P_.llt().matrixL();

  //set first column of sigma point matrix
  Xsig.col(0)  = x_;

  //set remaining sigma points
  for (int i = 0; i < n_x_; i++)
  {
    Xsig.col(i+1)     = x_ + sqrt(lambda_+n_x_) * A.col(i);
    Xsig.col(i+1+n_x_) = x_ - sqrt(lambda_+n_x_) * A.col(i);
  }

  //write result
  *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_aug_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_aug_+n_aug_) * L.col(i);
  }

  //write result
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd Xsig_aug, double delta_t) {

  //predict sigma points
  for (int i = 0; i < 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

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

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovariance() {

  //predicted state mean
  x_.fill(0.0);
  //iterate over sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  //if estimated velocity is negative, invert sign and add pi to the predicted yaw
  if (x_(2) < 0){
    x_(2) = -x_(2);
    x_(3) += M_PI;
  }
  //normalize yaw component to range [-pi, pi]
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;

  //predicted state covariance matrix
  P_.fill(0.0);
  //iterate over sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
  }
}

void UKF::PredictRadarMeasurement(MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out) {

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // extract values for better readability
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_radar_,n_z_radar_);
  R.fill(0.0);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  S = S + R;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

void UKF::PredictLidarMeasurement(MatrixXd* Zsig_out, VectorXd* z_out, MatrixXd* S_out) {

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);

  //transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    // measurement model
    Zsig(0,i) = Xsig_pred_(0,i); // p_x
    Zsig(1,i) = Xsig_pred_(1,i); // p_y
  }

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_laser_);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_laser_,n_z_laser_);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S = S + weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_laser_,n_z_laser_);
  R.fill(0.0);
  R(0,0) = std_laspx_*std_laspx_;
  R(1,1) = std_laspy_*std_laspy_;

  S = S + R;

  //write result
  *z_out = z_pred;
  *S_out = S;
  *Zsig_out = Zsig;
}

void UKF::UpdateState(VectorXd z, VectorXd z_pred, MatrixXd S, MatrixXd Zsig) {

  // infer measurement dimensions
  int n_z = z.size();

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 sigma points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  //if estimated velocity is negative, invert sign and add pi to the predicted yaw
  if (x_(2) < 0){
    x_(2) = -x_(2);
    x_(3) += M_PI;
  }
  //normalize yaw component to range [-pi, pi]
  while (x_(3)> M_PI) x_(3)-=2.*M_PI;
  while (x_(3)<-M_PI) x_(3)+=2.*M_PI;
  P_ = P_ - K * S * K.transpose();
}

double UKF::CalculateNIS(VectorXd z, VectorXd z_pred, MatrixXd S){
  VectorXd z_diff = z - z_pred;
  double e = z_diff.transpose() * S.inverse() * z_diff;
  return e;
}
