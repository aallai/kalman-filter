#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
  P_.setIdentity();
  F_.setIdentity();
  Q_.setIdentity();
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict(long long delta_T) {
    F_ << 1, 0, delta_T , 0,
          0, 1, 0, delta_T,
          0, 0, 1, 0,
          0, 0, 0, 1;

    Q_ << pow(delta_T, 4)/4* sig_a, 0, pow(delta_T, 3)/2*sig_a, 0,
         0, pow(delta_T, 4)/4* sig_a, 0, pow(delta_T, 3)/2*sig_a,
         pow(delta_T, 3)/2*sig_a, 0, pow(delta_T, 2) * sig_a, 0,
         0, pow(delta_T, 3)/2*sig_a, 0, pow(delta_T, 2) * sig_a;

    x_ = F_ * x_;
    P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
    update(z - H_ * x_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // Project predicted state into polar coordinates.
    double px = x_[0], py = x_[1], vx = x_[2], vy = x_[3];
    auto polar = VectorXd(3);
    polar(0) = sqrt((px*px) + py*py);
    polar(1) = atan2(py, px); // divide by zero
    polar(2) = (px*vx + py*vy) / sqrt(px*px+py*py); // divide by zero

    VectorXd innovation = z - polar;

    // Normalize angles.
    innovation(1) = atan2(sin(innovation[1]), cos(innovation[1]));

    update(innovation);
}

void KalmanFilter::update(const VectorXd &innovation) {
    auto S = R_ + H_ * P_ * H_.transpose();
    auto K = P_ * H_.transpose() * S.inverse();
    x_ += K * innovation;
    P_ -= K * H_ * P_;
}
