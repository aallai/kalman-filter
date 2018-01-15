#include "kalman_filter.h"
#include "tools.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {
  P_.setIdentity(4, 4);
  F_.setIdentity(4, 4);
  Q_.setIdentity(4, 4);
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict(double delta_T) {
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

    // Map predicted state into polar coordinates.
    double px = x_(0), py = x_(1), vx = x_(2), vy = x_(3);
    auto polar = VectorXd(3);
    polar(0) = sqrt((px*px) + py*py);
    polar(1) = atan2(py, px);
    polar(2) = (px*vx + py*vy) / dbz_guard(sqrt(px*px+py*py));

    VectorXd innovation = z - polar;

    // Normalize angles.
    innovation(1) = atan2(sin(innovation(1)), cos(innovation(1)));

    update(innovation);
}

void KalmanFilter::update(const VectorXd &innovation) {
    MatrixXd S = H_ * P_ * H_.transpose() + R_;
    MatrixXd K = P_ * H_.transpose() * S.inverse();

    x_ += K * innovation;
    P_ -= K * H_ * P_;
}
