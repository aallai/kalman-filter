#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

// Can still result in INF.
double dbz_guard(double x)
{
  if (x == 0.0)
      return 0.000000001;
  return x;
}

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  auto sum = VectorXd(4);
  sum << 0, 0, 0, 0;

  for (int i = 0; i < estimations.size(); i++) {
    VectorXd diff = estimations[i] - ground_truth[i];
    diff = diff.array() * diff.array();
    sum += diff;
  }

  return (sum.array() / estimations.size()).sqrt();
}

// Calculate the jacobian of the state-to-measurement function with respect to state.
MatrixXd Tools::CalculateJacobian(const VectorXd& x) {

  double px = x(0), py = x(1), vx = x(2), vy = x(3);

  auto J = MatrixXd(3, 4);

  J << px/dbz_guard(sqrt(px*px+py*py)), py/dbz_guard(sqrt(px*px+py*py)), 0, 0,
       -py/dbz_guard(px*px+py*py), px/dbz_guard(px*px+py*py), 0, 0,
       py*(vx*py-vy*px)/dbz_guard(pow(px*px+py*py, 3.0/2.0)), px*(vy*px-vx*py)/dbz_guard(pow(px*px+py*py, 3.0/2.0)),
       px/dbz_guard(sqrt(px*px+py*py)), dbz_guard(py/sqrt(px*px+py*py));

  return J;
}
