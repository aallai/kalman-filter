#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
          0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
          0, 0.0009, 0,
          0, 0, 0.09;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

VectorXd polar_to_cartesian(VectorXd polar)
{
    double rho = polar[0], phi = polar[1], rho_dot = polar[2];
    auto px = sqrt((rho*rho) / dbz_guard(1 - tan(phi))); // divide by zero
    auto py = sqrt((rho*rho) - (px * px));
    auto cartesian = VectorXd(4);
    cartesian << px, py, 0.0, 0.0;
    return cartesian;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
        TODO:
          * Initialize the state ekf_.x_ with the first measurement.
          * Create the covariance matrix.
          * Remember: you'll need to convert radar from polar to cartesian coordinates.
        */
        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
          ekf_.x_ = polar_to_cartesian(measurement_pack.raw_measurements_);
        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
          ekf_.x_ << measurement_pack.raw_measurements_(0), measurement_pack.raw_measurements_(1), 0, 0;
        }

        previous_timestamp_ = measurement_pack.timestamp_;

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     TODO:
       * Update the state transition matrix F according to the new elapsed time.
        - Time is measured in seconds.
       * Update the process noise covariance matrix.
       * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    ekf_.Predict((measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0);

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        ekf_.R_ = MatrixXd(3, 3);
        ekf_.R_ << R_radar_;
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        ekf_.R_ = MatrixXd(2, 2);
        ekf_.R_ << R_laser_;
        ekf_.H_ = MatrixXd(2, 4);
        ekf_.H_ << 1, 0, 0, 0,
                   0, 1, 0, 0;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;

    previous_timestamp_ = measurement_pack.timestamp_;
}
