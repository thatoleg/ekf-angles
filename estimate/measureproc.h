/**
Copyright 2011 Oleg Popov

This file is part of EKF-Angles.

EKF-Angles is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

EKF-Angles is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with EKF-Angles. If not, see http://www.gnu.org/licenses/.

**/

#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#ifndef MEASUREPROC_H_
#define MEASUREPROC_H_



#define VECTOR3 3
#define DoubleMAX_VALUE DBL_MAX

typedef struct triple {
  int triple[3];
} triple;


typedef struct MV3 { double x[3]; } MV3;
typedef struct MV2 { double x[2]; } MV2;
typedef struct MMAT3 { MV3 x[3]; } MMAT3;




const static double GRAVITATION_ACCELERATION_CONSTANT = 9.80665;
const static double DEFAULT_ACCEL_DISTURBANCE_THRESHOLD_MS2 = 2.941995; //0.3 * GRAVITATION_ACCELERATION_CONSTANT;

const static double DEFAULT_SIGMA_WINDOW_SECONDS = 0.2;

const static double DOUBLE_DEFAULT_ACCELETOMETER_NOISE = 0.0859355543230206;

const static double DEFAULT_GYRO_NOISE = 0.004;



static MV3 EPSILON_VECTOR ;

/*
 * On Samsung Galaxy S2 that lies horizontally on the balcony floor the mean
 * of gyroscope data is 0.00896271374747443 -0.00517417406954337
 * -0.00371192004660972 stdev of gyroscope data is 0.00344555661801865
 * 0.00303404060391456 0.00431348773767872
 *
 * so we set the processNoiseCovarianceQ to the stdev values
 *
 * BUT
 *
 * Simon Kalman filter book says that by adding some fictious noise to the
 * measurementnoise, we can compensate for modelling errors. Since our gyro
 * had bias and the model is nonlinear and sensitive to double precision
 * calculations accuracy - let's add some noise. Let this noise to be
 * estimated
 *
 * as triple compared to the actual gyro error
 */
// 1.618d*1.618d; since we are integrating over time, the error will be 10
// times lower
// but for unmodeled bias + numerical errors + modelling errors add x2.7
// multiplier
// so it will be 3d/10d
const static double FICTIOUS_NOISE_FACTOR = 2.7; // 1.0d;

static MV3 INITIAL_Z_ESTIMATE ;

//const long SKIP_AFTER = 1000;

/**
 * diagonal matrix constant
 */
static MMAT3 diagonalOne;

typedef struct kalman {

  MV3 measurementNoiseCovarianceR;

  MV3 processNoiseCovarianceQ;

  MV3 aVec;

  MV3 aVecGravityNormalized;

  MV3 currentCovariancePkEstimate;

  MV2 currentRollPitchEstimate;

  MV3 currentXEstimate;
  MV3 currentZEstimate;

  double idleTimeLength;

  bool isCurrentlyAccelDisturbanePresent;
  MMAT3 kalmanA;

  MV3 lastGoodZEstimate;

  MV3 lastInputAcceleration;

  MV3 lastInputAngularVelocity;
  MMAT3 pma;

  double sigmaWindowLengthSeconds;

  long long skipperCounter;

  MV3 wkhNormDt;

  /**
   * Angular velocity vector wx, wy, wz
   */
  MV3 wVec;

} KF;

void initKf(KF* pkf);
void processMeasurement(KF* pkf, double a1, double a2, double a3, double w1,
                        double w2, double w3, double dt);

void rotateVectorByA(KF* pkf, MV3* pv, MV3* pvout);
void initLibrary();
int currentSigmaValue(KF * pkf) ;
int isErrorVector(double * vec, int size);

#endif /* MEASUREPROC_H_ */
