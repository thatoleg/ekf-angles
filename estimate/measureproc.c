#include "measureproc.h"

/**
Copyright 2011 Oleg Popov

This file is part of EKF-Angles.

EKF-Angles is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

EKF-Angles is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with EKF-Angles. If not, see http://www.gnu.org/licenses/.

**/

bool Double_isNaN(double d) {
  return isnan(d);
}

double Math_abs (double d) {
  return fabs(d);
}

double Math_acos (double d) {
  return acos(d);
}

double Math_atan (double d) {
  return atan(d);
}

double Math_cos(double d) {
  return cos(d);
}

double Math_sin(double d) {
  return sin(d);
}

double Math_sqrt(double d) {
  return sqrt(d);
}


int isErrorVector(double * vec, int size) {
  int result = 0;
  int i;
  for (i = 0; i < size; i++) {
    if ( isnan(vec[i]) || !isfinite(vec[i])  ) {
      return 1;
    }
  }


  return result;

}


void updateSigmaValues(KF * pkf,double accelerationNorm, double dt) {

  if (Math_abs(accelerationNorm - GRAVITATION_ACCELERATION_CONSTANT) < DEFAULT_ACCEL_DISTURBANCE_THRESHOLD_MS2) {
    pkf->idleTimeLength += dt;
  } else {
    pkf->idleTimeLength = EPSILON_VECTOR.x[0]; //assume it is 100 microseconds in order to avoid divisionby zero in another routines
  }

}

/**
 * Assign the values from vector b to the corresponing values in vector a
 *
 * @param va
 * @param vb
 */
void v3AssignVector( MV3*  __restrict pva, MV3    vb) {
  pva->x[0] = vb.x[0];
  pva->x[1] = vb.x[1];
  pva->x[2] = vb.x[2];
}

void v3AssignMatrix(MMAT3 *  pma, MMAT3 mb) {
  pma->x[0] =mb.x[0];
  pma->x[1] =mb.x[1];
  pma->x[2] =mb.x[2];


}

void v3DivideScalar(  MV3 inV3vec, double scalar,
                      MV3 * __restrict pOutV3Vec) {
  pOutV3Vec->x[0] = inV3vec.x[0] / scalar;
  pOutV3Vec->x[1] = inV3vec.x[1] / scalar;
  pOutV3Vec->x[2] = inV3vec.x[2] / scalar;
}

/*
 * in - va, vb out - vc vc[i] = va[i] + vb[i]
 */
void v3EbeAdd(  MV3 va,   MV3 vb,   MV3 * __restrict pvc) {
  pvc->x[0] = va.x[0] + vb.x[0];
  pvc->x[1] = va.x[1] + vb.x[1];
  pvc->x[2] = va.x[2] + vb.x[2];

}

/*
 * in - va, vb out - vc vc[i] = va[i]/vb[i]
 */
void v3EbeDivide(  MV3 va,   MV3 vb,   MV3 *  __restrict pvc) {
  pvc->x[0] = va.x[0] / vb.x[0];
  pvc->x[1] = va.x[1] / vb.x[1];
  pvc->x[2] = va.x[2] / vb.x[2];

}


		


/*
 * in - va, vb out - vc vc[i] = va[i]*vb[i]
 */
void v3EbeMultiply(  MV3 va,   MV3 vb,   MV3 * __restrict pvc) {
  pvc->x[0] = va.x[0] * vb.x[0];
  pvc->x[1] = va.x[1] * vb.x[1];
  pvc->x[2] = va.x[2] * vb.x[2];

}

/*
 * in - va, vb out - vc vc[i] = va[i]/vb[i]
 */
void v3EbeSubstract(  MV3 va,   MV3 vb,   MV3 * __restrict pvc) {
  pvc->x[0] = va.x[0] - vb.x[0];
  pvc->x[1] = va.x[1] - vb.x[1];
  pvc->x[2] = va.x[2] - vb.x[2];

}

double v3Norm(  MV3 v3vec) {
  return Math_sqrt(v3vec.x[0] * v3vec.x[0] + v3vec.x[1] * v3vec.x[1] + v3vec.x[2] 
                   * v3vec.x[2]);
}

bool v3IsZeroNorm(  MV3 inV3vec) {
  return Math_abs(v3Norm(inV3vec)) < EPSILON_VECTOR.x[0];
}

bool v3IsZeroVector(  MV3 inV3vec) {
  return inV3vec.x[0] < EPSILON_VECTOR.x[0] || inV3vec.x[1] < EPSILON_VECTOR.x[1]
    || inV3vec.x[2] < EPSILON_VECTOR.x[2];
}

void v3MultiplyScalar(  MV3 inV3vec, double scalar,
                        MV3 * __restrict pOutV3Vec) {
  pOutV3Vec->x[0] = scalar * inV3vec.x[0];
  pOutV3Vec->x[1] = scalar * inV3vec.x[1];
  pOutV3Vec->x[2] = scalar * inV3vec.x[2];
}


                





void destroy() {
  // TODO Auto-generated method stub

}



void initKf(KF* pkf) {

  memset(pkf, 0, sizeof(KF));

  /*
   * On Samsung Galaxy S2 that lies horizontally on the balcony floor the mean
   * of accelerometer data is 0.155163486717285 -0.217262627632744
   * 10.715392549481 stdev of accelerometer data is 0.0211647814647762
   * 0.0311250002004864 0.0426692683090581
   *
   * so we set the measurementNoiseCovarianceR to the stdev valuescause sdf\
   */
  pkf->measurementNoiseCovarianceR.x[0] = DOUBLE_DEFAULT_ACCELETOMETER_NOISE;
  pkf->measurementNoiseCovarianceR.x[1] = DOUBLE_DEFAULT_ACCELETOMETER_NOISE;
  pkf->measurementNoiseCovarianceR.x[2] = DOUBLE_DEFAULT_ACCELETOMETER_NOISE;


  // control weighting
  // matrix

  pkf->processNoiseCovarianceQ.x[0] = FICTIOUS_NOISE_FACTOR * DEFAULT_GYRO_NOISE;
  pkf->processNoiseCovarianceQ.x[1] = FICTIOUS_NOISE_FACTOR * DEFAULT_GYRO_NOISE;
  pkf->processNoiseCovarianceQ.x[2] = FICTIOUS_NOISE_FACTOR * DEFAULT_GYRO_NOISE;


  // state weighting
  // matrix
  /**
   * Acceleration vector ax, ay, az
   */
  pkf->aVec.x[0] = pkf->aVec.x[1] =  pkf->aVec.x[2] = 0  ;
  /**
   * Acceleration vector ax, ay, az divided by earth gravity
   */
  pkf->aVecGravityNormalized.x[0] = pkf->aVecGravityNormalized.x[1] = pkf->aVecGravityNormalized.x[2] = 0;

  pkf->currentCovariancePkEstimate.x[0] = 922.0557/1000.0;
  pkf->currentCovariancePkEstimate.x[1] = 958.3959/1000.0;
  pkf->currentCovariancePkEstimate.x[2] = 977.7422/1000.0;


  // 2.7998851412952314E-4, 3.132630280222944E-4, 4.3783599103602866E-4 };

  pkf->currentRollPitchEstimate.x[0] =
    pkf->currentRollPitchEstimate.x[1] = 0;


  v3AssignVector(&(pkf->currentXEstimate), INITIAL_Z_ESTIMATE);
  v3AssignVector(&(pkf->currentZEstimate), INITIAL_Z_ESTIMATE);



  // let's set it to max value - because we need even a very bad estimate to
  // begin our calculations with
  pkf->idleTimeLength = DEFAULT_SIGMA_WINDOW_SECONDS*2;

  /*
   *
   *
   * A =
   *
   * 1.0000 0.0000 0.0000 -0.0000 1.0000 -0.0000 -0.0000 0.0000 1.0000
   *
   * >> Q
   *
   * Q =
   *
   * 0.0090 0 0 0 0.0079 0 0 0 0.0113
   *
   * >> Qt=Q*dt
   *
   * Qt =
   *
   * 1.0e-03 *
   *
   * 0.1355 0 0 0 0.1193 0 0 0 0.1697
   *
   * >> P = dlyap(A,Qt)
   *
   * P =
   *
   * 922.0557 -37.6811 16.4091 -37.6811 958.3959 8.9598 16.4091 8.9598
   * 977.7422
   */

  pkf->isCurrentlyAccelDisturbanePresent = false;


  v3AssignVector( &(pkf->lastGoodZEstimate), INITIAL_Z_ESTIMATE);


  pkf->sigmaWindowLengthSeconds = DEFAULT_SIGMA_WINDOW_SECONDS;



  pkf->skipperCounter = 0;

  /*
   * variable to hold updateCovarianceResult - in order to avoid excessive
   * memory allocation
   */






  v3AssignMatrix( &(pkf->kalmanA), diagonalOne );

}










/**
 * w1, w2, w3 - are the components of the angular velocity return by the
 * Gyro
 *
 * @return
 */
MMAT3 createKalmanA(KF* pkf,   MV3 angularVelocities,
                     double timeStep) {


  MV3 wkh;
  MMAT3 m;

  v3MultiplyScalar(angularVelocities, timeStep, &wkh);

  double wkhDtNorm = v3Norm(wkh);
  

  // TODO VERY IMPORTANT HACK - start
  // if w1=w2=w3=0 it will result in matrix that is NaN, while the right
  // valie
  // according to the equation would be matrix I in this case. If the
  // input is
  // a zero vector - just return a one matrix I

  // if (v3IsZeroVector(angularVelocities) )
  // {
  // return diagonalOne;
  // }

  // TODO VERY IMPORTANT HACK - end

  // calculate sin and cos of each element of wkh for later use, use
  // sincos where possible

  /*
   *
   * It looks like the matrix form of these calculations gives a very big
   * error even in matlab
   *
   * for example:
   *
   * dt = 0.0 A(0.0) = I
   *
   * but in reality it gives us Z =
   *
   * 1.0000 0.0073 -0.0034 -0.0073 1.0000 -0.0046 0.0034 0.0046 1.0000
   *
   * which is quite far from the correct value of matrix diag(A)= [1 1 1]
   *
   *
   * so - try to replace it with more inline solution from the following
   * URL http://mathworld.wolfram.com/RodriguesRotationFormula.html
   *
   * RealMatrix matrixOne = new Array2DRowRealMatrix(diagonalOne);
   * RealMatrix matrixSwk = new
   * Array2DRowRealMatrix(angularVelocitiesSwk); RealMatrix
   * matrixSwkSquared = matrixSwk.copy().multiply(matrixSwk.copy());
   *
   *
   * RealMatrix result = new Array2DRowRealMatrix(kalmanAk); result =
   * matrixOne.subtract(
   * matrixSwk.scalarMultiply(Math.sin(wkhNorm)*(1d/wkhNorm)).add(
   * matrixSwkSquared.scalarMultiply( (1d -
   * Math_cos(wkhNorm))*(1d/wkhNormSquared) ) ) );
   */

  /*
   * The following S in matlab 0 -0.00763581553474069 0.00610865233466029
   * 0.00763581553474069 0 0.0106901414692402 -0.00610865233466029
   * -0.0106901414692402 0
   *
   * produces the following R 1 1.10627245596999e-20 -8.85017951283045e-21
   * -1.10627245596999e-20 1 -1.54878139787915e-20 8.85017951283045e-21
   * 1.54878139787915e-20 1
   */

  double wx = angularVelocities.x[0];
  double wy = angularVelocities.x[1];
  double wz = angularVelocities.x[2];
  double cf = Math_cos(wkhDtNorm);
  double sf = Math_sin(wkhDtNorm);

  // column 1
  // kalmanA[0][0] = cf + wx*wx*(1-cf);
  // kalmanA[1][0] = wz*sf + wx*wy*(1-cf);
  // kalmanA[2][0] = -wy*sf + wx*wz*(1-cf);
  //
  // //column 2
  // kalmanA[0][1] = wx * wy * (1 - cf) - wz * sf;
  // kalmanA[1][1] = cf + wy * wy * (1 - cf);
  // kalmanA[2][1] = wx * sf + wy * wz * (1 - cf);
  //
  // //column 3
  // kalmanA[0][2]= wy*sf + wx*wy*(1-cf);
  // kalmanA[1][2] =-wx*sf+wy*wz*(1-cf);
  // kalmanA[2][2] =cf+wz*wz*(1-cf);

  /*
   * w =
   *
   * [ 0, -wz, wy] [ wz, 0, -wx] [ -wy, wx, 0]
   *
   * according to the redrigues formula
   *
   * e^(w*tet) =
   *
   * [ (wy^2 + wz^2)*(cos(tet) - 1) + 1, - wz*sin(tet) - wx*wy*(cos(tet) -
   * 1), wy*sin(tet) - wx*wz*(cos(tet) - 1)] [ wz*sin(tet) -
   * wx*wy*(cos(tet) - 1), (wx^2 + wz^2)*(cos(tet) - 1) + 1, - wx*sin(tet)
   * - wy*wz*(cos(tet) - 1)] [ - wy*sin(tet) - wx*wz*(cos(tet) - 1),
   * wx*sin(tet) - wy*wz*(cos(tet) - 1), (wx^2 + wy^2)*(cos(tet) - 1) + 1]
   */

  // first row
  m.x[0].x[0] = (wy * wy + wz * wz) * (cf - 1) + 1;
  m.x[0].x[1] = -wz * sf - wx * wy * (cf - 1);
  m.x[0].x[2] = wy * sf - wx * wz * (cf - 1);

  // second row
  m.x[1].x[0] = wz * sf - wx * wy * (cf - 1);
  m.x[1].x[1] = (wx * wx + wz * wz) * (cf - 1) + 1;
  m.x[1].x[2] = -wx * sf - wy * wz * (cf - 1);

  m.x[2].x[0] = -wy * sf - wx * wz * (cf - 1);
  m.x[2].x[1] = wx * sf - wy * wz * (cf - 1);
  m.x[2].x[2] = (wx * wx + wy * wy) * (cf - 1) + 1;

  // TODO - verify that it will be treated correctly as an output,
  // and columns will NOT be mixed with lines!!!

  return m;

}

int currentSigmaValue(KF * pkf) {

  // if (skipperCounter > SKIP_AFTER ) {
  // return 0;
  // }
  // else {
  // skipperCounter++;
  // }

  //TODO - reenable the switching architecture
  //return 0;
  return (pkf->idleTimeLength >= pkf->sigmaWindowLengthSeconds) ? 1 : 0;

}


/**
 * returns current estimate of the covariance P ( that is being propagated through time)
 * @return
 */

MV3* getP(KF* pkf) {
  return &(pkf->currentCovariancePkEstimate);
}


MV2* getCurrentRollPitchRadians(KF* pkf) {
  return &(pkf->currentRollPitchEstimate);
}

/*
 * (non-Javadoc)
 *
 * @see
 * com.opopov.boat.compass.utils.MeasurementProcessor#getCurrentXEstimate()
 */

MV3* getCurrentXEstimate(KF* pkf) {
  return &(pkf->currentXEstimate);
}


MV3* getCurrentZEstimate(KF* pkf) {
  return &(pkf->currentZEstimate);
}


void initLibrary() {
  EPSILON_VECTOR.x[0] = 0.0000001;
  EPSILON_VECTOR.x[1] = 0.0000001;
  EPSILON_VECTOR.x[2] = 0.0000001;

  INITIAL_Z_ESTIMATE.x[0] = 0.01;
  INITIAL_Z_ESTIMATE.x[1] = 0.01;
  INITIAL_Z_ESTIMATE.x[2] = 0.98;


  //const long SKIP_AFTER = 1000;

  /**
   * diagonal matrix constant
   */
  memset(&diagonalOne, 0, sizeof(diagonalOne));
  diagonalOne.x[0].x[0] = 1;
  diagonalOne.x[1].x[1] = 1;
  diagonalOne.x[2].x[2] = 1;


}


void initialize(KF* pkf) {

}

/**
 * Multiplies 3x3 matrix a by 3 elements vector v as follows ma00*v0 +
 * ma01*v1 + ma02*v2 ma10*v0 + ma11*v1 + ma12*v2 ma20*v0 + ma21*v1 + ma22*v2
 */
void m3MultiplyByv3(   MMAT3 ma,   MV3 v,   MV3 * __restrict pOutV) {
  pOutV->x[0] = ma.x[0].x[0] * v.x[0] + ma.x[0].x[1] * v.x[1] + ma.x[0].x[2] * v.x[2];
  pOutV->x[1] = ma.x[1].x[0] * v.x[0] + ma.x[1].x[1] * v.x[1] + ma.x[1].x[2] * v.x[2];
  pOutV->x[2] = ma.x[2].x[0] * v.x[0] + ma.x[2].x[1] * v.x[1] + ma.x[2].x[2] * v.x[2];

}

/**
 * Multiplies transposed A 3x3 matrix a by 3 elements vector v as follows
 *
 * ma00*v0 + ma10*v1 + ma20*v2 ma01*v0 + ma11*v1 + ma21*v2 ma02*v0 + ma12*v1
 * + ma22*v2
 */
void m3TransposeMultiplyByv3(  MMAT3  ma ,    MV3 v,
                               MV3 * __restrict pOutV) {
  pOutV->x[0] = ma.x[0].x[0] * v.x[0] + ma.x[1].x[0] * v.x[1] + ma.x[2].x[0] * v.x[2];
  pOutV->x[1] = ma.x[0].x[1] * v.x[0] + ma.x[1].x[1] * v.x[1] + ma.x[2].x[1] * v.x[2];
  pOutV->x[2] = ma.x[0].x[2] * v.x[0] + ma.x[1].x[2] * v.x[1] + ma.x[2].x[2] * v.x[2];

}

/**
 * will apply it to the manetic vector directly - to interpolate between
 * measurements

 void multiplyTransposedTransitionMatrixByVector(MV3 v3InVector,
 MV3 v3OutVector) {
 m3TransposeMultiplyByv3(kalmanA, v3InVector, v3OutVector);

 }
*/


void rotateVectorByA(KF* pkf, MV3* pv, MV3* pvout) {

  //first calculate the norm of the input vector
  MV3 rv = *pv;
  MV3 ov;
  double n = v3Norm(rv);
  
  //if the norm is zero - just exit to avoid division by zero
  if (Math_abs(n) < 0.001) {
    return; 
  }


  //now normalize the vector
  v3DivideScalar(rv, n, &ov);
     
  
  //now rotate by A
  MMAT3 m = pkf->pma;
  m3MultiplyByv3(m, ov, &rv);
  

  //now multiply the vector back by its norm
  v3MultiplyScalar(rv, n, &ov);
     


  //return the result
  *pvout = ov; 
  

}


void updateCovarianceP(KF* pkf,double dt) {
  MV3 pkResult;

  MV3 pkResultBuffer;

  MV3 pkResultDivisionBuffer;

  MV3 pkResultSquared;

  MV3 vpk, vq,vr;



  v3AssignVector( &vpk , pkf->currentCovariancePkEstimate);
  v3AssignVector( &vq ,pkf->processNoiseCovarianceQ);
  v3AssignVector(&vr , pkf->measurementNoiseCovarianceR);


  // vq.mapMultiplyToSelf(dt); //Qt = Q*dt
  // vr.mapMultiplyToSelf(dt);

  v3EbeAdd(vpk, vq, &pkResult);

  // check if we need to take second part into account
  if (  currentSigmaValue(pkf) != 0) {
    v3EbeMultiply(vpk, vpk, &pkResultSquared); // p^2
    v3EbeAdd(vpk, vr, &pkResultBuffer); // p + r
    v3EbeDivide(pkResultSquared, pkResultBuffer, &pkResultDivisionBuffer); // p^2/(p+r)

    // substract
    v3EbeSubstract(pkResult, pkResultDivisionBuffer, &pkResultBuffer);
  }
  else {
    //else just take the current result 
    pkResultBuffer = pkResult;
  }

  v3AssignVector(&(pkf->currentCovariancePkEstimate), pkResultBuffer);

  // correct to avoid infinity when diving by zero
  ///if (v3IsZeroVector(pkf->currentCovariancePkEstimate)) {
  ///	v3MultiplyScalar(EPSILON_VECTOR, 2.0, &(pkf->currentCovariancePkEstimate));
  ///}



}

void updateEstimates(KF* pkf) {



  // now calculate x
  // TODO somehow check Z estimate for zero here

  // if any part of Z is zero at this step - just skip updating X at this
  // stage
  if (v3IsZeroNorm(pkf->currentZEstimate)) {
    v3AssignVector(&(pkf->currentZEstimate), pkf->lastGoodZEstimate);


  } else {
    // FIXME - this is not in the model but i will skip norming Z vector
    // to ||Z|| in
    // case where current sigma is 0. because if it is 0, then we just
    // integrating over w1 w2 w3 values, and no need to norm to adjust
    // ot the gravity vector
    // if (currentSigmaValue() == 0) {
    // currentXEstimate[0] = currentZEstimate[0];
    // currentXEstimate[1] = currentZEstimate[1];
    // currentXEstimate[2] = currentZEstimate[2];
    // }
    // else {

    // FIXME - not sure why, but i think that Z needs to be normed to
    // itself too
    // since it will grow disproportionaly otherwise

    // it probably should be normed if only some value of Y was added
    // to it, namely, when sigma !=0

    double zNorm = v3Norm(pkf->currentZEstimate);

    MV3 zResult;
    v3DivideScalar(pkf->currentZEstimate, zNorm, &zResult);
			
    //no norming in this case - nothing has been added to Z at the previous step
			
    pkf->currentZEstimate =  zResult;
    pkf->currentXEstimate =  zResult;
                        
                                
    //v3AssignVector(&(pkf->currentXEstimate), pkf->currentZEstimate);

    /*
      if (Double_isNaN(currentXEstimate[0])
      || Double_isNaN(currentXEstimate[1])
      || Double_isNaN(currentXEstimate[2])) {
      throw new IllegalStateException(
      "X estimate is NaN: currentValue ZEstimate "
      + Arrays.toString(currentZEstimate)
      + " z norm="
      + zNorm
      + " current covariance is :"
      + Arrays.toString(currentCovariancePkEstimate)
      + " " + " lastInputAcceleration "
      + Arrays.toString(lastInputAcceleration)
      + " lastInputAngularVelocity "
      + Arrays.toString(lastInputAngularVelocity)
      + " ma " + ma.toString() + " mz "


      );

      }

    */

    v3AssignVector( &(pkf->lastGoodZEstimate), pkf->currentZEstimate);

    // }
  }

  double roll = Math_atan(pkf->currentXEstimate.x[1] / pkf->currentXEstimate.x[2]);
  double pitch = Math_acos(pkf->currentXEstimate.x[2] / Math_cos(roll));

  // todo - fixme, use acos and asin below
  pkf->currentRollPitchEstimate.x[0] = roll;
  pkf->currentRollPitchEstimate.x[1] = pitch;

}


/*
 * (non-Javadoc)
 *
 * @see
 * com.opopov.boat.compass.utils.MeasurementProcessor#processMeasurement
 * (double, double, double, double, double, double, double)
 */

void processMeasurement(KF* pkf, double a1, double a2, double a3, double w1,
			double w2, double w3, double dt) {

  // VERIFY THE CORRECTNESS OF ARGUMETS
  //math c will fail to calc sin and cos 
  //if the value of w1/w2/w3 is greater than pi/2 or less than -pi/2
  //here we are going to make an assumption that if it is greater than PI/2
  // than the angular velocity is just PI/2
            
  //but for now, for the simplicity - we'll just skip updating anything if
  //w1/w2/w3 are outsise the PI/2 range

  double angularVelocityLimit = (3.14159265358979323846 / 2.0) - 0.0001;
  if (
      Math_abs(w1) > angularVelocityLimit ||
      Math_abs(w2) > angularVelocityLimit ||
      Math_abs(w3) > angularVelocityLimit  ) {
    //just skip everything
    return;

  }


  pkf->wVec.x[0] = w1;
  pkf->wVec.x[1] = w2;
  pkf->wVec.x[2] = w3;

  pkf->aVec.x[0] = a1;
  pkf->aVec.x[1] = a2;
  pkf->aVec.x[2] = a3;

  double aNormOriginal = v3Norm(pkf->aVec);
  updateSigmaValues(pkf, aNormOriginal, dt);

  pkf->lastInputAcceleration.x[0] = a1;
  pkf->lastInputAcceleration.x[1] = a2;
  pkf->lastInputAcceleration.x[2] = a3;

  pkf->lastInputAngularVelocity.x[0] = w1;
  pkf->lastInputAngularVelocity.x[1] = w2;
  pkf->lastInputAngularVelocity.x[2] = w3;

  // in our model - y = y/g u=u/g where u is acceleration

  v3DivideScalar(pkf->aVec, GRAVITATION_ACCELERATION_CONSTANT,
                 &(pkf->aVecGravityNormalized));

  // FIXME HACK - start
  // it may be weird, but in cases where sigma == 0 (too much acceleration
  // to use the accelerometer as
  // a gravity sensor AND when angular acceleration is zero w1=w2=w3=0,
  // AND
  // we have not had any viable value of Z before - we should just skip
  // this measurement
  // completely AND just wait until we can initialize the system properly
  // if (currentSigmaValue() == 0 && w1 == 0.0d && w2 == 0.0d && w3 ==
  // 0.0d &&
  // currentZEstimate[0] == 0.0d && currentZEstimate[1] == 0.0d &&
  // currentZEstimate[2] == 0.0d) {
  // System.out.println("Waiting to initialize properly...");
  // return;
  // }

  // FIXME HACK - end

  MMAT3 kalmanA = createKalmanA(pkf, pkf->wVec, dt);

  v3AssignMatrix(&(pkf->pma), kalmanA);


  MV3 * pvz = &(pkf->currentZEstimate);;

  // check if we need to calculate the second part of the equation in
  // brackets
  // since if sigma is zero - then the whole second part is zero
  // meaning that we do not take the measurement of the acceleration
  // at this step into account at all
  MV3 secondPart;
  MV3 zEstimateBuffer;
  if (currentSigmaValue(pkf) != 0) {
    MV3 * pvy = &(pkf->aVecGravityNormalized);
    MV3 * pvpk = &(pkf->currentCovariancePkEstimate);
    MV3 * pvr = &(pkf->measurementNoiseCovarianceR);
    // secondPart = vpk.ebeDivide(
    // vpk.copy().add(vr)
    // ).ebeMultiply(vy.subtract(currentZEstimate)
    // );
    // vr.mapMultiplyToSelf(dt);
    // (a/b)*c

    //in order to avoid gc - have a buffer

    MV3 measurementResultBuffer;
    MV3 processMeasurementKalmanGainResultBuffer;

    MV3 processMeasurementKalmanGain;


    v3EbeSubstract(*pvy, *pvz, &measurementResultBuffer); // yk - zestk ,
    // aliased as c
    v3EbeAdd(*pvpk, *pvr, &processMeasurementKalmanGainResultBuffer);// pk
    // +r,
    // aliased
    // as b


    v3EbeDivide(*pvpk, processMeasurementKalmanGainResultBuffer,
                &processMeasurementKalmanGain);
    // kalman gain is now ready and stored in the
    // processMeasurementKalmanGain variable
    // multiply the difference between current and previous measument by
    // kalman gain
    v3EbeMultiply(measurementResultBuffer,
                  processMeasurementKalmanGain, &secondPart);

    // zest + KalmanGain*(y-zest), where KalmanGain = pk/(pk+r)
    v3EbeAdd(pkf->currentZEstimate, secondPart, &zEstimateBuffer);

    // update z itself
    pkf->currentZEstimate =  zEstimateBuffer;
  }

  // A*[zestk + sigma*KalmanGain*(yk - zestk)]
  m3MultiplyByv3(kalmanA, pkf->currentZEstimate, &zEstimateBuffer);


  // update z itself
  pkf->currentZEstimate =  zEstimateBuffer;

  updateEstimates(pkf);

  updateCovarianceP(pkf,dt);


}



