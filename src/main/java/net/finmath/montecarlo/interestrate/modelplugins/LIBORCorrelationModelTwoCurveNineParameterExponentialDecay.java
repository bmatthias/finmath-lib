/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 20.05.2006
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.time.TimeDiscretizationInterface;

/**
 * Simple correlation model given by R, where R is a factor reduced matrix
 * (see {@link net.finmath.functions.LinearAlgebra#factorReduction(double[][], int)}) created from the
 * \( n \) Eigenvectors of \( \tilde{R} \) belonging to the \( n \) largest non-negative Eigenvalues,
 * where \( \tilde{R} = \tilde{\rho}_{i,j} \) and
 * \[ \tilde{\rho}_{i,j} = b + (1-b) * \exp(-a |T_{i} - T_{j}| - c \max(T_{i},T_{j}))
 *
 * @see net.finmath.functions.LinearAlgebra#factorReduction(double[][], int)
 *
 * @author Christian Fries
 */
public class LIBORCorrelationModelTwoCurveNineParameterExponentialDecay extends AbstractLIBORCorrelationModelCombinedCurves {

    public LIBORCorrelationModelTwoCurveNineParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, boolean isCalibrateable) {
        this(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, new double[] {
                0.1, 0.1, 0.1, 0.1, 0.1, 0.1
        }, isCalibrateable);
    }

    public LIBORCorrelationModelTwoCurveNineParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, double[] calibrationParameter, boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, calibrationParameter, isCalibrateable);
    }

    @Override
    protected double[][] generateCorrelationMatrix() {
        int timeSteps = liborPeriodDiscretization.getNumberOfTimeSteps();
        double[][] correlationMatrix = new double[2 * timeSteps][2 * timeSteps];
        for (int row = 0; row < timeSteps; row++) {
            for (int col = row; col < timeSteps; col++) {
                // Exponentially decreasing instantaneous correlation
                double T1 = liborPeriodDiscretization.getTime(row);
                double T2 = liborPeriodDiscretization.getTime(col);

                double correlation1 = fixedParameter[1] + (1.0 - fixedParameter[1]) * Math.exp(-fixedParameter[0] * Math.abs(T1 - T2) - fixedParameter[2] * Math.max(T1, T2));
                double correlation2 = calibrationParameter[1] + (1.0 - calibrationParameter[1]) * Math.exp(-calibrationParameter[0] * Math.abs(T1 - T2) - calibrationParameter[2] * Math.max(T1, T2));
                double correlation3 = calibrationParameter[3] * Math.exp(-calibrationParameter[4] * Math.abs(T1 - T2));
                double correlation4 = calibrationParameter[3] * Math.exp(-calibrationParameter[5] * Math.abs(T1 - T2));

                correlationMatrix[row][col] = correlation1;
                correlationMatrix[col][row] = correlation1;

                correlationMatrix[row + timeSteps][col + timeSteps] = correlation2;
                correlationMatrix[col + timeSteps][row + timeSteps] = correlation2;

                correlationMatrix[row + timeSteps][col] = correlation3;
                correlationMatrix[col][row + timeSteps] = correlation3;

                correlationMatrix[row][col + timeSteps] = correlation4;
                correlationMatrix[col + timeSteps][row] = correlation4;
            }
            correlationMatrix[row][row] = 1.0;
            correlationMatrix[row + timeSteps][row + timeSteps] = 1.0;
        }
        return correlationMatrix;
    }

    @Override
    void adjustParameters() {
        calibrationParameter[0] = Math.max(calibrationParameter[0], 0.0);
        calibrationParameter[1] = Math.min(Math.max(calibrationParameter[1], 0.0), 1.0);
        calibrationParameter[2] = Math.max(calibrationParameter[2], 0.0);
        calibrationParameter[3] = Math.min(Math.max(calibrationParameter[3], 0.0), 1.0);
        calibrationParameter[4] = Math.max(calibrationParameter[4], 0.0);
        calibrationParameter[5] = Math.max(calibrationParameter[5], 0.0);
    }
}
