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
public class LIBORCorrelationModelSeperateCurveFiveParameterExponentialDecay extends AbstractLIBORCorrelationModelSeperateCurves {

    public LIBORCorrelationModelSeperateCurveFiveParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, boolean isCalibrateable) {
        this(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, new double[] {
                0.1, 0.1, 0.1, 0.1, 0.1
        }, isCalibrateable);
    }

    public LIBORCorrelationModelSeperateCurveFiveParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, double[] calibrationParameter, boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, calibrationParameter, isCalibrateable);
    }

    @Override
    protected void generateCorrelationMatrices(int timeSteps, double[][] firstCurveCorrelationMatrix, double[][] secondCurveCorrelationMatrix, double[][] curvesCorrelationMatrix) {
        for (int row = 0; row < timeSteps; row++) {
            for (int col = row + 1; col < timeSteps; col++) {
                // Exponentially decreasing instantaneous correlation
                double T1 = liborPeriodDiscretization.getTime(row);
                double T2 = liborPeriodDiscretization.getTime(col);

                double correlation1 = Math.exp(-fixedParameter[0] * Math.abs(T1 - T2));
                double correlation2 = Math.exp(-calibrationParameter[0] * Math.abs(T1 - T2));
                double correlation3 = calibrationParameter[4] + (calibrationParameter[1] - calibrationParameter[4]) * Math.exp(-calibrationParameter[2] * Math.abs(T1 - T2));
                double correlation4 = calibrationParameter[4] + (calibrationParameter[1] - calibrationParameter[4]) * Math.exp(-calibrationParameter[3] * Math.abs(T1 - T2));

                firstCurveCorrelationMatrix[row][col] = correlation1;
                firstCurveCorrelationMatrix[col][row] = correlation1;

                secondCurveCorrelationMatrix[row][col] = correlation2;
                secondCurveCorrelationMatrix[col][row] = correlation2;

                curvesCorrelationMatrix[row][col] = correlation3;
                curvesCorrelationMatrix[col][row] = correlation4;
            }
            firstCurveCorrelationMatrix[row][row] = 1.0;
            secondCurveCorrelationMatrix[row][row] = 1.0;
        }
    }

    @Override
    void adjustParameters() {
        calibrationParameter[0] = Math.max(calibrationParameter[0], 0.0);
        calibrationParameter[1] = Math.min(Math.max(calibrationParameter[1], 0.0), 1.0);
        calibrationParameter[2] = Math.max(calibrationParameter[2], 0.0);
        calibrationParameter[3] = Math.min(Math.max(calibrationParameter[3], 0.0), 1.0);
        calibrationParameter[4] = Math.max(calibrationParameter[4], 0.0);
    }
}
