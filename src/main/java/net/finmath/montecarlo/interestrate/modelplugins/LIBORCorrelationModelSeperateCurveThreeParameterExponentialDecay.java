/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 20.05.2006
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.MatrixUtils;
import net.finmath.functions.LinearAlgebra;
import net.finmath.time.TimeDiscretizationInterface;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.jblas.DoubleMatrix;

import java.util.stream.IntStream;

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
public class LIBORCorrelationModelSeperateCurveThreeParameterExponentialDecay extends AbstractLIBORCorrelationModelSeperateCurves {

    public LIBORCorrelationModelSeperateCurveThreeParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, boolean isCalibrateable) {
        this(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, new double[] { 0.1, 0.1, 0.1 }, isCalibrateable);
    }

    public LIBORCorrelationModelSeperateCurveThreeParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, double[] calibrationParameter, boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, calibrationParameter, isCalibrateable);
    }

    @Override
    synchronized void initialize() {
        int numberOfFactors = getNumberOfFactors();

        int timeSteps = liborPeriodDiscretization.getNumberOfTimeSteps();
        double[][] firstCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] secondCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] thirdCurveCorrelationMatrix = new double[timeSteps][timeSteps];

        generateCorrelationMatrices(timeSteps, firstCurveCorrelationMatrix, secondCurveCorrelationMatrix, thirdCurveCorrelationMatrix);

        double[][] firstFactorMatrix = LinearAlgebra.factorReduction(firstCurveCorrelationMatrix, numberOfFactors / 2);
        DoubleMatrix secondFactorMatrix = new DoubleMatrix(LinearAlgebra.factorReduction(secondCurveCorrelationMatrix, numberOfFactors / 2));

        double sqrt = Math.sqrt(1.0 - calibrationParameter[1] * calibrationParameter[1]);
        DoubleMatrix rho = DoubleMatrix.diag(new DoubleMatrix(IntStream.range(0, numberOfFactors / 2).mapToDouble((i) -> calibrationParameter[1]).toArray()));
        DoubleMatrix sqrtRho = DoubleMatrix.diag(new DoubleMatrix(IntStream.range(0, numberOfFactors / 2).mapToDouble((i) -> sqrt).toArray()));

        DoubleMatrix thirdFactorMatrix = secondFactorMatrix.mmul(rho);
        secondFactorMatrix = secondFactorMatrix.mmul(sqrtRho);

        DoubleMatrix M = new DoubleMatrix(2 * timeSteps, numberOfFactors);
        MatrixUtils.setSubMatrix(M, firstFactorMatrix, 0, 0);
        MatrixUtils.setSubMatrix(M, secondFactorMatrix.toArray2(), timeSteps, (int) Math.ceil(0.5 * numberOfFactors));
        MatrixUtils.setSubMatrix(M, thirdFactorMatrix.toArray2(), timeSteps, 0);

        this.factorMatrix = M.toArray2();
        this.correlationMatrix = new DoubleMatrix(this.factorMatrix).mmul(new DoubleMatrix(this.factorMatrix).transpose()).toArray2();
    }

    @Override
    void generateCorrelationMatrices(int dim, double[][] firstCurveCorrelationMatrix, double[][] secondCurveCorrelationMatrix, double[][] curvesCorrelationMatrix) {
        int timeSteps = getLiborPeriodDiscretization().getNumberOfTimeSteps();
        for (int row = 0; row < timeSteps; row++) {
            for (int col = row; col < timeSteps; col++) {
                // Exponentially decreasing instantaneous correlation
                double T1 = liborPeriodDiscretization.getTime(row);
                double T2 = liborPeriodDiscretization.getTime(col);

                double correlation1 = Math.exp(-fixedParameter[0] * Math.abs(T1 - T2));
                double correlation2 = Math.exp(-calibrationParameter[0] * Math.abs(T1 - T2));

                firstCurveCorrelationMatrix[row][col] = correlation1;
                firstCurveCorrelationMatrix[col][row] = correlation1;

                secondCurveCorrelationMatrix[row][col] = correlation2;
                secondCurveCorrelationMatrix[col][row] = correlation2;

            }
            firstCurveCorrelationMatrix[row][row] = 1.0;
            secondCurveCorrelationMatrix[row][row] = 1.0;
        }
    }

    @Override
    void adjustParameters() {
        calibrationParameter[0] = Math.max(calibrationParameter[0], 0.0);
        calibrationParameter[1] = Math.min(Math.max(calibrationParameter[1], 0.0), 1.0);
    }
}
