package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.functions.LinearAlgebra;
import net.finmath.time.TimeDiscretizationInterface;
import org.apache.commons.math3.linear.RealMatrix;
import org.jblas.DoubleMatrix;

public abstract class AbstractLIBORCorrelationModelSeperateCurves extends AbstractMultiCurveLIBORCorrelationModel {

    public AbstractLIBORCorrelationModelSeperateCurves(TimeDiscretizationInterface timeDiscretization,
                                                       TimeDiscretizationInterface liborPeriodDiscretization,
                                                       int numberOfFactors,
                                                       double[] fixedParameter,
                                                       double[] calibrationParameter,
                                                       boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, calibrationParameter, isCalibrateable);
    }

    @Override
    synchronized void initialize() {

        int timeSteps = liborPeriodDiscretization.getNumberOfTimeSteps();
        double[][] firstCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] secondCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] curvesCorrelationMatrix = new double[timeSteps][timeSteps];

        generateCorrelationMatrices(
                timeSteps, firstCurveCorrelationMatrix, secondCurveCorrelationMatrix, curvesCorrelationMatrix
        );

        DoubleMatrix factorMatrix = LinearAlgebra.getFactorMatrixFromBlockCorrelationMatrix(
                getNumberOfFactors(), firstCurveCorrelationMatrix, secondCurveCorrelationMatrix, curvesCorrelationMatrix
        );
        this.factorMatrix = factorMatrix.toArray2();
        this.correlationMatrix = factorMatrix.mmul(factorMatrix.transpose()).toArray2();
    }

    abstract void generateCorrelationMatrices(
            int dim, double[][] firstCurveCorrelationMatrix,
            double[][] secondCurveCorrelationMatrix, double[][] curvesCorrelationMatrix
    );
}
