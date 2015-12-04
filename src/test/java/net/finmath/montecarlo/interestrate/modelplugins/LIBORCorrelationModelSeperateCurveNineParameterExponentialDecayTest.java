package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.MatrixUtils;
import net.finmath.functions.LinearAlgebra;
import net.finmath.time.TimeDiscretization;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class LIBORCorrelationModelSeperateCurveNineParameterExponentialDecayTest {

    //TODO: Eigenvalues change slightly, so matrices are not equal - check if factorization can be improved.
    @Test
    public void testThatCorrelationMatrixIsUnchangedWhenNumberOfFactorsEqualsNumberOfTimeSteps() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelSeperateCurveNineParameterExponentialDecay model =
                new LIBORCorrelationModelSeperateCurveNineParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, 2 * timeSteps, new double[]{ 0.1, 0.5, 0.1 },
                        new double[]{ 0.1, 0.1, 0.1, 0.1, 0.7, 1.0, 0.3 }, true
                );

        double[][] firstCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] secondCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] curvesCorrelationMatrix = new double[timeSteps][timeSteps];

        model.generateCorrelationMatrices(timeSteps, firstCurveCorrelationMatrix, secondCurveCorrelationMatrix, curvesCorrelationMatrix);

        RealMatrix originalCorrelationMatrix = new Array2DRowRealMatrix(2 * timeSteps, 2 * timeSteps);
        originalCorrelationMatrix.setSubMatrix(firstCurveCorrelationMatrix, 0, 0);
        originalCorrelationMatrix.setSubMatrix(secondCurveCorrelationMatrix, timeSteps, timeSteps);
        originalCorrelationMatrix.setSubMatrix(curvesCorrelationMatrix, 0, timeSteps);
        originalCorrelationMatrix.setSubMatrix(new Array2DRowRealMatrix(curvesCorrelationMatrix).transpose().getData(), timeSteps, 0);

        model.initialize();

        assertTrue(MatrixUtils.matricesAreEqual(model.correlationMatrix, originalCorrelationMatrix.getData(), 0.1, 0));
    }

    //TODO: Not a good test. Matrix has negative eigenvalues. The model should not allow this.
    @Test
    public void testThatCorrelationMatrixIsUnchangedWhenNumberOfFactorsEqualsNumberOfTimeStepsAndFirstBlockIsSingular() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelSeperateCurveNineParameterExponentialDecay model =
                new LIBORCorrelationModelSeperateCurveNineParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, 2 * timeSteps, new double[]{ 0.0, 1.0, 0.1 },
                        new double[]{ 0.1, 0.1, 0.1, 0.1, 0.7, 1.0, 0.3 }, true
                );

        double[][] firstCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] secondCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] curvesCorrelationMatrix = new double[timeSteps][timeSteps];

        model.generateCorrelationMatrices(timeSteps, firstCurveCorrelationMatrix, secondCurveCorrelationMatrix, curvesCorrelationMatrix);

        RealMatrix originalCorrelationMatrix = new Array2DRowRealMatrix(2 * timeSteps, 2 * timeSteps);
        originalCorrelationMatrix.setSubMatrix(firstCurveCorrelationMatrix, 0, 0);
        originalCorrelationMatrix.setSubMatrix(secondCurveCorrelationMatrix, timeSteps, timeSteps);
        originalCorrelationMatrix.setSubMatrix(curvesCorrelationMatrix, 0, timeSteps);
        originalCorrelationMatrix.setSubMatrix(new Array2DRowRealMatrix(curvesCorrelationMatrix).transpose().getData(), timeSteps, 0);

        model.initialize();

        assertTrue(MatrixUtils.matricesAreEqual(model.correlationMatrix, originalCorrelationMatrix.getData(), 0.1, 0));
    }

    @Test
    public void testThatEigenvalueSumIsEqualAfterFactorReduction() {
        double dt	= 1;
        double end	= 4;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelSeperateCurveNineParameterExponentialDecay model =
                new LIBORCorrelationModelSeperateCurveNineParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, timeSteps, new double[]{ 0.1, 0.5, 0.1 },
                        new double[]{ 0.1, 0.1, 0.1, 0.1, 0.7, 1.0, 0.3 }, true
                );

        double[][] firstCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] secondCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] curvesCorrelationMatrix = new double[timeSteps][timeSteps];

        model.generateCorrelationMatrices(timeSteps, firstCurveCorrelationMatrix, secondCurveCorrelationMatrix, curvesCorrelationMatrix);

        RealMatrix originalCorrelationMatrix = new Array2DRowRealMatrix(2 * timeSteps, 2 * timeSteps);
        originalCorrelationMatrix.setSubMatrix(firstCurveCorrelationMatrix, 0, 0);
        originalCorrelationMatrix.setSubMatrix(secondCurveCorrelationMatrix, timeSteps, timeSteps);
        originalCorrelationMatrix.setSubMatrix(curvesCorrelationMatrix, 0, timeSteps);
        originalCorrelationMatrix.setSubMatrix(new Array2DRowRealMatrix(curvesCorrelationMatrix).transpose().getData(), timeSteps, 0);

        model.initialize();

        DoubleMatrix eigenvalues1 = Eigen.symmetricEigenvalues(new DoubleMatrix(originalCorrelationMatrix.getData()));
        DoubleMatrix eigenvalues2 = Eigen.symmetricEigenvalues(new DoubleMatrix(model.correlationMatrix));

        //Check if correlation is similar after factor reduction
        double eigenvalueSum1 = 0.0;
        double eigenvalueSum2 = 0.0;
        for (int i = 0; i < eigenvalues1.length; i++) {
            eigenvalueSum1 += eigenvalues1.get(i);
            eigenvalueSum2 += eigenvalues2.get(i);
        }

        assertTrue(Math.abs(eigenvalueSum1 / eigenvalueSum2 - 1.0) < 0.01);
    }
}