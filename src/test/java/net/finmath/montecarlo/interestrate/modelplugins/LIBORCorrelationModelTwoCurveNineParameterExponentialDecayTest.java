package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.time.TimeDiscretization;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class LIBORCorrelationModelTwoCurveNineParameterExponentialDecayTest {
    @Test
    public void testThatCorrelationMatrixIsUnchangedWhenNumberOfFactorsEqualsNumberOfTimeSteps() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelTwoCurveNineParameterExponentialDecay model =
                new LIBORCorrelationModelTwoCurveNineParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, 2 * timeSteps, new double[]{ 0.1, 0.2, 0.1 },
                        new double[]{ 0.1, 0.1, 0.1, 0.1, 0.7, 1.0, 0.3 }, false
                );

        double[][] originalCorrelationMatrix = model.generateCorrelationMatrix();

        //Check if correlation is the same
        boolean matricesAreEqual = true;
        for (int row = 0; row < originalCorrelationMatrix.length; row++) {
            for (int col = 0; col < originalCorrelationMatrix[row].length; col++) {
                matricesAreEqual &= Math.abs(originalCorrelationMatrix[row][col] / model.getCorrelation(0, row, col) - 1.0) < 0.1;
            }
        }

        assertTrue(matricesAreEqual);
    }

    //very inaccurate test - might fail due to small changes
    @Test
    public void testThatCorrelationMatrixIsSimilarAfterFactorReduction() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelTwoCurveNineParameterExponentialDecay model =
                new LIBORCorrelationModelTwoCurveNineParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, 2 * timeSteps - 8, new double[]{ 0.1, 0.2, 0.1 },
                        new double[]{ 0.1, 0.1, 0.1, 0.1, 0.7, 1.0, 0.3 }, false
                );

        double[][] originalCorrelationMatrix = model.generateCorrelationMatrix();

        //Check if correlation is similar after factor reduction
        int differentEntries = 0;
        for (int row = 0; row < originalCorrelationMatrix.length; row++) {
            for (int col = 0; col < originalCorrelationMatrix[row].length; col++) {
                if (Math.abs(originalCorrelationMatrix[row][col] / model.getCorrelation(0, row, col) - 1.0) > 0.1) {
                    differentEntries++;
                }
            }
        }

        assertTrue(differentEntries < 0.2 * 4 * timeSteps * timeSteps);
    }
}