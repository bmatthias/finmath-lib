package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.MatrixUtils;
import net.finmath.time.TimeDiscretization;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class LIBORCorrelationModelTwoCurveTenParameterExponentialDecayTest {
    @Test
    public void testThatCorrelationMatrixIsUnchangedWhenNumberOfFactorsEqualsNumberOfTimeSteps() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelTwoCurveTenParameterExponentialDecay model =
                new LIBORCorrelationModelTwoCurveTenParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, 2 * timeSteps, new double[]{ 0.1, 0.2, 0.1 },
                        new double[]{ 0.1, 0.1, 0.1, 0.1, 0.7, 1.0, 0.3, 0.0 }, false
                );

        double[][] originalCorrelationMatrix = model.generateCorrelationMatrix();
        model.initialize();

        assertTrue(MatrixUtils.frobeniusNorm(originalCorrelationMatrix, model.correlationMatrix) < 1e-5);
    }

    //very inaccurate test - might fail due to small changes
    @Test
    public void testThatCorrelationMatrixIsSimilarAfterFactorReduction() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelTwoCurveTenParameterExponentialDecay model =
                new LIBORCorrelationModelTwoCurveTenParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, 2 * timeSteps - 8, new double[]{ 0.1, 0.2, 0.1 },
                        new double[]{ 0.1, 0.1, 0.1, 0.1, 0.7, 1.0, 0.3 }, false
                );

        double[][] originalCorrelationMatrix = model.generateCorrelationMatrix();
        model.initialize();

        //Check if correlation is similar after factor reduction
        assertTrue(MatrixUtils.frobeniusNorm(originalCorrelationMatrix, model.correlationMatrix) < 1);
    }
}