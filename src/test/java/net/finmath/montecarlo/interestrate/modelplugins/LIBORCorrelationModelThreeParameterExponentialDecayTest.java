package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.time.TimeDiscretization;
import org.junit.Test;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class LIBORCorrelationModelThreeParameterExponentialDecayTest {

    @Test
    public void testThatCorrelationMatrixIsUnchangedWhenNumberOfFactorsEqualsNumberOfTimeSteps() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelThreeParameterExponentialDecay model =
                new LIBORCorrelationModelThreeParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, timeSteps, 0.1, 0.3, 0.2, false
                );

        double[][] originalCorrelationMatrix = model.generateCorrelationMatrix();

        //Check if correlation is the same
        boolean matricesDiffer = false;
        for (int row = 0; row < originalCorrelationMatrix.length; row++) {
            for (int col = 0; col < originalCorrelationMatrix[row].length; col++) {
                matricesDiffer |= Math.abs(originalCorrelationMatrix[row][col] / model.getCorrelation(0, row, col) - 1.0) > 0.1;
            }
        }

        assertFalse(matricesDiffer);
    }

    //very inaccurate test - might fail due to small changes
    @Test
    public void testThatCorrelationMatrixIsSimilarAfterFactorReduction() {
        double dt	= 0.5;
        double end	= 20.0;
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (end / dt), dt);
        int timeSteps = timeDiscretization.getNumberOfTimeSteps();

        LIBORCorrelationModelThreeParameterExponentialDecay model =
                new LIBORCorrelationModelThreeParameterExponentialDecay(
                        timeDiscretization,	timeDiscretization, timeSteps - 4, 0.1, 0.3, 0.2, false
                );

        double[][] originalCorrelationMatrix = model.generateCorrelationMatrix();

        //Check if correlation is the same
        int differentEntries = 0;
        for (int row = 0; row < originalCorrelationMatrix.length; row++) {
            for (int col = 0; col < originalCorrelationMatrix[row].length; col++) {
                if (Math.abs(originalCorrelationMatrix[row][col] / model.getCorrelation(0, row, col) - 1.0) > 0.1) {
                    differentEntries++;
                }
            }
        }

        assertTrue(differentEntries < 0.1 * timeSteps * timeSteps);
    }
}