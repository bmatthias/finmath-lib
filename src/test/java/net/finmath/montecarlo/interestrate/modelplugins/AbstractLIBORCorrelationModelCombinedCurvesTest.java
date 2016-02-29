package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.MatrixUtils;
import net.finmath.time.TimeDiscretization;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.junit.Before;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class AbstractLIBORCorrelationModelCombinedCurvesTest {

    private AbstractLIBORCorrelationModelCombinedCurves model;
    private RealMatrix factorMatrix;
    private int numberOfFactors = 2;

    @Before
    public void setUp() throws Exception {
        double sqrt = 1.0 / Math.sqrt(2.0);

        //generates a correlation matrix with eigenvalues 0, 1 and 2
        factorMatrix = new Array2DRowRealMatrix(new double[][] {
                { sqrt, sqrt }, { 1.0, 0.0 }, { 0.0, 1.0 }
        });

        TimeDiscretization timeDiscretization = new TimeDiscretization(1.0, 1.5, 2.0, 3.0);
        model = new AbstractLIBORCorrelationModelCombinedCurves(
                timeDiscretization, timeDiscretization, numberOfFactors, null, null, false
        ) {
            @Override
            double[][] generateCorrelationMatrix() {
                RealMatrix factorMatrix = AbstractLIBORCorrelationModelCombinedCurvesTest.this.factorMatrix;
                return factorMatrix.multiply(factorMatrix.transpose()).getData();
            }

            @Override
            void adjustParameters() {}
        };

        model.initialize();
    }

    @Test
    public void testThatEigenvaluesAreCorrectAfterFactorization() {
        double[][] modelFactorMatrix = new Array2DRowRealMatrix(model.factorMatrix).getData();

        //columns of factor matrix are normed and multiplied with square root of eigenvalues
        List<Double> eigenvalues = new ArrayList<>();
        for (int i = 0; i < numberOfFactors; i++) {
            double eigenvalue = 0.0;
            for (int j = 0; j < 3; j++) {
                eigenvalue += modelFactorMatrix[j][i] * modelFactorMatrix[j][i];
            }
            eigenvalues.add(eigenvalue);
        }
        Collections.sort(eigenvalues);

        assertTrue(Math.abs(eigenvalues.get(0) - 1.0) < 1e-14);
        assertTrue(Math.abs(eigenvalues.get(1) - 2.0) < 1e-14);
    }

    @Test
    public void testThatFactorMatrixGeneratesTheSameCorrelationMatrix() {
        double[][] correlationMatrix = this.factorMatrix.multiply(this.factorMatrix.transpose()).getData();
        RealMatrix factorMatrix = new Array2DRowRealMatrix(model.factorMatrix);
        double[][] modelCorrelationMatrix = factorMatrix.multiply(factorMatrix.transpose()).getData();

        assertTrue(MatrixUtils.frobeniusNorm(correlationMatrix, modelCorrelationMatrix) < 1e-10);
    }

    @Test
    public void testThatCorrelationMatrixIsCorrectAfterFactorization() {
        double[][] correlationMatrix = factorMatrix.multiply(factorMatrix.transpose()).getData();

        assertTrue(MatrixUtils.frobeniusNorm(correlationMatrix, model.correlationMatrix) < 1e-10);
    }
}