package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.MatrixUtils;
import net.finmath.time.TimeDiscretization;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

@RunWith(Parameterized.class)
public class AbstractLIBORCorrelationModelSeperateCurvesTest {

    private AbstractLIBORCorrelationModelSeperateCurves model;
    private double[][] correlationMatrix;
    private int numberOfFactors;
    private int dim;

    @Parameterized.Parameters(name = "{index}: dim={0}, a={1}, b={2}, c={3}")
    public static Collection<Object[]> data() {
        return Arrays.asList(new Object[][]{
                {12, 0.3, 0.4, 0.1},
                {12, 0.3, 1.0, 0.1}, //creates a matrix with singular lower right block
                {12, 1.0, 0.3, 0.1}, //creates a matrix with singular upper left block
                {12, 1.0, 1.0, 0.33} //creates a singular matrix
        });
    }

    public AbstractLIBORCorrelationModelSeperateCurvesTest(int dim, double a, double b, double c) {
        this.dim = dim;
        this.correlationMatrix = new double[dim][];

        for (int i = 0; i < correlationMatrix.length; i++) {
            double[] row = new double[correlationMatrix.length];
            for (int j = 0; j < correlationMatrix.length; j++) {
                if (i >= correlationMatrix.length / 2 && j >= correlationMatrix.length / 2) {
                    row[j] = b + (1.0 - b) * Math.exp(-c * Math.abs(i - j));
                } else {
                    row[j] = a + (1.0 - a) * Math.exp(-c * Math.abs(i - j));
                }
            }
            correlationMatrix[i] = row;
        }

        //this makes the test matrix positive semi-definite,
        //because negative eigenvalues are dismissed during factorization
        init(dim);
        correlationMatrix = model.correlationMatrix;
    }

    @Test
    public void testThatEigenvalueSumIsEqualAfterFactorReduction() {
        init(4);

        DoubleMatrix eigenvalues1 = Eigen.symmetricEigenvalues(new DoubleMatrix(correlationMatrix));
        double[][] modelFactorMatrix = new Array2DRowRealMatrix(model.factorMatrix).getData();

        //columns of factor matrix should be normed and multiplied with square root of eigenvalues
        List<Double> eigenvalues2 = new ArrayList<>();
        for (int i = 0; i < numberOfFactors; i++) {
            double eigenvalue = 0.0;
            for (int j = 0; j < modelFactorMatrix.length; j++) {
                eigenvalue += modelFactorMatrix[j][i] * modelFactorMatrix[j][i];
            }
            eigenvalues2.add(eigenvalue);
        }

        double eigenvalueSum1 = 0.0;
        double eigenvalueSum2 = 0.0;
        for (int i = 0; i < eigenvalues1.length; i++) {
            eigenvalueSum1 += eigenvalues1.get(i);
            if(i < eigenvalues2.size()) eigenvalueSum2 += eigenvalues2.get(i);
        }

        assertTrue(Math.abs(eigenvalueSum1 / eigenvalueSum2 - 1.0) < 0.01);
    }

    @Test
    public void testThatEigenvaluesAreEqualWithoutFactorReduction() {
        init(dim);

        DoubleMatrix eigenvalues1 = Eigen.symmetricEigenvalues(new DoubleMatrix(correlationMatrix));
        DoubleMatrix eigenvalues2 = Eigen.symmetricEigenvalues(new DoubleMatrix(model.correlationMatrix));

        //columns of factor matrix should be normed and multiplied with square root of eigenvalues
        for (int i = 0; i < numberOfFactors; i++) {
            assertEquals(eigenvalues1.get(i), eigenvalues2.get(i), 1e-10);
        }
    }

    @Test
    public void testThatFactorMatrixGeneratesTheSameCorrelationMatrix() {
        init(dim);

        RealMatrix factorMatrix = new Array2DRowRealMatrix(model.factorMatrix);
        double[][] modelCorrelationMatrix = factorMatrix.multiply(factorMatrix.transpose()).getData();

        assertTrue(MatrixUtils.frobeniusNorm(modelCorrelationMatrix, model.correlationMatrix) < 1e-10);
        assertTrue(MatrixUtils.frobeniusNorm(modelCorrelationMatrix, correlationMatrix) < 1e-10);
    }

    @Test
    public void testThatCorrelationMatrixIsCorrectAfterFactorization() {
        init(dim);

        assertTrue(MatrixUtils.frobeniusNorm(model.correlationMatrix, correlationMatrix) < 1e-10);
    }

    private void init(int numberOfFactors) {
        this.numberOfFactors = numberOfFactors;

        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, correlationMatrix.length / 2, 1.0);
        model = new AbstractLIBORCorrelationModelSeperateCurves(
                timeDiscretization, timeDiscretization, numberOfFactors, null, null, false
        ) {
            @Override
            void generateCorrelationMatrices(int dim, double[][] firstCurveCorrelationMatrix, double[][] secondCurveCorrelationMatrix, double[][] curvesCorrelationMatrix) {
                for (int i = 0; i < dim; i++) {
                    for (int j = 0; j < dim; j++) {
                        firstCurveCorrelationMatrix[i][j] = AbstractLIBORCorrelationModelSeperateCurvesTest.this.correlationMatrix[i][j];
                        secondCurveCorrelationMatrix[i][j] = AbstractLIBORCorrelationModelSeperateCurvesTest.this.correlationMatrix[i + dim][j + dim];
                        curvesCorrelationMatrix[i][j] = AbstractLIBORCorrelationModelSeperateCurvesTest.this.correlationMatrix[i][j+ dim];
                    }
                }
            }

            @Override
            void adjustParameters() {}
        };

        model.initialize();
    }
}