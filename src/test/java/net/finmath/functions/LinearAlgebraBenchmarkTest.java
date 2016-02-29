package net.finmath.functions;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.EigenvalueDecomposition;
import net.finmath.MatrixUtils;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RRQRDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.random.MersenneTwister;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.commons.math3.random.SynchronizedRandomGenerator;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

import static org.junit.Assert.assertTrue;

/*
* Includes simple benchmark tests for the various linear algebra libraries. Note that the times required
* for copying the arrays (i.e. creating the matrix objects from the plain Java arrays) are included in the results.
*/
public class LinearAlgebraBenchmarkTest {

    @Test
    public void testMatrixMultiplication() throws Exception {
        long totalTimeMath = 0;
        long totalTimeJBlas = 0;
        long totalTimeTrivial = 0;

        RandomGenerator randomGenerator = new SynchronizedRandomGenerator(new MersenneTwister(3271));
        ExecutorService executorService = Executors.newFixedThreadPool(10);
        List<Future<Double[]>> results = new ArrayList<>();
        for(int dim1 = 10; dim1 <= 100; dim1 += 10) {
            final int finalDim = dim1;
            results.add(executorService.submit(() -> {
                Double[] result = new Double[] { 0.0, 0.0, 0.0 };

                int dim2 =  (finalDim * 8) / 10;
                double[][] matrix1 = new double[finalDim][dim2];
                double[][] matrix2 = new double[dim2][finalDim];

                for (int n = 0; n < 10000; n++) {
                    for (int i = 0; i < finalDim; i++) {
                        for (int j = 0; j < dim2; j++) {
                            matrix1[i][j] = randomGenerator.nextDouble();
                            matrix2[j][i] = randomGenerator.nextDouble();
                        }
                    }

                    long time = System.currentTimeMillis();
                    Array2DRowRealMatrix mathResult = new Array2DRowRealMatrix(matrix1).multiply(new Array2DRowRealMatrix(matrix2));
                    result[0] += System.currentTimeMillis() - time;

                    time = System.currentTimeMillis();
                    DoubleMatrix jblasResult = new DoubleMatrix(matrix1).mmul(new DoubleMatrix(matrix2));
                    result[1] += System.currentTimeMillis() - time;

                    time = System.currentTimeMillis();
                    double[][] javaResult = new double[finalDim][finalDim];
                    for (int k = 0; k < finalDim; k++) {
                        for (int l = 0; l < finalDim; l++) {
                            for (int m = 0; m < dim2; m++) {
                                javaResult[k][l] += matrix1[k][m] * matrix2[m][l];
                            }
                        }
                    }
                    result[2] += System.currentTimeMillis() - time;

                    //Check if the results are the same and in particular check for
                    //concurrency issues with JBlas that would disqualify it for use in multi threaded tasks.
                    assertTrue(MatrixUtils.frobeniusNorm(mathResult.getDataRef(), jblasResult.toArray2()) < 1e-10);
                    assertTrue(MatrixUtils.frobeniusNorm(mathResult.getDataRef(), javaResult) < 1e-10);
                }

                return result;
            }));
        }

        for (Future<Double[]> futureResult : results) {
            Double[] result = futureResult.get();

            totalTimeMath += result[0];
            totalTimeJBlas += result[1];
            totalTimeTrivial += result[2];
        }

        System.out.println("Math:\t" + totalTimeMath);
        System.out.println("JBlas:\t" + totalTimeJBlas);
        System.out.println("Java:\t" + totalTimeTrivial);
    }

    @Test
    public void testSymmetricEigenvalueDecomposition() throws Exception {
        long totalTimeMath = 0;
        long totalTimeJBlas = 0;
        long totalTimeColt = 0;

        RandomGenerator randomGenerator = new SynchronizedRandomGenerator(new MersenneTwister(3271));
        ExecutorService executorService = Executors.newFixedThreadPool(10);
        List<Future<Double[]>> results = new ArrayList<>();
        for(int dim = 10; dim <= 100; dim += 10) {
            final int finalDim = dim;
            results.add(executorService.submit(() -> {
                Double[] result = new Double[] { 0.0, 0.0, 0.0 };

                double[][] matrix = new double[finalDim][finalDim];
                for (int n = 0; n < 1000; n++) {
                    for (int i = 0; i < finalDim; i++) {
                        for (int j = 0; j < finalDim; j++) {
                            matrix[i][j] = matrix[j][i] = randomGenerator.nextDouble();
                        }
                    }

                    long time = System.currentTimeMillis();
                    new EigenDecomposition(new Array2DRowRealMatrix(matrix));
                    result[0] += System.currentTimeMillis() - time;

                    time = System.currentTimeMillis();
                    Eigen.symmetricEigenvalues(new DoubleMatrix(matrix));
                    result[1] += System.currentTimeMillis() - time;

                    time = System.currentTimeMillis();
                    new EigenvalueDecomposition(new DenseDoubleMatrix2D(matrix));
                    result[2] += System.currentTimeMillis() - time;
                }
                return result;
            }));
        }

        for (Future<Double[]> futureResult : results) {
            Double[] result = futureResult.get();

            totalTimeMath += result[0];
            totalTimeJBlas += result[1];
            totalTimeColt += result[2];
        }

        System.out.println("Math:\t" + totalTimeMath);
        System.out.println("JBlas:\t" + totalTimeJBlas);
        System.out.println("Colt:\t" + totalTimeColt);
    }

    @Test
    public void testFactorReduction() throws Exception {
        long totalTimeMath = 0;
        long totalTimeJBlas = 0;

        RandomGenerator randomGenerator = new SynchronizedRandomGenerator(new MersenneTwister(3271));
        ExecutorService executorService = Executors.newFixedThreadPool(10);
        List<Future<Double[]>> results = new ArrayList<>();
        for(int dim = 10; dim <= 80; dim += 10) {
            final int finalDim = dim;
            results.add(executorService.submit(() -> {
                Double[] result = new Double[] { 0.0, 0.0 };

                double[][] firstMatrix = new double[finalDim][finalDim];
                double[][] secondMatrix = new double[finalDim][finalDim];
                double[][] thirdMatrix = new double[finalDim][finalDim];
                for (int n = 0; n < 1000; n++) {
                    for (int i = 0; i < finalDim; i++) {
                        for (int j = 0; j < finalDim; j++) {
                            firstMatrix[i][j] = firstMatrix[j][i] = randomGenerator.nextDouble();
                            secondMatrix[i][j] = secondMatrix[j][i] = randomGenerator.nextDouble();
                            thirdMatrix[i][j] = thirdMatrix[j][i] = randomGenerator.nextDouble();
                        }
                        firstMatrix[i][i] = secondMatrix[i][i] = 1.0;
                    }

                    long time = System.currentTimeMillis();
                    getFactorMatrixFromBlockCorrelationMatrix(finalDim, firstMatrix, secondMatrix, thirdMatrix);
                    result[0] += System.currentTimeMillis() - time;

                    time = System.currentTimeMillis();
                    LinearAlgebra.getFactorMatrixFromBlockCorrelationMatrix(finalDim, firstMatrix, secondMatrix, thirdMatrix);
                    result[1] += System.currentTimeMillis() - time;
                }
                return result;
            }));
        }

        for (Future<Double[]> futureResult : results) {
            Double[] result = futureResult.get();

            totalTimeMath += result[0];
            totalTimeJBlas += result[1];
        }

        System.out.println("Math:\t" + totalTimeMath);
        System.out.println("JBlas:\t" + totalTimeJBlas);
    }

    public static RealMatrix getFactorMatrixFromBlockCorrelationMatrix(int numberOfFactors,
                                                                       double[][] firstCurveCorrelationMatrix,
                                                                       double[][] secondCurveCorrelationMatrix,
                                                                       double[][] curvesCorrelationMatrix) {
        int dimension = firstCurveCorrelationMatrix.length;
        RealMatrix factorMatrix;
        RRQRDecomposition rrqrDecomposition = new RRQRDecomposition(new Array2DRowRealMatrix(firstCurveCorrelationMatrix));
        if (rrqrDecomposition.getRank(1e-10) == dimension) {
            RealMatrix inverseOfFirstMatrix = rrqrDecomposition.getSolver().getInverse();
            RealMatrix C = new Array2DRowRealMatrix(curvesCorrelationMatrix);
            RealMatrix inverseATimesC = inverseOfFirstMatrix.multiply(C);
            Array2DRowRealMatrix B = (Array2DRowRealMatrix)new Array2DRowRealMatrix(secondCurveCorrelationMatrix).subtract(C.transpose().multiply(inverseATimesC));

            double[][] identityMatrix = MatrixUtils.createIdentityMatrix(dimension);
            Array2DRowRealMatrix N = new Array2DRowRealMatrix(2 * dimension, 2 * dimension);
            N.setSubMatrix(identityMatrix, 0, 0);
            N.setSubMatrix(identityMatrix, dimension, dimension);
            N.setSubMatrix(inverseATimesC instanceof Array2DRowRealMatrix ?
                            ((Array2DRowRealMatrix)inverseATimesC).getDataRef() : inverseATimesC.getData(), 0, dimension);

            /*
             * Perform a factor decomposition (and reduction if numberOfFactors < correlationMatrix.columns())
             */
            LinearAlgebra.fixSymmetry(B.getDataRef());

            double[][] firstFactorMatrix = LinearAlgebra.factorReduction(firstCurveCorrelationMatrix, (int) Math.ceil(0.5 * numberOfFactors));
            double[][] secondFactorMatrix = LinearAlgebra.factorReduction(B.getDataRef(), numberOfFactors / 2);

            RealMatrix M = new Array2DRowRealMatrix(2 * dimension, numberOfFactors);
            M.setSubMatrix(firstFactorMatrix, 0, 0);
            M.setSubMatrix(secondFactorMatrix, dimension, (int) Math.ceil(0.5 * numberOfFactors));

            double[][] fM = N.transpose().multiply(M).getData();
            LinearAlgebra.normalizeRows(fM);

            double[][] reducedCorrelationMatrix = new DoubleMatrix(fM).mmul(new DoubleMatrix(fM).transpose()).toArray2();
            factorMatrix = new Array2DRowRealMatrix(LinearAlgebra.getFactorMatrix(reducedCorrelationMatrix, numberOfFactors));
        } else {
            rrqrDecomposition = new RRQRDecomposition(new Array2DRowRealMatrix(secondCurveCorrelationMatrix));
            if (rrqrDecomposition.getRank(1e-10) == dimension) {
                RealMatrix inverseOfSecondMatrix = rrqrDecomposition.getSolver().getInverse();
                RealMatrix C = new Array2DRowRealMatrix(curvesCorrelationMatrix);
                RealMatrix cTimesInverseD = C.multiply(inverseOfSecondMatrix);
                Array2DRowRealMatrix B = (Array2DRowRealMatrix)new Array2DRowRealMatrix(firstCurveCorrelationMatrix).subtract(cTimesInverseD.multiply(C.transpose()));

                double[][] identityMatrix = MatrixUtils.createIdentityMatrix(dimension);
                Array2DRowRealMatrix N = new Array2DRowRealMatrix(2 * dimension, 2 * dimension);
                N.setSubMatrix(identityMatrix, 0, 0);
                N.setSubMatrix(identityMatrix, dimension, dimension);
                N.setSubMatrix(cTimesInverseD instanceof Array2DRowRealMatrix ?
                        ((Array2DRowRealMatrix)cTimesInverseD).getDataRef() : cTimesInverseD.getData(), 0, dimension);

                /*
                 * Perform a factor decomposition (and reduction if numberOfFactors < correlationMatrix.columns())
                 */
                LinearAlgebra.fixSymmetry(B.getDataRef());

                double[][] firstFactorMatrix = LinearAlgebra.factorReduction(B.getDataRef(), (int) Math.ceil(0.5 * numberOfFactors));
                double[][] secondFactorMatrix = LinearAlgebra.factorReduction(secondCurveCorrelationMatrix, numberOfFactors / 2);

                RealMatrix M = new Array2DRowRealMatrix(2 * dimension, numberOfFactors);
                M.setSubMatrix(firstFactorMatrix, 0, 0);
                M.setSubMatrix(secondFactorMatrix, dimension, (int) Math.ceil(0.5 * numberOfFactors));

                double[][] fM = N.multiply(M).getData();
                LinearAlgebra.normalizeRows(fM);

                double[][] reducedCorrelationMatrix = new DoubleMatrix(fM).mmul(new DoubleMatrix(fM).transpose()).toArray2();
                factorMatrix = new Array2DRowRealMatrix(LinearAlgebra.getFactorMatrix(reducedCorrelationMatrix, numberOfFactors));
            } else {
                Array2DRowRealMatrix correlationMatrix = new Array2DRowRealMatrix(2 * dimension, 2 * dimension);
                correlationMatrix.setSubMatrix(firstCurveCorrelationMatrix, 0, 0);
                correlationMatrix.setSubMatrix(secondCurveCorrelationMatrix, dimension, dimension);
                correlationMatrix.setSubMatrix(curvesCorrelationMatrix, 0, dimension);
                correlationMatrix.setSubMatrix(new Array2DRowRealMatrix(curvesCorrelationMatrix).transpose().getData(), dimension, 0);

                factorMatrix = new Array2DRowRealMatrix(LinearAlgebra.factorReduction(correlationMatrix.getDataRef(), numberOfFactors));
            }
        }
        return factorMatrix;
    }
}