package net.finmath;

import org.jblas.DoubleMatrix;

public class MatrixUtils {
    public static boolean checkSymmetry(double[][] matrix) {
        for (int row = 0; row < matrix.length; row++) {
            for (int col = 0; col < matrix[row].length; col++) {
                if(Math.abs(matrix[row][col] - matrix[col][row]) > 1e-14) return false;
            }
        }
        return true;
    }

    public static double frobeniusNorm(double[][] A, double[][] B) {
        return new DoubleMatrix(A).distance2(new DoubleMatrix(B));
    }

    public static String printMatrix(double[][] matrix) {
        StringBuilder stringBuilder = new StringBuilder();
        for (double[] row : matrix) {
            for (double value : row) {
                stringBuilder.append(value).append("\t");
            }
            stringBuilder.append("\n");
        }
        return stringBuilder.toString();
    }

    public static double[][] createIdentityMatrix(int dimension) {
        double[][] identityMatrix = new double[dimension][dimension];
        for(int i = 0; i < dimension; i++) {
            identityMatrix[i][i] = 1.0;
        }
        return identityMatrix;
    }

    public static void setSubMatrix(DoubleMatrix matrix, double[][] subMatrix, int fromRow, int fromColumn) {
        for(int row = 0; row < subMatrix.length; row++) {
            for(int column = 0; column < subMatrix[row].length; column++) {
                matrix.put(fromRow + row, fromColumn + column, subMatrix[row][column]);
            }
        }
    }
}
