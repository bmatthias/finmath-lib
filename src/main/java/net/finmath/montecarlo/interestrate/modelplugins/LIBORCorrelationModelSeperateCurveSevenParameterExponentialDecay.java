/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 20.05.2006
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.functions.LinearAlgebra;
import net.finmath.time.TimeDiscretizationInterface;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;

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
public class LIBORCorrelationModelSeperateCurveSevenParameterExponentialDecay extends LIBORCorrelationModel {

	private int		numberOfFactors;
	private double[] fixedParameter;
    private double[] calibrationParameter;
	private final boolean isCalibrateable;

    private transient double[][]	correlationMatrix;
	private transient double[][]	factorMatrix;

    public LIBORCorrelationModelSeperateCurveSevenParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, boolean isCalibrateable) {
        this(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, new double[] {
                0.1, 0.1, 0.1, 0.1
        }, isCalibrateable);
    }

	public LIBORCorrelationModelSeperateCurveSevenParameterExponentialDecay(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors, double[] fixedParameter, double[] calibrationParameter, boolean isCalibrateable) {
		super(timeDiscretization, liborPeriodDiscretization, isCalibrateable);

		this.numberOfFactors = numberOfFactors;
        this.fixedParameter = fixedParameter;
        this.isCalibrateable = isCalibrateable;
        this.calibrationParameter = calibrationParameter;
	}

	@Override
	public double[] getParameter() {
		return calibrationParameter;
	}

	@Override
	public void setParameter(double[] parameter) {
        this.calibrationParameter = parameter;

        calibrationParameter[0] = Math.max(calibrationParameter[0], 0.0);
        calibrationParameter[1] = Math.min(Math.max(calibrationParameter[1], 0.0), 1.0);
        calibrationParameter[2] = Math.max(calibrationParameter[2], 0.0);

		factorMatrix = null;
		correlationMatrix = null;
	}

	@Override
    public double	getFactorLoading(int timeIndex, int factor, int component) {
		if(factorMatrix == null) initialize();

		return factorMatrix[component][factor];
	}

	@Override
    public double	getCorrelation(int timeIndex, int component1, int component2) {
		if(correlationMatrix == null) initialize();

		return correlationMatrix[component1][component2];
	}

	@Override
    public int getNumberOfFactors() {
		return numberOfFactors;
	}

	private synchronized void initialize() {
        /*
         * Create instantaneous correlation matrix
		 */
        int timeSteps = liborPeriodDiscretization.getNumberOfTimeSteps();
        double[][] firstCurveCorrelationMatrix = new double[timeSteps][timeSteps];
        double[][] secondCurveCorrelationMatrix = new double[timeSteps][timeSteps];

        generateCorrelationMatrices(timeSteps, firstCurveCorrelationMatrix, secondCurveCorrelationMatrix);

        double[][] firstFactorMatrix = LinearAlgebra.factorReduction(firstCurveCorrelationMatrix, numberOfFactors / 2);
        double[][] secondFactorMatrix = LinearAlgebra.factorReduction(secondCurveCorrelationMatrix, numberOfFactors / 2);

        RealMatrix factorMatrix = new Array2DRowRealMatrix(2 * timeSteps, numberOfFactors);
        factorMatrix.setSubMatrix(firstFactorMatrix, 0, 0);
        factorMatrix.setSubMatrix(secondFactorMatrix, timeSteps, numberOfFactors / 2);

        this.factorMatrix = factorMatrix.getData();

        for (int i = 0; i < timeSteps; i++) {
            this.factorMatrix[i + timeSteps][0] = calibrationParameter[3];
        }

        this.correlationMatrix = new Array2DRowRealMatrix(this.factorMatrix)
                .multiply(new Array2DRowRealMatrix(this.factorMatrix).transpose()).getData();
    }

    protected void generateCorrelationMatrices(int timeSteps, double[][] firstCurveCorrelationMatrix, double[][] secondCurveCorrelationMatrix) {
        for (int row = 0; row < timeSteps; row++) {
            for (int col = row; col < timeSteps; col++) {
                // Exponentially decreasing instantaneous correlation
                double T1 = liborPeriodDiscretization.getTime(row);
                double T2 = liborPeriodDiscretization.getTime(col);

                double correlation1 = fixedParameter[1]
                        + (1 - fixedParameter[1]) * Math.exp(-fixedParameter[0] * Math.abs(T1 - T2) - fixedParameter[2] * Math.max(T1, T2));
                double correlation2 = calibrationParameter[1] - calibrationParameter[3] * calibrationParameter[3]
                        + (1 - calibrationParameter[1]) * Math.exp(-calibrationParameter[0] * Math.abs(T1 - T2) - calibrationParameter[2] * Math.max(T1, T2));

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
	public Object clone() {
		initialize();

		LIBORCorrelationModelSeperateCurveSevenParameterExponentialDecay newModel = new LIBORCorrelationModelSeperateCurveSevenParameterExponentialDecay(
				super.getTimeDiscretization(),
				super.getLiborPeriodDiscretization(),
				numberOfFactors,
                fixedParameter, calibrationParameter,
                isCalibrateable);

		newModel.correlationMatrix	= this.correlationMatrix;
		newModel.factorMatrix		= this.factorMatrix;
		
		return newModel;
	}
}
