package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.functions.LinearAlgebra;
import net.finmath.time.TimeDiscretizationInterface;

public abstract class AbstractLIBORCorrelationModelCombinedCurves extends AbstractMultiCurveLIBORCorrelationModel {

    public AbstractLIBORCorrelationModelCombinedCurves(TimeDiscretizationInterface timeDiscretization,
                                                       TimeDiscretizationInterface liborPeriodDiscretization,
                                                       int numberOfFactors,
                                                       double[] fixedParameter,
                                                       double[] calibrationParameter,
                                                       boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter, calibrationParameter, isCalibrateable);
    }

    @Override
    synchronized void initialize() {

        correlationMatrix = generateCorrelationMatrix();
        factorMatrix = LinearAlgebra.factorReduction(correlationMatrix, getNumberOfFactors());

        for (int component1 = 0; component1 < factorMatrix.length; component1++) {
            for (int component2 = component1 + 1; component2 < factorMatrix.length; component2++) {
                double correlation = 0.0;
                for (int factor = 0; factor < factorMatrix[component1].length; factor++) {
                    correlation += factorMatrix[component1][factor] * factorMatrix[component2][factor];
                }
                correlationMatrix[component1][component2] = correlation;
                correlationMatrix[component2][component1] = correlation;
            }
            correlationMatrix[component1][component1] = 1.0;
        }
    }

    abstract double[][] generateCorrelationMatrix();
}
