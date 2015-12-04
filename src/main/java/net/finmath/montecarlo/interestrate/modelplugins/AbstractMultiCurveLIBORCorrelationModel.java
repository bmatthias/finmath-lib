package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.time.TimeDiscretizationInterface;

public abstract class AbstractMultiCurveLIBORCorrelationModel extends LIBORCorrelationModel {

    private int		numberOfFactors;
    protected final double[] fixedParameter;
    protected double[] calibrationParameter;

    protected transient double[][]	correlationMatrix;
    protected transient double[][]	factorMatrix;

    public AbstractMultiCurveLIBORCorrelationModel(TimeDiscretizationInterface timeDiscretization,
                                                   TimeDiscretizationInterface liborPeriodDiscretization,
                                                   int numberOfFactors,
                                                   double[] fixedParameter,
                                                   double[] calibrationParameter,
                                                   boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, isCalibrateable);
        this.numberOfFactors = numberOfFactors;
        this.fixedParameter = fixedParameter;
        this.calibrationParameter = calibrationParameter;
    }

    @Override
    public double[] getParameter() {
        return calibrationParameter;
    }

    @Override
    public void setParameter(double[] parameter) {
        if(!isCalibrateable) return;

        this.calibrationParameter = parameter;

        adjustParameters();

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

    abstract void initialize();
    abstract void adjustParameters();

    @Override
    public Object clone() {
        initialize();

        try {
            AbstractMultiCurveLIBORCorrelationModel newModel = this.getClass().getConstructor(
                    TimeDiscretizationInterface.class,
                    TimeDiscretizationInterface.class,
                    int.class, double[].class, double[].class, boolean.class).newInstance(
                    timeDiscretization, liborPeriodDiscretization, numberOfFactors, fixedParameter,
                    calibrationParameter, isCalibrateable
            );

            newModel.correlationMatrix	= this.correlationMatrix;
            newModel.factorMatrix		= this.factorMatrix;

            return newModel;
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }
}
