/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 15.12.2007
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.Arrays;

/**
 * A covariance model build from a volatility model implementing
 * <code>LIBORVolatilityModel</code> and a correlation model
 * implementing <code>LIBORCorrelationModel</code>.
 * 
 * <p>
 * The model parameters are given by the concatenation of the
 * parameters of the <code>LIBORVolatilityModel</code> and
 * the parameters of the <code>LIBORCorrelationModel</code>,
 * in this ordering
 * </p>
 * 
 * @author Christian Fries
 */
public class LIBORCovarianceModelFromVolatilityAndCorrelation extends AbstractLIBORCovarianceModelParametric {

	private final LIBORVolatilityModel	volatilityModel;
	private final LIBORCorrelationModel	correlationModel;
	
	public LIBORCovarianceModelFromVolatilityAndCorrelation(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, LIBORVolatilityModel volatilityModel, LIBORCorrelationModel correlationModel) {
		super(timeDiscretization, liborPeriodDiscretization, correlationModel.getNumberOfFactors());

		this.volatilityModel = volatilityModel;
		this.correlationModel = correlationModel;
	}

	@Override
    public RandomVariableInterface[] getFactorLoading(int timeIndex, int component, RandomVariableInterface[] realizationAtTimeIndex) {
		RandomVariableInterface[] factorLoading = new RandomVariableInterface[correlationModel.getNumberOfFactors()];

		RandomVariableInterface volatility	= volatilityModel.getVolatility(timeIndex, component);
		for (int factorIndex = 0; factorIndex < factorLoading.length; factorIndex++) {
			factorLoading[factorIndex] = volatility.mult(correlationModel.getFactorLoading(timeIndex, factorIndex, component));
		}
		
		return factorLoading;
	}
	
	@Override
    public RandomVariableInterface getFactorLoadingPseudoInverse(int timeIndex, int component, int factor, RandomVariableInterface[] realizationAtTimeIndex) {
		// Note that we assume that the correlation model getFactorLoading gives orthonormal vectors
		RandomVariableInterface factorLoadingPseudoInverse = volatilityModel.getVolatility(timeIndex, component).invert()
                .mult(correlationModel.getFactorLoading(timeIndex, factor, component));

        // @todo numberOfComponents should be stored as a member?!
        int numberOfComponents = getLiborPeriodDiscretization().getNumberOfTimeSteps();
        
        double factorWeight = 0.0;
        for(int componentIndex=0; componentIndex<numberOfComponents; componentIndex++) {
            double factorElement = correlationModel.getFactorLoading(timeIndex, factor, componentIndex);            
            factorWeight +=  factorElement*factorElement;                                                                                                                 
        }

        factorLoadingPseudoInverse = factorLoadingPseudoInverse.mult(1/factorWeight);

        return factorLoadingPseudoInverse;		
	}

    /* (non-Javadoc)
     * @see net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel#getCovariance(int, int, int)
     */
    @Override
    public RandomVariableInterface getCovariance(int timeIndex, int component1, int component2, RandomVariableInterface[] realizationAtTimeIndex) {
        RandomVariableInterface covariance = new RandomVariable(0.0, correlationModel.getCorrelation(timeIndex, component1, component2));
        covariance = covariance.mult(volatilityModel.getVolatility(timeIndex, component1))
                .mult(volatilityModel.getVolatility(timeIndex, component2));

        return covariance;
    }

	@Override
	public double[] getParameter() {
		double[] volatilityParameter	= volatilityModel.getParameter();
		double[] correlationParameter	= correlationModel.getParameter();
		
		int parameterLength = 0;
		parameterLength += volatilityModel.isCalibrateable() ? volatilityParameter.length : 0;
		parameterLength += correlationModel.isCalibrateable() ? correlationParameter.length : 0;
		
		double[] parameter = new double[parameterLength];

		int parameterIndex = 0;
		if(volatilityModel.isCalibrateable()) {
			System.arraycopy(volatilityParameter, 0, parameter, parameterIndex, volatilityParameter.length);
			parameterIndex += volatilityParameter.length;
		}
		if(correlationModel.isCalibrateable()) {
			System.arraycopy(correlationParameter, 0, parameter, parameterIndex, correlationParameter.length);
		}

		return parameter;
	}

	@Override
	public void setParameter(double[] parameter) {
		double[] volatilityParameter = volatilityModel.getParameter();
		double[] correlationParameter = correlationModel.getParameter();

		int parameterIndex = 0;
		if(volatilityModel.isCalibrateable()) {
			double[] newVolatilityParameter = new double[volatilityParameter.length];
			System.arraycopy(parameter, parameterIndex, newVolatilityParameter, 0, newVolatilityParameter.length);
			parameterIndex += newVolatilityParameter.length;
			if(!Arrays.equals(newVolatilityParameter, volatilityModel.getParameter()))
				volatilityModel.setParameter(newVolatilityParameter);
		}
		if(correlationModel.isCalibrateable()) {
			double[] newCorrelationParameter = new double[correlationParameter.length];
			System.arraycopy(parameter, parameterIndex, newCorrelationParameter, 0, newCorrelationParameter.length);
			if(!Arrays.equals(newCorrelationParameter, correlationModel.getParameter()))
				correlationModel.setParameter(newCorrelationParameter);
		}


        //IMPORTANT: Otherwise the optimizer may not know that the parameters have been adjusted,
        //leading to bad calibration results!
        if (!Arrays.equals(getParameter(), parameter)) {
            System.arraycopy(getParameter(), 0, parameter, 0, parameter.length);
        }
	}

	@Override
	public Object clone() {
		return new LIBORCovarianceModelFromVolatilityAndCorrelation(
				this.getTimeDiscretization(),
				this.getLiborPeriodDiscretization(),
				(LIBORVolatilityModel)volatilityModel.clone(), (LIBORCorrelationModel)correlationModel.clone());
	}

	public LIBORVolatilityModel getVolatilityModel() {
		return volatilityModel;
	}

	public LIBORCorrelationModel getCorrelationModel() {
		return correlationModel;
	}
}
