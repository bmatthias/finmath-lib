/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 08.08.2005
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretizationInterface;

/**
 * @author Christian Fries
 */
public class LIBORVolatilityModelTwoCurveEightParameterExponentialForm extends LIBORVolatilityModel {

    private double[] fixedParameter;
    private double[] calibrationParameter;
    private boolean isCalibrateable;

    /**
     * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
     * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
     * @param fixedParameter
     */
    public LIBORVolatilityModelTwoCurveEightParameterExponentialForm(
            TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization,
            double[] fixedParameter, boolean isCalibrateable) {
        this(timeDiscretization, liborPeriodDiscretization, fixedParameter, new double[] {
                0.1, 0.1, 0.1, 0.2
        }, isCalibrateable);
    }

    /**
     * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
     * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
     * @param fixedParameter
     * @param calibrationParameter The parameters.
     */
    public LIBORVolatilityModelTwoCurveEightParameterExponentialForm(
            TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization,
            double[] fixedParameter, double[] calibrationParameter, boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, isCalibrateable);
        this.fixedParameter = fixedParameter;
        this.calibrationParameter = calibrationParameter;
        this.isCalibrateable = isCalibrateable;
    }

    @Override
    public double[] getParameter() {
        if(!isCalibrateable) return null;

        return calibrationParameter;
    }

    @Override
    public void setParameter(double[] parameter) {
        if(!isCalibrateable) return;

        this.calibrationParameter = parameter;
    }

    /* (non-Javadoc)
     * @see net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel#getVolatility(int, int)
     */
    @Override
    public RandomVariableInterface getVolatility(int timeIndex, int liborIndex) {
        // Create a very simple volatility model here
        double time = getTimeDiscretization().getTime(timeIndex);
        int timeSteps = getLiborPeriodDiscretization().getNumberOfTimeSteps();
        double timeToMaturity;
        double volatilityInstantaneous;
        double maturity = getLiborPeriodDiscretization().getTime(liborIndex <= timeSteps ? liborIndex : liborIndex - timeSteps);
        timeToMaturity = maturity - time;
        if (timeToMaturity <= 0) {
            volatilityInstantaneous = 0.0;   // This forward rate is already fixed, no volatility
        } else if (liborIndex < timeSteps) {
            volatilityInstantaneous = (fixedParameter[0] + fixedParameter[1] * timeToMaturity) * Math.exp(-fixedParameter[2] * timeToMaturity) + fixedParameter[3];
        } else {
            volatilityInstantaneous = (calibrationParameter[0] + calibrationParameter[1] * timeToMaturity) * Math.exp(-calibrationParameter[2] * timeToMaturity) + calibrationParameter[3];
        }
        if (volatilityInstantaneous < 0.0) volatilityInstantaneous = Math.max(volatilityInstantaneous, 0.0);
        return new RandomVariable(time, volatilityInstantaneous);
    }

	@Override
	public Object clone() {
		return new LIBORVolatilityModelTwoCurveEightParameterExponentialForm(
				super.getTimeDiscretization(),
				super.getLiborPeriodDiscretization(),
                fixedParameter, calibrationParameter, isCalibrateable
        );
	}
}
