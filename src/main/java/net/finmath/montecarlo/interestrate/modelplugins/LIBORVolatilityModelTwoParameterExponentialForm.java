/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 08.08.2005
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.montecarlo.RandomVariable;
import net.finmath.time.TimeDiscretizationInterface;

/**
 * Implements the volatility model &sigma;<sub>i</sub>(t<sub>j</sub>) = a * exp(-b (T<sub>i</sub>-t<sub>j</sub>))
 * 
 * @author Christian Fries
 */
public class LIBORVolatilityModelTwoParameterExponentialForm extends LIBORVolatilityModel {

    private double a;
    private double b;

    /**
     * Creates the volatility model &sigma;<sub>i</sub>(t<sub>j</sub>) = a * exp(-b (T<sub>i</sub>-t<sub>j</sub>))
     * 
     * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
     * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
     * @param a The parameter a: an initial volatility level.
     * @param b The parameter b: exponential decay of the volatility.
     */
    public LIBORVolatilityModelTwoParameterExponentialForm(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, double a, double b) {
    	this(timeDiscretization, liborPeriodDiscretization, a, b, false);
    }

    /**
     * Creates the volatility model &sigma;<sub>i</sub>(t<sub>j</sub>) = a * exp(-b (T<sub>i</sub>-t<sub>j</sub>))
     * 
     * @param timeDiscretization The simulation time discretization t<sub>j</sub>.
     * @param liborPeriodDiscretization The period time discretization T<sub>i</sub>.
     * @param a The parameter a: an initial volatility level.
     * @param b The parameter b: exponential decay of the volatility.
     * @param isCalibrateable Set this to true, if the parameters are available for calibration.
     */
    public LIBORVolatilityModelTwoParameterExponentialForm(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, double a, double b, boolean isCalibrateable) {
        super(timeDiscretization, liborPeriodDiscretization, isCalibrateable);
        this.a = a;
        this.b = b;
    }
    
	@Override
	public double[] getParameter() {
		double[] parameter = new double[2];
		parameter[0] = a;
		parameter[1] = b;

		return parameter;
	}

	@Override
	public void setParameter(double[] parameter) {
        if(!isCalibrateable) return;

		this.a = parameter[0];
        this.b = parameter[1];
	}

    /* (non-Javadoc)
     * @see net.finmath.montecarlo.interestrate.modelplugins.LIBORVolatilityModel#getVolatility(int, int)
     */
    @Override
    public RandomVariable getVolatility(int timeIndex, int liborIndex) {
        // Create a very simple volatility model here
        double time             = getTimeDiscretization().getTime(timeIndex);
        double maturity         = getLiborPeriodDiscretization().getTime(liborIndex);
        double timeToMaturity   = maturity-time;

        double volatilityInstanteaneous; 
        if(timeToMaturity <= 0)
        {
            volatilityInstanteaneous = 0;   // This forward rate is already fixed, no volatility
        }
        else
        {
            volatilityInstanteaneous = a * Math.exp(-b * timeToMaturity);
        }

        return new RandomVariable(getTimeDiscretization().getTime(timeIndex), volatilityInstanteaneous);
    }

	@Override
	public Object clone() {
		return new LIBORVolatilityModelTwoParameterExponentialForm(
				getTimeDiscretization(),
				getLiborPeriodDiscretization(),
				a,
				b,
				isCalibrateable
				);
	}
}
