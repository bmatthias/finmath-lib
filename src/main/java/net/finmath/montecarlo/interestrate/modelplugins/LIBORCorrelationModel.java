/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 20.05.2006
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.time.TimeDiscretizationInterface;


/**
 * Abstract base class and interface description of a correlation model
 * (as it is used in {@link LIBORCovarianceModelFromVolatilityAndCorrelation}).
 * 
 * Derive from this class and implement the <code>getFactorLoading</code> method.
 * You have to call the constructor of this class to set the time
 * discretizations.
 * 
 * @author Christian Fries
 */
public abstract class LIBORCorrelationModel {
    final TimeDiscretizationInterface	timeDiscretization;
    final TimeDiscretizationInterface	liborPeriodDiscretization;
    final boolean isCalibrateable;
	
	public LIBORCorrelationModel(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, boolean isCalibrateable) {
		super();
		this.timeDiscretization = timeDiscretization;
		this.liborPeriodDiscretization = liborPeriodDiscretization;
        this.isCalibrateable = isCalibrateable;
    }

    public abstract double[]	getParameter();
    public abstract void		setParameter(double[] parameter);

    public abstract	double	getFactorLoading(int timeIndex, int factor, int component);
	public abstract	double	getCorrelation(int timeIndex, int component1, int component2);
	public abstract int		getNumberOfFactors();

	/**
	 * @return Returns the liborPeriodDiscretization.
	 */
	public TimeDiscretizationInterface getLiborPeriodDiscretization() {
		return liborPeriodDiscretization;
	}

	/**
	 * @return Returns the timeDiscretization.
	 */
	public TimeDiscretizationInterface getTimeDiscretization() {
		return timeDiscretization;
	}

	@Override
    public abstract Object clone();

    public boolean isCalibrateable() {
        return isCalibrateable;
    }
}
