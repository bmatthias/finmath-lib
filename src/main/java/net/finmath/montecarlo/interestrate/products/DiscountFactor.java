/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 09.02.2004
 */
package net.finmath.montecarlo.interestrate.products;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;

public class DiscountFactor extends AbstractLIBORMonteCarloProduct {

	private final double	maturity;

	public DiscountFactor(double maturity) {
		super();
		this.maturity = maturity;
	}

    @Override
    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		return new RandomVariable(evaluationTime, model.getModel().getDiscountCurve().getDiscountFactor(maturity));
    }
}
