/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 09.02.2004
 */
package net.finmath.montecarlo.interestrate.products;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.MultiCurveLIBORMarketModel;
import net.finmath.stochastic.RandomVariableInterface;

public class Forward extends AbstractLIBORMonteCarloProduct {

	private final double	maturity;

	public Forward(double maturity) {
		super();
		this.maturity = maturity;
	}

    @Override
    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {
		int timeIndex = model.getTimeIndex(maturity);
        int liborIndex = model.getLiborPeriodIndex(maturity);

		RandomVariableInterface	forward					= model.getModel() instanceof MultiCurveLIBORMarketModel ?
                ((MultiCurveLIBORMarketModel)model.getModel()).getForward(timeIndex, liborIndex) : model.getLIBOR(timeIndex, liborIndex);
		RandomVariableInterface	numeraire				= model.getNumeraire(model.getTime(timeIndex + 1));
		RandomVariableInterface	monteCarloProbabilities	= model.getMonteCarloWeights(timeIndex + 1);

		forward = forward.div(numeraire).mult(monteCarloProbabilities);

		RandomVariableInterface	numeraireAtValuationTime				= model.getNumeraire(evaluationTime);
		RandomVariableInterface	monteCarloProbabilitiesAtValuationTime	= model.getMonteCarloWeights(evaluationTime);
		forward = forward.mult(numeraireAtValuationTime).div(monteCarloProbabilitiesAtValuationTime);

		return forward;
    }
}
