/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 09.02.2004
 */
package net.finmath.montecarlo.interestrate.products;

import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.stochastic.RandomVariableInterface;

public class LIBOR extends AbstractLIBORMonteCarloProduct {

	private final double	maturity;
	private final double	periodLength;

	public LIBOR(double maturity, double periodLength) {
		super();
		this.maturity = maturity;
		this.periodLength = periodLength;
	}

    @Override
    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) throws CalculationException {

		RandomVariableInterface	libor					= model.getLIBOR(maturity, maturity, maturity + periodLength);
		RandomVariableInterface	numeraire				= model.getNumeraire(maturity + periodLength);
		RandomVariableInterface	monteCarloProbabilities	= model.getMonteCarloWeights(model.getTimeIndex(maturity + periodLength));

		libor = libor.div(numeraire).mult(monteCarloProbabilities);

		RandomVariableInterface	numeraireAtValuationTime				= model.getNumeraire(evaluationTime);
		RandomVariableInterface	monteCarloProbabilitiesAtValuationTime	= model.getMonteCarloWeights(evaluationTime);
		libor = libor.mult(numeraireAtValuationTime).div(monteCarloProbabilitiesAtValuationTime);

		return libor;
    }
}
