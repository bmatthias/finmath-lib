/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 17.05.2007
 * Created on 30.03.2014
 */
package net.finmath.montecarlo.interestrate.products;

import net.finmath.functions.AnalyticFormulas;
import net.finmath.integration.SimpsonRealIntegrator;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.MultiCurveLIBORMarketModel;
import net.finmath.stochastic.RandomVariableInterface;

import java.util.Arrays;
import java.util.function.DoubleUnaryOperator;

public class CapletAnalyticApproximation extends AbstractLIBORMonteCarloProduct {

    public enum MultiCurveApproximation {
        INTEGRATEDVARIANCE,
        FENTON_WILKINSON,
        HO_SCHWARTZ_YEH
    }

    private final double    strike;
    private final double	maturity;
    private ValueUnit   valueUnit;

    public CapletAnalyticApproximation(double strike, double maturity, ValueUnit valueUnit) {
        super();
        this.strike = strike;
        this.maturity	= maturity;
        this.valueUnit	= valueUnit;
    }

    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model, MultiCurveApproximation multiCurveApproximation, ValueUnit valueUnit) {
        return getValues(evaluationTime, model.getModel(), multiCurveApproximation, valueUnit);
    }

    @Override
    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) {
        return getValue(evaluationTime, model, MultiCurveApproximation.FENTON_WILKINSON, this.valueUnit);
    }

    public RandomVariableInterface getValues(double evaluationTime, LIBORMarketModelInterface model, MultiCurveApproximation multiCurveApproximation, ValueUnit valueUnit) {
        if(evaluationTime > 0) throw new RuntimeException("Forward start evaluation currently not supported.");

        int maturityIndex  = model.getLiborPeriodIndex(maturity);
        double periodLength = model.getLiborPeriodDiscretization().getTimeStep(maturityIndex);

        double integratedVariance;
        if (model instanceof MultiCurveLIBORMarketModel) {
            switch (multiCurveApproximation) {
                case FENTON_WILKINSON:
                    integratedVariance = calculateIntegratedCovarianceFentonWilkinson((MultiCurveLIBORMarketModel) model);
                    break;
                case HO_SCHWARTZ_YEH:
                    integratedVariance = calculateIntegratedCovarianceHoSchwartzYeh((MultiCurveLIBORMarketModel) model);
                    break;
                default:
                    integratedVariance = calculateIntegratedCovarianceMultiCurve((MultiCurveLIBORMarketModel) model);
                    break;
            }
        } else {
            integratedVariance = model.getIntegratedLIBORCovariance()[maturityIndex - 1][maturityIndex][maturityIndex];
        }

        if(valueUnit == ValueUnit.INTEGRATEDVARIANCE) return new RandomVariable(evaluationTime, integratedVariance);

        double volatility		= Math.sqrt(integratedVariance / maturity);

        if(valueUnit == ValueUnit.VOLATILITY) return new RandomVariable(evaluationTime, volatility);

        double valueCaplet = AnalyticFormulas.blackModelCapletValue(
                model.getForwardRateCurve().getForward(model.getAnalyticModel(), maturity),
                volatility, maturity, strike, periodLength,
                model.getDiscountCurve().getDiscountFactor(model.getAnalyticModel(), maturity + periodLength)
        );
        return new RandomVariable(evaluationTime, valueCaplet);
    }

    private double calculateIntegratedCovarianceMultiCurve(MultiCurveLIBORMarketModel model) {
        int maturityIndex = model.getLiborPeriodIndex(maturity);

        double timeStep = model.getLiborPeriodDiscretization().getTimeStep(maturityIndex);

        RandomVariableInterface[] initialState = model.getInitialState();
        double forward = model.applyStateSpaceTransform(maturityIndex, initialState[maturityIndex]).getAverage();
        double spread = model.applyStateSpaceTransform(maturityIndex + model.getNumberOfComponents() / 2, initialState[maturityIndex + model.getNumberOfComponents() / 2]).getAverage();

        double[] integratedLIBORCovariance = new double[3];

        for (int timeIndex = 0; timeIndex < model.getTimeIndex(maturity); timeIndex++) {
            double dt = model.getTime(timeIndex + 1) - model.getTime(timeIndex);
            RandomVariableInterface[] forwardFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex, null);
            RandomVariableInterface[] spreadFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex + model.getNumberOfComponents() / 2, null);

            for(int factorIndex = 0; factorIndex < model.getNumberOfFactors(); factorIndex++) {
                integratedLIBORCovariance[0] +=
                        forwardFactorLoading[factorIndex].get(0) * forwardFactorLoading[factorIndex].get(0) * dt;

                integratedLIBORCovariance[1] +=
                        spreadFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;

                integratedLIBORCovariance[2] +=
                        2.0 * forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
            }
        }

        double coeff1 = 1.0, coeff2 = 1.0, coeff3 = 1.0;
        if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
            coeff1 = (1.0 + timeStep * spread) * (1.0 + timeStep * spread);
            coeff2 = (1.0 + timeStep * forward) * (1.0 + timeStep * forward);
            coeff3 = (1.0 + timeStep * spread) * (1.0 + timeStep * forward);
        }

        double integratedSwapRateVariance = coeff1 * forward * forward * integratedLIBORCovariance[0]
                + coeff2 * spread * spread * integratedLIBORCovariance[1]
                + coeff3 * forward * spread * integratedLIBORCovariance[2];

        return integratedSwapRateVariance  / (strike * strike);
    }

    private double calculateIntegratedCovarianceFentonWilkinson(MultiCurveLIBORMarketModel model) {
        int maturityIndex = model.getLiborPeriodIndex(maturity);

        RandomVariableInterface[] realization = model.getInitialState();
        double[] integratedCovariance = new double[3];
        for (int timeIndex = 0; timeIndex < model.getTimeIndex(maturity); timeIndex++) {
            double dt = model.getTime(timeIndex + 1) - model.getTime(timeIndex);
            RandomVariableInterface[] forwardFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex, null);
            RandomVariableInterface[] spreadFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex + model.getNumberOfComponents() / 2, null);

            for(int factorIndex = 0; factorIndex < model.getNumberOfFactors(); factorIndex++) {
                integratedCovariance[0] += forwardFactorLoading[factorIndex].get(0) * forwardFactorLoading[factorIndex].get(0) * dt;
                integratedCovariance[1] += spreadFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
                integratedCovariance[2] += forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
            }
        }

        double gaussianMeanOfForward = realization[maturityIndex].getAverage() - 0.5 * integratedCovariance[0];
        double gaussianMeanOfSpread = realization[maturityIndex + model.getNumberOfComponents() / 2].getAverage() - 0.5 * integratedCovariance[1];

        double u1 = Math.exp(gaussianMeanOfForward + 0.5 * integratedCovariance[0])
                + Math.exp(gaussianMeanOfSpread + 0.5 * integratedCovariance[1]);
        double u2 = Math.exp(2.0 * (gaussianMeanOfForward + integratedCovariance[0]))
                + Math.exp(2.0 * (gaussianMeanOfSpread + integratedCovariance[1]))
                + 2.0 * Math.exp(gaussianMeanOfForward + gaussianMeanOfSpread +
                0.5 * integratedCovariance[0] + 0.5 * integratedCovariance[1] + integratedCovariance[2]);

        if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
            double timeStep = model.getLiborPeriodDiscretization().getTimeStep(maturityIndex);
            double gaussianMeanOfProduct = Math.log(timeStep) + gaussianMeanOfForward + gaussianMeanOfSpread;
            double varianceOfProduct = integratedCovariance[0] + integratedCovariance[1] + 2.0 * integratedCovariance[2];

            u1 += Math.exp(gaussianMeanOfProduct + 0.5 * varianceOfProduct);
            u2 += Math.exp(2.0 * (gaussianMeanOfProduct + varianceOfProduct))
                    + 2.0 * Math.exp(gaussianMeanOfForward + gaussianMeanOfProduct +
                    0.5 * varianceOfProduct + 0.5 * integratedCovariance[1] + 0.0)
                    + 2.0 * Math.exp(gaussianMeanOfForward + gaussianMeanOfProduct +
                    0.5 * integratedCovariance[0] + 0.5 * varianceOfProduct + 0.0);
        }

        return Math.log(u2) - 2.0 * Math.log(u1);
    }

    private double calculateIntegratedCovarianceHoSchwartzYeh(MultiCurveLIBORMarketModel model) {
        int maturityIndex = model.getLiborPeriodIndex(maturity);

        RandomVariableInterface[] realization = model.getInitialState();
        int numberOfComponents = model.getNumberOfComponents() / 2;
        double[] integratedCovariance = new double[3];
        for (int timeIndex = 0; timeIndex < model.getTimeIndex(maturity); timeIndex++) {
            double dt = model.getTime(timeIndex + 1) - model.getTime(timeIndex);
            RandomVariableInterface[] forwardFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex, null);
            RandomVariableInterface[] spreadFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex + numberOfComponents, null);

            for(int factorIndex = 0; factorIndex < model.getNumberOfFactors(); factorIndex++) {
                integratedCovariance[0] += forwardFactorLoading[factorIndex].get(0) * forwardFactorLoading[factorIndex].get(0) * dt;
                integratedCovariance[1] += spreadFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
                integratedCovariance[2] += forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
            }
        }

        double gaussianMeanOfForward = realization[maturityIndex].getAverage() - 0.5 * integratedCovariance[0];
        double gaussianMeanOfSpread = realization[maturityIndex + numberOfComponents].getAverage() - 0.5 * integratedCovariance[1];

        double[] mu = new double[]{ gaussianMeanOfForward, gaussianMeanOfSpread };
        Arrays.sort(mu);

        double mw = mu[0] - mu[1];
        double variancew = integratedCovariance[0] + integratedCovariance[1] - 2.0 * integratedCovariance[2];
        double variancey1 = mu[1] == gaussianMeanOfForward ? integratedCovariance[0] : integratedCovariance[1];
        double sigmaw = Math.sqrt(variancew);
        double eta = -mw / sigmaw;
        double sqrtTwoPi = Math.sqrt(2.0 * Math.PI);

        DoubleUnaryOperator fw = (w) -> Math.exp(-(w - mw) * (w - mw) / (2.0 * variancew)) / (sqrtTwoPi * sigmaw);

        double i0 = integrate((v) -> Math.exp(-0.5 * (Math.log(v) - eta) * (Math.log(v) - eta)) / sqrtTwoPi);
        double i1 = integrate((v) -> Math.log(1.0 + v) * (fw.applyAsDouble(Math.log(v)) + fw.applyAsDouble(-Math.log(v))));
        double i2 = integrate((v) -> (fw.applyAsDouble(Math.log(v)) + fw.applyAsDouble(-Math.log(v))) / (1.0 + 1.0 / v));
        double i3 = integrate((v) -> {
            double log = Math.log(1.0 + v);
            return log * log * (fw.applyAsDouble(Math.log(v)) + fw.applyAsDouble(-Math.log(v)));
        });
        double i4 = integrate((v) -> -Math.log(1.0 + v) * Math.log(v) * fw.applyAsDouble(-Math.log(v)));

        double a0 = sigmaw / sqrtTwoPi * Math.exp(-0.5 * eta * eta) + mw * i0;
        double g1 = a0 + i1;
        double g2 = i3 + 2.0 * i4 + variancew * i0 + mw * a0;

        return variancey1 - g1 * g1 + g2 - 2.0 * variancey1 * (i2 + i0);
    }

    private double integrate(DoubleUnaryOperator function) {
        return new SimpsonRealIntegrator(1e-15, 1.0, 100, true).integrate(v -> function.applyAsDouble(v) / v);
    }
}
