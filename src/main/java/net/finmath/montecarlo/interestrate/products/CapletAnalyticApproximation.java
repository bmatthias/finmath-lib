/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 17.05.2007
 * Created on 30.03.2014
 */
package net.finmath.montecarlo.interestrate.products;

import net.finmath.functions.AnalyticFormulas;
import net.finmath.functions.NormalDistribution;
import net.finmath.integration.SimpsonRealIntegrator;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulationInterface;
import net.finmath.montecarlo.interestrate.MultiCurveLIBORMarketModel;
import net.finmath.montecarlo.interestrate.ShiftedLIBORMarketModelInterface;
import net.finmath.stochastic.RandomVariableInterface;
import org.apache.commons.math3.util.FastMath;

import java.util.function.DoubleUnaryOperator;

public class CapletAnalyticApproximation extends AbstractLIBORMonteCarloProduct {

    public enum MultiCurveApproximation {
        DRIFT_FREEZE,
        INTEGRATED_EXPECTATION,
        LEVY,
        JU,
        HO
    }

    private final double    strike;
    private final double	maturity;
    private ValueUnit   valueUnit;
    private MultiCurveApproximation multiCurveApproximation;

    public CapletAnalyticApproximation(double strike, double maturity, ValueUnit valueUnit) {
        this(strike, maturity, valueUnit, MultiCurveApproximation.INTEGRATED_EXPECTATION);
    }

    public CapletAnalyticApproximation(double strike, double maturity, ValueUnit valueUnit, MultiCurveApproximation multiCurveApproximation) {
        super();
        this.strike = strike;
        this.maturity	= maturity;
        this.valueUnit	= valueUnit;
        this.multiCurveApproximation = multiCurveApproximation;
    }

    @Override
    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) {
        return getValue(evaluationTime, model, this.multiCurveApproximation, this.valueUnit);
    }

    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model, MultiCurveApproximation multiCurveApproximation, ValueUnit valueUnit) {
        return getValues(evaluationTime, model.getModel(), multiCurveApproximation, valueUnit);
    }

    public RandomVariableInterface getValues(double evaluationTime, LIBORMarketModelInterface model, MultiCurveApproximation multiCurveApproximation, ValueUnit valueUnit) {
        if(evaluationTime > 0) throw new RuntimeException("Forward start evaluation currently not supported.");

        int maturityIndex  = model.getLiborPeriodIndex(maturity);
        double periodLength = model.getLiborPeriodDiscretization().getTimeStep(maturityIndex);

        double integratedVariance;
        if (model instanceof MultiCurveLIBORMarketModel && ((MultiCurveLIBORMarketModel) model).isMultiCurve()) {
            switch (multiCurveApproximation) {
                case DRIFT_FREEZE:
                case INTEGRATED_EXPECTATION:
                    integratedVariance = calculateIntegratedVarianceMultiCurve((MultiCurveLIBORMarketModel) model);
                    break;
                case LEVY:
                case JU:
                    double[] params = calculateIntegratedCovarianceLevy((MultiCurveLIBORMarketModel) model);
                    return getValue(evaluationTime, model, valueUnit, params);
                case HO:
                    params = calculateIntegratedCovarianceHoSchwartzYeh((MultiCurveLIBORMarketModel) model);
                    return getValue(evaluationTime, model, valueUnit, params);
                default:
                    throw new UnsupportedOperationException("Approximation method" + multiCurveApproximation +  "not implemented.");
            }
        } else {
            integratedVariance = model.getIntegratedLIBORCovariance()[model.getTimeDiscretization().getTimeIndex(maturity) - 1][maturityIndex][maturityIndex];
        }

        if(valueUnit == ValueUnit.INTEGRATEDVARIANCE) return new RandomVariable(evaluationTime, integratedVariance);

        double volatility		= Math.sqrt(integratedVariance / maturity);

        if(valueUnit == ValueUnit.VOLATILITY) return new RandomVariable(evaluationTime, volatility);

        double shift = (model instanceof ShiftedLIBORMarketModelInterface) ? ((ShiftedLIBORMarketModelInterface)model).getLIBORShift(maturityIndex) : 0.0;
        double forward = model.getForwardRateCurve().getForward(model.getAnalyticModel(), maturity);

        double valueCaplet = AnalyticFormulas.blackModelCapletValue(
                forward + shift, volatility, maturity, strike + shift, periodLength,
                model.getDiscountCurve().getDiscountFactor(model.getAnalyticModel(), maturity + periodLength)
        );
        return new RandomVariable(evaluationTime, valueCaplet);
    }

    private RandomVariableInterface getValue(double evaluationTime, LIBORMarketModelInterface model, ValueUnit valueUnit, double[] params) {
        int maturityIndex  = model.getLiborPeriodIndex(maturity);
        double periodLength = model.getLiborPeriodDiscretization().getTimeStep(maturityIndex);

        double shift = (model instanceof ShiftedLIBORMarketModelInterface) ? ((ShiftedLIBORMarketModelInterface)model).getLIBORShift(maturityIndex) : 0.0;
        double forward = model.getForwardRateCurve().getForward(model.getAnalyticModel(), maturity);

        double blackStrike = strike + shift;

        double payoffUnit = periodLength * model.getDiscountCurve().getDiscountFactor(model.getAnalyticModel(), maturity + periodLength);

        double valueCaplet = optionValue(forward + shift, blackStrike, params[0], params[1], payoffUnit) + payoffUnit * params[2];

        if(valueUnit == ValueUnit.VALUE) return new RandomVariable(evaluationTime, valueCaplet);
        double vola = AnalyticFormulas.blackScholesOptionImpliedVolatility(forward + shift, maturity, blackStrike, payoffUnit, valueCaplet);

        if (valueUnit == ValueUnit.VOLATILITY) {
            return new RandomVariable(evaluationTime, vola);
        } else {
            return new RandomVariable(evaluationTime, vola * vola * maturity);
        }
    }

    public static double optionValue(double forward, double strike, double mean, double variance, double payoffUnit) {
        double dPlus = (mean - Math.log(strike) + variance) / Math.sqrt(variance);
        double dMinus = dPlus - Math.sqrt(variance);
        double valueCaplet = payoffUnit * (forward * NormalDistribution.cumulativeDistribution(dPlus) - strike * NormalDistribution.cumulativeDistribution(dMinus));
        return valueCaplet;
    }

    public double calculateIntegratedVarianceMultiCurve(MultiCurveLIBORMarketModel model) {
        int maturityIndex = model.getLiborPeriodIndex(maturity);
        int spreadIndex = maturityIndex + model.getNumberOfLibors();

        double timeStep = model.getLiborPeriodDiscretization().getTimeStep(maturityIndex);

        RandomVariableInterface[] initialState = model.getInitialState();
        double shiftedForward = model.applyStateSpaceTransform(maturityIndex, initialState[maturityIndex]).getAverage();
        double shiftedSpread = model.applyStateSpaceTransform(spreadIndex, initialState[spreadIndex]).getAverage();
        double forward = shiftedForward - model.getShiftParameter(maturityIndex);

        double[] integratedLIBORCovariance = new double[3];

        for (int timeIndex = 0; timeIndex < model.getTimeIndex(maturity); timeIndex++) {
            double dt = model.getTime(timeIndex + 1) - model.getTime(timeIndex);
            RandomVariableInterface[] forwardFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex, model.getInitialValue());
            RandomVariableInterface[] spreadFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, spreadIndex, model.getInitialValue());

            for(int factorIndex = 0; factorIndex < model.getNumberOfFactors(); factorIndex++) {
                integratedLIBORCovariance[0] +=
                        forwardFactorLoading[factorIndex].get(0) * forwardFactorLoading[factorIndex].get(0) * dt;

                integratedLIBORCovariance[1] +=
                        spreadFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;

                integratedLIBORCovariance[2] +=
                        forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;

                if (multiCurveApproximation == MultiCurveApproximation.INTEGRATED_EXPECTATION &&
                        model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
                    shiftedSpread -= timeStep / (1 + timeStep * forward) * shiftedForward * shiftedSpread
                            * forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
                }
            }
        }

        double coeff1 = 1.0, coeff2 = 1.0, coeff3 = 1.0;
        if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
            coeff1 = (1.0 + timeStep * shiftedSpread) * (1.0 + timeStep * shiftedSpread);
            coeff2 = (1.0 + timeStep * shiftedForward) * (1.0 + timeStep * shiftedForward);
            coeff3 = (1.0 + timeStep * shiftedSpread) * (1.0 + timeStep * shiftedForward);
        }

        double variance = coeff1 * shiftedForward * shiftedForward * integratedLIBORCovariance[0]
                + coeff2 * shiftedSpread * shiftedSpread * integratedLIBORCovariance[1]
                + 2.0 * coeff3 * shiftedForward * shiftedSpread * integratedLIBORCovariance[2];

        double shiftedLibor = shiftedForward + shiftedSpread;
        if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
            shiftedLibor += timeStep * shiftedForward * shiftedSpread;
        }

        return variance  / (shiftedLibor * shiftedLibor);
    }

    //It may be inefficient to get the integrated covariance for all components from the model -
    //however, since we are most likely reusing the model and the integrated covariance is cached in the model,
    //in most cases we will actually save some calculations.
    private double[] calculateIntegratedCovarianceLevy(MultiCurveLIBORMarketModel model) {
        int maturityIndex = model.getLiborPeriodIndex(maturity);
        int spreadIndex = maturityIndex + model.getNumberOfComponents() / 2;
        double timeStep = model.getLiborPeriodDiscretization().getTimeStep(maturityIndex);

        RandomVariableInterface[] realization = model.getInitialState();

        double shiftedForward = model.applyStateSpaceTransform(maturityIndex, realization[maturityIndex]).getAverage();
        double shiftedSpread = model.applyStateSpaceTransform(spreadIndex, realization[spreadIndex]).getAverage();

        double[][] integratedCovariance = new double[2][2];
        for (int timeIndex = 0; timeIndex < model.getTimeIndex(maturity); timeIndex++) {
            double dt = model.getTime(timeIndex + 1) - model.getTime(timeIndex);
            RandomVariableInterface[] forwardFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, maturityIndex, null);
            RandomVariableInterface[] spreadFactorLoading = model.getCovarianceModel().getFactorLoading(timeIndex, spreadIndex, null);

            for(int factorIndex = 0; factorIndex < model.getNumberOfFactors(); factorIndex++) {
                integratedCovariance[0][0] += forwardFactorLoading[factorIndex].get(0) * forwardFactorLoading[factorIndex].get(0) * dt;
                integratedCovariance[1][1] += spreadFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
                integratedCovariance[0][1] += forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
                integratedCovariance[1][0] += forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;

                if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
                    shiftedSpread -= timeStep / (1 + timeStep * shiftedForward) * shiftedForward * shiftedSpread *
                            forwardFactorLoading[factorIndex].get(0) * spreadFactorLoading[factorIndex].get(0) * dt;
                }
            }
        }

        double[] mu = new double[] { shiftedForward, shiftedSpread };
        double[] moments;
        if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
            moments = getMomentsLevyMult(mu, integratedCovariance, timeStep);
        } else {
            moments = getMomentsLevy(mu, integratedCovariance);
        }

        if (multiCurveApproximation == MultiCurveApproximation.JU) {
            double adjustment = getTaylorExpansionByJuAdjustment(mu, moments, integratedCovariance, strike + model.getLIBORShift(maturityIndex));
            return new double[]{ moments[0], moments[1], adjustment };
        } else {
            return moments;
        }
    }

    public static double[] getMomentsLevy(double[] expMeans, double[][] integratedCovariance) {
        double shiftedForward = expMeans[0];
        double shiftedSpread = expMeans[1];

        double u1 = shiftedForward + shiftedSpread;
        double u2 = shiftedForward * shiftedForward * Math.exp(integratedCovariance[0][0])
                + shiftedSpread * shiftedSpread * Math.exp(integratedCovariance[1][1])
                + 2.0 * shiftedForward * shiftedSpread * Math.exp(integratedCovariance[0][1]);

        double mean = 2.0 * Math.log(u1) - 0.5 * Math.log(u2);
        double variance = Math.log(u2) - 2.0 * Math.log(u1);

        return new double[] { mean, variance, 0.0 };
    }

    public static double[] getMomentsLevyMult(double[] expMeans, double[][] integratedCovariance, double timeStep) {
        double shiftedForward = expMeans[0];
        double shiftedSpread = expMeans[1];

        double u1 = shiftedForward + shiftedSpread;
        double u2 = shiftedForward * shiftedForward * Math.exp(integratedCovariance[0][0])
                + shiftedSpread * shiftedSpread * Math.exp(integratedCovariance[1][1])
                + 2.0 * shiftedForward * shiftedSpread * Math.exp(integratedCovariance[0][1]);

        double varianceOfProduct = integratedCovariance[0][0] + integratedCovariance[1][1] + 2.0 * integratedCovariance[0][1];
        double product = timeStep * shiftedForward * shiftedSpread;

        u1 += product;
        u2 += product * product * Math.exp(varianceOfProduct)
                + 2.0 * shiftedSpread * product * Math.exp(0.5 * varianceOfProduct + 0.5 * integratedCovariance[1][1] + (integratedCovariance[1][1] + integratedCovariance[0][1]) / Math.sqrt( varianceOfProduct * integratedCovariance[1][1]))
                + 2.0 * shiftedForward * product * Math.exp(0.5 * integratedCovariance[0][0] + 0.5 * varianceOfProduct + (integratedCovariance[0][0] + integratedCovariance[0][1]) / Math.sqrt(varianceOfProduct *integratedCovariance[0][0]));

        double mean = 2.0 * Math.log(u1) - 0.5 * Math.log(u2);
        double variance = Math.log(u2) - 2.0 * Math.log(u1);

        return new double[] { mean, variance, 0.0 };
    }

    public static double getTaylorExpansionByJuAdjustment(double[] expMeans, double[] momentsLevy, double[][] integratedCovariance, double strike) {
        double u20 = expMeans[0] * expMeans[0] + expMeans[1] * expMeans[1] + 2.0 * expMeans[0] * expMeans[1];
        double u21 = expMeans[0] * expMeans[0] * integratedCovariance[0][0]
                + expMeans[1] * expMeans[1] * integratedCovariance[1][1]
                + 2.0 * expMeans[0] * expMeans[1] * integratedCovariance[0][1];
        double u22 = expMeans[0] * expMeans[0] * integratedCovariance[0][0] * integratedCovariance[0][0]
                + expMeans[1] * expMeans[1] * integratedCovariance[1][1] * integratedCovariance[1][1]
                + 2.0 * expMeans[0] * expMeans[1] * integratedCovariance[0][1] * integratedCovariance[0][1];
        double u23 = expMeans[0] * expMeans[0] * integratedCovariance[0][0] * integratedCovariance[0][0] * integratedCovariance[0][0]
                + expMeans[1] * expMeans[1] * integratedCovariance[1][1] * integratedCovariance[1][1] * integratedCovariance[1][1]
                + 2.0 * expMeans[0] * expMeans[1] * integratedCovariance[0][1] * integratedCovariance[0][1] * integratedCovariance[0][1];

        double a1 = -0.5 * u21 / u20;
        double a2 = 2.0 * a1 * a1 - 0.5 * u22 / u20;
        double a3 = 6.0 * a1 * a2 - 4 * a1 * a1 * a1 - 0.5 * u23 / u20;

        double A0 = expMeans[0] + expMeans[1];

        double[] sums = sums(expMeans, integratedCovariance);

        double b1 = 0.5 / (A0 * A0 * A0) * sums[0];
        double b2 = a1 * a1 - 0.5 * a2;

        double c1 = -a1 * b1;
        double c2 = 1.0 / 144.0 / (A0 * A0 * A0 * A0) * (72.0 * sums[1] + 18.0 * u21 * u22 + 24.0 * sums[2]);
        double c3 = 1.0 / 48.0 / (A0 * A0 * A0) * (24.0 * sums[3] + 8.0 * sums[4]);
        double c4 = a1 * a2 - 2.0 / 3.0 * a1 * a1 * a1 - a3 / 6.0;

        //double d1 = 0.5 * (6.0 * a1 * a1 - a2 - 4.0 * b1 + 2.0 * b2) - 1.0 / 6.0 * (120.0 * a1 * a1 * a1 - a3 + 6.0 * (24.0 * c1 - 6.0 * c2 + 2.0 * c3 - c4));
        double d2 = 0.5 * (10.0 * a1 * a1 + a2 - 6.0 * b1 + 2.0 * b2) - (128.0 / 3.0 * a1 * a1 * a1 - 1.0 / 6.0 * a3 + 2.0 * a1 * b1 - a1 * b2 + 50.0 * c1 - 11.0 * c2 + 3.0 * c3 - c4);
        double d3 = 2.0 * a1 * a1 - b1 - 1.0 / 3.0 * (88.0 * a1 * a1 * a1 + 3.0 * a1 * (5.0 * b1 - 2.0 * b2) + 3.0 * (35.0 * c1 - 6.0 * c2 + c3));
        double d4 = -20.0 / 3.0 * a1 * a1 * a1 + a1 * (-4.0 * b1 + b2) - 10.0 * c1 + c2;

        double z1 = d2 - d3 + d4;
        double z2 = d3 - d4;
        double z3 = d4;

        double mean = momentsLevy[0] - 0.5 * momentsLevy[1];
        double variance = momentsLevy[1];

        double y = Math.log(strike);
        double std1 = normalDensity(y, mean, variance);
        double std2 = std1 * (mean - y) / variance;
        double std3 = std1 * ((mean - y) * (mean - y) / variance / variance - 1.0 / variance);

        return strike * (z1 * std1 + z2 * std2 + z3 * std3);
    }

    private static double normalDensity(double evaluationPoint, double mean, double variance) {
        double y = evaluationPoint - mean;
        return FastMath.exp(- 0.5 * y * y / variance) / FastMath.sqrt(2.0 * Math.PI * variance);
    }

    private static double[] sums(double[] processes, double[][] covariances) {
        double[] sums = new double[5];
        for(int i = 0; i < processes.length; i++) {
            for(int j = 0; j < processes.length; j++) {
                for(int k = 0; k < processes.length; k++) {
                    double productIJK = processes[i] * processes[j] * processes[k];
                    sums[0] += productIJK * covariances[i][k] * covariances[j][k];
                    for(int l = 0; l < processes.length; l++) {
                        sums[1] += productIJK * processes[l] * covariances[i][l] * covariances[j][k] * covariances[k][l];
                        sums[2] += productIJK * processes[l] * covariances[i][l] * covariances[j][l] * covariances[k][l];
                    }
                    sums[3] += productIJK * covariances[i][k] * covariances[j][k] * covariances[j][k];
                    sums[4] += productIJK * covariances[i][j] * covariances[i][k] * covariances[j][k];
                }
            }
        }
        return sums;
    }


    private double[] calculateIntegratedCovarianceHoSchwartzYeh(MultiCurveLIBORMarketModel model) {
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

        return getMomentsHo(mu, integratedCovariance);
    }

    public static double[] getMomentsHo(double[] mu, double[] integratedCovariance) {
        //sort them such that the bigger value is mu[0]
        if (Math.abs(mu[1]) > Math.abs(mu[0])) {
            double temp = mu[1];
            mu[1] = mu[0];
            mu[0] = temp;

            temp = integratedCovariance[1];
            integratedCovariance[1] = integratedCovariance[0];
            integratedCovariance[0] = temp;
        }

        double mw = mu[1] - mu[0];
        double variancew = integratedCovariance[0] + integratedCovariance[1] - 2.0 * integratedCovariance[2];
        double sigmaw = Math.sqrt(variancew);
        double eta = -mw / sigmaw;
        double sqrtTwoPi = Math.sqrt(2.0 * Math.PI);

        DoubleUnaryOperator fw = (w) -> Math.exp(-(w - mw) * (w - mw) / (2.0 * variancew)) / Math.sqrt(Math.PI * variancew);

        double i0 = integrate((v) -> Math.exp(-0.5 * (Math.log(v) - eta) * (Math.log(v) - eta)) / sqrtTwoPi);
        double i1 = integrate((v) -> Math.log(1.0 + v) * (fw.applyAsDouble(Math.log(v)) + fw.applyAsDouble(-Math.log(v))));
        double i2 = integrate((v) -> (fw.applyAsDouble(Math.log(v)) - fw.applyAsDouble(-Math.log(v))) / (1.0 + 1.0 / v));
        double i3 = integrate((v) -> {
            double log = Math.log(1.0 + v);
            return log * log * (fw.applyAsDouble(Math.log(v)) + fw.applyAsDouble(-Math.log(v)));
        });
        double i4 = integrate((v) -> -Math.log(1.0 + v) * Math.log(v) * fw.applyAsDouble(-Math.log(v)));

        double a0 = sigmaw / sqrtTwoPi * Math.exp(-0.5 * eta * eta) + mw * i0;
        double g1 = a0 + i1;
        double g2 = i3 + 2.0 * i4 + variancew * i0 + mw * a0;

        double mean = mu[0] + g1;
        double variance = integratedCovariance[0] - g1 * g1 - 2.0 * integratedCovariance[0] * (i2 + i0) + g2;
        return new double[] { mean, variance, 0.0 };
    }

    private static double integrate(DoubleUnaryOperator function) {
        return new SimpsonRealIntegrator(-3.3, 3.3, 1000, true).integrate(x -> {
            double u = Math.exp(- 2.0 * Math.sinh(x));
            double v = 1.0 / (1.0 + u);
            return 2.0 * function.applyAsDouble(v) * v * u * Math.cosh(x);
        });
    }
}
