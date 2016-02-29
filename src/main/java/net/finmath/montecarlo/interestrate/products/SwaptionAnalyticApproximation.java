/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 17.05.2007
 * Created on 30.03.2014
 */
package net.finmath.montecarlo.interestrate.products;

import java.lang.ref.WeakReference;
import java.util.HashMap;
import java.util.Map;

import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.*;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.*;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

/**
 * This class implements an analytic swaption valuation formula under
 * a LIBOR market model. The algorithm implemented here is the
 * OIS discounting version of the algorithm described in
 * ISBN 0470047224 (see {@link net.finmath.montecarlo.interestrate.products.SwaptionSingleCurveAnalyticApproximation}).
 *
 * The approximation assumes that the forward rates (LIBOR) follow a
 * log normal model and that the model provides the integrated
 * instantaneous covariance of the log-forward rates.
 *
 * The getValue method calculates the approximated integrated instantaneous variance of the swap rate,
 * using the approximation
 * \[
 * 	\frac{d log(S(t))}{d log(L(t))} \approx \frac{d log(S(0))}{d log(L(0))} = : w.
 * \]
 *
 * Since \( L \) is a vector, \( w \) is a gradient (vector). The class then approximates
 * the Black volatility of a swaption via
 * \[
 * 	\sigma_S^{2} T := \sum_{i,j} w_{i} \gamma_{i,j} w_{j}
 * \]
 * where \( (\gamma_{i,j})_{i,j = 1,...,m} \) is the covariance matrix of the forward rates.
 *
 * The valuation can be performed in terms of value or implied Black volatility.
 *
 *
 * @author Christian Fries
 * @date 17.05.2007.
 */
public class SwaptionAnalyticApproximation extends AbstractLIBORMonteCarloProduct {

    public enum MultiCurveApproximation {
        SINGLE_CURVE,
        DRIFT_FREEZE,
        INTEGRATED_EXPECTATION,
    }

    private final double      swaprate;
    private final double[]    swapTenor;       // Vector of swap tenor (period start and end dates). Start of first period is the option maturity.
    private ValueUnit   valueUnit;
    private MultiCurveApproximation multiCurveApproximation;

    private Map<String, double[]>						cachedLogSwaprateDerivative;
    private WeakReference<TimeDiscretizationInterface>	cachedLogSwaprateDerivativeTimeDiscretization;
    private WeakReference<DiscountCurveInterface>		cachedLogSwaprateDerivativeDiscountCurve;
    private WeakReference<ForwardCurveInterface>		cachedLogSwaprateDerivativeForwardCurve;
    private Object					cachedLogSwaprateDerivativeLock = new Object();

    /**
     * Create an analytic swaption approximation product for
     * log normal forward rate model.
     *
     * Note: It is implicitly assumed that swapTenor.getTime(0) is the exercise date (no forward starting).
     *
     * @param swaprate The strike swap rate of the swaption.
     * @param swapTenor The swap tenor in doubles.
     */
    public SwaptionAnalyticApproximation(double swaprate, TimeDiscretizationInterface swapTenor) {
        this(swaprate, swapTenor.getAsDoubleArray(), ValueUnit.VALUE);
    }

    public SwaptionAnalyticApproximation(double swaprate, double[] swapTenor, ValueUnit valueUnit) {
        this(swaprate, swapTenor, valueUnit, MultiCurveApproximation.DRIFT_FREEZE);
    }

    /**
     * Create an analytic swaption approximation product for
     * log normal forward rate model.
     *
     * Note: It is implicitly assumed that swapTenor[0] is the exercise date (no forward starting).
     *
     * @param swaprate The strike swap rate of the swaption.
     * @param swapTenor The swap tenor in doubles.
     * @param valueUnit The unit of the quantity returned by the getValues method.
     */
    public SwaptionAnalyticApproximation(double swaprate, double[] swapTenor, ValueUnit valueUnit, MultiCurveApproximation multiCurveApproximation) {
        super();
        this.swaprate	= swaprate;
        this.swapTenor	= swapTenor;
        this.valueUnit	= valueUnit;
        this.multiCurveApproximation = multiCurveApproximation;
    }

    @Override
    public RandomVariableInterface getValue(double evaluationTime, LIBORModelMonteCarloSimulationInterface model) {
        return getValues(evaluationTime, model.getModel(), this.valueUnit, this.multiCurveApproximation);
    }

    /**
     * Calculates the approximated integrated instantaneous variance of the swap rate,
     * using the approximation d log(S(t))/d log(L(t)) = d log(S(0))/d log(L(0)).
     *
     * @param evaluationTime Time at which the product is evaluated.
     * @param model A model implementing the LIBORModelMonteCarloSimulationInterface
     * @return Depending on the value of value unit, the method returns either
     * the approximated integrated instantaneous variance of the swap rate (ValueUnit.INTEGRATEDVARIANCE)
     * or the value using the Black formula (ValueUnit.VALUE).
     * @TODO make initial values an arg and use evaluation time.
     */
    public RandomVariableInterface getValues(double evaluationTime, LIBORMarketModelInterface model) {
        return getValues(evaluationTime, model, this.valueUnit, this.multiCurveApproximation);
    }

    /**
     * Calculates the approximated integrated instantaneous variance of the swap rate,
     * using the approximation d log(S(t))/d log(L(t)) = d log(S(0))/d log(L(0)).
     *
     * @param evaluationTime Time at which the product is evaluated.
     * @param model A model implementing the LIBORModelMonteCarloSimulationInterface
     * @return Depending on the value of value unit, the method returns either
     * the approximated integrated instantaneous variance of the swap rate (ValueUnit.INTEGRATEDVARIANCE)
     * or the value using the Black formula (ValueUnit.VALUE).
     * @TODO make initial values an arg and use evaluation time.
     */
    public RandomVariableInterface getValues(double evaluationTime, LIBORMarketModelInterface model, ValueUnit valueUnit, MultiCurveApproximation multiCurveApproximation) {
        if(evaluationTime > 0) throw new RuntimeException("Forward start evaluation currently not supported.");

        double swapStart    = swapTenor[0];
        double swapEnd      = swapTenor[swapTenor.length-1];

        int swapStartIndex  = model.getLiborPeriodIndex(swapStart);
        int swapEndIndex    = model.getLiborPeriodIndex(swapEnd);
        int optionMaturityIndex = model.getCovarianceModel().getTimeDiscretization().getTimeIndex(swapStart)-1;

        double integratedSwapRateVariance;
        if (model instanceof MultiCurveLIBORMarketModel && !(multiCurveApproximation == MultiCurveApproximation.SINGLE_CURVE)) {
            integratedSwapRateVariance = calculateIntegratedSwapRateVarianceMultiCurve((MultiCurveLIBORMarketModel) model, swapStartIndex, swapEndIndex, optionMaturityIndex, multiCurveApproximation);
        } else {
            integratedSwapRateVariance = calculateIntegratedSwapRateVarianceSingleCurve(model, swapStartIndex, swapEndIndex, optionMaturityIndex);
        }

        // Return integratedSwapRateVariance if requested
        if(valueUnit == ValueUnit.INTEGRATEDVARIANCE) return new RandomVariable(evaluationTime, integratedSwapRateVariance);

        double volatility		= Math.sqrt(integratedSwapRateVariance / swapStart);

        //if(model instanceof MultiCurveLIBORMarketModel) System.out.println(volatility);

        // Return integratedSwapRateVariance if requested
        if(valueUnit == ValueUnit.VOLATILITY) return new RandomVariable(evaluationTime, volatility);

        // Use black formula for swaption to calculate the price
        TimeDiscretizationInterface tenor = new TimeDiscretization(swapTenor);

        double shift = model instanceof ShiftedLIBORMarketModelInterface ?
                net.finmath.marketdata.products.Swap.getShift(tenor, tenor, (ShiftedLIBORMarketModelInterface)model) : 0.0;

        double parSwaprate		= net.finmath.marketdata.products.Swap.getForwardSwapRate(tenor, tenor, model.getForwardRateCurve(), model.getDiscountCurve());
        double swapAnnuity      = net.finmath.marketdata.products.SwapAnnuity.getSwapAnnuity(tenor, model.getDiscountCurve());

        double optionMaturity	= swapStart;

        double valueSwaption = AnalyticFormulas.blackModelSwaptionValue(parSwaprate + shift, volatility, optionMaturity, swaprate + shift, swapAnnuity);
        return new RandomVariable(evaluationTime, valueSwaption);
    }

    private double calculateIntegratedSwapRateVarianceMultiCurve(MultiCurveLIBORMarketModel model, int swapStartIndex, int swapEndIndex, int optionMaturityIndex, MultiCurveApproximation multiCurveApproximation) {
        // Get the integrated libor covariance from the model
        double[][][][]	integratedLIBORCovariance = model.getIntegratedLIBORCovariances();

        TimeDiscretizationInterface liborPeriodDiscretization = model.getLiborPeriodDiscretization();

        TimeDiscretizationInterface tenor = new TimeDiscretization(swapTenor);
        double shift = net.finmath.marketdata.products.Swap.getShift(tenor, tenor, model);

        RandomVariableInterface[] realization = model.getInitialState();
        for(int i = 0; i < realization.length; i++) {
            realization[i] = model.applyStateSpaceTransform(i, realization[i]);
        }
        RandomVariableInterface[][] realizations = new RandomVariableInterface[optionMaturityIndex + 1][];
        double[] swapRates = new double[optionMaturityIndex + 1];
        double[][] swapCovarianceWeights = new double[optionMaturityIndex + 1][];
        if (multiCurveApproximation == MultiCurveApproximation.INTEGRATED_EXPECTATION) {
            evolveProcessesAndWeights(optionMaturityIndex, swapStartIndex, swapEndIndex, realization, model,
                    realizations, swapRates, swapCovarianceWeights);
        } else {
            realizations[0] = realization;
            swapRates[0] = swaprate + shift;
            swapCovarianceWeights[0] = getLogSwaprateDerivative(model).get("values");
        }

        int timeEndIndex = multiCurveApproximation == MultiCurveApproximation.INTEGRATED_EXPECTATION ? optionMaturityIndex : 0;

        // Calculate integrated swap rate covariance
        double integratedSwapRateVariance = 0.0;
        for (int timeIndex = 0; timeIndex <= timeEndIndex; timeIndex++) {
            for(int componentIndex1 = swapStartIndex; componentIndex1 < swapEndIndex; componentIndex1++) {
                for(int componentIndex2 = componentIndex1; componentIndex2 < swapEndIndex; componentIndex2++) {
                    double timeStep1 = liborPeriodDiscretization.getTimeStep(componentIndex1);
                    double timeStep2 = liborPeriodDiscretization.getTimeStep(componentIndex2);

                    double forwardOfComponent1 = realizations[timeIndex][componentIndex1].getAverage();
                    double forwardOfComponent2 = realizations[timeIndex][componentIndex2].getAverage();

                    double spread1 = 0.0, spread2 = 0.0;
                    if (model.isMultiCurve()) {
                        int numberOfComponents = model.getNumberOfComponents() / 2;
                        spread1 = realizations[timeIndex][componentIndex1 + numberOfComponents].getAverage();
                        spread2 = realizations[timeIndex][componentIndex2 + numberOfComponents].getAverage();
                    }

                    forwardOfComponent1 -= model.getShiftParameter(componentIndex1);
                    forwardOfComponent2 -= model.getShiftParameter(componentIndex2);
                    spread1 -= model.getShiftParameter(componentIndex1 + model.getNumberOfComponents() / 2);
                    spread2 -= model.getShiftParameter(componentIndex2 + model.getNumberOfComponents() / 2);

                    if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MMARTINGALE) {
                        double alpha = model.getCorrelationFactor();

                        double coeff1 = spread1 * Math.pow(forwardOfComponent1, alpha - 1.0);
                        double coeff2 = spread2 * Math.pow(forwardOfComponent2, alpha - 1.0);

                        double[] coefficients = new double[5];

                        coefficients[0] = 1.0;
                        coefficients[1] = timeStep1 * coeff1 * forwardOfComponent1; //spread1 * Math.pow(forwardOfComponent1, alpha);
                        coefficients[2] = timeStep2 * coeff2 * forwardOfComponent2; //spread2 * Math.pow(forwardOfComponent2, alpha);

                        coefficients[3] = coefficients[1] * coefficients[2] + coeff1 * coeff2 * alpha * alpha;

                        coefficients[1] += coeff1 * alpha;
                        coefficients[2] += coeff2 * alpha;

                        coefficients[4] = coeff1 * coeff2 * (
                                1 + timeStep1 * forwardOfComponent1 + timeStep2 * forwardOfComponent2 +
                                        timeStep1 * forwardOfComponent1 * timeStep2 * forwardOfComponent2
                        );

                        for (int i = 0; i < 5; i++) {
                            double increment;
                            if (multiCurveApproximation == MultiCurveApproximation.INTEGRATED_EXPECTATION) {
                                increment = integratedLIBORCovariance[i][timeIndex][componentIndex1][componentIndex2] -
                                        (timeIndex > 0 ? integratedLIBORCovariance[i][timeIndex - 1][componentIndex1][componentIndex2] : 0.0);
                            } else {
                                increment = integratedLIBORCovariance[i][optionMaturityIndex][componentIndex1][componentIndex2];
                            }

                            integratedSwapRateVariance += (componentIndex1 == componentIndex2 ? 1.0 : 2.0) *
                                    forwardOfComponent1 * forwardOfComponent2 *
                                    swapCovarianceWeights[timeIndex][componentIndex1-swapStartIndex] * swapCovarianceWeights[timeIndex][componentIndex2-swapStartIndex] *
                                    coefficients[i] * increment;
                        }
                    } else {
                        double coeff1 = 1.0, coeff2 = 1.0, coeff3 = 1.0;
                        if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
                            coeff1 = (1.0 + timeStep1 * spread1) * (1.0 + timeStep2 * spread2);
                            coeff2 = (1.0 + timeStep1 * forwardOfComponent1) * (1.0 + timeStep2 * forwardOfComponent2);
                            coeff3 = (1.0 + timeStep1 * spread1) * (1.0 + timeStep2 * forwardOfComponent2);
                        }

                        double[] increment = new double[3];
                        for (int i = 0; i < 3; i++) {
                            if (multiCurveApproximation == MultiCurveApproximation.INTEGRATED_EXPECTATION) {
                                increment[i] = integratedLIBORCovariance[i][timeIndex][componentIndex1][componentIndex2] -
                                        (timeIndex > 0 ? integratedLIBORCovariance[i][timeIndex - 1][componentIndex1][componentIndex2] : 0.0);
                            } else {
                                increment[i] = integratedLIBORCovariance[i][optionMaturityIndex][componentIndex1][componentIndex2];
                            }
                        }

                        forwardOfComponent1 += model.getShiftParameter(componentIndex1);
                        forwardOfComponent2 += model.getShiftParameter(componentIndex2);
                        spread1 += model.getShiftParameter(componentIndex1 + model.getNumberOfComponents() / 2);
                        spread2 += model.getShiftParameter(componentIndex2 + model.getNumberOfComponents() / 2);

                        integratedSwapRateVariance += (componentIndex1 == componentIndex2 ? 1.0 : 2.0)
                                * swapCovarianceWeights[timeIndex][componentIndex1-swapStartIndex] * swapCovarianceWeights[timeIndex][componentIndex2-swapStartIndex]
                                * (coeff1 * forwardOfComponent1 * forwardOfComponent2 * increment[0]
                                + coeff2 * spread1 * spread2 * increment[1]
                                + coeff3 * forwardOfComponent1 * spread2 * increment[2]);
                    }
                }
            }
        }

        return integratedSwapRateVariance;
    }

    private void evolveProcessesAndWeights(int optionMaturityIndex, int swapStartIndex, int swapEndIndex,
                                           RandomVariableInterface[] realization, MultiCurveLIBORMarketModel model,
                                           RandomVariableInterface[][] realizations, double[] swapRates, double[][] swapCovarianceWeights) {
        TimeDiscretization tenor = new TimeDiscretization(swapTenor);
        double shift = net.finmath.marketdata.products.Swap.getShift(tenor, tenor, model);

        for (int timeIndex = 0; timeIndex <= optionMaturityIndex; timeIndex++) {
            realizations[timeIndex] = new RandomVariableInterface[model.getNumberOfComponents()];
            double dt = model.getTime(timeIndex + 1) - model.getTime(timeIndex);

            TimeDiscretizationInterface liborPeriodDiscretization = model.getLiborPeriodDiscretization();

            double[] forwardsAtTimeIndex = new double[liborPeriodDiscretization.getNumberOfTimes()];
            double[] liborsAtTimeIndex = new double[liborPeriodDiscretization.getNumberOfTimes()];
            for (int k = 0; k < liborPeriodDiscretization.getNumberOfTimeSteps(); k++) {
                forwardsAtTimeIndex[k] = realization[k].getAverage() - model.getShiftParameter(k);
                double spread = 0.0;
                if (model.isMultiCurve()) {
                    int numberOfComponents = model.getNumberOfComponents() / 2;
                    spread = realization[k + numberOfComponents].getAverage() - model.getShiftParameter(k + numberOfComponents);
                }
                liborsAtTimeIndex[k] = forwardsAtTimeIndex[k] + spread
                        + (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.ADDITIVE ?
                        0.0 : liborPeriodDiscretization.getTimeStep(k) * forwardsAtTimeIndex[k] * spread);
            }

            ForwardCurve liborCurve = ForwardCurve.createForwardCurveFromForwards(
                    model.getForwardRateCurve().getName(),
                    model.getDiscountCurve().getName(),
                    liborPeriodDiscretization.getAsDoubleArray(),
                    liborsAtTimeIndex,
                    liborPeriodDiscretization.getTimeStep(0)
            );
            ForwardCurveInterface forwardCurve = ForwardCurve.createForwardCurveFromForwards(
                    model.getRiskFreeCurve().getName(),
                    model.getDiscountCurve().getName(),
                    liborPeriodDiscretization.getAsDoubleArray(),
                    forwardsAtTimeIndex,
                    liborPeriodDiscretization.getTimeStep(0)
            );
            DiscountCurveInterface discountCurve = new DiscountCurveFromForwardCurve(model.getDiscountCurve().getName(), forwardCurve);

            double shiftedSwapRate	= net.finmath.marketdata.products.Swap.getForwardSwapRate(tenor, tenor, liborCurve, discountCurve) + shift;
            swapRates[timeIndex] = shiftedSwapRate;

            Map<String, double[]> logSwapRateDerivative = getLogSwaprateDerivative(liborPeriodDiscretization, discountCurve, liborCurve, model);
            swapCovarianceWeights[timeIndex] = logSwapRateDerivative.get("values");
            double[] swapAnnuities = logSwapRateDerivative.get("swapAnnuities");

            double[] drift = getDriftUnderSwapMeasure(timeIndex, swapStartIndex, swapEndIndex, realization, model, discountCurve, swapAnnuities[0]);
            for (int component = 0; component < (model.isMultiCurve() ? model.getNumberOfComponents() : model.getNumberOfLibors()); component++) {
                realizations[timeIndex][component] = realization[component];
                realization[component] = realization[component].mult(Math.exp((drift[component] - 0.5 * model.getCovarianceModel().getCovariance(timeIndex, component, component, realization).getAverage()) * dt));
            }
        }
    }

    private double[] getDriftUnderSwapMeasure(int timeIndex, int swapStartIndex, int swapEndIndex,
                                              RandomVariableInterface[] realizationAtTimeIndex,
                                              MultiCurveLIBORMarketModel model,
                                              DiscountCurveInterface discountCurve, double annuity) {
        double[] drift = new double[model.getNumberOfComponents()];

        TimeDiscretizationInterface liborPeriodDiscretization = model.getLiborPeriodDiscretization();
        int		firstLiborIndex		= model.getLiborPeriodIndex(model.getTime(timeIndex)) + 1;
        if(firstLiborIndex<0) firstLiborIndex = -firstLiborIndex;


        for (int k = Math.max(firstLiborIndex, swapStartIndex); k < swapEndIndex; k++) {
            for (int j = swapStartIndex; j < swapEndIndex; j++) {
                double timeStepJ = liborPeriodDiscretization.getTimeStep(j);
                int start = Math.min(j+1, k+1);
                int end = Math.max(k, j);
                for (int i = start; i <= end; i++) {
                    double timeStepI = liborPeriodDiscretization.getTimeStep(i);
                    double shift = model.getShiftParameter(i);
                    double forwardI = realizationAtTimeIndex[i].getAverage() - shift;
                    double discountedForward = timeStepI * (forwardI + shift) / (1.0 * timeStepI * forwardI);
                    RandomVariableInterface covarianceForward = model.getCovarianceModel().getCovariance(timeIndex, k, i, realizationAtTimeIndex);

                    drift[k] += covarianceForward.mult(discountedForward).getAverage();

                    if (model.isMultiCurve()) {
                        int numberOfComponents = model.getNumberOfComponents() / 2;
                        RandomVariableInterface covarianceSpread = model.getCovarianceModel().getCovariance(timeIndex, k + numberOfComponents, i, realizationAtTimeIndex);
                        drift[k + numberOfComponents] += covarianceSpread.mult(discountedForward).getAverage();
                    }
                }
                double discount = (k >= j ? 1.0 : -1.0) * timeStepJ * discountCurve.getDiscountFactor(null, liborPeriodDiscretization.getTime(j + 1)) / annuity;
                drift[k] = drift[k] * discount;
                if (model.isMultiCurve()) {
                    int numberOfComponents = model.getNumberOfComponents() / 2;
                    drift[k + numberOfComponents] = drift[k + numberOfComponents] * discount;
                }
            }

            //Adjust the drift in the multiplicative model
            if (model.getMultiCurveModel() == MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE) {
                double timeStep = liborPeriodDiscretization.getTimeStep(k);
                double shiftedForwardK = realizationAtTimeIndex[k].getAverage();
                double forward = shiftedForwardK - model.getShiftParameter(k);

                if (model.isMultiCurve()) {
                    int numberOfComponents = model.getNumberOfComponents() / 2;
                    double shiftedSpreadK = realizationAtTimeIndex[k + numberOfComponents].getAverage();
                    drift[k + numberOfComponents] -= timeStep / (1 + timeStep * forward) * shiftedForwardK * shiftedSpreadK *
                            model.getCovarianceModel().getCovariance(timeIndex, k, k + numberOfComponents, realizationAtTimeIndex).getAverage();
                }
            }
        }

        return drift;
    }

    private double calculateIntegratedSwapRateVarianceSingleCurve(LIBORMarketModelInterface model, int swapStartIndex, int swapEndIndex, int optionMaturityIndex) {
        Map<String, double[]>  logSwaprateDerivative  = getLogSwaprateDerivative(model);
        double[] swapCovarianceWeights = logSwaprateDerivative.get("values");

        // Get the integrated libor covariance from the model
        double[][]	integratedLIBORCovariance = model.getIntegratedLIBORCovariance()[optionMaturityIndex];

        // Calculate integrated swap rate covariance
        double integratedSwapRateVariance = 0.0;
        for(int componentIndex1 = swapStartIndex; componentIndex1 < swapEndIndex; componentIndex1++) {
            // Sum the libor cross terms (use symmetry)
            for(int componentIndex2 = componentIndex1+1; componentIndex2 < swapEndIndex; componentIndex2++) {
                double forward1 = model.getForwardRateCurve().getForward(model.getAnalyticModel(), model.getLiborPeriodDiscretization().getTime(componentIndex1));
                double forward2 = model.getForwardRateCurve().getForward(model.getAnalyticModel(), model.getLiborPeriodDiscretization().getTime(componentIndex2));

                if (model instanceof ShiftedLIBORMarketModelInterface) {
                    forward1 += ((ShiftedLIBORMarketModelInterface) model).getShiftParameter(componentIndex1);
                    forward2 += ((ShiftedLIBORMarketModelInterface) model).getShiftParameter(componentIndex2);
                }

                integratedSwapRateVariance += (componentIndex1 == componentIndex2 ? 1.0 : 2.0) *
                        swapCovarianceWeights[componentIndex1-swapStartIndex] * swapCovarianceWeights[componentIndex2-swapStartIndex]
                        * forward1 * forward2 * integratedLIBORCovariance[componentIndex1][componentIndex2];
            }
        }
        return integratedSwapRateVariance;
    }

    /**
     * This function calculate the partial derivative <i>d log(S) / d log(L<sub>k</sub>)</i> for
     * a given swap rate with respect to a vector of forward rates (on a given forward rate tenor).
     *
     * It also returns some useful other quantities like the corresponding discout factors and swap annuities.
     *
     * @param model The LIBOR market model.
     * @return A map containing the partial derivatives (key "value"), the discount factors (key "discountFactors") and the annuities (key "annuities") as vectors of double[] (indexed by forward rate tenor index starting at swap start)
     */
    public Map<String, double[]> getLogSwaprateDerivative(LIBORMarketModelInterface model) {
        ForwardCurveInterface forwardCurve = model.getForwardRateCurve();
        DiscountCurveInterface discountCurve = model.getDiscountCurve();
        TimeDiscretizationInterface liborPeriodDiscretization = model.getLiborPeriodDiscretization();

        return getLogSwaprateDerivative(liborPeriodDiscretization, discountCurve, forwardCurve, model);
    }

    /**
     * This function calculate the partial derivative <i>d log(S) / d log(L<sub>k</sub>)</i> for
     * a given swap rate with respect to a vector of forward rates (on a given forward rate tenor).
     *
     * It also returns some useful other quantities like the corresponding discount factors and swap annuities.
     *
     * @param liborPeriodDiscretization The libor period discretization.
     * @param discountCurve The discount curve. If this parameter is null, the discount curve will be calculated from the forward curve.
     * @param forwardCurve The forward curve.
     * @param model The LIBOR market model.
     * @return A map containing the partial derivatives (key "value"), the discount factors (key "discountFactors") and the annuities (key "annuities") as vectors of double[] (indexed by forward rate tenor index starting at swap start)
     */
    public Map<String, double[]> getLogSwaprateDerivative(TimeDiscretizationInterface liborPeriodDiscretization, DiscountCurveInterface discountCurve, ForwardCurveInterface forwardCurve, LIBORMarketModelInterface model) {

        /*
		 * We cache the calculation of the log swaprate derivative. In a calibration this method might be called quite often with the same arguments.
		 */
        synchronized (cachedLogSwaprateDerivativeLock) {
            if(
                    cachedLogSwaprateDerivative != null
                            && liborPeriodDiscretization == cachedLogSwaprateDerivativeTimeDiscretization.get()
                            && discountCurve == cachedLogSwaprateDerivativeDiscountCurve.get()
                            && forwardCurve == cachedLogSwaprateDerivativeForwardCurve.get()
                    ) {
                return cachedLogSwaprateDerivative;
            }
            cachedLogSwaprateDerivativeTimeDiscretization = new WeakReference<>(liborPeriodDiscretization);
            cachedLogSwaprateDerivativeDiscountCurve = new WeakReference<>(discountCurve);
            cachedLogSwaprateDerivativeForwardCurve = new WeakReference<>(forwardCurve);

            AnalyticModelInterface analyticModel = model.getAnalyticModel();
            if(discountCurve == null) {
                discountCurve	= new DiscountCurveFromForwardCurve(forwardCurve.getName());
                analyticModel					= new AnalyticModel(new CurveInterface[] { forwardCurve, discountCurve });
            }

            double swapStart    = swapTenor[0];
            double swapEnd      = swapTenor[swapTenor.length-1];

            // Get the indices of the swap start and end on the forward rate tenor
            int swapStartIndex  = liborPeriodDiscretization.getTimeIndex(swapStart);
            int swapEndIndex    = liborPeriodDiscretization.getTimeIndex(swapEnd);

            // Precalculate forward rates and discount factors. Note: the swap contains swapEndIndex-swapStartIndex forward rates
            double[] forwardRates       = new double[swapEndIndex-swapStartIndex+1];
            double[] shiftedForwardRates       = new double[swapEndIndex-swapStartIndex+1];
            double[] discountFactors    = new double[swapEndIndex-swapStartIndex+1];

            // Calculate discount factor at swap start
            discountFactors[0] = discountCurve.getDiscountFactor(analyticModel, swapStart);

            // Calculate discount factors for swap period ends (used for swap annuity)
            TimeDiscretizationInterface tenor = new TimeDiscretization(swapTenor);
            for(int liborPeriodIndex = swapStartIndex; liborPeriodIndex < swapEndIndex; liborPeriodIndex++) {
                double libor = forwardCurve.getForward(analyticModel, liborPeriodDiscretization.getTime(liborPeriodIndex));
                double shift = model instanceof ShiftedLIBORMarketModelInterface ?
                        net.finmath.marketdata.products.Swap.getShift(tenor, tenor, (ShiftedLIBORMarketModelInterface)model) : 0.0;

                forwardRates[liborPeriodIndex-swapStartIndex]       = libor;
                shiftedForwardRates[liborPeriodIndex-swapStartIndex]= libor + shift;
                discountFactors[liborPeriodIndex-swapStartIndex+1]  = discountCurve.getDiscountFactor(analyticModel, liborPeriodDiscretization.getTime(liborPeriodIndex+1));
            }

            // Precalculate swap annuities
            double[]    swapAnnuities   = new double[swapTenor.length-1];
            double      swapAnnuity     = 0.0;
            for(int swapPeriodIndex = swapTenor.length-2; swapPeriodIndex >= 0; swapPeriodIndex--) {
                int periodEndIndex = liborPeriodDiscretization.getTimeIndex(swapTenor[swapPeriodIndex+1]);
                swapAnnuity += discountFactors[periodEndIndex-swapStartIndex] * (swapTenor[swapPeriodIndex+1]-swapTenor[swapPeriodIndex]);
                swapAnnuities[swapPeriodIndex] = swapAnnuity;
            }

            // Precalculate weights: The formula is take from ISBN 0470047224
            double[] swapCovarianceWeights = new double[swapEndIndex-swapStartIndex];

            double valueFloatLeg = 0.0;
            for(int liborPeriodIndex = swapStartIndex; liborPeriodIndex < swapEndIndex; liborPeriodIndex++) {
                double liborPeriodLength = liborPeriodDiscretization.getTimeStep(liborPeriodIndex);
                valueFloatLeg += shiftedForwardRates[liborPeriodIndex-swapStartIndex] * discountFactors[liborPeriodIndex-swapStartIndex+1] * liborPeriodLength;
            }

            int swapPeriodIndex = 0;
            double valueFloatLegUpToSwapStart = 0.0;
            for(int liborPeriodIndex = swapStartIndex; liborPeriodIndex < swapEndIndex; liborPeriodIndex++) {
                if(liborPeriodDiscretization.getTime(liborPeriodIndex) >= swapTenor[swapPeriodIndex+1]) swapPeriodIndex++;

                double libor				= forwardRates[liborPeriodIndex-swapStartIndex];
                double shiftedLibor			= shiftedForwardRates[liborPeriodIndex-swapStartIndex];
                double liborPeriodLength	= liborPeriodDiscretization.getTimeStep(liborPeriodIndex);

                valueFloatLegUpToSwapStart += shiftedLibor * discountFactors[liborPeriodIndex-swapStartIndex+1] * liborPeriodLength;

                double discountFactorAtPeriodEnd = discountCurve.getDiscountFactor(analyticModel, liborPeriodDiscretization.getTime(liborPeriodIndex+1));
                double derivativeFloatLeg	= (discountFactorAtPeriodEnd + (valueFloatLegUpToSwapStart - valueFloatLeg)  / (1.0 + libor * liborPeriodLength)) * liborPeriodLength / valueFloatLeg;
                double derivativeFixLeg		= swapAnnuities[swapPeriodIndex] / swapAnnuity * liborPeriodLength / (1.0 + libor * liborPeriodLength);

                swapCovarianceWeights[liborPeriodIndex-swapStartIndex] = derivativeFloatLeg + derivativeFixLeg;
            }

            // Return results
            Map<String, double[]> results = new HashMap<>();
            results.put("values",			swapCovarianceWeights);
            results.put("discountFactors",	discountFactors);
            results.put("swapAnnuities",	swapAnnuities);

            cachedLogSwaprateDerivative = results;

            return results;
        }
    }

    static public double[][][] getIntegratedLIBORCovariance(LIBORMarketModel model) {
        return model.getIntegratedLIBORCovariance();
    }
}
