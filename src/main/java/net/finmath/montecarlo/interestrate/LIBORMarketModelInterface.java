package net.finmath.montecarlo.interestrate;

import java.util.ArrayList;
import java.util.Map;

import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.DiscountCurveInterface;
import net.finmath.marketdata.model.curves.ForwardCurveInterface;
import net.finmath.marketdata.model.volatilities.AbstractSwaptionMarketData;
import net.finmath.marketdata.products.Swap;
import net.finmath.marketdata.products.SwapAnnuity;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.SwaptionAnalyticApproximation;
import net.finmath.montecarlo.interestrate.products.SwaptionSimple;
import net.finmath.montecarlo.model.AbstractModelInterface;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.RegularSchedule;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public interface LIBORMarketModelInterface extends AbstractModelInterface {

    RandomVariableInterface getLIBOR(int timeIndex, int liborIndex) throws CalculationException;

    /**
     * The tenor time discretization of the forward rate curve.
     *
     * @return The tenor time discretization of the forward rate curve.
     */
    TimeDiscretizationInterface getLiborPeriodDiscretization();

    /**
     * Get the number of LIBORs in the LIBOR discretization.
     *
     * @return The number of LIBORs in the LIBOR discretization
     */
    int getNumberOfLibors();

    /**
     * The period start corresponding to a given forward rate discretization index.
     *
     * @param timeIndex The index corresponding to a given time (interpretation is start of period)
     * @return The period start corresponding to a given forward rate discretization index.
     */
    double getLiborPeriod(int timeIndex);

    /**
     * Same as java.util.Arrays.binarySearch(liborPeriodDiscretization,time). Will return a negative value if the time is not found, but then -index-1 corresponds to the index of the smallest time greater than the given one.
     *
     * @param time The period start.
     * @return The index corresponding to a given time (interpretation is start of period)
     */
    int getLiborPeriodIndex(double time);

    /**
     * Return the associated analytic model, a collection of market date object like discount curve, forward curve
     * and volatility surfaces.
     *
     * @return The associated analytic model.
     */
    AnalyticModelInterface getAnalyticModel();

    /**
     * Return the discount curve associated the forwards.
     *
     * @return the discount curve associated the forwards.
     */
    DiscountCurveInterface getDiscountCurve();

    /**
     * Return the initial forward rate curve.
     *
     * @return the forward rate curve
     */
    ForwardCurveInterface getForwardRateCurve();

    /**
     * Return the covariance model.
     *
     * @return The covariance model.
     */
    AbstractLIBORCovarianceModel getCovarianceModel();

    /**
     * Create a new object implementing LIBORMarketModelInterface, using the new covariance model.
     *
     * @param calibrationCovarianceModel The new covariance model.
     * @return A new object implementing LIBORMarketModelInterface, using the new covariance model.
     */
    LIBORMarketModelInterface getCloneWithModifiedCovarianceModel(AbstractLIBORCovarianceModel calibrationCovarianceModel);

    /**
     * Create a new object implementing LIBORMarketModelInterface, using the new data.
     *
     * @param dataModified A map with values to be used in constructions (keys are identical to parameter names of the constructors).
     * @return A new object implementing LIBORMarketModelInterface, using the new data.
     * @throws net.finmath.exception.CalculationException Thrown if the valuation fails, specific cause may be available via the <code>cause()</code> method.
     */
    LIBORMarketModelInterface getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException;

    /**
     * Returns the integrated instantaneous log-forward rate covariance, i.e.,
     * <i>\int_0^t_i d log(L_j) d log(L_k) dt</i>.
     *
     * The array returned has the parametrization [i][j][k], i.e.,
     * <code>integratedLIBORCovariance[timeIndex][componentIndex1][componentIndex2]</code>.
     *
     * @return The integrated instantaneous log-LIBOR covariance.
     */
    double[][][] getIntegratedLIBORCovariance();

    public static class CalibrationItem {
        public final AbstractLIBORMonteCarloProduct calibrationProduct;
        public final double								calibrationTargetValue;
        public final double								calibrationWeight;

        public CalibrationItem(AbstractLIBORMonteCarloProduct calibrationProduct, double calibrationTargetValue, double calibrationWeight) {
            super();
            this.calibrationProduct		= calibrationProduct;
            this.calibrationTargetValue	= calibrationTargetValue;
            this.calibrationWeight		= calibrationWeight;
        }

        @Override
        public String toString() {
            return "CalibrationItem [calibrationProduct=" + calibrationProduct
                    + ", calibrationTargetValue=" + calibrationTargetValue
                    + ", calibrationWeight=" + calibrationWeight + "]";
        }
    }

    static CalibrationItem[] getCalibrationItems(TimeDiscretizationInterface liborPeriodDiscretization, ForwardCurveInterface forwardCurve, AbstractSwaptionMarketData swaptionMarketData, boolean isUseAnalyticApproximation) {
        if(swaptionMarketData == null) return null;

        TimeDiscretizationInterface	optionMaturities		= swaptionMarketData.getOptionMaturities();
        TimeDiscretizationInterface	tenor					= swaptionMarketData.getTenor();
        double						swapPeriodLength		= swaptionMarketData.getSwapPeriodLength();

        ArrayList<CalibrationItem> calibrationItems = new ArrayList<CalibrationItem>();
        for(int exerciseIndex=0; exerciseIndex<=optionMaturities.getNumberOfTimeSteps(); exerciseIndex++) {
            for(int tenorIndex=0; tenorIndex<=tenor.getNumberOfTimeSteps()-exerciseIndex; tenorIndex++) {

                // Create a swaption
                double exerciseDate	= optionMaturities.getTime(exerciseIndex);
                double swapLength	= tenor.getTime(tenorIndex);

                if(liborPeriodDiscretization.getTimeIndex(exerciseDate) < 0) continue;
                if(liborPeriodDiscretization.getTimeIndex(exerciseDate+swapLength) <= liborPeriodDiscretization.getTimeIndex(exerciseDate)) continue;

                int numberOfPeriods = (int)(swapLength / swapPeriodLength);

                double[] fixingDates      = new double[numberOfPeriods];
                double[] paymentDates     = new double[numberOfPeriods];
                double[] swapTenorTimes   = new double[numberOfPeriods+1];

                for(int periodStartIndex=0; periodStartIndex<numberOfPeriods; periodStartIndex++) {
                    fixingDates[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
                    paymentDates[periodStartIndex] = exerciseDate + (periodStartIndex+1) * swapPeriodLength;
                    swapTenorTimes[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
                }
                swapTenorTimes[numberOfPeriods] = exerciseDate + numberOfPeriods * swapPeriodLength;


                // Swaptions swap rate
                ScheduleInterface swapTenor = new RegularSchedule(new TimeDiscretization(swapTenorTimes));
                double swaprate = Swap.getForwardSwapRate(swapTenor, swapTenor, forwardCurve, null);

                // Set swap rates for each period
                double[] swaprates        = new double[numberOfPeriods];
                for(int periodStartIndex=0; periodStartIndex<numberOfPeriods; periodStartIndex++) {
                    swaprates[periodStartIndex] = swaprate;
                }

                if(isUseAnalyticApproximation) {
                    AbstractLIBORMonteCarloProduct swaption = new SwaptionAnalyticApproximation(swaprate, swapTenorTimes, SwaptionAnalyticApproximation.ValueUnit.VOLATILITY);
                    double impliedVolatility = swaptionMarketData.getVolatility(exerciseDate, swapLength, swaptionMarketData.getSwapPeriodLength(), swaprate);

                    calibrationItems.add(new CalibrationItem(swaption, impliedVolatility, 1.0));
                }
                else {
                    AbstractLIBORMonteCarloProduct swaption = new SwaptionSimple(swaprate, swapTenorTimes, SwaptionSimple.ValueUnit.VALUE);

                    double forwardSwaprate		= Swap.getForwardSwapRate(swapTenor, swapTenor, forwardCurve);
                    double swapAnnuity 			= SwapAnnuity.getSwapAnnuity(swapTenor, forwardCurve);
                    double impliedVolatility	= swaptionMarketData.getVolatility(exerciseDate, swapLength, swaptionMarketData.getSwapPeriodLength(), swaprate);

                    double targetValue = AnalyticFormulas.blackModelSwaptionValue(forwardSwaprate, impliedVolatility, exerciseDate, swaprate, swapAnnuity);

                    calibrationItems.add(new CalibrationItem(swaption, targetValue, 1.0));
                }
            }
        }

        return calibrationItems.toArray(new CalibrationItem[calibrationItems.size()]);
    }
}