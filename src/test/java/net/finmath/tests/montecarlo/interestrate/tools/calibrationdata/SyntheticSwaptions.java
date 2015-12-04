package net.finmath.tests.montecarlo.interestrate.tools.calibrationdata;

import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.*;
import net.finmath.marketdata.products.SwapAnnuity;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.SwaptionAnalyticApproximation;
import net.finmath.montecarlo.interestrate.products.SwaptionFactory;
import net.finmath.tests.montecarlo.interestrate.tools.CalibrationData;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.ArrayList;
import java.util.Arrays;

public class SyntheticSwaptions extends CalibrationData {
    public SyntheticSwaptions(TimeDiscretizationInterface timeDiscretization,
                              TimeDiscretizationInterface periodDiscretization,
                              AnalyticModelInterface analyticModel,
                              String secondCurveName, String firstCurveName,
                              LIBORMarketModelInterface.CalibrationItem[] firstCurveCalibrationItems,
                              LIBORMarketModelInterface.CalibrationItem[] secondCurveNameCalibrationItems,
                              LIBORMarketModelInterface.CalibrationItem[] firstCurveCalibrationItemsMonteCarlo,
                              LIBORMarketModelInterface.CalibrationItem[] secondCurveNameCalibrationItemsMonteCarlo) {
        super(timeDiscretization,
                periodDiscretization,
                analyticModel,
                secondCurveName,
                firstCurveName,
                firstCurveCalibrationItems,
                secondCurveNameCalibrationItems,
                firstCurveCalibrationItemsMonteCarlo,
                secondCurveNameCalibrationItemsMonteCarlo);
    }

    public static CalibrationData newInstance() {
        /*
		 * Create the libor tenor structure and the initial values
		 */
        double liborPeriodLength	= 0.5;
        double liborRateTimeHorzion	= 20.0;
        TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

        // Create the forward curve (initial value of the LIBOR market model)
        ForwardCurve euriborForwardCurve = ForwardCurve.createForwardCurveFromForwards(
                "EURIBOR"								/* name of the curve */,
                new double[]{0.5, 1.0, 5.0, 10.0, 40.0}	/* fixings of the forward */,
                new double[]{0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
                liborPeriodLength							/* tenor / period length */
        );

        // Create the discount curve

        ForwardCurveInterface eoniaForwardCurve = ForwardCurve.createForwardCurveFromForwards(
                "EONIA"								/* name of the curve */,
                "DiscountCurve",
                new double[]{0.5, 1.0, 5.0, 10.0, 40.0}	/* maturities */,
                new double[]{0.04, 0.04, 0.04, 0.04, 0.04}	/* forwards */,
                0.5
        );
        DiscountCurveInterface discountCurve = new DiscountCurveFromForwardCurve("DiscountCurve", eoniaForwardCurve);

        AnalyticModelInterface analyticModel = new AnalyticModel(new CurveInterface[] { discountCurve, eoniaForwardCurve, euriborForwardCurve });

		/*
		 * Create a set of calibration products.
		 */
        ArrayList<LIBORMarketModelInterface.CalibrationItem> eoniaCalibrationItems = new ArrayList<>();
        ArrayList<LIBORMarketModelInterface.CalibrationItem> euriborCalibrationItems = new ArrayList<>();
        ArrayList<LIBORMarketModelInterface.CalibrationItem> euriborCalibrationItemsMonteCarlo = new ArrayList<>();
        ArrayList<LIBORMarketModelInterface.CalibrationItem> eoniaCalibrationItemsMonteCarlo = new ArrayList<>();
        for (int exerciseIndex = 4; exerciseIndex <= liborPeriodDiscretization.getNumberOfTimeSteps() - 5; exerciseIndex+=4) {
            double exerciseDate = liborPeriodDiscretization.getTime(exerciseIndex);
            for (int numberOfPeriods = 1; numberOfPeriods < liborPeriodDiscretization.getNumberOfTimeSteps() - exerciseIndex - 5; numberOfPeriods+=4) {

                // Create swaptions

                double[]	swapTenor			= new double[numberOfPeriods + 1];
                double		swapPeriodLength	= 0.5;

                for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
                    swapTenor[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
                }
                swapTenor[numberOfPeriods] = exerciseDate + numberOfPeriods * swapPeriodLength;

                TimeDiscretizationInterface fixAndFloatTenor = new TimeDiscretization(swapTenor);

                // EONIA Swaptions swap rate
                double eoniaSwaprate = net.finmath.marketdata.products.Swap.getForwardSwapRate(
                        fixAndFloatTenor,
                        fixAndFloatTenor,
                        eoniaForwardCurve,
                        discountCurve
                );

                // EURIBOR Swaptions swap rate
                double euriborSwaprate = net.finmath.marketdata.products.Swap.getForwardSwapRate(
                        fixAndFloatTenor,
                        fixAndFloatTenor,
                        euriborForwardCurve,
                        discountCurve
                );

                // Set swap rates for each period
                double[] eoniaSwaprates = new double[numberOfPeriods];
                Arrays.fill(eoniaSwaprates, eoniaSwaprate);

                double[] euriborSwaprates = new double[numberOfPeriods];
                Arrays.fill(euriborSwaprates, euriborSwaprate);

                double swapAnnuity = SwapAnnuity.getSwapAnnuity(new TimeDiscretization(swapTenor), discountCurve);

                // This is just some swaption volatility used for testing, true market data should go here.
                double eoniaTargetValueVolatilty = 0.20 + 0.20 * Math.exp(-exerciseDate / 10.0) + 0.20 * Math.exp(-(exerciseDate+numberOfPeriods) / 10.0);
                double euriborTargetValueVolatilty = 0.20 + 0.20 * Math.exp(-exerciseDate / 10.0) + 0.20 * Math.exp(-(exerciseDate+numberOfPeriods) / 10.0);

                // Build our calibration products
                SwaptionAnalyticApproximation eoniaSwaptionAnalytic = new SwaptionAnalyticApproximation(
                        eoniaSwaprate, swapTenor, SwaptionAnalyticApproximation.ValueUnit.VOLATILITY
                );

                eoniaCalibrationItems.add(new LIBORMarketModelInterface.CalibrationItem(eoniaSwaptionAnalytic, eoniaTargetValueVolatilty, 1.0));

                SwaptionAnalyticApproximation euriborSwaptionAnalytic = new SwaptionAnalyticApproximation(
                        euriborSwaprate, swapTenor, SwaptionAnalyticApproximation.ValueUnit.VOLATILITY
                );

                euriborCalibrationItems.add(new LIBORMarketModelInterface.CalibrationItem(euriborSwaptionAnalytic, euriborTargetValueVolatilty, 1.0));

                AbstractLIBORMonteCarloProduct swaptionMonteCarlo = SwaptionFactory.createSwaption("SwaptionSimple", euriborSwaprate, fixAndFloatTenor, "VALUE");
                double targetValuePrice = AnalyticFormulas.blackModelSwaptionValue(euriborSwaprate, euriborTargetValueVolatilty, exerciseDate, euriborSwaprate, swapAnnuity);
                euriborCalibrationItemsMonteCarlo.add(new LIBORMarketModelInterface.CalibrationItem(swaptionMonteCarlo, targetValuePrice, 1.0));

                swaptionMonteCarlo = SwaptionFactory.createSwaption("SwaptionSimple", eoniaSwaprate, fixAndFloatTenor, "VALUE");
                targetValuePrice = AnalyticFormulas.blackModelSwaptionValue(eoniaSwaprate, eoniaTargetValueVolatilty, exerciseDate, eoniaSwaprate, swapAnnuity);
                eoniaCalibrationItemsMonteCarlo.add(new LIBORMarketModelInterface.CalibrationItem(swaptionMonteCarlo, targetValuePrice, 1.0));
            }
        }

        return new CalibrationData(liborPeriodDiscretization,
                liborPeriodDiscretization,
                analyticModel,
                "EONIA",
                "EURIBOR",
                eoniaCalibrationItems.toArray(new LIBORMarketModelInterface.CalibrationItem[eoniaCalibrationItems.size()]),
                euriborCalibrationItems.toArray(new LIBORMarketModelInterface.CalibrationItem[euriborCalibrationItems.size()]),
                eoniaCalibrationItemsMonteCarlo.toArray(new LIBORMarketModelInterface.CalibrationItem[eoniaCalibrationItemsMonteCarlo.size()]),
                euriborCalibrationItemsMonteCarlo.toArray(new LIBORMarketModelInterface.CalibrationItem[euriborCalibrationItemsMonteCarlo.size()]));
    }
}
