package net.finmath.tests.montecarlo.interestrate.tools.calibrationdata;

import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.*;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.products.*;
import net.finmath.tests.montecarlo.interestrate.tools.CalibrationData;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

public class SyntheticCaplets extends CalibrationData {
    public SyntheticCaplets(TimeDiscretizationInterface timeDiscretization,
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
        double liborPeriodLength	= 0.5;
        double liborRateTimeHorzion	= 20.0;
        TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);

        // Create the forward curve (initial value of the LIBOR market model)
        ForwardCurve euriborForwardCurve = ForwardCurve.createForwardCurveFromForwards(
                "EURIBOR",
                "DiscountCurve"             /* name of the curve */,
                new double[]{0.5, 1.0, 5.0, 10.0, 40.0}	/* fixings of the forward */,
                new double[]{0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
                liborPeriodLength							/* tenor / period length */
        );

        // Create the discount curve

        ForwardCurveInterface eoniaForwardCurve = ForwardCurve.createForwardCurveFromForwards(
                "EONIA"								/* name of the curve */,
                "DiscountCurve",
                new double[]{0.5, 1.0, 5.0, 10.0, 40.0}	/* maturities */,
                new double[]{0.025, 0.025, 0.025, 0.025, 0.025}	/* forwards */,
                liborPeriodLength
        );
        DiscountCurveInterface discountCurve = new DiscountCurveFromForwardCurve("DiscountCurve", eoniaForwardCurve);

        AnalyticModelInterface analyticModel = new AnalyticModel(new CurveInterface[] { discountCurve, eoniaForwardCurve, euriborForwardCurve });

        int numberOfTimeSteps = liborPeriodDiscretization.getNumberOfTimeSteps();

        LIBORMarketModelInterface.CalibrationItem[] eoniaCalibrationItemsAnalytic = new LIBORMarketModelInterface.CalibrationItem[numberOfTimeSteps - 1];
        LIBORMarketModelInterface.CalibrationItem[] euriborCalibrationItemsAnalytic = new LIBORMarketModelInterface.CalibrationItem[numberOfTimeSteps - 1];

        LIBORMarketModelInterface.CalibrationItem[] eoniaCalibrationItemsMonteCarlo = new LIBORMarketModelInterface.CalibrationItem[numberOfTimeSteps - 1];
        LIBORMarketModelInterface.CalibrationItem[] euriborCalibrationItemsMonteCarlo = new LIBORMarketModelInterface.CalibrationItem[numberOfTimeSteps - 1];

        for (int component = 1; component < numberOfTimeSteps; component++) {
            double maturity = liborPeriodDiscretization.getTime(component);
            double periodLength = liborPeriodDiscretization.getTimeStep(component);

            double eoniaStrike = eoniaForwardCurve.getForward(null, maturity);
            AbstractLIBORMonteCarloProduct eoniaCaplet = new Caplet(maturity, periodLength, eoniaStrike, false);
            AbstractLIBORMonteCarloProduct eoniaCapletAnalytic = new CapletAnalyticApproximation(eoniaStrike, maturity, CapletAnalyticApproximation.ValueUnit.VOLATILITY);

            double euriborStrike = euriborForwardCurve.getForward(null, maturity);
            AbstractLIBORMonteCarloProduct euriborCaplet = new Caplet(maturity, periodLength, euriborStrike, false);
            AbstractLIBORMonteCarloProduct euriborCapletAnalytic = new CapletAnalyticApproximation(euriborStrike, maturity, CapletAnalyticApproximation.ValueUnit.VOLATILITY);

            double eoniaVolatility = 0.35 + (0.45 - 0.35) * Math.exp(-(double)component / numberOfTimeSteps);
            double euriborVolatility = 0.30 + (0.40 - 0.30) * Math.exp(-(double)component / numberOfTimeSteps);

            eoniaCalibrationItemsAnalytic[component - 1] = new LIBORMarketModelInterface.CalibrationItem(eoniaCapletAnalytic, eoniaVolatility, 1.0);
            euriborCalibrationItemsAnalytic[component - 1] = new LIBORMarketModelInterface.CalibrationItem(euriborCapletAnalytic, euriborVolatility, 1.0);

            double discountFactor = discountCurve.getDiscountFactor(analyticModel, maturity + periodLength);
            double eoniaValue = AnalyticFormulas.blackModelCapletValue(eoniaStrike, eoniaVolatility, maturity, eoniaStrike, periodLength, discountFactor);
            double euriborValue = AnalyticFormulas.blackModelCapletValue(euriborStrike, euriborVolatility, maturity, euriborStrike, periodLength, discountFactor);

            eoniaCalibrationItemsMonteCarlo[component - 1] = new LIBORMarketModelInterface.CalibrationItem(eoniaCaplet, eoniaValue, 1.0);
            euriborCalibrationItemsMonteCarlo[component - 1] = new LIBORMarketModelInterface.CalibrationItem(euriborCaplet, euriborValue, 1.0);
        }

        return new CalibrationData(liborPeriodDiscretization,
                liborPeriodDiscretization,
                analyticModel,
                "EONIA",
                "EURIBOR",
                eoniaCalibrationItemsAnalytic,
                euriborCalibrationItemsAnalytic,
                eoniaCalibrationItemsMonteCarlo,
                euriborCalibrationItemsMonteCarlo);
    }
}
