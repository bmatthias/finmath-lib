/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christianfries.com.
 *
 * Created on 31.03.2014
 */

package net.finmath.montecarlo.interestrate.products;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.*;
import net.finmath.montecarlo.interestrate.*;
import net.finmath.montecarlo.interestrate.modelplugins.*;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.time.TimeDiscretization;
import org.junit.Assert;
import org.junit.Test;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.Locale;
import java.util.Map;

/**
 * @author Christian Fries
 *
 */
public class SwaptionMultiCurveAnalyticApproximationTest {

	private static DecimalFormat formatterMaturity	= new DecimalFormat("00.00", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterValue		= new DecimalFormat(" ##0.000%;-##0.000%", new DecimalFormatSymbols(Locale.ENGLISH));
	private static DecimalFormat formatterDeviation	= new DecimalFormat(" 0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));

	private final int numberOfPaths = 20000;
	private final int numberOfFactors = 5;

	@Test
	public void testMultiCurveModel() throws CalculationException {
		System.out.println("Runnning tests with a multi curve LIBOR Market Model");

		LIBORModelMonteCarloSimulationInterface liborMarketModel = createMultiCurveLIBORMarketModel(numberOfPaths, numberOfFactors);

		testModel(liborMarketModel);
	}

	public void testModel(LIBORModelMonteCarloSimulationInterface liborMarketModel) throws CalculationException {
		/*
		 * Value a swaption
		 */
		System.out.println("Swaption prices:\n");
		System.out.println("Maturity      Simulation (MC)   Analytic (MC)    Deviation (AMC-S)");

		double maxAbsDeviationSimulation	= 0.0;
		for (int maturityIndex = 1; maturityIndex <= liborMarketModel.getNumberOfLibors() - 10; maturityIndex++) {

			double exerciseDate = liborMarketModel.getLiborPeriod(maturityIndex);
			System.out.print(formatterMaturity.format(exerciseDate) + "          ");

			int numberOfPeriods = 5;

			// Create a swaption

			double[] fixingDates = new double[numberOfPeriods];
			double[] paymentDates = new double[numberOfPeriods];
			double[] swapTenor = new double[numberOfPeriods + 1];
			double swapPeriodLength = 0.5;

			for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
				fixingDates[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
				paymentDates[periodStartIndex] = exerciseDate + (periodStartIndex + 1) * swapPeriodLength;
				swapTenor[periodStartIndex] = exerciseDate + periodStartIndex * swapPeriodLength;
			}
			swapTenor[numberOfPeriods] = exerciseDate + numberOfPeriods * swapPeriodLength;

			// Swaptions swap rate
			double swaprate = 0.05;// getParSwaprate(liborMarketModel, swapTenor);

			// Set swap rates for each period
			double[] swaprates = new double[numberOfPeriods];
			for (int periodStartIndex = 0; periodStartIndex < numberOfPeriods; periodStartIndex++) {
				swaprates[periodStartIndex] = swaprate;
			}

			// Value with Monte Carlo
			Swaption swaptionMonteCarlo	= new Swaption(exerciseDate, fixingDates, paymentDates, swaprates);
			double valueSimulation = swaptionMonteCarlo.getValue(liborMarketModel);
			System.out.print(formatterValue.format(valueSimulation) + "           ");
			
			// Value analytic
			SwaptionAnalyticApproximation swaptionAnalyitc = new SwaptionAnalyticApproximation(swaprate, swapTenor, SwaptionAnalyticApproximation.ValueUnit.VALUE);
			double valueAnalytic = swaptionAnalyitc.getValue(liborMarketModel);
			System.out.print(formatterValue.format(valueAnalytic) + "          ");

			// Absolute deviation
			double deviation1 = (valueAnalytic - valueSimulation);
			System.out.print(formatterDeviation.format(deviation1) + "\n");

			maxAbsDeviationSimulation	= Math.max(maxAbsDeviationSimulation, Math.abs(deviation1));
		}

		System.out.println("__________________________________________________________________________________________\n");

		/*
		 * jUnit assertion: condition under which we consider this test successful
		 */
		Assert.assertEquals("Deviation", 0.0, maxAbsDeviationSimulation, 5E-3);
	}

	public static LIBORModelMonteCarloSimulationInterface createMultiCurveLIBORMarketModel(int numberOfPaths, int numberOfFactors) throws CalculationException {

		// Create the forward curve (initial value of the LIBOR market model)
		ForwardCurve forwardCurve = ForwardCurve.createForwardCurveFromForwards(
				"forwardCurve"								/* name of the curve */,
				new double[] {0.5 , 1.0 , 2.0 , 5.0 , 40.0}	/* fixings of the forward */,
				new double[] {0.05, 0.05, 0.05, 0.05, 0.05}	/* forwards */,
				0.5											/* tenor / period length */
				);

		// Create the discount curve
		DiscountCurveInterface discountCurve = new DiscountCurveFromForwardCurve(ForwardCurve.createForwardCurveFromForwards(
                "discountCurve"								/* name of the curve */,
                new double[]{0.5, 1.0, 2.0, 5.0, 40.0}	/* maturities */,
                new double[]{0.04, 0.04, 0.04, 0.04, 0.05}	/* zero rates */,
                0.5
        ));

		return createLIBORMarketModel(numberOfPaths, numberOfFactors, discountCurve, forwardCurve);
	}

	public static LIBORModelMonteCarloSimulationInterface createLIBORMarketModel(
			int numberOfPaths, int numberOfFactors, DiscountCurveInterface discountCurve, ForwardCurveInterface forwardCurve) throws CalculationException {

		/*
		 * Create the libor tenor structure and the initial values
		 */
		double liborPeriodLength	= 0.5;
		double liborRateTimeHorzion	= 20.0;
		TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, (int) (liborRateTimeHorzion / liborPeriodLength), liborPeriodLength);		
		
		/*
		 * Create a simulation time discretization
		 */
		double lastTime	= 20.0;
		double dt		= 0.5;

		TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, (int) (lastTime / dt), dt);

        LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTwoCurveEightParameterExponentialForm(
                timeDiscretization, liborPeriodDiscretization,
                new double[]{ 0.34, -0.024, 0.088, 0.21},
                new double[]{ 0.1, 0.1, 0.1, 0.2 },
                true);

        LIBORCorrelationModel correlationModel = new LIBORCorrelationModelTwoCurveTenParameterExponentialDecay(
                timeDiscretization, liborPeriodDiscretization, numberOfFactors,
                new double[]{ 0.0, 2.16, 1.89 },
                new double[]{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 },
                true);
		/*
		 * Combine volatility model and correlation model to a covariance model
		 */
		LIBORCovarianceModelFromVolatilityAndCorrelation covarianceModel =
				new LIBORCovarianceModelFromVolatilityAndCorrelation(
                        timeDiscretization,	liborPeriodDiscretization, volatilityModel, correlationModel);

		// BlendedLocalVolatlityModel (future extension)
		//		AbstractLIBORCovarianceModel covarianceModel2 = new BlendedLocalVolatlityModel(covarianceModel, 0.00, false);

		// Set model properties
		Map<String, String> properties = new HashMap<String, String>();

		// Choose the simulation measure
		properties.put("measure", LIBORMarketModel.Measure.SPOT.name());
		
		// Choose log normal model
		properties.put("stateSpace", LIBORMarketModel.StateSpace.LOGNORMAL.name());

		// Empty array of calibration items - hence, model will use given covariance
		LIBORMarketModel.CalibrationItem[] calibrationItems = new LIBORMarketModel.CalibrationItem[0];

		/*
		 * Create corresponding LIBOR Market Model
		 */
		LIBORMarketModelInterface liborMarketModel = new MultiCurveLIBORMarketModel(
				liborPeriodDiscretization, new AnalyticModel(new CurveInterface[] { forwardCurve, discountCurve }),
                forwardCurve.getName(), discountCurve.getName(), covarianceModel, calibrationItems, properties);

		ProcessEulerScheme process = new ProcessEulerScheme(
				new net.finmath.montecarlo.BrownianMotion(timeDiscretization,
						numberOfFactors, numberOfPaths, 3141 /* seed */));
		//		process.setScheme(ProcessEulerScheme.Scheme.PREDICTOR_CORRECTOR);

		return new LIBORModelMonteCarloSimulation(liborMarketModel, process);
	}
}
