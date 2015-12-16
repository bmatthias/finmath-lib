/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 10.02.2004
 */
package net.finmath.tests.montecarlo.interestrate;

import net.finmath.exception.CalculationException;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.*;
import net.finmath.modelling.ProductInterface;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.*;
import net.finmath.montecarlo.interestrate.LIBORMarketModel.Measure;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface.CalibrationItem;
import net.finmath.montecarlo.interestrate.modelplugins.*;
import net.finmath.montecarlo.interestrate.products.*;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.tests.montecarlo.interestrate.tools.CSVCreator;
import net.finmath.tests.montecarlo.interestrate.tools.CalibrationData;
import net.finmath.tests.montecarlo.interestrate.tools.TexCreator;
import net.finmath.tests.montecarlo.interestrate.tools.calibrationdata.*;
import net.finmath.time.TimeDiscretizationInterface;
import org.apache.commons.math3.primes.Primes;

import java.io.File;
import java.io.FileWriter;
import java.util.*;

/**
 * This class tests the LIBOR market model and products.
 *
 * @author Christian Fries
 */
public class MultiCurveLIBORMarketModelCalibrationTest {
    private static final String VALUE = "VALUE";
    private static final String VOLATILITY = "VOLATILITY";

    private int[] years = { 2010, 2011, 2012, 2013, 2014 };
    private int numberOfFactors = 20;
    private int numberOfPaths = 10000;
    private Integer numberOfParams = 17;
    private String measure = Measure.SPOT.name();
    private String multiCurveModel = MultiCurveLIBORMarketModel.MultiCurveModel.ADDITIVE.name();
    private boolean useAnalyticApproximation = true;
    private boolean useSeperateCorrelationModels = false;
    private String tenor = "1m";
    private String productType;
    private int numberOfSeeds = 1;

    private boolean calibrateSC = true;
    private boolean calibrateMC = true;
    private String stateSpace = LIBORMarketModel.StateSpace.LOGNORMAL.name();

    public MultiCurveLIBORMarketModelCalibrationTest() {}

    private MultiCurveLIBORMarketModelCalibrationTest(
            Integer numberOfFactors, Integer numberOfPaths, Measure measure,
            MultiCurveLIBORMarketModel.MultiCurveModel multiCurveModel, Boolean useAnalyticApproximation,
            Boolean useSeperateCorrelationModels, String tenor, Integer numberOfParams, String productType) {
        this.numberOfPaths = useAnalyticApproximation ? numberOfPaths : numberOfPaths / 10;
        this.numberOfParams = numberOfParams;
        this.numberOfFactors = numberOfFactors;
        this.productType = productType;
        this.measure = measure.name();
        this.multiCurveModel = multiCurveModel.name();
        this.useAnalyticApproximation = useAnalyticApproximation;
        this.useSeperateCorrelationModels = useSeperateCorrelationModels;
        this.tenor = tenor;
    }

    public static void main(String[] args) throws Exception {
        Integer[] numberOfParams = { 17 };
        Integer[] numberOfFactors = { 10 };
        Integer[] numberOfPaths = { 5000 };
        LIBORMarketModel.Measure[] measures = { Measure.SPOT };
        MultiCurveLIBORMarketModel.MultiCurveModel[] models = {
                MultiCurveLIBORMarketModel.MultiCurveModel.ADDITIVE,
                //MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE,
                //MultiCurveLIBORMarketModel.MultiCurveModel.MMARTINGALE
        };
        Boolean[] useAnalyticApproximation = { true };
        Boolean[] useSeperateCorrelationModels = { false };
        String[] tenors = new String[] { null };
        String[] swaptionTypes = new String[] { "synthCaplets" };

        Object[][] parameters = new Object[][]{
                numberOfFactors, numberOfPaths, measures, models, useAnalyticApproximation,
                useSeperateCorrelationModels, tenors, numberOfParams, swaptionTypes
        };

        Object[] parameterSet = new Object[9];
        run(0, parameters, parameterSet);
    }

    private static void run(int column, Object[][] parameters, Object[] parameterSet) {
        if (column == parameters.length) {
            MultiCurveLIBORMarketModelCalibrationTest test = null;
            try {
                test = new MultiCurveLIBORMarketModelCalibrationTest(
                        (Integer) parameterSet[0],
                        (Integer) parameterSet[1],
                        (Measure) parameterSet[2],
                        (MultiCurveLIBORMarketModel.MultiCurveModel) parameterSet[3],
                        (Boolean) parameterSet[4],
                        (Boolean) parameterSet[5],
                        (String) parameterSet[6],
                        (Integer) parameterSet[7],
                        (String) parameterSet[8]);
                switch (test.productType) {
                    case "synthSwaptions":
                        test.start(SyntheticSwaptions.newInstance());
                        break;
                    case "synthCaplets":
                        test.start(SyntheticCaplets.newInstance());
                        break;
                    default:
                        throw new IllegalArgumentException("Type not supported: " + test.productType);
                }
            } catch (Exception e) {
                e.printStackTrace();
            } finally {
                test = null;
                System.gc();
            }
        } else {
            for (int j = 0; j < parameters[column].length; j++) {
                parameterSet[column] = parameters[column][j];
                run(column + 1, parameters, parameterSet);
            }
        }
    }

    private void start(CalibrationData calibrationData) throws Exception {
        TimeDiscretizationInterface timeDiscretization = calibrationData.getTimeDiscretization();
        TimeDiscretizationInterface liborPeriodDiscretization = calibrationData.getPeriodDiscretization();
        AnalyticModelInterface analyticModel = calibrationData.getAnalyticModel();
        String eoniaCurveName = calibrationData.getFirstCurveName();
        String forwardCurveName = calibrationData.getSecondCurveName();
        CalibrationItem[] eoniaCalibrationItemsAnalytic = calibrationData.getFirstCurveCalibrationItems();
        CalibrationItem[] euriborCalibrationItemsAnalytic = calibrationData.getSecondCurveCalibrationItems();
        CalibrationItem[] eoniaCalibrationItemsMonteCarlo = calibrationData.getFirstCurveCalibrationItemsMonteCarlo();
        CalibrationItem[] euriborCalibrationItemsMonteCarlo = calibrationData.getSecondCurveCalibrationItemsMonteCarlo();

        final String baseName = numberOfSeeds + productType + measure + numberOfFactors + multiCurveModel + numberOfPaths +
                useAnalyticApproximation + tenor + useSeperateCorrelationModels + numberOfParams;

        System.out.println(baseName);

        Map<String, Object> scProperties = new HashMap<>();
        scProperties.put("measure", measure);
        scProperties.put("stateSpace", stateSpace);

        Map<String, Object> mcProperties = new HashMap<>();
        mcProperties.put("measure", measure);
        mcProperties.put("stateSpace", stateSpace);
        mcProperties.put("multiCurveModel", multiCurveModel);

		/*
		 * Create a LIBOR Market Model
		 */
        AbstractLIBORCovarianceModelParametric covarianceModelParametric =
                getSingleCurveCovarianceModel(timeDiscretization, liborPeriodDiscretization);

        /*covarianceModelParametric.setParameter(new double[]{
                0.1, -0.01, 0.05, 0.5, 0.1, 0.5, 0.1
//                0.0, 0.0, 0.10, 0.2, 0.0
//                0.20, 0.05, 0.10, 0.20, 0.10
        });*/

        ForwardCurveInterface forwardCurve = analyticModel.getForwardCurve(eoniaCurveName);
        DiscountCurveInterface discountCurve = analyticModel.getDiscountCurve(forwardCurve.getDiscountCurveName());
        LIBORMarketModelInterface scLiborMarketModelCalibrated = new LIBORMarketModel(
                liborPeriodDiscretization,
                forwardCurve, discountCurve, covarianceModelParametric,
                calibrateSC ? eoniaCalibrationItemsAnalytic : new CalibrationItem[0],
                scProperties);

        double[] scParams = ((AbstractLIBORCovarianceModelParametric) scLiborMarketModelCalibrated.getCovarianceModel()).getParameter();
        for (double p : scParams) System.out.println(p);
        System.out.println();

        covarianceModelParametric =
                getMultiCurveCovarianceModel(timeDiscretization, liborPeriodDiscretization, scParams);

        /*covarianceModelParametric.setParameter(new double[]{
//                0.0, 0.00, -0.04, 0.2, 2.0, 0.3, 1.0, 0.15, 0.03, 0.0, 0.55
//                0.20, 0.05, 0.10, 0.20, 0.10, 0.0, Double.MAX_VALUE, Double.MAX_VALUE, 0.0 //perfectly correlated spread in 13 param model
                0.0, 0.0, 0.0, 0.0, Double.MAX_VALUE, 0.0, Double.MAX_VALUE, 0.0 //zero vola spread in 15 param model
        });*/
        /*double[] spreadParams = new double[10];
        System.arraycopy(scParams, 0, spreadParams, 0, 7);
        spreadParams[7] = 1.0; spreadParams[8] = Double.MAX_VALUE; spreadParams[9] = Double.MAX_VALUE;
        covarianceModelParametric.setParameter(spreadParams);*/

        long start = System.currentTimeMillis();

        LIBORMarketModelInterface mcLiborMarketModelCalibrated = new MultiCurveLIBORMarketModel(
                liborPeriodDiscretization, analyticModel,
                forwardCurveName, eoniaCurveName, covarianceModelParametric,
                calibrateMC ? euriborCalibrationItemsAnalytic : new CalibrationItem[0],
                mcProperties);

        long duration = System.currentTimeMillis() - start;

        double[] mcParams = ((AbstractLIBORCovarianceModelParametric) mcLiborMarketModelCalibrated.getCovarianceModel()).getParameter();
        for (double p : mcParams) System.out.println(p);
        System.out.println();

		/*
		 * Test our calibration
		 */
        boolean isSpotMeasure = Measure.valueOf(measure) == Measure.SPOT;
        Map<String, Object> dataModified = new HashMap<>();
        dataModified.put("measure", isSpotMeasure ? Measure.TERMINAL.name() : Measure.SPOT.name());

        LIBORMarketModelInterface scModelSpot = isSpotMeasure ? scLiborMarketModelCalibrated : scLiborMarketModelCalibrated.getCloneWithModifiedData(dataModified);
        LIBORMarketModelInterface mcModelSpot = isSpotMeasure ? mcLiborMarketModelCalibrated : mcLiborMarketModelCalibrated.getCloneWithModifiedData(dataModified);

        LIBORMarketModelInterface scModelTerminal = !isSpotMeasure ? scLiborMarketModelCalibrated : scLiborMarketModelCalibrated.getCloneWithModifiedData(dataModified);
        LIBORMarketModelInterface mcModelTerminal = !isSpotMeasure ? mcLiborMarketModelCalibrated : mcLiborMarketModelCalibrated.getCloneWithModifiedData(dataModified);

        File outputFile = new File("results/" + baseName + ".csv");
        String fileName = baseName;
        for (int i = 0; outputFile.exists(); i++) {
            fileName = baseName + "_" + i;
            outputFile = new File("results/" + fileName + ".csv");
        }

        outputFile.getParentFile().mkdirs();
        FileWriter writer = new FileWriter(outputFile);
        try {
            int primeNumber = 1009;
            List<Integer> primes = new ArrayList<>();
            while (primes.size() < numberOfSeeds) {
                primes.add(primeNumber);
                primeNumber = Primes.nextPrime(primeNumber + 1);
            }

            double[][] results = new double[16][euriborCalibrationItemsAnalytic.length];

            LIBORModelMonteCarloSimulation mcSimSpot = null, mcSimTerminal = null, scSimSpot = null, scSimTerminal = null;
            ProcessEulerScheme mcProcessSpot = null, mcProcessTerminal = null, scProcessSpot = null, scProcessTerminal = null;
            for (Integer prime : primes) {
                System.out.println(prime);
                BrownianMotionInterface scBrownianMotion = new BrownianMotion(timeDiscretization, numberOfFactors, numberOfPaths, prime /* seed */);
                BrownianMotionInterface mcBrownianMotion = new BrownianMotion(timeDiscretization, 2 * numberOfFactors, numberOfPaths, prime /* seed */);

                scProcessSpot = new ProcessEulerScheme(scBrownianMotion);
                scProcessTerminal = new ProcessEulerScheme(scBrownianMotion);
                scSimSpot = new LIBORModelMonteCarloSimulation(scModelSpot, scProcessSpot);
                //scSimTerminal = new LIBORModelMonteCarloSimulation(scModelTerminal, scProcessTerminal);

                mcProcessSpot = new ProcessEulerScheme(mcBrownianMotion);
                mcProcessTerminal = new ProcessEulerScheme(mcBrownianMotion);
                mcSimSpot  = new LIBORModelMonteCarloSimulation(mcModelSpot, mcProcessSpot);
                //mcSimTerminal = new LIBORModelMonteCarloSimulation(mcModelTerminal, mcProcessTerminal);

                CSVCreator.printProcess(mcProcessSpot, mcModelSpot);
                CSVCreator.printProcess(scProcessSpot, scModelSpot);

                for (int i = 0; i < euriborCalibrationItemsAnalytic.length; i++) {
                    results[0][i] += euriborCalibrationItemsAnalytic[i].calibrationTargetValue;
                    results[1][i] += eoniaCalibrationItemsAnalytic[i].calibrationTargetValue;

                    results[2][i] += getOptionValue(euriborCalibrationItemsAnalytic[i].calibrationProduct, mcSimSpot, VOLATILITY);
                    results[3][i] += getOptionValue(eoniaCalibrationItemsAnalytic[i].calibrationProduct, scSimSpot, VOLATILITY);

                    //results[0][i] += results[2][i];
                    //results[1][i] += results[3][i];

                    results[4][i] += getOptionValue(euriborCalibrationItemsAnalytic[i].calibrationProduct, mcSimSpot, VALUE);
                    results[5][i] += getOptionValue(eoniaCalibrationItemsAnalytic[i].calibrationProduct, scSimSpot, VALUE);

                    results[6][i] += euriborCalibrationItemsMonteCarlo[i].calibrationTargetValue;
                    results[7][i] += eoniaCalibrationItemsMonteCarlo[i].calibrationTargetValue;

                    results[8][i] += getOptionValue(euriborCalibrationItemsMonteCarlo[i].calibrationProduct, mcSimSpot, VALUE);
                    results[9][i] += getOptionValue(eoniaCalibrationItemsMonteCarlo[i].calibrationProduct, scSimSpot, VALUE);

                    results[10][i] += getOptionValue(euriborCalibrationItemsMonteCarlo[i].calibrationProduct, mcSimSpot, VOLATILITY);
                    results[11][i] += getOptionValue(eoniaCalibrationItemsMonteCarlo[i].calibrationProduct, scSimSpot, VOLATILITY);

                    results[12][i] += getOptionValue(euriborCalibrationItemsMonteCarlo[i].calibrationProduct, mcSimTerminal, VALUE);
                    results[13][i] += getOptionValue(eoniaCalibrationItemsMonteCarlo[i].calibrationProduct, scSimTerminal, VALUE);

                    results[14][i] += getOptionValue(euriborCalibrationItemsMonteCarlo[i].calibrationProduct, mcSimTerminal, VOLATILITY);
                    results[15][i] += getOptionValue(eoniaCalibrationItemsMonteCarlo[i].calibrationProduct, scSimTerminal, VOLATILITY);
                }
            }

            double[] deviations = new double[36];

            for (int calibrationItemIndex = 0; calibrationItemIndex < results[0].length; calibrationItemIndex++) {
                for (int resultIndex = 0; resultIndex < results.length; resultIndex++) {
                    results[resultIndex][calibrationItemIndex] /= primes.size();
                }

                CSVCreator.writeValues(writer, results, deviations, calibrationItemIndex);
            }

            writer.flush();
            writer.close();

            writer = new FileWriter(new File("results/" + fileName + "_errs.csv"));
            CSVCreator.writeErrorStatistics(eoniaCalibrationItemsAnalytic, euriborCalibrationItemsAnalytic, writer, deviations);
            writer.flush();
            writer.close();

            writer = new FileWriter(new File("results/" + fileName + "_params.csv"));
            CSVCreator.writeParamsAndDuration(scParams, duration, mcParams, writer);
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            writer.flush();
            writer.close();
        }

        TexCreator.createPlot(fileName,
                Math.max((int) Math.floor(euriborCalibrationItemsAnalytic[euriborCalibrationItemsAnalytic.length - 1].calibrationTargetValue * 10 - 1) * 10, 0),
                (int) Math.floor(euriborCalibrationItemsAnalytic[0].calibrationTargetValue * 10 + 2) * 10,
                numberOfPaths, numberOfParams, numberOfFactors,
                multiCurveModel, useSeperateCorrelationModels
        );

        System.out.println("Calibration time:\t" + duration);
        System.out.println("__________________________________________________________________________________________\n");
    }

    private AbstractLIBORCovarianceModelParametric getMultiCurveCovarianceModel(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, double[] scParams) {
        double[] fixedVolParams = new double[4];
        System.arraycopy(scParams, 0, fixedVolParams, 0, 4);

        double[] fixedCorParams = new double[scParams.length - 4];
        System.arraycopy(scParams, 4, fixedCorParams, 0, scParams.length - 4);

        LIBORCorrelationModel correlationModel;
        if (MultiCurveLIBORMarketModel.MultiCurveModel.MMARTINGALE.name().equals(multiCurveModel.toUpperCase())) {
            correlationModel = new LIBORCorrelationModelIndependentCurveSevenParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true);
        } else if (numberOfParams == 13) {
            correlationModel = useSeperateCorrelationModels ?
                    new LIBORCorrelationModelSeperateCurveFiveParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true) :
                    new LIBORCorrelationModelTwoCurveFiveParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true);
        } else if (numberOfParams == 15) {
            correlationModel = new LIBORCorrelationModelSeperateCurve7ParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true);
        } else if (numberOfParams == 16) {
            correlationModel = new LIBORCorrelationModelSeperateCurveSevenParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true);
        } else if (numberOfParams == 17) {
            correlationModel = useSeperateCorrelationModels ?
                    new LIBORCorrelationModelSeperateCurveNineParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true) :
                    new LIBORCorrelationModelTwoCurveNineParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true);
        } else if (numberOfParams == 18) {
            correlationModel = useSeperateCorrelationModels ?
                    new LIBORCorrelationModelSeperateCurveTenParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true) :
                    new LIBORCorrelationModelTwoCurveTenParameterExponentialDecay(timeDiscretization, liborPeriodDiscretization, 2 * numberOfFactors, fixedCorParams, true);
        } else {
            throw new IllegalArgumentException("" + numberOfParams + " params not supported.");
        }

        LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTwoCurveEightParameterExponentialForm(timeDiscretization, liborPeriodDiscretization, fixedVolParams, true);
        return new LIBORCovarianceModelFromVolatilityAndCorrelation(timeDiscretization, liborPeriodDiscretization, volatilityModel, correlationModel);
    }

    private AbstractLIBORCovarianceModelParametric getSingleCurveCovarianceModel(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization) {
        AbstractLIBORCovarianceModelParametric covarianceModelParametric;
        if (numberOfParams == 13) {
            covarianceModelParametric = new LIBORCovarianceModelExponentialForm5Param(timeDiscretization, liborPeriodDiscretization, numberOfFactors);
        } else {
            covarianceModelParametric = new LIBORCovarianceModelExponentialForm7Param(timeDiscretization, liborPeriodDiscretization, numberOfFactors);
        }
        return covarianceModelParametric;
    }

    private double getOptionValue(AbstractLIBORMonteCarloProduct product,
                                  LIBORModelMonteCarloSimulationInterface simulation,
                                  String valueUnit) throws CalculationException {
        if (simulation == null) {
            return 0.0;
        } else if (product instanceof SwaptionAnalyticApproximation) {
            return ((SwaptionAnalyticApproximation) product).getValues(0.0, simulation.getModel(), ProductInterface.ValueUnit.valueOf(valueUnit)).getAverage();
        } else if (product instanceof SwaptionSimple) {
            return ((SwaptionSimple) product).getValue(0.0, simulation, ProductInterface.ValueUnit.valueOf(valueUnit)).getAverage();
        } else if (product instanceof Caplet) {
            return ((Caplet) product).getValue(0.0, simulation, ProductInterface.ValueUnit.valueOf(valueUnit)).getAverage();
        } else if (product instanceof CapletAnalyticApproximation) {
            return ((CapletAnalyticApproximation) product).getValue(0.0, simulation,
                    CapletAnalyticApproximation.MultiCurveApproximation.FENTON_WILKINSON,
                    ProductInterface.ValueUnit.valueOf(valueUnit)).getAverage();
        } else {
            throw new IllegalArgumentException("Product must be a Caplet or Swaption!");
        }
    }
}