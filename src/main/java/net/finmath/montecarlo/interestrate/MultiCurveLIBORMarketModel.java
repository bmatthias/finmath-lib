package net.finmath.montecarlo.interestrate;

import net.finmath.exception.CalculationException;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.marketdata.model.curves.*;
import net.finmath.marketdata.model.volatilities.AbstractSwaptionMarketData;
import net.finmath.marketdata.products.Swap;
import net.finmath.marketdata.products.SwapAnnuity;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModel;
import net.finmath.montecarlo.interestrate.modelplugins.AbstractLIBORCovarianceModelParametric;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.interestrate.products.SwaptionAnalyticApproximation;
import net.finmath.montecarlo.interestrate.products.SwaptionSimple;
import net.finmath.montecarlo.model.AbstractModel;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.RegularSchedule;
import net.finmath.time.ScheduleInterface;
import net.finmath.time.TimeDiscretization;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class MultiCurveLIBORMarketModel extends AbstractModel implements LIBORMarketModelInterface {

    public enum Measure				{ SPOT, TERMINAL }
    public enum StateSpace			{ NORMAL, LOGNORMAL }
    public enum MultiCurveModel     { ADDITIVE, MULTIPLICATIVE, MMARTINGALE}

    private final TimeDiscretizationInterface liborPeriodDiscretization;
    private ForwardCurveInterface forwardRateCurve;
    private ForwardCurveInterface riskFreeCurve;
    private DiscountCurveInterface discountCurve;
    private AbstractLIBORCovarianceModel covarianceModel;
    private AnalyticModelInterface curveModel;

    private Measure	measure	= Measure.SPOT;
    private StateSpace stateSpace = StateSpace.LOGNORMAL;
    private MultiCurveModel multiCurveModel = MultiCurveModel.MULTIPLICATIVE;

    private double[][][][] integratedLIBORCovariance;
    private Double[][] c;
    private double liborCap = Double.POSITIVE_INFINITY;

    public MultiCurveLIBORMarketModel(TimeDiscretizationInterface liborPeriodDiscretization,
                                      AnalyticModelInterface curveModel,
                                      String liborCurveName,
                                      String riskFreeCurveName,
                                      AbstractLIBORCovarianceModel covarianceModel,
                                      CalibrationItem[] calibrationItems,
                                      Map<String, ?> properties) throws CalculationException {
        this.liborPeriodDiscretization = liborPeriodDiscretization;
        this.forwardRateCurve = curveModel.getForwardCurve(liborCurveName);
        this.riskFreeCurve = curveModel.getForwardCurve(riskFreeCurveName);
        this.discountCurve = curveModel.getDiscountCurve(riskFreeCurve.getDiscountCurveName());
        this.covarianceModel = covarianceModel;
        this.curveModel = curveModel;

        if(properties != null && properties.containsKey("measure"))	    measure		= Measure.valueOf((properties.get("measure").toString()).toUpperCase());
        if(properties != null && properties.containsKey("stateSpace"))	stateSpace	= StateSpace.valueOf((properties.get("stateSpace").toString()).toUpperCase());
        if(properties != null && properties.containsKey("multiCurveModel"))	multiCurveModel = MultiCurveModel.valueOf((properties.get("multiCurveModel").toString()).toUpperCase());

        Map<String,Object> calibrationParameters = null;
        if(properties != null && properties.containsKey("calibrationParameters"))	calibrationParameters	= (Map<String,Object>)properties.get("calibrationParameters");

        // Perform calibration, if data is given
        if(calibrationItems != null && calibrationItems.length > 0) {
            AbstractLIBORCovarianceModelParametric covarianceModelParametric = null;
            try {
                covarianceModelParametric = (AbstractLIBORCovarianceModelParametric)covarianceModel;
            }
            catch(Exception e) {
                throw new ClassCastException("Calibration is currently restricted to parametric covariance models (AbstractLIBORCovarianceModelParametric).");
            }

            // @TODO Should be more elegant. Convert array for constructor
            AbstractLIBORMonteCarloProduct[]	calibrationProducts		= new AbstractLIBORMonteCarloProduct[calibrationItems.length];
            double[]							calibrationTargetValues	= new double[calibrationItems.length];
            double[]							calibrationWeights		= new double[calibrationItems.length];
            for(int i=0; i<calibrationTargetValues.length; i++) {
                calibrationProducts[i]		= calibrationItems[i].calibrationProduct;
                calibrationTargetValues[i]	= calibrationItems[i].calibrationTargetValue;
                calibrationWeights[i]		= calibrationItems[i].calibrationWeight;
            }

            this.covarianceModel    = covarianceModelParametric.getCloneCalibrated(this, calibrationProducts, calibrationTargetValues, calibrationWeights, calibrationParameters);
        }
    }

    public MultiCurveLIBORMarketModel(TimeDiscretizationInterface liborPeriodDiscretization,
                                      AnalyticModelInterface curveModel,
                                      String liborCurveName,
                                      String riskFreeCurveName,
                                      AbstractLIBORCovarianceModel covarianceModel,
                                      AbstractSwaptionMarketData swaptionMarketData,
                                      Map<String, ?> properties) throws CalculationException {
        this(liborPeriodDiscretization,
                curveModel,
                liborCurveName,
                riskFreeCurveName,
                covarianceModel,
                LIBORMarketModelInterface.getCalibrationItems(
                        liborPeriodDiscretization,
                        curveModel.getForwardCurve(liborCurveName),
                        swaptionMarketData,
                        true
                ),
                properties);
    }

    public MultiCurveModel getMultiCurveModel() {
        return multiCurveModel;
    }

    public ForwardCurveInterface getRiskFreeCurve() {
        return riskFreeCurve;
    }

    @Override
    public RandomVariableInterface getLIBOR(int timeIndex, int liborIndex) throws CalculationException {
        return getProcessValue(timeIndex, liborIndex);
    }

    public RandomVariableInterface getForward(int timeIndex, int liborIndex) throws CalculationException {
        return getProcess().getProcessValue(timeIndex, liborIndex);
    }

    @Override
    public RandomVariableInterface getProcessValue(int timeIndex, int componentIndex) throws CalculationException {
        RandomVariableInterface forward = getProcess().getProcessValue(timeIndex, componentIndex);
        RandomVariableInterface spread = getProcess().getProcessValue(timeIndex, componentIndex + getNumberOfComponents() / 2);

        switch (multiCurveModel) {
            case ADDITIVE:
                return forward.add(spread);
            case MULTIPLICATIVE:
                return forward.add(spread).add(forward.mult(spread).mult(getTimeDiscretization().getTimeStep(timeIndex)));
            case MMARTINGALE:
                spread = spread.mult(forward.pow(getCorrelationFactor())).mult(getNormalizationFactor(componentIndex, timeIndex));
                return forward.add(spread).add(forward.mult(spread).mult(getTimeDiscretization().getTimeStep(timeIndex)));
            default:
                throw new IllegalArgumentException("Model type " + multiCurveModel + " not supported.");
        }
    }

    @Override
    public TimeDiscretizationInterface getLiborPeriodDiscretization() {
        return liborPeriodDiscretization;
    }

    @Override
    public int getNumberOfLibors() {
        return getNumberOfComponents() / 2;
    }

    @Override
    public double getLiborPeriod(int timeIndex) {
        if(timeIndex >= liborPeriodDiscretization.getNumberOfTimes() || timeIndex < 0) {
            throw new ArrayIndexOutOfBoundsException("Index for LIBOR period discretization out of bounds: " + timeIndex + ".");
        } else {
            return liborPeriodDiscretization.getTime(timeIndex);
        }
    }

    @Override
    public int getLiborPeriodIndex(double time) {
        return liborPeriodDiscretization.getTimeIndex(time);
    }

    @Override
    public AnalyticModelInterface getAnalyticModel() {
        return curveModel;
    }

    @Override
    public DiscountCurveInterface getDiscountCurve() {
        if(discountCurve == null) {
            discountCurve = new DiscountCurveFromForwardCurve(getForwardRateCurve());
        }

        return discountCurve;
    }

    @Override
    public ForwardCurveInterface getForwardRateCurve() {
        return forwardRateCurve;
    }

    @Override
    public AbstractLIBORCovarianceModel getCovarianceModel() {
        return covarianceModel;
    }

    @Override
    public LIBORMarketModelInterface getCloneWithModifiedCovarianceModel(AbstractLIBORCovarianceModel calibrationCovarianceModel) {
        MultiCurveLIBORMarketModel model = (MultiCurveLIBORMarketModel)this.clone();
        model.covarianceModel = calibrationCovarianceModel;
        return model;
    }

    @Override
    public LIBORMarketModelInterface getCloneWithModifiedData(Map<String, Object> dataModified) throws CalculationException {
        TimeDiscretizationInterface liborPeriodDiscretization	= this.liborPeriodDiscretization;
        AnalyticModelInterface analyticModel				= this.curveModel;
        ForwardCurveInterface forwardRateCurve			= this.forwardRateCurve;
        ForwardCurveInterface riskFreeCurve				= this.riskFreeCurve;
        AbstractLIBORCovarianceModel covarianceModel				= this.covarianceModel;
        AbstractSwaptionMarketData swaptionMarketData			    = null;		// No recalibration, unless new swaption data is specified
        Map<String, Object>				properties					= new HashMap<String, Object>();
        properties.put("measure",		measure.name());
        properties.put("stateSpace",	stateSpace.name());
        properties.put("multiCurveModel",	multiCurveModel.name());

        if(dataModified.containsKey("liborPeriodDiscretization")) {
            liborPeriodDiscretization = (TimeDiscretizationInterface)dataModified.get("liborPeriodDiscretization");
        }
        if(dataModified.containsKey("liborCurve")) {
            forwardRateCurve = (ForwardCurveInterface)dataModified.get("liborCurve");
        }
        if(dataModified.containsKey("riskFreeCurve")) {
            riskFreeCurve = (ForwardCurveInterface)dataModified.get("riskFreeCurve");
        }
        if(dataModified.containsKey("forwardRateShift")) {
            throw new RuntimeException("Forward rate shift clone currently disabled.");
        }
        if(dataModified.containsKey("covarianceModel")) {
            covarianceModel = (AbstractLIBORCovarianceModel)dataModified.get("covarianceModel");
        }
        if(dataModified.containsKey("swaptionMarketData")) {
            swaptionMarketData = (AbstractSwaptionMarketData)dataModified.get("swaptionMarketData");
        }
        if(dataModified.containsKey("measure")) {
            properties.put("measure", dataModified.get("measure").toString().toUpperCase());
        }

        MultiCurveLIBORMarketModel newModel = new MultiCurveLIBORMarketModel(
                liborPeriodDiscretization,
                analyticModel,
                forwardRateCurve.getName(),
                riskFreeCurve.getName(),
                covarianceModel,
                swaptionMarketData,
                properties
        );
        return newModel;
    }

    @Override
    public Object clone() {
        Map<String, Object> properties = new HashMap<String, Object>();
        properties.put("measure",		measure.name());
        properties.put("stateSpace",	stateSpace.name());
        properties.put("multiCurveModel",	multiCurveModel.name());
        try {
            return new MultiCurveLIBORMarketModel(
                    liborPeriodDiscretization,
                    curveModel,
                    forwardRateCurve.getName(),
                    riskFreeCurve.getName(),
                    covarianceModel,
                    new CalibrationItem[0],
                    properties);
        } catch (CalculationException e) {
            return null;
        }
    }

    @Override
    public synchronized double[][][] getIntegratedLIBORCovariance() {
        return getIntegratedLIBORCovariances()[0];
    }

    @Override
    public int getNumberOfComponents() {
        return liborPeriodDiscretization.getNumberOfTimeSteps() * 2;
    }

    @Override
    public RandomVariableInterface applyStateSpaceTransform(int componentIndex, RandomVariableInterface randomVariable) {
        RandomVariableInterface value = randomVariable;

        if(stateSpace == StateSpace.LOGNORMAL)	value = value.exp();

        if(!Double.isInfinite(liborCap)) value = value.cap(liborCap);

        return value;
    }

    @Override
    public RandomVariableInterface[] getInitialState() {

        RandomVariableInterface[] initialStateRandomVariable = new RandomVariableInterface[getNumberOfComponents()];
        for(int componentIndex = 0; componentIndex < getNumberOfComponents() / 2; componentIndex++) {
            double timeStep = liborPeriodDiscretization.getTimeStep(componentIndex);
            double libor = forwardRateCurve.getForward(curveModel, liborPeriodDiscretization.getTime(componentIndex));
            double forward = riskFreeCurve.getForward(curveModel, liborPeriodDiscretization.getTime(componentIndex));

            double spread = libor - forward;

            if (multiCurveModel == MultiCurveModel.MULTIPLICATIVE) {
                spread /= (1.0 + timeStep * forward);
            } else if (multiCurveModel == MultiCurveModel.MMARTINGALE) {
                spread /= (1.0 + timeStep * forward) * Math.pow(forward, getCorrelationFactor());
            }

            double initialForwardState = (stateSpace == StateSpace.LOGNORMAL) ? Math.log(Math.max(forward, 0)) : forward; //We expect positive interest rates
            double initialSpreadState = (stateSpace == StateSpace.LOGNORMAL) ? Math.log(Math.max(spread, 0)) : spread; //We expect positive spreads

            initialStateRandomVariable[componentIndex] = new RandomVariable(initialForwardState);
            initialStateRandomVariable[componentIndex + getNumberOfComponents() / 2] = new RandomVariable(initialSpreadState);
        }
        return initialStateRandomVariable;
    }

    @Override
    public RandomVariableInterface getNumeraire(double time) throws CalculationException {
        int timeIndex = getLiborPeriodIndex(time);

        if(timeIndex < 0) {
            // Interpolation of Numeraire: log linear interpolation.
            int upperIndex = -timeIndex-1;
            int lowerIndex = upperIndex-1;
            if(lowerIndex < 0) throw new IllegalArgumentException("Numeraire requested for time " + time + ". Unsupported");

            double alpha = (time-getLiborPeriod(lowerIndex)) / (getLiborPeriod(upperIndex) - getLiborPeriod(lowerIndex));
            RandomVariableInterface numeraire = getNumeraire(getLiborPeriod(upperIndex)).log().mult(alpha).add(getNumeraire(getLiborPeriod(lowerIndex)).log().mult(1.0-alpha)).exp();

			/*
			 * Adjust for discounting, i.e. funding or collateralization
			 */
            if(discountCurve != null) {
                // This includes a control for zero bonds
                double deterministicNumeraireAdjustment = numeraire.invert().getAverage() / discountCurve.getDiscountFactor(curveModel, time);
                numeraire = numeraire.mult(deterministicNumeraireAdjustment);
            }

            return numeraire;
        }

        // Calculate the numeraire, when time is part of liborPeriodDiscretization

        // Get the start and end of the product
        int firstLiborIndex, lastLiborIndex;

        if(measure == Measure.TERMINAL) {
            firstLiborIndex	= getLiborPeriodIndex(time);
            if(firstLiborIndex < 0) {
                throw new CalculationException("Simulation time discretization not part of forward rate tenor discretization.");
            }

            lastLiborIndex 	= liborPeriodDiscretization.getNumberOfTimeSteps()-1;
        } else if(measure == Measure.SPOT) {
            // Spot measure
            firstLiborIndex	= 0;
            lastLiborIndex	= getLiborPeriodIndex(time)-1;
        } else {
            throw new CalculationException("Numeraire not implemented for specified measure.");
        }

		/*
		 * Calculation of the numeraire
		 */

        // Initialize to 1.0
        RandomVariableInterface numeraire = getProcess().getBrownianMotion().getRandomVariableForConstant(1.0);

        // The product
        for(int liborIndex = firstLiborIndex; liborIndex <= lastLiborIndex; liborIndex++) {
            RandomVariableInterface forward = getForward(getTimeIndex(Math.min(time, liborPeriodDiscretization.getTime(liborIndex))), liborIndex);

            double periodLength = liborPeriodDiscretization.getTimeStep(liborIndex);

            if(measure == Measure.SPOT) {
                numeraire = numeraire.accrue(forward, periodLength);
            } else {
                numeraire = numeraire.discount(forward, periodLength);
            }
        }

		/*
		 * Adjust for discounting
		 */
        if(discountCurve != null) {
            double deterministicNumeraireAdjustment = numeraire.invert().getAverage() / discountCurve.getDiscountFactor(curveModel, time);
            numeraire = numeraire.mult(deterministicNumeraireAdjustment);
        }
        return numeraire;
    }

    @Override
    public RandomVariableInterface[] getDrift(int timeIndex, RandomVariableInterface[] realizationAtTimeIndex, RandomVariableInterface[] realizationPredictor) {
        double	time				= getTime(timeIndex);
        int		firstLiborIndex		= this.getLiborPeriodIndex(time) + 1;
        if(firstLiborIndex<0) firstLiborIndex = -firstLiborIndex;

        RandomVariableInterface zero	= new RandomVariable(time, 0.0);

        // Allocate drift vector and initialize to zero (will be used to sum up drift components)
        RandomVariableInterface[]	drift = new RandomVariableInterface[getNumberOfComponents()];
        for(int componentIndex=firstLiborIndex; componentIndex < getNumberOfComponents() / 2; componentIndex++) {
            drift[componentIndex] = zero;
            drift[componentIndex + getNumberOfComponents() / 2] = zero;
        }

        RandomVariableInterface[]	covarianceFactorSums	= new RandomVariableInterface[getNumberOfFactors()];
        for(int factorIndex=0; factorIndex<getNumberOfFactors(); factorIndex++) {
            covarianceFactorSums[factorIndex] = zero;
        }

        int k = multiCurveModel == MultiCurveModel.MULTIPLICATIVE ? 1 : 0;
        if(measure == Measure.SPOT) {
            // Calculate drift for the component componentIndex (starting at firstLiborIndex, others are zero)
            for(int componentIndex = firstLiborIndex; componentIndex < getNumberOfComponents() / 2; componentIndex++) {
                double periodLength	= liborPeriodDiscretization.getTimeStep(componentIndex);
                RandomVariableInterface forward			= realizationAtTimeIndex[componentIndex];
                RandomVariableInterface oneStepMeasureTransform = (getProcess().getBrownianMotion().getRandomVariableForConstant(periodLength)).discount(forward, periodLength);

                if(stateSpace == StateSpace.LOGNORMAL) oneStepMeasureTransform = oneStepMeasureTransform.mult(forward);

                int spreadComponentIndex = componentIndex + getNumberOfComponents() / 2;
                RandomVariableInterface[]	forwardFactorLoading   	= covarianceModel.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
                RandomVariableInterface[]	spreadFactorLoading   	= covarianceModel.getFactorLoading(timeIndex, spreadComponentIndex - k, realizationAtTimeIndex);
                for(int factorIndex=0; factorIndex < getNumberOfFactors(); factorIndex++) {
                    if (multiCurveModel == MultiCurveModel.ADDITIVE) {
                        covarianceFactorSums[factorIndex] = covarianceFactorSums[factorIndex].addProduct(oneStepMeasureTransform, forwardFactorLoading[factorIndex]);
                        drift[spreadComponentIndex] = drift[spreadComponentIndex].addProduct(covarianceFactorSums[factorIndex], spreadFactorLoading[factorIndex]);
                    } else {
                        drift[spreadComponentIndex] = drift[spreadComponentIndex].addProduct(covarianceFactorSums[factorIndex], spreadFactorLoading[factorIndex]);
                        covarianceFactorSums[factorIndex] = covarianceFactorSums[factorIndex].addProduct(oneStepMeasureTransform, forwardFactorLoading[factorIndex]);
                    }
                    drift[componentIndex] = drift[componentIndex].addProduct(covarianceFactorSums[factorIndex], forwardFactorLoading[factorIndex]);
                }
            }
        }
        else if(measure == Measure.TERMINAL) {
            // Calculate drift for the component componentIndex (starting at firstLiborIndex, others are zero)
            for(int componentIndex = getNumberOfComponents() / 2 - 1; componentIndex >= firstLiborIndex; componentIndex--) {
                double					periodLength	= liborPeriodDiscretization.getTimeStep(componentIndex);
                RandomVariableInterface forward			= realizationAtTimeIndex[componentIndex];
                RandomVariableInterface oneStepMeasureTransform = (getProcess().getBrownianMotion().getRandomVariableForConstant(periodLength)).discount(forward, periodLength);

                if(stateSpace == StateSpace.LOGNORMAL) oneStepMeasureTransform = oneStepMeasureTransform.mult(forward);

                int spreadComponentIndex = componentIndex + getNumberOfComponents() / 2;
                RandomVariableInterface[]	forwardFactorLoading   	= covarianceModel.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
                RandomVariableInterface[]	spreadFactorLoading   	= covarianceModel.getFactorLoading(timeIndex, spreadComponentIndex, realizationAtTimeIndex);
                for(int factorIndex=0; factorIndex < getNumberOfFactors(); factorIndex++) {
                    drift[componentIndex] = drift[componentIndex].addProduct(covarianceFactorSums[factorIndex], forwardFactorLoading[factorIndex]);
                    if (multiCurveModel == MultiCurveModel.ADDITIVE) {
                        drift[spreadComponentIndex] = drift[spreadComponentIndex].addProduct(covarianceFactorSums[factorIndex], spreadFactorLoading[factorIndex]);
                        covarianceFactorSums[factorIndex] = covarianceFactorSums[factorIndex].sub(oneStepMeasureTransform.mult(forwardFactorLoading[factorIndex]));
                    } else if (multiCurveModel == MultiCurveModel.MULTIPLICATIVE) {
                        covarianceFactorSums[factorIndex] = covarianceFactorSums[factorIndex].sub(oneStepMeasureTransform.mult(forwardFactorLoading[factorIndex]));
                        drift[spreadComponentIndex] = drift[spreadComponentIndex].addProduct(covarianceFactorSums[factorIndex], spreadFactorLoading[factorIndex]);
                    } else {
                        covarianceFactorSums[factorIndex] = covarianceFactorSums[factorIndex].sub(oneStepMeasureTransform.mult(forwardFactorLoading[factorIndex]));
                    }
                }
            }
        }

        if(stateSpace == StateSpace.LOGNORMAL) {
            // Drift adjustment for log-coordinate in each component
            for(int componentIndex=firstLiborIndex; componentIndex < getNumberOfComponents() / 2; componentIndex++) {
                int spreadComponentIndex = componentIndex + getNumberOfComponents() / 2;
                RandomVariableInterface forwardVariance		= covarianceModel.getCovariance(timeIndex, componentIndex, componentIndex, realizationAtTimeIndex);
                RandomVariableInterface spreadVariance		= covarianceModel.getCovariance(timeIndex, spreadComponentIndex, spreadComponentIndex, realizationAtTimeIndex);
                drift[componentIndex] = drift[componentIndex].addProduct(forwardVariance, -0.5);
                drift[spreadComponentIndex] = drift[spreadComponentIndex].addProduct(spreadVariance, -0.5);
            }
        }

        return drift;
    }

    @Override
    public RandomVariableInterface[] getFactorLoading(int timeIndex, int componentIndex, RandomVariableInterface[] realizationAtTimeIndex) {
        return covarianceModel.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
    }

    public double getCorrelationFactor() {
        if(!(covarianceModel instanceof AbstractLIBORCovarianceModelParametric)) throw new UnsupportedOperationException();
        double[] parameter = ((AbstractLIBORCovarianceModelParametric)covarianceModel).getParameter();
        return parameter[parameter.length - 1];
    }

    public double getNormalizationFactor(int componentIndex, int maxTimeIndex) {
        if (c == null) c = new Double[getNumberOfComponents()][liborPeriodDiscretization.getNumberOfTimeSteps()];
        if (c[componentIndex][maxTimeIndex] == null) {
            double alpha = getCorrelationFactor();

            RandomVariableInterface squareIntegratedVolatility = new RandomVariable(0.0);
            for (int timeIndex = 0; timeIndex <= maxTimeIndex; timeIndex++) {
                RandomVariableInterface[] factorLoadings = covarianceModel.getFactorLoading(timeIndex, componentIndex, null);
                for (RandomVariableInterface factorLoading : factorLoadings) {
                    double dt = liborPeriodDiscretization.getTimeStep(timeIndex);
                    squareIntegratedVolatility = squareIntegratedVolatility.addProduct(factorLoading, factorLoading.mult(dt));
                }
            }

            c[componentIndex][maxTimeIndex] = Math.exp(-0.5 * alpha * (alpha - 1.0) * squareIntegratedVolatility.getAverage());
        }
        return c[componentIndex][maxTimeIndex];
    }

    public synchronized double[][][][] getIntegratedLIBORCovariances() {
        if (integratedLIBORCovariance != null) return integratedLIBORCovariance;

        int numberOfCovariances = multiCurveModel == MultiCurveModel.MMARTINGALE ? 5 : 3;
        int numberOfComponents = getNumberOfComponents() / 2;
        int numberOfTimeSteps = getTimeDiscretization().getNumberOfTimeSteps();

        integratedLIBORCovariance = new double[numberOfCovariances][numberOfTimeSteps][numberOfComponents][numberOfComponents];

        for (int timeIndex = 0; timeIndex < numberOfTimeSteps; timeIndex++) {
            double dt = getTime(timeIndex + 1) - getTime(timeIndex);
            RandomVariableInterface[][] factorLoadings = new RandomVariableInterface[getNumberOfComponents()][];
            // Prefetch factor loadings
            for(int componentIndex = 0; componentIndex < getNumberOfComponents(); componentIndex++) {
                factorLoadings[componentIndex] = getCovarianceModel().getFactorLoading(timeIndex, componentIndex, null);
            }

            for (int componentIndex1 = 0; componentIndex1 < numberOfComponents; componentIndex1++) {
                RandomVariableInterface[] factorLoadingComponent1Forward = factorLoadings[componentIndex1];
                RandomVariableInterface[] factorLoadingComponent1Spread = factorLoadings[componentIndex1 + numberOfComponents];

                for (int componentIndex2 = componentIndex1; componentIndex2 < numberOfComponents; componentIndex2++) {
                    if(getLiborPeriod(componentIndex1) > getTime(timeIndex)) {
                        RandomVariableInterface[] factorLoadingComponent2Forward = factorLoadings[componentIndex2];
                        RandomVariableInterface[] factorLoadingComponent2Spread = factorLoadings[componentIndex2 + numberOfComponents];

                        if (multiCurveModel == MultiCurveModel.MMARTINGALE) {
                            double c1 = getNormalizationFactor(componentIndex1, timeIndex);
                            double c2 = getNormalizationFactor(componentIndex2, timeIndex);
                            for(int factorIndex = 0; factorIndex < getNumberOfFactors(); factorIndex++) {
                                double covarianceForward = factorLoadingComponent1Forward[factorIndex].get(0) * factorLoadingComponent2Forward[factorIndex].get(0) * dt;
                                double covarianceSpread = factorLoadingComponent1Spread[factorIndex].get(0) * factorLoadingComponent2Spread[factorIndex].get(0) * dt;

                                integratedLIBORCovariance[0][timeIndex][componentIndex1][componentIndex2] += covarianceForward;
                                integratedLIBORCovariance[1][timeIndex][componentIndex1][componentIndex2] += covarianceForward * c1;
                                integratedLIBORCovariance[2][timeIndex][componentIndex1][componentIndex2] += covarianceForward * c2;
                                integratedLIBORCovariance[3][timeIndex][componentIndex1][componentIndex2] += covarianceForward * c1 * c2;
                                integratedLIBORCovariance[4][timeIndex][componentIndex1][componentIndex2] += covarianceSpread * c1 * c2;
                            }
                        } else {
                            for(int factorIndex = 0; factorIndex < getNumberOfFactors(); factorIndex++) {
                                integratedLIBORCovariance[0][timeIndex][componentIndex1][componentIndex2] +=
                                        factorLoadingComponent1Forward[factorIndex].get(0) * factorLoadingComponent2Forward[factorIndex].get(0) * dt;

                                integratedLIBORCovariance[1][timeIndex][componentIndex1][componentIndex2] +=
                                        factorLoadingComponent1Spread[factorIndex].get(0) * factorLoadingComponent2Spread[factorIndex].get(0) * dt;

                                integratedLIBORCovariance[2][timeIndex][componentIndex1][componentIndex2] +=
                                        (factorLoadingComponent1Forward[factorIndex].get(0) * factorLoadingComponent2Spread[factorIndex].get(0) +
                                                factorLoadingComponent2Forward[factorIndex].get(0) * factorLoadingComponent1Spread[factorIndex].get(0)) * dt;
                            }
                        }

                        // Integrate over time (i.e. sum up).
                        if(timeIndex > 0) for (int i = 0; i < numberOfCovariances; i++) {
                            integratedLIBORCovariance[i][timeIndex][componentIndex1][componentIndex2] =
                                    integratedLIBORCovariance[i][timeIndex-1][componentIndex1][componentIndex2] +
                                            integratedLIBORCovariance[i][timeIndex][componentIndex1][componentIndex2];
                            integratedLIBORCovariance[i][timeIndex][componentIndex2][componentIndex1] =
                                    integratedLIBORCovariance[i][timeIndex][componentIndex1][componentIndex2];
                        }
                    }
                }
            }
        }

        return integratedLIBORCovariance;
    }
}
