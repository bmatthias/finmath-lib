/*
 * (c) Copyright Christian P. Fries, Germany. All rights reserved. Contact: email@christian-fries.de.
 *
 * Created on 20.05.2006
 */
package net.finmath.montecarlo.interestrate.modelplugins;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.BrownianMotionInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.LIBORModelMonteCarloSimulation;
import net.finmath.montecarlo.interestrate.products.AbstractLIBORMonteCarloProduct;
import net.finmath.montecarlo.process.ProcessEulerScheme;
import net.finmath.optimizer.LevenbergMarquardt;
import net.finmath.optimizer.SolverException;
import net.finmath.time.TimeDiscretizationInterface;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Base class for parametric covariance models, see also {@link AbstractLIBORCovarianceModel}.
 * 
 * Parametric models feature a parameter vector which can be inspected
 * and modified for calibration purposes.
 * 
 * The parameter vector may have zero length, which indicated that the model
 * is not calibrateable.
 * 
 * This class includes the implementation of a generic calibration algorithm.
 * If you provide an arbitrary list of calibration products, the class can return
 * a new instance where the parameters are chosen such that the (weighted) root-mean-square 
 * error of the difference of the value of the calibration products and given target
 * values is minimized.
 * 
 * @author Christian Fries
 * @date 20.05.2006
 * @date 23.02.2014
 * @version 1.1
 */
public abstract class AbstractLIBORCovarianceModelParametric extends AbstractLIBORCovarianceModel {

	private static final Logger logger = Logger.getLogger("net.finmath");

	/**
	 * Constructor consuming time discretizations, which are handled by the super class.
	 * 
	 * @param timeDiscretization The vector of simulation time discretization points.
	 * @param liborPeriodDiscretization The vector of tenor discretization points.
	 * @param numberOfFactors The number of factors to use (a factor reduction is performed)
	 */
	public AbstractLIBORCovarianceModelParametric(TimeDiscretizationInterface timeDiscretization, TimeDiscretizationInterface liborPeriodDiscretization, int numberOfFactors) {
		super(timeDiscretization, liborPeriodDiscretization, numberOfFactors);
	}

	/**
	 * Get the parameters of determining this parametric
	 * covariance model. The parameters are usually free parameters
	 * which may be used in calibration.
	 * 
	 * @return Parameter vector.
	 */
	public abstract double[]	getParameter();

	public abstract void		setParameter(double[] parameter);

	@Override
	public abstract Object clone();

	public AbstractLIBORCovarianceModelParametric getCloneWithModifiedParameters(double[] parameters) {
		AbstractLIBORCovarianceModelParametric calibrationCovarianceModel = (AbstractLIBORCovarianceModelParametric)AbstractLIBORCovarianceModelParametric.this.clone();
		calibrationCovarianceModel.setParameter(parameters);

		return calibrationCovarianceModel;
	}

	public AbstractLIBORCovarianceModelParametric getCloneCalibrated(final LIBORMarketModelInterface calibrationModel, final AbstractLIBORMonteCarloProduct[] calibrationProducts, double[] calibrationTargetValues, double[] calibrationWeights) throws CalculationException {
		return getCloneCalibrated(calibrationModel, calibrationProducts, calibrationTargetValues, calibrationWeights, null);
	}

	/**
	 * Performs a generic calibration of the parametric model by
	 * trying to match a given vector of calibration product to a given vector of target values
	 * using a given vector of weights.
	 * 
	 * Optional calibration parameters may be passed using the map calibrationParameters. The keys are (<code>String</code>s):
	 * <ul>
	 * 	<li><tt>brownianMotion</tt>: Under this key an object implementing {@link net.finmath.montecarlo.BrownianMotionInterface} may be provided. If so, this Brownian motion is used to build the valuation model.</li>
	 * 	<li><tt>maxIterations</tt>: Under this key an object of type Integer may be provided specifying the maximum number of iterations.</li>
	 * 	<li><tt>accuracy</tt>: Under this key an object of type Double may be provided specifying the desired accuracy. Note that this is understood in the sense that the solver will stop if the iteration does not improve by more than this number.</li>
	 * </ul>
	 * 
	 * @param calibrationModel The LIBOR market model to be used for calibrations (specifies forward curve and tenor discretization).
	 * @param calibrationProducts The array of calibration products.
	 * @param calibrationTargetValues The array of target values.
	 * @param calibrationWeights The array of weights.
	 * @param calibrationParameters A map of type Map&lt;String, Object&gt; specifying some (optional) calibration parameters.
	 * @return A new parametric model of the same type than <code>this</code> one, but with calibrated parameters.
	 * @throws CalculationException Thrown if calibration has failed.
	 */
	public AbstractLIBORCovarianceModelParametric getCloneCalibrated(final LIBORMarketModelInterface calibrationModel, final AbstractLIBORMonteCarloProduct[] calibrationProducts, final double[] calibrationTargetValues, double[] calibrationWeights, Map<String,Object> calibrationParameters) throws CalculationException {

		double[] initialParameters = this.getParameter();

		if(calibrationParameters == null) calibrationParameters = new HashMap<String,Object>();
		Integer numberOfPathsParameter	= (Integer)calibrationParameters.get("numberOfPaths");
		Integer seedParameter			= (Integer)calibrationParameters.get("seed");
		Integer maxIterationsParameter	= (Integer)calibrationParameters.get("maxIterations");
		Double	accuracyParameter		= (Double)calibrationParameters.get("accuracy");
		BrownianMotionInterface brownianMotionParameter	= (BrownianMotionInterface)calibrationParameters.get("brownianMotion");

		int numberOfPaths	= numberOfPathsParameter != null ? numberOfPathsParameter.intValue() : 2000;
		int seed			= seedParameter != null ? seedParameter.intValue() : 31415;
		int maxIterations	= maxIterationsParameter != null ? maxIterationsParameter.intValue() : 400;
		double accuracy		= accuracyParameter != null ? accuracyParameter.doubleValue() : 1E-7;

		final BrownianMotionInterface brownianMotion = brownianMotionParameter != null ? brownianMotionParameter : new BrownianMotion(getTimeDiscretization(), getNumberOfFactors(), numberOfPaths, seed);

		/*
		 * Create separate thread pools for values and derivatives, in order to prevent deadlocks.
		 * We allow for getParameter().length simultaneous calibration models.
		 * Note: In the case of a Monte-Carlo calibration, the memory requirement is that of
		 * one model with getParameter().length times the number of paths. In the case of an analytic calibration
		 * memory requirement is not the limiting factor.
		 * Note: It is important not to re-create the thread pools on every run, in order to avoid performance issues
		 * from threading overhead on machines with a large number of processors. This is especially a problem when the
		 * actual calculations are much faster than the creation of the threads, i.e. with analytic calibration.
		 */
		final ExecutorService executorValues = Executors.newFixedThreadPool(Math.min(Math.max(Runtime.getRuntime().availableProcessors(), 1), calibrationProducts.length));
		final ExecutorService executorDerivatives = Executors.newFixedThreadPool(Math.min(Math.max(Runtime.getRuntime().availableProcessors(), 1), getParameter().length));
		LevenbergMarquardt optimizer = new LevenbergMarquardt(initialParameters, calibrationTargetValues, maxIterations, executorDerivatives) {
			// Calculate model values for given parameters
			@Override
			public void setValues(double[] parameters, double[] values) throws SolverException {

				AbstractLIBORCovarianceModelParametric calibrationCovarianceModel = AbstractLIBORCovarianceModelParametric.this.getCloneWithModifiedParameters(parameters);

				// Create a LIBOR market model with the new covariance structure.
				LIBORMarketModelInterface model = calibrationModel.getCloneWithModifiedCovarianceModel(calibrationCovarianceModel);
				ProcessEulerScheme process = new ProcessEulerScheme(brownianMotion);
				final LIBORModelMonteCarloSimulation liborMarketModelMonteCarloSimulation =  new LIBORModelMonteCarloSimulation(model, process);

				ArrayList<Future<Double>> valueFutures = new ArrayList<Future<Double>>(calibrationProducts.length);

				for(int calibrationProductIndex=0; calibrationProductIndex<calibrationProducts.length; calibrationProductIndex++) {
					final int workerCalibrationProductIndex = calibrationProductIndex;
					Callable<Double> worker = () -> {
                        try {
                            return calibrationProducts[workerCalibrationProductIndex].getValue(liborMarketModelMonteCarloSimulation);
                        } catch (CalculationException e) {
							e.printStackTrace();
                            // We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
                            return calibrationTargetValues[workerCalibrationProductIndex];
                        } catch (Exception e) {
							e.printStackTrace();
                            // We do not signal exceptions to keep the solver working and automatically exclude non-working calibration products.
                            return calibrationTargetValues[workerCalibrationProductIndex];
                        }
                    };
                    Future<Double> valueFuture = executorValues.submit(worker);
                    valueFutures.add(calibrationProductIndex, valueFuture);
				}
				for(int calibrationProductIndex=0; calibrationProductIndex<calibrationProducts.length; calibrationProductIndex++) {
					try {
						double value = valueFutures.get(calibrationProductIndex).get();
						values[calibrationProductIndex] = value;
					} catch (InterruptedException e) {
						throw new SolverException(e);
					} catch (ExecutionException e) {
						throw new SolverException(e);
					}
				}
			}
		};

		// Set solver parameters
		optimizer.setWeights(calibrationWeights);
		optimizer.setErrorTolerance(accuracy);

		try {
			optimizer.run();
		} catch(SolverException e) {
			throw new CalculationException(e);
		} finally {
			executorValues.shutdown();
			executorDerivatives.shutdown();
		}

		// Get covariance model corresponding to the best parameter set.
		double[] bestParameters = optimizer.getBestFitParameters();
		AbstractLIBORCovarianceModelParametric calibrationCovarianceModel = this.getCloneWithModifiedParameters(bestParameters);

		// Diagnostic output
		if (logger.isLoggable(Level.FINE)) {
			logger.fine("The solver required " + optimizer.getIterations() + " iterations. The best fit parameters are:");

			String logString = "Best parameters:";
			for(int i=0; i<bestParameters.length; i++) {
				logString += "\tparameter["+i+"]: " + bestParameters[i];
			}
			logger.fine(logString);
		}

		return calibrationCovarianceModel;    	
	}

	@Override
	public String toString() {
		return "AbstractLIBORCovarianceModelParametric [getParameter()="
				+ Arrays.toString(getParameter()) + "]";
	}
}
