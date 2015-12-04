package net.finmath.montecarlo.interestrate;

import net.finmath.marketdata.model.AnalyticModel;
import net.finmath.marketdata.model.curves.ForwardCurve;
import net.finmath.montecarlo.BrownianMotion;
import net.finmath.montecarlo.RandomVariable;
import net.finmath.montecarlo.interestrate.modelplugins.*;
import net.finmath.montecarlo.process.AbstractProcess;
import net.finmath.stochastic.RandomVariableInterface;
import net.finmath.time.TimeDiscretization;
import org.junit.Before;
import org.junit.Test;
import org.mockito.Mock;
import org.mockito.MockitoAnnotations;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;

import static org.junit.Assert.*;
import static org.mockito.Matchers.any;
import static org.mockito.Matchers.anyDouble;
import static org.mockito.Matchers.anyInt;
import static org.mockito.Matchers.eq;
import static org.mockito.Matchers.isNull;
import static org.mockito.Mockito.*;

public class MultiCurveLIBORMarketModelTest {
    private MultiCurveLIBORMarketModel additiveModel;
    private MultiCurveLIBORMarketModel multiplicativeModel;
    private MultiCurveLIBORMarketModel martingaleModel;

    @Mock ForwardCurve riskFreeCurve;
    @Mock ForwardCurve liborCurve;
    @Mock AnalyticModel analyticModel;
    @Mock AbstractLIBORCovarianceModelParametric covarianceModel;
    @Mock AbstractProcess process;

    private TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, 0.25, 0.5, 0.75, 1.0);
    private TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, 0.5, 1.0);

    @Before
    public void setUp() throws Exception {
        MockitoAnnotations.initMocks(this);

        when(process.getTimeDiscretization()).thenReturn(timeDiscretization);
        when(analyticModel.getForwardCurve("ois")).thenReturn(riskFreeCurve);
        when(analyticModel.getForwardCurve("libor")).thenReturn(liborCurve);
        when(covarianceModel.getParameter()).thenReturn(new double[]{ 10.0, 0.3 });
        when(riskFreeCurve.getForward(analyticModel, 0.0)).thenReturn(0.01);
        when(liborCurve.getForward(analyticModel, 0.0)).thenReturn(0.02);
        when(riskFreeCurve.getForward(analyticModel, 0.5)).thenReturn(0.03);
        when(liborCurve.getForward(analyticModel, 0.5)).thenReturn(0.05);

        Map<String, Object> properties = new HashMap<>();
        properties.put("multiCurveModel", MultiCurveLIBORMarketModel.MultiCurveModel.ADDITIVE);

        additiveModel = new MultiCurveLIBORMarketModel(
                liborPeriodDiscretization, analyticModel, "libor", "ois",
                covarianceModel, new LIBORMarketModelInterface.CalibrationItem[0], properties
        );

        properties.put("multiCurveModel", MultiCurveLIBORMarketModel.MultiCurveModel.MULTIPLICATIVE);
        multiplicativeModel = new MultiCurveLIBORMarketModel(
                liborPeriodDiscretization, analyticModel, "libor", "ois",
                covarianceModel, new LIBORMarketModelInterface.CalibrationItem[0], properties
        );

        properties.put("multiCurveModel", MultiCurveLIBORMarketModel.MultiCurveModel.MMARTINGALE);
        martingaleModel = new MultiCurveLIBORMarketModel(
                liborPeriodDiscretization, analyticModel, "libor", "ois",
                covarianceModel, new LIBORMarketModelInterface.CalibrationItem[0], properties
        );

        additiveModel.setProcess(process);
        multiplicativeModel.setProcess(process);
        martingaleModel.setProcess(process);
    }

    @Test
    public void testGetInitialStateAdditive() throws Exception {
        RandomVariableInterface[] initialState = additiveModel.getInitialState();

        verify(liborCurve, times(2)).getForward(eq(analyticModel), anyDouble());
        verify(riskFreeCurve, times(2)).getForward(eq(analyticModel), anyDouble());

        assertEquals(initialState.length, 4);
        assertEquals(Math.log(0.01), initialState[0].getAverage(), 1e-15);
        assertEquals(Math.log(0.03), initialState[1].getAverage(), 1e-15);
        assertEquals(Math.log(0.01), initialState[2].getAverage(), 1e-15);
        assertEquals(Math.log(0.02), initialState[3].getAverage(), 1e-15);
    }

    @Test
    public void testGetInitialStateMultiplicative() throws Exception {
        RandomVariableInterface[] initialState = multiplicativeModel.getInitialState();

        verify(liborCurve, times(2)).getForward(eq(analyticModel), anyDouble());
        verify(riskFreeCurve, times(2)).getForward(eq(analyticModel), anyDouble());

        assertEquals(initialState.length, 4);
        assertEquals(Math.log(0.01), initialState[0].getAverage(), 1e-15);
        assertEquals(Math.log(0.03), initialState[1].getAverage(), 1e-15);
        assertEquals(Math.log(0.01 / 1.005), initialState[2].getAverage(), 1e-15);
        assertEquals(Math.log(0.02 / 1.015), initialState[3].getAverage(), 1e-15);
    }

    @Test
    public void testGetInitialStateMartingale() throws Exception {
        RandomVariableInterface[] initialState = martingaleModel.getInitialState();

        verify(liborCurve, times(2)).getForward(eq(analyticModel), anyDouble());
        verify(riskFreeCurve, times(2)).getForward(eq(analyticModel), anyDouble());
        verify(covarianceModel, times(2)).getParameter();

        assertEquals(initialState.length, 4);
        assertEquals(Math.log(0.01), initialState[0].getAverage(), 1e-15);
        assertEquals(Math.log(0.03), initialState[1].getAverage(), 1e-15);
        assertEquals(Math.log(0.01 / 1.005 / Math.pow(0.01, 0.3)), initialState[2].getAverage(), 1e-15);
        assertEquals(Math.log(0.02 / 1.015 / Math.pow(0.03, 0.3)), initialState[3].getAverage(), 1e-15);
    }

    @Test
    public void testGetIntegratedLIBORCovariancesAdditiveAndMultiplicative() throws Exception {
        when(covarianceModel.getFactorLoading(anyInt(), anyInt(), (RandomVariable[]) isNull()))
                .thenReturn(new RandomVariableInterface[]{
                        new RandomVariable(1), new RandomVariable(3), new RandomVariable(6)
                });

        //test that the third factor with value 6 is ignored
        when(process.getNumberOfFactors()).thenReturn(2);

        //dt = 1.0 for all t
        when(process.getTime(0)).thenReturn(0.0);
        when(process.getTime(1)).thenReturn(1.0);
        when(process.getTime(2)).thenReturn(2.0);
        when(process.getTime(3)).thenReturn(3.0);

        double[][][][] covariances = additiveModel.getIntegratedLIBORCovariances();

        assertEquals(covariances.length, 3);
        assertEquals(covariances[0].length, timeDiscretization.getNumberOfTimeSteps());
        assertEquals(covariances[0][0].length, liborPeriodDiscretization.getNumberOfTimeSteps());
        assertEquals(covariances[0][0][0].length, liborPeriodDiscretization.getNumberOfTimeSteps());

        assertEquals(10, (int)covariances[0][0][1][1]);

        //integrate for libor period > time TODO: Should it be >=?
        assertEquals(0, (int)covariances[0][1][0][0]);

        assertTrue(Arrays.deepEquals(covariances, multiplicativeModel.getIntegratedLIBORCovariances()));
    }

    @Test
    public void testGetDriftAdditive() throws Exception {
        double forward = 1.0;
        double periodLength = liborPeriodDiscretization.getTimeStep(0);
        double factorLoading1 = 1.0, factorLoading2 = 3.0;
        double covariance = 3.0;

        when(covarianceModel.getFactorLoading(anyInt(), anyInt(), any()))
                .thenReturn(new RandomVariableInterface[]{
                        new RandomVariable(factorLoading1), new RandomVariable(factorLoading2), new RandomVariable(6)
                });
        when(covarianceModel.getCovariance(anyInt(), anyInt(), anyInt(), any())).thenReturn(new RandomVariable(covariance));

        BrownianMotion brownianMotion = mock(BrownianMotion.class);
        when(process.getBrownianMotion()).thenReturn(brownianMotion);
        when(process.getNumberOfFactors()).thenReturn(2);
        when(brownianMotion.getRandomVariableForConstant(anyDouble())).thenReturn(new RandomVariable(periodLength));

        RandomVariableInterface[] realizationAtTimeIndex = new RandomVariableInterface[]{
                new RandomVariable(forward), new RandomVariable(forward)
        };

        RandomVariableInterface[] drift = additiveModel.getDrift(0, realizationAtTimeIndex, null);
        RandomVariableInterface[] driftSimple = calculateDriftSlow(0, realizationAtTimeIndex, additiveModel);

        assertNull(drift[0]);
        assertNull(drift[2]);
        assertTrue(drift[1].equals(driftSimple[1]));
        assertTrue(drift[3].equals(driftSimple[3]));

        //we are looking at only one period, so we can easily calculate the drift from the chosen mock values
        double driftForward = periodLength * forward / (1.0 + periodLength * forward) *
                (Math.pow(factorLoading1, 2) + Math.pow(factorLoading2, 2)) - covariance * 0.5;

        assertEquals(drift[1].get(0), driftForward, 1e-15);
        assertEquals(drift[3].get(0), driftForward, 1e-15);
    }

    @Test
    public void testGetDriftMultiplicative() throws Exception {
        double forward = 1.0;
        double periodLength = liborPeriodDiscretization.getTimeStep(0);
        double factorLoading1 = 1.0, factorLoading2 = 3.0;
        double covariance = 3.0;

        when(covarianceModel.getFactorLoading(anyInt(), anyInt(), any()))
                .thenReturn(new RandomVariableInterface[]{
                        new RandomVariable(factorLoading1), new RandomVariable(factorLoading2), new RandomVariable(6)
                });
        when(covarianceModel.getCovariance(anyInt(), anyInt(), anyInt(), any())).thenReturn(new RandomVariable(covariance));

        BrownianMotion brownianMotion = mock(BrownianMotion.class);
        when(process.getBrownianMotion()).thenReturn(brownianMotion);
        when(process.getNumberOfFactors()).thenReturn(2);
        when(brownianMotion.getRandomVariableForConstant(anyDouble())).thenReturn(new RandomVariable(periodLength));

        RandomVariableInterface[] realizationAtTimeIndex = new RandomVariableInterface[]{
                new RandomVariable(forward), new RandomVariable(forward)
        };

        RandomVariableInterface[] drift = multiplicativeModel.getDrift(0, realizationAtTimeIndex, null);
        RandomVariableInterface[] driftSimple = calculateDriftSlow(0, realizationAtTimeIndex, multiplicativeModel);

        assertNull(drift[0]);
        assertNull(drift[2]);
        assertTrue(drift[1].equals(driftSimple[1]));
        assertTrue(drift[3].equals(driftSimple[3]));

        //we are looking at only one period, so we can easily calculate the drift from the chosen mock values
        double driftForward = periodLength * forward / (1.0 + periodLength * forward) *
                (Math.pow(factorLoading1, 2) + Math.pow(factorLoading2, 2)) - covariance * 0.5;
        double driftSpread = 0.0 - covariance * 0.5;

        assertEquals(drift[1].get(0), driftForward, 1e-15);
        assertEquals(drift[3].get(0), driftSpread, 1e-15);
    }

    @Test
    public void testGetDriftAdditiveWithoutMock() throws Exception {
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, 0.5, 1.0, 1.5, 2.0);
        TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, 0.5, 1.0, 1.5, 2.0);

        AbstractMultiCurveLIBORCorrelationModel correlationModel = new LIBORCorrelationModelSeperateCurveFiveParameterExponentialDecay(
                timeDiscretization, liborPeriodDiscretization, 8, new double[]{ 0.1, 0.1, 0.1 }, false
        );
        LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTwoCurveEightParameterExponentialForm(
                timeDiscretization, liborPeriodDiscretization, new double[]{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }, false
        );
        covarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(
                timeDiscretization, liborPeriodDiscretization, volatilityModel, correlationModel
        );
        Map<String, Object> properties = new HashMap<>();
        properties.put("multiCurveModel", MultiCurveLIBORMarketModel.MultiCurveModel.ADDITIVE);
        MultiCurveLIBORMarketModel model = new MultiCurveLIBORMarketModel(
                liborPeriodDiscretization, analyticModel, "libor", "ois",
                covarianceModel, new LIBORMarketModelInterface.CalibrationItem[0], properties
        );
        model.setProcess(process);

        BrownianMotion brownianMotion = mock(BrownianMotion.class);
        when(process.getBrownianMotion()).thenReturn(brownianMotion);
        when(process.getNumberOfFactors()).thenReturn(8);
        when(process.getTimeDiscretization()).thenReturn(timeDiscretization);
        when(brownianMotion.getRandomVariableForConstant(anyDouble())).thenReturn(new RandomVariable(0.5));

        RandomVariableInterface[] realizationAtTimeIndex = new RandomVariableInterface[]{
                new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0),
                new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0)
        };

        RandomVariableInterface[] drift = model.getDrift(0, realizationAtTimeIndex, null);
        RandomVariableInterface[] driftSimple = calculateDriftSlow(0, realizationAtTimeIndex, model);

        for (int i = 0; i < drift.length; i++) {
            assertTrue((drift[i] == null && driftSimple[i] == null) || drift[i].equals(driftSimple[i]));
        }
    }

    @Test
    public void testGetDriftMultiplicativeWithoutMock() throws Exception {
        TimeDiscretization timeDiscretization = new TimeDiscretization(0.0, 0.5, 1.0, 1.5, 2.0);
        TimeDiscretization liborPeriodDiscretization = new TimeDiscretization(0.0, 0.5, 1.0, 1.5, 2.0);

        AbstractMultiCurveLIBORCorrelationModel correlationModel = new LIBORCorrelationModelSeperateCurveFiveParameterExponentialDecay(
                timeDiscretization, liborPeriodDiscretization, 8, new double[]{ 0.1, 0.1, 0.1 }, false
        );
        LIBORVolatilityModel volatilityModel = new LIBORVolatilityModelTwoCurveEightParameterExponentialForm(
                timeDiscretization, liborPeriodDiscretization, new double[]{ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 }, false
        );
        covarianceModel = new LIBORCovarianceModelFromVolatilityAndCorrelation(
                timeDiscretization, liborPeriodDiscretization, volatilityModel, correlationModel
        );
        Map<String, Object> properties = new HashMap<>();
        properties.put("multiCurveModel", MultiCurveLIBORMarketModel.MultiCurveModel.ADDITIVE);
        MultiCurveLIBORMarketModel model = new MultiCurveLIBORMarketModel(
                liborPeriodDiscretization, analyticModel, "libor", "ois",
                covarianceModel, new LIBORMarketModelInterface.CalibrationItem[0], properties
        );
        model.setProcess(process);

        BrownianMotion brownianMotion = mock(BrownianMotion.class);
        when(process.getBrownianMotion()).thenReturn(brownianMotion);
        when(process.getNumberOfFactors()).thenReturn(8);
        when(process.getTimeDiscretization()).thenReturn(timeDiscretization);
        when(brownianMotion.getRandomVariableForConstant(anyDouble())).thenReturn(new RandomVariable(0.5));

        RandomVariableInterface[] realizationAtTimeIndex = new RandomVariableInterface[]{
                new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0),
                new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0), new RandomVariable(1.0)
        };

        RandomVariableInterface[] drift = model.getDrift(0, realizationAtTimeIndex, null);
        RandomVariableInterface[] driftSimple = calculateDriftSlow(0, realizationAtTimeIndex, model);

        for (int i = 0; i < drift.length; i++) {
            assertTrue((drift[i] == null && driftSimple[i] == null) || drift[i].equals(driftSimple[i]));
        }
    }

    //this is a different, straightforward but slower method of calculating the drift
    //should give the same result as the implementation of the model
    //this is only for SPOT measure and log-normal process
    private RandomVariableInterface[] calculateDriftSlow(int timeIndex, RandomVariableInterface[] realizationAtTimeIndex, MultiCurveLIBORMarketModel model) {
        int numberOfComponents = model.getNumberOfComponents();
        int numberOfFactors = model.getNumberOfFactors();
        MultiCurveLIBORMarketModel.MultiCurveModel multiCurveModel = model.getMultiCurveModel();

        double time = model.getTime(timeIndex);
        int		firstLiborIndex		= model.getLiborPeriodIndex(time) + 1;
        if(firstLiborIndex<0) firstLiborIndex = -firstLiborIndex;

        RandomVariableInterface zero	= new RandomVariable(time, 0.0);

        // Allocate drift vector and initialize to zero (will be used to sum up drift components)
        RandomVariableInterface[]	drift = new RandomVariableInterface[numberOfComponents];
        for(int componentIndex = firstLiborIndex; componentIndex < numberOfComponents / 2; componentIndex++) {

            drift[componentIndex] = zero;
            drift[componentIndex + numberOfComponents / 2] = zero;

            int spreadComponentIndex = componentIndex + numberOfComponents / 2;
            RandomVariableInterface[]	forwardFactorLoading   	= covarianceModel.getFactorLoading(timeIndex, componentIndex, realizationAtTimeIndex);
            RandomVariableInterface[]	spreadFactorLoading   	= covarianceModel.getFactorLoading(timeIndex, spreadComponentIndex, realizationAtTimeIndex);

            for (int componentIndex2 = firstLiborIndex; componentIndex2 <= componentIndex; componentIndex2++) {
                double periodLength	= model.getLiborPeriodDiscretization().getTimeStep(componentIndex);
                RandomVariableInterface forward			= realizationAtTimeIndex[componentIndex];
                RandomVariableInterface oneStepMeasureTransform = (process.getBrownianMotion()
                        .getRandomVariableForConstant(periodLength))
                        .discount(forward, periodLength)
                        .mult(forward);

                RandomVariableInterface[] forwardFactorLoading2 = covarianceModel.getFactorLoading(timeIndex, componentIndex2, realizationAtTimeIndex);

                for(int factorIndex = 0; factorIndex < numberOfFactors; factorIndex++) {
                    drift[componentIndex] =  drift[componentIndex].addProduct(oneStepMeasureTransform.mult(forwardFactorLoading2[factorIndex]), forwardFactorLoading[factorIndex]);
                    if(componentIndex2 < componentIndex || multiCurveModel == MultiCurveLIBORMarketModel.MultiCurveModel.ADDITIVE)
                        drift[spreadComponentIndex] = drift[spreadComponentIndex].addProduct(oneStepMeasureTransform.mult(forwardFactorLoading2[factorIndex]), spreadFactorLoading[factorIndex]);
                }
            }
        }

        for(int componentIndex=firstLiborIndex; componentIndex < numberOfComponents / 2; componentIndex++) {
            int spreadComponentIndex = componentIndex + numberOfComponents / 2;
            RandomVariableInterface forwardVariance		= covarianceModel.getCovariance(timeIndex, componentIndex, componentIndex, realizationAtTimeIndex);
            RandomVariableInterface spreadVariance		= covarianceModel.getCovariance(timeIndex, spreadComponentIndex, spreadComponentIndex, realizationAtTimeIndex);
            drift[componentIndex] = drift[componentIndex].addProduct(forwardVariance, -0.5);
            drift[spreadComponentIndex] = drift[spreadComponentIndex].addProduct(spreadVariance, -0.5);
        }

        return drift;
    }
}