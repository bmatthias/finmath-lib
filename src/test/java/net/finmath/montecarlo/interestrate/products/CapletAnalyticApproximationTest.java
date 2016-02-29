package net.finmath.montecarlo.interestrate.products;

import cern.jet.random.engine.MersenneTwister64;
import net.finmath.functions.AnalyticFormulas;
import net.finmath.functions.NormalDistribution;
import org.junit.Test;

import java.util.concurrent.Callable;

import static org.junit.Assert.assertTrue;

public class CapletAnalyticApproximationTest {

    MersenneTwister64 mersenneTwister = new MersenneTwister64(3026);
    private double randomDrawSum(double correlation, double[] means, double[] covariances) {
        double stdNormalDist1 = NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());
        double stdNormalDist2 = NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());

        double correlatedNormalDist = correlation * stdNormalDist1 + Math.sqrt(1 - correlation * correlation) * stdNormalDist2;

        double normalDist1 = means[0] - 0.5 * covariances[0] + Math.sqrt(covariances[0]) * stdNormalDist1;
        double normalDist2 = means[1] - 0.5 * covariances[1] + Math.sqrt(covariances[1]) * correlatedNormalDist;

        return Math.exp(normalDist1) + Math.exp(normalDist2);
    }

    @Test
    public void testSimpleApproximation() throws Exception {
        int numberOfTests = 0;
        double numberOfSimulations = 1E5;
        double meanNormalizedErrorSimApprox = 0.0, meanNormalizedErrorBlackApprox = 0.0;
        for(int i = 1; i < 10; i++) {
            for(int j = 0; j < 10; j++) {
                for(int k = 0; k < 10; k++) {
                    double correlation = 0.1 * k;
                    double[] means = new double[] { -1.0 - 1.2 * i, -2.0 - 0.6 * i };
                    double[] covariances = new double[] { 0.1 + 0.5 * j, 0.6 * j, Math.sqrt(0.1 + 0.5 * j) * Math.sqrt(0.6 * j) * correlation };

                    double forward = Math.exp(means[0]);
                    double spread = Math.exp(means[1]);
                    double libor = forward + spread;
                    double strike = libor;
                    double[] moments = new double[] {
                            Math.log(libor),
                            (forward * forward * covariances[0] + spread * spread * covariances[1] + 2.0 * forward * spread * covariances[2]) / (libor * libor)
                    };

                    Callable<Double> randomDrawSimple = () -> {
                        double normalDist1 = moments[0] - 0.5 * moments[1] + Math.sqrt(moments[1]) * NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());
                        return Math.exp(normalDist1);
                    };

                    double probabilitySimApprox = 0.0, probabilitySum = 0.0;
                    for(int n = 0; n < numberOfSimulations; n++) {
                        probabilitySum += Math.max(randomDrawSum(correlation, means, covariances) - strike, 0.0);
                        probabilitySimApprox += Math.max(randomDrawSimple.call() - strike, 0.0);
                    }

                    probabilitySum /= numberOfSimulations;
                    probabilitySimApprox /= numberOfSimulations;
                    double probabilityBlackApprox = AnalyticFormulas.blackModelCapletValue(libor, Math.sqrt(moments[1]), 1.0, strike, 1.0, 1.0);

                    System.out.println("Sim Sum: \t\t" + probabilitySum);
                    System.out.println("Sim Approx.: \t" + probabilitySimApprox);
                    System.out.println("Black Approx.: \t" + probabilityBlackApprox);
                    System.out.println();

                    meanNormalizedErrorSimApprox += 2.0 * Math.abs(probabilitySum - probabilitySimApprox) / Math.abs(probabilitySum + probabilitySimApprox);
                    meanNormalizedErrorBlackApprox += 2.0 * Math.abs(probabilitySum - probabilityBlackApprox) / Math.abs(probabilitySum + probabilityBlackApprox);
                    numberOfTests++;
                }
            }
        }
        System.out.println("Mean Error Sim Approx.: \t\t" + (meanNormalizedErrorSimApprox /= numberOfTests));
        System.out.println("Mean Error Black  Approx.: \t\t" + (meanNormalizedErrorBlackApprox /= numberOfTests));

        assertTrue("The normalized mean error of the approximation in a Monte Carlo simulation is less than 5%.", meanNormalizedErrorSimApprox < 5e-2);
        assertTrue("The normalized mean error of the approximation by Black formula is less than 5%.", meanNormalizedErrorBlackApprox < 5e-2);
    }

    @Test
    public void testGetMomentsLevyJu() throws Exception {
        int numberOfTests = 0;
        double numberOfSimulations = 1E5;
        double meanNormalizedErrorSimLevy = 0.0, meanNormalizedErrorSimJu = 0.0, meanNormalizedErrorBlackLevy = 0.0, meanNormalizedErrorBlackJu = 0.0;
        for(int i = 1; i < 10; i++) {
            for(int j = 0; j < 10; j++) {
                for(int k = 0; k < 10; k++) {
                    for(int l = 10; l >= 10; l-=2) {
                        double correlation = 0.1 * k;
                        double[] means = new double[] { -1.0 - 1.2 * i, -2.0 - 0.6 * i };
                        double[] covariances = new double[] { 0.1 + 0.5 * j, 0.6 * j, Math.sqrt(0.1 + 0.5 * j) * Math.sqrt(0.6 * j) * correlation };
                        double strike = 0.1 * l * (Math.exp(means[0]) + Math.exp(means[1]));

                        double[][] covarianceMatrix = new double[][]{
                                { covariances[0], covariances[2] }, { covariances[2], covariances[1] }
                        };
                        double[] expMeans = new double[] { Math.exp(means[0]), Math.exp(means[1]) };
                        double[] momentsLevy = CapletAnalyticApproximation.getMomentsLevy(expMeans, covarianceMatrix);

                        Callable<Double> randomDrawLevy = () -> {
                            double normalDist1 = momentsLevy[0] + Math.sqrt(momentsLevy[1]) * NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());
                            return Math.exp(normalDist1);
                        };

                        double probabilitySimLevy = 0.0, probabilitySum = 0.0;
                        for(int n = 0; n < numberOfSimulations; n++) {
                            probabilitySum += Math.max(randomDrawSum(correlation, means, covariances) - strike, 0.0);
                            probabilitySimLevy += Math.max(randomDrawLevy.call() - strike, 0.0);
                        }

                        double forward = Math.exp(means[0]) + Math.exp(means[1]);
                        double adjustment = CapletAnalyticApproximation.getTaylorExpansionByJuAdjustment(expMeans, momentsLevy, covarianceMatrix, strike);
                        probabilitySum /= numberOfSimulations;
                        probabilitySimLevy /= numberOfSimulations;
                        double probabilitySimJu = probabilitySimLevy + adjustment;
                        double probabilityBlackLevy = CapletAnalyticApproximation.optionValue(forward, strike, momentsLevy[0], momentsLevy[1], 1.0);
                        double probabilityBlackJu = CapletAnalyticApproximation.optionValue(forward, strike, momentsLevy[0], momentsLevy[1], 1.0) + adjustment;

                        System.out.println("Sim Sum: \t\t" + probabilitySum);
                        System.out.println("Sim Levy: \t\t" + probabilitySimLevy);
                        System.out.println("Sim Ju: \t\t" + probabilitySimJu);
                        System.out.println("Black Levy: \t" + probabilityBlackLevy);
                        System.out.println("Black Ju: \t\t" + probabilityBlackJu);
                        System.out.println();

                        meanNormalizedErrorSimLevy += 2.0 * Math.abs(probabilitySum - probabilitySimLevy) / Math.abs(probabilitySum + probabilitySimLevy);
                        meanNormalizedErrorSimJu += 2.0 * Math.abs(probabilitySum - probabilitySimJu) / Math.abs(probabilitySum + probabilitySimJu);
                        meanNormalizedErrorBlackLevy += 2.0 * Math.abs(probabilitySum - probabilityBlackLevy) / Math.abs(probabilitySum + probabilityBlackLevy);
                        meanNormalizedErrorBlackJu += 2.0 * Math.abs(probabilitySum - probabilityBlackJu) / Math.abs(probabilitySum + probabilityBlackJu);
                        numberOfTests++;
                    }
                }
            }
        }
        System.out.println("Mean Error Sim Levy: \t\t" + (meanNormalizedErrorSimLevy /= numberOfTests));
        System.out.println("Mean Error Sim Ju: \t\t" + (meanNormalizedErrorSimJu /= numberOfTests));
        System.out.println("Mean Error Black Levy: \t\t" + (meanNormalizedErrorBlackLevy /= numberOfTests));
        System.out.println("Mean Error Black Ju: \t\t" + (meanNormalizedErrorBlackJu /= numberOfTests));

        assertTrue("The normalized mean error of the Levy approximation in a Monte Carlo simulation is less than 5%.", meanNormalizedErrorSimLevy < 5e-2);
        assertTrue("The normalized mean error of the Ju approximation in a Monte Carlo simulation is less than 5%.", meanNormalizedErrorSimJu < 5e-2);
        assertTrue("The normalized mean error of the Levy approximation by Black formula is less than 5%.", meanNormalizedErrorBlackLevy < 5e-2);
        assertTrue("The normalized mean error of the Ju approximation by Black formula is less than 5%.", meanNormalizedErrorBlackJu < 5e-2);
    }

    @Test
    public void testGetMomentsHo() throws Exception {
        for(int i = 1; i < 10; i++) {
            //Create two lognormals with parameters taken from the paper of Chia-Lu Ho
            double correlation = 0.0;
            double[] covariances = new double[] { 2.0, 2.0, 2.0 * correlation };
            double[] means = new double[] { -1.0 - i, -3.0 };
            double eta = i;// * Math.abs(means[1] - means[0]) / Math.sqrt(covariances[0] + covariances[1] + covariances[2]);
            double numberOfSimulations = 1E7;

            double probabilitySum = 0.0;
            for(int n = 0; n < numberOfSimulations; n++) {
                if(randomDrawSum(correlation, means, covariances) > eta) probabilitySum += 1;
            }

            double[] momentsHo = CapletAnalyticApproximation.getMomentsHo(new double[] { means[0] - 0.5 * covariances[0], means[1] - 0.5 * covariances[1] }, covariances);

            Callable<Double> randomDrawHo = () -> {
                double normalDist1 = momentsHo[0] + Math.sqrt(momentsHo[1]) * NormalDistribution.inverseCumulativeDistribution(mersenneTwister.nextDouble());
                return Math.exp(normalDist1);
            };

            double probabilityHo = 0.0;
            for(int n = 0; n < numberOfSimulations; n++) {
                if(randomDrawHo.call() > eta) probabilityHo += 1;
            }

            System.out.println("Sim Sum: \t" + probabilitySum / numberOfSimulations);
            System.out.println("Sim Ho: \t" + probabilityHo / numberOfSimulations);
        }
    }
}