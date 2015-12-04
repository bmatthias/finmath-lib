package net.finmath.tests.montecarlo.interestrate.tools;

import net.finmath.exception.CalculationException;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.montecarlo.interestrate.MultiCurveLIBORMarketModel;
import net.finmath.montecarlo.process.ProcessInterface;
import net.finmath.stochastic.RandomVariableInterface;

import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Locale;

public class CSVCreator {
    private static DecimalFormat formatterValue		= new DecimalFormat("##0.000;-##0.000", new DecimalFormatSymbols(Locale.ENGLISH));
    private static DecimalFormat formatterDeviation	= new DecimalFormat("0.00000E00;-0.00000E00", new DecimalFormatSymbols(Locale.ENGLISH));

    public static void writeParamsAndDuration(double[] scParams, long duration, double[] mcParams, FileWriter writer) throws IOException {
        writer.write("Calibration Time" + "\t" + "Parameter");
        for (double p : scParams) writer.write("\n" + duration + "\t" + p);
        for (double p : mcParams) writer.write("\n" + duration + "\t" + p);
    }

    public static void writeErrorStatistics(LIBORMarketModelInterface.CalibrationItem[] eoniaCalibrationItems, LIBORMarketModelInterface.CalibrationItem[] euriborCalibrationItems, FileWriter writer, double[] deviations) throws IOException {
        writer.write("Error Type"
                        + "\t" + "Vola Multi-Curve Analytic"
                        + "\t" + "Vola Multi-Curve Monte-Carlo (Spot)"
                        + "\t" + "Vola Multi-Curve Monte-Carlo (Terminal)"
                        + "\t" + "Vola Single-Curve Analytic"
                        + "\t" + "Vola Single-Curve Monte-Carlo (Spot)"
                        + "\t" + "Vola Single-Curve Monte-Carlo (Terminal)"
                        + "\t" + "Value Multi-Curve Analytic"
                        + "\t" + "Value Multi-Curve Monte-Carlo (Spot)"
                        + "\t" + "Value Multi-Curve Monte-Carlo (Terminal)"
                        + "\t" + "Value Single-Curve Analytic"
                        + "\t" + "Value Single-Curve Monte-Carlo (Spot)"
                        + "\t" + "Value Single-Curve Monte-Carlo (Terminal)"
                        + "\n" + "Max Deviation"
                        + "\t" + formatterDeviation.format(deviations[0])
                        + "\t" + formatterDeviation.format(deviations[1])
                        + "\t" + formatterDeviation.format(deviations[2])
                        + "\t" + formatterDeviation.format(deviations[3])
                        + "\t" + formatterDeviation.format(deviations[4])
                        + "\t" + formatterDeviation.format(deviations[5])
                        + "\t" + formatterDeviation.format(deviations[6])
                        + "\t" + formatterDeviation.format(deviations[7])
                        + "\t" + formatterDeviation.format(deviations[8])
                        + "\t" + formatterDeviation.format(deviations[9])
                        + "\t" + formatterDeviation.format(deviations[10])
                        + "\t" + formatterDeviation.format(deviations[11])
                        + "\n" + "Mean Deviation"
                        + "\t" + formatterDeviation.format(deviations[12] / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[13] / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[14] / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[15] / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[16] / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[17] / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[18] / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[19] / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[20] / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[21] / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[22] / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(deviations[23] / eoniaCalibrationItems.length)
                        + "\n" + "Mean Square Deviation"
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[24]) / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[25]) / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[26]) / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[27]) / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[28]) / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[29]) / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[30]) / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[31]) / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[32]) / euriborCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[33]) / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[34]) / eoniaCalibrationItems.length)
                        + "\t" + formatterDeviation.format(Math.sqrt(deviations[35]) / eoniaCalibrationItems.length)
        );
    }

    public static void writeValues(FileWriter writer, double[][] results, double[] deviations, int calibrationItemIndex) throws IOException {
        if (calibrationItemIndex == 0) {
            writer.write("Vola Target (Multi-Curve)" +
                            "\t" + "Vola Analytic" +
                            "\t" + "Vola Monte-Carlo (Spot)" +
                            "\t" + "Vola Monte-Carlo (Terminal)" +
                            "\t" + "Vola Target (Single Curve)" +
                            "\t" + "Vola Analytic (Single Curve)" +
                            "\t" + "Vola Monte-Carlo (Single Curve Spot)" +
                            "\t" + "Vola Monte-Carlo (Single Curve Terminal)" +
                            "\t" + "Value Target (Multi-Curve)" +
                            "\t" + "Value Analytic" +
                            "\t" + "Value Monte-Carlo (Spot)" +
                            "\t" + "Value Monte-Carlo (Terminal)" +
                            "\t" + "Value Target (Single-Curve)" +
                            "\t" + "Value Analytic (Single Curve)" +
                            "\t" + "Value Monte-Carlo (Single Curve Spot)" +
                            "\t" + "Value Monte-Carlo (Single Curve Terminal)" + "\n"
            );
        }
        writer.write(formatterValue.format(Math.min(Math.max(100 * results[0][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[2][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[10][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[14][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[1][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[3][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[11][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[15][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[6][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[4][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[8][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[12][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[7][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[5][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[9][calibrationItemIndex], 0), 999))
                        + "\t" + formatterValue.format(Math.min(Math.max(100 * results[13][calibrationItemIndex], 0), 999)) + "\n"
        );

        double deviationVolaMultiCurveAnalytic = results[2][calibrationItemIndex] - results[0][calibrationItemIndex];
        double deviationVolaMultiCurveMonteCarloSpot = results[10][calibrationItemIndex] - results[0][calibrationItemIndex];
        double deviationVolaMultiCurveMonteCarloTerminal = results[14][calibrationItemIndex] - results[0][calibrationItemIndex];

        double deviationVolaSingleCurveAnalytic = results[3][calibrationItemIndex] - results[1][calibrationItemIndex];
        double deviationVolaSingleCurveMonteCarloSpot = results[11][calibrationItemIndex] - results[1][calibrationItemIndex];
        double deviationVolaSingleCurveMonteCarloTerminal = results[15][calibrationItemIndex] - results[1][calibrationItemIndex];

        double deviationValueMultiCurveAnalytic = results[4][calibrationItemIndex] - results[6][calibrationItemIndex];
        double deviationValueMultiCurveMonteCarloSpot = results[8][calibrationItemIndex] - results[6][calibrationItemIndex];
        double deviationValueMultiCurveMonteCarloTerminal = results[12][calibrationItemIndex] - results[6][calibrationItemIndex];

        double deviationValueSingleCurveAnalytic = results[5][calibrationItemIndex] - results[7][calibrationItemIndex];
        double deviationValueSingleCurveMonteCarloSpot = results[9][calibrationItemIndex] - results[7][calibrationItemIndex];
        double deviationValueSingleCurveMonteCarloTerminal = results[13][calibrationItemIndex] - results[7][calibrationItemIndex];

        deviations[0] = Math.max(Math.abs(deviationVolaMultiCurveAnalytic), deviations[0]);
        deviations[1] = Math.max(Math.abs(deviationVolaMultiCurveMonteCarloSpot), deviations[1]);
        deviations[2] = Math.max(Math.abs(deviationVolaMultiCurveMonteCarloTerminal), deviations[2]);

        deviations[3] = Math.max(Math.abs(deviationVolaSingleCurveAnalytic), deviations[3]);
        deviations[4] = Math.max(Math.abs(deviationVolaSingleCurveMonteCarloSpot), deviations[4]);
        deviations[5] = Math.max(Math.abs(deviationVolaSingleCurveMonteCarloTerminal), deviations[5]);

        deviations[6] = Math.max(Math.abs(deviationValueMultiCurveAnalytic), deviations[6]);
        deviations[7] = Math.max(Math.abs(deviationValueMultiCurveMonteCarloSpot), deviations[7]);
        deviations[8] = Math.max(Math.abs(deviationValueMultiCurveMonteCarloTerminal), deviations[8]);

        deviations[9] = Math.max(Math.abs(deviationValueSingleCurveAnalytic), deviations[9]);
        deviations[10] = Math.max(Math.abs(deviationValueSingleCurveMonteCarloSpot), deviations[10]);
        deviations[11] = Math.max(Math.abs(deviationValueSingleCurveMonteCarloTerminal), deviations[11]);

        deviations[12] += Math.abs(deviationVolaMultiCurveAnalytic);
        deviations[13] += Math.abs(deviationVolaMultiCurveMonteCarloSpot);
        deviations[14] += Math.abs(deviationVolaMultiCurveMonteCarloTerminal);

        deviations[15] += Math.abs(deviationVolaSingleCurveAnalytic);
        deviations[16] += Math.abs(deviationVolaSingleCurveMonteCarloSpot);
        deviations[17] += Math.abs(deviationVolaSingleCurveMonteCarloTerminal);

        deviations[18] += Math.abs(deviationValueMultiCurveAnalytic);
        deviations[19] += Math.abs(deviationValueMultiCurveMonteCarloSpot);
        deviations[20] += Math.abs(deviationValueMultiCurveMonteCarloTerminal);

        deviations[21] += Math.abs(deviationValueSingleCurveAnalytic);
        deviations[22] += Math.abs(deviationValueSingleCurveMonteCarloSpot);
        deviations[23] += Math.abs(deviationValueSingleCurveMonteCarloTerminal);

        deviations[24] += deviationVolaMultiCurveAnalytic * deviationVolaMultiCurveAnalytic;
        deviations[25] += deviationVolaMultiCurveMonteCarloSpot * deviationVolaMultiCurveMonteCarloSpot;
        deviations[26] += deviationVolaMultiCurveMonteCarloTerminal * deviationVolaMultiCurveMonteCarloTerminal;

        deviations[27] += deviationVolaSingleCurveAnalytic * deviationVolaMultiCurveMonteCarloTerminal;
        deviations[28] += deviationVolaSingleCurveMonteCarloSpot * deviationVolaSingleCurveMonteCarloSpot;
        deviations[29] += deviationVolaSingleCurveMonteCarloTerminal * deviationVolaSingleCurveMonteCarloTerminal;

        deviations[30] += deviationValueMultiCurveAnalytic * deviationValueMultiCurveAnalytic;
        deviations[31] += deviationValueMultiCurveMonteCarloSpot * deviationValueMultiCurveMonteCarloSpot;
        deviations[32] += deviationValueMultiCurveMonteCarloTerminal * deviationValueMultiCurveMonteCarloTerminal;

        deviations[33] += deviationValueSingleCurveAnalytic * deviationValueSingleCurveAnalytic;
        deviations[34] += deviationValueSingleCurveMonteCarloSpot * deviationValueSingleCurveMonteCarloSpot;
        deviations[35] += deviationValueSingleCurveMonteCarloTerminal * deviationValueSingleCurveMonteCarloTerminal;
    }

    public static void printProcess(ProcessInterface process, LIBORMarketModelInterface model) throws CalculationException {
        int numberOfComponents = model.getNumberOfComponents();
        int numberOfTimes = process.getTimeDiscretization().getNumberOfTimes();
        if (model instanceof MultiCurveLIBORMarketModel) {
            for (int component = 0; component < numberOfComponents / 2; component++) {
                for (int timeIndex = 0; timeIndex < numberOfTimes &&
                        process.getTime(timeIndex) < model.getLiborPeriodDiscretization().getTime(component); timeIndex++) {
                    RandomVariableInterface forward = process.getProcessValue(timeIndex, component);
                    RandomVariableInterface spread = process.getProcessValue(timeIndex, component + numberOfComponents / 2);
                    try {
                        System.out.println(process.getTime(timeIndex) +
                                "\t" + forward.getAverage() +
                                "\t" + spread.getAverage() +
                                "\t" + ((MultiCurveLIBORMarketModel) model).getProcessValue(timeIndex, component).getAverage() +
                                "\t" + model.getNumeraire(process.getTime(timeIndex + 1)).getAverage() +
                                "\t" + ((MultiCurveLIBORMarketModel) model).getProcessValue(timeIndex, component)
                                .div(model.getNumeraire(process.getTime(timeIndex + 1)))
                                .div(model.getDiscountCurve().getDiscountFactor(process.getTime(timeIndex + 1))).getAverage());
                    } catch (CalculationException e) {
                        e.printStackTrace();
                    }
                }
                System.out.println();
            }
        } else {
            for (int timeIndex = 0; timeIndex < numberOfComponents; timeIndex++) {
                RandomVariableInterface forward = process.getProcessValue(timeIndex, timeIndex);
                try {
                    System.out.println(process.getTime(timeIndex) +
                            "\t" + forward.getAverage() +
                            "\t" + model.getNumeraire(process.getTime(timeIndex + 1)).getAverage() +
                            "\t" + forward
                            .div(model.getNumeraire(process.getTime(timeIndex + 1)))
                            .div(model.getDiscountCurve().getDiscountFactor(process.getTime(timeIndex + 1))).getAverage());
                } catch (CalculationException e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
