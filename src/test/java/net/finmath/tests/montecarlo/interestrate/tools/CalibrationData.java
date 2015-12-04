package net.finmath.tests.montecarlo.interestrate.tools;

import net.finmath.marketdata.model.AnalyticModelInterface;
import net.finmath.montecarlo.interestrate.LIBORMarketModelInterface;
import net.finmath.time.TimeDiscretizationInterface;

public class CalibrationData {
    private TimeDiscretizationInterface timeDiscretization;
    private  TimeDiscretizationInterface periodDiscretization;
    private AnalyticModelInterface analyticModel;
    private  String firstCurveName;
    private  String secondCurveName;
    private  LIBORMarketModelInterface.CalibrationItem[] firstCurveCalibrationItems;
    private  LIBORMarketModelInterface.CalibrationItem[] secondCurveCalibrationItems;
    private  LIBORMarketModelInterface.CalibrationItem[] firstCurveCalibrationItemsMonteCarlo;
    private  LIBORMarketModelInterface.CalibrationItem[] secondCurveCalibrationItemsMonteCarlo;

    public CalibrationData(TimeDiscretizationInterface timeDiscretization,
                           TimeDiscretizationInterface periodDiscretization,
                           AnalyticModelInterface analyticModel,
                           String firstCurveName,
                           String secondCurveName,
                           LIBORMarketModelInterface.CalibrationItem[] firstCurveCalibrationItems,
                           LIBORMarketModelInterface.CalibrationItem[] secondCurveCalibrationItems,
                           LIBORMarketModelInterface.CalibrationItem[] firstCurveCalibrationItemsMonteCarlo,
                           LIBORMarketModelInterface.CalibrationItem[] secondCurveCalibrationItemsMonteCarlo) {
        this.timeDiscretization = timeDiscretization;
        this.periodDiscretization = periodDiscretization;
        this.analyticModel = analyticModel;
        this.firstCurveName = firstCurveName;
        this.secondCurveName = secondCurveName;
        this.firstCurveCalibrationItems = firstCurveCalibrationItems;
        this.secondCurveCalibrationItems = secondCurveCalibrationItems;
        this.firstCurveCalibrationItemsMonteCarlo = firstCurveCalibrationItemsMonteCarlo;
        this.secondCurveCalibrationItemsMonteCarlo = secondCurveCalibrationItemsMonteCarlo;
    }

    public TimeDiscretizationInterface getTimeDiscretization() {
        return timeDiscretization;
    }

    public TimeDiscretizationInterface getPeriodDiscretization() {
        return periodDiscretization;
    }

    public AnalyticModelInterface getAnalyticModel() {
        return analyticModel;
    }

    public String getSecondCurveName() {
        return secondCurveName;
    }

    public String getFirstCurveName() {
        return firstCurveName;
    }

    public LIBORMarketModelInterface.CalibrationItem[] getFirstCurveCalibrationItems() {
        return firstCurveCalibrationItems;
    }

    public LIBORMarketModelInterface.CalibrationItem[] getSecondCurveCalibrationItems() {
        return secondCurveCalibrationItems;
    }

    public LIBORMarketModelInterface.CalibrationItem[] getFirstCurveCalibrationItemsMonteCarlo() {
        return firstCurveCalibrationItemsMonteCarlo;
    }

    public LIBORMarketModelInterface.CalibrationItem[] getSecondCurveCalibrationItemsMonteCarlo() {
        return secondCurveCalibrationItemsMonteCarlo;
    }
}
