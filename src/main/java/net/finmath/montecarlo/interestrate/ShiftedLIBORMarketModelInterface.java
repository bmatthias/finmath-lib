package net.finmath.montecarlo.interestrate;

public interface ShiftedLIBORMarketModelInterface extends LIBORMarketModelInterface {

    double getShiftParameter(int component);
    double getLIBORShift(int liborIndex);
}