package net.finmath.tests.montecarlo.interestrate.tools;

import au.com.bytecode.opencsv.CSVReader;

import java.io.*;

public class TexCreator {
    public static void createPlot(String fileName, int yMin, int yMax, int numberOfPaths,
                                  int numberOfParams, int numberOfFactors,
                                  String multiCurveModel,
                                  boolean useSeperateCorrelationModels) throws IOException {
        String caption = ""
                + numberOfPaths + " paths, "
                + numberOfParams + " params, "
                + numberOfFactors + " factors, "
                + multiCurveModel + " model, "
                + (useSeperateCorrelationModels ? "with " : "without ") + "LDLT-decomp";

        try (FileReader fileReader = new FileReader("templates/plot_template.tex");
             BufferedReader reader = new BufferedReader(fileReader);
             Writer writer = new FileWriter(fileName + ".tex")
        ) {
            reader.lines().forEach(
                    line -> {
                        try {
                            writer.write(line
                                            .replace("fileName", fileName)
                                            .replace("yMin", String.valueOf(yMin))
                                            .replace("yMax", String.valueOf(yMax))
                                            .replace("Caption", caption) + "\n"
                            );
                        } catch (IOException e) {
                            throw new RuntimeException(e);
                        }
                    }
            );
        }
    }

    public static void generateStats() throws IOException {
        for (File params : new File("csv").listFiles()) {
            if (params.getName().contains("_params")) {
                try (FileReader paramsFileReader = new FileReader(params);
                     FileReader fileReader = new FileReader("templates/stats_template.tex");
                     BufferedReader reader = new BufferedReader(fileReader);
                     Writer writer = new FileWriter(params.getAbsolutePath().replace("_params", "").replace("csv", "tex"))
                ) {
                    CSVReader paramsReader = new CSVReader(paramsFileReader);
                    paramsReader.readNext();
                    String calibrationTime = paramsReader.readNext()[0].split("\t")[0];

                    reader.lines().forEach(
                            line -> {
                                try {
                                    writer.write(line
                                                    .replace("calibrationTime", calibrationTime)
                                                    .replace("fileName", params.getName().replace("_params", "_errs")) + "\n"
                                    );
                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }
                            }
                    );
                }
            }
        }
    }

    public static void generatePlots() throws IOException {
        for (File params : new File("csv").listFiles()) {
            if (params.getName().contains("_params")) {
                try (FileReader paramsFileReader = new FileReader(params);
                     FileReader fileReader = new FileReader("templates/plots_template.tex");
                     BufferedReader reader = new BufferedReader(fileReader);
                     Writer writer = new FileWriter(params.getAbsolutePath().replace("_params", "").replace("csv", "tex"))
                ) {
                    CSVReader paramsReader = new CSVReader(paramsFileReader);
                    paramsReader.readNext();
                    String calibrationTime = paramsReader.readNext()[0].split("\t")[0];

                    reader.lines().forEach(
                            line -> {
                                try {
                                    writer.write(line
                                                    .replace("calibrationTime", calibrationTime)
                                                    .replace("fileName", params.getName().replace("_params", "_errs")) + "\n"
                                    );
                                } catch (IOException e) {
                                    throw new RuntimeException(e);
                                }
                            }
                    );
                }
            }
        }
    }

    public static void main(String[] args) throws IOException {
        generateStats();
        generatePlots();
    }
}

