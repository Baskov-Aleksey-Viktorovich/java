package laboratorna_02;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import net.objecthunter.exp4j.Expression;
import net.objecthunter.exp4j.ExpressionBuilder;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

public class FunctionPlotter extends JFrame {

    private JTextField functionField, startField, stopField, stepField;
    private JPanel chartPanel;

    public FunctionPlotter() {
        setTitle("GUIApplication");
        setSize(600, 500);
        setDefaultCloseOperation(EXIT_ON_CLOSE);
        setLayout(new BorderLayout());

        // Верхня панель
        JPanel inputPanel = new JPanel();
        inputPanel.setLayout(new GridLayout(2, 1));

        JPanel funcPanel = new JPanel();
        funcPanel.add(new JLabel("f(x):"));
        functionField = new JTextField("sin(x)/x", 20);
        funcPanel.add(functionField);
        inputPanel.add(funcPanel);

        JPanel rangePanel = new JPanel();
        rangePanel.add(new JLabel("Start:"));
        startField = new JTextField("-6", 5);
        rangePanel.add(startField);
        rangePanel.add(new JLabel("Stop:"));
        stopField = new JTextField("6", 5);
        rangePanel.add(stopField);
        rangePanel.add(new JLabel("Step:"));
        stepField = new JTextField("0.01", 5);
        rangePanel.add(stepField);

        JButton plotButton = new JButton("Plot");
        JButton exitButton = new JButton("Exit");

        rangePanel.add(plotButton);
        rangePanel.add(exitButton);

        inputPanel.add(rangePanel);

        add(inputPanel, BorderLayout.NORTH);

        // Панель графіка
        chartPanel = new JPanel(new BorderLayout());
        add(chartPanel, BorderLayout.CENTER);

        // Обробник кнопки "Plot"
        plotButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                plotFunction();
            }
        });

        // Обробник кнопки "Exit"
        exitButton.addActionListener(e -> System.exit(0));
    }

    private void plotFunction() {
        try {
            String funcText = functionField.getText();
            double start = Double.parseDouble(startField.getText());
            double stop = Double.parseDouble(stopField.getText());
            double step = Double.parseDouble(stepField.getText());

            Expression expression = new ExpressionBuilder(funcText)
                    .variables("x")
                    .build();

            XYSeries seriesFunc = new XYSeries("Function");
            XYSeries seriesDeriv = new XYSeries("Derivative");

            for (double x = start; x <= stop; x += step) {
                expression.setVariable("x", x);
                double y = expression.evaluate();
                seriesFunc.add(x, y);

                // чисельне наближення похідної
                expression.setVariable("x", x + 1e-5);
                double y2 = expression.evaluate();
                double derivative = (y2 - y) / 1e-5;
                seriesDeriv.add(x, derivative);
            }

            XYSeriesCollection dataset = new XYSeriesCollection();
            dataset.addSeries(seriesFunc);
            dataset.addSeries(seriesDeriv);

            JFreeChart chart = ChartFactory.createXYLineChart(
                    "Function and Derivative",
                    "X", "Y",
                    dataset,
                    PlotOrientation.VERTICAL,
                    true, true, false
            );

            chartPanel.removeAll();
            chartPanel.add(new ChartPanel(chart), BorderLayout.CENTER);
            chartPanel.validate();

        } catch (Exception ex) {
            JOptionPane.showMessageDialog(this, "Error: " + ex.getMessage());
        }
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            FunctionPlotter app = new FunctionPlotter();
            app.setVisible(true);
        });
    }
}
