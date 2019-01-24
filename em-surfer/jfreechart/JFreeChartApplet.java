import org.jfree.chart.JFreeChart;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.general.DefaultPieDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;


import java.applet.Applet;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.net.URL;

public class JFreeChartApplet extends Applet
{
	// Constants that determine which characters are used to separate points
	private final String POINT_SEPARATOR = ";";
	private final String XY_SEPARATOR = ",";
	// Holds the plotted values that are displayed on the chart
	private XYSeries points;
	// Reference to the actual chart that will be used
	JFreeChart chart;
	
	public void init()
	{
		// draw the points specified in the parameter passed
		String pointsDescription = getParameter("points");
		setPoints(pointsDescription);
	}

	public void stop()
	{
	}
	
	// To allow easy calling from javascript, this function will receive a series of 2D points of
	// the form "1,2;4,5;8,10" representing pairs of points. Each ";" determines the beginning of a new pair
	// and each point's X and Y coordinates are separated by a ","
	public void setPoints(String pointsDescription)
	{
		// Reset the points container
		points = new XYSeries("XYSeries");
		// Separate pairs of points
		String[] pointPairs = pointsDescription.split(POINT_SEPARATOR);
		for(String pointString : pointPairs)
		{
			String[] xyString = pointString.split(XY_SEPARATOR);
			// assuming that there are indeed 2, the first should be X and the second Y
			points.add(Double.parseDouble(xyString[0]), Double.parseDouble(xyString[1]));
		}
		// Add the series to your data set
		XYSeriesCollection dataset = new XYSeriesCollection();
		dataset.addSeries(points);
		// Generate the graph
		chart = ChartFactory.createXYLineChart(null, // Title
														"Zernike Descriptor Number", // x-axis Label
														"Value", // y-axis Label
														dataset, // Dataset
														PlotOrientation.VERTICAL, // Plot Orientation
														false, // Show Legend
														false, // Use tooltips
														false // Configure chart to generate URLs?
														);
		XYPlot chartPlot = (XYPlot)chart.getPlot();
		XYItemRenderer renderer = chartPlot.getRenderer();
		chartPlot.setBackgroundPaint(Color.white);
		renderer.setSeriesPaint(0, Color.black);
		chart.setAntiAlias(true);
	}

	public void paint(Graphics g)
	{
		if ( chart!=null )
		{
			chart.draw( (Graphics2D)g,getBounds()); //repaints the whole chart
		}
	}
 }
