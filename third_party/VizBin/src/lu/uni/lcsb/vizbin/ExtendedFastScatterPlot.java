package lu.uni.lcsb.vizbin;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.List;

import lu.uni.lcsb.vizbin.graphics.PointShape;

import org.apache.log4j.Logger;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.CrosshairState;
import org.jfree.chart.plot.FastScatterPlot;
import org.jfree.chart.plot.PlotRenderingInfo;
import org.jfree.ui.RectangleEdge;

public class ExtendedFastScatterPlot extends FastScatterPlot {
	/**
	 * Default class logger.
	 */
	@SuppressWarnings("unused")
	private final Logger										logger						= Logger.getLogger(ExtendedFastScatterPlot.class);

	/**
    *
    */
	private static final long	serialVersionUID	= 1L;

	
	private final static int	POLYGON_EDGE_SIZE	= 6;

	private int[]											sizes;
	private Paint[]										colors;
	private PointShape[]							shapes;
	// Double[] alpha;
	private List<Point2D>							polygon;

	public ExtendedFastScatterPlot(float[][] data, NumberAxis domainAxis, NumberAxis rangeAxis, int[] sizes, Paint[] colors, PointShape[] shapes,
			List<Point2D> polygon) {
		super(data, domainAxis, rangeAxis);
		this.sizes = sizes;
		this.colors = colors;
		this.shapes = shapes;
		this.polygon = polygon;

	}

	@Override
	public void render(Graphics2D g2, Rectangle2D dataArea, PlotRenderingInfo info, CrosshairState crosshairState) {
		ShapePrinter shapePrinter = new ShapePrinter(g2);
		if (this.getData() != null) {
			for (int i = 0; i < this.getData()[0].length; i++) {
				float x = this.getData()[0][i];
				float y = this.getData()[1][i];
				int size = this.sizes[i];
				int transX = (int) this.getDomainAxis().valueToJava2D(x, dataArea, RectangleEdge.BOTTOM);
				int transY = (int) this.getRangeAxis().valueToJava2D(y, dataArea, RectangleEdge.LEFT);

				Paint color = this.colors[i];
				PointShape shape = this.shapes[i];
				
				shapePrinter.drawShape(size, transX, transY, color, shape);
			}
		}
		g2.setPaint(Color.RED); // Set the color of the polygon
		g2.setStroke(new BasicStroke(2)); // Make the stroke of the polygon a bit
																			// thicker
		for (int i = 0; i < polygon.size() - 1; i++) {
			int x1 = (int) this.getDomainAxis().valueToJava2D(polygon.get(i).getX(), dataArea, RectangleEdge.BOTTOM);
			int y1 = (int) this.getRangeAxis().valueToJava2D(polygon.get(i).getY(), dataArea, RectangleEdge.LEFT);

			int x2 = (int) this.getDomainAxis().valueToJava2D(polygon.get(i + 1).getX(), dataArea, RectangleEdge.BOTTOM);
			int y2 = (int) this.getRangeAxis().valueToJava2D(polygon.get(i + 1).getY(), dataArea, RectangleEdge.LEFT);

			g2.drawLine(x1, y1, x2, y2);

			// Move the polygon edges just a slight bit such that they are nicely
			// centered with the edges of the line.
			x1 = x1 - (POLYGON_EDGE_SIZE / 2);
			y1 = y1 - (POLYGON_EDGE_SIZE / 2);
			x2 = x2 - (POLYGON_EDGE_SIZE / 2);
			y2 = y2 - (POLYGON_EDGE_SIZE / 2);

			g2.fillRect(x1, y1, POLYGON_EDGE_SIZE, POLYGON_EDGE_SIZE);
			g2.fillRect(x2, y2, POLYGON_EDGE_SIZE, POLYGON_EDGE_SIZE);
		}
		if (polygon.size() > 1) {
			int x1 = (int) this.getDomainAxis().valueToJava2D(polygon.get(0).getX(), dataArea, RectangleEdge.BOTTOM);
			int y1 = (int) this.getRangeAxis().valueToJava2D(polygon.get(0).getY(), dataArea, RectangleEdge.LEFT);

			int x2 = (int) this.getDomainAxis().valueToJava2D(polygon.get(polygon.size() - 1).getX(), dataArea, RectangleEdge.BOTTOM);
			int y2 = (int) this.getRangeAxis().valueToJava2D(polygon.get(polygon.size() - 1).getY(), dataArea, RectangleEdge.LEFT);

			g2.drawLine(x1, y1, x2, y2);

			// Move the polygon edges just a slight bit such that they are nicely
			// centered with the edges of the line.
			x1 = x1 - (POLYGON_EDGE_SIZE / 2);
			y1 = y1 - (POLYGON_EDGE_SIZE / 2);
			x2 = x2 - (POLYGON_EDGE_SIZE / 2);
			y2 = y2 - (POLYGON_EDGE_SIZE / 2);

			g2.fillRect(x1, y1, POLYGON_EDGE_SIZE, POLYGON_EDGE_SIZE);
			g2.fillRect(x2, y2, POLYGON_EDGE_SIZE, POLYGON_EDGE_SIZE);
		}
		if (polygon.size() == 1) { // Draw a filled rectancle based on the clicked
																// coordinates
			int x1 = (int) this.getDomainAxis().valueToJava2D(polygon.get(0).getX(), dataArea, RectangleEdge.BOTTOM);
			int y1 = (int) this.getRangeAxis().valueToJava2D(polygon.get(0).getY(), dataArea, RectangleEdge.LEFT);

			// Move the polygon edges just a slight bit such that they are nicely
			// centered with the edges of the line.
			x1 = x1 - (POLYGON_EDGE_SIZE / 2);
			y1 = y1 - (POLYGON_EDGE_SIZE / 2);

			g2.fillRect(x1, y1, POLYGON_EDGE_SIZE, POLYGON_EDGE_SIZE);
		}

	}

}
