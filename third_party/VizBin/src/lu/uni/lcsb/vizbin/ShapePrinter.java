package lu.uni.lcsb.vizbin;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Paint;
import java.awt.Stroke;

import lu.uni.lcsb.vizbin.graphics.PointColor;
import lu.uni.lcsb.vizbin.graphics.PointShape;
import lu.uni.lcsb.vizbin.graphics.StarPolygon;

public class ShapePrinter {

	/**
	 * Start angle of the star shape.
	 */
	private static final int		START_ANGLE	= 60;

	/**
	 * How many points should have a star.
	 */
	private static final int		STAR_POINTS	= 5;
	/**
	 * Bold stroke.
	 */
	private static final Stroke	BOLD_STROKE	= new BasicStroke(1.5f);

	private Graphics2D					g2;
	private Stroke							defaultStroke;

	public ShapePrinter(Graphics2D graphics) {
		g2 = graphics;
		defaultStroke = g2.getStroke();
	}

	protected void drawShape(int size, int transX, int transY, Paint color, PointShape shape) {
		g2.setPaint(color);
		int[] xCoords;
		int[] yCoords;
		switch (shape) {
			case RECTANGLE:
				g2.fillRect(transX, transY, size, size);
				break;
			case CIRCLE:
				g2.fillOval(transX, transY, size, size);
				break;
			case TRIANGLE_DOWN:
				xCoords = new int[] { transX - size, transX, transX + size };
				yCoords = new int[] { transY, transY - size, transY };
				g2.fillPolygon(xCoords, yCoords, xCoords.length);
				break;
			case TRIANGLE_UP:
				xCoords = new int[] { transX - size, transX, transX + size };
				yCoords = new int[] { transY, transY + size, transY };
				g2.fillPolygon(xCoords, yCoords, xCoords.length);
				break;
			case CROSS:
				g2.setStroke(BOLD_STROKE); // Make the stroke of the
																		// polygon a bit thicker
				g2.drawLine((int) transX - (size / 2), (int) transY - (size / 2), (int) transX + (size / 2), (int) transY + (size / 2));
				g2.drawLine((int) transX + (size / 2), (int) transY - (size / 2), (int) transX - (size / 2), (int) transY + (size / 2));
				g2.setStroke(defaultStroke); // Reset the stroke;not sure if
																			// needed;
				break;
			case PLUS:
				g2.setStroke(BOLD_STROKE); // Make the stroke of the
																		// polygon a bit thicker
				g2.drawLine((int) transX - (size / 2), (int) transY, (int) transX + (size / 2), (int) transY);
				g2.drawLine((int) transX, (int) transY - (size / 2), (int) transX, (int) transY + (size / 2));
				g2.setStroke(defaultStroke); // Reset the stroke;not sure if
				// needed;
				break;
			case STAR:
				StarPolygon star = new StarPolygon(transX, transY, size + 1, size / 2, STAR_POINTS, START_ANGLE);
				g2.fill(star);
				// The star shape is a particular shape and to discriminate it
				// from similarly colored "normal" points it gets a black boundary
				g2.setPaint(Color.BLACK);
				g2.setStroke(BOLD_STROKE); // Make the stroke of the
				// polygon a bit thicker
				g2.draw(star);
				g2.setStroke(defaultStroke); // Reset the stroke;not sure if
				// needed

				break;
			default:
				throw new InvalidArgumentException("Unknown shape: " + shape);
		}
	}

	protected void drawShape(int size, int transX, int transY, int labelId) {
		Color color = getColor(labelId);
		PointShape shape = getShape(labelId);
		drawShape(size, transX, transY, color, shape);
	}

	public static Color getColor(int id) {
		int size = PointColor.values().length;
		int base = id % size;
		return PointColor.values()[base].getColor();
	}

	public static PointShape getShape(int id) {
		int base = id % (PointShape.values().length - 1);
		return PointShape.values()[base];
	}

}
