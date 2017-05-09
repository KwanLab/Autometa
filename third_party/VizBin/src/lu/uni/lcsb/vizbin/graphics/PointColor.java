package lu.uni.lcsb.vizbin.graphics;

import java.awt.Color;

/**
 * Defines available colors of points for vizualization.
 * 
 * @author Piotr Gawron
 * 
 */
public enum PointColor {
	/**
	 * {@link Color#BLUE}.
	 */
	BLUE(Color.BLUE),
	/**
	 * {@link Color#RED}.
	 */
	RED(Color.RED),
	/**
	 * {@link Color#GREEN}.
	 */
	GREEN(Color.GREEN),
	/**
	 * {@link Color#ORANGE}.
	 */
	ORANGE(Color.ORANGE),
	/**
	 * {@link Color#BLACK}.
	 */
	BLACK(Color.BLACK);

	/**
	 * Color definition.
	 */
	private Color	color;

	/**
	 * Default constructor.
	 * 
	 * @param color
	 *          {@link #color}
	 */
	private PointColor(Color color) {
		this.color = color;
	}

	/**
	 * @return the color
	 * @see #color
	 */
	public Color getColor() {
		return color;
	}
}
