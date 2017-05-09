package lu.uni.lcsb.vizbin;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.util.HashMap;
import java.util.Map;

import javax.swing.JFrame;
import javax.swing.JPanel;

import lu.uni.lcsb.vizbin.graphics.PointShape;

public class LegendFrame extends javax.swing.JFrame {

	private Map<Integer, String>	labels		= new HashMap<Integer, String>();
	private boolean								star			= false;
	private String								starLabel	= "STAR :)";

	class DrawPane extends JPanel {
		private JFrame	frame;

		private DrawPane(JFrame frame) {
			this.frame = frame;
		}

		/**
		 * 
		 */
		private static final long	serialVersionUID	= 1L;

		public void paintComponent(Graphics g) {
			ShapePrinter shapePrinter = new ShapePrinter((Graphics2D) g);

			int width = 200;
			int y = 20;
			for (Integer integer : labels.keySet()) {
				y += 20;
				width = Math.max(width, g.getFontMetrics().stringWidth(labels.get(integer)));
			}
			if (isStar()) {
				width = Math.max(width, g.getFontMetrics().stringWidth(labels.get(getStarLabel())));
			}
			g.setColor(Color.WHITE);
			g.fillRect(0, 0, 100 + width, y + 50);
			frame.setSize(100 + width, y + 50);
			g.setColor(Color.BLACK);
			y = 20;
			for (Integer integer : labels.keySet()) {
				shapePrinter.drawShape(5, 20, y, integer);
				g.drawString(labels.get(integer), 30, y + 5);
				y += 20;
				width = Math.max(width, g.getFontMetrics().stringWidth(labels.get(integer)));
			}
			if (isStar()) {
				shapePrinter.drawShape(5, 20, y, ShapePrinter.getColor(0), PointShape.STAR);
				g.drawString(getStarLabel(), 30, y + 5);
			}

		}
	}

	/**
	 * 
	 */
	private static final long	serialVersionUID	= 1L;

	public LegendFrame() {
		super.setLocation(100, 100);
		super.setSize(320, 320);
		setContentPane(new DrawPane(this));
		setResizable(false);
		setTitle("Legend");
	}

	public void setLabel(Integer labelId, String labelName) {
		labels.put(labelId, labelName);
	}

	/**
	 * @return the star
	 * @see #star
	 */
	public boolean isStar() {
		return star;
	}

	/**
	 * @param star
	 *          the star to set
	 * @see #star
	 */
	public void setStar(boolean star) {
		this.star = star;
	}

	/**
	 * @return the starLabel
	 * @see #starLabel
	 */
	public String getStarLabel() {
		return starLabel;
	}

	/**
	 * @param starLabel the starLabel to set
	 * @see #starLabel
	 */
	public void setStarLabel(String starLabel) {
		this.starLabel = starLabel;
	}

}
