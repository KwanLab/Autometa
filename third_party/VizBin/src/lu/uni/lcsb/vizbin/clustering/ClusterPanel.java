package lu.uni.lcsb.vizbin.clustering;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.SwingUtilities;

import lu.uni.lcsb.vizbin.DataExporter;
import lu.uni.lcsb.vizbin.ExtendedFastScatterPlot;
import lu.uni.lcsb.vizbin.LegendFrame;
import lu.uni.lcsb.vizbin.ShapePrinter;
import lu.uni.lcsb.vizbin.data.Cluster;
import lu.uni.lcsb.vizbin.data.DataSet;
import lu.uni.lcsb.vizbin.data.Sequence;
import lu.uni.lcsb.vizbin.graphics.PointShape;

import org.apache.log4j.Logger;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.ChartRenderingInfo;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.ValueAxis;
import org.jfree.chart.plot.FastScatterPlot;
import org.jfree.ui.RectangleEdge;

/**
 * This class creates and manages {@link ChartPanel} that has informations about
 * clusters (in fact it's responsible for the visualization of the results of
 * the system).
 * 
 * @author Piotr Gawron
 * 
 */
public class ClusterPanel {
	/**
	 * Maximum alpha value.
	 */
	private static final int	MAX_ALPHA_VALUE						= 255;
	/**
	 * What size should be used when information about sequence length is missing.
	 */
	private static final int	DEFAULT_MIN_SEQUENCE_SIZE	= 4;
	/**
	 * Default class logger.
	 */
	private Logger						logger										= Logger.getLogger(ClusterPanel.class);
	/**
	 * Panel where the data is visualized.
	 */
	private ChartPanel				panel;
	/**
	 * Currently selected polygon on the {@link #panel}.
	 */
	private List<Point2D>			polygon										= new ArrayList<Point2D>();

	/**
	 * Dataset that is visualized.
	 */
	private DataSet						dataSet;

	/**
	 * Parent {@link JFrame} where the {@link #panel} is located.
	 */
	private JFrame						frame;

	/**
	 * Name of the file that is visualized.
	 */
	private String						filename;

	/**
	 * Frame with the legend about visualized data.
	 */
	private LegendFrame				legendFrame;

	/**
	 * Class that extends {@link MouseAdapter} and will allow to select polygon on
	 * the {@link MouseMarker#chart chart}.
	 * 
	 * @author Piotr Gawron
	 * 
	 */
	private final class MouseMarker extends MouseAdapter {
		/**
		 * Chart where we edit polygon.
		 */
		private final JFreeChart	chart;
		/**
		 * Panel with the {@link #chart}.
		 */
		private final ChartPanel	panel;
		/**
		 * Structure where the polygon is stored.
		 */
		private List<Point2D>			polygon;

		/**
		 * Default constructor.
		 * 
		 * @param panel
		 *          {@link #panel}
		 * @param polygon
		 *          {@link #polygon}
		 */
		public MouseMarker(ChartPanel panel, List<Point2D> polygon) {
			this.panel = panel;
			this.chart = panel.getChart();
			this.polygon = polygon;
			// this.plot.setDomainGridlinesVisible(false);
			// this.plot.setRangeGridlinesVisible(false);

		}

		@Override
		public void mouseClicked(MouseEvent e) {
			// Put edges of the polygon for selecting a set of points in 2D
			if (SwingUtilities.isLeftMouseButton(e)) {
				// Motivated by http://www.jfree.org/phpBB2/viewtopic.php?p=54140
				int mouseX = e.getX();
				int mouseY = e.getY();
				Point2D p = panel.translateScreenToJava2D(new Point(mouseX, mouseY));
				FastScatterPlot plot = (FastScatterPlot) chart.getPlot();
				ChartRenderingInfo info = panel.getChartRenderingInfo();
				Rectangle2D dataArea = info.getPlotInfo().getDataArea();

				ValueAxis domainAxis = plot.getDomainAxis();
				// RectangleEdge domainAxisEdge = plot.getDomainAxisEdge();
				ValueAxis rangeAxis = plot.getRangeAxis();
				// RectangleEdge rangeAxisEdge = plot.getRangeAxisEdge();
				double chartX = domainAxis.java2DToValue(p.getX(), dataArea, RectangleEdge.BOTTOM);
				double chartY = rangeAxis.java2DToValue(p.getY(), dataArea, RectangleEdge.LEFT);
				// DEBUG
				polygon.add(new Point2D.Double(chartX, chartY));

			}

		}
	}

	/**
	 * Default constructor.
	 * 
	 * @param ds
	 *          {@link ClusterPanel#dataSet}
	 * @param inFileName
	 *          {@link ClusterPanel#filename}
	 * @param parentFrame
	 *          {@link ClusterPanel#frame}
	 */
	public ClusterPanel(DataSet ds, String inFileName, JFrame parentFrame) {

		this.dataSet = ds;
		this.frame = parentFrame;
		this.filename = inFileName;
		ArrayList<Sequence> pointList = (ArrayList<Sequence>) ds.getSequences();
		int count = pointList.size();

		int minSize = DEFAULT_MIN_SEQUENCE_SIZE;
		// If there is additional length information provided use a smaller default
		// minimum size
		if ((double) 0 < ds.getMinSequenceLength()) {
			minSize = 2;
		}
		float[][] data = new float[2][count];
		int[] sizes = new int[count];
		Color[] colors = new Color[count];
		PointShape[] shapes = new PointShape[count];

		int counter = 0;
		float alpha = 0f;
		legendFrame = new LegendFrame();

		for (Sequence sequence : pointList) {
			data[0][counter] = (float) sequence.getLocation().getX();
			data[1][counter] = (float) sequence.getLocation().getY();
			// DEBUG
			// logger.debug(sequence.getLocation());
			colors[counter] = ShapePrinter.getColor(sequence.getLabelId());
			shapes[counter] = ShapePrinter.getShape(sequence.getLabelId());

			legendFrame.setLabel(sequence.getLabelId(), sequence.getLabelName());

			if (sequence.getCoverage() != null) {
				alpha = (sequence.getCoverage()).floatValue();
				alpha = alpha * alpha;
				colors[counter] = new Color(
						(int) colors[counter].getRed(), (int) colors[counter].getGreen(), (int) colors[counter].getBlue(), (int) (alpha * MAX_ALPHA_VALUE));
			} else // Either coverage OR GC content
			if (sequence.getGc() != null) {
				alpha = (sequence.getGc()).floatValue();
				alpha = alpha * alpha;
				colors[counter] = new Color(
						(int) colors[counter].getRed(), (int) colors[counter].getGreen(), (int) colors[counter].getBlue(), (int) (alpha * MAX_ALPHA_VALUE));
			}

			if (sequence.getMarker() != null && sequence.getMarker()) {
				shapes[counter] = PointShape.STAR;
			}
			// Use a minimum size
			sizes[counter] = minSize;
			// Adjust the size if additional length information is given
			if (sequence.getLength() != null) {
				double magnification = (sequence.getLength()).doubleValue() - ds.getMinSequenceLength();
				// magnification = magnification * magnification;
				sizes[counter] = (int) (sizes[counter] + (sizes[counter] * magnification));

			}
			counter++;
		}

		final NumberAxis domainAxis = new NumberAxis("");
		domainAxis.setAutoRangeIncludesZero(false);
		domainAxis.setTickMarksVisible(false);
		domainAxis.setTickLabelsVisible(false);

		final NumberAxis rangeAxis = new NumberAxis("");
		rangeAxis.setAutoRangeIncludesZero(false);
		rangeAxis.setTickMarksVisible(false);
		rangeAxis.setTickLabelsVisible(false);

		ExtendedFastScatterPlot plot = new ExtendedFastScatterPlot(data, domainAxis, rangeAxis, sizes, colors, shapes, polygon);
		plot.setDomainGridlinesVisible(false);
		plot.setRangeGridlinesVisible(false);
		plot.setDomainPannable(true);
		plot.setRangePannable(true);
		final JFreeChart chart = new JFreeChart("", plot);
		panel = new ChartPanel(chart);
		panel.addMouseListener(new MouseMarker(panel, polygon));
		panel.setMouseWheelEnabled(true); // Enable scroll-wheel support for zooming

		JPopupMenu popup = panel.getPopupMenu();

		JMenu selectionMenu = new JMenu("Selection");
		JMenuItem exportMenu = new JMenuItem("Export");
		exportMenu.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				Cluster selectedCluster = new ClusterFactory().createClusterFromPolygon(dataSet, polygon);
				if (selectedCluster.getElements().size() > 0) {
					int option;
					option = JOptionPane
							.showConfirmDialog(
									null, selectedCluster.getElements().size() + " sequences selected, export to file? \n\nPress Cancel to remove your selection.",
									"Export selected cluster", JOptionPane.YES_NO_CANCEL_OPTION);
					if (option == JOptionPane.YES_OPTION) {
						File outFile = DataExporter.exportCluster(frame, filename, selectedCluster.getElements());

						try {
							String png = outFile.getCanonicalPath() + ".png";
							int width = panel.getWidth();
							int height = panel.getHeight();
							BufferedImage image = new BufferedImage((int) width, (int) height, BufferedImage.TYPE_INT_ARGB);
							Graphics2D g = image.createGraphics();
							panel.paint(g);
							panel.printAll(g);
							File f = new File(png);
							// png is an image format (like gif or jpg)
							ImageIO.write(image, "png", f);
						} catch (IOException ex) {
							logger.error(ex, ex);
						}
					} else if (option == JOptionPane.CANCEL_OPTION) {
						polygon.clear();

						panel.repaint();
					}
				}
			}
		});
		JMenuItem clearMenu = new JMenuItem("Clear selection");
		clearMenu.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				polygon.clear();
				chart.fireChartChanged();
			}
		});

		selectionMenu.add(exportMenu);
		selectionMenu.add(clearMenu);

		popup.add(selectionMenu);

		JMenuItem legendMenu = new JMenuItem("Legend");

		legendMenu.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				LegendFrame legendFrame = getLegendFrame();
				legendFrame.setVisible(!legendFrame.isVisible());
			}
		});
		popup.add(legendMenu);

	}

	/**
	 * Returns {@link #panel}.
	 * 
	 * @return {@link #panel}
	 */
	public ChartPanel getChartPanel() {
		return panel;
	}

	/**
	 * @return the legendFrame
	 * @see #legendFrame
	 */
	public LegendFrame getLegendFrame() {
		return legendFrame;
	}

	/**
	 * @param legendFrame
	 *          the legendFrame to set
	 * @see #legendFrame
	 */
	public void setLegendFrame(LegendFrame legendFrame) {
		this.legendFrame = legendFrame;
	}
}
