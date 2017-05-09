package lu.uni.lcsb.vizbin;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JProgressBar;
import javax.swing.JTabbedPane;

/**
 * This immutable class defines all gui object parameters required for
 * computations.
 * 
 * @author Piotr Gawron
 * 
 */
public class ProcessGuiParameters {
	/**
	 * {@link JLabel} where the status is displayed.
	 */
	private JLabel				status;
	/**
	 * {@link JProgressBar} where the progress is displayed.
	 */
	private JProgressBar	progBar;
	/**
	 * {@link JTabbedPane} where the chart will be displayed.
	 */
	private JTabbedPane		tabPane;
	/**
	 * Frame of the program.
	 */
	private JFrame				parentFrame;
	/**
	 * TODO not used.
	 */
	private boolean				drawAxes;

	/**
	 * Default constructor.
	 * 
	 * @param status
	 *          {@link #status}
	 * @param progBar
	 *          {@link #progBar}
	 * @param tabPane
	 *          {@link #tabPane}
	 * @param parentFrame
	 *          {@link #parentFrame}
	 * @param drawAxes
	 *          {@link #drawAxes}
	 */
	public ProcessGuiParameters(JLabel status, JProgressBar progBar, JTabbedPane tabPane, JFrame parentFrame, boolean drawAxes) {
		this.status = status;
		this.progBar = progBar;
		this.tabPane = tabPane;
		this.parentFrame = parentFrame;
		this.drawAxes = drawAxes;
	}

	/**
	 * @return the status
	 * @see #status
	 */
	public JLabel getStatusLabel() {
		return status;
	}

	/**
	 * @return the progBar
	 * @see #progBar
	 */
	public JProgressBar getProgessBar() {
		return progBar;
	}

	/**
	 * @return the tabPane
	 * @see #tabPane
	 */
	public JTabbedPane getTabPane() {
		return tabPane;
	}

	/**
	 * @return the parentFrame
	 * @see #parentFrame
	 */
	public JFrame getParentFrame() {
		return parentFrame;
	}

	/**
	 * @return the drawAxes
	 * @see #drawAxes
	 */
	public boolean isDrawAxes() {
		return drawAxes;
	}

	/**
	 * @param drawAxes
	 *          the drawAxes to set
	 * @see #drawAxes
	 */
	public void setDrawAxes(boolean drawAxes) {
		this.drawAxes = drawAxes;
	}

}
