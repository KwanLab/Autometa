package lu.uni.lcsb.vizbin.clustering;

import java.awt.geom.Point2D;
import java.util.ArrayList;
import javax.swing.JPanel;

/**
 * 
 * @author <a href="mailto:valentin.plugaru.001@student.uni.lu">Valentin
 *         Plugaru</a>
 */
public interface ClusteringPlugin {
	/**
	 * Returns list of clusters obtained by the plugin.
	 * 
	 * @param pointList
	 *          list of {@link Point2D points} that should be clustered
	 * @return list of clusters obtained by the plugin
	 */
	ArrayList<ArrayList<Point2D>> getClusters(ArrayList<Point2D> pointList);

	/**
	 * Adds options of the algorithm on the {@link JPanel}.
	 * 
	 * @param optionsPanel
	 *          where the options should be presented
	 */
	void drawOptionsPanel(JPanel optionsPanel);

	/**
	 * Returns version of the plugin.
	 * 
	 * @return version of the plugin
	 */
	int getPluginModelVersion();
}
