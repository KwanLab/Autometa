package lu.uni.lcsb.vizbin;

import java.awt.geom.Point2D;
import java.util.ArrayList;

// import java.util.List;
import javax.swing.JLabel;
import javax.swing.JPanel;

import lu.uni.lcsb.vizbin.clustering.ClusteringPlugin;

/**
 * 
 * Sample class to show how a Clustering plugin looks like. The getClusters()
 * method simply splits the given points in two clusters.
 * 
 * @author <a href="mailto:valentin.plugaru.001@student.uni.lu">Valentin
 *         Plugaru</a>
 */
public class SamplePlugin implements ClusteringPlugin {

	public ArrayList<ArrayList<Point2D>> getClusters(ArrayList<Point2D> pointList) {
		ArrayList<ArrayList<Point2D>> clusters = new ArrayList<ArrayList<Point2D>>();
		ArrayList<Point2D> cluster = new ArrayList<Point2D>();
		for (int i = 0; i < pointList.size() / 2; i++) {
			cluster.add(pointList.get(i));
		}
		clusters.add(cluster);
		cluster = new ArrayList<Point2D>();
		for (int i = (pointList.size() / 2) + 1; i < pointList.size(); i++) {
			cluster.add(pointList.get(i));
		}
		clusters.add(cluster);
		return clusters;
	}

	public int getPluginModelVersion() {
		return 1;
	}

	@Override
	public void drawOptionsPanel(JPanel optionsPanel) {
		// create interface elements on optionsPanel
		JLabel label = new JLabel("This plugin has no configurable options.");
		optionsPanel.add(label);
	}
}
