package lu.uni.lcsb.vizbin.clustering;

import java.awt.Polygon;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

import lu.uni.lcsb.vizbin.data.Cluster;
import lu.uni.lcsb.vizbin.data.DataSet;
import lu.uni.lcsb.vizbin.data.Sequence;

import org.apache.log4j.Logger;

/**
 * Factory that creates clusters from dataset..
 * 
 * @author Piotr Gawron
 * 
 */
public class ClusterFactory {
	/**
	 * Default class logger.
	 */
	@SuppressWarnings("unused")
	private final Logger	logger	= Logger.getLogger(ClusterFactory.class);

	/**
	 * Creates a {@link Cluster} from a given dataset and a polygon.
	 * 
	 * @param dataSet
	 *          input dataset
	 * @param polygon
	 *          polygon from which cluster is created
	 * @return a {@link Cluster} from a given dataset and a polygon
	 */
	public Cluster createClusterFromPolygon(DataSet dataSet, List<Point2D> polygon) {
		List<Sequence> elements = new ArrayList<Sequence>();

		Polygon awtPolygon = new Polygon();

		for (Point2D point2d : polygon) {
			awtPolygon.addPoint((int) point2d.getX(), (int) point2d.getY());
		}

		for (Sequence sequence : dataSet.getSequences()) {
			if (awtPolygon.contains(sequence.getLocation())) {
				elements.add(sequence);
			}
		}

		Cluster result = new Cluster(elements);
		return result;
	}
}
