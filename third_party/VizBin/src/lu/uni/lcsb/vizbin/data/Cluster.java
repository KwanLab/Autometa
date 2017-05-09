package lu.uni.lcsb.vizbin.data;

import java.util.ArrayList;
import java.util.List;

/**
 * Defaines a cluster of {@link Sequence sequences}.
 * 
 * @author Piotr Gawron
 * 
 */
public class Cluster {
	/**
	 * List of sequences.
	 */
	private List<Sequence>	elements	= new ArrayList<Sequence>();

	/**
	 * Default constructor.
	 * 
	 * @param elements
	 *          {@link #elements sequences} to add
	 */
	public Cluster(List<Sequence> elements) {
		this.elements = elements;
	}

	/**
	 * @return the elements
	 * @see #elements
	 */
	public List<Sequence> getElements() {
		return elements;
	}

	/**
	 * @param elements
	 *          the elements to set
	 * @see #elements
	 */
	public void setElements(List<Sequence> elements) {
		this.elements = elements;
	}

}
