package lu.uni.lcsb.vizbin.data;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

/**
 * Class defining single input data set.
 * 
 * @author Piotr Gawron
 * 
 */
public class DataSet {
	/**
	 * Default class logger.
	 */
	@SuppressWarnings("unused")
	private final Logger		logger						= Logger.getLogger(DataSet.class);

	/**
	 * List of sequences.
	 */
	private List<Sequence>	sequences					= new ArrayList<Sequence>();

	/**
	 * Minimum sequence length.
	 * 
	 * @see Sequence#length
	 */
	private double					minSequenceLength	= 0;

	/**
	 * Adds sequence to dataset.
	 * 
	 * @param sequence
	 *          sequence to add
	 */
	public void addSequence(Sequence sequence) {
		sequences.add(sequence);
		if (sequence.getLength() != null) {
			// First sequence that is visited? -> Initialize to a meaningful value
			if ((double) 0 == this.minSequenceLength) {
				this.minSequenceLength = sequence.getLength();
			} else {
				this.minSequenceLength = Math.min(this.minSequenceLength, sequence.getLength());
			}
		}
	}

	/**
	 * Computes and return number of different labels.
	 * 
	 * @return number of different labels
	 */
	public Integer getLabelsCount() {
		Set<Integer> labels = new HashSet<Integer>();
		for (Sequence s : sequences) {
			labels.add(s.getLabelId());
		}
		return labels.size();
	}

	/**
	 * Returns number of sequences in a dataset.
	 * 
	 * @return number of sequences in a dataset
	 */
	public int getSize() {
		return sequences.size();
	}

	/**
	 * @return the maxSequenceLength
	 * @see #minSequenceLength
	 */
	public double getMinSequenceLength() {
		return minSequenceLength;
	}

	/**
	 * @param minSequenceLength
	 *          the minSequenceLength to set
	 * @see #minSequenceLength
	 */
	public void setMinSequenceLength(double minSequenceLength) {
		this.minSequenceLength = minSequenceLength;
	}

	/**
	 * @return the sequences
	 * @see #sequences
	 */
	public List<Sequence> getSequences() {
		return sequences;
	}

	/**
	 * @param sequences
	 *          the sequences to set
	 * @see #sequences
	 */
	public void setSequences(List<Sequence> sequences) {
		this.sequences = sequences;
	}
}
