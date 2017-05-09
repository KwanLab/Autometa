package lu.uni.lcsb.vizbin;

// Structure used in DatasSetUtils class to hold return values of createKmers
// function
public class KmersResult {
	private Integer[]	result;
	private int				numKmersRemoved;

	public KmersResult(Integer[] result, int numKmersRemoved) {
		this.result = result;
		this.numKmersRemoved = numKmersRemoved;
	}

	public Integer[] getResult() {
		return result;
	}

	public int getNumKmersRemoved() {
		return numKmersRemoved;
	}
}
