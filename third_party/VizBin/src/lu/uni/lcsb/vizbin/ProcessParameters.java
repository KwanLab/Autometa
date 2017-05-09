package lu.uni.lcsb.vizbin;

import lu.uni.lcsb.vizbin.pca.PcaType;

/**
 * This immutable class defines parameters for computations.
 * 
 * @author Piotr Gawron
 * 
 */
public class ProcessParameters {
	/**
	 * Input fasta file.
	 */
	private String	inputFastaFile;
	/**
	 * Input label file.
	 */
	private String	inputLabelFile;

	/**
	 * Input point file.
	 */
	private String	inputPointFile;
	/**
	 * Output point file (if defined the data will be saved here, not in the
	 * temporary file).
	 */
	private String	outputFile;

	/**
	 * Output k-mer debug file.
	 */
	private String	kmerDebugFile;

	/**
	 * Length of the k-mer.
	 */
	private Integer	kMerLength				= Config.DEFAULT_KMER_LENGTH;
	/**
	 * Minimum contig length.
	 */
	private Integer	contigLength			= Config.DEFAULT_CONTIG_LENGTH;
	/**
	 * Number of threads to use.
	 */
	private Integer	threads						= Config.DEFAULT_THREAD_NUM;

	/**
	 * Pca columns.
	 */
	private Integer	pcaColumns				= Config.DEFAULT_PCA_COLUMNS;

	/**
	 * Random seed.
	 */
	private Integer	seed							= Config.DEFAULT_SEED;

	/**
	 * Pca algorithm {@link PcaType type}.
	 */
	private PcaType	pcaAlgorithmType	= Config.DEFAULT_PCA_TYPE;

	/**
	 * Theta value.
	 */
	private Double	theta							= Config.DEFAULT_THETA;

	/**
	 * Perplexity.
	 */
	private Double	perplexity				= Config.DEFAULT_PERPLEXILITY;

	/**
	 * Should the k-mers and reverse k-mers be merged.
	 */
	private Boolean	merge							= Config.DEFAULT_MERGE;

	/**
	 * Extra information in logs.
	 */
	private Boolean	extendedLogs			= false;

	/**
	 * Default constructor that creates a copy of the {@link ProcessParameters}
	 * object.
	 * 
	 * @param original
	 *          original parameter object
	 */
	public ProcessParameters(ProcessParameters original) {
		this.inputFastaFile = original.inputFastaFile;
		this.inputPointFile = original.inputPointFile;
		this.inputLabelFile = original.inputLabelFile;
		this.outputFile = original.outputFile;
		this.kmerDebugFile = original.kmerDebugFile;
		this.kMerLength = original.kMerLength;
		this.contigLength = original.contigLength;
		this.threads = original.threads;
		this.pcaColumns = original.pcaColumns;
		this.seed = original.seed;
		this.pcaAlgorithmType = original.pcaAlgorithmType;
		this.theta = original.theta;
		this.perplexity = original.perplexity;
		this.merge = original.merge;
		this.extendedLogs = original.extendedLogs;
	}

	/**
	 * Default constructor.
	 */
	public ProcessParameters() {

	}

	/**
	 * @return the inputFastaFile
	 * @see #inputFastaFile
	 */
	public String getInputFastaFile() {
		return inputFastaFile;
	}

	/**
	 * @param inputFastaFile
	 *          the inputFastaFile to set
	 * @see #inputFastaFile
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters inputFastaFile(String inputFastaFile) {
		ProcessParameters result = new ProcessParameters(this);
		if ("".equals(inputFastaFile)) {
			result.inputFastaFile = null;
		} else {
			result.inputFastaFile = inputFastaFile;
		}
		return result;
	}

	/**
	 * @return the inputPointFile
	 * @see #inputPointFile
	 */
	public String getInputPointFile() {
		return inputPointFile;
	}

	/**
	 * @param inputPointFile
	 *          the inputPointFile to set
	 * @see #inputPointFile
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters inputPointFile(String inputPointFile) {
		ProcessParameters result = new ProcessParameters(this);
		if ("".equals(inputPointFile)) {
			result.inputPointFile = null;
		} else {
			result.inputPointFile = inputPointFile;
		}
		return result;
	}

	/**
	 * @return the outputFile
	 * @see #outputFile
	 */
	public String getOutputFile() {
		return outputFile;
	}

	/**
	 * @param outputFile
	 *          the outputFile to set
	 * @see #outputFile
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters outputFile(String outputFile) {
		ProcessParameters result = new ProcessParameters(this);
		if ("".equals(outputFile)) {
			result.outputFile = null;
		} else {
			result.outputFile = outputFile;
		}
		return result;
	}

	/**
	 * @return the kMerLength
	 * @see #kMerLength
	 */
	public Integer getkMerLength() {
		return kMerLength;
	}

	/**
	 * @param kMerLength
	 *          the kMerLength to set
	 * @see #kMerLength
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters kMerLength(Integer kMerLength) {
		ProcessParameters result = new ProcessParameters(this);
		result.kMerLength = kMerLength;
		return result;
	}

	/**
	 * @return the contigLength
	 * @see #contigLength
	 */
	public Integer getContigLength() {
		return contigLength;
	}

	/**
	 * @param contigLength
	 *          the contigLength to set
	 * @see #contigLength
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters contigLength(Integer contigLength) {
		ProcessParameters result = new ProcessParameters(this);
		result.contigLength = contigLength;
		return result;
	}

	/**
	 * @return the threads
	 * @see #threads
	 */
	public Integer getThreads() {
		return threads;
	}

	/**
	 * @param threads
	 *          the threads to set
	 * @see #threads
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters threads(Integer threads) {
		ProcessParameters result = new ProcessParameters(this);
		result.threads = threads;
		return result;
	}

	/**
	 * @return the inputLabelFile
	 * @see #inputLabelFile
	 */
	public String getInputLabelFile() {
		return inputLabelFile;
	}

	/**
	 * @param inputLabelFile
	 *          the inputLabelFile to set
	 * @see #inputLabelFile
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters inputLabelFile(String inputLabelFile) {
		ProcessParameters result = new ProcessParameters(this);
		if ("".equals(inputLabelFile)) {
			result.inputLabelFile = null;
		} else {
			result.inputLabelFile = inputLabelFile;
		}
		return result;
	}

	/**
	 * @return the kmerDebugFile
	 * @see #kmerDebugFile
	 */
	public String getKmerDebugFile() {
		return kmerDebugFile;
	}

	/**
	 * @param kmerDebugFile
	 *          the kmerDebugFile to set
	 * @see #kmerDebugFile
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters kmerDebugFile(String kmerDebugFile) {
		ProcessParameters result = new ProcessParameters(this);
		if ("".equals(kmerDebugFile)) {
			result.kmerDebugFile = null;
		} else {
			result.kmerDebugFile = kmerDebugFile;
		}
		return result;
	}

	/**
	 * @return the pcaColumns
	 * @see #pcaColumns
	 */
	public Integer getPcaColumns() {
		return pcaColumns;
	}

	/**
	 * @param pcaColumns
	 *          the pcaColumns to set
	 * @see #pcaColumns
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters pcaColumns(Integer pcaColumns) {
		ProcessParameters result = new ProcessParameters(this);
		result.pcaColumns = pcaColumns;
		return result;
	}

	/**
	 * @return the seed
	 * @see #seed
	 */
	public Integer getSeed() {
		return seed;
	}

	/**
	 * @param seed
	 *          the seed to set
	 * @see #seed
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters seed(Integer seed) {
		ProcessParameters result = new ProcessParameters(this);
		result.seed = seed;
		return result;
	}

	/**
	 * @return the pcaAlgorithmType
	 * @see #pcaAlgorithmType
	 */
	public PcaType getPcaAlgorithmType() {
		return pcaAlgorithmType;
	}

	/**
	 * @param pcaAlgorithmType
	 *          the pcaAlgorithmType to set
	 * @see #pcaAlgorithmType
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters pcaAlgorithmType(PcaType pcaAlgorithmType) {
		ProcessParameters result = new ProcessParameters(this);
		result.pcaAlgorithmType = pcaAlgorithmType;
		return result;
	}

	/**
	 * @return the theta
	 * @see #theta
	 */
	public Double getTheta() {
		return theta;
	}

	/**
	 * @param theta
	 *          the theta to set
	 * @see #theta
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters theta(Double theta) {
		ProcessParameters result = new ProcessParameters(this);
		result.theta = theta;
		return result;
	}

	/**
	 * @return the perplexity
	 * @see #perplexity
	 */
	public Double getPerplexity() {
		return perplexity;
	}

	/**
	 * @param perplexity
	 *          the perplexity to set
	 * @see #perplexity
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters perplexity(Double perplexity) {
		ProcessParameters result = new ProcessParameters(this);
		result.perplexity = perplexity;
		return result;
	}

	/**
	 * @return the merge
	 * @see #merge
	 */
	public Boolean getMerge() {
		return merge;
	}

	/**
	 * @param merge
	 *          the merge to set
	 * @see #merge
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters merge(Boolean merge) {
		ProcessParameters result = new ProcessParameters(this);
		result.merge = merge;
		return result;
	}

	/**
	 * @return the extendedLogs
	 * @see #extendedLogs
	 */
	public Boolean getExtendedLogs() {
		return extendedLogs;
	}

	/**
	 * @param extendedLogs
	 *          the extendedLogs to set
	 * @see #extendedLogs
	 * @return new {@link ProcessParameters} object with the parameter set
	 */
	public ProcessParameters extendedLogs(Boolean extendedLogs) {
		ProcessParameters result = new ProcessParameters(this);
		result.extendedLogs = extendedLogs;
		return result;
	}

}
