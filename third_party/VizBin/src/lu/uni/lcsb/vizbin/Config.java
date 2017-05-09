package lu.uni.lcsb.vizbin;

import lu.uni.lcsb.vizbin.pca.PcaType;

/**
 * Class with config parameters.
 * 
 * @author Piotr Gawron
 * 
 */
public final class Config {

	/**
	 * Default constructor for utility class. Prevents instatiation.
	 */
	private Config() {

	}

	/**
	 * Default minimum contig length.
	 */
	public static final Integer	DEFAULT_CONTIG_LENGTH	= 1000;
	/**
	 * Default number of threads.
	 */
	public static final Integer	DEFAULT_THREAD_NUM		= 1;
	/**
	 * Default k-mer length.
	 */
	public static final Integer	DEFAULT_KMER_LENGTH		= 5;
	/**
	 * Default number of columns to which
	 * {@link lu.uni.lcsb.vizbin.pca.IPrincipleComponentAnalysis} should reduce
	 * input data.
	 */
	public static final Integer	DEFAULT_PCA_COLUMNS		= 50;
	/**
	 * Default theta value used by tsne algorithm.
	 */
	public static final Double	DEFAULT_THETA					= 0.5;
	/**
	 * Should the kmers and "reversed and complementary kmers" be merged into one.
	 */
	public static final Boolean	DEFAULT_MERGE					= true;
	/**
	 * Default perplexity value used by tsne algorithm.
	 */
	public static final Double	DEFAULT_PERPLEXILITY	= 30.;
	/**
	 * Default seed for random function.
	 */
	public static final Integer	DEFAULT_SEED					= 0;
	/**
	 * Default {@link PcaType PCA} algorithm to be used.
	 */
	public static final PcaType	DEFAULT_PCA_TYPE			= PcaType.MTJ;

}
