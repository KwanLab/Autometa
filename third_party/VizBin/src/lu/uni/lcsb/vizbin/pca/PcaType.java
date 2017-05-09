package lu.uni.lcsb.vizbin.pca;

import org.apache.log4j.Logger;

/**
 * Defines implementation of the PCA algorithm that should be used by the client
 * side.
 * 
 * @author Piotr Gawron
 * 
 */
public enum PcaType {
	/**
	 * Mtj implementation.
	 * 
	 * @see PrincipleComponentAnalysisMtj
	 */
	MTJ("Mtj", PrincipleComponentAnalysisMtj.class),
	/**
	 * Mtj implementation with some optimisations.
	 */
	MTJ_OPTIMIZED("Mtj optimized", PrincipleComponentAnalysisMtj.class),
	/**
	 * EJML implementation with some optimisations.
	 * 
	 * @see PrincipleComponentAnalysisEJML
	 */
	EJML("EJML", PrincipleComponentAnalysisEJML.class);

	/**
	 * Name of the pca implementation that is human readable.
	 */
	private String																				commonName	= null;

	/**
	 * Default class logger.
	 */
	private final Logger																	logger			= Logger.getLogger(PcaType.class);
	/**
	 * Class with the implementation.
	 */
	private Class<? extends IPrincipleComponentAnalysis>	clazz;

	/**
	 * Default constructor.
	 * 
	 * @param name
	 *          {@link #commonName}
	 * @param clazz
	 *          {@link #clazz}
	 */
	private PcaType(String name, Class<? extends IPrincipleComponentAnalysis> clazz) {
		this.commonName = name;
		this.clazz = clazz;
	}

	/**
	 * Returns human readable name of the implementation.
	 * 
	 * @return {@link #commonName}
	 */
	public String getName() {
		return commonName;
	}

	/**
	 * Creates object that implements pca.
	 * 
	 * @return implementation of {@link IPrincipleComponentAnalysis}
	 */
	public IPrincipleComponentAnalysis getInstance() {
		try {
			return clazz.newInstance();
		} catch (InstantiationException | IllegalAccessException e) {
			logger.fatal("Problem with creating IPrincipleComponentAnalysis implementation: " + clazz, e);
			return null;
		}
	}
}
