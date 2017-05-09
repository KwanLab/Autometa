package lu.uni.lcsb.vizbin;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.MissingOptionException;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.log4j.Logger;

/**
 * This class is used for processing command line options.
 *
 * @author Piotr Gawron
 *
 */
public class CommandLineOptions {
	/**
	 * Default class logger.
	 */
	private final Logger				logger						= Logger.getLogger(CommandLineOptions.class);
	/**
	 * One letter option name in command line for setting number of threads.
	 */
	private static final String	THREAD_PARAM			= "t";

	/**
	 * One letter option name in command line for setting kmer length.
	 */
	private static final String	K_MER_PARAM				= "k";

	/**
	 * One letter option name in command line for setting sequence length
	 * filtering cutoff.
	 */
	private static final String	CUTOFF_PARAM			= "c";

	/**
	 * One letter option name in command line for setting name of output file.
	 */
	private static final String	OUTPUT_FILE_PARAM	= "o";

	/**
	 * One letter option name in command line for setting name of input file.
	 */
	private static final String	INPUT_FILE_PARAM	= "i";

	/**
	 * One letter option name in command line for setting preplexity parameter.
	 */
	private static final String	PERPLEXITY_PARAM	= "p";

	/**
	 * One letter option name in command line for setting number of PCA columns.
	 */
	private static final String	PCA_COLUMNS_PARAM	= "a";

	/**
	 * Header that should be printed in {@link CommandLineOptions#printHelp()}
	 * method.
	 */
	private static final String	HEADER						= "";
	/**
	 * Footer that should be printed in {@link CommandLineOptions#printHelp()}
	 * method.
	 */
	private static final String	FOOTER						= "";

	/**
	 * This object contains structure of the parameters.
	 */
	private Options							options;

	/**
	 * This object contains processed input data.
	 */
	private CommandLine					cmd;

	/**
	 * Was the input data provided to this class valid or not.
	 */
	private boolean							ok;

	/**
	 * Default constructor.
	 *
	 * @param args
	 *          parameters passed to the program
	 * @throws ParseException
	 *           thrown when there is a problem with parameters
	 */
	public CommandLineOptions(String[] args) throws ParseException {
		options = new Options();
		options.addOption(createOption(false, false, "h", "help", "print this help menu", null));
		options.addOption(createOption(true, true, INPUT_FILE_PARAM, "input", "input file in fasta format", "input-file"));
		options.addOption(createOption(true, true, OUTPUT_FILE_PARAM, "output", "output file containing coordinates", "output-file"));
		options.addOption(createOption(true, false, CUTOFF_PARAM, "cut-off", "minimum conting length [default=" + Config.DEFAULT_CONTIG_LENGTH + "]", "cut-off"));
		options.addOption(createOption(true, false, K_MER_PARAM, "k-mer", "k-mer length [default=" + Config.DEFAULT_KMER_LENGTH + "]", "length"));
		options.addOption(createOption(true, false, THREAD_PARAM, "thread", "number of threads [default=" + Config.DEFAULT_THREAD_NUM + "]", "number"));
		options.addOption(createOption(true, false, PERPLEXITY_PARAM, "perplexity", "perplexity parameter [default=" + Config.DEFAULT_PERPLEXILITY + "]", "number"));
		options.addOption(createOption(true, false, PCA_COLUMNS_PARAM, "pca", "number of PCA columns [default=" + Config.DEFAULT_PCA_COLUMNS + "]", "number"));

		CommandLineParser parser = new BasicParser();
		try {
			cmd = parser.parse(options, args);
			ok = true;
		} catch (MissingOptionException e) {
			ok = false;
		}

	}

	/**
	 * Creates {@link Option} object.
	 *
	 * @param arg
	 *          has the {@link Option} argument
	 * @param required
	 *          is the {@link Option} required
	 * @param abbreviation
	 *          abbreviation of the {@link Option}
	 * @param name
	 *          long name of the {@link Option}
	 * @param description
	 *          description of the {@link Option}
	 * @param paramName
	 *          what should be the name of argument
	 * @return {@link Option} created from input data
	 */
	private Option createOption(boolean arg, boolean required, String abbreviation, String name, String description, String paramName) {
		OptionBuilder.hasArg(arg);
		OptionBuilder.isRequired(required);
		OptionBuilder.withDescription(description);
		OptionBuilder.withLongOpt(name);
		if (paramName != null) {
			OptionBuilder.withArgName(paramName);
		} else {
			OptionBuilder.withArgName(null);
		}
		Option option = OptionBuilder.create(abbreviation);
		return option;
	}

	/**
	 *
	 * @return {@link #ok}
	 */
	public boolean isValid() {
		return ok;
	}

	/**
	 * Print help about command line (information how to pass arguments into it).
	 */
	public void printHelp() {
		HelpFormatter formatter = new HelpFormatter();
		System.out.println("Console usage: ");
		System.out.println("");
		formatter.printHelp("myapp", HEADER, options, FOOTER, true);

	}

	/**
	 * Create {@link ProcessParameters} object from the input data.
	 *
	 * @return {@link ProcessParameters} object from the input data
	 */
	public ProcessParameters getParameters() {
		ProcessParameters params = new ProcessParameters();
		if (cmd.getOptionValue(INPUT_FILE_PARAM) != null) {
			params = params.inputFastaFile(cmd.getOptionValue(INPUT_FILE_PARAM));
		}
		if (cmd.getOptionValue(OUTPUT_FILE_PARAM) != null) {
			params = params.outputFile(cmd.getOptionValue(OUTPUT_FILE_PARAM));
		}
		if (cmd.getOptionValue(CUTOFF_PARAM) != null) {
			try {
				Integer i = Integer.valueOf(cmd.getOptionValue(CUTOFF_PARAM));
				params = params.contigLength(i);
			} catch (NumberFormatException e) {
				logger.warn("Problem with parameter: " + cmd.getOptionValue(CUTOFF_PARAM) + ". Integer value expected. Using default.");
			}
		}
		if (cmd.getOptionValue(K_MER_PARAM) != null) {
			try {
				Integer i = Integer.valueOf(cmd.getOptionValue(K_MER_PARAM));
				params = params.kMerLength(i);
			} catch (NumberFormatException e) {
				logger.warn("Problem with parameter: " + cmd.getOptionValue(K_MER_PARAM) + ". Integer value expected. Using default.");
			}
		}
		if (cmd.getOptionValue(THREAD_PARAM) != null) {
			try {
				Integer i = Integer.valueOf(cmd.getOptionValue(THREAD_PARAM));
				params = params.threads(i);
			} catch (NumberFormatException e) {
				logger.warn("Problem with parameter: " + cmd.getOptionValue(THREAD_PARAM) + ". Integer value expected. Using default.");
			}
		}
		if (cmd.getOptionValue(PERPLEXITY_PARAM) != null) {
			try {
				Double d = Double.valueOf(cmd.getOptionValue(PERPLEXITY_PARAM));
				params = params.perplexity(d);
			} catch (NumberFormatException e) {
				logger.warn("Problem with parameter: " + cmd.getOptionValue(PERPLEXITY_PARAM) + ". Double value expected. Using default.");
			}
		}
		if (cmd.getOptionValue(PCA_COLUMNS_PARAM) != null) {
			try {
				Integer i = Integer.valueOf(cmd.getOptionValue(PCA_COLUMNS_PARAM));
				params = params.pcaColumns(i);
			} catch (NumberFormatException e) {
				logger.warn("Problem with parameter: " + cmd.getOptionValue(PCA_COLUMNS_PARAM) + ". Integer value expected. Using default.");
			}
		}
		return params;
	}
}
