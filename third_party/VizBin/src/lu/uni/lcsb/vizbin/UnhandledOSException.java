package lu.uni.lcsb.vizbin;

/**
 * Exception that should be thrown where the current operating system is not
 * supported by the tool.
 * 
 * @author Piotr Gawron
 * 
 */
public class UnhandledOSException extends Exception {

	/**
	 * 
	 */
	private static final long	serialVersionUID	= 1L;

	/**
	 * Default constructor.
	 * 
	 * @param message
	 *          exception meesage
	 */
	public UnhandledOSException(String message) {
		super(message);
	}

}
