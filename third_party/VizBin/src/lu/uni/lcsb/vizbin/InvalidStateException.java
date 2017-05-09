package lu.uni.lcsb.vizbin;

public class InvalidStateException extends RuntimeException {

	public InvalidStateException(String string) {
		super(string);
	}

	public InvalidStateException(Exception e) {
		super(e);
	}

	/**
	 * 
	 */
	private static final long	serialVersionUID	= 1L;

}
