package lu.uni.lcsb.vizbin;

public class InvalidArgumentException extends RuntimeException {

	public InvalidArgumentException(String string) {
		super(string);
	}

	public InvalidArgumentException(Exception e) {
		super(e);
	}

	/**
	 * 
	 */
	private static final long	serialVersionUID	= 1L;

}
