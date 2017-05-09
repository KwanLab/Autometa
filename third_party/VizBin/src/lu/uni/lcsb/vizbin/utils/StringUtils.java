package lu.uni.lcsb.vizbin.utils;

import lu.uni.lcsb.vizbin.InvalidArgumentException;

import org.apache.log4j.Logger;

public class StringUtils {

	/**
	 * Default constructor for utility class. Prevents instatiation.
	 */
	private StringUtils() {

	}

	/**
	 * Default class logger.
	 */
	@SuppressWarnings("unused")
	private final Logger	logger	= Logger.getLogger(StringUtils.class);

	public static String nextString(String current) {
		if (current.equals("Z") || current.length() == 0) {
			return null;
		}

		if (current.endsWith("Z")) {
			String result = nextString(current.substring(0, current.length() - 1));
			if (result == null) {
				return null;
			} else {
				return result + "A";
			}
		} else {
			int charValue = current.charAt(current.length() - 1);
			String next = String.valueOf((char) (charValue + 1));
			return current.substring(0, current.length() - 1) + next;
		}

	}

	public static String nextOverAlphabet(String current, String alphabet) {
		String lastLetter = alphabet.substring(alphabet.length() - 1, alphabet.length());

		if (current.equals(lastLetter) || current.length() == 0) {
			return null;
		}

		if (current.endsWith(lastLetter)) {
			String result = nextOverAlphabet(current.substring(0, current.length() - 1), alphabet);
			if (result == null) {
				return null;
			} else {
				return result + alphabet.substring(0, 1);
			}
		} else {
			String currentLetter = current.substring(current.length() - 1, current.length());
			int index = alphabet.indexOf(currentLetter);

			String nextLetter = alphabet.substring(index + 1, index + 2);

			return current.substring(0, current.length() - 1) + nextLetter;
		}

	}

	public static String reverseDnaSequence(String sequence) {
		String result = "";
		for (int i = sequence.length() - 1; i >= 0; i--) {
			Character ch = sequence.charAt(i);
			Character nextChar = null;
			if (ch == 'A') {
				nextChar = 'T';
			} else if (ch == 'C') {
				nextChar = 'G';
			} else if (ch == 'G') {
				nextChar = 'C';
			} else if (ch == 'T') {
				nextChar = 'A';
			}
			if (nextChar == null) {
				throw new InvalidArgumentException("Invalid input sequence");
			}
			result += nextChar;
		}
		return result;
	}

}
