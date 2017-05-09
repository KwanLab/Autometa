package lu.uni.lcsb.vizbin;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;

/**
 * This class represents zip workspace with all data computed from input file.
 * 
 * @author Piotr Gawron
 * 
 */
public class ZipProject {
	/**
	 * Size of the buffer that will be used to read files.
	 */
	private static final int		READ_BUFFER_SIZE	= 1024;

	/**
	 * Default class logger.
	 */
	@SuppressWarnings("unused")
	private final Logger				logger						= Logger.getLogger(ZipProject.class);

	/**
	 * Filename in zip archive where {@link #dataInput} is stored.
	 */
	private static final String	DATA_FILE_NAME		= "data.fa";
	/**
	 * Filename in zip archive where {@link #labelInput} is stored.
	 */
	private static final String	LABEL_FILE_NAME		= "label.txt";
	/**
	 * Filename in zip archive where {@link #pointInput} is stored.
	 */
	private static final String	POINT_FILE_NAME		= "point.txt";
	/**
	 * Filename in zip archive where {@link #logInput} is stored.
	 */
	private static final String	LOG_FILE_NAME			= "log.txt";

	/**
	 * Filename with the input data to be included in the zip file.
	 */
	private String							dataInput;
	/**
	 * Filename with the labels of input sequences to be included in the zip file.
	 */
	private String							labelInput;
	/**
	 * Filename with the points for sequence to be included in the zip file.
	 */
	private String							pointInput;
	/**
	 * Filename with the log file to be included in the zip file.
	 */
	private String							logInput;

	/**
	 * Default constructor from zip file.
	 * 
	 * @param zipFileName
	 *          zip archive where the workspace is stored
	 * @throws IOException
	 *           thrown when there is a problem with input file
	 */
	public ZipProject(String zipFileName) throws IOException {
		ZipInputStream zis = new ZipInputStream(new FileInputStream(zipFileName));
		ZipEntry ze = zis.getNextEntry();
		byte[] buffer = new byte[READ_BUFFER_SIZE];
		while (ze != null) {
			String fileName = ze.getName();

			String tmpFile = File.createTempFile("temp", Long.toString(System.nanoTime())).getCanonicalPath() + ".txt";

			if (fileName.equals(DATA_FILE_NAME)) {
				dataInput = tmpFile;
			} else if (fileName.equals(LABEL_FILE_NAME)) {
				labelInput = tmpFile;
			} else if (fileName.equals(POINT_FILE_NAME)) {
				pointInput = tmpFile;
			} else if (fileName.equals(LOG_FILE_NAME)) {
				logInput = tmpFile;
			} else {
				zis.close();
				throw new InvalidArgumentException("Zip file contains unknown file: " + fileName);
			}

			File newFile = new File(tmpFile);

			// create all non exists folders
			// else you will hit FileNotFoundException for compressed folder
			new File(newFile.getParent()).mkdirs();

			FileOutputStream fos = new FileOutputStream(newFile);

			int len;
			while ((len = zis.read(buffer)) > 0) {
				fos.write(buffer, 0, len);
			}

			fos.close();
			ze = zis.getNextEntry();
			newFile.deleteOnExit();
		}
		zis.closeEntry();
		zis.close();
	}

	/**
	 * Default constructor that creates zip from input files.
	 * 
	 * @param dataFileName
	 *          file with sequences
	 * @param labelFileName
	 *          file with labels
	 * @param pointFileName
	 *          file with computed points
	 * @param logFileName
	 *          log file
	 */
	public ZipProject(String dataFileName, String labelFileName, String pointFileName, String logFileName) {
		if (dataFileName == null) {
			throw new InvalidArgumentException("Data file name cannot be null");
		}
		if (pointFileName == null) {
			throw new InvalidArgumentException("Point file name cannot be null");
		}

		if (!fileExists(dataFileName)) {
			throw new InvalidArgumentException("Data file (" + dataFileName + ") doesn't exist.");
		}

		if (labelFileName != null && !fileExists(labelFileName)) {
			if (!labelFileName.equals("")) {
				throw new InvalidArgumentException("Label file (" + labelFileName + ") doesn't exist.");
			} else {
				labelFileName = null;
			}
		}

		if (!fileExists(pointFileName)) {
			throw new InvalidArgumentException("Point file (" + pointFileName + ") doesn't exist.");
		}

		if (logFileName != null && !fileExists(logFileName)) {
			if (!logFileName.equals("")) {
				throw new InvalidArgumentException("Log file (" + logFileName + ") doesn't exist.");
			} else {
				logFileName = null;
			}
		}

		dataInput = dataFileName;
		labelInput = labelFileName;
		pointInput = pointFileName;
		logInput = logFileName;
	}

	/**
	 * Checks if file exists.
	 * 
	 * @param fileName
	 *          file to check
	 * @return <code>true</code> if the file in parameter exists
	 */
	private boolean fileExists(String fileName) {
		return new File(fileName).exists();
	}

	/**
	 * Returns {@link InputStream} with input data (sequences).
	 * 
	 * @return {@link InputStream} with input data (sequences)
	 */
	public InputStream getDataInput() {
		if (dataInput == null) {
			return null;
		} else {
			try {
				return new FileInputStream(dataInput);
			} catch (FileNotFoundException e) {
				throw new InvalidStateException("File that should be there doesn't exist anymore...: " + dataInput);
			}
		}
	}

	/**
	 * Returns {@link InputStream} with labels or null if such file is not
	 * available.
	 * 
	 * @return {@link InputStream} with labels
	 */
	public InputStream getLabelInput() {
		if (labelInput == null) {
			return null;
		} else {
			try {
				return new FileInputStream(labelInput);
			} catch (FileNotFoundException e) {
				throw new InvalidStateException("File that should be there doesn't exist anymore...: " + labelInput);
			}
		}
	}

	/**
	 * Returns {@link InputStream} with points.
	 * 
	 * @return {@link InputStream} with points
	 */
	public InputStream getPointInput() {
		if (pointInput == null) {
			return null;
		} else {
			try {
				return new FileInputStream(pointInput);
			} catch (FileNotFoundException e) {
				throw new InvalidStateException("File that should be there doesn't exist anymore...: " + pointInput);
			}
		}
	}

	/**
	 * Returns {@link InputStream} with log or null if such file is not available.
	 * 
	 * @return {@link InputStream} with log
	 */
	public InputStream getLogInput() {
		if (logInput == null) {
			return null;
		} else {
			try {
				return new FileInputStream(logInput);
			} catch (FileNotFoundException e) {
				throw new InvalidStateException("File that should be there doesn't exist anymore...: " + logInput);
			}
		}
	}

	/**
	 * Saves project to file.
	 * 
	 * @param file
	 *          where the project should be saved
	 * @throws IOException
	 *           thrown when there is a problem with files (input or output)
	 */
	public void saveTo(String file) throws IOException {
		FileOutputStream fos = new FileOutputStream(file);
		ZipOutputStream zos = new ZipOutputStream(fos);

		InputStream di = getDataInput();
		ZipEntry diEntry = new ZipEntry(DATA_FILE_NAME);
		zos.putNextEntry(diEntry);
		IOUtils.copy(di, zos);
		zos.closeEntry();

		InputStream li = getLabelInput();
		if (li != null) {
			ZipEntry liEntry = new ZipEntry(LABEL_FILE_NAME);
			zos.putNextEntry(liEntry);
			IOUtils.copy(li, zos);
			zos.closeEntry();
		}

		InputStream pi = getPointInput();
		ZipEntry piEntry = new ZipEntry(POINT_FILE_NAME);
		zos.putNextEntry(piEntry);
		IOUtils.copy(pi, zos);
		zos.closeEntry();

		InputStream logi = getLogInput();
		if (logi != null) {
			ZipEntry liEntry = new ZipEntry(LOG_FILE_NAME);
			zos.putNextEntry(liEntry);
			IOUtils.copy(logi, zos);
			zos.closeEntry();
		}
		zos.close();
	}

	/**
	 * @return the dataInput
	 * @see #dataInput
	 */
	public String getDataInputFile() {
		return dataInput;
	}

	/**
	 * @return the labelInput
	 * @see #labelInput
	 */
	public String getLabelInputFile() {
		return labelInput;
	}

	/**
	 * @return the pointInput
	 * @see #pointInput
	 */
	public String getPointInputFile() {
		return pointInput;
	}

	/**
	 * @return the logInput
	 * @see #logInput
	 */
	public String getLogInputFile() {
		return logInput;
	}

}
