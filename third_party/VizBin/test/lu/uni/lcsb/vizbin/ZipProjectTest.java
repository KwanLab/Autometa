package lu.uni.lcsb.vizbin;

import static org.junit.Assert.assertTrue;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;

import org.apache.commons.io.IOUtils;
import org.apache.log4j.Logger;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class ZipProjectTest {
	Logger	logger	= Logger.getLogger(ZipProjectTest.class);

	@Before
	public void setUp() throws Exception {
	}

	private String createTemp() throws IOException {
		String fileName = File.createTempFile("temp", Long.toString(System.nanoTime())).getCanonicalPath() + ".txt";
		PrintWriter writer = new PrintWriter(fileName, "UTF-8");
		writer.println(Math.random());
		writer.println("The second line");
		writer.close();
		return fileName;
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void testCreateFromFiles() throws Exception {
		try {
			String file1 = createTemp();
			String file2 = createTemp();
			String file3 = createTemp();
			String file4 = createTemp();

			String fileName = File.createTempFile("temp", Long.toString(System.nanoTime())).getCanonicalPath() + ".zip";

			ZipProject zipProject = new ZipProject(file1, file2, file3, file4);

			zipProject.saveTo(fileName);

			File zip = new File(fileName);
			assertTrue(zip.exists());

			new File(file1).delete();
			new File(file2).delete();
			new File(file3).delete();
			new File(file4).delete();
			new File(fileName).delete();
		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}

	@Test
	public void testCreateFromZip() throws Exception {
		try {
			ZipProject zipProject = new ZipProject("testFiles/zip/ZipProjectTest.zip");
			String fileName = File.createTempFile("temp", Long.toString(System.nanoTime())).getCanonicalPath() + ".zip";

			zipProject.saveTo(fileName);

			ZipProject zipProject2 = new ZipProject(fileName);

			assertTrue(IOUtils.contentEquals(zipProject.getDataInput(), zipProject2.getDataInput()));
			assertTrue(IOUtils.contentEquals(zipProject.getLabelInput(), zipProject2.getLabelInput()));
			assertTrue(IOUtils.contentEquals(zipProject.getLogInput(), zipProject2.getLogInput()));
			assertTrue(IOUtils.contentEquals(zipProject.getPointInput(), zipProject2.getPointInput()));

			new File(fileName).delete();
		} catch (Exception e) {
			e.printStackTrace();
			throw e;
		}
	}

	boolean compareFiles(String filename1, String filename2) throws IOException {
		File file1 = new File(filename1);
		File file2 = new File(filename2);
		if (file1.length() != file2.length()) {
			return false;
		}

		try (InputStream in1 = new BufferedInputStream(new FileInputStream(file1)); InputStream in2 = new BufferedInputStream(new FileInputStream(file2));) {

			int value1, value2;
			do {
				// since we're buffered read() isn't expensive
				value1 = in1.read();
				value2 = in2.read();
				if (value1 != value2) {
					return false;
				}
			} while (value1 >= 0);

			// since we already checked that the file sizes are equal
			// if we're here we reached the end of both files without a mismatch
			return true;
		}
	}

}
