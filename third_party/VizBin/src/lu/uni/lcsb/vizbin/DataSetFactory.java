package lu.uni.lcsb.vizbin;

import java.awt.geom.Point2D;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import javax.swing.JOptionPane;

import lu.uni.lcsb.vizbin.data.DataSet;
import lu.uni.lcsb.vizbin.data.Sequence;

import org.apache.log4j.Logger;

public class DataSetFactory {
	/**
	 * Default constructor for utility class. Prevents instatiation.
	 */
	private DataSetFactory() {

	}

	/**
	 * Default class logger.
	 */
	private static Logger	logger	= Logger.getLogger(DataSetFactory.class);

	public static DataSet createDataSetFromPointFile(String fileName, String labelFileName, double scale, boolean log) throws IOException,
			InvalidMetaFileException {
		InputStream labelIS = null;
		if (labelFileName != null) {
			if (new File(labelFileName).exists()) {
				labelIS = new FileInputStream(new File(labelFileName));
			}
		}
		return createDataSetFromPointFile(new FileInputStream(fileName), labelIS, scale, log);
	}

	public static DataSet createDataSetFromPointFile(InputStream is, InputStream labelIS, double scale, boolean log) throws IOException,
			InvalidMetaFileException {
		Double minCoverage = null;
		Double maxCoverage = null;
		Double maxGc = null;

		DataSet result = new DataSet();
		int i = 0;
		Double minLength = null;

		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		try {
			String line = br.readLine();
			while (line != null) {
				Sequence sequence = new Sequence();
				String coordinates[] = line.split(",");
				if (coordinates.length == 2) {
					double x = Double.parseDouble(coordinates[0]);
					double y = Double.parseDouble(coordinates[1]);
					sequence.setLocation(new Point2D.Double(x * scale, y * scale));
					sequence.setId("Seq" + (i++));
					result.addSequence(sequence);
				}
				line = br.readLine();
			}
		} finally {
			br.close();
		}
		if (labelIS != null) {

			Integer labelColumn = null;
			Integer coverageColumn = null;
			Integer markerColumn = null;
			Integer gcColumn = null;
			Integer lengthColumn = null;

			int columns = 0;

			Map<String, Integer> labels = new HashMap<String, Integer>();
			br = new BufferedReader(new InputStreamReader(labelIS));
			try {

				String headerLine = br.readLine();
				String header[] = headerLine.split(",");
				columns = header.length;
				int colId = 0;
				for (String string : header) {
					if ("coverage".equalsIgnoreCase(string)) {
						coverageColumn = colId;
						maxCoverage = 0.0;
						minCoverage = Double.MAX_VALUE;
					} else if ("gc".equalsIgnoreCase(string)) {
						gcColumn = colId;
						maxGc = 0.0;
					} else if ("label".equalsIgnoreCase(string)) {
						labelColumn = colId;
					} else if ("isMarker".equalsIgnoreCase(string)) {
						markerColumn = colId;
					} else if ("length".equalsIgnoreCase(string)) {
						lengthColumn = colId;
						minLength = Double.MAX_VALUE;
					} else {
						throw new InvalidMetaFileException("Unknonw column header: " + string);
					}
					colId++;
				}

				String line = br.readLine();
				int id = 0;
				while (line != null && id < result.getSequences().size()) {

					Sequence sequence = result.getSequences().get(id);

					String label = "unknown";

					String row = line;
					String data[] = row.split(",");
					if (data.length != columns) {
						throw new InvalidMetaFileException("Line: \"" + row + "\" cotains invalid number of fields: " + columns + ", but " + data.length + " found.");
					}
					if (labelColumn != null) {
						label = data[labelColumn];
					}
					if (gcColumn != null) {
						String string = data[gcColumn];
						try {
							Double val = Double.valueOf(string);
							sequence.setGc(val);
							maxGc = Math.max(val, maxGc);
						} catch (NumberFormatException e) {
							throw new InvalidMetaFileException("Problem with parsing gc value: " + string + ". Double expected.");
						}
					}
					if (markerColumn != null) {
						String string = data[markerColumn];
						if ("0".equals(string)) {
							sequence.setMarker(false);
						} else if ("1".equals(string)) {
							sequence.setMarker(true);
						} else {
							throw new InvalidMetaFileException("Problem with parsing marker value: " + string + ". {0,1} expected.");
						}
					}
					if (coverageColumn != null) {
						String string = data[coverageColumn];
						try {
							Double val = Double.valueOf(string);
							if (log) {
								val = Math.log(val);
							}
							sequence.setCoverage(val);
							maxCoverage = Math.max(val, maxCoverage);
							if (val > 0) {
								minCoverage = Math.min(val, minCoverage);
							}
						} catch (NumberFormatException e) {
							throw new InvalidMetaFileException("Problem with parsing coverage value: " + string + ". Double expected.");
						}
					}
					if (lengthColumn != null) {
						String string = data[lengthColumn];
						try {
							Double val = Double.valueOf(string);
							if (log) {
								val = Math.log(val);
							}
							sequence.setLength(val);
							minLength = Math.min(val, minLength);
						} catch (NumberFormatException e) {
							throw new InvalidMetaFileException("Problem with parsing length value: " + string + ". Double expected.");
						}
					}

					Integer labelId = labels.get(label);
					if (labelId == null) {
						labelId = labels.size();
						labels.put(label, labelId);
					}
					sequence.setLabelId(labelId);
					sequence.setLabelName(label);
					line = br.readLine();
					id++;
				}
			} finally {
				br.close();
			}
		}
		for (Sequence sequence : result.getSequences()) {
			if (maxCoverage != null) {
				sequence.setCoverage(Math.max(minCoverage, sequence.getCoverage()) / maxCoverage);
			}
			if (maxGc != null) {
				sequence.setGc(sequence.getGc() / maxGc);
			}
		}

		// TODO Insert an assert() or so to make sure that IF length information is
		// provided is must be > 0
		if (minLength != null) {
			// Set the minimum length of the entire dataset
			result.setMinSequenceLength(minLength.doubleValue());
		}

		return result;
	}

	public static DataSet createDataSetFromFastaFile(String fileName, String filteredSequencesFileName, String labelFileName, String pointsFileName,
			Integer contigLen, boolean log, ProcessGuiParameters guiParameters) throws IOException, InvalidMetaFileException {

		Collection<Integer> filteredSequences = filterSequences(
				new FileInputStream(fileName), new FileOutputStream(filteredSequencesFileName), contigLen, guiParameters);

		FileInputStream sequencesFis = new FileInputStream(filteredSequencesFileName);

		FileInputStream labelFis = null;
		if (labelFileName != null) {
			labelFis = new FileInputStream(labelFileName);
		}

		FileInputStream pointsFis = null;
		if (pointsFileName != null) {
			pointsFis = new FileInputStream(pointsFileName);
		}

		return createDataSetFromFastaFile(sequencesFis, labelFis, pointsFis, log, filteredSequences);

	}

	public static int getLabelID(ArrayList<String> labels, String label) {
		for (int i = 0; i < labels.size(); i++) {
			if (labels.get(i).equals(label))
				return i;
		}
		return -1;
	}

	/**
	 * Filters input seuences based on contig length.
	 * 
	 * @param is
	 *          input stream from which data is read
	 * @param os
	 *          output stream where the filtered data is saved
	 * @param contigLen
	 *          minimum conting length for sequences that should be saved
	 * @return set of integers defining which sequences were skipped (0-based)
	 * @throws IOException
	 *           thrown when there is a problem with input/output stream
	 */
	public static Set<Integer> filterSequences(InputStream is, OutputStream os, Integer contigLen, ProcessGuiParameters guiParameters) throws IOException {
		Set<Integer> filteresSequences = new HashSet<Integer>();

		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(os));
		String line = br.readLine();
		String label = "";
		StringBuilder dna = new StringBuilder();
		int numShortContigs = 0, numContigs = 0;

		int contigId = 0;
		while (line != null) {
			label = line;
			line = "";
			do {
				dna.append(line);
				line = br.readLine();
			} while (line != null && (!line.startsWith(">")));

			if (dna.toString().length() >= contigLen) {
				bw.write(label + "\n");
				bw.write(dna.toString() + "\n");
				numContigs++;
			} else {
				numShortContigs++;
				filteresSequences.add(contigId);
			}
			dna = new StringBuilder();
			contigId++;
		}

		if (numShortContigs > 0) {
			String message = "Number of contigs longer or equal then " + contigLen + ": " + numContigs + "\n" + "Number of contigs shorter then " + contigLen + ": "
					+ numShortContigs;
			if (guiParameters != null) {
				JOptionPane.showMessageDialog(null, message);
			} else {
				logger.info(message);
			}
		}

		br.close();
		bw.close();
		return filteresSequences;
	}

	public static DataSet createDataSetFromFastaFile(InputStream is, InputStream labelIS, InputStream pointsIS, boolean log,
			Collection<Integer> filteredSequences) throws IOException, InvalidMetaFileException {

		String invalidCoverage = "";
		DataSet result = new DataSet();
		Double maxCoverage = null;
		Double minCoverage = null;
		Double maxGc = null;
		Double minLength = null;
		String label = "";
		ArrayList<String> labels = new ArrayList<String>();

		int numSequences = 0;
		int numPoints = 0;
		int pstart, pend, id;

		BufferedReader br = new BufferedReader(new InputStreamReader(is));
		BufferedReader brLabels = null;

		if (labelIS != null)
			brLabels = new BufferedReader(new InputStreamReader(labelIS));
		BufferedReader brPoints = null;
		if (pointsIS != null)
			brPoints = new BufferedReader(new InputStreamReader(pointsIS));

		try {
			String line;
			Sequence sequence = null;
			line = br.readLine();

			Integer labelColumn = null;
			Integer coverageColumn = null;
			Integer markerColumn = null;
			Integer gcColumn = null;
			Integer lengthColumn = null;

			int columns = 0;

			if (brLabels != null) {
				String headerLine = brLabels.readLine();
				String header[] = headerLine.split(",");
				columns = header.length;
				int colId = 0;
				for (String string : header) {
					if ("coverage".equalsIgnoreCase(string)) {
						coverageColumn = colId;
						maxCoverage = 0.0;
						minCoverage = Double.MAX_VALUE;
					} else if ("gc".equalsIgnoreCase(string)) {
						gcColumn = colId;
						maxGc = 0.0;
					} else if ("label".equalsIgnoreCase(string)) {
						labelColumn = colId;
					} else if ("isMarker".equalsIgnoreCase(string)) {
						markerColumn = colId;
					} else if ("length".equalsIgnoreCase(string)) {
						lengthColumn = colId;
						minLength = Double.MAX_VALUE;
					} else {
						throw new InvalidMetaFileException("Unknonw column header: " + string);
					}
					colId++;
				}
			}

			if (gcColumn != null && coverageColumn != null) {
				throw new InvalidMetaFileException("Either provide 'gc' or 'coverage' values!");
			}
			Integer labelSequenceId = 0;

			while (line != null) {
				sequence = new Sequence();
				sequence.setId("Seq" + (numSequences++));
				if (brLabels != null) {
					String row = brLabels.readLine();
					while (filteredSequences.contains(labelSequenceId)) {
						labelSequenceId++;
						row = brLabels.readLine();
					}
					String data[] = row.split(",");
					if (data.length != columns) {
						throw new InvalidMetaFileException("Line: \"" + row + "\" cotains invalid number of fields: " + columns + ", but " + data.length + " found.");
					}
					if (labelColumn != null) {
						label = data[labelColumn];
					}
					if (gcColumn != null) {
						String string = data[gcColumn];
						try {
							Double val = Double.valueOf(string);
							sequence.setGc(val);
							maxGc = Math.max(val, maxGc);
						} catch (NumberFormatException e) {
							throw new InvalidMetaFileException("Problem with parsing gc value: " + string + ". Double expected.");
						}
					}
					if (markerColumn != null) {
						String string = data[markerColumn];
						if ("0".equals(string)) {
							sequence.setMarker(false);
						} else if ("1".equals(string)) {
							sequence.setMarker(true);
						} else {
							throw new InvalidMetaFileException("Problem with parsing marker value: " + string + ". {0,1} expected.");
						}
					}
					if (coverageColumn != null) {
						String string = data[coverageColumn];
						try {
							Double val = Double.valueOf(string);
							if (log) {
								val = Math.log(val);
								if (val < 0) {
									invalidCoverage = string;
								}
							}
							sequence.setCoverage(val);
							maxCoverage = Math.max(val, maxCoverage);
							if (val > 0) {
								minCoverage = Math.min(val, minCoverage);
							}
						} catch (NumberFormatException e) {
							throw new InvalidMetaFileException("Problem with parsing coverage value: " + string + ". Double expected.");
						}
					}
					if (lengthColumn != null) {
						String string = data[lengthColumn];
						try {
							Double val = Double.valueOf(string);
							if (log) {
								val = Math.log(val);
							}
							sequence.setLength(val);
							minLength = Math.min(val, minLength);
						} catch (NumberFormatException e) {
							throw new InvalidMetaFileException("Problem with parsing length value: " + string + ". Double expected.");
						}
					}
				}
				if (brLabels == null || labelColumn == null) {
					// find and parse description if in `
					// description=" [...] " `format
					if (line.toLowerCase().contains("description=")) {
						pstart = line.indexOf("description=") + "description=".length();
						pstart = line.indexOf('"', pstart) + 1;
						pend = line.indexOf('"', pstart);
						label = line.substring(pstart, pend);
					}
				}

				// System.out.println("LABEL: "+label);
				sequence.setLabelName(label);
				id = getLabelID(labels, label);
				if (id == -1) {
					labels.add(label);
					sequence.setLabelId(labels.size() - 1);
				} else {
					sequence.setLabelId(id);
				}

				sequence.setDna(br.readLine());
				result.addSequence(sequence);

				line = br.readLine();
			}

			if (brPoints != null) {
				while (brPoints.readLine() != null) {
					numPoints++;
				}
			}
		} finally {
			br.close();
			if (brLabels != null)
				brLabels.close();
			if (brPoints != null)
				brPoints.close();
		}
		if (pointsIS != null) {
			if (numSequences != numPoints) {
				logger.debug("Points file does not match fasta file!");
				return null;
			}
		}

		if (result.getSize() == 0) {
			logger.debug("Zero sequences loaded from the file! Check file and parameters.");
			return null;
		}

		for (Sequence sequence : result.getSequences()) {
			if (maxCoverage != null) {
				sequence.setCoverage(Math.max(minCoverage, sequence.getCoverage()) / maxCoverage);
			}
			if (maxGc != null) {
				sequence.setGc(sequence.getGc() / maxGc);
			}
		}
		if (!invalidCoverage.equals("")) {
			JOptionPane.showMessageDialog(null, "WARNING: Unexpected coverage value of " + invalidCoverage
					+ " found.\nWill take the minimum strictly positive input coverage value as default for now.\nPlease verify your input annotation file.");

		}

		// TODO Insert an assert() or so to make sure that IF length information is
		// provided is must be > 0
		if (minLength != null) {
			// Set the minimum length of the entire dataset
			result.setMinSequenceLength(minLength.doubleValue());
		}

		return result;
	}

	public static void saveToPcaFile(DataSet dataSet, String fileName, int numThreads, double theta, double perplexity, int seed) throws FileNotFoundException,
			UnsupportedEncodingException {
		PrintWriter writer = new PrintWriter(fileName, "UTF-8");

		writer.println(dataSet.getSequences().size());
		writer.println(dataSet.getSequences().get(0).getPcaVector().length);
		writer.println(theta);
		writer.println(perplexity);
		writer.println(numThreads);
		writer.println(seed);

		for (Sequence sequence : dataSet.getSequences()) {
			double vector[] = sequence.getPcaVector();
			for (int i = 0; i < vector.length; i++)
				writer.print(vector[i] + " ");
			writer.println();
		}
		writer.close();
	}

}
