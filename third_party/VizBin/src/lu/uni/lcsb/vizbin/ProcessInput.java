package lu.uni.lcsb.vizbin;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.concurrent.atomic.AtomicInteger;

import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import lu.uni.lcsb.vizbin.clustering.ClusterPanel;
import lu.uni.lcsb.vizbin.data.DataSet;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;

/**
 * 
 * @author <a href="mailto:valentin.plugaru.001@student.uni.lu">Valentin
 *         Plugaru</a>
 */
public class ProcessInput extends ObjectWithProperties {

	/**
	 * Maximum value for {@link #progressVal} progress bar.
	 */
	private static final int			MAX_PROGRESS_VALUE		= 100;

	public final static String		FINISHED_PROPERTY			= "FINISHED";

	protected static final String	POINTS_FILE_PROPERTY	= "POINTS_FILE";

	private String								name;

	/**
	 * Default class logger.
	 */
	private final Logger					logger								= Logger.getLogger(ProcessInput.class.getName());

	private String								filteredSequencesFile	= "filteredSequences.fa";
	private String								pointsfile;

	private File									tsneCmd;

	private DataSet								dataSetOrig						= null;
	private DataSet								dataSet								= null;

	private volatile Integer			progressVal;

	private Boolean								processEnded					= true;

	private ProcessParameters			parameters;
	private ProcessGuiParameters	guiParameters;

	ProcessInput(ProcessParameters parameters, ProcessGuiParameters guiParams, File tsneCmd) throws IOException {
		logger.debug("Init of ProcessInput");
		this.parameters = parameters;
		this.guiParameters = guiParams;
		if (parameters.getInputPointFile() != null) {
			setPointsfile(parameters.getInputPointFile());
		}
		this.tsneCmd = tsneCmd;
		this.progressVal = 0;
	}

	void updateStatus(String message) {
		updateStatus(message, 0);
	}

	void updateStatus(final String message, int amount) {
		logger.debug(message);
		progressVal += amount;
		if (progressVal > MAX_PROGRESS_VALUE) {
			progressVal = 100;
		}
		if (progressVal < 0) {
			progressVal = 0;
		}
		final AtomicInteger value = new AtomicInteger(progressVal);

		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				if (guiParameters != null) {
					guiParameters.getStatusLabel().setText(message);
					guiParameters.getProgessBar().setValue(value.get());
				} else {
					logger.debug("[PROGRESS BAR] " + message + " (" + value.get() + ")");
				}
			}
		});

		// Do some work and update value.
	}

	public void doProcess() {
		Thread thread = new Thread() {
			public void run() {

				setProcessEnded(false);
				String localPath = System.getProperty("java.io.tmpdir");
				logger.debug(localPath);
				try {
					File myTempDir = File.createTempFile("map", "", new File(localPath));
					myTempDir.delete();
					myTempDir.mkdirs();

					String directory = myTempDir.getAbsolutePath();
					filteredSequencesFile = directory + "/" + filteredSequencesFile;

					logger.debug("Loading data from file.\nContig length treshold: " + parameters.getContigLength());
					updateStatus("Loading fasta file: " + parameters.getInputFastaFile(), 5);
					dataSet = DataSetFactory.createDataSetFromFastaFile(
							parameters.getInputFastaFile(), filteredSequencesFile, parameters.getInputLabelFile(), parameters.getInputPointFile(),
							parameters.getContigLength(), parameters.getExtendedLogs(), guiParameters);
					if (dataSet == null) {
						JOptionPane.showMessageDialog(null, "Error during loading data from given file! Check the logs.", "alert", JOptionPane.ERROR_MESSAGE);
						updateStatus("", -MAX_PROGRESS_VALUE);
						setProcessEnded(true);
						return;
					}
					// for (int i=0;i<dataSet.getSequences().size(); i++)
					// System.out.println(dataSet.getSequences().get(i).getLabelName());
					updateStatus("DataSet loaded (" + dataSet.getSequences().size() + " sequences)");

					// If points file is provided, no calculations are needed
					if (parameters.getInputPointFile() == null) {
						updateStatus("Creating kmers (k=" + parameters.getkMerLength() + ", merge = " + parameters.getMerge() + ")");
						DataSetUtils.createKmers(dataSet, parameters.getkMerLength(), parameters.getMerge(), guiParameters);
						updateStatus("Normalizing vectors...", 5);
						DataSetUtils.normalizeDescVectors(dataSet, parameters.getkMerLength());
						updateStatus("Clr normalization...", 5);
						DataSetUtils.createClrData(dataSet);
						if (parameters.getKmerDebugFile() != null) {
							DataSetUtils.saveClrData(dataSet, parameters.getKmerDebugFile());
						}
						updateStatus("Running PCA... (" + parameters.getPcaAlgorithmType().getName() + ")", 5);
						DataSetUtils.computePca(dataSet, parameters.getPcaColumns(), parameters.getPcaAlgorithmType());
						updateStatus("Running T-SNE...", 15);
						DataSetUtils.runTsneAndPutResultsToDir(
								dataSet, parameters.getThreads(), directory, parameters.getTheta(), parameters.getPerplexity(), parameters.getSeed(), tsneCmd, guiParameters);

						// progressVal = progBar.getValue();
						File f = new File(directory + "/points.txt");
						if (f.exists() && !f.isDirectory()) {
							updateStatus("Points created."); // 90% progress up
							// to// here
						} else {
							throw new FileNotFoundException("points.txt file not found. Probably bh_tsne binaries execution failed.");
						}

						setPointsfile(directory + "/points.txt");

					}

					File flabels = null;
					FileInputStream labelsIS = null;
					if (parameters.getInputLabelFile() != null) {
						try {
							flabels = new File(parameters.getInputLabelFile());
							labelsIS = new FileInputStream(flabels);
						} catch (Exception e) {
							e.printStackTrace();
						}
					}

					dataSetOrig = dataSet;
					int scale = 10; // Scaling is needed since AWT.Polygon() requires
													// int-coordinates for the polygon vertices.
					// Scaling by a factor of ten allows to zoom in and still get
					// meaningful int-coordinates form double points.
					// TODO: Refactor such that this is done internally in
					// ClusterFactory.createClusterFromPolygon()
					dataSet = DataSetFactory.createDataSetFromPointFile(new FileInputStream(pointsfile), labelsIS, scale, parameters.getExtendedLogs());
					updateStatus("Done.", 100); // 100% progress, make sure
					// progress bar is at 100

					// add label ID and name from initial dataset
					for (int i = 0; i < dataSet.getSequences().size(); i++) {
						dataSet.getSequences().get(i).setLabelId(dataSetOrig.getSequences().get(i).getLabelId());
						dataSet.getSequences().get(i).setLabelName(dataSetOrig.getSequences().get(i).getLabelName());
					}

					DataSetUtils.setDataSet(dataSet);
					DataSetUtils.setIsDataSetCreated(true);

					if (guiParameters != null) {
						ClusterPanel cp = new ClusterPanel(dataSet, filteredSequencesFile, guiParameters.getParentFrame());
						guiParameters.getTabPane().setComponentAt(1, cp.getChartPanel());
						// NotificationCenter.addObserver(cp);
						// Show data points
						JFrame frame = new JFrame("Cluster " + name);
						DataSetUtils.setDrawingFrame(frame);
						// frame.setDefaultCloseOperation (JFrame.EXIT_ON_CLOSE);
						ClusterPanel cpPopOut = new ClusterPanel(dataSet, filteredSequencesFile, guiParameters.getParentFrame());
						// NotificationCenter.addObserver(cpPopOut);
						// frame.getContentPane().add(cpPopOut);

						frame.setSize(800, 600);
						JScrollPane scrPane = new JScrollPane(JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
						scrPane.getViewport().add(cpPopOut.getChartPanel());
						frame.setContentPane(scrPane);
						frame.setVisible(true);
					}

				} catch (OutOfMemoryError e) {
					showMessageDialog("Error! Java machine ran out of memmory.\n" + "Check input file size, or increase java heap size.\n"
							+ "Application will now restart.");
					e.printStackTrace();
					restartApplication();
				} catch (InvalidMetaFileException e) {
					logger.error(e.getMessage(), e);
					updateStatus("Error! Check the logs.");
					showMessageDialog(e.getMessage(), "Label file error", JOptionPane.ERROR_MESSAGE);
				} catch (Exception e) {
					e.printStackTrace();
					updateStatus("Error! Check the logs.");
					logger.error(e.getMessage(), e);
				}
				setProcessEnded(true);
			}

			private void showMessageDialog(String message, String title, int errorMessage) {
				if (guiParameters != null) {
					JOptionPane.showMessageDialog(guiParameters.getParentFrame(), message, title, errorMessage);
				} else {
					logger.fatal("[" + title + "][" + errorMessage + "] " + message);
				}
			}

			private void showMessageDialog(String message) {
				if (guiParameters != null) {
					JOptionPane.showMessageDialog(guiParameters.getParentFrame(), message);
				} else {
					logger.info(message);
				}

			}
		};
		thread.start();
	}

	public void restartApplication() {
		final String javaBin = System.getProperty("java.home") + File.separator + "bin" + File.separator + "java";
		File currentJar = null;
		currentJar = new File(ProcessInput.class
				.getProtectionDomain().getCodeSource().getLocation().toString().replace("build/classes", "dist/VizBin-dist.jar ").replace("file:", ""));

		logger.error("Jar: " + currentJar);

		/* Build command: java -jar application.jar */
		final ArrayList<String> command = new ArrayList<String>();
		command.add(javaBin);
		command.add("-jar");
		command.add(currentJar.getPath());

		logger.error("Command: " + command);

		final ProcessBuilder builder = new ProcessBuilder(command);

		try {
			builder.start();
		} catch (IOException e) {
			e.printStackTrace();
		}
		System.exit(0);
	}

	public DataSet getDataSet() {
		return dataSet;
	}

	/**
	 * @param processEnded
	 *          the processEnded to set
	 * @see #processEnded
	 */
	private void setProcessEnded(Boolean processEnded) {
		Boolean oldValue = this.processEnded;
		this.processEnded = processEnded;
		firePropertyChange(FINISHED_PROPERTY, oldValue, processEnded);
	}

	/**
	 * @return the processEnded
	 * @see #processEnded
	 */
	public Boolean getProcessEnded() {
		return processEnded;
	}

	/**
	 * @return the progressVal
	 * @see #progressVal
	 */
	public Integer getProgressVal() {
		return progressVal;
	}

	/**
	 * @return the inpointsfile
	 * @see #pointsfile
	 */
	public String getInpointsfile() {
		return pointsfile;
	}

	/**
	 * @param inpointsfile
	 *          the inpointsfile to set
	 * @throws IOException
	 * @see #pointsfile
	 */
	public void setPointsfile(String inpointsfile) throws IOException {
		String oldValue = this.pointsfile;
		this.pointsfile = inpointsfile;
		firePropertyChange(POINTS_FILE_PROPERTY, oldValue, inpointsfile);
		if (parameters.getOutputFile() != null) {
			FileUtils.copyFile(new File(pointsfile), new File(parameters.getOutputFile()));
		}

	}

	/**
	 * @return the name
	 * @see #name
	 */
	public String getName() {
		return name;
	}

	/**
	 * @param name
	 *          the name to set
	 * @see #name
	 */
	public void setName(String name) {
		this.name = name;
	}

}
