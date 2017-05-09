package lu.uni.lcsb.vizbin;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
// import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
// import javax.swing.JPanel;
import javax.swing.JProgressBar;

import lu.uni.lcsb.vizbin.data.Sequence;

/**
 * 
 * 
 * @author <a href="mailto:valentin.plugaru.001@student.uni.lu">Valentin
 *         Plugaru</a>
 */
public final class DataExporter {

	/**
	 * Default constructor for utility class. Prevents instatiation.
	 */
	private DataExporter() {

	}

	/**
	 * Export list of sequences into a file.
	 *  
	 * @param parentFrame
	 * @param inFileName
	 * @param sequenceList
	 * @return
	 */
	public static File exportCluster(JFrame parentFrame, String inFileName, List<Sequence> sequenceList) {
		File outFile = null;
		File inFile = new File(inFileName);
		int pos;
		Thread progressThread;
		final JDialog dlg = new JDialog(parentFrame, "Sequence export progress", true);

		if (!sequenceList.isEmpty()) { // only try to export if there's something to
																		// export
			outFile = getSelectedFile();
			if (outFile != null) { // only try to export if the user chose an output
															// file
				BufferedReader br = null;
				BufferedWriter bw = null;
				try {
					// draw a progress dialog
					JProgressBar dpb = new JProgressBar(0, sequenceList.size());
					dlg.add(BorderLayout.CENTER, dpb);
					dlg.add(BorderLayout.NORTH, new JLabel("Progress..."));
					dlg.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
					dlg.setSize(300, 75);
					dlg.setLocationRelativeTo(parentFrame);
					progressThread = new Thread(new Runnable() {
						public void run() {
							dlg.setVisible(true);
						}
					});
					progressThread.start();

					br = new BufferedReader(new InputStreamReader(new FileInputStream(inFile)));
					bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(outFile)));
					String line = br.readLine();
					int seqCounter = -1;
					int exportedCount = 0;
					boolean addHeaderLine;

					int[] lookupTable = new int[sequenceList.size()];
					for (int i = 0; i < sequenceList.size(); i++) {
						// for some reason, the sequence ID is in the 'Seq##' format
						// so we need to retrieve only the ##
						lookupTable[i] = Integer.parseInt(sequenceList.get(i).getId().split("Seq")[1]);
					}

					while (line != null) {
						if (line.startsWith(">")) { // if we've encountered a new sequence
																				// start
							seqCounter += 1;
							// pos = seqPosition(sequenceList, seqCounter);
							pos = Arrays.binarySearch(lookupTable, seqCounter);
							if (pos >= 0) { // if this sequence is in the list to export
								System.out.println("Exporting sequence: " + sequenceList.get(pos).getId() + " " + sequenceList.get(pos).getLabelName());
								addHeaderLine = true;
								while (line != null && (addHeaderLine || !line.startsWith(">"))) {
									addHeaderLine = false;
									bw.write(line + System.getProperty("line.separator"));
									line = br.readLine();
								}
								dpb.setValue(++exportedCount);
							} else {
								line = br.readLine();
							}
						} else {
							line = br.readLine();
						}
					}
				} catch (IOException e) {
					e.printStackTrace();
				} finally {
					try {
						br.close();
						bw.close();
						dlg.setVisible(false);
						dlg.dispose();
					} catch (IOException e) {
						e.printStackTrace();
					}
				}
			}
		}
		return outFile;
	}

	private static File getSelectedFile() {
		final JFileChooser fc = new JFileChooser();
		fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		fc.setMultiSelectionEnabled(false);
		int returnVal = fc.showSaveDialog(null);
		if (fc.getSelectedFile() == null) {
			return null;
		} else {
			if (fc.getSelectedFile().exists()) {
				returnVal = JOptionPane.showConfirmDialog(
						null, "The file you specified already exists, do you wish to overwrite it?", "Confirm overwrite", JOptionPane.YES_NO_OPTION);
				if (returnVal != JOptionPane.YES_OPTION) {
					return null;
				} else {
					return fc.getSelectedFile();
				}
			} else {
				return fc.getSelectedFile();
			}
		}
	}

}
