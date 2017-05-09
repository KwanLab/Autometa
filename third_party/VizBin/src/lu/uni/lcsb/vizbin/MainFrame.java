package lu.uni.lcsb.vizbin;

import java.awt.Container;
import java.awt.Desktop;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.InputEvent;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;
import javax.swing.DefaultComboBoxModel;
import javax.swing.DefaultListModel;
import javax.swing.GroupLayout;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu.Separator;
import javax.swing.KeyStroke;
import javax.swing.filechooser.FileNameExtensionFilter;

import lu.uni.lcsb.vizbin.clustering.ClusterPanel;
import lu.uni.lcsb.vizbin.pca.PcaType;

import org.apache.log4j.Appender;
import org.apache.log4j.FileAppender;
import org.apache.log4j.Logger;

/**
 * 
 * @author <a href="mailto:valentin.plugaru.001@student.uni.lu">Valentin
 *         Plugaru</a>
 */
public class MainFrame extends javax.swing.JFrame {
	private static final FileNameExtensionFilter	ZIP_FILTER							= new FileNameExtensionFilter("zip files", "zip");
	private static final FileNameExtensionFilter	TEXT_FILTER							= new FileNameExtensionFilter("text files", "txt");

	/**
	 * 
	 */
	private static final long											serialVersionUID				= 1L;

	private static final String										SAVEABLE_PROPERTY				= "custom:saveable";
	private static final String										KMER_DATA_FILE_PROPERTY	= "custom:k-merDataFile";
	private static final String										WORKSPACE_NAME_PROPERTY	= "custom:workspaceName";

	/**
	 * Property defining if the project can be saved or not.
	 */
	private Boolean																saveable								= false;

	/**
	 * Property describing where the k-mer frequencies should be saved (for debug
	 * purposes). If null then this data won't be saved.
	 */
	private String																kmerDataFile						= null;

	/**
	 * Property describing the name of workspace.
	 */
	private String																workspaceName						= null;

	private Boolean																defLog									= true;
	private Boolean																moreOpionsVisible				= false;
	private Settings															settings								= null;

	/**
	 * Data file from which results where computed.
	 */
	private String																indatafile							= null;
	/**
	 * Labels file from which results where computed.
	 */
	private String																inlabelsfile						= null;
	/**
	 * Input file from which results where computed.
	 */
	private String																inpointsfile						= null;

	private File																	lastOpenPath						= null;

	private ProcessInput													processor								= null;

	private PcaType																pcaType									= PcaType.EJML;

	/**
	 * Default class logger.
	 */
	private final Logger													logger									= Logger.getLogger(MainFrame.class.getName());

	@SuppressWarnings({ "unchecked", "rawtypes" })
	public MainFrame() {
		logger.debug("Init of Main application frame");

		PropertyChangeListener propertyChangeLogger = new PropertyChangeListener() {
			@Override
			public void propertyChange(final PropertyChangeEvent arg0) {
				logger.debug("Property changed: " + arg0.getPropertyName() + ". Old: " + arg0.getOldValue() + " New: " + arg0.getNewValue());
			}

		};
		addPropertyChangeListener(propertyChangeLogger);

		initComponents();
		// set default values to the GUI
		this.formatfieldContigLen.setText(Integer.toString(Config.DEFAULT_CONTIG_LENGTH));
		this.formatfieldNumThreads.setText(Integer.toString(Config.DEFAULT_THREAD_NUM));
		this.formatfieldKmer.setText(Integer.toString(Config.DEFAULT_KMER_LENGTH));
		this.formatfieldPca.setText(Integer.toString(Config.DEFAULT_PCA_COLUMNS));
		this.formatfieldTheta.setText(Double.toString(Config.DEFAULT_THETA));
		this.formatfieldPerplexity.setText(Double.toString(Config.DEFAULT_PERPLEXILITY));
		this.formatfieldSeed.setText(Integer.toString(Config.DEFAULT_SEED));
		this.comboboxMerge.setModel(new DefaultComboBoxModel(new String[] { "Yes", "No" }));
		this.comboboxLog.setModel(new DefaultComboBoxModel(new String[] { "Yes", "No" }));
		this.comboboxPca.setModel(new DefaultComboBoxModel(new String[] { PcaType.MTJ.getName(), PcaType.EJML.getName() }));
		this.comboboxPca.setSelectedIndex(0);
		if (Config.DEFAULT_MERGE) {
			this.comboboxMerge.setSelectedIndex(0);
		} else {
			this.comboboxMerge.setSelectedIndex(1);
		}

		if (defLog) {
			this.comboboxLog.setSelectedIndex(0);
		} else {
			this.comboboxLog.setSelectedIndex(1);
		}

		// for easy debugging, pre-set input file selector:
		// this.textfield_file.setText("/Users/cedric.laczny/Documents/phd/projects/BINNING/publication/VizBin_-_Application_Note/data/EssentialGenes.fa");

		// this.textfield_file.setText("/Users/cedric.laczny/Documents/phd/projects/BINNING/publication/VizBin_-_Application_Note/revision_01/data/DaVis_testdat.fa");
		// this.textfield_points_file.setText("/Users/cedric.laczny/Documents/phd/projects/BINNING/publication/VizBin_-_Application_Note/revision_01/data/DaVis_testdat.points.txt");
		// this.textfield_labels.setText("/Users/cedric.laczny/Documents/phd/projects/BINNING/publication/VizBin_-_Application_Note/revision_01/data/DaVis_testdat.loglength.ann");
	}

	void setSettings(Settings settings) {
		this.settings = settings;
	}

	@SuppressWarnings({ "rawtypes" })
	// <editor-fold defaultstate="collapsed"
	// desc="Generated Code">//GEN-BEGIN:initComponents
	private void initComponents() {
		java.awt.GridBagConstraints gridBagConstraints;

		jMenuItem1 = new javax.swing.JMenuItem();
		tabpanel = new javax.swing.JTabbedPane();
		tabPanelMain = new javax.swing.JPanel();
		labelFile = new javax.swing.JLabel();
		labelContigLen = new javax.swing.JLabel();
		labelNumThreads = new javax.swing.JLabel();
		labelPointsFile = new javax.swing.JLabel();
		labelLabels = new javax.swing.JLabel();
		labelKmer = new javax.swing.JLabel();
		labelMerge = new javax.swing.JLabel();
		labelPcaDimensions = new javax.swing.JLabel();
		labelTheta = new javax.swing.JLabel();
		labelPerplexity = new javax.swing.JLabel();
		labelSeed = new javax.swing.JLabel();
		labelPcaLibrary = new javax.swing.JLabel();
		labelLog = new javax.swing.JLabel();
		textfieldFile = new javax.swing.JTextField();
		textfieldPointsFile = new javax.swing.JTextField();
		textfieldLabels = new javax.swing.JTextField();
		formatfieldKmer = new javax.swing.JFormattedTextField();
		comboboxMerge = new javax.swing.JComboBox();
		comboboxLog = new javax.swing.JComboBox();
		comboboxPca = new javax.swing.JComboBox();
		formatfieldContigLen = new javax.swing.JFormattedTextField();
		formatfieldNumThreads = new javax.swing.JFormattedTextField();
		formatfieldPca = new javax.swing.JFormattedTextField();
		formatfieldTheta = new javax.swing.JFormattedTextField();
		formatfieldPerplexity = new javax.swing.JFormattedTextField();
		formatfieldSeed = new javax.swing.JFormattedTextField();
		buttonFile = new javax.swing.JButton();
		buttonPointsFile = new javax.swing.JButton();
		buttonLabels = new javax.swing.JButton();
		buttonMoreOptions = new javax.swing.JButton();
		buttonProcess = new javax.swing.JButton();
		tabPanelVis = new javax.swing.JPanel();
		tabPanelPlugins = new javax.swing.JPanel();
		jScrollPane1 = new javax.swing.JScrollPane();
		pluginList = new javax.swing.JList();
		buttonReloadPlugins = new javax.swing.JButton();
		buttonLoadPlugin = new javax.swing.JButton();
		filler1 = new javax.swing.Box.Filler(new java.awt.Dimension(0, 0), new java.awt.Dimension(0, 0), new java.awt.Dimension(32767, 0));
		panelPluginOptions = new javax.swing.JPanel();
		labelPluginList = new javax.swing.JLabel();
		labelPluginOpts = new javax.swing.JLabel();
		statuspanel = new javax.swing.JPanel();
		labelStatus = new javax.swing.JLabel();
		progBar = new javax.swing.JProgressBar();

		jMenuItem1.setText("jMenuItem1");

		setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
		getContentPane().setLayout(new java.awt.GridBagLayout());

		tabPanelMain.setName("tab_panel_main"); // NOI18N
		tabPanelMain.setLayout(new java.awt.GridBagLayout());

		labelFile.setText("File to visualise:");
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 0;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelFile, gridBagConstraints);

		labelContigLen.setText("Minimal contig length:");
		labelContigLen.setVisible(true);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 1;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelContigLen, gridBagConstraints);

		labelNumThreads.setText("Number of threads:");
		labelNumThreads.setVisible(true);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 2;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelNumThreads, gridBagConstraints);

		labelPointsFile.setText("Point file (optional):");
		labelPointsFile.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 4;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelPointsFile, gridBagConstraints);

		labelLabels.setText("Annotation file (optional):");
		labelLabels.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 5;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelLabels, gridBagConstraints);

		labelKmer.setText("Kmer length:");
		labelKmer.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 6;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelKmer, gridBagConstraints);

		labelMerge.setText("Merge rev compl:");
		labelMerge.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 7;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelMerge, gridBagConstraints);

		labelPcaDimensions.setText("Intermediate dimensions:");
		labelPcaDimensions.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 8;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelPcaDimensions, gridBagConstraints);

		labelTheta.setText("Theta:");
		labelTheta.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 9;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelTheta, gridBagConstraints);

		labelPerplexity.setText("Perplexity:");
		labelPerplexity.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 10;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelPerplexity, gridBagConstraints);

		labelSeed.setText("Seed:");
		labelSeed.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 11;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelSeed, gridBagConstraints);

		labelPcaLibrary.setText("PCA library:");
		labelPcaLibrary.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 12;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelPcaLibrary, gridBagConstraints);

		textfieldFile.setMinimumSize(new java.awt.Dimension(4, 7));
		textfieldFile.setPreferredSize(new java.awt.Dimension(140, 15));
		textfieldFile.addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent evt) {
				textfieldFileMouseClicked(evt);
			}
		});
		labelLog.setText("Take logarithm of coverage & length?");
		labelLog.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 13;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.EAST;
		tabPanelMain.add(labelLog, gridBagConstraints);

		textfieldFile.setMinimumSize(new java.awt.Dimension(4, 7));
		textfieldFile.setPreferredSize(new java.awt.Dimension(140, 15));
		textfieldFile.addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent evt) {
				textfieldFileMouseClicked(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 0;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.ipadx = 136;
		gridBagConstraints.ipady = 12;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(textfieldFile, gridBagConstraints);

		formatfieldContigLen.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter()));
		formatfieldContigLen.setPreferredSize(new java.awt.Dimension(140, 29));
		formatfieldContigLen.setVisible(true);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 1;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(formatfieldContigLen, gridBagConstraints);

		formatfieldNumThreads.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter()));
		formatfieldNumThreads.setPreferredSize(new java.awt.Dimension(140, 29));
		formatfieldNumThreads.setVisible(true);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 2;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(formatfieldNumThreads, gridBagConstraints);

		buttonMoreOptions.setText("Show additional options");
		buttonMoreOptions.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				buttonMoreOptionsActionPerformed(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 3;
		tabPanelMain.add(buttonMoreOptions, gridBagConstraints);

		textfieldPointsFile.setVisible(false);
		textfieldPointsFile.setMinimumSize(new java.awt.Dimension(4, 7));
		textfieldPointsFile.setPreferredSize(new java.awt.Dimension(140, 15));
		textfieldPointsFile.addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent evt) {
				textfieldPointsFileMouseClicked(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 4;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.ipadx = 136;
		gridBagConstraints.ipady = 12;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(textfieldPointsFile, gridBagConstraints);

		textfieldLabels.setVisible(false);
		textfieldLabels.setMinimumSize(new java.awt.Dimension(4, 7));
		textfieldLabels.setPreferredSize(new java.awt.Dimension(140, 15));
		textfieldLabels.addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent evt) {
				textfieldLabelsMouseClicked(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 5;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.ipadx = 136;
		gridBagConstraints.ipady = 12;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(textfieldLabels, gridBagConstraints);

		formatfieldKmer.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter()));
		formatfieldKmer.setPreferredSize(new java.awt.Dimension(140, 29));
		formatfieldKmer.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 6;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(formatfieldKmer, gridBagConstraints);

		comboboxMerge.setPreferredSize(new java.awt.Dimension(140, 24));
		comboboxMerge.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 7;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.ipadx = 69;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(comboboxMerge, gridBagConstraints);

		formatfieldPca.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter()));
		formatfieldPca.setPreferredSize(new java.awt.Dimension(140, 29));
		formatfieldPca.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 8;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(formatfieldPca, gridBagConstraints);

		formatfieldTheta.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter()));
		formatfieldTheta.setPreferredSize(new java.awt.Dimension(140, 29));
		formatfieldTheta.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 9;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(formatfieldTheta, gridBagConstraints);

		formatfieldPerplexity.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter()));
		formatfieldPerplexity.setPreferredSize(new java.awt.Dimension(140, 29));
		formatfieldPerplexity.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 10;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(formatfieldPerplexity, gridBagConstraints);

		formatfieldSeed.setFormatterFactory(new javax.swing.text.DefaultFormatterFactory(new javax.swing.text.NumberFormatter()));
		formatfieldSeed.setPreferredSize(new java.awt.Dimension(140, 29));
		formatfieldSeed.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 11;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(formatfieldSeed, gridBagConstraints);

		comboboxPca.setPreferredSize(new java.awt.Dimension(140, 24));
		comboboxPca.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 12;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.ipadx = 69;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(comboboxPca, gridBagConstraints);

		comboboxLog.setPreferredSize(new java.awt.Dimension(140, 24));
		comboboxLog.setVisible(false);
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 13;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.ipadx = 69;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelMain.add(comboboxLog, gridBagConstraints);

		buttonFile.setText("Choose...");
		buttonFile.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				buttonFileActionPerformed(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 2;
		gridBagConstraints.gridy = 0;
		tabPanelMain.add(buttonFile, gridBagConstraints);

		buttonPointsFile.setVisible(false);
		buttonPointsFile.setText("Choose...");
		buttonPointsFile.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				buttonPointsFileActionPerformed(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 2;
		gridBagConstraints.gridy = 4;
		tabPanelMain.add(buttonPointsFile, gridBagConstraints);

		buttonLabels.setVisible(false);
		buttonLabels.setText("Choose...");
		buttonLabels.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				buttonLabelsActionPerformed(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 2;
		gridBagConstraints.gridy = 5;
		tabPanelMain.add(buttonLabels, gridBagConstraints);

		buttonProcess.setText("Start");
		buttonProcess.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				buttonProcessActionPerformed(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 2;
		gridBagConstraints.gridy = 7;
		tabPanelMain.add(buttonProcess, gridBagConstraints);

		tabpanel.addTab("Main", tabPanelMain);
		tabPanelMain.getAccessibleContext().setAccessibleName("Main");

		javax.swing.GroupLayout tabPanelVisLayout = new javax.swing.GroupLayout(tabPanelVis);
		tabPanelVis.setLayout(tabPanelVisLayout);
		tabPanelVisLayout.setHorizontalGroup(tabPanelVisLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addGap(0, 625, Short.MAX_VALUE));
		tabPanelVisLayout.setVerticalGroup(tabPanelVisLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING).addGap(0, 300, Short.MAX_VALUE));

		tabpanel.addTab("Visualisation", tabPanelVis);
		tabPanelVis.getAccessibleContext().setAccessibleName("Visualisation");

		tabPanelPlugins.setLayout(new java.awt.GridBagLayout());

		pluginList.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
		jScrollPane1.setViewportView(pluginList);

		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 1;
		gridBagConstraints.gridwidth = 2;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		gridBagConstraints.weightx = 0.25;
		gridBagConstraints.weighty = 1.0;
		tabPanelPlugins.add(jScrollPane1, gridBagConstraints);

		buttonReloadPlugins.setText("Reload list");
		buttonReloadPlugins.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				buttonReloadPluginsActionPerformed(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 2;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHWEST;
		tabPanelPlugins.add(buttonReloadPlugins, gridBagConstraints);

		buttonLoadPlugin.setText("Load plugin");
		buttonLoadPlugin.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				buttonLoadPluginActionPerformed(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 1;
		gridBagConstraints.gridy = 2;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.NORTHEAST;
		tabPanelPlugins.add(buttonLoadPlugin, gridBagConstraints);
		tabPanelPlugins.add(filler1, new java.awt.GridBagConstraints());

		panelPluginOptions.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));

		GroupLayout panelPluginOptionsLayout = new GroupLayout(panelPluginOptions);
		panelPluginOptions.setLayout(panelPluginOptionsLayout);
		panelPluginOptionsLayout.setHorizontalGroup(panelPluginOptionsLayout.createParallelGroup(GroupLayout.Alignment.LEADING).addGap(0, 360, Short.MAX_VALUE));
		panelPluginOptionsLayout.setVerticalGroup(panelPluginOptionsLayout.createParallelGroup(GroupLayout.Alignment.LEADING).addGap(0, 256, Short.MAX_VALUE));

		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 2;
		gridBagConstraints.gridy = 1;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.weightx = 0.75;
		gridBagConstraints.weighty = 1.0;
		tabPanelPlugins.add(panelPluginOptions, gridBagConstraints);

		labelPluginList.setText("Plugin list:");
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 0;
		gridBagConstraints.gridwidth = 2;
		tabPanelPlugins.add(labelPluginList, gridBagConstraints);

		labelPluginOpts.setText("Plugin options:");
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 2;
		gridBagConstraints.gridy = 0;
		tabPanelPlugins.add(labelPluginOpts, gridBagConstraints);

		/*
		 * Commented out - plugins are not used for now
		 */
		// tabpanel.addTab("Plugins", tab_panel_plugins);
		// tab_panel_plugins.getAccessibleContext().setAccessibleName("Plugins");

		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 0;
		gridBagConstraints.gridwidth = 11;
		gridBagConstraints.gridheight = 11;
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.weightx = 1.0;
		gridBagConstraints.weighty = 1.0;
		getContentPane().add(tabpanel, gridBagConstraints);
		tabpanel.getAccessibleContext().setAccessibleName("Main");

		statuspanel.setBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)));
		statuspanel.setLayout(new java.awt.GridBagLayout());

		labelStatus.setToolTipText("Double click this status bar to open the application's log file.");
		labelStatus.setBorder(javax.swing.BorderFactory.createBevelBorder(javax.swing.border.BevelBorder.RAISED));
		labelStatus.setMinimumSize(new java.awt.Dimension(4, 20));
		labelStatus.setPreferredSize(new java.awt.Dimension(20, 20));
		labelStatus.addMouseListener(new java.awt.event.MouseAdapter() {
			public void mouseClicked(java.awt.event.MouseEvent evt) {
				labelStatusMouseClicked(evt);
			}
		});
		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.fill = java.awt.GridBagConstraints.BOTH;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.LINE_START;
		gridBagConstraints.weightx = 0.5;
		statuspanel.add(labelStatus, gridBagConstraints);
		statuspanel.add(progBar, new java.awt.GridBagConstraints());

		gridBagConstraints = new java.awt.GridBagConstraints();
		gridBagConstraints.gridx = 0;
		gridBagConstraints.gridy = 11;
		gridBagConstraints.gridwidth = 11;
		gridBagConstraints.fill = java.awt.GridBagConstraints.HORIZONTAL;
		gridBagConstraints.anchor = java.awt.GridBagConstraints.PAGE_END;
		getContentPane().add(statuspanel, gridBagConstraints);

		setJMenuBar(createMenu());

		pack();
	}

	/**
	 * This method creates menu for the main form.
	 * 
	 * @return menu
	 */
	protected JMenuBar createMenu() {

		final JMenuBar menu = new JMenuBar();
		JMenu menuFile = new JMenu();
		JMenuItem menuFileReinitialize = new JMenuItem();
		Separator menuFileSep1 = new Separator();
		JMenuItem menuFileExit = new JMenuItem();
		JMenu menuOptions = new JMenu();
		JMenuItem menuOptionsExtVis = new JMenuItem();
		// menu_options_drawaxes = new javax.swing.JCheckBoxMenuItem();
		JMenuItem menuOptionsPlotToPng = new JMenuItem();
		JMenu menuAbout = new JMenu();
		JMenuItem menuOptionsShowlog = new JMenuItem();

		menuFile.setMnemonic('F');
		menuFile.setText("File");

		menuFileReinitialize.setMnemonic('r');
		menuFileReinitialize.setText("Reinitialize settings");
		menuFileReinitialize.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				menuFileReinitializeActionPerformed(evt);
			}
		});
		menuFile.add(menuFileReinitialize);

		final JMenuItem saveProject = new JMenuItem();
		saveProject.setText("Save workspace");
		saveProject.setMnemonic('s');
		saveProject.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				saveProject();
			}
		});
		saveProject.setEnabled(saveable);
		saveProject.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_S, InputEvent.CTRL_MASK));
		addPropertyChangeListener(new PropertyChangeListener() {
			@Override
			public void propertyChange(PropertyChangeEvent evt) {
				if (SAVEABLE_PROPERTY.equals(evt.getPropertyName())) {
					if ((Boolean) evt.getNewValue()) {
						saveProject.setEnabled(true);
					} else {
						saveProject.setEnabled(false);
					}
				}
			}
		});

		menuFile.add(saveProject);

		final JMenuItem openProject = new JMenuItem();
		openProject.setText("Open workspace");
		openProject.setMnemonic('o');
		openProject.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				openProject();
			}
		});
		openProject.setEnabled(true);
		openProject.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_O, InputEvent.CTRL_MASK));

		menuFile.add(openProject);

		menuFile.add(menuFileSep1);

		menuFileExit.setMnemonic('x');
		menuFileExit.setText("Exit");
		menuFileExit.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				menuFileExitActionPerformed(evt);
			}
		});
		menuFile.add(menuFileExit);

		menu.add(menuFile);

		menuOptions.setMnemonic('o');
		menuOptions.setText("Options");

		menuOptionsExtVis.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_O, java.awt.event.InputEvent.CTRL_MASK));
		menuOptionsExtVis.setText("Create external visualisation window");
		menuOptionsExtVis.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				popOutVisWindow(evt);
			}
		});
		menuOptions.add(menuOptionsExtVis);

		// menu_options_drawaxes.setText("Draw plot axes");
		/*
		 * Commented - not use for now. Axes drawing needs to be fixed.
		 */
		// menu_options.add(menu_options_drawaxes);

		menuOptionsPlotToPng.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_E, java.awt.event.InputEvent.CTRL_MASK));
		menuOptionsPlotToPng.setText("Export plot to PNG");
		menuOptionsPlotToPng.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				exportVisToFile(evt);
			}
		});
		menuOptions.add(menuOptionsPlotToPng);

		JMenu debugMenu = new JMenu();
		debugMenu.setText("Debug");
		debugMenu.setMnemonic('d');
		final JCheckBoxMenuItem kmerDataMenu = new JCheckBoxMenuItem();
		kmerDataMenu.setText("K-mer data");
		kmerDataMenu.setMnemonic('k');
		kmerDataMenu.setSelected(false);
		kmerDataMenu.addActionListener(new ActionListener() {
			@Override
			public void actionPerformed(ActionEvent e) {
				if (kmerDataMenu.isSelected()) {
					File output = getSelectedOutputFile(TEXT_FILTER);
					if (output == null) {
						kmerDataMenu.setSelected(false);
						setKmerDataFile(null);
					} else {
						try {
							setKmerDataFile(output.getCanonicalPath());
						} catch (IOException e1) {
							logger.error(e1, e1);
							JOptionPane.showMessageDialog(menu.getComponent(), e1.getMessage(), "Problem with file", JOptionPane.ERROR_MESSAGE);
						}
					}
				} else {
					setKmerDataFile(null);
				}
			}
		});
		debugMenu.add(kmerDataMenu);

		menuOptions.add(debugMenu);

		menu.add(menuOptions);

		menuAbout.setMnemonic('A');
		menuAbout.setText("About");

		menuOptionsShowlog.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_L, java.awt.event.InputEvent.CTRL_MASK));
		menuOptionsShowlog.setText("Show application log");
		menuOptionsShowlog.addActionListener(new java.awt.event.ActionListener() {
			public void actionPerformed(java.awt.event.ActionEvent evt) {
				menuOptionsShowlogActionPerformed(evt);
			}
		});
		menuAbout.add(menuOptionsShowlog);

		menu.add(menuAbout);
		return menu;
	}

	protected void saveProject() {
		File output = getSelectedOutputFile(ZIP_FILTER);
		if (output != null && !output.equals("")) {
			Appender app = Logger.getRootLogger().getAppender("R");
			String logFile = null;
			if (app instanceof FileAppender && app != null) {
				logFile = ((FileAppender) app).getFile();
			} else {
				logger.error("Cannot find log appender");
			}
			ZipProject zp = new ZipProject(indatafile, inlabelsfile, inpointsfile, logFile);
			try {
				if (new File(output.getCanonicalPath()).exists()) {
					new File(output.getCanonicalPath()).delete();
				}
				zp.saveTo(output.getCanonicalPath());
				JOptionPane.showMessageDialog(this, "File saved successfully", "Info", JOptionPane.INFORMATION_MESSAGE);
			} catch (IOException e) {
				logger.error(e, e);
				JOptionPane.showMessageDialog(this, e.getMessage(), "Problem with file", JOptionPane.ERROR_MESSAGE);
			}
		}
	}

	protected void openProject() {
		String input = getSelectedInputFile(ZIP_FILTER);
		if (input != null && !input.equals("")) {
			ZipProject zip;
			try {
				zip = new ZipProject(input);
				this.textfieldFile.setText(zip.getDataInputFile());
				if (zip.getLabelInputFile() != null) {
					this.textfieldLabels.setText(zip.getLabelInputFile());
				}
				this.textfieldPointsFile.setText(zip.getPointInputFile());
				File f = new File(input);
				setWorkspaceName(f.getName());
			} catch (IOException e) {
				logger.error(e, e);
				JOptionPane.showMessageDialog(this, e.getMessage(), "Problem with file", JOptionPane.ERROR_MESSAGE);
			}
		}
	}

	private void buttonProcessActionPerformed(java.awt.event.ActionEvent evt) {
		Integer contigLen = Config.DEFAULT_CONTIG_LENGTH;
		Integer numThreads = Config.DEFAULT_THREAD_NUM;
		Integer kmer = Config.DEFAULT_KMER_LENGTH;
		Integer pca = Config.DEFAULT_PCA_COLUMNS;
		Double theta = Config.DEFAULT_THETA;
		Double perplexity = Config.DEFAULT_PERPLEXILITY;
		Integer seed = Config.DEFAULT_SEED;
		Boolean merge = Config.DEFAULT_MERGE;
		Boolean log = defLog;

		indatafile = this.textfieldFile.getText();
		if (indatafile.isEmpty()) {
			JOptionPane.showMessageDialog(null, "You must specify an input Fasta file!");
			return;
		}

		try {
			contigLen = Integer.parseInt(this.formatfieldContigLen.getText().replaceAll(",", ""));
			if (contigLen < 0) {
				throw new NumberFormatException(this.formatfieldContigLen.getText());
			}
		} catch (NumberFormatException e) {
			logger.warn("Invalid minimal contig length value: " + this.formatfieldContigLen.getText() + ". Using value: " + contigLen);
		}

		try {
			numThreads = Integer.parseInt(this.formatfieldNumThreads.getText());

		} catch (NumberFormatException e) {
			logger.warn("Invalid numThreads value: " + this.formatfieldNumThreads.getText().replaceAll(",", "") + ". Using value: " + numThreads);
		}

		inpointsfile = this.textfieldPointsFile.getText();
		inlabelsfile = this.textfieldLabels.getText();
		try {
			kmer = Integer.parseInt(this.formatfieldKmer.getText());
		} catch (NumberFormatException e) {
			logger.warn("Invalid kmer value: " + this.formatfieldKmer.getText() + ". Using value: " + kmer);
		}

		merge = ((String) this.comboboxMerge.getSelectedItem()).equals("Yes");

		log = ((String) this.comboboxLog.getSelectedItem()).equals("Yes");

		pcaType = null;
		for (PcaType type : PcaType.values()) {
			if (((String) this.comboboxPca.getSelectedItem()).equals(type.getName())) {
				pcaType = type;
			}
		}
		if (pcaType == null) {
			logger.warn("Invalid PCA type: " + this.comboboxPca.getSelectedItem() + ".");
			pcaType = PcaType.EJML;
		}

		try {
			pca = Integer.parseInt(this.formatfieldPca.getText());
		} catch (NumberFormatException e) {
			logger.warn("Invalid PCA value: " + this.formatfieldPca.getText() + ". Using value: " + pca);
		}

		try {
			theta = Double.parseDouble(this.formatfieldTheta.getText());
		} catch (NumberFormatException e) {
			logger.warn("Invalid theta value: " + this.formatfieldTheta.getText() + ". Using value: " + theta);
		}

		try {
			perplexity = Double.parseDouble(this.formatfieldPerplexity.getText());
		} catch (NumberFormatException e) {
			logger.warn("Invalid perplexity value: " + this.formatfieldPerplexity.getText() + ". Using value: " + perplexity);
		}

		try {
			seed = Integer.parseInt(this.formatfieldSeed.getText());
		} catch (NumberFormatException e) {
			logger.warn("Invalid seed value: " + this.formatfieldSeed.getText() + ". Using value: " + seed);
		}

		if (processor == null || processor.getProcessEnded() == true) {
			try {
				ProcessParameters params = new ProcessParameters().inputFastaFile(indatafile).//
				contigLength(contigLen).//
						threads(numThreads).//
						inputPointFile(inpointsfile).//
						inputLabelFile(inlabelsfile).//
						kMerLength(kmer).//
						merge(merge).//
						pcaColumns(pca).//
						theta(theta).//
						perplexity(perplexity).//
						seed(seed).//
						pcaAlgorithmType(pcaType).//
						kmerDebugFile(kmerDataFile).//
						extendedLogs(log);

				ProcessGuiParameters guiParams = new ProcessGuiParameters(this.labelStatus, this.progBar, this.tabpanel, this, false);
				processor = new ProcessInput(params, guiParams, settings.getBinFile());

				processor.setName(getWorkspaceName());

				processor.addPropertyChangeListener(new PropertyChangeListener() {
					@Override
					public void propertyChange(PropertyChangeEvent evt) {
						if (ProcessInput.FINISHED_PROPERTY.equals(evt.getPropertyName())) {
							boolean saveable = false;
							if ((Boolean) evt.getNewValue() && processor.getProgressVal() == 100) {
								saveable = true;
							}
							setSaveable(saveable);
						}
					}
				});
				processor.addPropertyChangeListener(new PropertyChangeListener() {
					@Override
					public void propertyChange(PropertyChangeEvent evt) {
						if (ProcessInput.POINTS_FILE_PROPERTY.equals(evt.getPropertyName())) {
							setInpointsfile((String) evt.getNewValue());
						}
					}
				});

				processor.doProcess();
			} catch (Exception e) {
				logger.error("Problem with initializing compuatations...", e);
			}
		}
	}

	private String getSelectedInputFile(FileNameExtensionFilter filter) {
		final JFileChooser fc = new JFileChooser();
		if (filter != null) {
			fc.setFileFilter(filter);
		}
		fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		fc.setMultiSelectionEnabled(false);
		if (lastOpenPath != null)
			fc.setCurrentDirectory(lastOpenPath);
		fc.showOpenDialog(this.getParent());
		if (fc.getSelectedFile() == null) {
			return "";
		} else {
			if (!fc.getSelectedFile().exists()) {
				JOptionPane.showMessageDialog(null, "The file you specified doesn't exist!");
				return "";
			} else {
				lastOpenPath = fc.getSelectedFile();
				return fc.getSelectedFile().toString();
			}
		}
	}

	private File getSelectedOutputFile(FileNameExtensionFilter filter) {
		final JFileChooser fc = new JFileChooser();
		if (filter != null) {
			fc.setFileFilter(filter);
		}
		fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
		fc.setMultiSelectionEnabled(false);
		if (lastOpenPath != null) {
			fc.setCurrentDirectory(lastOpenPath);
		}
		int returnVal = fc.showSaveDialog(null);
		if (fc.getSelectedFile() == null) {
			return null;
		} else {
			lastOpenPath = fc.getSelectedFile();
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

	private void buttonFileActionPerformed(java.awt.event.ActionEvent evt) {
		this.textfieldFile.setText(getSelectedInputFile(null));
		setSaveable(false);
		File f = new File(this.textfieldFile.getText());
		setWorkspaceName(f.getName());

	}

	private void textfieldFileMouseClicked(java.awt.event.MouseEvent evt) {
		this.textfieldFile.setText(getSelectedInputFile(null));
		setSaveable(false);
		File f = new File(this.textfieldFile.getText());
		setWorkspaceName(f.getName());
	}

	private void buttonPointsFileActionPerformed(java.awt.event.ActionEvent evt) {
		this.textfieldPointsFile.setText(getSelectedInputFile(null));
		setSaveable(false);
	}

	private void textfieldPointsFileMouseClicked(java.awt.event.MouseEvent evt) {
		this.textfieldPointsFile.setText(getSelectedInputFile(null));
		setSaveable(false);
	}

	private void buttonLabelsActionPerformed(java.awt.event.ActionEvent evt) {
		this.textfieldLabels.setText(getSelectedInputFile(null));
		setSaveable(false);
	}

	private void buttonMoreOptionsActionPerformed(ActionEvent evt) {
		if (moreOpionsVisible) {
			labelPointsFile.setVisible(false);
			labelLabels.setVisible(false);
			labelKmer.setVisible(false);
			labelMerge.setVisible(false);
			labelPcaDimensions.setVisible(false);
			labelTheta.setVisible(false);
			labelPerplexity.setVisible(false);
			labelSeed.setVisible(false);
			labelPcaLibrary.setVisible(false);
			labelLog.setVisible(false);
			textfieldLabels.setVisible(false);
			textfieldPointsFile.setVisible(false);
			formatfieldKmer.setVisible(false);
			formatfieldPca.setVisible(false);
			formatfieldPerplexity.setVisible(false);
			formatfieldSeed.setVisible(false);
			formatfieldTheta.setVisible(false);
			comboboxMerge.setVisible(false);
			comboboxPca.setVisible(false);
			comboboxLog.setVisible(false);
			buttonMoreOptions.setText("Show additional options");
			moreOpionsVisible = false;
			buttonLabels.setVisible(false);
			buttonPointsFile.setVisible(false);
			repaint();
		} else {
			labelPointsFile.setVisible(true);
			labelLabels.setVisible(true);
			labelKmer.setVisible(true);
			labelMerge.setVisible(true);
			labelPcaDimensions.setVisible(true);
			labelTheta.setVisible(true);
			labelPerplexity.setVisible(true);
			labelSeed.setVisible(true);
			labelPcaLibrary.setVisible(true);
			labelLog.setVisible(true);
			textfieldLabels.setVisible(true);
			textfieldPointsFile.setVisible(true);
			formatfieldKmer.setVisible(true);
			formatfieldPca.setVisible(true);
			formatfieldPerplexity.setVisible(true);
			formatfieldSeed.setVisible(true);
			formatfieldTheta.setVisible(true);
			comboboxMerge.setVisible(true);
			comboboxPca.setVisible(true);
			comboboxLog.setVisible(true);
			buttonMoreOptions.setText("Hide additional options");
			moreOpionsVisible = true;
			buttonLabels.setVisible(true);
			buttonPointsFile.setVisible(true);
			pack();
			repaint();
		}
	}

	private void textfieldLabelsMouseClicked(java.awt.event.MouseEvent evt) {
		this.textfieldLabels.setText(getSelectedInputFile(null));
		setSaveable(false);
	}

	private void menuFileExitActionPerformed(java.awt.event.ActionEvent evt) {
		System.exit(0);
	}

	private void menuFileReinitializeActionPerformed(java.awt.event.ActionEvent evt) {
		Integer option = JOptionPane.showConfirmDialog(this, "Are you sure you want to reinitialize the settings?\n"
				+ "This will recreate the config file and redeploy the OS-specific TSNE application.", "Confirm reinitialization", JOptionPane.YES_NO_OPTION);
		if (option == JOptionPane.YES_OPTION) {
			settings.createSettings();
			settings.extractTSNEBin();
		}
	}

	private void labelStatusMouseClicked(java.awt.event.MouseEvent evt) {
		if (evt.getClickCount() == 2) {
			openLogs();
		}
	}

	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void buttonReloadPluginsActionPerformed(java.awt.event.ActionEvent evt) {
		DefaultListModel pluginListModel = new DefaultListModel();
		ArrayList<String> plugins = settings.getPluginUtils().listPlugins();
		for (int i = 0; i < plugins.size(); i++) {
			pluginListModel.addElement(plugins.get(i));
		}
		pluginList.setModel(pluginListModel);
	}

	private void buttonLoadPluginActionPerformed(java.awt.event.ActionEvent evt) {
		PluginUtils.PLUGINSTATE pluginState;
		if (pluginList.getSelectedIndex() != -1) {
			settings.getPluginUtils().setParentOptionsPanel(panelPluginOptions);
			settings.getPluginUtils().loadPlugin((String) pluginList.getSelectedValue());
			pluginState = settings.getPluginUtils().getPluginState();
			if (pluginState == PluginUtils.PLUGINSTATE.CANNOTLOAD) {
				JOptionPane.showMessageDialog(this, "The plugin could not be loaded.");
			} else if (pluginState == PluginUtils.PLUGINSTATE.NOTVALID) {
				JOptionPane.showMessageDialog(this, "The specified plugin doesn't appear to be a valid plugin.");
			} else if (pluginState == PluginUtils.PLUGINSTATE.NOTSUPPORTED) {
				JOptionPane.showMessageDialog(this, "The plugin's version is not supported!");
			} else {
				// perform action if plugin is valid
				return;
			}
		}

	}

	private void popOutVisWindow(java.awt.event.ActionEvent evt) {
		if (DataSetUtils.isIsDataSetCreated()) {
			JFrame frame = new JFrame("Visualisation");
			ClusterPanel panel = new ClusterPanel(DataSetUtils.getDataSet(), indatafile, this);
			frame.getContentPane().add(panel.getChartPanel());
			frame.pack();
			frame.setVisible(true);
			frame.setSize(800, 600);
		} else {
			JOptionPane.showMessageDialog(null, "No data set available!\nPlease process a data file first.", "Missing dataset", JOptionPane.WARNING_MESSAGE);
		}
	}

	private void menuOptionsShowlogActionPerformed(java.awt.event.ActionEvent evt) {
		openLogs();
	}

	private void openLogs() {
		try {
			Desktop.getDesktop().open(new File(((FileAppender) Logger.getRootLogger().getAppender("R")).getFile()));
		} catch (IOException e) {
			logger.warn("Failed to open log file. Trying to edit it...");
			try {
				Desktop.getDesktop().edit(new File(((FileAppender) Logger.getRootLogger().getAppender("R")).getFile()));
			} catch (IOException e1) {
				logger.error(e1, e1);
			}
		}
	}

	private void exportVisToFile(java.awt.event.ActionEvent evt) {
		if (DataSetUtils.isIsDataSetCreated()) {
			Container contentPane = DataSetUtils.getDrawingFrame().getContentPane();
			BufferedImage bufferedImage = new BufferedImage(contentPane.getWidth(), contentPane.getHeight(), BufferedImage.TYPE_INT_RGB);
			Graphics2D graphics2DContext = bufferedImage.createGraphics();
			contentPane.printAll(graphics2DContext);
			graphics2DContext.dispose();
			try {
				File outFile = getSelectedOutputFile(null);
				if (outFile != null)
					ImageIO.write(bufferedImage, "png", outFile);
			} catch (IOException ex) {
				ex.printStackTrace();
			}
		} else {
			JOptionPane.showMessageDialog(null, "No data set available!\nPlease process a data file first.", "Missing dataset", JOptionPane.WARNING_MESSAGE);
		}
	}

	/**
	 * @param args
	 *          the command line arguments
	 */
	public static void main(String args[]) {
		/* Set the Nimbus look and feel */
		// <editor-fold defaultstate="collapsed"
		// desc=" Look and feel setting code (optional) ">
		/*
		 * If Nimbus (introduced in Java SE 6) is not available, stay with the
		 * default look and feel. For details see
		 * http://download.oracle.com/javase/tutorial/uiswing/lookandfeel/plaf.html
		 */
		try {
			for (javax.swing.UIManager.LookAndFeelInfo info : javax.swing.UIManager.getInstalledLookAndFeels()) {
				if ("Nimbus".equals(info.getName())) {
					javax.swing.UIManager.setLookAndFeel(info.getClassName());
					break;
				}
			}
		} catch (ClassNotFoundException ex) {
			java.util.logging.Logger.getLogger(MainFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
		} catch (InstantiationException ex) {
			java.util.logging.Logger.getLogger(MainFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
		} catch (IllegalAccessException ex) {
			java.util.logging.Logger.getLogger(MainFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
		} catch (javax.swing.UnsupportedLookAndFeelException ex) {
			java.util.logging.Logger.getLogger(MainFrame.class.getName()).log(java.util.logging.Level.SEVERE, null, ex);
		}
		// </editor-fold>

		/* Create and display the form */
		java.awt.EventQueue.invokeLater(new Runnable() {
			public void run() {
				new MainFrame().setVisible(true);
			}
		});

	}

	/**
	 * @return the saveable
	 * @see #saveable
	 */
	public Boolean getSaveable() {
		return saveable;
	}

	/**
	 * @param saveable
	 *          the saveable to set
	 * @see #saveable
	 */
	public void setSaveable(Boolean saveable) {
		Boolean oldValue = this.saveable;
		this.saveable = saveable;
		firePropertyChange(SAVEABLE_PROPERTY, oldValue, saveable);
	}

	/**
	 * @return the inpointsfile
	 * @see #inpointsfile
	 */
	public String getInpointsfile() {
		return inpointsfile;
	}

	/**
	 * @param inpointsfile
	 *          the inpointsfile to set
	 * @see #inpointsfile
	 */
	private void setInpointsfile(String inpointsfile) {
		this.inpointsfile = inpointsfile;
	}

	/**
	 * @return the kmerDataFile
	 * @see #kmerDataFile
	 */
	public String getKmerDataFile() {
		return kmerDataFile;
	}

	/**
	 * @param kmerDataFile
	 *          the kmerDataFile to set
	 * @see #kmerDataFile
	 */
	public void setKmerDataFile(String kmerDataFile) {
		String oldValue = this.kmerDataFile;
		this.kmerDataFile = kmerDataFile;
		firePropertyChange(KMER_DATA_FILE_PROPERTY, oldValue, kmerDataFile);
	}

	/**
	 * @return the workspaceName
	 * @see #workspaceName
	 */
	public String getWorkspaceName() {
		return workspaceName;
	}

	/**
	 * @param workspaceName
	 *          the workspaceName to set
	 * @see #workspaceName
	 */
	public void setWorkspaceName(String workspaceName) {
		String oldValue = this.workspaceName;
		this.workspaceName = workspaceName;
		firePropertyChange(WORKSPACE_NAME_PROPERTY, oldValue, workspaceName);
	}

	// Variables declaration - do not modify//GEN-BEGIN:variables
	private javax.swing.JButton							buttonFile;
	private javax.swing.JButton							buttonPointsFile;
	private javax.swing.JButton							buttonLabels;
	private javax.swing.JButton							buttonMoreOptions;
	private javax.swing.JButton							buttonLoadPlugin;
	private javax.swing.JButton							buttonProcess;
	private javax.swing.JButton							buttonReloadPlugins;
	@SuppressWarnings("rawtypes")
	private javax.swing.JComboBox						comboboxMerge;
	@SuppressWarnings("rawtypes")
	private javax.swing.JComboBox						comboboxLog;
	@SuppressWarnings("rawtypes")
	private javax.swing.JComboBox						comboboxPca;
	private javax.swing.Box.Filler					filler1;
	private javax.swing.JFormattedTextField	formatfieldContigLen;
	private javax.swing.JFormattedTextField	formatfieldNumThreads;
	private javax.swing.JFormattedTextField	formatfieldKmer;
	private javax.swing.JFormattedTextField	formatfieldPca;
	private javax.swing.JFormattedTextField	formatfieldPerplexity;
	private javax.swing.JFormattedTextField	formatfieldSeed;
	private javax.swing.JFormattedTextField	formatfieldTheta;
	private javax.swing.JMenuItem						jMenuItem1;
	private javax.swing.JScrollPane					jScrollPane1;
	private javax.swing.JLabel							labelContigLen;
	private javax.swing.JLabel							labelNumThreads;
	private javax.swing.JLabel							labelFile;
	private javax.swing.JLabel							labelPointsFile;
	private javax.swing.JLabel							labelKmer;
	private javax.swing.JLabel							labelLabels;
	private javax.swing.JLabel							labelMerge;
	private javax.swing.JLabel							labelPcaDimensions;
	private javax.swing.JLabel							labelPerplexity;
	private javax.swing.JLabel							labelSeed;
	private javax.swing.JLabel							labelPcaLibrary;
	private javax.swing.JLabel							labelLog;
	private javax.swing.JLabel							labelPluginList;
	private javax.swing.JLabel							labelPluginOpts;
	private javax.swing.JLabel							labelStatus;
	private javax.swing.JLabel							labelTheta;
	// private javax.swing.JCheckBoxMenuItem menu_options_drawaxes;
	private javax.swing.JPanel							panelPluginOptions;
	@SuppressWarnings("rawtypes")
	private javax.swing.JList								pluginList;
	private javax.swing.JProgressBar				progBar;
	private javax.swing.JPanel							statuspanel;
	private javax.swing.JPanel							tabPanelMain;
	private javax.swing.JPanel							tabPanelPlugins;
	private javax.swing.JPanel							tabPanelVis;
	private javax.swing.JTabbedPane					tabpanel;
	private javax.swing.JTextField					textfieldFile;
	private javax.swing.JTextField					textfieldPointsFile;
	private javax.swing.JTextField					textfieldLabels;

	// End of variables declaration//GEN-END:variables

}
