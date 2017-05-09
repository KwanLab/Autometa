package lu.uni.lcsb.vizbin;

// import java.io.FileInputStream;
// import java.io.FileNotFoundException;
// import java.io.IOException;
// import java.util.Properties;

import javax.swing.JOptionPane;

import org.apache.log4j.Logger;

import java.awt.GraphicsEnvironment;

// import org.apache.log4j.PropertyConfigurator;

/**
 *
 * @author <a href="mailto:valentin.plugaru.001@student.uni.lu">Valentin
 *         Plugaru</a>
 *
 */
public class Main {

	/**
	 * Default constructor for utility class. Prevents instatiation.
	 */
	private Main() {

	}

	/**
	 * Default class logger.
	 */
	private static Logger		logger		= null;
	private static Settings	settings	= null;

	public static void main(String[] args) {
		settings = new Settings();
		logger = Logger.getLogger(Main.class);

		try {
			CommandLineOptions clo = new CommandLineOptions(args);

			if (clo.isValid()) {
				logger.debug("Running command line...");
				if (!settings.settingsExist()) {
					settings.createSettings();
					settings.extractTSNEBin();
				} else {
					settings.loadSettings();
				}
				ProcessParameters params = clo.getParameters();
				ProcessInput process = new ProcessInput(params, null, settings.getBinFile());
				process.doProcess();

			} else {
				clo.printHelp();
				if(!GraphicsEnvironment.isHeadless()){
					MainFrame mframe = new MainFrame();
					mframe.setVisible(true);
					mframe.setSettings(settings);

					if (!settings.settingsExist()) {
						JOptionPane.showMessageDialog(mframe, "Application settings not found, they will be created for you now.\n" + "Configuration file: "
								+ settings.getSettingsFile().toString());
						settings.createSettings();
						settings.extractTSNEBin();
					} else {
						settings.loadSettings();
					}
				}
			}
		} catch (Exception e) {
			logger.error(e, e);
			e.printStackTrace();
		}
	}

}
