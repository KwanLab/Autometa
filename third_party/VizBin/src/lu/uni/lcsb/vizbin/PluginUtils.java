package lu.uni.lcsb.vizbin;

import java.awt.geom.Point2D;
import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLClassLoader;
import java.util.ArrayList;

import javax.swing.JPanel;

import lu.uni.lcsb.vizbin.clustering.ClusteringPlugin;

/**
 * 
 * @author <a href="mailto:valentin.plugaru.001@student.uni.lu">Valentin
 *         Plugaru</a>
 */
public class PluginUtils {

	public final static int	PLUGINMODELVERSION	= 1;

	public enum PLUGINSTATE {
		VALID, CANNOTLOAD, NOTVALID, NOTSUPPORTED
	};

	private ClusteringPlugin	plugin										= null;
	private File							pluginPath;

	private JPanel						parentPluginOptionsPanel	= null;

	public PluginUtils(File pluginPath) {
		this.pluginPath = pluginPath;
	}

	public void setParentOptionsPanel(JPanel panel) {
		this.parentPluginOptionsPanel = panel;
	}

	public ArrayList<String> listPlugins() {
		ArrayList<String> pluginlist = new ArrayList<String>();
		File[] filelist = pluginPath.listFiles();
		for (int i = 0; i < filelist.length; i++) {
			if (filelist[i].getName().endsWith(".class")) {
				pluginlist.add(filelist[i].getName().replace(".class", ""));
			}
		}
		return pluginlist;
	}

	@SuppressWarnings("rawtypes")
	public void loadPlugin(String className) {
		URLClassLoader pluginClassLoader;
		Class pluginClass;
		try {
			URL url = null;
			try {
				url = new URL("file://" + pluginPath.getAbsolutePath() + System.getProperty("file.separator"));
			} catch (MalformedURLException ex) {
				ex.printStackTrace();
			}
			pluginClassLoader = new URLClassLoader(new URL[] { url });
			pluginClass = pluginClassLoader.loadClass(className);
			plugin = (ClusteringPlugin) pluginClass.newInstance();
			plugin.drawOptionsPanel(parentPluginOptionsPanel);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public PLUGINSTATE getPluginState() {
		int pluginModelVersion;
		if (plugin == null) {
			return PLUGINSTATE.CANNOTLOAD;
		} else {
			try {
				pluginModelVersion = plugin.getPluginModelVersion();
			} catch (Exception e) {
				e.printStackTrace();
				return PLUGINSTATE.NOTVALID;
			}

			if (pluginModelVersion == PluginUtils.PLUGINMODELVERSION) {
				return PLUGINSTATE.VALID;
			} else {
				return PLUGINSTATE.NOTSUPPORTED;
			}
		}
	}

	public ArrayList<ArrayList<Point2D>> getClusters(ArrayList<Point2D> pointList) {
		if (plugin == null)
			return null;
		else
			return ((ClusteringPlugin) plugin).getClusters(pointList);
	}
}