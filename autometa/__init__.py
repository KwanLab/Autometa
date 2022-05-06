from importlib_metadata import entry_points
import pkg_resources

dist = pkg_resources.get_distribution("autometa")
__version__ = dist.version
console_scripts = dist.get_entry_map()["console_scripts"].keys()
