#!/usr/bin/env python

import pkg_resources

__version__ = pkg_resources.get_distribution("autometa").version
console_scripts = (
    pkg_resources.get_distribution("autometa").get_entry_map()["console_scripts"].keys()
)
