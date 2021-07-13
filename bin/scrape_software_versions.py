#!/usr/bin/env python
# MIT License
# Copyright (c) 2018 nf-core
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import print_function
import os

# TODO nf-core: Add additional regexes for new tools in process get_software_versions
regexes = {
    "autometa": ["v_pipeline.txt", r"(\S+)"],
    "Nextflow": ["v_nextflow.txt", r"(\S+)"],
}
results = OrderedDict()
results["autometa"] = '<span style="color:#999999;">N/A</span>'
results["Nextflow"] = '<span style="color:#999999;">N/A</span>'
results = {}
version_files = [x for x in os.listdir(".") if x.endswith(".version.txt")]
for version_file in version_files:


# Search each file using its regex
for k, v in regexes.items():
    try:
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))
    except IOError:
        results[k] = False
    software = version_file.replace(".version.txt", "")
    if software == "pipeline":
        software = "nf-core/autometa"

    with open(version_file) as fin:
        version = fin.read().strip()
    results[software] = version

# Dump to YAML
print(
    """
id: 'software_versions'
section_name: 'autometa Software Versions'
section_href: 'https://github.com/autometa'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
"""
)
for k, v in sorted(results.items()):
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v))
print("    </dl>")

# Write out regexes as csv file:
with open("software_versions.tsv", "w") as f:
    for k, v in sorted(results.items()):
        f.write("{}\t{}\n".format(k, v))
