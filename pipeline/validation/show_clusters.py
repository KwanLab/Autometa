#!/usr/bin/env python2.7

from bokeh.plotting import figure, show
from bokeh.io import output_file, show
from bokeh.models import ColumnDataSource, HoverTool, BoxSelectTool, LassoSelectTool, TapTool
from bokeh.models import Spacer, NumeralTickFormatter, PrintfTickFormatter, Title, CustomJS, Button
from bokeh.palettes import d3, Spectral6
from bokeh.transform import linear_cmap, factor_cmap, CategoricalColorMapper
from bokeh.layouts import widgetbox, column, row, layout
from bokeh.models.widgets import RangeSlider, Select, Button, DataTable, DateFormatter, TableColumn
import pandas as pd
import numpy as np
import sys, colorsys, argparse
from os.path import dirname, join

# TODO: Run as bokeh server to visualize ML recruitment in Real Time

def restricted_float(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r is not in the brightness range. Brightness must be between 0 and 1. " % (x,))
    return (x)

#create color spectrum of resolution N and brightness I, return as list of decimal RGB value tuples
def generate_color_range(N, I):
    HSV_tuples = [ (x*1.0/N, 0.5, I) for x in range(N) ]
    RGB_tuples = map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples)
    conversion = []
    for RGB_tuple in RGB_tuples:
        conversion.append((int(RGB_tuple[0]*255), int(RGB_tuple[1]*255), int(RGB_tuple[2]*255)))
    hex_colors = [ rgb_to_hex(RGB_tuple) for RGB_tuple in conversion ]
    return hex_colors, conversion

# convert RGB tuple to hexadecimal code
def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

# convert hexadecimal to RGB tuple
def hex_to_dec(hex):
    red = ''.join(hex.strip('#')[0:2])
    green = ''.join(hex.strip('#')[2:4])
    blue = ''.join(hex.strip('#')[4:6])
    return (int(red, 16), int(green, 16), int(blue,16))

def contig_selection(source=None, source2=None):
    return CustomJS(args=dict(source=source, source2=source2),code="""
    var bk_obj = cb_obj;
    var selectedIndices = cb_obj.attributes.callback.attributes.args.source.attributes.selected["1d"].indices;
    var data = cb_obj.attributes.callback.attributes.args.source.attributes.data;
    var contigClusters = []
    var selectedContigs = []
    for (i = 0; i < selectedIndices.length; i++) {
        selectedContigs.push(data.contig[selectedIndices[i]]);
        contigClusters.push(data.cluster[selectedIndices[i]]);
    }
    source2.data = {
    "contig":selectedContigs,
    "cluster":contigClusters
    }
""")

parser = argparse.ArgumentParser(description="Builds a scatter plot of clusters from AutoMeta output",
    epilog="Output is generated as an html file. Open the file in your browser to view the clustering.")
parser.add_argument('-i', metavar='<autometa_table>', help='Table generated from AutoMeta output (ML_recruitment_output.tab)', required=True)
parser.add_argument('-o', metavar='<outfile>', help='prefix to the generated html file', default="autometa_clustering", required=False)
parser.add_argument('-b', metavar='<brightness>', help='Adjusts the color scheme hue.', default=.70, type=restricted_float)
parser.add_argument("-legend", action='store_true', help='Adds legend to the scatter plot', required=False)
parser.add_argument('-plot_title', metavar='<plot_title>', help='text to place as the title of the scatter plot', default="AutoMeta Clustered Contigs", required=False)

args = vars(parser.parse_args())
infile = args['i']
outfile = args['o']
brightness = args['b']
plot_title = args['plot_title']

outfile_name = outfile + '.html'
output_file(outfile_name, title=outfile + ' bokeh plot')
fh = open(infile, "r")
df = pd.read_table(fh, delimiter="\t")
fh.close()

source = ColumnDataSource(df, name="clusterTable")
contigsList = ColumnDataSource(data=dict())


crcolor, crRGBs = generate_color_range(len(df.cluster.unique()),brightness) # produce spectrum

cluster_plot = figure(title=None,
                      toolbar_location="right",
                      plot_height=550,
                      plot_width=550,
                      y_axis_location= "left",
                      x_axis_location= 'above',
                      logo=None
                     )

cluster_plot.xaxis.axis_label = "bh-tsne-x"
cluster_plot.yaxis.axis_label = "bh-tsne-y"
cluster_plot.yaxis.axis_label_text_align = "left"
cluster_plot.xaxis.axis_label_text_align = "right"
cluster_plot.add_layout(Title(text=plot_title, align="left"), "above")
cluster_plot.add_tools(HoverTool(tooltips=[("Contig", "@contig"),
                                           ("Cluster", "@cluster"),
                                           ("Coverage", "@cov %"),
                                           ("Length", "@length bp"),
                                           ("Protein Families", "@single_copy_PFAMs"),
                                           ("# Single Copies", "@num_single_copies"),
                                           ("GC%", "@gc %"),
                                           ("Kingdom", "@kingdom"),
                                           ("Phylum", "@phylum"),
                                           ("Class", "@class"),
                                           ("Order", "@order"),
                                           ("Family", "@family"),
                                           ("Genus", "@genus"),
                                           ("Species", "@species"),
                                           ("Tax ID", "@taxid"),
                                          ],
                                 formatters={'Coverage': 'printf',
                                             'Length': 'printf',
                                             'GC%': 'printf'
                                            }
                                ),
                       BoxSelectTool(callback=contig_selection(source=source, source2=contigsList)),
                       LassoSelectTool(callback=contig_selection(source=source, source2=contigsList)),
                       TapTool(),
                      )

if args['legend']:
    cluster_plot.scatter(x='bh_tsne_x',
                        y='bh_tsne_y',
                        source=source,
                        legend='cluster',
                        color='cluster',
                        fill_color=factor_cmap('cluster',
                                               palette=crcolor,
                                               factors=list(df.cluster.unique())))

    # TODO: Need to figure out how to place legend outside of plot while keeping same spatial layout with histograms
    cluster_plot.legend.location = (0,0)
    cluster_plot.legend.label_text_align = 'left'
    #from bokeh.models import Legend
    #legend = Legend(items=[(list(df.cluster.unique()), cluster_plot)], location=(0, -30))
    #cluster_plot.add_layout(legend, 'right')

else:
    cluster_plot.scatter(x='bh_tsne_x',
                        y='bh_tsne_y',
                        source=source,
                        color='cluster',
                        fill_color=factor_cmap('cluster',
                                               palette=crcolor,
                                               factors=list(df.cluster.unique())))



gc_hover = HoverTool(tooltips=[("Contig", "@contig"),
                               ("GC Percentage", "@gc %"),
                               ("Cluster", "@cluster"),
                               ("Length", "@length bp")
                              ],
                     formatters={'GC Percentage' : 'printf', "Length" : "printf"},
                     mode="hline",
                     line_policy='interp',
                     point_policy='follow_mouse')


max_gc = df.gc.max()*1.1
gc_plot = figure(title=None,
                 y_range=cluster_plot.y_range,
                 y_axis_location = None,
                 x_range = (0, max_gc),
                 plot_height=cluster_plot.plot_height,
                 plot_width=200,
                 toolbar_location=None
                )

gc_plot.segment(x0=0,
                y0='bh_tsne_y',
                x1='gc',
                y1='bh_tsne_y',
                line_width=1,
                line_cap = 'square',
                line_color = factor_cmap('cluster',
                                         factors = list(df.cluster.unique()),
                                         palette=crcolor),
                source = source,
                name = 'gc%')

gc_plot.add_layout(Title(text="GC% vs. Contig",align="center"), "above")
gc_plot.ygrid.grid_line_color = None
gc_plot.xaxis.major_label_orientation = np.pi/4
gc_plot.xaxis[0].formatter = PrintfTickFormatter(format="%f %%")



# Adding hover tool
length_hover = HoverTool(tooltips=[("Contig", "@contig"),
                                   ("Length", "@length bp"),
                                   ("Cluster", "@cluster"),
                                   ("GC%", "@gc")
                                   ],
                         formatters={"Length" : "printf", 'GC%' : 'printf'},
                         line_policy='interp',
                         point_policy='follow_mouse',
                         show_arrow=True,
                         mode='vline')


max_length = df.length.max()*1.1
length_plot = figure(
                    x_range=cluster_plot.x_range,
                    y_range=(0,max_length),
                    plot_height=200,
                    plot_width=cluster_plot.plot_width+20,
                    x_axis_location = None,
                    y_axis_location = "right",
                    toolbar_location=None
                   )


length_plot.segment(x0='bh_tsne_x',
                    y0=0,
                    x1='bh_tsne_x',
                    y1='length',
                    line_width=1,
                    line_cap = 'round',
                    line_color = factor_cmap('cluster',
                                             factors = list(df.cluster.unique()),
                                             palette=crcolor,
                                            ),
                    source = source,
                    name = 'length',
                    )


length_plot.add_layout(Title(text="Length vs. Contig",align="center"), "left")
length_plot.xgrid.grid_line_color = None
length_plot.yaxis.formatter = PrintfTickFormatter(format="%f bp")
length_plot.yaxis.major_label_orientation = np.pi/4

button = Button(label="Save selected contigs to list", button_type="success")
button.callback = CustomJS(args=dict(source=contigsList),
                           code=open(join(dirname(__file__), "download_contigs_list.js")).read())
cluster_layout = layout(children=[
    column(
        row(cluster_plot, gc_plot),
        row(length_plot, button))
    ])

show(cluster_layout)
