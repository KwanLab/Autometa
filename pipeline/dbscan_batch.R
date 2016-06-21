# Program to run dbscan over a range of eps values
# USAGE: Rscript dbscan_batch.r <input contig table> <start eps> <end eps> [<eps step> (default = 0.1)]

library(dbscan)
library(docopt)

# args is a named list
args <- docopt('USAGE: dbscan_batch.R <input_contig_table> <start_eps> <end_eps> [eps_step]')

if ("eps_step" %in% names(args)) {
	args[["eps_step"]] = 0.1
}

args[["start_eps"]] = as.numeric(args[["start_eps"]])
args[["end_eps"]] = as.numeric(args[["end_eps"]])

# Load data from table
data = read.table(args[["input_contig_table"]], header=TRUE)

eps_values = seq(args[["start_eps"]], args[["end_eps"]], by=args[["eps_step"]])

for (i in eps_values) {
	# First determine what the output file will be called
	output_filename = paste(args[["input_contig_table"]], "eps", toString(i), sep="_")
	d <- data.frame(data$vizbin_x, data$vizbin_y)
	db <- dbscan(data, eps=i, minPts=3)
	output_table <- data.frame(data, db.cluster=db$cluster)
	write.table(output_table, file=output_filename, sep="\t", quote=FALSE, row.names=FALSE)
}
