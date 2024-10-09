#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

rmarkdown::render(
  input=args[[1]],
  params=list(
    bins_path=args[[2]],
    assembly_to_locus_path=args[[2]],
    assembly_report_path=args[[3]],
    genus=FALSE
  ),
  knit_root_dir=getwd(),
  output_dir=getwd(),
  output_file="mock_data_report_by_assembly.html"
)
rmarkdown::render(
  input=args[[1]],
  params=list(
    bins_path= args[[2]],
    assembly_to_locus_path = args[[2]],
    assembly_report_path = args[[3]],
    genus=TRUE
  ),
  knit_root_dir=getwd(),
  output_dir=getwd(),
  output_file="mock_data_report_by_genus.html"
)
