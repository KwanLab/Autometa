#!/usr/bin/env nextflow

nextflow.enable.dsl=2


process SINGLE_RSYNC_DOWNLOAD {
  tag "Downloading: ${filename}"
  storeDir params.rsync_storedir
  cpus = 1

  input:
    val rsync_url //path to single file
    val filename

  output:
    path "${filename}"

  """
  rsync -a \
       --quiet \
       --checksum \
       ${rsync_url} \
       .
  """
}



process SINGLE_WGET_DOWNLOAD {
  tag "Downloading: ${filename}"
  storeDir params.rsync_storedir
  cpus = 1

  input:
    val wget_url //path to single file
    val wget_md5_url //path to single file's md5 
    val filename

  output:
    path "${filename}"

  """
  wget ${wget_url}
  # wget ${wget_md5_url} && md5sum -c --ignore-missing *.md5
  # 2nd wget commented until fix: https://github.com/KwanLab/Autometa/issues/149
  """
}
