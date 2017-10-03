![Logo]
# Lowest Common Ancestor (LCA)

The lowest common ancestor refers to the ORF most closely associated with other given ORFs parsed from a BLAST query that is lowest on the tree of life. LCA uses sparse table creation with a range minimum query algorithm to retrieve the lowest common ORF. The ancestor from which all descendants are related is returned. The LCA algorithm is deterministic and will return the LCA, even if the queried ORFs are highly divergent. As the divergence between ORF ancestry becomes more extreme, the LCA will near the root of the tree.

## Getting Started

LCA has been exclusively tested in Linux and MacOS X environments.

LCA is contained in the [Autometa repository][AutometaRepo].

### Prerequisites

python >= 2.7

install [pip](https://packaging.python.org/tutorials/installing-packages/ "python package pip homepage")

python module dependencies.
pip install these modules:

```bash
pip install numpy
pip install tqdm
pip install cython
```

### Installing

####First setup the lca functions

You will need to navigate to the directory in which you have placed lca_functions.pyx and setup_lca_functions.py

`$ python setup_lca_functions.py build_ext --inplace`

####If the setup worked properly
you will be able to access the LCA usage information

`$ python lca.py --help`

This will display:
```
usage: lca.py [-h] [-v] [-f bitscore filter] [-fail_info]
              {database_files,database_directory} ... BLAST output

Script to find the Lowest Common Ancestor for each ORF from a BLAST output table

positional arguments:
  BLAST output          Path to BLAST output file.

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbose         Indicate verbose progress reporting
  -f bitscore filter    Filter to parse percentage of top BLAST hits based on
                        bitscore.
  -fail_info            Writes out files with failure taxid/orf information

Method for Accession of Database Files:
  {database_files,database_directory}
                        Additional help may be provided by specifying
                        (database_files|database_directory) -h
    database_files      Accesses database files provided by individually
                        listing nodes.dmp, names.dmp and accession2taxid
    database_directory  Accesses database files from directory provided

The LCA analysis output will be directed to run_taxonomy.py

NOTE:
LCA analysis will produce best results when database files are up to date.
Database files can be automatically updated before performing LCA analysis by specifying:
"lca.py [-f] [-v] [-fail_info] database_directory <path to database directory> -update <BLAST output>"

Up to date versions of nodes.dmp and names.dmp may be found at:
ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

prot.accession2taxid is updated weekly and may be found at:
ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz
```
## Initial database setup

LCA has a built in feature of automatically downloading and preparing NCBI's database information.

If you wish to have LCA set up your database information, simply pass the update flag when performing your first LCA.

`$ python lca.py database_directory <path to where you want to keep the database files> -update <BLAST output file>`

LCA will navigate to the specified directory, download and extract the database files before performing LCA.

Additional help information may be found for the database selection by specifying help after selecting database_directory or database_files

`$ python lca.py database_directory -h` or `python lca.py database_files -h`

### Additional Features

LCA has additional optional arguments to provide greater flexibility

1. Bitscore Filter (-f)
1. Verbose  (-v)
1. Failure Information  (-fail_info)
___


**Bitscore filter** is the designated value for ORFs to parse from the BLAST results. The default filters 90% and above of the top BLAST bitscore for each ORF. For reduced resolution, you can specify lower values as decimals to parse more hits from the BLAST results. This may result in a higher final LCA for each ORF being queried.

_Example application of filter:_

`$ python lca.py database_directory <path to database directory> -f 0.4 <BLAST file>`

**Verbose** is the flag for printing LCA progress to the terminal. The default has verbose turned off. Simply specifying the `-v` flag will initiate printing to the terminal.

_Example application of verbose:_

`$ python lca.py database_directory <path to database directory> -v <BLAST file>`

**Failure Information** is a flag that writes out a file "output_filename_failed_orfs.lca" if any ORFs fail during the LCA analysis.

_Example application of failure information tracking:_
`$ python lca.py database_directory <path to database directory> -fail_info <BLAST file>`

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [numpy](http://www.numpy.org/) - Used for Sparse table and Range Minimum Query algorithm
* [tqdm](https://pypi.python.org/pypi/tqdm) - Used for progress reporting
* [cython](http://cython.org/) - Used to cythonize LCA functions

## Contributing

Please read [CONTRIBUTING.md](path to file for instructions for making pull requests) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

Could place something here if we plan on keeping updated

## Authors

* **Jason C. Kwan** - *Tree creation*
* **Evan R. Rees** - *RMQ and Sparse table algorithm as well as framework*

## License

This project is licensed under the >>NEED TO INSERT LICENSE TITLE HERE>> License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

##### Special Thanks to:

* _Miguel Pignatelli_
* _Ian J. Miller_


### Algorithms

Tree of Life Creation:
The tree of life is constructed from the node.dmp file from [NCBI's taxonomy database](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/ "NCBI's taxonomy database directory") located in the taxdump.tar.gz compressed file. Branches stemming from the root are constructed until the entire tree has been built. Paths between each node by are built by traversing the tree using an Eulerian tour method. During the Eulerian tour, features of each tax ID are stored for sparse table creation.
___

Sparse Table:
The depth of each taxonomic ID in relation to the rest of the tree is used to efficiently store the tree in memory for quick lookup of taxonomic information. This algorithm employs dynamic programming to assess each tax ID in a range of other tax IDs starting from a range of the tax ID and it's closest relative and increasing to a range from the tax ID and it's furthest relative.
___

Range Minimum Query (RMQ):
Following the generation of the sparse table, a list of lists populated by respective ORFs from the BLAST query is given to the RMQ algorithm to determine the LCA. The RMQ algorithm utilizes the generated tree of tax IDs, sparse table and features of each tax ID. i.e. depth and location within the tree of tax IDs. The RMQ algorithm will look at the ORFs in pairs reducing the LCA to a final lowest common ancestor. Upon receiving the ORF list input, the RMQ algorithm will look at ORF pairs, determine the tax IDs between the two and return the closest tax ID to the root. Each ORF pair has an array of tax IDs linking the relation between the two. The array of tax IDs between the two ORFs is investigated for a lowest common ancestor. An LCA is returned and subsequent RMQ is performed between the returned LCA and the next ORF until a final LCA is returned. As more divergent ORFs are introduced the LCA will get higher until the lowest common ancestor is the root.



[AutometaRepo]:(place autometa repository link here "Autometa's Repository")
[Logo]:(logo "Kwan Lab LCA")
