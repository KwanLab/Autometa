# Autometa Python API

## Running modules

Many of the Autometa modules may be run standalone.

Simply pass in the `-m` flag when calling a script to signify to python you are
running the script as an Autometa *module*.

I.e. `python -m autometa.common.kmers -h`

:::{note}
Autometa has many *entrypoints* available that are utilized by the {ref}`autometa-nextflow-workflow` and {ref}`autometa-bash-workflow`. If you have installed autometa,
all of these entrypoints will be available to you.

If you would like to get a better understanding of each entrypoint, we recommend reading the {ref}`step-by-step-tutorial` section.
:::

## Using Autometa's Python API

Autometa's classes and functions are available after installation.
To access these, do the same as importing any other python library.

### Examples

#### Samtools wrapper

To incorporate a call to `samtools sort` inside of your python code using the Autometa samtools wrapper.

```python
from autometa.common.external import samtools

# To see samtools.sort parameters try the commented command below:
# samtools.sort?

# Run samtools sort command in ipython interpreter
samtools.sort(sam="<path/to/alignment.sam>", out="<path/to/output/alignment.bam>", cpus=4)
```

#### Metagenome Description

Here is an example to easily assess your metagenome's characteristics using Autometa's Metagenome class

```python
from autometa.common.metagenome import Metagenome

# To see input parameters, instance attributes and methods
# Metagenome?

# Create a metagenome instance
mg = Metagenome(assembly="/path/to/metagenome.fasta")

# To see available methods (ignore any elements in the list with a double underscore)
dir(mg)

# Get pandas dataframe of metagenome details.
metagenome_df = mg.describe()

metagenome_df.to_csv("path/to/metagenome_description.tsv", sep='\t', index=True, header=True)
```

#### k-mer frequency counting, normalization, embedding

To quickly perform a k-mer frequency counting, normalization and embedding pipeline...

```python
from autometa.common import kmers

# Count kmers
counts = kmers.count(
    assembly="/path/to/metagenome.fasta",
    size=5
)

# Normalize kmers
norm_df = kmers.normalize(
    df=counts,
    method="ilr"
)

# Embed kmers
embed_df = kmers.embed(
    norm_df,
    pca_dimensions=50,
    embed_dimensions=3,
    method="densmap"
)
```
