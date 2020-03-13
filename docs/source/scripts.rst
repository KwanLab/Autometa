Scripts
=======================

Main
----

Main pipeline script calls::

    autometa.py </path/to/metagenome.config>

binning
-------

- recursive_dbscan.py

common
------

- coverage.py
- kmers.py
- mag.py
- markers.py
- metagenome.py
- utilities.py

external (within common):
^^^^^^^^^^^^^^^^^^^^^^^^^

    - bedtools.py
    - bowtie.py
    - diamond.py
    - hmmer.py
    - prodigal.py
    - samtools.py
    - work_queue.py

config
------

- default.config
- metagenome.config
- databases.py
- environ.py
- project.py
- user.py

datasets
--------

- simulated.py (To be written)
- synthetic.py (To be written)

taxonomy
--------

- lca.py
- majority_vote.py
- ncbi.py

validation
----------

- cluster_process.py
