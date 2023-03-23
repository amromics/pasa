
Pasa
====

Illumina short-read assemblies for microbial genomes are normally suffered from fragmentation due to repetitive 
elements that require long reads or close reference genomes to resolve. While long-read sequencing data is not 
always available, reference-based methods are usually subject to the bias problem which can be severe for highly
diversified strains.To tackle this issue, we use the pangenome graph built from a large set of multiple related genomes to guide 
the scaffolding process of a bacterial strain genome assembly. We propose Pasa, a graph-based algorithm framework 
that can utilize the pangenome graph and the assembly graph information to improve scaffolding quality.

----------
Installing
----------

Requirements:

* Python 3.6 or greater
* numpy
* scipy
* numba
* Panta

---------------
How to use Pasa
---------------

An example of using PASA is provided in ...

For more realistic examples and Python scripts to reproduce the results
in our paper are available at GitHub: https://github.com/....

-------
License
-------

The Pasa package is 3-clause BSD licensed.

This code was tested on 
Python 3.6, 3.7; numpy version 1.19.2; scipy version 1.5.3; numba version 0.52.0 
