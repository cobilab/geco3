<p align="left">
     <img src="https://github.com/cobilab/geco3/blob/master/geco3.png?raw=true" alt="GeCo3" height="140"/>
</p>
<a href="https://anaconda.org/bioconda/geco3"> <img src="https://anaconda.org/bioconda/geco3/badges/downloads.svg" /> </a>
<a href="https://conda.anaconda.org/bioconda"> <img src="https://anaconda.org/bioconda/geco3/badges/installer/conda.svg" /> </a>
<a href="https://anaconda.org/bioconda/geco3"> <img src="https://anaconda.org/bioconda/geco3/badges/platforms.svg" /> </a>

GeCo3 is a DNA compressor that uses a neural network to do the mixing of experts.

Installation
---
Installation is available with conda.
```
conda install -y -c bioconda geco3
```

Build
----
GeCo3 depends on `make` and `gcc`.

```
cd src
make
```
**NOTE:** The default compilation options use micro-architecture dependent instructions, because GeCo3 greatly benefits from vector instructions. This means that due to the different floating point accuracies (e.g. fused multiply–add) the compressed file might not decompress with a binary compressed in a different computer or with a different compiler version or options. To ensure the file decompresses, use binaries (GeCo3 and GeDe3) that were compiled in the same environment.

Example usage
----
Reference free:
```
# Compression of file BuEb, using level 1, learning rate 0.06 and 8 hidden nodes
./GeCo3 -l 1 -lr 0.06 -hs 8 BuEb

# Decompression
./GeDe3 BuEd.co
```

Referential:
```
# Compression of file GG_C4 using as reference file PT_C4.
./GeCo3 -rm 20:500:1:35:0.95/3:100:0.95 -rm 13:200:1:1:0.95/0:0:0 -rm 10:10:0:0:0.95/0:0:0 -lr 0.03 -hs 64 -r PT_C4 GG_C4

# Decompression
./GeDe3 -r PT_C4 GG_C4.co
```

A complete description of the parameters can be read by invoking:
```
./GeCo3 -h
./GeDe3 -h
```

Export mixer to other compressors
----
Check the instructions and code in: https://github.com/cobilab/nn-expert-mixer

Citation
----
If you use GeCo3, please cite:
* Milton Silva, Diogo Pratas, Armando J Pinho, **"Efficient DNA sequence compression with neural networks"**, GigaScience, Volume 9, Issue 11, November 2020, giaa119, https://doi.org/10.1093/gigascience/giaa119
