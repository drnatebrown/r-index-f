# R-Index-F
<!--- ```console
  ____            ___           _                     _____ 
 |  _ \          |_ _|_ __   __| | _____  __         |  ___|
 | |_) |  _____   | || '_ \ / _` |/ _ \ \/ /  _____  | |_   
 |  _ <  |_____|  | || | | | (_| |  __/>  <  |_____| |  _|  
 |_| \_\         |___|_| |_|\__,_|\___/_/\_\         |_|    
                                                            
```
-->

R-Index-F Library for String Indexing

Implemented and adapted from original work by Takaaki Nishimoto and Yasuo Tabei [1].

This library uses a simplified approach which follows the theory of the original paper. We store intervals consisting of full BWT-runs rather than sub-runs, representing a maximal interval mapping, and a custom block compression [2]. Joint work with Travis Gagie and Massimiliano Rossi. To reproduce experiments shown in RLBWT Tricks [2], accession codes for SARS-CoV-2 genomes from the [Covid-19 Data Portal](https://www.covid19dataportal.org/) are listed in the example data.

Efficiently performs decompression and count queries using interval mapping of BWT-runs.

*Current Version:* 0.2.0

# Example
### Download and Compile

```console
git clone https://github.com/drnatebrown/r-index-f.git
cd r-index-f

mkdir build && cd build
cmake ..
make
```

### Build
Builds the data structure on the example fasta file given, creating [filename].rif as output. The -f flag specifies we read in a fasta format. Other build flags affect the BWT build and are described in [Big-BWT](https://github.com/alshai/Big-BWT.git).
```console
python3 rif ../data/example_fasta/example.fasta -f
```

If using row splitting from [r-permute](https://github.com/drnatebrown/r-permute) to bound LF to $O(1)$ and $O(r)$-space, use the `-d` option. First, copy the output of r-permute to rename using the split parameter.
```console
cp <FASTA>.d_col <FASTA>.<d>_col
python3 rif ../data/example_fasta/example.fasta -f -d <d>
```

### Queries
The data structure should be imported and loaded as decribed in r-index-f.hpp once built, and supports LF computation needed to perform count queries. An example command prints the count query for a pattern to stdout, assuming the table was built using default settings.
```console
./test/src/count_query
```


# External Dependencies

* [Big-BWT](https://github.com/alshai/Big-BWT.git)
    * [gSACA-K](https://github.com/felipelouza/gsa-is.git)
    * [malloc_count](https://github.com/bingmann/malloc_count)
* [pfp_thresholds](https://github.com/maxrossi91/pfp-thresholds)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [divufsort](https://github.com/simongog/libdivsufsort) 
* [Google Benchmark](https://github.com/google/benchmark.git)
    * [Google Test](https://github.com/google/googletest)

# Authors

### Implementation:

* [Nathaniel Brown](https://github.com/drnatebrown)
* [Massimiliano Rossi](https://github.com/maxrossi91)

### Theory
* Nathaniel Brown
* Travis Gagie
* Massimiliano Rossi

# Citation
Please cite the original paper by Nishimoto and Tabei [1] if you refer only to their data structure

If you use the implementation in an academic setting, or the style of our data structure, please cite both the former as well as RLBWT Tricks [2].
 
# References

[1] Nishimoto, T., & Tabei, Y. (2020). Optimal-Time Queries on BWT-runs Compressed Indexes. arXiv preprint arXiv:2006.05104.  
[2] Brown, N.K., Gagie, T., & Rossi, M. (2022). RLBWT Tricks. arXiv preprint arXiv:2112.04271.
