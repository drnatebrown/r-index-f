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

*Current Version:* 0.1.1

# Example
### Download and Compile

```console
git clone https://github.com/drnatebrown/r-index-f.git
cd r-index-f

mkdir build && cd build
cmake ..
make
```

### Example

```console
[build script in progress]
```

# External Dependencies

* [Big-BWT](https://github.com/alshai/Big-BWT.git)
    * [gSACA-K](https://github.com/felipelouza/gsa-is.git)
    * [malloc_count](https://github.com/bingmann/malloc_count)
* [sdsl-lite](https://github.com/simongog/sdsl-lite)
    * [divufsort](https://github.com/simongog/libdivsufsort) 
* [Google Benchmark](https://github.com/google/benchmark.git)
    * [Google Test](https://github.com/google/googletest)

# Authors

### Implementation:

* [Nathaniel Brown](https://github.com/oma219)
* [Massimiliano Rossi](https://github.com/maxrossi91)

### Theory
* Nathaniel Brown
* Travis Gagie
* Massimiliano Rossi

# Citation
Please cite the original paper by Nishimoto and Tabei [1] if you refer only to their data structure

If you use the implementation in an academic setting, please cite both the former as well as RLBWT Tricks [2].
 
# References

[1] Nishimoto, T., & Tabei, Y. (2020). Optimal-Time Queries on BWT-runs Compressed Indexes. arXiv preprint arXiv:2006.05104.  
[2] Brown, N.K., Gagie, T., & Rossi, M. (2022). RLBWT Tricks. arXiv preprint arXiv:2112.04271.
