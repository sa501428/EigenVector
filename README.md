**POSSUM - PCA Of Sparse, SUper Massive Matrices**

R and C/C++ functions to compute a few leading eigenvectors of the correlation matrices of large sparse matrices. It is mostly aimed at HiC contact matrices but can also be applied to general matrices. Currently we assume that the matrices are symmetric but this condition may be lifted in the future.
  
  
**References**  
* Fine-mapping of nuclear compartments using ultra-deep Hi-C shows that active promoter and enhancer elements localize in the active A compartment even when adjacent sequences do not (https://www.biorxiv.org/content/10.1101/2021.10.03.462599v2)


See [C++ build, tests, and resolution-independent A/B orientation](C++/README.md)
for Power Method and Lanczos usage, including automatic bundled GC/RefSeq phasing for hg19, hg38, mm9, and mm10,
assembly overrides, confidence diagnostics, and optional activity tracks.
