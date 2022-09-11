# What is this?

C code implementation of splitting or merging algorithms for determining the number of clusters.

It is intended to be called in the R language.

# Requirement

- `mclust` package (for em-based method)
- `kohonen` package (for smlsom)
- `gtools` package (for smlsom)

# Direcotries
- `src` directory contains three directories:
  - `smlsom`
  - `two_means`
  - `em_based`

- `two_means` and `em_based` contain the following algorithms
  - two_means:
    - x-means: Pelleg, Dan, and Andrew W. Moore. "X-means: Extending k-means with efficient estimation of the number of clusters." The 17th International Conference on Machine Learning (2000): 727-734.
    - g-means: Hamerly, Greg, and Charles Elkan. "Learning the k in k-means." Advances in neural information processing systems 16 (2003).
    - dip-means: Kalogeratos, Argyris, and Aristidis Likas. "Dip-means: an incremental clustering method for estimating the number of clusters." Advances in neural information processing systems 25 (2012).
  - em_based:
    - pg-means: Feng, Yu, and Greg Hamerly. "PG-means: learning the number of clusters in data." Advances in neural information processing systems 19 (2006).
    - mml-em: Figueiredo, Mario A. T., and Anil K. Jain. "Unsupervised learning of finite mixture models." IEEE Transactions on pattern analysis and machine intelligence 24.3 (2002): 381-396.

# Compile for using in R
- Compiling as follows and move all *.so files to `lib` directory
- In your src/two_means
  - x-means:   `R CMD SHLIB -c xmeans.c twomeans.c list.c utils.c`
  - g-means:   `R CMD SHLIB -c gmeans.c twomeans.c list.c utils.c`
  - dip-means: `R CMD SHLIB -c dipmeans.c twomeans.c diptest.c dip.c list.c utils.c`
- In your src/em_based
  - pg-means:  `R CMD SHLIB -c pgmeans.c em_gmm.c kstest.c utils.c`
  - mml-em:    `R CMD SHLIB -c mixt4.c em_gmm.c utils.c`
- In your src/smlsom `R CMD SHLIB -c *.c -o smlsom.so`

- Use in your R as: `source("[algorithm name].R")`

# Note:
 The C implementation of smlsom is experimental code and very hard to read.
 "https://github.com/lipryou/smlsom" is a reorganized version in C++ and R for easier maintenance and readability.
 Note that while it is much easier to maintain, the execution speed is slower than the C code.
