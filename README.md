# GeneralizedSchurAlgorithm

This package aims at providing backward stable version [1] of generalized Schur
algorithm for indefinite matrices. This package is under development.

### Build Status and Code Coverage

-  Build status: [![Build Status][build-img]][build-link]
-  Code coverage: [![Coveralls][ca-img]][ca-link] [![Codecov][cc-img]][cc-link]

[build-img]:  https://travis-ci.org/KTH-AC/GeneralizedSchurAlgorithm.jl.svg?branch=master
[build-link]: https://travis-ci.org/KTH-AC/GeneralizedSchurAlgorithm.jl
[ca-img]: https://coveralls.io/repos/github/KTH-AC/GeneralizedSchurAlgorithm.jl/badge.svg?branch=master
[ca-link]: https://coveralls.io/github/KTH-AC/GeneralizedSchurAlgorithm.jl?branch=master
[cc-img]: https://codecov.io/gh/KTH-AC/GeneralizedSchurAlgorithm.jl/branch/master/graph/badge.svg
[cc-link]: https://codecov.io/gh/KTH-AC/GeneralizedSchurAlgorithm.jl

### Description

This package aims at providing backward stable version[1] of generalized Schur
algorithm for indefinite matrices. In particular the goal is to provide some of
the functionality regarding Toeplitz matrices from the
[SLICOT Library](https://github.com/KTH-AC/slicot) [2].

### References
-  [1] [Kailath, T. and Sayed, A. Fast Reliable Algorithms for Matrices with
        Structure. SIAM Publications, Philadelphia,
        1999](http://epubs.siam.org/doi/book/10.1137/1.9781611971354)
-  [2] [Kressner, D. and Van Dooren, P. Factorizations and linear system solvers
        for matrices with Toeplitz structure. SLICOT Working Note 2000-2,
        2000](www.icm.tu-bs.de/NICONET/REPORTS/SLWN2000-2.ps.gz)
