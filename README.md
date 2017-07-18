# GeneralizedSchurAlgorithm

[![Unix][unix-img]][unix-link]
[![Coveralls][ca-img]][ca-link]
[![Codecov][cc-img]][cc-link]

This package aims at providing backward stable version [1] of generalized Schur
algorithm for indefinite matrices. This package is under development.

[unix-img]: https://img.shields.io/travis/neveritt/GeneralizedSchurAlgorithm.jl/master.svg?label=unix
[unix-link]: https://travis-ci.org/neveritt/GeneralizedSchurAlgorithm.jl
[ca-img]: https://img.shields.io/coveralls/neveritt/GeneralizedSchurAlgorithm.jl/master.svg?label=coveralls
[ca-link]: https://coveralls.io/github/neveritt/GeneralizedSchurAlgorithm.jl?branch=master
[cc-img]: https://img.shields.io/codecov/c/github/neveritt/GeneralizedSchurAlgorithm.jl/master.svg?label=codecov
[cc-link]: https://codecov.io/gh/neveritt/GeneralizedSchurAlgorithm.jl?branch=master

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
