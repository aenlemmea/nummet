# nummet

Implementation of various algorithms for Numerical Methods class, in Fortran.

Mainly exposed as a library. Support for CMake and [FPM](https://fpm.fortran-lang.org/) is planned.

Done in a Test Driven Development fashion using [test-drive](https://github.com/fortran-lang/test-drive)

## Contents

1. Interpolation:
    1. [ ] Newton's Forward
    2. [ ] Newton's Backward
    3. [ ] Newton's Divided Difference
    4. [ ] Lagrange's
2. DiffyQ:
    1. [ ] Range-Kutta 2nd Order (Modified Euler)
    2. [ ] Euler's Method
3. Algebraic Equations:
    1. [ ] Bisection 
    2. [ ] Newton-Rhapson
    3. [ ] Regula Falsi
4. TBD

## Running Tests

```bash
cd nummet
fpm build
fpm test
```

## TODO

- [ ] Switch to using ISO types.
