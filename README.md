# nummet

Implementation of various algorithms for Numerical Methods class, in Fortran.

Mainly exposed as a library. Support for CMake and [FPM](https://fpm.fortran-lang.org/) is planned.

Done in a Test Driven Development fashion using [test-drive](https://github.com/fortran-lang/test-drive)

## Contents

1. Interpolation:
    1. [x] Newton's Forward
    2. [x] Newton's Backward
    3. [x] Newton's Divided Difference
    4. [x] Lagrange's
2. DiffyQ:
    1. [x] Range-Kutta 2 (Modified Euler)
    2. [x] Euler's Method
    3. [x] Range-Kutta 4
    4. [ ] Milne's Predictor Corrector 
    7. [ ] Taylor Series Method
3. Algebraic Equations:
    1. [x] Bisection 
    2. [ ] Newton-Rhapson
    3. [x] Regula Falsi
4. System of Linear Equations
    1. [ ] Gauss Elimination
    2. [ ] Gaus Jordan
    3. [ ] Gauss Seidel
    4. [x] LU Factorisation [Not Implementing. Use `stdlib_linalg`->`solve_lu`]
5. Finite Differences BVP
    1. [ ] 2nd Order Solver
6. TBD

> [!NOTE]
> Please check `iter_zero` after using `set_first_interval`. Only continue if `iter_zero` is false. If it is true then that indicates that the first interval search was not successfuly. You can override this behaviour by manually specifying the first interval while defining an instance of `algeb_prob` and setting `iter_zero` to false and setting correct `[a,b]` for it. When you do want to use `set_first_interval`, you can ignore `iter_zero` by not providing it in the constructor and setting any random `a,b`

> [!CAUTION]
> DiffyQ's return value is currently not precision controlled. It is on TODO 
for now. The output is correct.

Note: DiffyQ's support first order derivatives only. The `diffyq_prob` type should be used as a derived type to set a problem using `init_prob`. The `eqn_interface` provides a general interface to represent the diffyq. 

## Running Tests

```bash
cd nummet
fpm build
fpm test
```

## TODO

- [ ] GNUPlot Integration ([ogpf](https://github.com/kookma/ogpf))
- [ ] Switch to using ISO types.
- [ ] Look into real(dp) for diffyQ
- [ ] Use `select type` to extract the `solve` implementations in the diffyQ modules as a single function[^1].

[^1]: I am aware of the code repetition in the diffyq modules (the driver function is very much the same). The only reason they are kept that way is the function call changes in each. A simple `select type` would do the trick is what I am expecting but in an ideal case some metaprogramming way would be nice.


## Notes

- Newton's Backward: p + i - 1
- Newton's Forward:  p - i + 1