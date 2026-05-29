Cubic Hermite Spline
====================

![CI](https://github.com/Galdeano/CubicHermiteSpline/actions/workflows/ci.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/Galdeano/CubicHermiteSpline/badge.svg?branch=master&service=github)](https://coveralls.io/github/Galdeano/CubicHermiteSpline?branch=master)

A library to compute cubic Hermite spline.


## Compilation

In order to compile:

```sh
$ mkdir build && cd build
$ cmake ..
$ make 
```

Unitary test can be executed with:

```sh
$ ctest
```
## Usage

Then, you can launch a demo:

  * simple demo:

```sh
$ src/CHSplineDemo 
```

  * or with parameters:

```sh
$ src/CHSplineDemo t0 t1 p0 p1 v0 v1 Np 
```
Where:
* **t0** : Initial time
* **t1** : Final time
* **p0** : Initial position
* **p1** : Final position
* **v0** : Initial velocity
* **v1** : Final velocity
* **Np** : Number of point

Like:

```sh
$ src/CHSplineDemo 0.0 5.0 -1.0 5.0 -1.0 0.0 100 
```

## Monotonicity Preservation (Fritsch-Carlson)

By default, Cubic Hermite Splines can overshoot or oscillate when there are sharp steps in the input points. To prevent this, the library includes an implementation of the **Fritsch-Carlson algorithm** for monotonicity-preserving cubic interpolation.

If the input knot points are monotonic, this algorithm adjusts the derivatives so that the interpolated spline is guaranteed to also remain strictly monotonic, completely avoiding overshoot or undershoot.

### Usage

This feature is fully backward-compatible and opt-in. You can activate it by invoking the `initDerivativeMonotonicFritschCarlson()` builder after initializing your spline coordinates:

```cpp
#include "CHSpline/CHSpline.h"

// Define knot points (monotonic)
std::vector<double> ti = {0.0, 1.0, 2.0, 3.0};
std::vector<double> pi = {0.0, 1.0, 1.1, 2.1};
std::vector<double> vi = {0.0, 0.0}; // boundary conditions

CHSpline::Spline spline;
spline.initSpline(ti, pi, vi);

// Opt-in to monotonicity-preserving derivatives
spline.initDerivativeMonotonicFritschCarlson();

// The spline is now guaranteed not to overshoot [1.0, 1.1] on the second segment
double val = spline.evalSpline(1.5); // returns a value strictly in [1.0, 1.1]
```


