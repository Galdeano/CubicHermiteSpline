# Cubic Hermite Spline (CHSpline)

![CI](https://github.com/Galdeano/CubicHermiteSpline/actions/workflows/ci.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/Galdeano/CubicHermiteSpline/badge.svg?branch=master&service=github)](https://coveralls.io/github/Galdeano/CubicHermiteSpline?branch=master)

A lightweight, high-performance, and mathematically robust C++ library to compute **Cubic Hermite Spline** interpolations.

This library is designed for robotics, computer graphics, path planning, and trajectory generation where smooth, first-derivative-continuous curves are required.

---

## Key Features

* **High Performance**:
  * **Cache Locality (AoS Layout)**: Spline knots are stored in a contiguous Array-of-Structures (`std::vector<Knot>`) layout, optimizing CPU L1/L2 data cache hit ratios during evaluation loops.
  * **Binary Search Evaluation**: Segment lookup is optimized to $O(\log n)$ using `std::upper_bound`, making it extremely fast even for splines with thousands of knots.
* **API Safety & Cleanliness**:
  * **Namespace Protected**: Fully encapsulated inside the `CHSpline` namespace to avoid global naming collisions.
  * **Pass-by-Reference**: Minimal memory copies achieved by passing vectors by `const std::vector<double>&`.
  * **Error Boundary Safety**: Strict runtime exception safety checks prevent evaluations on uninitialized states.
* **Monotonicity Preservation**:
  * Out-of-the-box opt-in implementation of the **Fritsch-Carlson algorithm** to avoid curve overshoots/undershoots and eliminate spurious oscillations.
* **Modern Tooling**:
  * Compatible with modern compiler toolchains (requires CMake 3.5+ and standard C++11 or higher).

---

## Quick Start Integration

You can easily integrate `CHSpline` into your own C++ code. Below is a complete quick-start example:

```cpp
#include <iostream>
#include <vector>
#include "CHSpline/CHSpline.h"

int main()
{
  // Define time, position, and velocity constraints
  std::vector<double> time = {0.0, 1.5, 3.0};
  std::vector<double> pos  = {0.0, 4.5, 2.0};
  std::vector<double> vel  = {0.0, 0.0}; // Boundary velocities

  // Initialize the spline
  CHSpline::Spline spline;
  if (spline.initSpline(time, pos, vel))
  {
    // Auto-calculate intermediary derivatives using Catmull-Rom construction
    spline.initDerivativeCatmullRom();

    // Evaluate position at time t = 1.0
    double p = spline.evalSpline(1.0);
    std::cout << "Position at t = 1.0: " << p << std::endl;
  }
  return 0;
}
```

---

## Compilation & Testing

To compile the library and the demo executable locally, run:

```sh
$ mkdir build && cd build
$ cmake ..
$ cmake --build .
```

To run the unit test suite (requires Boost Test):

```sh
$ ctest --output-on-failure
```

---

## Running the Demo Application

The compilation builds a demo executable `CHSplineDemo` under `src/`. You can launch it to output spline evaluation points:

* **Simple Demo Mode**:
  ```sh
  $ src/CHSplineDemo
  ```

* **Custom Trajectory Mode**:
  ```sh
  $ src/CHSplineDemo t0 t1 p0 p1 v0 v1 Np
  ```
  Where:
  * `t0` / `t1` : Initial / Final time boundaries
  * `p0` / `p1` : Initial / Final position boundaries
  * `v0` / `v1` : Initial / Final velocity boundaries
  * `Np`        : Number of evaluation steps to print

  *Example:*
  ```sh
  $ src/CHSplineDemo 0.0 5.0 -1.0 5.0 -1.0 0.0 100
  ```

---

## Monotonicity Preservation (Fritsch-Carlson)

By default, Cubic Hermite Splines can overshoot or oscillate when there are sharp steps in the input points. To prevent this, the library includes a state-of-the-art implementation of the **Fritsch-Carlson algorithm** for monotonicity-preserving cubic interpolation.

If the input knot points are strictly monotonic, this algorithm adjusts the derivatives so that the interpolated spline is guaranteed to also remain strictly monotonic, completely avoiding overshoot or undershoot.

### Usage

This feature is fully backward-compatible and opt-in. You can activate it by invoking the `initDerivativeMonotonicFritschCarlson()` builder after initializing your spline coordinates:

```cpp
#include "CHSpline/CHSpline.h"

// Define knot points containing a sharp step
std::vector<double> ti = {0.0, 1.0, 2.0, 3.0};
std::vector<double> pi = {0.0, 1.0, 1.1, 2.1}; // Sharp step between t=1 and t=2
std::vector<double> vi = {0.0, 0.0}; 

CHSpline::Spline spline;
spline.initSpline(ti, pi, vi);

// Opt-in to monotonicity-preserving derivatives
spline.initDerivativeMonotonicFritschCarlson();

// The spline is now mathematically guaranteed to remain within [1.0, 1.1] 
// on the second segment, eliminating overshoot!
double val = spline.evalSpline(1.5); 
```

---

## License

This library is distributed under the BSD 3-Clause License. See [COPYING](COPYING) for details.
