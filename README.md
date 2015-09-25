Cubic Hermite Spline
====================

A library to compute cubic Hermite spline.

####Master branch build status: 
![](https://travis-ci.org/Galdeano/CubicHermiteSpline.svg?branch=master)

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

