#ifndef SPLINE_H
#define SPLINE_H

#include <vector>
#include <deque>
#include <cmath>

class Spline
{
 public:
  /**
     \brief Constructor for Cubic Hermite splines implementation with first
     derivative end conditions
  */
  Spline();

  /**
     \brief Destructor
  */
  ~Spline();

  /**
     \brief Initalise spline coefficients for two knot
  */
  bool initSpline(double t0, double t1, double p0, double p1, double v0,
                  double v1);

  /**
     \brief Initalise spline coefficients
  */
  bool initSpline(std::vector<double> ti, std::vector<double> pi,
                  std::vector<double> vi);

  /**
     \brief Derivatives calculation using Catmull-Rom Spline construction for
     minimal accelerations
  */
  bool initDerivativeCatmullRom();

  /**
     \brief Derivatives set to zero for avoiding overshot
  */
  bool initDerivativezero();

  /**
     \brief Derivatives set to vector values
  */
  bool initDerivatives(std::vector<double> vi);

  /**
     \brief Add a node
  */
  bool addNode(double ti, double pi, double vi);

  /**
     \brief Evaluate spline for one value
  */
  double evalSpline(double t) const;

  /**
     \brief Evaluate spline for a vector of values
  */
  bool evalVectorSpline(std::vector<double> t,
                        std::vector<double>& output) const;
  std::vector<double> evalVectorSpline(std::vector<double> t) const;

  /**
     \brief Accessors
  */
  const std::vector<double>& getTime() const
  {
    return t_;
  }
  const std::vector<double>& getPosition() const
  {
    return p_;
  }
  const std::vector<double>& getVelocity() const
  {
    return v_;
  }

 private:
  std::vector<double> t_, p_, v_;
};

#endif
