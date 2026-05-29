// -*- C++ -*-
/*!
 * @file  Spline.cpp
 * @brief Cubic Hermite splines implementation with first derivative end
 *conditions
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <stdexcept>
#include "CHSpline/CHSpline.h"

namespace CHSpline
{
void Spline::syncKnots()
{
  knots_.clear();
  knots_.reserve(t_.size());
  for (std::size_t i = 0; i < t_.size(); ++i)
  {
    knots_.push_back({t_[i], p_[i], v_[i]});
  }
}

Spline::Spline()
{
  t_.push_back(0.0);
  t_.push_back(1.0);
  p_.push_back(0.0);
  p_.push_back(1.0);
  v_.push_back(0.0);
  v_.push_back(0.0);
  syncKnots();
}

Spline::Spline(double ti0, double ti1, double pi0, double pi1, double vi0,
               double vi1)
{
  t_.push_back(ti0);
  t_.push_back(ti1);
  p_.push_back(pi0);
  p_.push_back(pi1);
  v_.push_back(vi0);
  v_.push_back(vi1);
  syncKnots();
}


bool Spline::initSpline(double ti0, double ti1, double pi0, double pi1,
                        double vi0, double vi1)
{
  t_.clear();
  p_.clear();
  v_.clear();

  t_.push_back(ti0);
  t_.push_back(ti1);
  p_.push_back(pi0);
  p_.push_back(pi1);
  v_.push_back(vi0);
  v_.push_back(vi1);

  syncKnots();
  return true;
}

bool Spline::initSpline(const std::vector<double>& ti, const std::vector<double>& pi,
                        const std::vector<double>& vi)
{
  // check vector sizes
  if ((ti.size() < 2) || (pi.size() < 2) || (vi.size() < 2))
  {
    std::cerr << "Spline initialisation: Error in vector sizes (too small)"
              << std::endl;
    return false;
  }
  if (ti.size() != pi.size())
  {
    std::cerr << "Spline initialisation: Different vector sizes for knot time "
                 "and position" << std::endl;
    return false;
  }
  if ((vi.size() != ti.size()) && (vi.size() != 2))
  {
    std::cerr << "Spline initialisation: Error in velocity vector size"
              << std::endl;
    return false;
  }

  // t is sorted
  for (std::size_t i = 1; i < ti.size(); i++)
  {
    if ((ti[i] - ti[i - 1]) <= 0)
    {
      std::cerr << "Spline initialisation: Error t vector is not sorted"
                << std::endl;
      return false;
    }
  }

  t_.clear();
  p_.clear();
  v_.clear();

  t_ = ti;
  p_ = pi;

  // if more than two knot and only boundaries conditions for velocities
  // use Catmull-Rom Spline construction: V_i=0.5*(p_(i+1)-p_(i+1))
  if ((vi.size() == 2) && (t_.size() != 2))
  {
    v_.push_back(vi[0]);
    for (std::size_t i = 1; i < (t_.size() - 1); i++)
    {
      v_.push_back(0.5 * (((p_[i] - p_[i - 1]) * (t_[i + 1] - t_[i])) /
                              ((t_[i] - t_[i - 1]) * (t_[i + 1] - t_[i - 1])) +
                          ((p_[i + 1] - p_[i]) * (t_[i] - t_[i - 1])) /
                              ((t_[i + 1] - t_[i]) * (t_[i + 1] - t_[i - 1]))));
    }
    v_.push_back(vi[1]);
  }
  else
  {
    v_ = vi;
  }

  syncKnots();
  return true;
}

bool Spline::initDerivativeCatmullRom()
{
  if (t_.size() > 2)
  {
    double vFront = v_.front();
    double vBack = v_.back();
    v_.clear();
    v_.push_back(vFront);
    for (std::size_t i = 1; i < (t_.size() - 1); i++)
    {
      v_.push_back(0.5 * (((p_[i] - p_[i - 1]) * (t_[i + 1] - t_[i])) /
                              ((t_[i] - t_[i - 1]) * (t_[i + 1] - t_[i - 1])) +
                          ((p_[i + 1] - p_[i]) * (t_[i] - t_[i - 1])) /
                              ((t_[i + 1] - t_[i]) * (t_[i + 1] - t_[i - 1]))));
    }
    v_.push_back(vBack);
    syncKnots();
    return true;
  }
  return false;
}

bool Spline::initDerivativeZero()
{
  if (t_.size() > 2)
  {
    double vFront = v_.front();
    double vBack = v_.back();
    v_.clear();
    v_.push_back(vFront);
    for (std::size_t i = 1; i < (t_.size() - 1); i++)
    {
      v_.push_back(0.0);
    }
    v_.push_back(vBack);
    syncKnots();
    return true;
  }
  return false;
}

bool Spline::initDerivativeMonotonicFritschCarlson()
{
  if (t_.size() <= 2)
  {
    return false;
  }

  // 1. Initialize velocities using Catmull-Rom as a starting point
  initDerivativeCatmullRom();

  std::size_t n = t_.size();
  std::vector<double> d(n - 1); // secant slopes

  // 2. Compute secant slopes
  for (std::size_t i = 0; i < n - 1; ++i)
  {
    d[i] = (p_[i + 1] - p_[i]) / (t_[i + 1] - t_[i]);
  }

  // 3. Apply Fritsch-Carlson monotonicity filter
  for (std::size_t i = 0; i < n - 1; ++i)
  {
    double secant = d[i];
    if (std::abs(secant) < 1e-9) // flat segment
    {
      v_[i] = 0.0;
      v_[i + 1] = 0.0;
    }
    else
    {
      // Enforce sign agreement
      if (v_[i] * secant < 0) v_[i] = 0.0;
      if (v_[i + 1] * secant < 0) v_[i + 1] = 0.0;

      double alpha = v_[i] / secant;
      double beta = v_[i + 1] / secant;
      double sumSquare = alpha * alpha + beta * beta;

      if (sumSquare > 9.0)
      {
        double tau = 3.0 / std::sqrt(sumSquare);
        v_[i] = tau * alpha * secant;
        v_[i + 1] = tau * beta * secant;
      }
    }
  }

  syncKnots();
  return true;
}

bool Spline::initDerivatives(const std::vector<double>& vi)
{
  if (vi.size() == t_.size())
  {
    v_ = vi;
    syncKnots();
    return true;
  }
  return false;
}

bool Spline::addNode(double ti, double pi, double vi)
{

  if ((t_.back()) < ti)
  {
    t_.push_back(ti);
    p_.push_back(pi);
    v_.push_back(vi);
    syncKnots();
    return true;
  }
  else
  {
    return false;
  }
}

double Spline::evalSpline(double te) const
{
  if (knots_.empty())
  {
    throw std::runtime_error("Spline is not initialized");
  }
  if (te <= knots_.front().time)
  {
    return knots_.front().position;
  }
  if (te >= knots_.back().time)
  {
    return knots_.back().position;
  }

  // Find the right knot using binary search (O(log n))
  auto it = std::upper_bound(knots_.begin(), knots_.end(), te,
                             [](double val, const Knot& knot) {
                               return val < knot.time;
                             });
  std::size_t noSpline = static_cast<std::size_t>(std::distance(knots_.begin(), it) - 1);

  double tn = (te - knots_[noSpline].time) /
              (knots_[noSpline + 1].time - knots_[noSpline].time);  // normalized time

  double h1 = 2 * tn * tn * tn - 3 * tn * tn + 1;  // calculate basis function 1
  double h2 = -2 * tn * tn * tn + 3 * tn * tn;     // calculate basis function 2
  double h3 = tn * tn * tn - 2 * tn * tn + tn;     // calculate basis function 3
  double h4 = tn * tn * tn - tn * tn;              // calculate basis function 4
  return h1 * knots_[noSpline].position +      // multiply and sum all functions
         h2 * knots_[noSpline + 1].position +  // together to build the interpolated
         h3 * knots_[noSpline].velocity *
             (knots_[noSpline + 1].time - knots_[noSpline].time) +  // point along the curve.
         h4 * knots_[noSpline + 1].velocity * (knots_[noSpline + 1].time - knots_[noSpline].time);
}

bool Spline::evalVectorSpline(const std::vector<double>& t,
                              std::vector<double>& output) const
{
  output.clear();
  for (std::size_t i = 0; i < t.size(); i++)
  {
    output.push_back(evalSpline(t[i]));
  }
  return true;
}

std::vector<double> Spline::evalVectorSpline(const std::vector<double>& t) const
{
  std::vector<double> output;
  for (std::size_t i = 0; i < t.size(); i++)
  {
    output.push_back(evalSpline(t[i]));
  }
  return output;
}
} // namespace CHSpline
