// -*- C++ -*-
/*!
 * @file  Spline.cpp
 * @brief Cubic Hermite splines implementation with first derivative end
 *conditions
 * $Date$
 *
 * $Id$
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include "CHSpline/CHSpline.h"

Spline::Spline()
{
  t_.push_back(0.0);
  t_.push_back(1.0);
  p_.push_back(0.0);
  p_.push_back(1.0);
  v_.push_back(0.0);
  v_.push_back(0.0);
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
}

Spline::~Spline()
{
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

  return true;
}

bool Spline::initSpline(std::vector<double> ti, std::vector<double> pi,
                        std::vector<double> vi)
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
  for (std::vector<double>::size_type i = 1; i < ti.size(); i++)
  {
    if ((ti[i] - ti[i - 1]) < 0)
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
    for (std::vector<double>::size_type i = 1; i < (t_.size() - 1); i++)
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
    for (std::vector<double>::size_type i = 1; i < (t_.size() - 1); i++)
    {
      v_.push_back(0.5 * (((p_[i] - p_[i - 1]) * (t_[i + 1] - t_[i])) /
                              ((t_[i] - t_[i - 1]) * (t_[i + 1] - t_[i - 1])) +
                          ((p_[i + 1] - p_[i]) * (t_[i] - t_[i - 1])) /
                              ((t_[i + 1] - t_[i]) * (t_[i + 1] - t_[i - 1]))));
    }
    v_.push_back(vBack);
    return true;
  }
  return false;
}

bool Spline::initDerivativezero()
{
  if (t_.size() > 2)
  {
    double vFront = v_.front();
    double vBack = v_.back();
    v_.clear();
    v_.push_back(vFront);
    for (std::vector<double>::size_type i = 1; i < (t_.size() - 1); i++)
    {
      v_.push_back(0.0);
    }
    v_.push_back(vBack);
    return true;
  }
  return false;
}

bool Spline::initDerivatives(std::vector<double> vi)
{
  if (vi.size() == t_.size())
  {
    v_ = vi;
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
    return true;
  }
  else
  {
    return false;
  }
}

double Spline::evalSpline(double te) const
{
  if (te <= t_.front())
  {
    return p_.front();
  }
  if (te >= t_.back())
  {
    return p_.back();
  }

  // Find the right knot
  std::vector<double>::size_type noSpline = 0;
  for (std::vector<double>::size_type i = 0; i < t_.size(); ++i)
  {
    if (!(te < t_[i]))
    {
      noSpline = i;
    }
  }

  double tn = (te - t_[noSpline]) /
              (t_[noSpline + 1] - t_[noSpline]);  // normalized time

  double h1 = 2 * tn * tn * tn - 3 * tn * tn + 1;  // calculate basis function 1
  double h2 = -2 * tn * tn * tn + 3 * tn * tn;     // calculate basis function 2
  double h3 = tn * tn * tn - 2 * tn * tn + tn;     // calculate basis function 3
  double h4 = tn * tn * tn - tn * tn;              // calculate basis function 4
  return h1 * p_[noSpline] +      // multiply and sum all functions
         h2 * p_[noSpline + 1] +  // together to build the interpolated
         h3 * v_[noSpline] *
             (t_[noSpline + 1] - t_[noSpline]) +  // point along the curve.
         h4 * v_[noSpline + 1] * (t_[noSpline + 1] - t_[noSpline]);
}

bool Spline::evalVectorSpline(std::vector<double> t,
                              std::vector<double>& output) const
{
  output.clear();
  for (std::vector<double>::size_type i = 0; i < t.size(); i++)
  {
    output.push_back(evalSpline(t[i]));
  }
  return true;
}

std::vector<double> Spline::evalVectorSpline(std::vector<double> t) const
{
  std::vector<double> output;
  for (std::vector<double>::size_type i = 0; i < t.size(); i++)
  {
    output.push_back(evalSpline(t[i]));
  }
  return output;
}
