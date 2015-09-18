
#define BOOST_TEST_MODULE TestsCHSpline
#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <sstream>
#include <map>
#include <vector> 
#include "CHSpline.h"


BOOST_AUTO_TEST_SUITE(TestsCHSpline)

// ------------- Tests Follow --------------
BOOST_AUTO_TEST_CASE(DefaultConstructors)
{

  Spline Sp;
  Spline Sp2;
  

  BOOST_CHECK(Sp.tAcces() == Sp2.tAcces());
  BOOST_CHECK(Sp.pAcces() == Sp2.pAcces());
  BOOST_CHECK(Sp.vAcces() == Sp2.vAcces());
  
}

BOOST_AUTO_TEST_CASE(DefaultConstructors_VectorSize)
{
  Spline Sp;
  
  BOOST_CHECK_EQUAL(Sp.tAcces().size(), 2);
  BOOST_CHECK_EQUAL(Sp.pAcces().size(), 2);
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), 2);
}

BOOST_AUTO_TEST_CASE(InitialiserDouble_VectorSizeCheck)
{
  Spline Sp;
  
  double t0=0.0, t1=2.0, p0=1.0, p1=5.0, v0=-1.0, v1=4.0; 
  
  Sp.initSpline(t0,t1,p0,p1,v0,v1);
  
  BOOST_CHECK_EQUAL(Sp.tAcces().size(), 2);
  BOOST_CHECK_EQUAL(Sp.pAcces().size(), 2);
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), 2);
}


BOOST_AUTO_TEST_CASE(InitialiserDouble_VectorContentCheck)
{
  Spline Sp;
  
  double t0=0.0, t1=2.0, p0=1.0, p1=5.0, v0=-1.0, v1=4.0; 
  
  Sp.initSpline(t0,t1,p0,p1,v0,v1);
  
  BOOST_CHECK_EQUAL(Sp.tAcces().front(), t0);  
  BOOST_CHECK_EQUAL(Sp.tAcces().back(), t1);
  BOOST_CHECK_EQUAL(Sp.pAcces().front(), p0);  
  BOOST_CHECK_EQUAL(Sp.pAcces().back(), p1);
  BOOST_CHECK_EQUAL(Sp.vAcces().front(), v0);  
  BOOST_CHECK_EQUAL(Sp.vAcces().back(), v1);

}

BOOST_AUTO_TEST_CASE(InitialiserVector_VectorContentCheck)
{
  Spline Sp;
  
  double t0=0.0, t1=2.2, t2=4.0, t3=6.2, p0=1.1, p1=5.0, p2=2.0, p3=5.4, v0=-1.0, v1=4.2, v2=3.0, v3=2.2; 
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();
  ti.push_back(t0);
  ti.push_back(t1);
  ti.push_back(t2);
  ti.push_back(t3);
  
  pi.push_back(p0);
  pi.push_back(p1);
  pi.push_back(p2);
  pi.push_back(p3);
  
  vi.push_back(v0);
  vi.push_back(v1);
  vi.push_back(v2);
  vi.push_back(v3);
  
  Sp.initSpline(ti,pi,vi);
  
  BOOST_CHECK_EQUAL(Sp.tAcces().size(), ti.size());  
  for( int i = 0; i < ti.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.tAcces()[i], ti[i]);
  }
  
  BOOST_CHECK_EQUAL(Sp.pAcces().size(), pi.size());  
  for( int i = 0; i < pi.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.pAcces()[i], pi[i]);
  }
  
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), vi.size());  
  for( int i = 0; i < vi.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.vAcces()[i], vi[i]);
  }
  
}

BOOST_AUTO_TEST_CASE(InitialiserVector_EmptyVector)
{
  Spline Sp;
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();

  BOOST_CHECK_EQUAL(Sp.initSpline(ti,pi,vi), false);
}
  
BOOST_AUTO_TEST_CASE(InitialiserVector_ShortVector)
{
  Spline Sp;
  
  double t0=0.0, t1=2.2, t2=4.0, t3=6.2, p0=1.1, p1=5.0, p2=2.0, p3=5.4, v0=-1.0, v1=4.2, v2=3.0, v3=2.2; 
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();
  
  ti.push_back(t0);
  pi.push_back(p0);
  vi.push_back(v0);
  BOOST_CHECK_EQUAL(Sp.initSpline(ti,pi,vi), false);
}
    
BOOST_AUTO_TEST_CASE(InitialiserVector_DifferentVectorSize)
{
  Spline Sp;
  
  double t0=0.0, t1=2.2, t2=4.0, t3=6.2, p0=1.1, p1=5.0, p2=2.0, p3=5.4, v0=-1.0, v1=4.2, v2=3.0, v3=2.2; 
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();
  ti.push_back(t0);
  ti.push_back(t1);
  ti.push_back(t2);
  
  pi.push_back(p0);
  pi.push_back(p1);
  pi.push_back(p2);
  pi.push_back(p3);
  
  vi.push_back(v0);
  vi.push_back(v1);
  vi.push_back(v2);
  vi.push_back(v3);
  BOOST_CHECK_EQUAL(Sp.initSpline(ti,pi,vi), false);
  
  ti.push_back(t3);
  pi.pop_back();
  BOOST_CHECK_EQUAL(Sp.initSpline(ti,pi,vi), false);
  
  pi.push_back(p3);
  vi.pop_back();
  BOOST_CHECK_EQUAL(Sp.initSpline(ti,pi,vi), false);
}

BOOST_AUTO_TEST_CASE(InitialiserVector_ChronologicalTime)
{
  Spline Sp;
  
  double t0=0.0, t1=2.2, t2=4.0, p0=1.1, p1=5.0, p2=2.0, v0=-1.0, v1=4.2, v2=3.0; 
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();
  ti.push_back(t0);
  ti.push_back(t2);
  ti.push_back(t1);
  
  pi.push_back(p0);
  pi.push_back(p1);
  pi.push_back(p2);
  
  vi.push_back(v0);
  vi.push_back(v1);
  vi.push_back(v2);
  BOOST_CHECK_EQUAL(Sp.initSpline(ti,pi,vi), false);
  
  ti.clear();
  ti.push_back(t1);
  ti.push_back(t0);
  ti.push_back(t2);
  BOOST_CHECK_EQUAL(Sp.initSpline(ti,pi,vi), false);
}

BOOST_AUTO_TEST_CASE(Initialiser_Consistency)
{

  Spline Sp;
  Spline Sp2;
  
  double t0=0.1, t1=2.2, p0=-1.1, p1=-5.4, v0=-1.0, v1=-4.4; 
  
  Sp.initSpline(t0,t1,p0,p1,v0,v1);
  
  
  std::vector<double> ti, pi, vi;
  ti.push_back(t0);
  ti.push_back(t1);
  pi.push_back(p0);
  pi.push_back(p1);
  vi.push_back(v0);
  vi.push_back(v1);
  Sp2.initSpline(ti,pi,vi);
  
  BOOST_CHECK(Sp.tAcces() == Sp2.tAcces());
  BOOST_CHECK(Sp.pAcces() == Sp2.pAcces());
  BOOST_CHECK(Sp.vAcces() == Sp2.vAcces());

}

BOOST_AUTO_TEST_CASE(Initialiser_ResetOldVectors)
{
  Spline Sp;
  
  double t0=1.0, t1=2.0, t2=3.0, t3=4.0, p0=5.0, p1=6.0, p2=7.0, p3=8.0, v0=9.0, v1=10.0, v2=11.0, v3=12.0; 
  
  Sp.initSpline(t0,t1,p0,p1,v0,v1);
  
  Sp.initSpline(t1,t2,p1,p2,v1,v2);

  
  BOOST_CHECK_EQUAL(Sp.tAcces().size(), 2);
  BOOST_CHECK_EQUAL(Sp.pAcces().size(), 2);
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), 2);
  
  BOOST_CHECK_EQUAL(Sp.tAcces().front(), t1);  
  BOOST_CHECK_EQUAL(Sp.tAcces().back(), t2);
  BOOST_CHECK_EQUAL(Sp.pAcces().front(), p1);  
  BOOST_CHECK_EQUAL(Sp.pAcces().back(), p2);
  BOOST_CHECK_EQUAL(Sp.vAcces().front(), v1);  
  BOOST_CHECK_EQUAL(Sp.vAcces().back(), v2);
  
  std::vector<double> ti, pi, vi;
  ti.push_back(t0);
  ti.push_back(t1);
  pi.push_back(p0);
  pi.push_back(p1);
  vi.push_back(v0);
  vi.push_back(v1);
  Sp.initSpline(ti,pi,vi);
  
  BOOST_CHECK_EQUAL(Sp.tAcces().size(), ti.size());
  BOOST_CHECK_EQUAL(Sp.pAcces().size(), pi.size());
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), vi.size());
  
  for( int i = 0; i < ti.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.tAcces()[i], ti[i]);
  }
  for( int i = 0; i < pi.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.pAcces()[i], pi[i]);
  }
  for( int i = 0; i < vi.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.vAcces()[i], vi[i]);
  }

  ti.push_back(t2);
  ti.push_back(t3);
  pi.push_back(p2);
  pi.push_back(p3);
  vi.push_back(v2);
  vi.push_back(v3);
  Sp.initSpline(ti,pi,vi);
  
  BOOST_CHECK_EQUAL(Sp.tAcces().size(), ti.size());
  BOOST_CHECK_EQUAL(Sp.pAcces().size(), pi.size());
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), vi.size());
  
  for( int i = 0; i < ti.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.tAcces()[i], ti[i]);
  }
  for( int i = 0; i < pi.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.pAcces()[i], pi[i]);
  }
  for( int i = 0; i < vi.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.vAcces()[i], vi[i]);
  }

}

BOOST_AUTO_TEST_CASE(Derivative_Manipulation)
{
  Spline Sp;
  
  double t0=0.0, t1=2.2, t2=4.0, t3=6.2, p0=1.1, p1=5.0, p2=2.0, p3=5.4, v0=-1.0, v1=4.2, v2=3.0, v3=2.2; 
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();
  ti.push_back(t0);
  ti.push_back(t1);
  ti.push_back(t2);
  ti.push_back(t3);
  
  pi.push_back(p0);
  pi.push_back(p1);
  pi.push_back(p2);
  pi.push_back(p3);
  
  vi.push_back(v0);
  vi.push_back(v1);
  vi.push_back(v2);
  vi.push_back(v3);
  
  Sp.initSpline(ti,pi,vi);
  
  Sp.initDerivativeCatmullRom();
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), vi.size());  
  BOOST_CHECK_EQUAL(Sp.vAcces().front(), v0);  
  BOOST_CHECK_EQUAL(Sp.vAcces().back(), v3);

  Sp.initDerivativezero();
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), vi.size());  
  BOOST_CHECK_EQUAL(Sp.vAcces().front(), v0);  
  BOOST_CHECK_EQUAL(Sp.vAcces().back(), v3);
  for( int i = 1; i < vi.size()-1; ++i )
  {
    BOOST_CHECK_EQUAL(Sp.vAcces()[i], 0.0);
  }

  Sp.initDerivatives(vi);
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), vi.size());  
  for( int i = 0; i < vi.size(); ++i )
  {
    BOOST_CHECK_EQUAL(Sp.vAcces()[i], vi[i]);
  }
  
}

BOOST_AUTO_TEST_CASE(AddNode_VectorSize)
{
  Spline Sp;
  
  Sp.addNode(5.0,2.0,1.0);
  
  BOOST_CHECK_EQUAL(Sp.tAcces().size(), 3);
  BOOST_CHECK_EQUAL(Sp.pAcces().size(), 3);
  BOOST_CHECK_EQUAL(Sp.vAcces().size(), 3);
}

BOOST_AUTO_TEST_CASE(AddNode_Chronology)
{
  Spline Sp;
  
  double t0=0.0, t1=2.0, p0=1.0, p1=5.0, v0=0.0, v1=1.0; 
  
  Sp.initSpline(t0,t1,p0,p1,v0,v1);
  
  double invalidChronologicTime=1.0, validChronologicTime=5.0;
  
  BOOST_CHECK_EQUAL(Sp.addNode(invalidChronologicTime,0.0,0.0), false);
  BOOST_CHECK_EQUAL(Sp.addNode(validChronologicTime,0.0,0.0), true);
}

BOOST_AUTO_TEST_CASE(EvalSpline_MiddleValue)
{
  Spline Sp;
  
  double t0=0.0, t1=2.0, p0=1.0, p1=5.0, v0=0.0, v1=0.0; 
  
  Sp.initSpline(t0,t1,p0,p1,v0,v1);
  
  BOOST_CHECK_EQUAL(Sp.evalSpline((t0+t1)/2), (p0+p1)/2);
}

BOOST_AUTO_TEST_CASE(EvalSpline_Boundaries)
{
  Spline Sp;
  
  double t0=0.0, t1=2.0, t2=4.0, p0=1.0, p1=5.0, p2=4.0, v0=0.0, v1=-1.0, v2=2.0; 
  
  Sp.initSpline(t0,t1,p0,p1,v0,v1);

  Sp.addNode(t2,p2,v2);
  
  BOOST_CHECK_EQUAL(Sp.evalSpline(t0), p0);
  BOOST_CHECK_EQUAL(Sp.evalSpline(t1), p1);
  BOOST_CHECK_EQUAL(Sp.evalSpline(t2), p2);
}

BOOST_AUTO_TEST_CASE(EvalVectorSplineParam_Boundaries)
{
  Spline Sp;
  
  double t0=0.0, t1=2.2, t2=4.0, t3=6.2, p0=1.1, p1=5.0, p2=2.0, p3=5.4, v0=-1.0, v1=4.2, v2=3.0, v3=2.2; 
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();
  ti.push_back(t0);
  ti.push_back(t1);
  ti.push_back(t2);
  ti.push_back(t3);
  
  pi.push_back(p0);
  pi.push_back(p1);
  pi.push_back(p2);
  pi.push_back(p3);
  
  vi.push_back(v0);
  vi.push_back(v1);
  vi.push_back(v2);
  vi.push_back(v3);
  
  Sp.initSpline(ti,pi,vi);
  
  std::vector<double> eval;
  Sp.evalVectorSpline(ti, eval);
  
  BOOST_CHECK_EQUAL(eval.size(), pi.size()); 
  for( int i = 0; i < ti.size(); ++i )
  {
    BOOST_CHECK_EQUAL(eval[i], pi[i]);
  }
}

BOOST_AUTO_TEST_CASE(EvalVectorSplineReturn_Boundaries2)
{
  Spline Sp;
  
  double t0=0.0, t1=2.2, t2=4.0, t3=6.2, p0=1.1, p1=5.0, p2=2.0, p3=5.4, v0=-1.0, v1=4.2, v2=3.0, v3=2.2; 
  
  std::vector<double> ti, pi, vi;
  ti.clear();
  pi.clear();
  vi.clear();
  ti.push_back(t0);
  ti.push_back(t1);
  ti.push_back(t2);
  ti.push_back(t3);
  
  pi.push_back(p0);
  pi.push_back(p1);
  pi.push_back(p2);
  pi.push_back(p3);
  
  vi.push_back(v0);
  vi.push_back(v1);
  vi.push_back(v2);
  vi.push_back(v3);
  
  Sp.initSpline(ti,pi,vi);
  
  std::vector<double> eval = Sp.evalVectorSpline(ti);
  
  BOOST_CHECK_EQUAL(eval.size(), pi.size()); 
  for( int i = 0; i < ti.size(); ++i )
  {
    BOOST_CHECK_EQUAL(eval[i], pi[i]);
  }
}

BOOST_AUTO_TEST_SUITE_END( )
