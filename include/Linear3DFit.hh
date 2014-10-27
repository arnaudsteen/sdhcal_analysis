#ifndef LINEAR3DFIT_HH
#define LINEAR3DFIT_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "ThreeVector.hh"
/*!
 * Class Linear3DFit.
 * Fit a 3D line in 3D space. The line equation is the following one :
 *       x = a*z + b
 *       y = c*z + d
 *       z = t
 */

class Linear3DFit {

protected :
  std::vector<ThreeVector> positions;
  std::vector<int> clSize;
  float chi2;
  float params[4];
  float paramsError[4];
  void ComputeChi2();

public :
  Linear3DFit( std::vector<ThreeVector> pos , std::vector<int> ClustersSize);
  virtual ~Linear3DFit();
  void Fit();
  inline float GetChi2(){ ComputeChi2(); return chi2; }
  ThreeVector VectorFromRealLine( const ThreeVector& );
  ThreeVector GetNormaleProjection( const ThreeVector& );
  inline float* GetFitParameters(){ return params; }
  inline float* GetFitParError(){ return paramsError; }
  inline std::vector<int> GetClusterSize(){return clSize;}
};
#endif
