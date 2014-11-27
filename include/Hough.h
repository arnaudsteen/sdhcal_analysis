#ifndef HOUGH_HH
#define HOUGH_HH

#include <iostream>
#include <vector>
#include <math.h>
#include <algorithm>
#include "Cluster.h"
#include "Track.h"
#include "HoughPoint.h"
#include "Distance.h"

static const unsigned int cutValue=6;

typedef struct{
  std::vector<Cluster*> clusters;
  int theta;
  int rho;
}HoughBin;

class Hough
{
 public:
  Hough();
  ~Hough();
  void Init(std::vector<Cluster*>& clVec);
  void ComputeHoughTransform();
  std::vector<HoughBin> getHoughSpace(std::vector<Cluster*> &clVec, bool ZX);
  void RemoveIsolatedClusters(std::vector<Cluster*>& clVec);
  void RemoveTrackedClusters(std::vector<HoughBin>& hBinVec);
  inline void setClusters(std::vector<Cluster*>& clVec){clusters=clVec;}
  inline void addTrack(Track* patatrack){tracks.push_back(patatrack);}
  inline std::vector<Cluster*>& getClusters(){return clusters;}
  inline std::vector<Track*>& ReturnTracks(){return tracks;}
 private:
  std::vector<Cluster*> clusters;
  std::vector<Track*> tracks;
};

class HoughBinFunction
{
 public:
  HoughBinFunction(){;}
  ~HoughBinFunction(){;}
  static bool SortBySize(HoughBin hBin1, HoughBin hBin2){return hBin1.clusters.size()>hBin2.clusters.size();}
  static bool SmallBin(HoughBin hBin){return hBin.clusters.size()<cutValue;}
  static bool TrackInLayer(HoughBin hBin){return hBin.theta>=48&&hBin.theta<=52;}
};

#endif
