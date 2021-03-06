#include "Distance.h"
#include <cmath>
#include "marlin/VerbosityLevels.h"

ThreeVector Distance::VectorProduct(ThreeVector vec1, ThreeVector vec2)
{
  ThreeVector vec( vec1.y()*vec2.z()-vec1.z()*vec2.y(),
		   vec1.z()*vec2.x()-vec1.x()*vec2.z(),
		   vec1.x()*vec2.y()-vec1.y()*vec2.x() );
  return vec;
}

float Distance::VectorNorm(ThreeVector vec)
{
  return sqrt(vec.x()*vec.x() + vec.y()*vec.y() + vec.z()*vec.z());
}

//----------------------------------------------------------------------------------

DistanceBetweenTwoHits::DistanceBetweenTwoHits() : Distance() , hit1(0), hit2(0)
{}

void DistanceBetweenTwoHits::Init(EVENT::CalorimeterHit* it1,EVENT::CalorimeterHit* it2)
{
  hit1=it1;
  hit2=it2;
}

float DistanceBetweenTwoHits::CalculateDistance()
{
  ThreeVector AB(hit1->getPosition()[0]-hit2->getPosition()[0],
		 hit1->getPosition()[1]-hit2->getPosition()[1],
		 hit1->getPosition()[2]-hit2->getPosition()[2]);
  return AB.mag();
}

//----------------------------------------------------------------------------------

DistanceBetweenTwoClusters::DistanceBetweenTwoClusters() : Distance() , cluster1(0), cluster2(0)
{}

void DistanceBetweenTwoClusters::Init(Cluster* c1,Cluster* c2)
{
  cluster1=c1;
  cluster2=c2;
}

float DistanceBetweenTwoClusters::CalculateDistance()
{
  ThreeVector AB(cluster1->getClusterPosition().x()-cluster2->getClusterPosition().x(),
		 cluster1->getClusterPosition().y()-cluster2->getClusterPosition().y(),
		 cluster1->getClusterPosition().z()-cluster2->getClusterPosition().z());
  return AB.mag();
}

//----------------------------------------------------------------------------------

DistanceBetweenOneHitAndOneCluster::DistanceBetweenOneHitAndOneCluster() : Distance() , hit(0), cluster(0)
{}

void DistanceBetweenOneHitAndOneCluster::Init(EVENT::CalorimeterHit* it,Cluster* cl)
{
  hit=it;
  cluster=cl;
}

float DistanceBetweenOneHitAndOneCluster::CalculateDistance()
{
  ThreeVector AB(hit->getPosition()[0]-cluster->getClusterPosition().x(),
		 hit->getPosition()[1]-cluster->getClusterPosition().y(),
		 hit->getPosition()[2]-cluster->getClusterPosition().z());
  return AB.mag();
}

//----------------------------------------------------------------------------------

DistanceBetweenOneClusterAndOneTrack::DistanceBetweenOneClusterAndOneTrack() : Distance(), trackParams(0), Nx(0), Ny(0), u(0), B(0), normU(0)
{
  /*
    cluster C(x_c,y_c,z_c)
    track T : x_t = param[0] + param[1]*z_t => plan equation; normal vector Nx(-1,0,param[1])
    track T : y_t = param[2] + param[3]*z_t => plan equation; normal vector Ny(0,-1,param[3])
    track T orientation vector u  : u = Nx * Ny
    d(C,T) = || vec(BC) * u || / || u || where B is a point from the track
  */
}

void DistanceBetweenOneClusterAndOneTrack::Init(std::vector<float> &params)
{
  trackParams=params;
  //Nx,Ny : plans containing the track
  Nx = ThreeVector(-1., 0., trackParams[1]);
  Ny = ThreeVector(0., -1., trackParams[3]);
  //u : track orientation vector
  u=Nx.cross(Ny);
  normU=u.mag();
  //B : a track point 
  B = ThreeVector(trackParams[0],trackParams[2],0.);
}

float DistanceBetweenOneClusterAndOneTrack::CalculateDistance(Cluster* cluster)
{
  ThreeVector C(cluster->getClusterPosition().x(),cluster->getClusterPosition().y(),cluster->getClusterPosition().z());
  ThreeVector BC(B.x()-C.x(),B.y()-C.y(),B.z()-C.z());
  if(normU>0)
    return (BC.cross(u)).mag()/normU;
  else{
    std::cout << "ORIENTATION VECTOR u IS NULL => should return exception (cluster and track)" << std::endl;
    std::cout << "trackParams[1] = " << trackParams[1] 
    	      << "\t trackParams[3] = " << trackParams[3] << std::endl;
    return 0;
  }
}

//----------------------------------------------------------------------------------

DistanceBetweenOneHitAndOneTrack::DistanceBetweenOneHitAndOneTrack() : Distance(), trackParams(0), Nx(0), Ny(0), u(0), B(0), normU(0)
{  
  /*
    hit H(x_c,y_c,z_c)
    track T : x_t = param[0] + param[1]*z_t => plan equation; normal vector Nx(-1,0,param[1])
    track T : y_t = param[2] + param[3]*z_t => plan equation; normal vector Ny(0,-1,param[3])
    track T orientation vector u  : u = Nx * Ny
    d(H,T) = || vec(BH) * u || / || u || where B is a point from the track
  */
}

void DistanceBetweenOneHitAndOneTrack::Init(std::vector<float> &params)
{
  trackParams=params;
  //Nx,Ny : plans containing the track
  Nx = ThreeVector(-1., 0., trackParams[1]);
  Ny = ThreeVector(0., -1., trackParams[3]);
  //u : track orientation vector
  u=Nx.cross(Ny);
  normU=u.mag();
  //B : a track point 
  B = ThreeVector(trackParams[0],trackParams[2],0.);
}

float DistanceBetweenOneHitAndOneTrack::CalculateDistance(EVENT::CalorimeterHit* hit)
{
  ThreeVector H(hit->getPosition()[0],hit->getPosition()[1],hit->getPosition()[2]);
  ThreeVector BH(B.x()-H.x(),B.y()-H.y(),B.z()-H.z());
  if(normU>0)
    return (BH.cross(u)).mag()/normU;
  else{
    streamlog_out( MESSAGE ) << "ORIENTATION VECTOR u IS NULL => should return exception (hit and track)" << std::endl;
    streamlog_out( MESSAGE ) << "DistanceBetweenOneHitAndOneTrack::CalculateDistance==>\t"
			     << trackParams[1] << "\t"
			     << trackParams[3] << std::endl;
    return 0;
  }
}

//----------------------------------------------------------------------------------

