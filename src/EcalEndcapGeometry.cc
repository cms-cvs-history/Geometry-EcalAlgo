#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <CLHEP/Geometry/Point3D.h>
#include <CLHEP/Geometry/Plane3D.h>

EcalEndcapGeometry::EcalEndcapGeometry() :
  _nnmods(0 ),_nncrys(0)
{}

EcalEndcapGeometry::~EcalEndcapGeometry() {}

void 
EcalEndcapGeometry::initialize()
{
  zeP=0.;
  zeN=0.;
  unsigned nP=0;
  unsigned nN=0;
  m_nref = 0 ;
  for( CaloSubdetectorGeometry::CellCont::const_iterator i ( cellGeometries().begin() ) ; 
       i != cellGeometries().end(); ++i )
  {
//      addCrystalToZGridmap(i->first,dynamic_cast<const TruncatedPyramid*>(i->second));
     float z=dynamic_cast<const TruncatedPyramid*>(i->second)->getPosition(0.).z();
     if(z>0.)
     {
	zeP+=z;
	++nP;
     }
     else
     {
	zeN+=z;
	++nN;
     }
     const EEDetId myId ( i->first ) ;
     const unsigned int ix ( myId.ix() ) ;
     const unsigned int iy ( myId.iy() ) ;
     if( abs( ix ) > m_nref ) m_nref = abs( ix ) ;
     if( abs( iy ) > m_nref ) m_nref = abs( iy ) ;
  }
  zeP/=(float)nP;
  zeN/=(float)nN;

  m_href = 0 ;
  for( CaloSubdetectorGeometry::CellCont::const_iterator i ( cellGeometries().begin() ) ; 
       i != cellGeometries().end(); ++i )
  {
     const EEDetId myId ( i->first ) ;
     const unsigned int ix ( myId.ix() ) ;
     const unsigned int iy ( myId.iy() ) ;
     if( ix == m_nref ) 
     {
	const GlobalPoint p ( dynamic_cast<const TruncatedPyramid*>
			      (i->second)->getPosition(0.)  ) ;
	const float x ( p.x()*fabs(zeP/p.z()) ) ;
	if( m_href < x ) m_href = x ;
     }
     if( iy == m_nref ) 
     {
	const GlobalPoint p ( dynamic_cast<const TruncatedPyramid*>
			      (i->second)->getPosition(0.)  ) ;
	const float x ( p.y()*fabs(zeP/p.z()) ) ;
	if( m_href < x ) m_href = x ;
     }
  }
  m_href = m_href*(1.*m_nref)/( 1.*m_nref - 1. ) ;
  m_wref = m_href/(0.5*m_nref) ;
}


unsigned int 
EcalEndcapGeometry::index( float x ) const
{
   int i ( 1 + (int)floor( ( x + m_href )/m_wref ) ) ;
   if( 1      > i )
   {
      i = 1 ;
   }
   else
   {
      if( m_nref < (unsigned int)i ) i = m_nref ;
   }
   return i ;
}

EEDetId 
EcalEndcapGeometry::gId( float x, 
			 float y, 
			 float z ) const
{
   const double       fac ( fabs(zeP/z) ) ;
   const unsigned int ix  ( index( x*fac ) ) ; 
   const unsigned int iy  ( index( y*fac ) ) ; 
   const unsigned int iz  ( z>0 ? 1 : -1 ) ;

   if( EEDetId::validDetId( ix, iy, iz ) ) 
   {
      return EEDetId( ix, iy, iz ) ; // first try is on target
   }
   else // try nearby coordinates, spiraling out from center
   {
      for( unsigned int i ( 1 ) ; i != 6 ; ++i )
      {
	 for( unsigned int k ( 0 ) ; k != 8 ; ++k )
	 {
	    const int jx ( 0 == k || 4 == k || 5 == k ? +i :
			   ( 1 == k || 5 < k ? -i : 0 ) ) ;
	    const int jy ( 2 == k || 4 == k || 6 == k ? +i :
			   ( 3 == k || 5 == k || 7 == k ? -i : 0 ) ) ;
	    if( EEDetId::validDetId( ix + jx, iy + jy, iz ) ) 
	    {
	       return EEDetId( ix + jx, iy + jy, iz ) ;
	    }
	 }
      }
   }
   return EEDetId() ; // nowhere near any crystal
}


// Get closest cell, etc...
DetId 
EcalEndcapGeometry::getClosestCell( const GlobalPoint& r ) const 
{
   try
   {
      EEDetId mycellID ( gId( r.x(), r.y(), r.z() ) ) ; // educated guess

      if( EEDetId::validDetId( mycellID.ix(), 
			       mycellID.iy(),
			       mycellID.zside() ) )
      {
	 // now get points in convenient ordering

	 HepPoint3D  A;
	 HepPoint3D  B;
	 HepPoint3D  C;
	 HepPoint3D  point(r.x(),r.y(),r.z());
	 // D.K. : equation of plane : AA*x+BB*y+CC*z+DD=0;
	 // finding equation for each edge
	 
	 // ================================================================
	 double x,y,z;
	 unsigned offset=0;
	 int zsign=1;
	 //================================================================
	 std::vector<double> SS;
      
	 // compute the distance of the point with respect of the 4 crystal lateral planes
	 const GlobalPoint& myPosition=getGeometry(mycellID)->getPosition();
	 
	 x=myPosition.x();
	 y=myPosition.y();
	 z=myPosition.z();
	 
	 offset=0;
	 // This will disappear when Andre has applied his fix
	 zsign=1;
	 
	 if(z>0)
	 {
	    if(x>0&&y>0)
	       offset=1;
	    else  if(x<0&&y>0)
	       offset=2;
	    else if(x>0&&y<0)
	       offset=0;
	    else if (x<0&&y<0)
	       offset=3;
	    zsign=1;
	 }
	 else
	 {
	    if(x>0&&y>0)
	       offset=3;
	    else if(x<0&&y>0)
	       offset=2;
	    else if(x>0&&y<0)
	       offset=0;
	    else if(x<0&&y<0)
	       offset=1;
	    zsign=-1;
	 }
	 std::vector<GlobalPoint> corners;
	 corners.clear();
	 corners.resize(8);
	 for(unsigned ic=0;ic<4;++ic)
	 {
	    corners[ic]=getGeometry(mycellID)->getCorners()[(unsigned)((zsign*ic+offset)%4)];
	    corners[4+ic]=getGeometry(mycellID)->getCorners()[(unsigned)(4+(zsign*ic+offset)%4)];
	 }
	 
	 for (short i=0; i < 4 ; ++i)
	 {
	    A = HepPoint3D(corners[i%4].x(),corners[i%4].y(),corners[i%4].z());
	    B = HepPoint3D(corners[(i+1)%4].x(),corners[(i+1)%4].y(),corners[(i+1)%4].z());
	    C = HepPoint3D(corners[4+(i+1)%4].x(),corners[4+(i+1)%4].y(),corners[4+(i+1)%4].z());
	    HepPlane3D plane(A,B,C);
	    plane.normalize();
	    double distance = plane.distance(point);
	    if (corners[0].z()<0.) distance=-distance;
	    SS.push_back(distance);
	 }
	 
	 // Only one move in necessary direction
	 
	 const bool yout ( 0 > SS[0]*SS[2] ) ;
	 const bool xout ( 0 > SS[1]*SS[3] ) ;
	 
	 if( yout || xout )
	 {
	    const int ydel ( !yout ? 0 :  ( 0 < SS[0] ? -1 : 1 ) ) ;
	    const int xdel ( !xout ? 0 :  ( 0 < SS[1] ? -1 : 1 ) ) ;
	    const unsigned int ix ( mycellID.ix() + xdel ) ;
	    const unsigned int iy ( mycellID.iy() + ydel ) ;
	    const unsigned int iz ( mycellID.zside()     ) ;
	    if( EEDetId::validDetId( ix, iy, iz ) ) 
	       mycellID = EEDetId( ix, iy, iz ) ;
	 }
  
	 return mycellID;
      }
      return DetId(0);
   }
   catch ( cms::Exception &e ) 
   { 
      return DetId(0);
   }
}
