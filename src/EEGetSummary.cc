#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "Geometry/EcalAlgo/interface/EcalEndcapGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "FWCore/Utilities/interface/Exception.h"
#include <CLHEP/Geometry/Point3D.h>
#include <CLHEP/Geometry/Plane3D.h>

#include <Math/Transform3D.h>
#include <Math/EulerAngles.h>

typedef CaloCellGeometry::CCGFloat CCGFloat ;
typedef CaloCellGeometry::Pt3D     Pt3D     ;
typedef CaloCellGeometry::Pt3DVec  Pt3DVec  ;
typedef HepGeom::Plane3D<CCGFloat> Pl3D     ;

typedef CaloCellGeometry::Tr3D     Tr3D     ;

typedef HepGeom::Vector3D<double> DVec3D ;
typedef HepGeom::Plane3D<double>  DPlane3D ;
typedef HepGeom::Point3D<double>  DPt3D ;


void
EcalEndcapGeometry::getSummary( CaloSubdetectorGeometry::TrVec&  tVec ,
				CaloSubdetectorGeometry::IVec&   iVec ,   
				CaloSubdetectorGeometry::DimVec& dVec   )  const
{
//hack for producing default EE geometry in db format. 

   dVec.reserve( numberOfShapes()*numberOfParametersPerShape() ) ;

   // std::cout<<"nshape is "<<numberOfShapes()<<std::endl; //always 1 for EE
   const double L  ( 2*parVecVec().front()[0] ) ;
   const double a  ( 2*parVecVec().front()[3] ) ;
   const double A  ( 2*parVecVec().front()[7] ) ;


   dVec.push_back( L/2. ) ;
   dVec.push_back( atan( (A-a)/(sqrt(2.)*L) ) ) ;
   dVec.push_back( M_PI/4. ) ;
   dVec.push_back( a/2.   ) ;
   dVec.push_back( a/2.   ) ;
   dVec.push_back( a/2.   ) ;
   dVec.push_back( 0.      ) ;
   dVec.push_back( A/2.   ) ;
   dVec.push_back( A/2.   ) ;
   dVec.push_back( A/2.   ) ;
   dVec.push_back( 0.      ) ;

//OK Not because below part
 //    for( ParVecVec::const_iterator ivv ( parVecVec().begin() ) ; ivv != parVecVec().end() ; ++ivv )
//     {
//        const ParVec& pv ( *ivv ) ;
//        for( ParVec::const_iterator iv ( pv.begin() ) ; iv != pv.end() ; ++iv )
//        {
//           dVec.push_back( *iv ) ;
//        }
//     }

   const unsigned int nC ( m_validIds.size() ) ; //#of crystals. 

   tVec = TrVec( nC*6, 0. ) ;  //#position pars 3 for Translation, 3 for rotation
   iVec = IVec(  nC  , 0  ) ;  //# DetIDs //
   unsigned int it ( 0 ) ;
   for( uint32_t i ( 0 ) ; i != nC ; ++i, it+=6 )
   {
      iVec[ i      ] = 0 ;  // why it's 0 for segmentation geometry? 
//      iVec[ i + nC ] = 1 ;  



      // first we get the existing crystal cell geometry object
      const CaloCellGeometry* ptr ( cellGeomPtr( i ) ) ;
      assert( 0 != ptr ) ;
      const TruncatedPyramid* tptr ( dynamic_cast<const TruncatedPyramid*>(ptr) ) ;

      // get local coord of shape 1
      Pt3DVec lc     ( 8, Pt3D(0.,0.,0.) ) ; //8 of 3D points specified below, initialize to 0
      Pt3D    lFront ( 0,0,0 ) ;
      TruncatedPyramid::localCorners( lc, &dVec.front(), lFront ) ; //local corners for pyramaid shape from first element of dVec, lFront point is changed to center of front. 
      const Pt3D  lBack   ( 0.25*(lc[4]+lc[5]+lc[6]+lc[7]) ) ;
      const DPt3D dlCorn  ( lc[0].x(), lc[0].y(), lc[0].z() ) ;
      const DPt3D dlFront ( lFront.x(), lFront.y(), lFront.z() ) ;
      const DPt3D dlBack  ( lBack.x(), lBack.y(), lBack.z() ) ;

      const double dFtoC ( ( dlCorn-dlFront ).mag() ) ;
      const double dFtoB ( ( dlBack-dlFront ).mag() ) ;
      const double aCtoB ( ( dlCorn-dlFront ).angle( dlBack-dlFront ) ) ;

      // the front face center
      const GlobalPoint& p ( ptr->getPosition() ) ;
      const DPt3D        dgFront ( p.x(), p.y(), p.z() ) ;
      
      // the corners of the existing crystals in global coordinates
      const CaloCellGeometry::CornersVec& co ( ptr->getCorners() ) ;
      const DPt3D gCorn ( co[0].x(), co[0].y(), co[0].z() ) ;
      DPt3D dgCorn ( gCorn ) ;

      // the centers of all faces must lie along the axis of existing crystal
      const double fr1ToBk1 ( sqrt( L*L + (A-a)*(A-a)/2. ) ) ;
      const GlobalPoint pBack ( p + fr1ToBk1*tptr->axis() ) ;
      const DPt3D gBack ( pBack.x(), pBack.y(), pBack.z() ) ;
      const DPt3D dgBack ( dgFront + dFtoB*( gBack-dgFront ).unit() ) ;

      //-----------------------------
      // this section is to deal with a precision check in Tr3D's constructor
      // that has to do with precision. It prints an error if the angle
      // between two vectors defined by the initial 3 points is preserved
      // to within 1 micro-radian. Hence this angle is checked and
      // a rotation made to correct it before passing the points
      // to the constructor.

      double aCtoBg ( ( dgCorn-dgFront).angle( dgBack-dgFront ) ) ;
      double dangle ( aCtoB - aCtoBg ) ;

      if( 1.e-6 < fabs(dangle) )//guard against precision problems
      {
	 const DPlane3D dgPl ( dgFront, dgCorn, dgBack ) ;
	 const DPt3D    dp2  ( dgFront + dgPl.normal().unit() ) ;

	 dgCorn = ( dgFront + HepGeom::Rotate3D( -dangle, dgFront, dp2 )*
		    DVec3D( dgCorn - dgFront ) ) ;

	 aCtoBg = (dgCorn-dgFront).angle( dgBack-dgFront ) ;
	 dangle = ( aCtoB - aCtoBg ) ;
	 if( 1.e-6<dangle) std::cout<<"After Fix Dangle="<<dangle<<std::endl;
      }
       //-----------------------------

      // the transform
//      Tr3D tr ;
//      const GlobalPoint& gp ( ptr->getPosition() ) ; 
//      tr = HepGeom::Translate3D( gp.x(), gp.y(), gp.z() ) ;

      const Tr3D tr ( dlFront, dlBack, dlCorn,
		      dgFront, dgBack, dgCorn ) ;

      // translation components
      const CLHEP::Hep3Vector  tt ( tr.getTranslation() ) ;
      tVec[ it + 0 ] = tt.x() ;
      tVec[ it + 1 ] = tt.y() ;
      tVec[ it + 2 ] = tt.z() ;
//      std::cout<<"i="<<i<<", Front translation = "<<tt<<std::endl;

      // rotation components
      const CLHEP::HepRotation rr ( tr.getRotation() ) ;
      const ROOT::Math::Transform3D rot ( rr.xx(), rr.xy(), rr.xz(), tt.x(),
					  rr.yx(), rr.yy(), rr.yz(), tt.y(),
					  rr.zx(), rr.zy(), rr.zz(), tt.z()  ) ;

      ROOT::Math::EulerAngles ea ;
      rot.GetRotation( ea ) ;
//      std::cout<<"Front rotation = "<<ea<<std::endl;
      tVec[ it + 3 ] = ea.Phi() ;
      tVec[ it + 4 ] = ea.Theta() ;
      tVec[ it + 5 ] = ea.Psi() ;

   }
}
