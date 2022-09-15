#pragma once
#include "open3d/Open3D.h"
//#include <KyBase/KyBaseExportDll.h>
//#include <KyBase/KyVARRAY.h>
//#include <KyBase/KyTMatrix.h>


class KyGeoICP {
public:
	// Specify the source and target data sets.
	void SetSource(const KyP3Ds& source);
	void SetTarget(const KyP3Ds& target);
	void Set(const KyP3Ds& source, const KyP3Ds& target);

	// Get a rigid transformation matrix
	//	- an iterative closest points algorithm
	bool GetRigidTM(KyTMatrix& tm, const int& iter = 50, const double& tol = 0.01, bool errorMetric = true,bool bMappingFixX=false,bool bMappingFixY=false,bool bMappingFixZ=false);

	const KyTMatrix& GetRigidTM() { return m_RigidTM; }

	double GetError() const { return m_Error; }
protected:

	void PreProcess();
	// Calc. translation t and rotation R 
	//		that minimizes the sum of the squared error
	KyTMatrix CalcRigidTM(bool bMappingFixX = false, bool bMappingFixY = false, bool bMappingFixZ = false) const;
	// Transform a source points set 
	// using a rigid transformation matrix
	bool TransformSource(const KyTMatrix& tm);
	// Find unit vectors 
	// which is perpendicular to this on and to each other.
	void Perpendiculars(const KyV3D& v0, 
						KyV3D& v1, KyV3D& v2, double& theta) const;
	
	double EstiError(bool mode =  true);
public:
	KyGeoICP();
	~KyGeoICP();

private:
	KyP3Ds m_Source;
	KyP3Ds m_Target;
	double m_Error;

	KyTMatrix m_RigidTM;
};