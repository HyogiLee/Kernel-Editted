#pragma once
#include "open3d/Open3D.h"
#include "KyMath.h"
//#include <KyBase/KyBaseExportDll.h>
//#include <KyBase/KyVARRAY.h>
//#include <KyBase/KyTMatrix.h>


class KyGeoICP {
public:
	// Specify the source and target data sets.
	void SetSource(const std::vector<Eigen::Vector3d>& source);
	void SetTarget(const std::vector<Eigen::Vector3d>& target);
	void Set(const std::vector<Eigen::Vector3d>& source, const std::vector<Eigen::Vector3d>& target);

	// Get a rigid transformation matrix
	//	- an iterative closest points algorithm
	bool GetRigidTM(Eigen::Matrix4d& tm, const int& iter = 50, const double& tol = 0.01, bool errorMetric = true,bool bMappingFixX=false,bool bMappingFixY=false,bool bMappingFixZ=false);

	const Eigen::Matrix4d& GetRigidTM() { return m_RigidTM; }

	double GetError() const { return m_Error; }
protected:

	void PreProcess();
	// Calc. translation t and rotation R 
	//		that minimizes the sum of the squared error
	Eigen::Matrix4d  CalcRigidTM(bool bMappingFixX = false, bool bMappingFixY = false, bool bMappingFixZ = false) const;
	// Transform a source points set 
	// using a rigid transformation matrix
	bool TransformSource(const Eigen::Matrix4d& tm);
	// Find unit vectors 
	// which is perpendicular to this on and to each other.
	void Perpendiculars(const Eigen::Vector3d& v0,
		Eigen::Vector3d& v1, Eigen::Vector3d& v2, double& theta) const;
	
	double EstiError(bool mode =  true);
public:
	KyGeoICP();
	~KyGeoICP();

private:
	std::vector<Eigen::Vector3d> m_Source;
	std::vector<Eigen::Vector3d> m_Target;
	double m_Error=0.0;

	Eigen::Matrix4d  m_RigidTM;
};