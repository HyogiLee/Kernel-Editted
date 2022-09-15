#pragma once
#include "open3d/Open3D.h"
/*
enum KyPrimEnum {
    KyPrim_LINE2 = 0,
	KyPrim_LINE3 = 1,
    KyPrim_CIRCLE2 = 2,
	KyPrim_CIRCLE3 = 3,
    KyPrim_ELLIPSE2 = 4,
	KyPrim_ELLIPSE3 = 5,
	KyPrim_PLANE = 6,
	KyPrim_CYLINDER = 7,
	KyPrim_CONE = 8,
	KyPrim_SPHERE = 9,
	KyPrim_UNDEF = 9999
};
*/

class KyPrimFit
{
public:
    bool MakePlane(const std::vector<Eigen::Vector3d>& arrPts,
                          Eigen::Vector3d& ptOrg,
                          Eigen::Vector3d& vNor);

};


class KyPrincipalAxes
{
public:
	Eigen::Vector3d GetCenter()	const { return m_Center; }
	Eigen::Vector3d GetXAxis()	const { return m_XAxis; }
	Eigen::Vector3d GetYAxis()	const { return m_YAxis; }
	Eigen::Vector3d GetZAxis()	const { return m_ZAxis; }

    void Set(const std::vector<Eigen::Vector3d> pts);

public:
    KyPrincipalAxes(const std::vector<Eigen::Vector3d> pts);
	KyPrincipalAxes();

	~KyPrincipalAxes();
protected:
	void Init0();
	void Clear();
private:
	Eigen::Vector3d m_Center;
	Eigen::Vector3d m_XAxis;
	Eigen::Vector3d m_YAxis;
	Eigen::Vector3d m_ZAxis;

	// a matrix of the eigenvalue problem
	double** m_EigenvalueProblem;
	// for efficiency reasons parts of the eigenvalue problem are computed separately
	double** m_EigenvalueProblemDiag;
	double** m_Eigenvectors;
	double*  m_Eigenvalues;
};
 
