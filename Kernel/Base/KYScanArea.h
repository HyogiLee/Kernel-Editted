#pragma once
#include <vector>
//#include "SGAttrib.h"
#include "SGScanArea.h"
#include "open3d/Open3D.h"
enum EKYScanAnalysisType
{
	SCAN_ANALYSIS_FLATNESS,
	SCAN_ANALYSIS_CYLINDER,
	SCAN_ANALYSIS_CONE,
	SCAN_ANALYSIS_RING,
	SCAN_ANALYSIS_SPHERE,
	SCAN_ANALYSIS_UNDEFINED = 999,
};

class CSGObjectOwner {
public:
    
};


class CKYScanArea : CSGObjectOwner 
{
public:
	CKYScanArea();
	virtual ~CKYScanArea();

public:
	EKYScanAnalysisType m_AnalysisType;
    Eigen::Vector3d m_ptMin, m_ptMax;
    int m_Nx;  // Nz
    int m_Ny;  // Nr

public:
	void SetAnalysisType(EKYScanAnalysisType type) { m_AnalysisType = type; }
	EKYScanAnalysisType GetAnalysisType() { return m_AnalysisType; }

	void SetMinPt(Eigen::Vector3d pt) { m_ptMin = pt; }
    void SetMaxPt(Eigen::Vector3d pt) { m_ptMax = pt; }
    Eigen::Vector3d GetMinPt() { return m_ptMin; }
    Eigen::Vector3d GetMaxPt() { return m_ptMax; }

	void SetNx(int nx) { m_Nx = nx; }
	void SetNy(int ny) { m_Ny = ny; }
	int GetNx() { return m_Nx; }
	int GetNy() { return m_Ny; }

	int GetRefN();
	std::vector<Eigen::Vector3d> GetPartialPts(int nStep = 100);
    static CKYScanArea* MakeArea(std::vector<Eigen::Vector3d>* pts,double vRefSize);

    protected:
	CSGScanAreaAttrib* m_pScanAreaAttrib;

public:	
	CSGScanAreaAttrib* GetScanAreaAttrib() { return m_pScanAreaAttrib; }
};