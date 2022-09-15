#include "KYScanArea.h"


//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif


CKYScanArea::CKYScanArea()
{
    m_pScanAreaAttrib = new CSGScanAreaAttrib;
    m_pScanAreaAttrib->SetAttribOwner(this);

	m_Nx = m_Ny = 0;
	m_AnalysisType = SCAN_ANALYSIS_UNDEFINED;
}

CKYScanArea::~CKYScanArea() 
{
	if (m_Nx == 0)
		return;
}

std::vector<Eigen::Vector3d> CKYScanArea::GetPartialPts(int nStep /* =100 */)
{
	std::vector<Eigen::Vector3d> pts;
    CSGScanAreaAttrib::GetPartialScanPts(GetScanAreaAttrib(), &pts, nStep);

	return pts;
}

CKYScanArea* CKYScanArea::MakeArea(std::vector<Eigen::Vector3d>* pts, double vRefSize) {
    CKYScanArea* pArea = new CKYScanArea;
    CSGScanAreaAttrib* pAttrib = pArea->GetScanAreaAttrib();
    pAttrib->MakeGrid(pts, vRefSize);

    return pArea;
}

//���� �м��� ��, ��ҿ� ���Ե� �ּ� ����Ʈ ����
int CKYScanArea::GetRefN()
{
	int NMax = max(m_Nx, m_Ny);

	double v = GetScanAreaAttrib()->GetNPts() / (NMax * 60.0) + 0.5;

	int nRef = (int)pow(v, 0.4);

	if (nRef < 1)
		nRef = 1;

	return nRef;
}

