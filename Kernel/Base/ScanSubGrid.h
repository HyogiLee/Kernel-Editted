#pragma once
#include <vector>
//#include <open3d/Open3D.h>
#include "V3D.h"
class CScanSubGrid
{
public:
	CScanSubGrid() {}
	~CScanSubGrid();
	bool m_bFlag;
	bool m_bFlag2;

private:
    std::vector<Eigen::Vector3d> m_Pts;
    V3D m_PtsBoundingBox;

public:
    void AddPoint(Eigen::Vector3d* pSP) { m_Pts.push_back(*pSP); }
    std::vector<Eigen::Vector3d>* GetPts() { return &m_Pts; }
    void SetPts(std::vector<Eigen::Vector3d>* pts) { m_Pts = *pts; }

	V3D GetPtsBoundingBox() { return m_PtsBoundingBox; }
	void SetPtsBoundingBox();
	void SetbFlagPts(bool flag);
};
