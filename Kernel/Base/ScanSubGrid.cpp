#include "ScanSubGrid.h"

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif

CScanSubGrid::~CScanSubGrid()
{
    m_Pts.clear();
}

void CScanSubGrid::SetbFlagPts(bool flag)
{
	//for(int p=0; p<(int)m_Pts.size(); p++)
	//	m_Pts[p]->m_bFlag = flag;
}

void CScanSubGrid::SetPtsBoundingBox()
{
	int nStep = 1;
    int size = m_Pts.size();
	if (m_Pts.size() > 1000) 
		nStep = (int)(size / 1000.);

	//m_PtsBoundingBox.Init();
	for (int p = 0; p < size; p += nStep)
        m_PtsBoundingBox.Update(m_Pts[p].x(), m_Pts[p].y(), m_Pts[p].z());
}
