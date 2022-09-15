#include "ScanGrid.h"

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif


CScanGrid::~CScanGrid()
{
	DeleteGridData();
}

void CScanGrid::DeleteGridData()
{
	if (!m_pSubGrid)
		return;

	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
				delete m_pSubGrid[i][j][k];
			delete[] m_pSubGrid[i][j];
		}
		delete[] m_pSubGrid[i];
	}
	delete[] m_pSubGrid;
	m_pSubGrid = NULL;
}


void CScanGrid::MakeScanSubGrid(int nSubMax, int nz)
{
	DeleteGridData();

	double lenX = m_GridBox.Dx();
	double lenY = m_GridBox.Dy();
	double lenZ = m_GridBox.Dz();

	if (lenX >= lenY && lenX >= lenZ)
	{
		m_NX = nSubMax;
		m_NY = max((int)ceil(nSubMax * (lenY / lenX)), 1);
		m_NZ = max((int)ceil(nSubMax * (lenZ / lenX)), 1);
	}
	else if (lenY >= lenZ)
	{
		m_NX = max((int)ceil(nSubMax * (lenX / lenY)), 1);
		m_NY = nSubMax;
		m_NZ = max((int)ceil(nSubMax * (lenZ / lenY)), 1);
	}
	else
	{
		m_NX = max((int)ceil(nSubMax * (lenX / lenZ)), 1);
		m_NY = max((int)ceil(nSubMax * (lenY / lenZ)), 1);
		m_NZ = nSubMax;
	}
	if (nz > 0)
		m_NZ = nz;

	m_pSubGrid = new CScanSubGrid***[m_NX];
	for (int i = 0; i < m_NX; i++)
	{
		m_pSubGrid[i] = new CScanSubGrid**[m_NY];
		for (int j = 0; j < m_NY; j++)
		{
			m_pSubGrid[i][j] = new CScanSubGrid*[m_NZ];
			for (int k = 0; k < m_NZ; k++)
			{
				m_pSubGrid[i][j][k] = new CScanSubGrid;
			}
		}
	}
}

void CScanGrid::AddToSubGrid(Eigen::Vector3d* pSP,int nSubMax)
{
	if (!m_pSubGrid)
		MakeScanSubGrid(nSubMax);

	int nSx, nSy, nSz;
	GetPos(pSP,nSx,nSy,nSz);

	m_pSubGrid[nSx][nSy][nSz]->AddPoint(pSP);
}

void CScanGrid::GetPos(Eigen::Vector3d* pt, int& nx, int& ny, int& nz) 
{
	nx = (int)floor((pt->x() - m_GridBox.Min.x()) / m_GridBox.Dx() * m_NX);
	ny = (int)floor((pt->y() - m_GridBox.Min.y()) / m_GridBox.Dy() * m_NY);
	nz = (int)floor((pt->z() - m_GridBox.Min.z()) / m_GridBox.Dz() * m_NZ);

	if (nx < 0) nx = 0;
	else if (nx >= m_NX - 1) nx = m_NX - 1;
	if (ny < 0) ny = 0;
	else if (ny >= m_NY - 1) ny = m_NY - 1;
	if (nz < 0) nz = 0;
	else if (nz >= m_NZ - 1) nz = m_NZ - 1;
}


void CScanGrid::UpdateNumberOfGridPts()
{
	m_NPts = 0;

	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
				m_NPts += (int)m_pSubGrid[i][j][k]->GetPts()->size();
		}
	}
}

void CScanGrid::SetPtsBoundingBox()
{
	/*m_PtsBoundingBox.Init();*/

	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
			{
				m_pSubGrid[i][j][k]->SetPtsBoundingBox();
				m_PtsBoundingBox.Union(m_pSubGrid[i][j][k]->GetPtsBoundingBox());
			}
		}
	}
}

void CScanGrid::SetbFlagPts(bool flag)
{
	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
			{
				m_pSubGrid[i][j][k]->SetbFlagPts(flag);
			}
		}
	}
}

void CScanGrid::InitGridPtsList()
{
	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
			{
				m_pSubGrid[i][j][k]->GetPts()->clear();
			}
		}
	}
}