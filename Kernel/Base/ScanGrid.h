#pragma once
#include <vector>
#include "V3D.h"
#include "ScanSubGrid.h"

class CScanGrid
{
public:
    CScanGrid() {
        m_pSubGrid = nullptr; 
	m_NX = m_NY = m_NZ = 0;
	}
	~CScanGrid();

	bool m_bFlag;

private:
	CScanSubGrid**** m_pSubGrid;

	V3D m_GridBox;
    V3D m_PtsBoundingBox;

	short int m_NX, m_NY, m_NZ;

	int m_NPts;

public:
	void DeleteGridData();
	void MakeScanSubGrid(int nSub, int nz = 0);
	void AddToSubGrid(Eigen::Vector3d* pSP,int nMaxSub);
    void GetPos(Eigen::Vector3d* pt, int& nx, int& ny, int& nz);

	void SetGridBox(V3D box) { m_GridBox = box; }
	void SetSubGrid(CScanSubGrid****& pSubGrid) { m_pSubGrid = pSubGrid; }
	void SetN(int nx, int ny, int nz) { m_NX = nx; m_NY = ny; m_NZ = nz; }

	int GetNX() { return m_NX; }
	int GetNY() { return m_NY; }
	int GetNZ() { return m_NZ; }
	int GetNPts() { return m_NPts; }

	CScanSubGrid* GetGrid(int x, int y, int z) { return m_pSubGrid[x][y][z]; }
    V3D GetGridBox() { return m_GridBox; }

	void SetPtsBoundingBox(V3D box) { m_PtsBoundingBox = box; }
    V3D GetPtsBoundingBox() { return m_PtsBoundingBox; }
	void SetPtsBoundingBox();

	void UpdateNumberOfGridPts();
	
	void SetbFlagPts(bool flag);
	void InitGridPtsList();
};
