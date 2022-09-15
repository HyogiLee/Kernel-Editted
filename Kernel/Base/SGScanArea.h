#pragma once
#include <tuple>
#include <vector>
#include "open3d/Open3D.h"
#include "V3D.h"
#include "ScanGrid.h"

enum class EPointSetDisplayMode
{
	ePointSetDisplayUndef = -1,
	ePointSetDisplayAlwaysAll = 0, 
	ePointSetDisplayLOD = 1,
};
enum class EPointSetUpdateType
{
	ePointSetNotInitialized = 0,
	ePointSetInitialized = 1,
	ePointSetUpdatePicking = 2,
	ePointSetUpdateAll = 3,
};

class CSGObjectOwner;
class CSGScanAreaAttrib
{
public:
	CSGScanAreaAttrib();
	virtual ~CSGScanAreaAttrib();


private:
	int m_NPts;
    CSGObjectOwner* m_pOwner;
	CScanGrid**** m_pGrid;
    V3D m_BV;
	V3D m_Box;  //���� ������ ���� ��� �ڽ�
    V3D m_lmtBox;  //������ ���� �뵵, ���ο������� ����� ����Ʈ�� ����
	int m_NX, m_NY, m_NZ, m_NMax;
	double m_MaxLen;
	int m_NSubMax;


public:
	int GetNX() { return m_NX; }
	int GetNY() { return m_NY; }
	int GetNZ() { return m_NZ; }

	int GetNPts() { return m_NPts; }
	void SetNPts(int nPts) { m_NPts = nPts; }
    void SetBoundingBox(V3D box) { m_Box = box; }
    void SetLimitedBoundingBox(V3D box) { m_lmtBox = box; }
	void SetMaxLength(double vLen) { m_MaxLen; }
    V3D GetBoundBox() { return m_Box; }
    V3D GetLimitedBoundBox() { return m_lmtBox; }
	double GetMaxLength() { return m_MaxLen; }
	CScanGrid* GetGrid(int x, int y, int z) { 
		if (!m_pGrid)
			return NULL;
		return m_pGrid[x][y][z]; }

	//GRID
	
	void MakeGrid(vector<Eigen::Vector3d>* ptsScan, double vSize=0);

	//�߰��� ���� ���Ͽ� ���� ���� ������Ʈ
    void UpdateGrid(vector<Eigen::Vector3d>* ptsScan,std::vector<float>& rgbs);
    void CalcBoundingBox(vector<Eigen::Vector3d>* ptsScan);
	void MakeMainGrid(double vSize=0);
	void MakeSubGrid();
	void DeleteGrid();
	void InitGridPtsList();
    void MakeScanPoint(vector<Eigen::Vector3d>* ptsScan,
                           std::vector<float>& rgbs);
    void SetScanPtsToGrid(vector<Eigen::Vector3d>* pScanPts);
    void GetPos(Eigen::Vector3d* pt, int& nx, int& ny, int& nz);
	void AddToGrid();
    void AddToGrid(Eigen::Vector3d* pSP);
	void UpdateGridPtsData();
    void SetAttribOwner(CSGObjectOwner* pOwner) { m_pOwner = pOwner; }
	Eigen::Vector3d GetCenter();

	//After Moving
	void ClearWOPtsAll();
	void CalcBoundingBoxUsingAllPts();
    void CalcBoundingBoxUsingPts(std::vector<Eigen::Vector3d*>& pScanPts);
    V3D GetBoundingBox(CSGScanAreaAttrib* pArea, int nStep);
    static void GetPartialScanPts(CSGScanAreaAttrib* pArea,
                                              vector<Eigen::Vector3d>* pts,
                                              int nStep); 
	static void GetPartialGridScanPts(CScanGrid* pGrid,
                                      vector<Eigen::Vector3d>* pts,
                                      int nStep);
    public:	
	BOOL CalcBV();

// for glsl shader
private:
//	std::vector<CScanPoint*> m_allPts;
	EPointSetUpdateType m_eNeedToUpdate;
	EPointSetDisplayMode m_eDisplayMode; 

public:
//	std::vector<CScanPoint*>& GetAllScanPoints() { return m_allPts; }
	void SetUpdateFlag(EPointSetUpdateType type);
	void ClearUpdateFlag() { m_eNeedToUpdate = EPointSetUpdateType::ePointSetInitialized; }
	bool IsUpdated() { return m_eNeedToUpdate > EPointSetUpdateType::ePointSetInitialized; }
	EPointSetUpdateType GetUpdateFlag() { return m_eNeedToUpdate; }

	void SetDisplayMode(EPointSetDisplayMode mode) { m_eDisplayMode = mode; }
	void ClearDisplayMode() { m_eDisplayMode = EPointSetDisplayMode::ePointSetDisplayUndef; }
	EPointSetDisplayMode GetDisplayMode() { return m_eDisplayMode; }

	int m_initialStep = 0;

};
