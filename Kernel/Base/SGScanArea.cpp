#include "SGScanArea.h"
#include <chrono>

//#ifdef _DEBUG
//#define new DEBUG_NEW
//#endif


CSGScanAreaAttrib::CSGScanAreaAttrib()
{
	m_pGrid = NULL;
	m_eNeedToUpdate = EPointSetUpdateType::ePointSetNotInitialized;
	m_eDisplayMode = EPointSetDisplayMode::ePointSetDisplayAlwaysAll;

	m_NX = m_NY = m_NZ = 0;
	m_NPts = 0;
    m_lmtBox.Min = m_lmtBox.Max = {0, 0, 0};
}

CSGScanAreaAttrib::~CSGScanAreaAttrib()
{
	DeleteGrid();
}

void CSGScanAreaAttrib::DeleteGrid()
{
	if (!m_pGrid)
		return;

	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
				delete m_pGrid[i][j][k];
			delete[] m_pGrid[i][j];
		}
		delete[] m_pGrid[i];
	}

	delete[] m_pGrid;
	m_pGrid = NULL;
}

void CSGScanAreaAttrib::InitGridPtsList()
{
	if (!m_pGrid)
		return;

	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
				m_pGrid[i][j][k]->InitGridPtsList();
		}
	}
}

V3D CSGScanAreaAttrib::GetBoundingBox(CSGScanAreaAttrib* pArea, int nStep) 
{
    V3D mb;
    if (pArea->GetNPts() == 0) return mb;

    if (nStep < 1) nStep = 1;

    int NX = pArea->GetNX();
    int NY = pArea->GetNY();
    int NZ = pArea->GetNZ();

    int nCnt = 0;
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                CScanGrid* pGrid = pArea->GetGrid(i, j, k);

                int nx = pGrid->GetNX();
                int ny = pGrid->GetNY();
                int nz = pGrid->GetNZ();

                for (int ii = 0; ii < nx; ii++) {
                    for (int jj = 0; jj < ny; jj++) {
                        for (int kk = 0; kk < nz; kk++) {
                            CScanSubGrid* pSubGrid = pGrid->GetGrid(ii, jj, kk);

                            for (int p = 0; p < (int)pSubGrid->GetPts()->size();
                                 p++) {
                                if (nCnt % nStep == 0)
                                    mb.Update(pSubGrid->GetPts()->at(p).x(),
                                              pSubGrid->GetPts()->at(p).y(),
                                              pSubGrid->GetPts()->at(p).z());
                                nCnt++;
                            }
                        }
                    }
                }
            }
        }
    }

    return mb;
}

BOOL CSGScanAreaAttrib::CalcBV()
{
	/*if (m_bValidBV)
		return TRUE;*/

	//m_BV.Init();

	int nStep = int(m_NPts / 100000);
	if (nStep < 1) nStep = 1;

	vector<Eigen::Vector3d> pts;
	m_BV = GetBoundingBox(this, nStep);

	//m_bValidBV = TRUE;
	return TRUE;
}

void CSGScanAreaAttrib::MakeGrid(vector<Eigen::Vector3d>* ptsScan, double vSize)
{
	if (m_pGrid)
		return;

	m_NPts = (int)ptsScan->size();

	//���ڽ� ũ�� ���
	CalcBoundingBox(ptsScan);

	//�׸��� ũ�� ���� �� �޸� �Ҵ�
	MakeMainGrid(vSize);

	//��ĵ ����Ʈ ���� �� �׸��忡 �ֱ�
	std::vector<float> rgbs;
	MakeScanPoint(ptsScan, rgbs);

	//���� ���� ������Ʈ
	UpdateGridPtsData();

	CalcBV();
}

void CSGScanAreaAttrib::UpdateGrid(vector<Eigen::Vector3d>* ptsScan, std::vector<float>& rgbs) 
{
	
	m_NPts = (int)ptsScan->size();

	V3D box0;
	CalcBoundingBox(ptsScan);
    m_Box.Update(box0.Min.x(), box0.Min.y(), box0.Min.z());
    m_Box.Update(box0.Max.x(), box0.Max.y(), box0.Max.z());

	if (m_pGrid)
	{
		for (int i = 0; i < m_NX; i++)
		{
			for (int j = 0; j < m_NY; j++)
			{
				for (int k = 0; k < m_NZ; k++)
					delete m_pGrid[i][j][k];

				delete[] m_pGrid[i][j];
			}
			delete[] m_pGrid[i];
		}
		delete[] m_pGrid;
		m_pGrid = NULL;
	}

	MakeMainGrid();

	AddToGrid();

	MakeScanPoint(ptsScan, rgbs);

	UpdateGridPtsData();
}

void CSGScanAreaAttrib::MakeScanPoint(vector<Eigen::Vector3d>* ptsScan, std::vector<float>& rgbs) 
{
    int size = ptsScan->size();
    for (int i = 0; i < size; i++)
	{
        Eigen::Vector3d pSP = ptsScan->at(i);
		AddToGrid(&pSP);

		//if (!rgbs.empty())
		//{
// 			unsigned int rgb;
// 			if (rgbs[i] < 1)
// 				rgb = *reinterpret_cast<int*>(&rgbs[i]);
// 			else
// 				rgb = (int)rgbs[i];
// 			byte r, g, b;
// 			r = (rgb >> 16) & 0x0000ff;
// 			g = (rgb >> 8) & 0x0000ff;
// 			b = (rgb) & 0x0000ff;
// 
// //			if (abs(r - g) < 2 && abs(g - b) < 2)
// //				g = b = 0;
			//pSP.SetRGB(rgbs[i]);
		//}
	}
}

void CSGScanAreaAttrib::SetScanPtsToGrid(vector<Eigen::Vector3d>* pScanPts) 
{
    int size = (int)pScanPts->size();
	for (int i = 0; i < size; i++)
		AddToGrid(&pScanPts->at(i));
}

void CSGScanAreaAttrib::AddToGrid()
{
	int NX = this->GetNX();
	int NY = this->GetNY();
	int NZ = this->GetNZ();

	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NY; j++)
		{
			for (int k = 0; k < NZ; k++)
			{
				CScanGrid* pGrid = GetGrid(i, j, k);

				int nx = pGrid->GetNX();
				int ny = pGrid->GetNY();
				int nz = pGrid->GetNZ();

				for (int ii = 0; ii < nx; ii++)
				{
					for (int jj = 0; jj < ny; jj++)
					{
						for (int kk = 0; kk < nz; kk++)
						{
							CScanSubGrid* pSubGrid = pGrid->GetGrid(ii, jj, kk);

							for (int p = 0; p < (int)pSubGrid->GetPts()->size(); p++)
							{
								AddToGrid(&pSubGrid->GetPts()->at(p));
							}
						}
					}
				}
			}
		}
	}
}

void CSGScanAreaAttrib::AddToGrid(Eigen::Vector3d *pSP)
{
	int nx, ny, nz;
	GetPos(pSP, nx, ny, nz);

	m_pGrid[nx][ny][nz]->AddToSubGrid(pSP, m_NSubMax);
}

void CSGScanAreaAttrib::GetPos(Eigen::Vector3d *pt, int& nx, int& ny, int& nz) 
{
	nx = (int)floor((pt->x() - m_Box.Min.x()) / m_Box.Dx() * m_NX);
	ny = (int)floor((pt->y() - m_Box.Min.y()) / m_Box.Dy() * m_NY);
	nz = (int)floor((pt->z() - m_Box.Min.z()) / m_Box.Dz() * m_NZ);

	if (nx < 0) nx = 0;
	else if (nx >= m_NX - 1) nx = m_NX - 1;
	if (ny < 0) ny = 0;
	else if (ny >= m_NY - 1) ny = m_NY - 1;
	if (nz < 0) nz = 0;
	else if (nz >= m_NZ - 1) nz = m_NZ - 1;
}

void CSGScanAreaAttrib::MakeSubGrid()
{
	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
			{
				m_pGrid[i][j][k]->MakeScanSubGrid(m_NSubMax);
			}
		}
	}
}

void CSGScanAreaAttrib::MakeMainGrid(double vSize)
{
	m_NMax = (int)(pow(m_NPts / (27 * 500), 0.333));
	if (m_NMax < 1) m_NMax = 1;

	double vWidth = m_MaxLen / m_NMax;

	double lenX = m_Box.Dx();
	double lenY = m_Box.Dy();
	double lenZ = m_Box.Dz();

	//������ ������ ��ü ���ڰ� ���� �����ȴ�. �� ��� �ִ� ���ڸ� Ű���ش�.
	double nx = (m_MaxLen / max(lenX, vWidth));
	double ny = (m_MaxLen / max(lenY, vWidth));
	double nz = (m_MaxLen / max(lenZ, vWidth));

	m_NMax = int(m_NMax * pow(nx * ny * nz, 0.444));

	//Segmentation�� Ȱ���ϱ� ���� ũ�⸦ �����Ͽ� ���� ������ ����
	if (vSize > 0)
		m_NMax = (int)(m_MaxLen / vSize);
	if (m_NMax < 1) m_NMax = 1;


	if (lenX >= lenY && lenX >= lenZ)
	{
		m_NX = m_NMax;
		m_NY = max((int)ceil(m_NMax * (lenY / lenX)), 1);
		m_NZ = max((int)ceil(m_NMax * (lenZ / lenX)), 1);
	}
	else if (lenY >= lenZ)
	{
		m_NX = max((int)ceil(m_NMax * (lenX / lenY)), 1);
		m_NY = m_NMax;
		m_NZ = max((int)ceil(m_NMax * (lenZ / lenY)), 1);
	}
	else
	{
		m_NX = max((int)ceil(m_NMax * (lenX / lenZ)), 1);
		m_NY = max((int)ceil(m_NMax * (lenY / lenZ)), 1);
		m_NZ = m_NMax;
	}

	int NGrid = m_NX * m_NY * m_NZ;
	m_NSubMax = (int)pow(m_NPts / (NGrid * 500), 0.333);
	if (m_NSubMax < 1)
		m_NSubMax = 1;

	// 	CString strMsg;
	// 	strMsg.Format(_T("m_NMax: %d, m_NSubMax: %d"), m_NMax, m_NSubMax);
	// 	UtilMsg::UpdateReport(strMsg);

	m_pGrid = new CScanGrid * **[m_NX];
	double dX = lenX / m_NX;
	double dY = lenY / m_NY;
	double dZ = lenZ / m_NZ;

	double vCurX = m_Box.Min.x();
	double vCurY = m_Box.Min.y();
	double vCurZ = m_Box.Min.z();

	for (int i = 0; i < m_NX; i++)
	{
		m_pGrid[i] = new CScanGrid * *[m_NY];
		for (int j = 0; j < m_NY; j++)
		{
			m_pGrid[i][j] = new CScanGrid * [m_NZ];
			for (int k = 0; k < m_NZ; k++)
			{
				m_pGrid[i][j][k] = new CScanGrid;

				V3D box;
                box.Min = {vCurX, vCurY, vCurZ};
                box.Max = {vCurX + dX, vCurY + dY, vCurZ + dZ};
				m_pGrid[i][j][k]->SetGridBox(box);
				//				m_pGrid[i][j][k]->SetPos(i, j, k);

				vCurZ += dZ;
			}

			vCurZ = m_Box.Min.z();
			vCurY += dY;
		}
		vCurY = m_Box.Min.y();
		vCurX += dX;
	}
}


void CSGScanAreaAttrib::CalcBoundingBox(vector<Eigen::Vector3d>* ptsScan) 
{
	//m_Box.Init();
    int size = ptsScan->size();
    double val = (double)size / 10000.0 + 1;
	int nRef = int(sqrt(val));

	for (int i = 0; i < size; i += nRef)
	   m_Box.Update(ptsScan->at(i).x(), ptsScan->at(i).y(), ptsScan->at(i).z());


	m_MaxLen = max(m_Box.Dx(), max(m_Box.Dy(), m_Box.Dz()));
	m_Box.Expand((float)(m_MaxLen / 1000));
}

void CSGScanAreaAttrib::GetPartialScanPts(CSGScanAreaAttrib* pArea,vector<Eigen::Vector3d>* pts,int nStep) 
{
    if (pArea->GetNPts() == 0) return;

    if (nStep < 1) nStep = 1;

    int NX = pArea->GetNX();
    int NY = pArea->GetNY();
    int NZ = pArea->GetNZ();

    int nSize = (int)(pArea->GetNPts() / (double)nStep + 0.9999);
    pts->resize(nSize);

    int nCnt = 0;
    int nP = 0;
    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            for (int k = 0; k < NZ; k++) {
                CScanGrid* pGrid = pArea->GetGrid(i, j, k);

                int nx = pGrid->GetNX();
                int ny = pGrid->GetNY();
                int nz = pGrid->GetNZ();

                for (int ii = 0; ii < nx; ii++) {
                    for (int jj = 0; jj < ny; jj++) {
                        for (int kk = 0; kk < nz; kk++) {
                            CScanSubGrid* pSubGrid = pGrid->GetGrid(ii, jj, kk);

                            for (int p = 0; p < (int)pSubGrid->GetPts()->size();
                                 p++) {
                                if (nCnt % nStep == 0) {
                                    pts->at(nP) = pSubGrid->GetPts()->at(p);
                                    nP++;

                                    if (nP == nSize) return;
                                }
                                nCnt++;
                            }
                        }
                    }
                }
            }
        }
    }

    if (nP != nSize) pts->pop_back();
}

void CSGScanAreaAttrib::GetPartialGridScanPts(CScanGrid* pGrid,
                                              vector<Eigen::Vector3d>* pts,
                                              int nStep) 
{
    int nx = pGrid->GetNX();
    int ny = pGrid->GetNY();
    int nz = pGrid->GetNZ();

    int nSize = (int)(pGrid->GetNPts() / (double)nStep + 0.9999);
    pts->resize(nSize);

    int nCnt = 0;
    int nP = 0;
    for (int ii = 0; ii < nx; ii++) {
        for (int jj = 0; jj < ny; jj++) {
            for (int kk = 0; kk < nz; kk++) {
                CScanSubGrid* pSubGrid = pGrid->GetGrid(ii, jj, kk);

                for (int p = 0; p < (int)pSubGrid->GetPts()->size(); p++) {
                    if (nCnt % nStep == 0) {
                        pts->at(nP) = pSubGrid->GetPts()->at(p);
                        nP++;
                    }
                    nCnt++;
                }
            }
        }
    }
}


void CSGScanAreaAttrib::CalcBoundingBoxUsingAllPts()
{
	//m_Box.Init();
	double val = GetNPts() / 10000.0 + 1;
	int nStep = int(sqrt(val));
	if (nStep < 1) nStep = 1;

	vector<Eigen::Vector3d> pts;
	GetPartialScanPts(this, &pts, nStep);
	for (int i = 0; i < (int)pts.size(); i++)
		m_Box.Update(pts[i]);

	m_MaxLen = max(m_Box.Dx(), max(m_Box.Dy(), m_Box.Dz()));
	m_Box.Expand((float)(m_MaxLen / 1000));
}

void CSGScanAreaAttrib::UpdateGridPtsData()
{
	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
			{
				m_pGrid[i][j][k]->UpdateNumberOfGridPts();
				m_pGrid[i][j][k]->SetPtsBoundingBox();
			}
		}
	}
}

void CSGScanAreaAttrib::ClearWOPtsAll()
{
	for (int i = 0; i < m_NX; i++)
	{
		for (int j = 0; j < m_NY; j++)
		{
			for (int k = 0; k < m_NZ; k++)
			{
				CScanGrid* pGrid = this->GetGrid(i, j, k);

				int nx = pGrid->GetNX();
				int ny = pGrid->GetNY();
				int nz = pGrid->GetNZ();

				for (int ii = 0; ii < nx; ii++)
				{
					for (int jj = 0; jj < ny; jj++)
					{
						for (int kk = 0; kk < nz; kk++)
						{
							CScanSubGrid* pSubGrid = pGrid->GetGrid(ii, jj, kk);
							pSubGrid->GetPts()->clear();//���� ����Ʈ�� ������ �ʱ� ���ؼ� ����Ʈ�� �ʱ�ȭ ��
						}
					}
				}

				delete m_pGrid[i][j][k];
			}
			delete[] m_pGrid[i][j];
		}
		delete[] m_pGrid[i];
	}
	delete[] m_pGrid;
	m_pGrid = NULL;
}


Eigen::Vector3d CSGScanAreaAttrib::GetCenter() 
{
	Eigen::Vector3d ptCenter = {0, 0, 0};
	if (m_NPts == 0)
		return ptCenter;

	int NX = this->GetNX();
	int NY = this->GetNY();
	int NZ = this->GetNZ();

	for (int i = 0; i < NX; i++)
	{
		for (int j = 0; j < NY; j++)
		{
			for (int k = 0; k < NZ; k++)
			{
				CScanGrid* pGrid = GetGrid(i, j, k);

				int nx = pGrid->GetNX();
				int ny = pGrid->GetNY();
				int nz = pGrid->GetNZ();

				for (int ii = 0; ii < nx; ii++)
				{
					for (int jj = 0; jj < ny; jj++)
					{
						for (int kk = 0; kk < nz; kk++)
						{
							CScanSubGrid* pSubGrid = pGrid->GetGrid(ii, jj, kk);

							for (int p = 0; p < (int)pSubGrid->GetPts()->size(); p++)
								ptCenter += pSubGrid->GetPts()->at(p);
						}
					}
				}
			}
		}
	}

	ptCenter /= m_NPts;
	return ptCenter;
}

void CSGScanAreaAttrib::SetUpdateFlag(EPointSetUpdateType type)
{
	if (m_eNeedToUpdate > type)
		return;
	else
		m_eNeedToUpdate = type;
}