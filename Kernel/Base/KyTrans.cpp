#include "KyTrans.h"
#include "KyMatching.h"
#include "KyPrimFit.h"
//#include "V3D.h" 

#define OFFSET_ORG 139.96
#define OFFSET_BLOCK 600.37

 

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

KyTrans::KyTrans(void)
{
	m_bMappingFixX = m_bMappingFixY = m_bMappingFixZ = false;
	m_vRmsError = 0;
}

KyTrans::~KyTrans(void)
{
}

void KyTrans::UpdatePoint(std::vector<Eigen::Vector3d>& pts, Eigen::Matrix4d& tMat)
{
	Eigen::Matrix3d R = tMat.block<3, 3>(0, 0);
	for (size_t i = 0; i < pts.size(); i++) {
		//pts[i] = tMat * pts[i];
		pts[i] = R * pts[i];
	}
}

void KyTrans::UpdatePoint(std::vector<Eigen::Vector3f>& pts, Eigen::Matrix4f& tMat)
{
	Eigen::Matrix3f R = tMat.block<3, 3>(0, 0);
	for (size_t i = 0; i < pts.size(); i++) {
		//pts[i] = tMat * pts[i];
		pts[i] = R * pts[i];
	}
}

Eigen::Matrix4d KyTrans::MovePtsNtoN(std::vector<Eigen::Vector3d>& source, std::vector<Eigen::Vector3d>& target,bool bSrcMove)
{
	m_tMat.setIdentity();
	if(source.size() != target.size())
		return m_tMat;
	if(source.size() < 1)
		return m_tMat;

	m_errDx.clear();
	m_errDy.clear();
	m_errDz.clear();

	m_icp.SetSource(source);
	m_icp.SetTarget(target);

	m_icp.GetRigidTM(m_tMat,100,0.001,true,m_bMappingFixX,m_bMappingFixY,m_bMappingFixZ);

	std::vector<Eigen::Vector3d> newSource = source;
	for (size_t i = 0; i < source.size(); i++) {
		//newSource[i] = m_tMat * source[i];
		Eigen::Matrix3d R = m_tMat.block<3, 3>(0, 0);
		//pts[i] = tMat * pts[i];
		newSource[i] = R * source[i];
		Eigen::Vector3d TranslationVec = m_tMat.block<3, 1>(0, 3);
		newSource[i] += TranslationVec;
	}

	if(bSrcMove)
		source = newSource;

	m_vRmsError = 0;
	for(size_t i=0; i<source.size(); i++)
	{
		m_errDx.push_back(target[i].x() - newSource[i].x());
		m_errDy.push_back(target[i].y() - newSource[i].y());
		m_errDz.push_back(target[i].z() - newSource[i].z());
		
		//double dist = target[i].Dist(newSource[i]);
		double dist = KyMath::Dist(target[i], newSource[i]);
		m_vRmsError += dist * dist;
	}
	m_vRmsError /= source.size();
	m_vRmsError = std::sqrt(m_vRmsError);
	/*Translation Matrix*/
	return m_tMat;
}


double KyTrans::GetRmsError()
{
	double eTol = 0.;
	int N = int((m_errDx.size() + m_errDy.size() + m_errDz.size())/3.);
	if(N == 0)
		return -1.;

	for(unsigned int i=0; i < m_errDx.size(); i++)
		eTol += m_errDx[i] * m_errDx[i];
	for(unsigned int i=0; i < m_errDy.size(); i++)
		eTol += m_errDy[i] * m_errDy[i];
	for(unsigned int i=0; i < m_errDz.size(); i++)
		eTol += m_errDz[i] * m_errDz[i];

	return sqrt(eTol/N);
}

double KyTrans::GetPipeRmsError()
{
	double eTol = 0.;
	int N = int((m_errDx.size() + m_errDy.size()) / 2.);// +m_errDz.size());
	if(N == 0)
		return -1.;

	for(unsigned int i=0; i<m_errDx.size(); i++)
		eTol += m_errDx[i]*m_errDx[i];
	for(unsigned int i=0; i<m_errDy.size(); i++)
		eTol += m_errDy[i]*m_errDy[i];
// 	for(unsigned int i=0; i<m_errDz.size(); i++)
// 		eTol += m_errDz[i]*m_errDz[i];

	return sqrt(eTol/N);
}


void KyTrans::MovePtsNtoM(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget,std::vector<int>& indexMsur,std::vector<int>& indexCad,double vMatchingError, bool bOnlyRotateZ)
{
	if(ptsSource.size() == 0 || ptsTarget.size() == 0)
		return;

	KyAutoMatch matchAlgo;
	int i;
	Eigen::Vector3d pt;


	//set target pts data
	int nDest = int(ptsTarget.size());
	matchAlgo.AllocBaseData(nDest);

	double point[3];
	for (i = 0;i < nDest; i++) {
		point[0] = ptsTarget[i].x();
		point[1] = ptsTarget[i].y();
		point[2] = ptsTarget[i].z();
		matchAlgo.SetBaseData(i, i+1, point);
	}

	//set source pts data
	int nSrc =  int(ptsSource.size());
	matchAlgo.AllocMeasures(nSrc);

	for (i = 0;i < nSrc; i++) {
		point[0] = ptsSource[i].x();
		point[1] = ptsSource[i].y();
		point[2] = ptsSource[i].z();
		matchAlgo.SetMeasures(i, i+1, point);
	}	

	double filterBound = vMatchingError;
	matchAlgo.SetDeltaFiltered(filterBound, true);
	matchAlgo.with_Filter_doit(bOnlyRotateZ);

	int nMapped =  matchAlgo.GetNFiltered();
	double **TransformedPts = matchAlgo.GetTransformed();
	int *mIndex = matchAlgo.GetIdxM();
	int *cIndex = matchAlgo.GetIdxO();

	ptsSource.clear();
	for(i=0; i<nMapped; i++)
	{
		ptsSource.push_back(Eigen::Vector3d(TransformedPts[i][0],TransformedPts[i][1],TransformedPts[i][2]));
		indexMsur.push_back(mIndex[i]-1);
		indexCad.push_back(cIndex[i]-1);
	}
}

//nMaxDelete: ����� �ִ� ����Ʈ ����
//bDeleteUniform: ������ ���� ����Ʈ�� �ʹ� �����̿� �ִ� ����Ʈ�� ������ ����
//���� �������� �ʾ���
Eigen::Matrix4d KyTrans::MovePtsNtoMRemoveLargeError(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget, double vRMSTol, int nMaxDelete, bool bDeleteUniform)
{
	Eigen::Matrix4d tMat = MovePtsNtoM(ptsSource, ptsTarget, 50, false);
	V3D box;
	
	KyMath::GetBoundingBox(ptsTarget, box);
	double vUniformTol = box.MaxWidth() / 10;

	std::vector<Eigen::Vector3d> ptsDeleted;
	do
	{
		int nPrevSize = (int)ptsDeleted.size();
		DeleteMaxRMSPoint(ptsDeleted, vUniformTol, vRMSTol, bDeleteUniform);

		//���� ����Ʈ�� ���� ��� ����
		if(nPrevSize == (int)ptsDeleted.size())
			break;

		//�ִ� ����Ʈ ���� ��ŭ ���� ��� �۾� ����
		if ((int)ptsDeleted.size() == nMaxDelete)
			break;
	} while (1);

	if (ptsDeleted.empty())
		return tMat;

	//���� ����Ʈ�� �̿��ؼ� �ٽ� ����
	tMat = MovePtsNtoN(m_ptsNewSource, m_ptsNewTarget, false);

	return tMat;
}

void KyTrans::DeleteMaxRMSPoint(std::vector<Eigen::Vector3d>& ptsDeleted, double vUniformTol, double vRMSTol, bool bDeleteUniform)
{
	int posMax = -1;
	double vMaxErr = -1e10;
	for (int i = 0; i < (int)m_errDx.size(); i++)
	{
		double vRMS = sqrt(m_errDx[i] * m_errDx[i] + m_errDy[i] * m_errDy[i] + m_errDz[i] * m_errDz[i]);
		if (vRMS > vRMSTol && vRMS > vMaxErr)
		{
			vMaxErr = vRMS;
			posMax = i;
		}
	}

	if (posMax == -1)
		return;
	
	bool bOK = true;
	if (bDeleteUniform)
	{
		for (int j = 0; j < (int)ptsDeleted.size(); j++)
		{
			double dist = KyMath::Dist(m_ptsNewTarget[posMax], ptsDeleted[j]);
			if (dist < vUniformTol)
			{
				bOK = false;
				break;
			}
		}
	}
	if (!bOK)
		return;

	ptsDeleted.push_back(m_ptsNewTarget[posMax]);

	m_errDx.erase(m_errDx.begin() + posMax);
	m_errDy.erase(m_errDy.begin() + posMax);
	m_errDz.erase(m_errDz.begin() + posMax);
	m_ptsNewSource.erase(m_ptsNewSource.begin() + posMax);
	m_ptsNewTarget.erase(m_ptsNewTarget.begin() + posMax);

	return;
}

//nMinPts: ������ �����ߴٰ� �Ǵ��ϴ� �ּ� ����Ʈ ����
//vAllowedRMS: �־��� RMS���� ���� ��� �����ߴٰ� �Ǵ�, �ٸ� �ʱ� ��ġ�� ���ؼ� ���̻� �������� ����
//vAllowedDist: �� ����Ʈ ���� �Ÿ��� �־��� �Ÿ� ������ ���� ��쿡�� ���տ� ���
Eigen::Matrix4d KyTrans::MovePtsNtoMUsingMultiInitPosition(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget, 
														   double& vResultRMS, int& nResultMatchedPts, 
														   int nMinPts/* =5 */, double vAllowedRMS/* =5 */, 
														   double vAllowedDist /* = 20 */, bool bSrcMove /* = true */)
{
	std::vector<Eigen::Vector3d> ptsRefOrg;
	GetReferenceCenter(ptsSource, ptsRefOrg);
	Eigen::Vector3d ptTagetCenter = KyMath::GetCenter(ptsTarget);

	double vMinRMS = 1e10;
	Eigen::Matrix4d tMinMat;
	double vRMS;
	nResultMatchedPts = 0;
	vResultRMS = 1e10;

	for (int i = 0; i < (int)ptsRefOrg.size(); i++)
	{
		double px, py, pz;//���� ����
		px = py = pz = 0.0;
		double stx, sty, stz;//���� ����
		stx = sty = stz = 45;

		Eigen::Matrix4d tMat;

		while (px < 360.0) {
			py = 0.0;
			while (py < 360.0) {
				pz = 0.0;
				while (pz < 360.0) {
					if (i == 6 && px == 90 && py == 225 && pz == 135)
						int a = 10;
					//�߽��̵�
					tMat.setIdentity();
					//tMat.Translate(ptTagetCenter - ptsRefOrg[i]);
					//tMat = tMat * (ptTagetCenter - ptsRefOrg[i]);
					Eigen::Vector3d temp = ptTagetCenter - ptsRefOrg[i];
					tMat.block<3, 3>(0, 0) <<temp.x(),0,0,0,temp.y(),0,0,0,temp.z();
					std::vector<Eigen::Vector3d> pts = ptsSource;
					KyTrans::MovePts(tMat, pts);

					int nMatchedPts;
					Eigen::Matrix4d mat = MovePtsNtoMUsingMultiInitPosition(px, py, pz, pts, ptsTarget, vAllowedDist, nMinPts, vRMS, nMatchedPts);

					if (vRMS < vMinRMS && nMatchedPts >= nMinPts)
					{
						tMinMat = mat * tMat;
						vMinRMS = vRMS;
						vResultRMS = vMinRMS;
						nResultMatchedPts = nMatchedPts;

						if (vRMS < vAllowedRMS)
							return tMinMat;
					}
					pz += stz;
				}
				py += sty;
			}
			px += stx;
		}
	}

	return tMinMat;
}

Eigen::Matrix4d KyTrans::MovePtsNtoMUsingMultiInitPosition(double px, double py, double pz, 
														   std::vector<Eigen::Vector3d>& ptsSource, 
														   std::vector<Eigen::Vector3d>& ptsInitTarget, 
														   double vAllowedDist, int nMinPt, double& vRMS, int& nMatchedPts)
{
	Eigen::Vector3d ptCenter = KyMath::GetCenter(ptsSource);

	//1.�ʱ� ������ŭ ȸ���ϱ�
	Eigen::Affine3d tRotX, tRotY, tRotZ;
	Eigen::Transform<double, 4, Eigen::Affine> trans;
	Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(px), Eigen::Vector3d(1, 0, 0)));
	Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(py), Eigen::Vector3d(0, 1, 0)));
	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(pz), Eigen::Vector3d(0, 0, 1)));
	Eigen::Affine3d r = rx * ry * rz;
	Eigen::Matrix4d tRot = r.matrix();
	//tRotX.RotateX(ptCenter, px / 180 * KY_PI);
	//tRotY.RotateY(ptCenter, py / 180 * KY_PI);
	//tRotZ.RotateZ(ptCenter, pz / 180 * KY_PI);
	//Eigen::Matrix4d tRot = tRotX * tRotY * tRotZ;

	KyTrans::MovePts(tRot, ptsSource);
	Eigen::Matrix4d tFinalMat = tRot;
	std::vector<Eigen::Vector3d> ptsLastSource, ptsLastTarget;

	double vRMS1, vRMS0;

	vRMS0 = vRMS = 1e10;
	int nCnt = 0;
	double vAllowedDistTol;
	double vAllowedRMSTol = 0.1;
	do
	{
		if (nCnt < 3)
			vAllowedDistTol = 500;
		else if (nCnt < 5)
			vAllowedDistTol = 100;
		else if (nCnt < 7)
			vAllowedDistTol = vAllowedDist * 3;
		else
			vAllowedDistTol = vAllowedDist;

		//�־��� ������ ������ �ּ� ��ġ�� ������ ã��
		std::vector<Eigen::Vector3d> ptsNewTarget, ptsNewSource;
		FindClosestPts(ptsSource, ptsInitTarget, ptsNewSource, ptsNewTarget, vAllowedDistTol);

		if((int)ptsNewSource.size() < nMinPt)
			break;

		KyTrans trans;
		Eigen::Matrix4d tMat = trans.MovePtsNtoN(ptsNewSource, ptsNewTarget);
		vRMS1 = trans.GetRmsError();
		vRMS = vRMS1;
		nMatchedPts = (int)ptsNewSource.size();

		if (vRMS1 < vAllowedRMSTol || (nCnt > 7 && fabs(vRMS1 - vRMS0) < vRMS1 / 1000))
		{
			tFinalMat = tMat * tFinalMat;
			ptsLastTarget = ptsNewTarget;
			ptsLastSource = ptsNewSource;

			//���� ����Ʈ�� ������ �� �ִ� ���� ���з� ����
			KyMath::DeleteSamePts(ptsLastTarget);
			if ((int)ptsLastTarget.size() < nMinPt)
			{
				vRMS = 1e10;
				break;
			}
			else
				int a = 10;

			break;
		}
		else if (vRMS1 > vRMS0*1.1)
			break;

		ptsLastTarget = ptsNewTarget;
		ptsLastSource = ptsNewSource;

		MovePts(tMat, ptsSource);//����Ʈ�� ���հ���� �̿��ؼ� �̵����Ѽ� ���� �ܰ� ����
		tFinalMat = tMat * tFinalMat;
		nCnt++;

		vRMS0 = vRMS1;
	} while (nCnt < 20);
	
	m_ptsNewTarget.clear();
	m_ptsNewSource.clear();
	for (size_t i = 0; i < ptsLastSource.size(); i++)
	{
		m_ptsNewTarget.push_back(ptsLastTarget[i]);
		m_ptsNewSource.push_back(ptsLastSource[i]);
	}
	m_errDx.clear();
	m_errDy.clear();
	m_errDz.clear();
	
	for (size_t i = 0; i < m_ptsNewSource.size(); i++)
	{
		m_errDx.push_back(m_ptsNewTarget[i].x() - m_ptsNewSource[i].x());
		m_errDy.push_back(m_ptsNewTarget[i].y() - m_ptsNewSource[i].y());
		m_errDz.push_back(m_ptsNewTarget[i].z() - m_ptsNewSource[i].z());
	}

	return tFinalMat;
}

void KyTrans::FindClosestPts(std::vector<Eigen::Vector3d>& ptsSource0, std::vector<Eigen::Vector3d>& ptsTarget0, std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget, double vAllowedDist)
{
	for (int i = 0; i < (int)ptsSource0.size(); i++)
	{
		double vMin = 1e10;
		Eigen::Vector3d ptTarget;
		for (int j = 0; j < (int)ptsTarget0.size(); j++)
		{
			double dis = KyMath::Dist(ptsSource0[i], ptsTarget0[j]);
			if (dis < vAllowedDist && dis < vMin)
			{
				vMin = dis;
				ptTarget = ptsTarget0[j];
			}
		}

		if (vMin < vAllowedDist)
		{
			ptsSource.push_back(ptsSource0[i]);
			ptsTarget.push_back(ptTarget);
		}
	}
}

void KyTrans::GetReferenceCenter(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsRefOrg)
{
	V3D box;
	for (int i = 0; i < ptsSource.size(); i++) {
		box.Update(ptsSource[i]);
	}
	ptsRefOrg.push_back(box.Center());

	double vMove = box.MaxWidth() / 2;
	Eigen::Vector3d pt;
	pt = box.Min;
	Eigen::Vector3d ptCenter1 = box.Center() + (pt - box.Center()).normalized() * vMove;
	pt.x() = box.Max.x();
	Eigen::Vector3d ptCenter2 = box.Center() + (pt - box.Center()).normalized() * vMove;
	pt.y() = box.Max.y();
	Eigen::Vector3d ptCenter3 = box.Center() + (pt - box.Center()).normalized() * vMove;
	pt.x() = box.Min.x();
	Eigen::Vector3d ptCenter4 = box.Center() + (pt - box.Center()).normalized() * vMove;

	pt = box.Min; 
	pt.z() = box.Max.z();
	Eigen::Vector3d ptCenter5 = box.Center() + (pt - box.Center()).normalized() * vMove * 0.75;
	pt.x() = box.Max.x();
	Eigen::Vector3d ptCenter6 = box.Center() + (pt - box.Center()).normalized() * vMove * 0.75;
	pt.y() = box.Max.y();
	Eigen::Vector3d ptCenter7 = box.Center() + (pt - box.Center()).normalized() * vMove * 0.75;
	pt.x() = box.Min.x();
	Eigen::Vector3d ptCenter8 = box.Center() + (pt - box.Center()).normalized() * vMove * 0.75;

	ptsRefOrg.push_back(ptCenter1);
	ptsRefOrg.push_back(ptCenter2);
	ptsRefOrg.push_back(ptCenter3);
	ptsRefOrg.push_back(ptCenter4);
	ptsRefOrg.push_back(ptCenter5);
	ptsRefOrg.push_back(ptCenter6);
	ptsRefOrg.push_back(ptCenter7);
	ptsRefOrg.push_back(ptCenter8);
}

Eigen::Matrix4d KyTrans::MovePtsNtoM(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget,double vMatchingError,bool bSrcMove, bool bOnlyRotateZ)
{
	//m_tMat.MakeIdentity();
	m_tMat.setIdentity();
	if(ptsSource.size() < 1 || ptsTarget.size() < 1)
		return m_tMat;

	std::vector<int> indexMsur;
	std::vector<int> indexCad;

	//01. N to M ����
	std::vector<Eigen::Vector3d> ptsInitSource = ptsSource;
	MovePtsNtoM(ptsSource, ptsTarget,indexMsur,indexCad,vMatchingError, bOnlyRotateZ);

	//�������� ����
	m_NtoMIdxMsr = indexMsur;
	m_NtoMIdxCad = indexCad;

	//02. ���յ� �����͸��� ������� ����
	m_ptsNewSource.clear();
	m_ptsNewTarget.clear();

	int nCnt = 0;
	int theJ;
	for(size_t i=0; i<ptsInitSource.size(); i++)
	{
		bool bWork = false;
		for(size_t j=0; j<indexMsur.size(); j++)
		{
			int n = indexMsur[j];
			if(nCnt == indexMsur[j])
			{
				theJ = (int)j;
				bWork = true;
				break;
			}
		}

		if(bWork)
		{
			Eigen::Vector3d p1 = ptsInitSource[i];
			int n = indexCad[theJ];
			Eigen::Vector3d p2 = ptsTarget[n];
			m_ptsNewSource.push_back(ptsInitSource[i]);
			m_ptsNewTarget.push_back(ptsTarget[n]);
		}
		nCnt++;
	}

	if(m_ptsNewSource.size() == 0)
		return m_tMat;

	//03. ������� ���ĵ� �����͸� �̿��ؼ� ��ȯ��� ���
	m_tMat = MovePtsNtoN(m_ptsNewSource, m_ptsNewTarget,true);

	m_errDx.clear();
	m_errDy.clear();
	m_errDz.clear();
	for(size_t i=0; i< m_ptsNewSource.size(); i++)
	{
		m_errDx.push_back(m_ptsNewTarget[i].x() - m_ptsNewSource[i].x());
		m_errDy.push_back(m_ptsNewTarget[i].y() - m_ptsNewSource[i].y());
		m_errDz.push_back(m_ptsNewTarget[i].z() - m_ptsNewSource[i].z());
	}

	return m_tMat;
}


Eigen::Matrix4d KyTrans::MovePtsToXYPlane(std::vector<Eigen::Vector3d>& source,bool bMove)
{
	Eigen::Vector3d ptrOrg,vPlaneNor;
	KyPrimFit KyPrimFitTemp;
	bMove = KyPrimFitTemp.MakePlane(source, ptrOrg, vPlaneNor);

	m_tMat = GetMatrixPtsToXYPlane(vPlaneNor,ptrOrg);

	if(bMove)
		UpdatePoint(source,m_tMat);

	return m_tMat;
}


Eigen::Matrix4d KyTrans::GetMatrixPtsToXYPlane(Eigen::Vector3d vPlaneNor, Eigen::Vector3d ptOrg)
{
	Eigen::Affine3d tRotY, tRotZ;

	Eigen::Matrix4d transM,rotY,rotX;
	//1. Y-axis ��ġ
	double yAng = -atan2(vPlaneNor.x(),vPlaneNor.z());
	if(yAng > M_PI_2)yAng -= M_PI;
	if(yAng < -M_PI_2)yAng += M_PI;
	//rotY.RotateY(yAng);
	//vPlaneNor = rotY*vPlaneNor;
	Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(yAng), Eigen::Vector3d(0, 1, 0)));

	//2. Z-axis ��ġ
	double xAng = atan2(vPlaneNor.y(),vPlaneNor.z());
	if(xAng > M_PI_2)xAng -= M_PI;
	if(xAng < -M_PI_2)xAng += M_PI;
	//vPlaneNor = rotX*vPlaneNor;
	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(xAng), Eigen::Vector3d(0, 0, 1)));
	Eigen::Affine3d r = ry * rz;
	Eigen::Matrix4d tRot = r.matrix();
	vPlaneNor = tRot.block<3,3>(0,0) * vPlaneNor;
	//3. ���� �̵�
	//transM.SetT(-ptOrg);
	//rotX.RotateX(xAng);
	Eigen::Affine3d t(Eigen::Translation3d(-ptOrg));
	Eigen::Matrix4d m = (t * r).matrix();
	return m;
	//return rotX*rotY*transM;
}

int KyTrans::GetAverageError(double& dErrX,double& dErrY,double& dErrZ,bool bRMS)
{
	dErrX = dErrY = dErrZ = 0;
	
	for(size_t i=0; i<m_errDx.size(); i++)
	{
		if(bRMS)
			dErrX += fabs(m_errDx[i]);
		else
			dErrX += (m_errDx[i]);
	}
	if(m_errDx.size() > 0)
		dErrX /= m_errDx.size();

	for(size_t i=0; i<m_errDy.size(); i++)
	{
		if(bRMS)
			dErrY += fabs(m_errDy[i]);
		else
			dErrY += (m_errDy[i]);
	}
	if(m_errDy.size() > 0)
		dErrY /= m_errDy.size();

	for(size_t i=0; i<m_errDz.size(); i++)
	{
		if(bRMS)
			dErrZ += fabs(m_errDz[i]);
		else
			dErrZ += (m_errDz[i]);
	}
	if(m_errDz.size() > 0)
		dErrZ /= m_errDz.size();

	return (int)m_errDx.size();
}

Eigen::Matrix4d KyTrans::MoveToZAxis(Eigen::Vector3d vPlaneNor, Eigen::Vector3d ptOrg)
{
	Eigen::Matrix4d transM,rotY,rotX;
		
	//1. Y-axis ��ġ
	double yAng = -atan2(vPlaneNor.x(),vPlaneNor.z());
	if(yAng > M_PI_2)yAng -= M_PI;
	if(yAng < -M_PI_2)yAng += M_PI;
	//rotY.RotateY(yAng);
	//vPlaneNor = rotY*vPlaneNor;
	Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(yAng), Eigen::Vector3d(0, 1, 0)));

	//2. Z-axis ��ġ
	double xAng = atan2(vPlaneNor.y(),vPlaneNor.z());
	if(xAng > M_PI_2)xAng -= M_PI;
	if(xAng < -M_PI_2)xAng += M_PI;
	//vPlaneNor = rotX*vPlaneNor;
	Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(xAng), Eigen::Vector3d(1, 0, 0)));
	Eigen::Affine3d r = ry * rx;
	Eigen::Matrix4d tRot = r.matrix();
	vPlaneNor = tRot.block<3, 3>(0, 0) * vPlaneNor;

	//3. ���� �̵�
	//transM.SetT(-ptOrg);
	//rotX.RotateX(xAng);
	Eigen::Affine3d t(Eigen::Translation3d(-ptOrg));
	Eigen::Matrix4d m = (t * r).matrix();
	return m;

	//return rotX*rotY*transM;
}


Eigen::Matrix4d KyTrans::MoveToXZAxis(Eigen::Vector3d ptOrg, Eigen::Vector3d vX, Eigen::Vector3d vZ)
{
	Eigen::Matrix4d tMatZAxis = MoveToZAxis(vZ, ptOrg);

	Eigen::Vector3d ptZ = tMatZAxis.block<3, 3>(0, 0) * (ptOrg + vZ);
	if (ptZ.z() < 0)
	{
		//Eigen::Matrix4d rotX;
		//rotX.RotateX(KY_PI);
		//tMatZAxis = rotX * tMatZAxis;	
		Eigen::Affine3d rx =
			Eigen::Affine3d(Eigen::AngleAxisd(M_PI, Eigen::Vector3d(1, 0, 0)));
		tMatZAxis = rx.matrix() * tMatZAxis;
	}
	Eigen::Vector3d ptX = tMatZAxis.block<3, 3>(0, 0) * (ptOrg + vX);
	double vAng = -atan2(ptX.y(), ptX.x());
	//Eigen::Matrix4d rotZ;
	//rotZ.RotateZ(vAng);

	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(vAng), Eigen::Vector3d(0, 0, 1)));

	//return rotZ*tMatZAxis;
	return rz.matrix()*tMatZAxis;
}

bool KyTrans::IsValidMatrix(Eigen::Matrix4d tMat)
{
	double v[16];
	//tMat.Get4x4(v);
	for (int i = 0; i < 16; i++)
	{
		v[i] = tMat(i);
	}
	double sum = 0;
	for(int i=0; i<12; i++)
	{
		if(fabs(v[i]) > 1.0 + 1.0E-6)
			return false;
		if(v[i] != v[i])
			return false;
	}

	if(fabs(v[3]) > 1.0E-6)
		return false;
	if(fabs(v[7]) > 1.0E-6)
		return false;
	if(fabs(v[11]) > 1.0E-6)
		return false;

	if(fabs(v[15])-1 > 1.0E-6)
		return false;

	for(int i=0; i<16; i++)
	{
		if(fabs(v[i]) > 1e15)
			return false;
	}

	return true;
}

//bOrg=true : ù��° ������ ���� ����
//bOrg=false: ù��° ������ ���� �ι�° ������ ���� �߽�
Eigen::Matrix4d KyTrans::MovePtsBy5PtsCord(std::vector<Eigen::Vector3d> ptsCord,bool bOrg)
{
	Eigen::Vector3d ptOrg, ptX, ptXY;
	Eigen::Vector3d pt1,pt2;

	// ����, x�� ����
	pt1 = ptsCord[0];
	pt2 = ptsCord[1];

	// ��� ������ ���� ����Ʈ �߰�
	std::vector<Eigen::Vector3d> ptsPlane;
	for(size_t i=2; i<ptsCord.size(); i++)
		ptsPlane.push_back(ptsCord[i]);

	// z�� ����
	Eigen::Vector3d ptPlaneOrg,vZDir;
	KyPrimFit kypf;
	kypf.MakePlane(ptsPlane, ptPlaneOrg, vZDir);
	//KyPrimFit::MakePlane(ptsPlane,ptPlaneOrg,vZDir);
	if(vZDir.z() < 0.) vZDir *= -1;

	// ����Ʈ ��� ����
	Eigen::Vector3d ptP1,ptP2;
	DistPtToPlane(ptPlaneOrg,vZDir,pt1,ptP1);
	DistPtToPlane(ptPlaneOrg,vZDir,pt2,ptP2);

	if(bOrg)
		ptOrg = ptP1;
	else
		ptOrg = (ptP2+ptP1)/2.;
	ptX = ptOrg + (ptP2-ptP1);

	//ptXY = ptOrg + vZDir.Cross(ptX);
	ptXY = ptOrg + vZDir.cross(ptX);
	Eigen::Matrix4d tMat = MovePtsBy3PtsCord(ptOrg,ptX,ptXY);

	//
	Eigen::Vector3d ptNewX = tMat.block<3, 3>(0, 0) *pt2;
	Eigen::Matrix4d tRotZ;
	if(ptNewX.x() < 0.)
	{
		//tRotZ.RotateZ(KY_PI);
		//tMat = tRotZ*tMat;
		Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(M_PI), Eigen::Vector3d(0, 0, 1)));
		tMat = rz.matrix() * tMat;
	}


	return tMat;
}


Eigen::Matrix4d KyTrans::MovePtsBy5PtsCordOFD(std::vector<Eigen::Vector3d> ptsCord,bool bIsY)
{
	Eigen::Matrix4d tmat;

	if(ptsCord.size()<5)
	{
		return tmat;
	}


	std::vector<Eigen::Vector3d> plane;
	for(size_t i=2; i<ptsCord.size(); i++)
		plane.push_back(ptsCord[i]);

	Eigen::Vector3d planeorg, planenor;
	KyPrimFit pf;
	pf.MakePlane(plane,planeorg,planenor);


	tmat=MovePtsBy5PtsCord(ptsCord,true);

	if(bIsY)
	{
		//Eigen::Matrix4d zRot;
		//zRot.Rotate(ptsCord[0],planenor,-M_PI/2.);
		//tmat = zRot*tmat;
		/*Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(-M_PI/2), Eigen::Vector3d(1, 0, 0)));
		Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(-M_PI / 2), Eigen::Vector3d(0, 1, 0)));
		Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(-M_PI / 2), Eigen::Vector3d(0, 0, 1)));*/
		//Eigen::Affine3d r = rx * ry * rz;
		Eigen::AngleAxisd r(-M_PI_2, planenor);
		Eigen::Affine3d t(Eigen::Translation3d(ptsCord.at(0)));
		Eigen::Matrix4d m = (t * r).matrix();
		tmat = m * tmat;
	}
	//������ ���� ��ġ
	Eigen::Vector3d org = ptsCord[0];
	DistPtToPlane(planeorg,planenor,org,org);
	//Eigen::Matrix4d zTrans;
	//zTrans.Translate(0.,0.,org.z()-ptsCord[0].z()+OFFSET_ORG); //������ ���� ��ġ
	//tmat = zTrans*tmat;
	Eigen::Affine3d t(Eigen::Translation3d(0., 0., org.z() - ptsCord[0].z() + OFFSET_ORG));
	Eigen::Matrix4d m = t.matrix();
	tmat = m * tmat;
	return tmat;
}



Eigen::Matrix4d KyTrans::MovePtsBy5PtsCordXZ(std::vector<Eigen::Vector3d> ptsCord,bool bOrg)
{
	Eigen::Matrix4d tMat = MovePtsBy5PtsCord(ptsCord,bOrg);

// 	KyP3D ptOrg;
// 	if(bOrg)
// 		ptOrg = tMat*ptsCord[ ];
// 	else
// 		ptOrg = tMat*(ptsCord[ ]+ptsCord[1])/2.;
	//Eigen::Matrix4d xRot;
	//xRot.RotateX(KY_PI/2.);
	//tMat = xRot * tMat;
	Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(M_PI/2), Eigen::Vector3d(0, 0, 1)));
	tMat = rx.matrix() * tMat;

	return tMat;
}

Eigen::Matrix4d KyTrans::GetMatrix3PtsCoord(Eigen::Vector3d ptOrg, Eigen::Vector3d ptX, Eigen::Vector3d ptXY)
{
	Eigen::Vector3d vX, vY, vZ;
	//vX = (ptX - ptOrg).UnitV();
	//vY = (ptXY - ptOrg).UnitV();
	//vZ = vX.Cross(vY).UnitV();
	//vY = vZ.Cross(vX);
	vX = (ptX - ptOrg).normalized();
	vY = (ptXY - ptOrg).normalized();
	vZ = (vX.cross(vY)).normalized();
	vY = vZ.cross(vX);

	Eigen::Matrix4d tMat;
	//tMat.SetX(vX);
	//tMat.SetY(vY);
	//tMat.SetZ(vZ);
	tMat.setIdentity();
	tMat.block<3, 1>(0, 0) = vX;
	tMat.block<3, 1>(0, 1) = vY;
	tMat.block<3, 1>(0, 2) = vZ;

	Eigen::Vector3d v0 = ptOrg;
	//tMat.SetT(v0);
	tMat.block<3, 1>(0, 3) = v0;

	//return tMat.InverseM();
	return tMat.inverse();
}

Eigen::Matrix4d KyTrans::MovePtsBy3PtsCord(Eigen::Vector3d ptOrg, Eigen::Vector3d ptX, Eigen::Vector3d ptXY)
{

	Eigen::Matrix4d tMat,tMatXYPlane,tMatOrg,tMatXAxis;

	//1. xy ������� �̵���Ű�� ���
	std::vector<Eigen::Vector3d> arrPts;
	arrPts.push_back(ptOrg);
	arrPts.push_back(ptX);
	arrPts.push_back(ptXY);

	tMatXYPlane = MovePtsToXYPlane(arrPts);

	//2. ���� �̵�
	//tMatOrg.SetT(-arrPts[0]);
	Eigen::Affine3d t(Eigen::Translation3d(-arrPts[0]));
	tMatOrg = t.matrix();

	MovePts(tMatOrg,arrPts);
	
	//3. X�� ��ġ
	double ang = atan2(arrPts[1].y(),arrPts[1].x());
	//tMatXAxis.RotateZ(-ang);
	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(-ang), Eigen::Vector3d(0, 0, 1)));
	tMat = rz.matrix() * tMat;

	tMat = tMatXAxis*tMatOrg*tMatXYPlane;
	
	//4. Y���� -�� ��� +�� ����
	Eigen::Vector3d ptY = tMat.block<3, 3>(0, 0) *ptXY;
	Eigen::Matrix4d tRotX;
	if(ptY.y() < 0.)
	{
 		//tRotX.RotateX(KY_PI);
 		//tMat = tRotX*tMat;
		Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(M_PI), Eigen::Vector3d(0, 0, 1)));
		tMat = rx.matrix() * tMat;
	}
	
	return tMat;
}

Eigen::Matrix4d KyTrans::MovePtsBy3PtsCordPipe(Eigen::Vector3d ptOrg, Eigen::Vector3d ptX, Eigen::Vector3d ptXY)
{
	Eigen::Matrix4d tMat,tMatXYPlane,tMatOrg,tMatXAxis;
	//1. xy ������� �̵���Ű�� ���
	std::vector<Eigen::Vector3d> arrPts;
	arrPts.push_back(ptOrg);
	arrPts.push_back(ptX);
	arrPts.push_back(ptXY);

	tMatXYPlane = MovePtsToXYPlane(arrPts);

	//2. ���� �̵�
	//tMatOrg.SetT(-arrPts[0]);
	Eigen::Affine3d t(Eigen::Translation3d(-arrPts[0]));
	tMatOrg = t.matrix();

	MovePts(tMatOrg,arrPts);

	//3. X�� ��ġ
	double ang = atan2(arrPts[1].y(),arrPts[1].x());
	//tMatXAxis.RotateZ(-ang);
	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(-ang), Eigen::Vector3d(0, 0, 1)));
	tMat = rz.matrix() * tMat;

	tMat = tMatXAxis*tMatOrg*tMatXYPlane;

	//4. Y���� +�� ��� -�� ����
	Eigen::Vector3d ptY = tMat.block<3, 3>(0, 0) *ptXY;

	Eigen::Matrix4d tRotX;
	if(ptY.y() > 0.)
	{
		//tRotX.RotateX(KY_PI);
		//tMat = tRotX*tMat;
		Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(M_PI), Eigen::Vector3d(0, 0, 1)));
		tMat = rx.matrix() * tMat;
	}

	return tMat;
}

//ptSt2->ptEd2 ���� ptSt1->ptEd1�࿡ ��ġ
Eigen::Matrix4d KyTrans::MovePtsBy4Pts(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2)
{
	Eigen::Matrix4d tMatOrg,tMatZ,tMatY,tMatX,tMat;

	//���� ��ġ
	//tMatOrg.SetT(-(ptSt2-ptSt1));
	Eigen::Affine3d t(Eigen::Translation3d(-(ptSt2 - ptSt1)));
	tMatOrg = t.matrix();

	ptSt2 = tMatOrg.block<3, 3>(0, 0) *ptSt2;
	ptEd2 = tMatOrg.block<3, 3>(0, 0) *ptEd2;

	//z�� ȸ��
	double ang1,ang2,ang;
	ang1 = atan2(ptSt1.y()-ptEd1.y(),ptSt1.x() -ptEd1.x());
	ang2 = atan2(ptSt2.y() -ptEd2.y(),ptSt2.x() -ptEd2.x());
	ang = ang1 - ang2;
	//tMatZ.RotateZ(ptSt1,ang);
	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(ang), ptSt1));
	tMatZ = rz.matrix();

	ptSt2 = tMatZ.block<3, 3>(0, 0) *ptSt2;
	ptEd2 = tMatZ.block<3, 3>(0, 0) *ptEd2;

	//y�� ȸ��
	ang1 = atan2(ptSt1.x() -ptEd1.x(),ptSt1.z() -ptEd1.z());
	ang2 = atan2(ptSt2.x() -ptEd2.x(),ptSt2.z() -ptEd2.z());
	ang = ang1 - ang2;
	//tMatY.RotateY(ptSt1,ang);
	Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(ang), ptSt1));
	tMatY = ry.matrix();

	ptSt2 = tMatY.block<3, 3>(0, 0) *ptSt2; 
	ptEd2 = tMatY.block<3, 3>(0, 0) *ptEd2;

	//x�� ȸ��
	ang1 = atan2(ptSt1.z() -ptEd1.z(),ptSt1.y() -ptEd1.y());
	ang2 = atan2(ptSt2.z() -ptEd2.z(),ptSt2.y() -ptEd2.y());
	ang = ang1 - ang2;
	//tMatX.RotateX(ptSt1,ang);
	Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(ang), ptSt1));
	tMatY = rx.matrix();


	ptSt2 = tMatX.block<3, 3>(0, 0) *ptSt2;
	ptEd2 = tMatX.block<3, 3>(0, 0) *ptEd2;

	tMat = tMatX*tMatY*tMatZ*tMatOrg;

	return tMat;
}

//ptSt2->ptEd2 ���� ptSt1->ptEd1�࿡ ��ġ : Z�࿡ ���� ȸ���� ����
Eigen::Matrix4d KyTrans::MovePtsBy4PtsWithZAxis(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2)
{
	Eigen::Matrix4d tMatOrg,tMatZ,tMat;

	//���� ��ġ
	//tMatOrg.SetT(-(ptSt2-ptSt1));
	Eigen::Affine3d t(Eigen::Translation3d(-(ptSt2 - ptSt1)));
	tMatOrg = t.matrix();


	//z�� ȸ��
	double ang1,ang2,ang;
	ang1 = atan2(ptSt1.y()-ptEd1.y(),ptSt1.x() -ptEd1.x());
	ang2 = atan2(ptSt2.y() -ptEd2.y(),ptSt2.x() -ptEd2.x());
	ang = ang1 - ang2;
	//tMatZ.RotateZ(ptSt1,ang);
	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(ang), ptSt1));
	tMatZ = rz.matrix();

	tMat = tMatZ*tMatOrg;

	return tMat;
}

//ptSt2->ptEd2 ���� ptSt1->ptEd1�࿡ ��ġ : Y�࿡ ���� ȸ���� ����
Eigen::Matrix4d KyTrans::MovePtsBy4PtsWithYAxis(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2)
{
	Eigen::Matrix4d tMatOrg,tMatZ,tMat;

	//���� ��ġ	
	//tMatOrg.SetT(-(ptSt2-ptSt1));
	Eigen::Affine3d t(Eigen::Translation3d(-(ptSt2 - ptSt1)));
	tMatOrg = t.matrix();


	//y�� ȸ��
	double ang1,ang2,ang;
	ang1 = atan2(ptSt1.x()-ptEd1.x(),ptSt1.z() -ptEd1.z());
	ang2 = atan2(ptSt2.x() -ptEd2.x(),ptSt2.z() -ptEd2.z());
	ang = ang1 - ang2;
	//tMatZ.RotateY(ptSt1,ang);
	Eigen::Affine3d ry = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(ang), ptSt1));
	tMatZ = ry.matrix();

	tMat = tMatZ*tMatOrg;

	return tMat;
}

//ptSt2->ptEd2 ���� ptSt1->ptEd1�࿡ ��ġ : X�࿡ ���� ȸ���� ����
Eigen::Matrix4d KyTrans::MovePtsBy4PtsWithXAxis(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2)
{
	Eigen::Matrix4d tMatOrg,tMatZ,tMat;

	//���� ��ġ
	//tMatOrg.SetT(-(ptSt2-ptSt1));
	Eigen::Affine3d t(Eigen::Translation3d(-(ptSt2 - ptSt1)));
	tMatOrg = t.matrix();


	//y�� ȸ��
	double ang1,ang2,ang;
	ang1 = atan2(ptSt1.z()-ptEd1.z(),ptSt1.y() -ptEd1.y());
	ang2 = atan2(ptSt2.z() -ptEd2.z(),ptSt2.y() -ptEd2.y());
	ang = ang1 - ang2;
	//tMatZ.RotateX(ptSt1,ang);
	Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(ang), ptSt1));
	tMatZ = rz.matrix();

	tMat = tMatZ*tMatOrg;

	return tMat;
}

void KyTrans::MovePts(Eigen::Matrix4d tMat, std::vector<Eigen::Vector3d>& arrPts)
{
	for(size_t i=0; i<arrPts.size(); i++)
	{		
		//arrPts[i] = tMat*arrPts[i];
		arrPts[i] = (tMat.block<3, 3>(0, 0)) *arrPts[i];
	}
}

void KyTrans::MovePts(Eigen::Matrix4d tMat, std::vector<std::vector<Eigen::Vector3d>>& arrPtss)
{
	for (size_t i = 0; i < arrPtss.size(); i++)
	{
		for(size_t j=0; j<arrPtss[i].size(); j++)
			//arrPtss[i][j] = tMat * arrPtss[i][j];
			arrPtss[i][j] = (tMat.block<3, 3>(0, 0)) * arrPtss[i][j];
	}
}


void KyTrans::MovePts(Eigen::Matrix4f tMat, std::vector<Eigen::Vector3f>& arrPts)
{
	for (size_t i = 0; i < arrPts.size(); i++)
	{
		//arrPts[i] = tMat*arrPts[i];
		arrPts[i] = (tMat.block<3, 3>(0, 0)) *arrPts[i];
	}
}

//Pt1, ptNor1 -> X��
//Pt2 + ptNor2 * vLen2 = pt3 -> xy��� ���� ��
Eigen::Matrix4d KyTrans::GetMatrixToPipeForH(Eigen::Vector3d Pt1, Eigen::Vector3d ptNor1, Eigen::Vector3d Pt2, Eigen::Vector3d ptNor2,double vLen1,double vLen2)
{
	Eigen::Matrix4d tMat;
	//tMat.MakeIdentity();
	tMat.setIdentity();

	Eigen::Vector3d Pt1_2 = Pt1 + ptNor1 * vLen1;
	Eigen::Vector3d Pt2_2 = Pt2 + ptNor2 * vLen2;

	Eigen::Vector3d ptOrg,ptX,ptZ;
	ptOrg = Pt1;
	ptX = Pt1_2;
	ptZ = Pt2_2;

	tMat = MovePtsBy3PtsCord(ptOrg,ptX,ptZ);
	if(!IsValidMatrix(tMat) || tMat.isIdentity())
	{
		//tMat.MakeIdentity();
		tMat.setIdentity();

		return tMat;
	}

	ptX = tMat.block<3, 3>(0, 0) * ptX;
	ptZ = tMat.block<3, 3>(0, 0) * ptZ;
	if(ptZ.y() < 0.0)
	{
		//Eigen::Matrix4d tRot;
		//tRot.RotateX(KY_PI);
		//tMat = tRot * tMat;
		Eigen::Affine3d rz = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(M_PI), Eigen::Vector3d(0, 0, 1)));
		tMat = rz.matrix() * tMat;

	}

	return tMat;
}

//Pt1, ptNor1 -> X��
//Pt2, ptNor2 -> XY ���� ����
Eigen::Matrix4d KyTrans::GetMatrixToPipeForB(Eigen::Vector3d Pt1, Eigen::Vector3d ptNor1, Eigen::Vector3d Pt2, Eigen::Vector3d ptNor2,double vLen1,double vLen2)
{
	Eigen::Vector3d Pt1_2 = Pt1 + ptNor1 * vLen1;
	Eigen::Vector3d Pt2_2 = Pt2 + ptNor2 * vLen2;

	Eigen::Vector3d iPt2 = Pt2;
	Eigen::Vector3d iPt2_2 = Pt2_2;

	Eigen::Matrix4d tMat = GetMatrixToPipeForH(Pt1,ptNor1,Pt2,ptNor2,vLen1,vLen2);

	Pt1   = tMat.block<3, 3>(0, 0) * Pt1;
	Pt1_2 = tMat.block<3, 3>(0, 0) * Pt1_2;
	Pt2   = tMat.block<3, 3>(0, 0) * Pt2;
	Pt2_2 = tMat.block<3, 3>(0, 0) * Pt2_2;

	double dAng = M_PI/180 * 0.001;
	Eigen::Matrix4d tRot;
	//tRot.RotateX(dAng);
	Eigen::Affine3d rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(dAng), Eigen::Vector3d(1, 0, 0)));
	tRot = rx.matrix();

	int nCnt = 0;
	int minCnt;
	double minV = 1e10;
	do 
	{
		if(fabs(Pt2.z() - Pt2_2.z()) < minV)
		{
			minV = fabs(Pt2.z() - Pt2_2.z());
			minCnt = nCnt;
		}

		if(dAng * nCnt > M_PI * 2.)
			break;

		nCnt++;

		Pt2   = tRot.block<3, 3>(0, 0) * Pt2;
		Pt2_2 = tRot.block<3, 3>(0, 0) * Pt2_2;
	} while (1);

	//tRot.MakeIdentity();
	//tRot.RotateX(dAng*minCnt);
	tRot.setIdentity();
	rx = Eigen::Affine3d(Eigen::AngleAxisd(KyMath::Radian(dAng * minCnt), Eigen::Vector3d(1, 0, 0)));
	tRot = rx.matrix();


	tMat  = tRot * tMat;

	iPt2   = tMat.block<3, 3>(0, 0) * iPt2;
	iPt2_2 = tMat.block<3, 3>(0, 0) * iPt2_2;
	return tMat;
}

//ptThe���� ��� ���� ������Ų ��(intP) ������ �Ÿ� ���
double KyTrans::DistPtToPlane(Eigen::Vector3d ptOrgPlane, Eigen::Vector3d vNorPlane, Eigen::Vector3d ptThe, Eigen::Vector3d&intP)
{
	double d = ptOrgPlane.x() *vNorPlane.x() + ptOrgPlane.y() *vNorPlane.y() + ptOrgPlane.z() *vNorPlane.z();

	return DistPtToPlane(vNorPlane.x(), vNorPlane.y(), vNorPlane.z(), d, vNorPlane, ptThe, intP);
}


//theP���� vec�������� ���(a,b,c,d)�� ������ �� �Ÿ��� ���
double KyTrans::DistPtToPlane(double a, double b, double c, double d, Eigen::Vector3d vec, Eigen::Vector3d theP, Eigen::Vector3d&intP)
{
	double u, v, w;//���� ����

	u = vec.x();
	v = vec.y();
	w = vec.z();

	intP.x() = intP.y() = intP.z() = 1e21;//�ʱ�ȭ
									//���� ������ ������ ���� �������� ����� �� ����.
	intP.setZero();
	//intP.SetZero();
	if (fabs(a*u + b*v + c*w) < 1e-10)
		return 0.;

	//find t
	double t = ((d - (a*theP.x() + b*theP.y() + c*theP.z())) / (a*u + b*v + c*w));

	intP.x() = u*t + theP.x();
	intP.y() = v*t + theP.y();
	intP.z() = w*t + theP.z();

	//return intP.Dist(theP);
	return  KyMath::Dist(intP, theP);
}
