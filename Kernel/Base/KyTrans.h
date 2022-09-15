#pragma once
#include "KyGeoICP.h"
#include "open3d/Open3D.h"
//#include <KyBase/KyBaseExportDll.h>
//#include <KyBase/KyMatrix.h>

enum CORD_TRANS_TYPE
{
	TRANS_TYPE_1POINT = 1001,
	TRANS_TYPE_2POINT = 1002,
	TRANS_TYPE_3POINT = 1003,
	TRANS_TYPE_4POINT = 1004,
	TRANS_TYPE_O1_X2  = 1005,
	TRANS_TYPE_X12_Y3 = 1006,
	TRANS_TYPE_X12_Z3 = 1007,
	TRANS_TYPE_ANG    = 1008,
	TRANS_TYPE_DIST   = 1009,

	TRANS_TYPE_MOVE   = 1101,
	TRANS_TYPE_ROTATE = 1102,
	TRANS_TYPE_P2     = 1103,

	TRANS_TYPE_AUTO1  = 1201,
	TRANS_TYPE_AUTO2  = 1202,

	TRANS_TYPE_NONE   = 9999
};

class KyTrans
{
public:
	KyTrans(void);
	~KyTrans(void);


public:
	KyGeoICP m_icp;
	Eigen::Matrix4d	m_tMat;
	bool m_bMappingFixX;
	bool m_bMappingFixY;
	bool m_bMappingFixZ;

	std::vector<double> m_errDx;
	std::vector<double> m_errDy;
	std::vector<double> m_errDz;

	std::vector<int> m_NtoMIdxMsr;
	std::vector<int> m_NtoMIdxCad;
	double m_vRmsError;

	std::vector<Eigen::Vector3d> m_ptsNewTarget;
	std::vector<Eigen::Vector3d> m_ptsNewSource;

public:
	bool IsValidMatrix(Eigen::Matrix4d tMat);

	void UpdatePoint(std::vector<Eigen::Vector3d>& pts, Eigen::Matrix4d& tMat);
	void UpdatePoint(std::vector<Eigen::Vector3f>& pts, Eigen::Matrix4f& tMat);


public:

	double GetRmsError();
	double GetPipeRmsError();
	int GetAverageError(double& dErrX,double& dErrY,double& dErrZ,bool bRMS=true);

	//01 Registration
	//01-1 N to N
	Eigen::Matrix4d MovePtsNtoN(std::vector<Eigen::Vector3d>& source, std::vector<Eigen::Vector3d>& target,bool bSrcMove=true);
	//01-2 N to M
	void MovePtsNtoM(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget,std::vector<int>& indexMsur,std::vector<int>& indexCad,double vMatchingError, bool bOnlyRotateZ= false);
	Eigen::Matrix4d MovePtsNtoM(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget,double vMatchingError=500.,bool bSrcMove=true,bool bOnlyRotateZ=false);
	Eigen::Matrix4d MovePtsNtoMRemoveLargeError(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget, double vRMSTol, int nMaxDelete, bool bDeleteUniform);
	void DeleteMaxRMSPoint(std::vector<Eigen::Vector3d>& ptsDeleted, double vUniformTol, double vRMSTol, bool bDeleteUniform);

	//�ʱⰪ�� �߽ɻӸ� �ƴ϶� �߰��� 8���� �� �̵����� ���鼭 ���� ��ġ�� ã��
	Eigen::Matrix4d MovePtsNtoMUsingMultiInitPosition(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget, double& vResultRMS, int& nResultMatchedPts, int nMinPts=5, double vAllowedRMS=5, double vAllowedDist = 10, bool bSrcMove = true);
	Eigen::Matrix4d MovePtsNtoMUsingMultiInitPosition(double px, double py, double pz, std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget, double vAllowedDist, int nMinPt, double& vRMS, int& nMatchedPts);
	void FindClosestPts(std::vector<Eigen::Vector3d>& ptsSource0, std::vector<Eigen::Vector3d>& ptsTarget0, std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsTarget, double vAllowedDist);
	void GetReferenceCenter(std::vector<Eigen::Vector3d>& ptsSource, std::vector<Eigen::Vector3d>& ptsRefOrg);

	//02 Move
	//02-1 Move to XY Plane
	Eigen::Matrix4d GetMatrixPtsToXYPlane(Eigen::Vector3d vPlaneNor, Eigen::Vector3d ptOrg);
	Eigen::Matrix4d MovePtsToXYPlane(std::vector<Eigen::Vector3d>& source,bool bMove=true);

	Eigen::Matrix4d MoveToZAxis(Eigen::Vector3d vPlaneNor, Eigen::Vector3d ptOrg);

	Eigen::Matrix4d MoveToXZAxis(Eigen::Vector3d ptOrg, Eigen::Vector3d vX, Eigen::Vector3d vZ);

	//
	Eigen::Matrix4d MovePtsBy5PtsCord(std::vector<Eigen::Vector3d> ptsCord,bool bOrg=true);
	Eigen::Matrix4d MovePtsBy5PtsCordXZ(std::vector<Eigen::Vector3d> ptsCord,bool bOrg=true);
	Eigen::Matrix4d MovePtsBy3PtsCord(Eigen::Vector3d ptOrg, Eigen::Vector3d ptX, Eigen::Vector3d ptXY);
	Eigen::Matrix4d MovePtsBy3PtsCordPipe(Eigen::Vector3d ptOrg, Eigen::Vector3d ptX, Eigen::Vector3d ptXY);

	Eigen::Matrix4d MovePtsBy4Pts(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2);
	Eigen::Matrix4d MovePtsBy4PtsWithXAxis(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2);
	Eigen::Matrix4d MovePtsBy4PtsWithYAxis(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2);
	Eigen::Matrix4d MovePtsBy4PtsWithZAxis(Eigen::Vector3d ptSt1, Eigen::Vector3d ptEd1, Eigen::Vector3d ptSt2, Eigen::Vector3d ptEd2);

	Eigen::Matrix4d MovePtsBy5PtsCordOFD(std::vector<Eigen::Vector3d> ptsCord,bool bIsY);
	
	Eigen::Matrix4d GetMatrixToPipeForH(Eigen::Vector3d Pt1, Eigen::Vector3d ptNor1, Eigen::Vector3d Pt2, Eigen::Vector3d ptNor2,double vLen1,double vLen2);
	Eigen::Matrix4d GetMatrixToPipeForB(Eigen::Vector3d Pt1, Eigen::Vector3d ptNor1, Eigen::Vector3d Pt2, Eigen::Vector3d ptNor2,double vLen1,double vLen2);

	Eigen::Matrix4d GetMatrix3PtsCoord(Eigen::Vector3d ptOrg, Eigen::Vector3d ptX, Eigen::Vector3d ptXY);
	
	static double DistPtToPlane(Eigen::Vector3d ptOrgPlane, Eigen::Vector3d vNorPlane, Eigen::Vector3d ptThe, Eigen::Vector3d&intP);
	static double DistPtToPlane(double a, double b, double c, double d, Eigen::Vector3d vec, Eigen::Vector3d theP, Eigen::Vector3d&intP);

public:
	static void MovePts(Eigen::Matrix4d tMat, std::vector<Eigen::Vector3d>& arrPts);
	static void MovePts(Eigen::Matrix4d tMat, std::vector<Eigen::Vector3f>& arrPts);
	static void MovePts(Eigen::Matrix4d tMat, std::vector<std::vector<Eigen::Vector3d>>& arrPtss);
};
