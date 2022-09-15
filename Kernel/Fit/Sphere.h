#pragma once

#include <vector>
#include "open3d/Open3D.h"
#include "../Base/Statistics.h"
#include "../Base/KYScanArea.h"
#include "../Base/KyMath.h"
using namespace std;
#define KySQR(x) ((x) * (x))

class CKyScanSegmentData {
public:
    CKyScanSegmentData() {
        m_Center = {0, 0, 0};
        m_Normal = {0, 0, 0};
        m_ptSt = {0, 0, 0};
        m_ptEd = {0, 0, 0};
        m_Radius = 0;
        m_Flag = false;
        m_RefDistance = 0;
        m_NPts = 0;
        m_vErr = 0;
    }
    ~CKyScanSegmentData() {}

public:
    Eigen::Vector3d m_Center;
    Eigen::Vector3d m_Normal;
    double m_Radius;

    int m_NPts;
    double m_vErr;

    Eigen::Vector3d m_ptSt, m_ptEd;

    double m_RefDistance;

    bool m_Flag;
};

class Sphere{
public:
    void FindSphereWithUniformSegmentation(std::shared_ptr<open3d::geometry::PointCloud> pts,
                                           vector<Eigen::Vector3d>& ptsCenter,
                                           double vFindingRadius1,
                                           double vFindingRadius2,
                                           double vGiveRadius,
                                           bool bCreate = false,
                                           double vStdTol = 4);
    
    void FindSphereWithUniformSegmentation(vector<Eigen::Vector3d>* pts,
                                           double vMinR,
                                           double vMaxR,
                                           vector<Eigen::Vector3d>& ptsCenter,
                                           vector<double>& vRadiuss,
                                           vector<double>& vRMSs,
                                           double vGivenRadius,
                                           bool bCreate,
                                           double vStdTol);

    void FindSphereWithUniformSegmentation(vector<Eigen::Vector3d>* pts,
                                           double vMinR,
                                           double vMaxR,
                                           vector<Eigen::Vector3d>& ptsCenter,
                                           std::vector<double>& vRadiuss,
                                           std::vector<double>& vRMSs,
                                           double vGivenRadius,
                                           double vStdTol);

private:
    double SqDist(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const {
        return (KySQR(v1.x() - v2.x()) + KySQR(v1.y() - v2.y()) + KySQR(v1.z() - v2.z()));
    }

    double SqDistXY(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const {
        return (KySQR(v1.x() - v2.x()) + KySQR(v1.y() - v2.y()));
    }

    double SqDistYZ(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const {
        return (KySQR(v1.y() - v2.y()) + KySQR(v1.z() - v2.z()));
    }

    double SqDistZX(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const {
        return (KySQR(v1.z() - v2.z()) + KySQR(v1.x() - v2.x()));
    }

    double Dist(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const { return sqrt(SqDist(v1,v2)); }
    double DistXY(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const { return sqrt(SqDistXY(v1, v2)); }
    double DistYZ(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const { return sqrt(SqDistYZ(v1, v2)); }
    double DistZX(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const { return sqrt(SqDistZX(v1, v2)); }

    double SqSize(const Eigen::Vector3d& v) const {return (KySQR(v.x()) + KySQR(v.y()) + KySQR(v.z()));}
    double Size(const Eigen::Vector3d& v) const {return sqrt(SqSize(v)); }
    

private:
    double GetRadiusOfSphere(vector<Eigen::Vector3d>& pts,
                                     Eigen::Vector3d ptCenter);
    double Dist(Eigen::Vector3d& v1, Eigen::Vector3d& v2); 
    
    Eigen::Vector3d GetCenter(vector<Eigen::Vector3d> arrPt);
    
    void FittingSphere(vector<Eigen::Vector3d> pts,
                       Eigen::Vector3d& ptCenter,
                       double& vRad);
    
    void FittingSphereGivenRadius(vector<Eigen::Vector3d> pts,
                                  const double vGivenRadius,
                                  Eigen::Vector3d& ptCenter);
    
    void MakeSphereData(vector<Eigen::Vector3d>& pts,
                        Eigen::Vector3d& ptCenter,
                        double& vRadius,
                        double& vErrStd,
                        double vRmvRatio=0.4,
                        double vGivenRadius=0,
                        bool bMessage=false);

    void ReMakeSphereData(vector<Eigen::Vector3d> pts,
                          vector<Eigen::Vector3d>& ptsNew,
                          double vRemoveRatio,
                          double& vRadius,
                          Eigen::Vector3d& ptCenter,
                          double vGivenRadius);

    void MergeSimilarSpheres(
            std::vector<CKyScanSegmentData*>& pSegmentDatas,
            int nPtsTol,
            vector<Eigen::Vector3d>& ptCenters,
            std::vector<double>& vRadiuss,
            std::vector<double>& vRMSs);

    void FindInitSpheres(CKYScanArea* pArea,
                         double vMinR,
                         double vMaxR,
                         int nPtTol,
                         double vErrTol,
                         vector<Eigen::Vector3d>& ptCenters,
                         std::vector<double>& vRadiuss);
    
    void FindSpheres(CKYScanArea* pArea,
                     vector<Eigen::Vector3d>* ptsAll,
                     double vMinR,
                     double vMaxR,
                     vector<Eigen::Vector3d>& ptCenters,
                     vector<double>& vRadiuss,
                     vector<double>& vRMSs,
                     double vGivenRadius,
                     double vStdTol);

    void FindPtsOfSpheres(vector<Eigen::Vector3d>* ptsAll,
                         vector<Eigen::Vector3d> ptCenters,
                         std::vector<double> vRadiuss,
                         int nPtTol,
                         vector<vector<Eigen::Vector3d>>& ptsSphere);

    bool MakeFinalSpheres(vector<vector<Eigen::Vector3d>>& ptsSphere,
                          double vMinR,
                          double vMaxR,
                          int nPtTol,
                          double vStdTol,
                          vector<Eigen::Vector3d>& ptCenters,
                          std::vector<double>& vRadiuss,
                          std::vector<double>& vRMSs,
                          double vGivenRadius=0);

    bool IsSphereWithArea(vector<Eigen::Vector3d>& ptsAll,
                          Eigen::Vector3d ptCenter,
                          double vRadius,
                          double vRatioTol);

    bool IsSphereWithPartSphere(vector<Eigen::Vector3d>& ptsAll);
    
    double GetPlaneRMS(vector<Eigen::Vector3d>& pts, double vRatio);
};

class CKyDisPtComparer {
public:
    bool operator()(std::pair<double, Eigen::Vector3d>& it1,
                    std::pair<double, Eigen::Vector3d>& it2) {
        return fabs(it1.first) < fabs(it2.first);
    }
};