#include "Sphere.h"
#include <time.h>

//#include "open3d/Open3D.h"

#define SQUARE(x) ((x) * (x))

// 1.0E-1
const double TOL_ZERO_1 = 1.0E-1;


void Sphere::FittingSphere(vector<Eigen::Vector3d> pts,
                           Eigen::Vector3d& ptCenter,
                           double& vRad) {
    int N = (int)pts.size();
    double x, y, z;
    double Sx, Sy, Sz;
    double Sxx, Syy, Szz, Sxy, Sxz, Syz;    
    double Sxxx, Syyy, Szzz, Sxyy, Sxzz, Sxxy, Sxxz, Syyz, Syzz;
    x = y = z = Sx = Sy = Sz = 0;
    Sxx = Syy = Szz = Sxy = Sxz = Syz = 0;
    Sxxx = Syyy = Szzz = Sxyy = Sxzz = Sxxy = Sxxz = Syyz = Syzz = 0;
    for (int i = 0; i < N; i++) {
        x = pts[i].x();
        y = pts[i].y();
        z = pts[i].z();

        Sx += x;
        Sy += y;
        Sz += z;

        Sxx += (x * x);
        Syy += (y * y);
        Szz += (z * z);
        Sxy += (x * y);
        Sxz += (x * z);
        Syz += (y * z);

		Sxxx += (x * x * x);
        Syyy += (y * y * y);
        Szzz += (z * z * z);
        Sxyy += (x * y * y);
        Sxzz += (x * z * z);
        Sxxy += (x * x * y);
        Sxxz += (x * x * z);
        Syyz += (y * y * z);
        Syzz += (y * z * z);
    }

	double A1 = Sxx + Syy + Szz;
    double a = 2 * Sx * Sx - 2 * (double)N * Sxx;
    double b = 2 * Sx * Sy - 2 * (double)N * Sxy;
    double c = 2 * Sx * Sz - 2 * (double)N * Sxz;
    double d = -N * (Sxxx + Sxyy + Sxzz) + A1 * Sx;
    double e = 2 * Sx * Sy - 2 * (double)N * Sxy;
    double f = 2 * Sy * Sy - 2 * (double)N * Syy;
    double g = 2 * Sy * Sz - 2 * (double)N * Syz;
    double h = -N * (Sxxy + Syyy + Syzz) + A1 * Sy;
    double j = 2 * Sx * Sz - 2 * (double)N * Sxz;
    double k = 2 * Sy * Sz - 2 * (double)N * Syz;
    double l = 2 * Sz * Sz - 2 * (double)N * Szz;
    double m = -N * (Sxxz + Syyz + Szzz) + A1 * Sz;
    double delta = a * (f * l - g * k) - e * (b * l - c * k) + j * (b * g - c * f);
    double xc = (d * (f * l - g * k) - h * (b * l - c * k) + m * (b * g - c * f)) / delta;
    double yc = (a * (h * l - m * g) - e * (d * l - m * c) + j * (d * g - h * c)) / delta;
    double zc = (a * (f * m - h * k) - e * (b * m - d * k) + j * (b * h - d * f)) / delta;
    double R = sqrt(xc * xc + yc * yc + zc * zc +(A1 - 2 * (xc * Sx + yc * Sy + zc * Sz)) / N);

    vRad = R;
    ptCenter.x() = xc;
    ptCenter.y() = yc;
    ptCenter.z() = zc;
}

void Sphere::FittingSphereGivenRadius(vector<Eigen::Vector3d> pts,
                                      const double vGivenRadius,
                                      Eigen::Vector3d& ptCenter) {
    Eigen::Vector3d ptCenter0, ptCenterPrev;
    double vRadius0;

    FittingSphere(pts, ptCenter0, vRadius0);

  	int nCnt = 0;
    double vTol = vRadius0 / 10000;
    int size = (int)pts.size();

    do {
        ptCenterPrev = ptCenter0;
        vector<Eigen::Vector3d> ptsCenter;

        for (int i = 0; i < size; i++) 
        {
            //Eigen::Vector3d ptNew = pts->points_[i] +(ptCenter0 - pts->points_[i]).UnitV() * vGivenRadius;
            Eigen::Vector3d ptNew = pts[i] + (ptCenter0 - pts[i]).normalized() * vGivenRadius;
            ptsCenter.push_back(ptNew);
        }

        ptCenter0 = GetCenter(ptsCenter);
        
        double dis = Dist(ptCenter0, ptCenterPrev);
        if (dis < vTol) break;

        nCnt++;
        if (nCnt > 100) break;
    } while (1);
}

Eigen::Vector3d Sphere::GetCenter(vector<Eigen::Vector3d> arrPt) 
{
    Eigen::Vector3d cen(0.f, 0.f, 0.f);

    if (arrPt.size() == 0) return cen;

    int size = arrPt.size();

    for (size_t i = 0; i < size; i++) {
        cen += arrPt[i];
    }

    return cen / float(arrPt.size());
}

double Sphere::Dist(Eigen::Vector3d& v1, Eigen::Vector3d& v2) 
{
    return sqrt(SQUARE(v1.x() - v2.x()) + 
                SQUARE(v1.y() - v2.y()) +
                SQUARE(v1.z() - v2.z()));
}

void Sphere::FindSphereWithUniformSegmentation(std::shared_ptr<open3d::geometry::PointCloud> pcd,
                                               vector<Eigen::Vector3d> &ptsCenter,
                                               double vFindingRadius1,
                                               double vFindingRadius2,
                                               double vGiveRadius,
                                               bool bCreate,
                                               double vStdTol) {
    double vMinR = min(vFindingRadius1, vFindingRadius2);
    double vMaxR = max(vFindingRadius1, vFindingRadius2);

    if (vMinR < 1) vMinR = 1;
    if (pcd == nullptr) return;
    if (pcd->points_.size() == 0) return;

    vector<Eigen::Vector3d>*pts = &pcd->points_;

    std::vector<double> vRadiuss0;
    std::vector<double> vRMS0;

    FindSphereWithUniformSegmentation(pts, vMinR, vMaxR, ptsCenter,
                                      vRadiuss0, vRMS0, vGiveRadius, bCreate,
                                      vStdTol);
}

void Sphere::FindSphereWithUniformSegmentation(
        vector<Eigen::Vector3d>* pts,
        double vMinR,
        double vMaxR,
        vector<Eigen::Vector3d>& ptsCenter,
        vector<double>& vRadiuss,
        vector<double>& vRMSs,
        double vGivenRadius,
        bool bCreate,
        double vStdTol) 
{
    if (pts->size() < 10) return;

    clock_t t = clock();

    FindSphereWithUniformSegmentation(pts, vMinR, vMaxR, ptsCenter, vRadiuss,
                                      vRMSs, vGivenRadius, vStdTol);

    if (ptsCenter.empty()) return;

    //kykernel implemented drawing part at here

    t = abs(t - clock());
    open3d::utility::LogInfo("elapsed time - {} sec", ((float)t / CLOCKS_PER_SEC));
    int i = 0;
    for(auto pt : ptsCenter) {
        open3d::utility::LogInfo("center xyz : {} {} {}, radius : {}, rms : {}", pt.x(), pt.y(), pt.z(), vRadiuss[i], vRMSs[i]);
        i++;
    }
    
}

void Sphere::FindSphereWithUniformSegmentation(
        vector<Eigen::Vector3d>* pts,
        double vMinR,
        double vMaxR,
        vector<Eigen::Vector3d>& ptsCenter,
        std::vector<double>& vRadiuss,
        std::vector<double>& vRMSs,
        double vGivenRadius,
        double vStdTol) {

    V3D box;
    //pts.GetBoundingBox(box);

    int size = pts->size();
    box.Min = {0, 0, 0};
    box.Max = {0, 0, 0};

    for(int i=0; i<size; i++) {
        box.Update(pts->at(i));
    }

    double vRefSize = vMinR * 2.1;
    if (vRefSize < 5) vRefSize = 5;

    if (vRefSize > box.MaxWidth() / 2.01)  //�ּ� 2ĭ�� ����� ������
        vRefSize = box.MaxWidth() / 2.01;
    if (vRefSize < box.MaxWidth() / 100) vRefSize = box.MaxWidth() / 100;

    //���� ������ �ڷ� ������ ���� �����
    CKYScanArea* pArea = CKYScanArea::MakeArea(pts, vRefSize);

    //
    FindSpheres(pArea, pts, vMinR, vMaxR, ptsCenter, vRadiuss, vRMSs,
                vGivenRadius, vStdTol);

    delete pArea;
}

double Sphere::GetRadiusOfSphere(vector<Eigen::Vector3d>& pts, Eigen::Vector3d ptCenter) 
{
    double vRadius = 0;
    int size = pts.size();
    for (int i = 0; i < size; i++) {
        vRadius += Dist(ptCenter, pts[i]);
    }

    vRadius /= pts.size();
    return vRadius;
}


void Sphere::MakeSphereData(vector<Eigen::Vector3d>& pts,
                            Eigen::Vector3d& ptCenter,
                                   double& vRadius,
                                   double& vErrStd,
                                   double vRmvRatio,
                                   double vGivenRadius,
                                   bool bMessage) {
    vRadius = 0;
    ptCenter = {0,0,0};
    if (pts.size() < 4) return;
    
    ptCenter = GetCenter(pts);
    vRadius = GetRadiusOfSphere(pts, ptCenter);

    vector<Eigen::Vector3d> ptsNew;
    if (vGivenRadius < TOL_ZERO_1)
        FittingSphere(pts, ptCenter, vRadius);
    else {
        FittingSphereGivenRadius(pts, vGivenRadius, ptCenter);
        vRadius = vGivenRadius;
    }

    double vRadius0 = vRadius;
    if (vRadius > 1e8) {
        vRadius = 0;
        ptCenter = {0, 0, 0};
        vErrStd = 100000;
        return;
    }

    if (vRmvRatio >= 0.01) {
        int N = (int)(vRmvRatio / 0.1);
        if (N > 4) N = 4;
        if (N < 1) N = 1;
        double vRatio = vRmvRatio / N;
        for (int i = 0; i < N; i++) {
            if (i == 0) {
                ReMakeSphereData(pts, ptsNew, vRatio, vRadius, ptCenter, vGivenRadius);
            } else{
                ReMakeSphereData(ptsNew, ptsNew, vRatio, vRadius, ptCenter,vGivenRadius);
            }
                
        }
    } else
        ptsNew = pts;

    std::vector<double> vErrs;
    vErrs.resize(ptsNew.size());
    for (int i = 0; i < (int)ptsNew.size(); i++) {
        double len = Dist(ptCenter,ptsNew[i]);
        vErrs[i] = vRadius - len;
    }

    if (vRadius > 1e8) {
        vRadius = 0;
        ptCenter = {0, 0, 0};
        vErrStd = 100000;
        return;
    }

    Statistics stat;
    double vAve;
    stat.CalcStatistics(vErrs, vErrStd, vAve);
    //open3d::utility::LogInfo("cx:{}, cy:{}, cz:{}, radius:{}, standard error:{}", ptCenter.x(),ptCenter.y(), ptCenter.z(), vRadius, vErrStd);
}

void Sphere::ReMakeSphereData(vector<Eigen::Vector3d> pts,
                              vector<Eigen::Vector3d>& ptsNew,
                                     double vRemoveRatio,
                                     double& vRadius,
                              Eigen::Vector3d& ptCenter,
                                     double vGivenRadius) {
    int size = pts.size();
    if (size < 4) return;

    ptsNew.clear();

    std::vector<std::pair<double, Eigen::Vector3d>> ptsDis;
    ptsDis.resize(size);
    for (int i = 0; i < size; i++) {
        Eigen::Vector3d pt = ptCenter + (pts[i] - ptCenter).normalized() * vRadius;
        double vErr = Dist(pts[i] ,pt);
        ptsDis[i] = std::make_pair(vErr, pts[i]);
    }

    struct {
        bool operator()(std::pair<double, Eigen::Vector3d>& it1,
                        std::pair<double, Eigen::Vector3d>& it2) const {
            return fabs(it1.first) < fabs(it2.first);
        }
    } disPtCmp;

    std::sort(ptsDis.begin(), ptsDis.end(), disPtCmp);

    int Nr = int(pts.size() * (1 - vRemoveRatio));
    if (Nr < 3) Nr = 3;

    ptsNew.resize(Nr);
    for (int i = 0; i < Nr; i++) ptsNew[i] = ptsDis[i].second;

    if (vGivenRadius < TOL_ZERO_1)
        FittingSphere(ptsNew, ptCenter, vRadius);
    else {
        FittingSphereGivenRadius(ptsNew, vGivenRadius, ptCenter);
        vRadius = vGivenRadius;
    }
}

void Sphere::FindInitSpheres(CKYScanArea* pArea,
                                          double vMinR,
                                          double vMaxR,
                                          int nPtTol,
                                          double vErrTol,
                                          vector<Eigen::Vector3d>& ptCenters,
                                          std::vector<double>& vRadiuss) 
{
    vMinR *= 0.9;
    vMaxR *= 1.1;
    nPtTol = int(nPtTol * 0.5);

    int NX = pArea->GetScanAreaAttrib()->GetNX();
    int NY = pArea->GetScanAreaAttrib()->GetNY();
    int NZ = pArea->GetScanAreaAttrib()->GetNZ();

    std::vector<CKyScanSegmentData*> pSegmentDatas;
    int nFile = 0;
    //#pragma omp parallel
    {
        //#pragma omp for  schedule(static)
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {
                for (int k = 0; k < NZ; k++) {
                    CScanGrid* pGrid =
                            pArea->GetScanAreaAttrib()->GetGrid(i, j, k);

                    int nStep = 1;

                    nStep = (int)(pGrid->GetNPts() / 1000.);
                    if (nStep < 1) nStep = 1;

                    vector<Eigen::Vector3d> pts;
                    CSGScanAreaAttrib::GetPartialGridScanPts(pGrid, &pts, nStep);

                    if (pts.size() < nPtTol / nStep) continue;

                    Eigen::Vector3d ptCenter;
                    double vRadius;
                    double vErrStd;
                    MakeSphereData(pts, ptCenter, vRadius, vErrStd, 0.3);  

                    bool bOK = false;

                    if (vErrStd < vErrTol / 3 && vRadius > vMinR * 0.8 &&
                        vRadius < vMaxR * 1.2)
                        bOK = true;
                    else if (vErrStd < vErrTol / 2 && vRadius > vMinR * 0.9 &&
                             vRadius < vMaxR * 1.1)
                        bOK = true;
                    else if (vErrStd < vErrTol && vRadius > vMinR &&
                             vRadius < vMaxR)
                        bOK = true;

                    if (bOK) {
                        double vPlaneRMS = GetPlaneRMS(pts, 0.3);
                        if (vPlaneRMS < vErrStd) continue;

                        CKyScanSegmentData* pData = new CKyScanSegmentData;
                        pData->m_Center = ptCenter;
                        pData->m_Radius = vRadius;
                        pData->m_NPts = (int)pts.size() * nStep;
                        pData->m_vErr = vErrStd;
                        pSegmentDatas.push_back(pData);
                    }
                }
            }
        }
    }

    if (pSegmentDatas.empty()) return;

    std::vector<double> vRMSs;
    MergeSimilarSpheres(pSegmentDatas, nPtTol, ptCenters, vRadiuss, vRMSs);
}

double Sphere::GetPlaneRMS(vector<Eigen::Vector3d>& pts, double vRatio) {
    Eigen::Vector3d ptOrg, vNormal;
    vector<V3D> plane_pts;
    KyMath math;

    math.MakePlane(pts, ptOrg, vNormal);

    std::vector<std::pair<double, Eigen::Vector3d>> ptsDis;
    ptsDis.resize(pts.size());
    Eigen::Vector3d ptInt;
    for (int i = 0; i < (int)pts.size(); i++) {
        double dis = math.DistPtToPlane(ptOrg, vNormal, pts[i], ptInt);
        if (Dist(pts[i], ptOrg) < Dist(pts[i] , ptOrg + vNormal)) 
            dis *= -1;
        ptsDis[i] = std::make_pair(dis, pts[i]);
    }

    std::sort(ptsDis.begin(), ptsDis.end(), CKyDisPtComparer());

    int Nr = int(pts.size() * (1 - vRatio));
    if (Nr < 3) Nr = 3;

    vector<Eigen::Vector3d> ptsNew;
    ptsNew.resize(Nr);
    for (int i = 0; i < Nr; i++) ptsNew[i] = ptsDis[i].second;

    math.MakePlane(ptsNew, ptOrg, vNormal);
    std::vector<double> vErrs;
    for (int i = 0; i < (int)ptsNew.size(); i++) {
        double dis = math.DistPtToPlane(ptOrg, vNormal, ptsNew[i], ptInt);
        if (Dist(ptsNew[i], ptInt) < Dist(ptsNew[i] , ptInt + vNormal)) dis *= -1;
        vErrs.push_back(dis);
    }

    double vStd, vAve;
    Statistics stat;
    stat.CalcStatistics(vErrs, vStd, vAve);
    return vStd;
}

bool Sphere::IsSphereWithArea(vector<Eigen::Vector3d>& ptsAll,
                              Eigen::Vector3d ptCenter,
                              double vRadius,
                              double vRatioTol)
{
    int N1 = 10;
    int N2 = 20;
    bool** bMesh = new bool*[N1];

    for (int i = 0; i < N1; i++) bMesh[i] = new bool[N2];

    for (int i = 0; i < N1; i++) {
        for (int j = 0; j < N2; j++) bMesh[i][j] = false;
    }

    double vErrTol = vRadius / 10;
    for (int i = 0; i < (int)ptsAll.size(); i++) {
        Eigen::Vector3d pt = ptsAll[i];

        pt -= ptCenter;

        double vErr = fabs(vRadius - Size(pt));
        if (vErr > vErrTol) continue;

        double ang = atan2(pt.z(), pt.y());
        if (ang < 0) ang += 2 * EIGEN_PI;
        double v2 = ang / (2 * EIGEN_PI) * N2;

        Eigen::AngleAxisd rollAngle(-ang, Eigen::Vector3d::UnitX());
        Eigen::AngleAxisd pitchAngle(0, Eigen::Vector3d::UnitY());
        Eigen::AngleAxisd yawAngle(0, Eigen::Vector3d::UnitZ());
        Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;

        Eigen::Matrix3d tRotX = q.matrix();
        //KyTMatrix tRotX;
        //tRotX.RotateX(-ang);
        Eigen::Vector3d pt0 = tRotX * pt;
        ang = atan2(pt0.y(), pt0.x());
        double v1 = ang / (EIGEN_PI)*N1;

        int n1 = int(v1);
        int n2 = int(v2);

        bMesh[n1][n2] = true;
    }

    double vTotalArea = 0;
    double vArea = 0;
    for (int i = 0; i < N1; i++) {
        double area;
        if (i == 0 || i == 9)
            area = 0.0150663;
        else if (i == 1 || i == 8)
            area = 0.0438402;
        else if (i == 2 || i == 7)
            area = 0.068412;
        else if (i == 3 || i == 6)
            area = 0.0862433;
        else if (i == 4 || i == 5)
            area = 0.095589;

        for (int j = 0; j < N2; j++) {
            if (bMesh[i][j]) {
                vArea += area;
            }
            vTotalArea += area;
        }
    }

    double vRatio = vArea / vTotalArea;

    if (vRatio > vRatioTol) return true;

    return false;
}

bool Sphere::MakeFinalSpheres(vector<vector<Eigen::Vector3d>>& ptsSphere,
                              double vMinR,
                              double vMaxR,
                              int nPtTol,
                              double vStdTol,
                              vector<Eigen::Vector3d>& ptCenters,
                              std::vector<double>& vRadiuss,
                              std::vector<double>& vRMSs,
                              double vGivenRadius)
{
    vMinR *= 0.98;
    vMaxR *= 1.02;

    std::vector<CKyScanSegmentData*> pSegmentDatas;
    ptCenters.clear();
    vRadiuss.clear();
    vRMSs.clear();

    double vRatioTol = 0.15;

    for (int i = 0; i < (int)ptsSphere.size(); i++) {
        if (ptsSphere[i].size() < nPtTol) continue;

        Eigen::Vector3d ptCenter;
        double vRadius;
        double vErrRMS;
        MakeSphereData(ptsSphere[i], ptCenter, vRadius, vErrRMS, 0.4, vGivenRadius,false);

        if (vErrRMS < vStdTol && vRadius > vMinR && vRadius < vMaxR) {
            if (!IsSphereWithArea(ptsSphere[i], ptCenter, vRadius, vRatioTol))
                continue;

            if (vErrRMS > vStdTol / 2 && !IsSphereWithPartSphere(ptsSphere[i]))
                continue;

            CKyScanSegmentData* pData = new CKyScanSegmentData;
            pData->m_Center = ptCenter;
            if (vGivenRadius == 0)
                pData->m_Radius = vRadius;
            else
                pData->m_Radius = vGivenRadius;

            pData->m_NPts = (int)ptsSphere[i].size();
            pData->m_vErr = vErrRMS;

            pSegmentDatas.push_back(pData);
        }
    }

    if (pSegmentDatas.empty()) return true;

    int nSeg = (int)pSegmentDatas.size();
    MergeSimilarSpheres(pSegmentDatas, nPtTol, ptCenters, vRadiuss, vRMSs);

    if (nSeg > (int)ptCenters.size()) return false;

    return true;
}

bool Sphere::IsSphereWithPartSphere(vector<Eigen::Vector3d>& ptsAll) {
    vector<Eigen::Vector3d> ptsPart1;
    vector<Eigen::Vector3d> ptsPart2;
    vector<Eigen::Vector3d> ptsAll2;

    //�ұ����ϰ� ����Ʈ�� ���ؼ� �Ǵ�
    int n = (int)(ptsAll.size() / 3);
    int ns = (int)(ptsAll.size() / 5);
    for (int i = 0; i < (int)ptsAll.size(); i++) {
        if (i < n && i % 3 == 0)
            ptsPart2.push_back(ptsAll[i]);
        else if (i > n && i < n * 2 && i % 7 == 0)
            ptsPart2.push_back(ptsAll[i]);
        else if (i > n * 2 && i % 11 == 0)
            ptsPart2.push_back(ptsAll[i]);

        if (i > ns && i % 4 == 0) ptsPart1.push_back(ptsAll[i]);

        if (i % 3 == 0) ptsAll2.push_back(ptsAll[i]);
    }

    if (ptsPart1.size() < 10 || ptsPart2.size() < 10) return true;

    //�������� �־��� ��� �ٽ� ������ ����� ���ϸ� ������ Ŀ�� �� ����
    //�ٽ� ������ ����� ��
    Eigen::Vector3d ptCenter;
    double vRadius, vErrStd;
    MakeSphereData(ptsAll2, ptCenter, vRadius, vErrStd, 0.4);

    double vErrTol = vRadius;
    double vTol = vRadius / 10;

    Eigen::Vector3d ptCenterNew;
    double vRadiusNew;

    // 33% ���ø��� ����� ���� ����� �ſ� �����ϸ� ���� �Ǵ�
    MakeSphereData(ptsPart1, ptCenterNew, vRadiusNew, vErrStd, 0.4);
    double dis = Dist(ptCenter,ptCenterNew);
    if (fabs(vRadius - vRadiusNew) < vTol / 2 && dis < vTol / 2) return true;

    //�б����ϰ� ���ø��� ����� ���̰� ���ϸ� ���� �ƴ�
    MakeSphereData(ptsPart2, ptCenterNew, vRadiusNew, vErrStd,0.4);
    dis = Dist(ptCenter,ptCenterNew);
    if (fabs(vRadius - vRadiusNew) > vTol || dis > vTol) return false;

    if (fabs(vRadius - vRadiusNew) + dis > vTol * 1.5) return false;

    return true;
}

void Sphere::FindSpheres(CKYScanArea* pArea,
                                      vector<Eigen::Vector3d>* ptsAll,
                                      double vMinR,
                                      double vMaxR,
                                      vector<Eigen::Vector3d>& ptCenters,
                                      vector<double>& vRadiuss,
                                      vector<double>& vRMSs,
                                      double vGivenRadius,
                                      double vStdTol) {
    vector<Eigen::Vector3d> ptCenters0;
    vector<double> vRadiuss0;

    int NX = pArea->GetScanAreaAttrib()->GetNX();
    int NY = pArea->GetScanAreaAttrib()->GetNY();
    int NZ = pArea->GetScanAreaAttrib()->GetNZ();
    int NPts = pArea->GetScanAreaAttrib()->GetNPts();

    int nPtTol = (int)(NPts / (NX * NY * NZ * 1.3));
    if (nPtTol < 100) nPtTol = 100;
    if (nPtTol > 3000) nPtTol = 3000;
    FindInitSpheres(pArea, vMinR, vMaxR, nPtTol, vStdTol * 1.5, ptCenters0, vRadiuss0);
    if (ptCenters0.empty()) return;

    vector<vector<Eigen::Vector3d>> ptsSphere;
    FindPtsOfSpheres(ptsAll, ptCenters0, vRadiuss0, nPtTol, ptsSphere);

    if (!MakeFinalSpheres(ptsSphere, vMinR, vMaxR, nPtTol, vStdTol, ptCenters,
                          vRadiuss, vRMSs, vGivenRadius)) {
        //�������� ���� �������� ��찡 �߻��ϸ�, �� �� ��Ȯ�ϰ� �ϱ� ���ؼ� �ٽ�
        //����Ʈ�� ã�Ƽ� ������
        vector<vector<Eigen::Vector3d>> ptsSphere;
        FindPtsOfSpheres(ptsAll, ptCenters, vRadiuss, nPtTol, ptsSphere);

        MakeFinalSpheres(ptsSphere, vMinR, vMaxR, nPtTol, vStdTol, ptCenters,
                         vRadiuss, vRMSs, vGivenRadius);
    }
}

void Sphere::MergeSimilarSpheres(
        std::vector<CKyScanSegmentData*>& pSegmentDatas,
        int nPtsTol,
        vector<Eigen::Vector3d>& ptCenters,
        std::vector<double>& vRadiuss,
        std::vector<double>& vRMSs) {
    if (pSegmentDatas.empty()) return;

    double vAveRadius = 0;
    for (int i = 0; i < (int)pSegmentDatas.size() - 1; i++)
        vAveRadius += pSegmentDatas[i]->m_Radius;
    vAveRadius /= pSegmentDatas.size();

    for (int i = 0; i < (int)pSegmentDatas.size() - 1; i++) {
        CKyScanSegmentData* pData1 = pSegmentDatas[i];
        for (int j = i + 1; j < (int)pSegmentDatas.size(); j++) {
            CKyScanSegmentData* pData2 = pSegmentDatas[j];

            if (Dist(pData1->m_Center, pData2->m_Center) < vAveRadius) {
                if (fabs(pData1->m_Radius - pData2->m_Radius) <
                    vAveRadius / 4) {
                    if (pData1->m_vErr < pData1->m_vErr / 3)
                        ;
                    else if (pData2->m_vErr < pData1->m_vErr / 3) {
                        pData1->m_Center = pData2->m_Center;
                        pData1->m_Radius = pData2->m_Radius;
                        pData1->m_NPts = pData2->m_NPts;
                        pData1->m_vErr = pData2->m_vErr;
                    } else {
                        pData1->m_Center = (pData1->m_Center * pData1->m_NPts +
                                            pData2->m_Center * pData2->m_NPts) /
                                           (pData1->m_NPts + pData2->m_NPts);
                        pData1->m_Radius = (pData1->m_Radius * pData1->m_NPts +
                                            pData2->m_Radius * pData2->m_NPts) /
                                           (pData1->m_NPts + pData2->m_NPts);
                        pData1->m_vErr = (pData1->m_vErr * pData1->m_NPts +
                                          pData2->m_vErr * pData2->m_NPts) /
                                         (pData1->m_NPts + pData2->m_NPts);
                        pData1->m_NPts = pData1->m_NPts + pData2->m_NPts;
                    }

                    delete pSegmentDatas[j];
                    pSegmentDatas.erase(pSegmentDatas.begin() + j);
                    j--;
                }
            }
        }
    }
    for (int i = 0; i < (int)pSegmentDatas.size() - 1; i++) {
        CKyScanSegmentData* pData1 = pSegmentDatas[i];
        for (int j = i + 1; j < (int)pSegmentDatas.size(); j++) {
            CKyScanSegmentData* pData2 = pSegmentDatas[j];

            if (Dist(pData1->m_Center, pData2->m_Center) < vAveRadius) {
                if (fabs(pData1->m_Radius - pData2->m_Radius) <
                    vAveRadius / 4) {
                    if (pData1->m_vErr < pData1->m_vErr / 3)
                        ;
                    else if (pData2->m_vErr < pData1->m_vErr / 3) {
                        pData1->m_Center = pData2->m_Center;
                        pData1->m_Radius = pData2->m_Radius;
                        pData1->m_NPts = pData2->m_NPts;
                        pData1->m_vErr = pData2->m_vErr;
                    } else {
                        pData1->m_Center = (pData1->m_Center * pData1->m_NPts +
                                            pData2->m_Center * pData2->m_NPts) /
                                           (pData1->m_NPts + pData2->m_NPts);
                        pData1->m_Radius = (pData1->m_Radius * pData1->m_NPts +
                                            pData2->m_Radius * pData2->m_NPts) /
                                           (pData1->m_NPts + pData2->m_NPts);
                        pData1->m_vErr = (pData1->m_vErr * pData1->m_NPts +
                                          pData2->m_vErr * pData2->m_NPts) /
                                         (pData1->m_NPts + pData2->m_NPts);
                        pData1->m_NPts = pData1->m_NPts + pData2->m_NPts;
                    }

                    delete pSegmentDatas[j];
                    pSegmentDatas.erase(pSegmentDatas.begin() + j);
                    j--;
                }
            }
        }
    }

    for (int i = 0; i < (int)pSegmentDatas.size(); i++) {
        CKyScanSegmentData* pData1 = pSegmentDatas[i];
        if (pData1->m_NPts > nPtsTol) {
            ptCenters.push_back(pData1->m_Center);
            vRadiuss.push_back(pData1->m_Radius);
            vRMSs.push_back(pData1->m_vErr);
        }
    }

    for (int i = 0; i < (int)pSegmentDatas.size(); i++) delete pSegmentDatas[i];
}

void Sphere::FindPtsOfSpheres(vector<Eigen::Vector3d>* ptsAll,
                              vector<Eigen::Vector3d> ptCenters,
                                           std::vector<double> vRadiuss,
                                           int nPtTol,
                              vector<vector<Eigen::Vector3d>>& ptsSphere) 
{
    for (int i = 0; i < ptCenters.size(); i++) {
        vector<Eigen::Vector3d> ptsF;
        ptsSphere.push_back(ptsF);
    }

    int nStep = 1;
    // 	nStep = (int)(ptsAll.size() / 100000.);
    // 	if (nStep < 1) nStep = 1;

    std::vector<int> nPtsRange1;
    std::vector<int> nPtsRange2;
    for (int i = 0; i < (int)ptCenters.size(); i++) {
        nPtsRange1.push_back(0);
        nPtsRange2.push_back(1);
    }

    for (int i = 0; i < (int)ptsAll->size(); i++) {
        if (i % nStep != 0) continue;

        for (int j = 0; j < (int)ptCenters.size(); j++) {
            double dis = fabs(Dist(ptsAll->at(i) ,ptCenters[j]) - vRadiuss[j]);
            if (dis < vRadiuss[j] / 5) ptsSphere[j].push_back(ptsAll->at(i));
        }
    }

    for (int i = 0; i < (int)ptsSphere.size(); i++) {
        if (ptsSphere[i].size() < nPtTol) {
            ptsSphere.erase(ptsSphere.begin() + i);
            i--;
        }
    }
}
