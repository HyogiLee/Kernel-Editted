#include "KyMath.h"
#include "Kyprimfit.h"
void KyMath::MakePlane(std::vector<Eigen::Vector3d>& arrPts,
                       Eigen::Vector3d& ptOrg,
                       Eigen::Vector3d& vNor) 
{
    KyPrimFit ky_fp;
    ky_fp.MakePlane(arrPts, ptOrg, vNor);
}

double KyMath::DistPtToPlane(Eigen::Vector3d ptOrgPlane,
                             Eigen::Vector3d vNorPlane,
                             Eigen::Vector3d ptThe,
                             Eigen::Vector3d& intP) {
    double d = ptOrgPlane.x() * vNorPlane.x() + ptOrgPlane.y() * vNorPlane.y() +
               ptOrgPlane.z() * vNorPlane.z();

    return DistPtToPlane(vNorPlane.x(), vNorPlane.y(), vNorPlane.z(), d, vNorPlane, ptThe, intP);
}

void KyMath::GetBoundingBox(std::vector<Eigen::Vector3d>& pts, V3D& box)
{
    for (size_t i = 0; i < pts.size(); i++)
    {
        box.Update(pts[i]);
    }
}

Eigen::Vector3d KyMath::GetCenter(std::vector<Eigen::Vector3d>& pts)
{
    Eigen::Vector3d cen;

    if (pts.size() == 0)
        return cen;

    for (size_t i = 0; i < pts.size(); i++)
    {
        cen += pts[i];
    }

    return cen / double(pts.size());
    return Eigen::Vector3d();
}

double KyMath::DistPtToPlane(double a,
                             double b,
                             double c,
                             double d,
                             Eigen::Vector3d vec,
                             Eigen::Vector3d theP,
                             Eigen::Vector3d& intP) {
    double u, v, w;

    u = vec.x();
    v = vec.y();
    w = vec.z();

    intP = {-9.999E10, -9.999E10, -9.999E10};
    //���� ������ ������ ���� �������� ����� �� ����.

    if (fabs(a * u + b * v + c * w) < 1e-10) return 0.;

    // find t
    double t = ((d - (a * theP.x() + b * theP.y() + c * theP.z())) /
                (a * u + b * v + c * w));

    intP.x() = u * t + theP.x();
    intP.y() = v * t + theP.y();
    intP.z() = w * t + theP.z();

    return Dist(intP, theP);
}

void KyMath::DeleteSamePts(std::vector<Eigen::Vector3d>& pts, double vTol)
{
    Eigen::Vector3d pt1, pt2;
    for (int i = 0; i < (int)pts.size() - 1; i++)
    {
        pt1 = pts[i];
        for (int j = i + 1; j < (int)pts.size(); j++)
        {
            pt2 = pts[j];
            //if (pt1.Dist(pt2) < vTol)
            if (KyMath::Dist(pt1, pt2) < vTol)
            {
                pts.erase(pts.begin() + j);
                j--;
            }
        }
    }
}