#pragma once

#include "open3d/Open3D.h"
#include <minmax.h>
#include "V3D.h"
#define KySQR(x) ((x) * (x))
#define _USE_MATH_DEFINES
#include <cmath>
#include <math.h>
class KyMath
{
public:
    KyMath(void){};
    ~KyMath(void){};

	public:
        double Sq(double x) { return x * x; }


        static double SqDist(const Eigen::Vector3d& v1,
                      const Eigen::Vector3d& v2) {
            return (KySQR(v1.x() - v2.x()) + KySQR(v1.y() - v2.y()) +
                    KySQR(v1.z() - v2.z()));
        }

        double SqDistXY(const Eigen::Vector3d& v1,
                        const Eigen::Vector3d& v2) const {
            return (KySQR(v1.x() - v2.x()) + KySQR(v1.y() - v2.y()));
        }

        double SqDistYZ(const Eigen::Vector3d& v1,
                        const Eigen::Vector3d& v2) const {
            return (KySQR(v1.y() - v2.y()) + KySQR(v1.z() - v2.z()));
        }

        double SqDistZX(const Eigen::Vector3d& v1,
                        const Eigen::Vector3d& v2) const {
            return (KySQR(v1.z() - v2.z()) + KySQR(v1.x() - v2.x()));
        }

        static double Dist(const Eigen::Vector3d& v1,
                    const Eigen::Vector3d& v2) {
            return sqrt(SqDist(v1, v2));
        }
        double DistXY(const Eigen::Vector3d& v1,
                      const Eigen::Vector3d& v2) const {
            return sqrt(SqDistXY(v1, v2));
        }
        double DistYZ(const Eigen::Vector3d& v1,
                      const Eigen::Vector3d& v2) const {
            return sqrt(SqDistYZ(v1, v2));
        }
        double DistZX(const Eigen::Vector3d& v1,
                      const Eigen::Vector3d& v2) const {
            return sqrt(SqDistZX(v1, v2));
        }

        double SqSize(const Eigen::Vector3d& v) const {
            return (KySQR(v.x()) + KySQR(v.y()) + KySQR(v.z()));
        }
        double Size(const Eigen::Vector3d& v) const { return sqrt(SqSize(v)); }
    
    public:
        void MakePlane(std::vector<Eigen::Vector3d>& arrPts,
                          Eigen::Vector3d& ptOrg,
                          Eigen::Vector3d& vNor);

        double DistPtToPlane(double a,
                         double b,
                         double c,
                         double d,
                         Eigen::Vector3d vec,
                         Eigen::Vector3d theP,
                         Eigen::Vector3d& intP);

        double DistPtToPlane(Eigen::Vector3d ptOrgPlane,
                         Eigen::Vector3d vNorPlane,
                         Eigen::Vector3d ptThe,
                         Eigen::Vector3d& intP);
        
        static void GetBoundingBox(std::vector<Eigen::Vector3d>& pts, V3D& box);
        static Eigen::Vector3d GetCenter(std::vector<Eigen::Vector3d>& pts);
        //static double Radian(double degree) { return degree * M_PI / 180.0; }
        static double Radian(double degree) { return degree * 3.14159265358979323846 / 180.0; }
        static void DeleteSamePts(std::vector<Eigen::Vector3d>& pts, double vTol = 1.0E-1);
};

class KyP4DsComparer
{
public:
    bool operator()(std::pair<double, std::vector<Eigen::Vector3d>>& rhs, std::pair<double, std::vector<Eigen::Vector3d>>& lhs)
	{
		return (rhs.first < lhs.first);
	}
};