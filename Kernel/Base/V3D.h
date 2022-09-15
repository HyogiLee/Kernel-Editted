#pragma once
#include "open3d/Open3D.h"

class V3D:Eigen::Vector3d 
{

public:
    V3D() {
        Min = {10000, 10000, 10000};
        Max = {-10000, -10000, -10000};
    };
    Eigen::Vector3d Min;
    Eigen::Vector3d Max;
    Eigen::Vector3d Center(void) const { return (Max + Min) / double(2); }
    double Cx() const { return (Max.x() + Min.x()) / double(2); }
    double Cy() const { return (Max.y() + Min.y()) / double(2); }
    double Cz() const { return (Max.z() + Min.z()) / double(2); }
    double Dx() const { return Max.x() - Min.x(); }
    double Dy() const { return Max.y() - Min.y(); }
    double Dz() const { return Max.z() - Min.z(); }

    double MaxWidth() const {
        double w;
        if (Dx() >= Dy() && Dx() >= Dz())
            w = Dx();
        else if (Dy() >= Dz())
            w = Dy();
        else
            w = Dz();
        return w;
    }

    double MinWidth() const {
        double w;
        if (Dx() <= Dy() && Dx() <= Dz())
            w = Dx();
        else if (Dy() <= Dz())
            w = Dy();
        else
            w = Dz();
        return w;
    }

    void Update(Eigen::Vector3d vec) { 
        Update(vec.x(), vec.y(), vec.z());
    }

    void Update(double x, double y, double z) {

        if (x < Min.x())
            Min.x() = x;
        else if (x > Max.x())
            Max.x() = x;

        if (y < Min.y())
            Min.y() = y;
        else if (y > Max.y())
            Max.y() = y;

        if (z < Min.z())
            Min.z() = z;
        else if (z > Max.z())
            Max.z() = z;
    }

    void Union(const V3D& mb) 
    {
        Min.x() = min(Min.x(), mb.Min.x());
        Min.y() = min(Min.y(), mb.Min.y());
        Min.z() = min(Min.z(), mb.Min.z());
        Max.x() = min(Max.x(), mb.Max.x());
        Max.y() = min(Max.y(), mb.Max.y());
        Max.z() = min(Max.z(), mb.Max.z());
   }

    void Expand(double d) {
       Min.x() -= d;
       Min.y() -= d;
       Min.z() -= d;
       Max.x() += d;
       Max.y() += d;
       Max.z() += d;
   }
   //
   void ExpandX(double d) {
       Min.x() -= d;
       Max.x() += d;
   }
   void ExpandY(double d) {
       Min.y() -= d;
       Max.y() += d;
   }
   void ExpandZ(double d) {
       Min.z() -= d;
       Max.z() += d;
   }
};