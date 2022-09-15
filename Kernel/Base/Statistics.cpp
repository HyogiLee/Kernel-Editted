#include "Statistics.h"

Statistics::Statistics(void) {}

void Statistics::CalcStatistics(std::vector<double>& arrV,
                                  double& vStd,
                                  double& vAve) {
    double minDis, maxDis, vMinMax;

    CalcStatistics(arrV, vStd, vAve, minDis, maxDis, vMinMax);
}

void Statistics::CalcStatistics(std::vector<double>& arrV,
                                  double& vStd,
                                  double& vAve,
                                  double& minDis,
                                  double& maxDis,
                                  double& vMinMax) {
    vAve = 0.;
    vStd = 0.;
    vMinMax = 0.;
    minDis = 0.;
    maxDis = 0.;

    if (arrV.size() < 1) return;

    for (size_t i = 0; i < arrV.size(); i++) vAve += arrV[i];

    vAve /= arrV.size();

    if (arrV.size() < 2) return;

    for (size_t i = 0; i < arrV.size(); i++)
        vStd += (vAve - arrV[i]) * (vAve - arrV[i]);
    vStd = sqrt(vStd / (arrV.size() - 1));

    maxDis = -1e10;
    minDis = 1e10;

    for (size_t i = 0; i < arrV.size(); i++) {
        if (arrV[i] > maxDis) maxDis = arrV[i];
        if (arrV[i] < minDis) minDis = arrV[i];
    }

    vMinMax = maxDis - minDis;
}
