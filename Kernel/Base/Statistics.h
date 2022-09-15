#pragma once
#include <vector>

class Statistics {
public:
    Statistics(void);
    void CalcStatistics(std::vector<double>& arrV, double& vStd, double& vAve);
    void CalcStatistics(std::vector<double>& arrV, double& vStd, double& vAve, double& minDis, double& maxDis, double& vMinMax);
};