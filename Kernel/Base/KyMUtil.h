#pragma once
#include "open3d/Open3D.h"

    
/** Math Utility
*/
class KyMUtil {
public:
	static int Jacobi(double **a, int n, double *w, double **v);
};
