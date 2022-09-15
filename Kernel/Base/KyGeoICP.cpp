#include "KyGeoICP.h"
#include "KyMUtil.h"
#include "Kyprimfit.h"
//#include <KyBase/KyGeoICP.h>
//#include <KyBase/KyMUtil.h>
//#include <KyBase/KyPrimFit.h>


void KyGeoICP::SetSource(const std::vector<Eigen::Vector3d>& source)
{
	m_Source = source;
}

void KyGeoICP::SetTarget(const std::vector<Eigen::Vector3d>& target)
{
	m_Target = target;
}

void KyGeoICP::Set(const std::vector<Eigen::Vector3d>& source, const std::vector<Eigen::Vector3d>& target)
{
	SetSource(source);
	SetTarget(target);
}

void KyGeoICP::PreProcess()
{
/*	double theta = 0.;
	double beta = 0.;
	double gamma = 0.;

	KyPrincipalAxes sAxis(m_Source);
	KyPrincipalAxes tAxis(m_Target);

	sAxis.GetZAxis().Angle(tAxis.GetZAxis(), theta, beta, gamma);
	
	KyTMatrix tm;
	tm.RotateX(sAxis.GetCenter(), theta);
	tm.RotateY(sAxis.GetCenter(), beta);
	tm.RotateZ(sAxis.GetCenter(), gamma);
	tM.Translate(m_Target.GetCenteer()-sAxis.GetCenter());
	TransformSource(tm);
	m_Source.EstiMB();
*/

}

// classical iterative closest points algorithm
bool KyGeoICP::GetRigidTM(Eigen::Matrix4d& tm, const int& iter, const double& tol, bool errorMetric,bool bMappingFixX, bool bMappingFixY,bool bMappingFixZ)
{
	double error0 = EstiError(errorMetric); // previous error
	double error1 = 0; // current error
	bool improved = false;
	//m_RigidTM.MakeIdentity();
	m_RigidTM.Identity();
	for (int i = 0; i < iter; i++) {
		//tm.MakeIdentity();
		tm.Identity();
		tm = CalcRigidTM(bMappingFixX,bMappingFixY,bMappingFixZ);	
		TransformSource(tm);
		error1 = EstiError(errorMetric);
		if (error0 < error1) {// 
			if (i > 0)
				improved = true;
			//Eigen::Matrix4d inverseM = tm.InverseM();
			Eigen::Matrix4d inverseM = tm.inverse();
			TransformSource(inverseM);
			break;
		}

		m_RigidTM = tm*m_RigidTM;
		
		//if (KyIsZero(error1, tol)) {
		if (fabs(error1) < tol){
		improved = true;
			error0 = error1;
			break;
		}
		//if (KyIsLE(error1, error0, tol)) { // not improved
		if( (error1 - error0) <=tol ){
			improved = true;
		}
		error0 = error1;
	}

	tm = m_RigidTM;
	m_Error = error0;
	
	return improved;
}

// Refer vtkLanderMarkTransform
Eigen::Matrix4d KyGeoICP::CalcRigidTM(bool bMappingFixX/* = false*/,	bool bMappingFixY /*= false*/,	bool bMappingFixZ /*= false*/) const
{ 
	Eigen::Matrix4d tm;
	//tm.MakeIdentity();
	tm.setIdentity();
	int count = (int)m_Source.size();
	if(count == 0)
		return tm;

	Eigen::Vector3d sCp, tCp;
	//sCp = tCp = KY_ORIGIN;
	sCp = tCp = { 0,0,0 };
	int i;
	for (i = 0; i < count; i++) {
		sCp += m_Source[i];
		tCp += m_Target[i];
	}
	
	sCp /= count;
	tCp /= count;

	// -- if only one point, stop right here
	if (count == 1)	{
		//tm.SetT(tCp - sCp);
		tm.block<3, 1>(0, 3) = (tCp - sCp);
		return tm;
	}

	// -- build the 3x3 matrix M --
	double M[3][3], KY[3][3];
	unsigned int j;
	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			KY[i][j] = M[i][j] = 0.;
		}
	}

	//KyStdArray<Eigen::Vector3d>::const_iterator sit;
	//KyStdArray<Eigen::Vector3d>::const_iterator tit;
	std::vector<Eigen::Vector3d>::const_iterator sit;
	std::vector<Eigen::Vector3d>::const_iterator tit;
	Eigen::Vector3d a, b;
	for (sit = m_Source.begin(), tit = m_Target.begin(); 
		 sit != m_Source.end() && tit != m_Target.end(); 
		 sit++, tit++) {		
		// get the origin-centred point (a) in the source set
		a = (*sit); 
		a -= sCp;
		// get the origin-centred point (b) in the target set
		b = (*tit);
		b -= tCp;

		//���� ���� : ������ ����
		if(bMappingFixX)
			a.x() = b.x() = 0.;
		if(bMappingFixY)
			a.y() = b.y() = 0.;
		if(bMappingFixZ)
			a.z() = b.z() = 0.;

		// accumulate the products a*T(b) into the matrix M
		for (i = 0; i < 3; i++) {
			for (j = 0; j < 3; j++) {
				M[i][j] += a[i]*b[j];
			}
		}
	}

	// -- build the 4x4 matrix N --
	double  Ndata[4][4];
	double* N[4];
	for (i = 0; i < 4; i++) {
		N[i] = Ndata[i];
		N[i][0]=N[i][1]=N[i][2]=N[i][3]=0.0; // fill N with zeros
	 }
	// on-diagonal elements
	N[0][0] = M[0][0]+M[1][1]+M[2][2];
	N[1][1] = M[0][0]-M[1][1]-M[2][2];
	N[2][2] = -M[0][0]+M[1][1]-M[2][2];
	N[3][3] = -M[0][0]-M[1][1]+M[2][2];
	// off-diagonal elements
	N[0][1] = N[1][0] = M[1][2]-M[2][1];
	N[0][2] = N[2][0] = M[2][0]-M[0][2];
	N[0][3] = N[3][0] = M[0][1]-M[1][0];

	N[1][2] = N[2][1] = M[0][1]+M[1][0];
	N[1][3] = N[3][1] = M[2][0]+M[0][2];
	N[2][3] = N[3][2] = M[1][2]+M[2][1];

	// -- eigen-decompose N (is symmetric) --
	double eigenvectorData[4][4];
	double *eigenvectors[4],eigenvalues[4];

	eigenvectors[0] = eigenvectorData[0];
	eigenvectors[1] = eigenvectorData[1];
	eigenvectors[2] = eigenvectorData[2];
	eigenvectors[3] = eigenvectorData[3];

	KyMUtil::Jacobi(N,4,eigenvalues,eigenvectors);

	// the eigenvector with the largest eigenvalue is the quaternion we want
	// (they are sorted in decreasing order for us by JacobiN)
	double w, x, y, z;

	// first: if points are collinear, choose the quaternion that 
	// results in the smallest rotation.
	if (eigenvalues[0] == eigenvalues[1] || count == 2) {
		Eigen::Vector3d s0, t0, s1, t1;
		s0 = m_Source[0];
		s1 = m_Source[1];
		t0 = m_Target[0];
		t1 = m_Target[1];
     
		Eigen::Vector3d ds = s1 - s0;
		Eigen::Vector3d dt = t1 - t0;
		//double rs = ds.Size();
		//double rt = dt.Size();
		double rs = KyMath::Dist(ds, { 0,0,0 });
		double rt = KyMath::Dist(dt, { 0,0,0 });

		ds.normalized();
		dt.normalized();


		  // take dot & cross product
		w = ds[0]*dt[0] + ds[1]*dt[1] + ds[2]*dt[2];
		x = ds[1]*dt[2] - ds[2]*dt[1];
		y = ds[2]*dt[0] - ds[0]*dt[2];
		z = ds[0]*dt[1] - ds[1]*dt[0];

	    double r = sqrt(x*x + y*y + z*z);
        double theta = atan2(r,w);

		// construct quaternion
		w = cos(theta/2);
		//if (!KyIsZero(r)) {
		if(!(fabs(r)<1.0E-6)){
			r = sin(theta/2)/r;
			x = x*r;
			y = y*r;
			z = z*r;
		} else {// rotation by 180 degrees: special case			
			// rotate around a vector perpendicular to ds
			//vtkMath::Perpendiculars(ds,dt,0,0);
			r = sin(theta/2);
			x = dt.x()*r;
			y = dt.y()*r;
			z = dt.z()*r;
		}
	} else { // points are not collinear
		w = eigenvectors[0][0];
		x = eigenvectors[1][0];
		y = eigenvectors[2][0];
		z = eigenvectors[3][0];
	}
	// convert quaternion to a rotation matrix
	double ww = w*w;
	double wx = w*x;
	double wy = w*y;
	double wz = w*z;

	double xx = x*x;
	double yy = y*y;
	double zz = z*z;

	double xy = x*y;
	double xz = x*z;
	double yz = y*z;

	Eigen::Vector3d rx, ry, rz;
	//rx.Set(ww + xx - yy - zz, 2.0 * (wz + xy), 2.0 * (-wy + xz));
    //ry.Set(2.0*(-wz + xy), ww - xx + yy - zz, 2.0*(wx + yz));
	//rz.Set(2.0*(wy + xz), 2.0*(-wx + yz), ww - xx - yy + zz);
	rx = { ww + xx - yy - zz, 2.0 * (wz + xy), 2.0 * (-wy + xz) };
	ry = { 2.0 * (-wz + xy), ww - xx + yy - zz, 2.0 * (wx + yz) };
	rz = { 2.0 * (wy + xz), 2.0 * (-wx + yz), ww - xx - yy + zz };
	
	//tm.SetX(rx);
	//tm.SetY(ry);
	//tm.SetZ(rz);
	tm.setIdentity();
	tm.block<3, 1>(0, 0) = rx;
	tm.block<3, 1>(0, 1) = ry;
	tm.block<3, 1>(0, 2) = rz;
	// the translation is given by the difference in the transformed source
	// centroid and the target centroid
	//tm.SetT(tCp - tm*sCp);
	tm.block<3, 1>(0, 3) = tCp - tm * sCp;

	return tm;
}


// Find unit vectors which is perpendicular to this on and to
// each other.
void KyGeoICP::Perpendiculars(const Eigen::Vector3d& v0,
	Eigen::Vector3d& v1, Eigen::Vector3d& v2, double& theta) const
{
	int dx, dy, dz;
	double x2 = v0.x() * v0.x();
	double y2 = v0.y() * v0.y();
	double z2 = v0.z() * v0.z();
	double r = sqrt(x2 + y2 + z2);

	// transpose the vector to avoid divide-by-zero error
	if (x2 > y2 && x2 > z2) {
		dx = 0; dy = 1; dz = 2;
	} else if (y2 > z2) {
		dx = 1; dy = 2; dz = 0;
	} else {
		dx = 2; dy = 0; dz = 1;
	}

	double a = v0[dx]/r;
	double b = v0[dy]/r;
	double c = v0[dz]/r;

	double tmp = sqrt(a*a+c*c);

	//if ( KyIsNotZero(theta) ) {
	if ( fabs(theta) >= 1.0E-6) {
		double sintheta = sin(theta);
		double costheta = cos(theta);
		if (&v1) {
			v1[dx] = (c*costheta - a*b*sintheta)/tmp;
			v1[dy] = sintheta*tmp;
			v1[dz] = (-a*costheta - b*c*sintheta)/tmp;
		}

		if (&v2) {
			v2[dx] = (-c*sintheta - a*b*costheta)/tmp;
			v2[dy] = costheta*tmp;
			v2[dz] = (a*sintheta - b*c*costheta)/tmp;
		}
	} else {
		if (&v1) {
			v1[dx] = c/tmp;
			v1[dy] = 0;
			v1[dz] = -a/tmp;
		}

		if (&v2) {
			v2[dx] = -a*b/tmp;
			v2[dy] = tmp;
			v2[dz] = -b*c/tmp;
		}
	}      

}

bool KyGeoICP::TransformSource(const Eigen::Matrix4d& tm)
{	
	//KyStdArray<Eigen::Vector3d>::iterator sit;
	std::vector<Eigen::Vector3d>::iterator sit;
	for (sit = m_Source.begin(); sit != m_Source.end(); sit++) {		
		(*sit) =  tm*(*sit);
	}
	return true;
}


double KyGeoICP::EstiError(bool mode)
{
	//KyStdArray<Eigen::Vector3d>::const_iterator sit;
	//KyStdArray<Eigen::Vector3d>::const_iterator tit;
	std::vector<Eigen::Vector3d>::const_iterator sit;
	std::vector<Eigen::Vector3d>::const_iterator tit;

	double totalDist = 0.;
	double meanDist = 0.;
	for (sit = m_Source.begin(), tit = m_Target.begin(); 
		 sit != m_Source.end() && tit != m_Target.end(); 
		 sit++, tit++) {		
		if (mode) // RMS
			//totalDist += (*sit).Dist(*tit) * (*sit).Dist(*tit); 
			totalDist += KyMath::Dist((*sit), (*tit)) * KyMath::Dist((*sit), (*tit));
		else
			totalDist += sqrt(KyMath::Dist((*sit), (*tit)));
	}

	if (mode) // RMS
		meanDist = sqrt(totalDist/m_Source.size());
	else
		meanDist = totalDist/m_Source.size();
   
	return meanDist;
}
	
KyGeoICP::KyGeoICP()
{
	;
}
	
KyGeoICP::~KyGeoICP()
{
	
	;
}
