#include "Kyprimfit.h"
#include "KyMUtil.h"

bool KyPrimFit::MakePlane(const std::vector<Eigen::Vector3d>& arrPts,
                          Eigen::Vector3d& org,
                          Eigen::Vector3d& norm) {
    KyPrincipalAxes prinsAxis(arrPts);
	org = prinsAxis.GetCenter();
	
	norm = prinsAxis.GetZAxis();
    norm.normalize(); 

	return true;	
}

void KyPrincipalAxes::Set(const std::vector<Eigen::Vector3d> pts)
{
	int i, j;
	
	// reset the intermediate arrays
    m_Center = {0, 0, 0};

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
			m_EigenvalueProblem[i][j] = 0.;

	for (i = 0; i < 3; i++)
		for (j = 0; j < 3; j++)
		    m_EigenvalueProblemDiag[i][j] = 0.;
 
    // compute the center
	for (i = 0; i < (int)pts.size(); i++) {
		m_Center += pts[i];
	}
 	m_Center /= pts.size();

   // create the eigenvalue-problem
   // using the symmetry of the result matrix
	for (i = 0; i < 3; i++)
		for (j = i; j < 3; j++)
			m_EigenvalueProblem[i][j] = -m_Center[i]*pts.size()*m_Center[j];
 
	Eigen::Vector3d x; 
	for (int m = 0; m < (int)pts.size(); m++) {
		x = pts[m];
		for(i = 0; i < 3;i++)
			for(j = i; j < 3;j++)
				m_EigenvalueProblemDiag[i][j] +=  x[i]*x[j];
	}

	for (i = 0; i < 3; i++)
		for (j = i; j <3; j++)
			m_EigenvalueProblem[i][j] += m_EigenvalueProblemDiag[i][j];

	for (i = 0; i < 3; i++)
		for (j = 0; j < i; j++)
			m_EigenvalueProblem[i][j] = m_EigenvalueProblem[j][i];

	KyMUtil::Jacobi(m_EigenvalueProblem, 3, m_Eigenvalues, m_Eigenvectors);
 
	// update Axes
	for (i = 0; i < 3; i++) {
		m_XAxis[i] = m_Eigenvectors[i][0];
		m_YAxis[i] = m_Eigenvectors[i][1];
		m_ZAxis[i] = m_Eigenvectors[i][2];
	}
}

void KyPrincipalAxes::Init0()
{
    m_Center = {0, 0, 0};
    m_XAxis = {1, 0, 0};
    m_YAxis = {0, 1, 0};
    m_ZAxis = {0, 0, 1};

	m_EigenvalueProblem = new double*[3];
	m_EigenvalueProblemDiag = new double*[3];
	m_Eigenvectors = new double*[3];
	for (int i = 0; i < 3; i++) {
		m_EigenvalueProblem[i] = new double[3];
		m_EigenvalueProblemDiag[i] = new double[3];
		m_Eigenvectors[i] = new double[3];
	}
	m_Eigenvalues = new double[3];
}

void KyPrincipalAxes::Clear()
{
	for (int i = 0; i < 3; i++) {
	   delete [] m_EigenvalueProblem[i];
	   delete [] m_EigenvalueProblemDiag[i];
	   delete [] m_Eigenvectors[i];
	}
	delete [] m_EigenvalueProblem;
	delete [] m_EigenvalueProblemDiag;
	delete [] m_Eigenvectors;
	delete [] m_Eigenvalues;
}


KyPrincipalAxes::KyPrincipalAxes(const std::vector<Eigen::Vector3d> pts)
{
	Init0();
	Set(pts);
}

KyPrincipalAxes::KyPrincipalAxes()
{
	Init0();
}

KyPrincipalAxes::~KyPrincipalAxes()
{
	Clear();
}
