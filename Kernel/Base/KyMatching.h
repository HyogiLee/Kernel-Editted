#pragma once
#include "open3d/Open3D.h"
//#include <KyBase/KyBaseExportDll.h>

// Degree of constraints
#define LVL_LOWER	2
#define LVL_MIDDLE	5
#define LVL_UPPER	10

#define NUM_OF_DATA  10000



class KyAutoMatch
{
public:
	KyAutoMatch();
	~KyAutoMatch();

	void Init();

	void with_Constraint_doit();
	void with_Filter_doit(bool bRotateZ=true);

	void AllocBaseData(int sz);
	void AllocMeasures(int sz);
	void AllocConstraints(int sz);

	void SetBaseData(int i, int idx, double p[3]);
	void SetMeasures(int i, int idx, double p[3]);
	void SetConstraint(int i, int idx, int level);
	void SetDeltaFiltered(double d, bool b) {m_delta=d; m_bFiltered=b;}

	// -- result
	void GetBaseData(int i, int &idx, double p[3]) const;
	void GetMeasures(int i, int &idx, double p[3]) const;
	double** GetTransformed() const { return m_oTransformed; }
	int* GetIdxM() const { return m_mIndex; }
	int* GetIdxO() const { return m_oIndex; }
	int  GetNFiltered() const { return m_nfiltered; }

private:
	double m_delta;
	bool m_bFiltered;
	int m_nfiltered;

	int m_nBase, m_nMeasured, m_nConstraints;
	int *m_bIndex, *m_mIndex, *m_oIndex;
	int *m_cIndex, *m_cLevel;
	double **m_bData,**m_mData,**m_oTransformed;
};


 
struct Closest_Pair {
  double First[3];
  double Base[3];
  double dist;
  struct Closest_Pair *next;
};

void	with_Constraint(int nNormal,int *normalIndex,double **normalXYZ,
						 int nmeasured,int *measuredIndex,double **measuredXYZ,
						 int nconstraint,int *constraintIndex,int *constraintLevel,
						 int *nfiltered,double **transformedXYZ,int *filteringindex,
						 double *distData,double delta, bool filteredYN);
void with_Filtering(int nNormal,int *normalIndex,double **normalXYZ,
					int nmeasured,int *measuredIndex,double **measuredXYZ,
					int *nfiltered,double **transformedXYZ,int *filteringindex,
					double *distData,double delta, bool filteredYN, bool bOnlyRotateZ =false);

void with_Filtering_native(int nNormal,int *normalIndex,double **normalXYZ,
					int nmeasured,int *measuredIndex,double **measuredXYZ,
					int *nfiltered,double **transformedXYZ,int *filteringindex,
					double *distData,double delta, bool filteredYN);

double find_motion(double **bData,int nB,double **mRotatedData,int nM);
double find_motion_final(double **bData,int nB,double **mRotatedData,int nM,double cbx,double cby,double cbz);
double find_motion_without_link(double **bData,int nB,double **mRotatedData,int *nM,double *mDist, int *mCorres);
void change_InitialPosition(double x,double y,double z,int n,double **mData,double **mR);
void createTable(double **bD,int nB,double **mR,int nM,double *mD,int *mC, double delta, bool useFileter);
void updateCorres(double **bD,int nB,double **mR,int nM,double *mD,int *mC, double delta, bool useFileter);
void moveBackBaseData(double **bD,int nB,double x,double y,double z);
void deleteFilteredData(int *measuredIndex,int *nfiltered,int *filteringIndex,double **oTransformed,double *dist,bool filteredYN);
void sort_dist(double **mR,int nM,double *mD,int *mC,int left, int right);
void swap(double **mR,double *mD,int i, int j);

double Assign_Check_Final(double set1[NUM_OF_DATA][3],struct Closest_Pair *Head);
struct Closest_Pair *Make_List_Of_Pair(int n_set);
void Initialize(double Mat[4][4]);
void eigsrt(double d[5], double v[5][5], int n);
void jacobi(double a[5][5], int n, double d[5], double v[5][5], int *nrot);
void Find_Closest_Pair_List(double set1[NUM_OF_DATA][3],int n_set1,double base[NUM_OF_DATA][3],int n_base,double MuP[3],double MuX[3],struct Closest_Pair* Head);
void Update(double Rot[3][3],double p[3],struct Closest_Pair *Head);
void Calc_Q(double Q[4][4],double Covariance[3][3]);
void Calc_Covariance(double Covariance[3][3],struct Closest_Pair *Head,double MuP[3],double MuX[3]);
void Rotation_Matrix(double Rot[3][3], double quaternion[4]);
void Assign_Points(struct Closest_Pair* Current,double Points[3],double Min_Pair[3],double Min_Dist);
void Make_Translation_Vector(double qt[3],double Rot[3][3],double MuP[3],double MuX[3]);
void Out_Result(struct Closest_Pair *Head,int n_set1,double cbx,double cby,double cbz,double **mR);
void Out_Result_without_link(struct Closest_Pair *Head,int n_set1,double **mR,double *mD);
void Out_Trans_Rot(double Rot[3][3],double qt[3]);
void Free_List(struct Closest_Pair *Head);

int findMax(double *dist,int nf);
double Find_Min_Dist_Pair(double Points[3],double Min_Pair[3],int n_base,double base[NUM_OF_DATA][3]);
void Assign_Points(struct Closest_Pair* Current,double Points[3],double Min_Pair[3],double Min_Dist);
