#include <math.h>
#include "KyMatching.h"
//#include <KyBase/KyMatching.h>

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif


KyAutoMatch::KyAutoMatch()
{
	m_nBase = 0 ;
	m_nMeasured = 0 ;
	m_nConstraints = 0;
	m_nfiltered = 0;

	Init();
}


KyAutoMatch::~KyAutoMatch()
{
	Init();
}

void KyAutoMatch::Init()
{
	int i;
	if ( m_nBase > 0 ) {
		for ( i = 0 ; i < m_nBase ; i++ )
			free(m_bData[i]);
		free(m_bData);
		free(m_bIndex);
		m_nBase = 0;
	}
	if ( m_nMeasured > 0 ) {
		for ( i = 0 ; i < m_nMeasured ; i++ ) {
			free(m_mData[i]);
			free(m_oTransformed[i]);
		}
		free(m_mData);
		free(m_oTransformed);
		free(m_mIndex);
		free(m_oIndex);
		m_nMeasured = 0;
	}
	if ( m_nConstraints > 0 ) {
		free(m_cIndex);
		free(m_cLevel);
		m_nConstraints = 0;
	}

	m_bIndex = NULL;
	m_mIndex = NULL;
	m_oIndex = NULL;
	m_cIndex = NULL;
	m_cLevel = NULL;

	m_bData = NULL;
	m_mData = NULL;
	m_oTransformed = NULL;
}


void KyAutoMatch::AllocBaseData(int sz)
{
	if ( sz < 1 ) return;

	m_nBase = sz;

	m_bData = (double**) malloc(sizeof(double*)*(sz));
	m_bIndex= (int*) malloc(sizeof(int)*(sz));
	for(int i=0;i<sz;i++) {
		m_bData[i]=(double*) malloc(sizeof(double)*3);
	}
}

void KyAutoMatch::AllocMeasures(int sz)
{
	if ( sz < 1 ) return;

	m_nMeasured = sz;

	m_mData = (double**) malloc(sizeof(double*)*(sz));
	m_oTransformed = (double **) malloc(sizeof(double*)*(sz));
	m_mIndex= (int*) malloc(sizeof(int)*(sz));
	m_oIndex= (int*) malloc(sizeof(int)*(sz));
	for(int i=0;i<sz;i++) {
		m_mData[i]=(double*) malloc(sizeof(double)*3);
		m_oTransformed[i] = (double*) malloc(sizeof(double)*3);
	}
}

void KyAutoMatch::AllocConstraints(int sz)
{
	if ( sz < 1 ) return;

	m_nConstraints = sz;

	m_cIndex= (int*) malloc(sizeof(int)*(sz));
	m_cLevel= (int*) malloc(sizeof(int)*(sz));
}

void KyAutoMatch::GetBaseData(int i, int &idx, double p[3]) const
{ 
	idx = m_bIndex[i]; 
	p[0] = m_bData[i][0]; 
	p[1] = m_bData[i][1]; 
	p[2] = m_bData[i][2]; 
}

void KyAutoMatch::GetMeasures(int i, int &idx, double p[3]) const
{ 
	idx = m_mIndex[i]; 
	p[0] = m_mData[i][0]; 
	p[1] = m_mData[i][1]; 
	p[2] = m_mData[i][2]; 
}

void KyAutoMatch::SetBaseData(int i, int idx, double p[3])
{
	if ( i < 0 || i >= m_nBase ) return;

	for ( int j = 0 ; j < 3 ; j++ )
		m_bData[i][j] = p[j];
	m_bIndex[i] = idx;
}

void KyAutoMatch::SetMeasures(int i, int idx, double p[3])
{
	if ( i < 0 || i >= m_nMeasured ) return;

	for ( int j = 0 ; j < 3 ; j++ )
		m_mData[i][j] = p[j];
	m_mIndex[i] = idx;
}

void KyAutoMatch::SetConstraint(int i, int idx, int level)
{
	int j;
	if ( i < 0 || i >= m_nConstraints ) return;

	for ( j = 0 ; j < m_nBase ; j++ )
		if ( m_bIndex[j] == idx ) break;
	m_cIndex[i] = j;
	m_cLevel[i] = level;
}

void KyAutoMatch::with_Constraint_doit()
{
	// Function call for constrained matching

	// nBase: ���� ����Ÿ�� ����
	// bIndex: ���� ����Ÿ�� �ε���
	// bData: ���� ����Ÿ ����Ʈ
	// nMeasured: ���� ����Ÿ�� ����
	// mIndex: ���� ����Ÿ�� �ε���
	// mData: ���� ����Ÿ ����Ʈ
	// nconstraint: constraint�� ����
	// cIndex: ���������� �ο��Ǵ� ������ �ε���
	// cLevel: ���� ������ ����
	// nfiltered: ���յ� ������ ���� (���͸��� �� ������ ���� ����)
	// oTransformed: ���յ� ���� ��. ���ο� ��ȯ ������ �Լ������� ���� �޴´�
	// filteringIndex: �� �������� ��Ī�Ǵ� ���԰��� �ε���
	// dist: �� ���� ���鰣�� �ּ� �Ÿ�
	// filteredYN: ���͸� ����


	int indexMax;
	double subdelta;

	double *dist = (double *)malloc(sizeof(double)*(m_nMeasured));
	double *distBackup = (double *)malloc(sizeof(double)*(m_nMeasured));

/*	with_Constraint(m_nBase,m_bIndex,m_bData,
			m_nMeasured,m_mIndex,m_mData,m_nConstraints,m_cIndex,m_cLevel,
			&m_nfiltered,m_oTransformed,m_oIndex,dist,m_delta,m_bFiltered);*/

	int i,i2;

	if(m_bFiltered) {

		with_Constraint(m_nBase,m_bIndex,m_bData,
			m_nMeasured,m_mIndex,m_mData,m_nConstraints,m_cIndex,m_cLevel,
			&m_nfiltered,m_oTransformed,m_oIndex,dist,m_delta,false);
		
		for(i=0;i<m_nMeasured;i++)
			distBackup[i] = dist[i];

		/* find the max dist */

		indexMax = findMax(dist,m_nfiltered);
		if (indexMax == -1)
			return;

		i2 = indexMax;

		while(dist[indexMax] > m_delta) {
			subdelta = distBackup[i2]-1.0e-6;
			for(i=0;i<m_nMeasured;i++)
				m_mIndex[i] = i;
			with_Constraint(m_nBase,m_bIndex,m_bData,
			m_nMeasured,m_mIndex,m_mData,m_nConstraints,m_cIndex,m_cLevel,
			&m_nfiltered,m_oTransformed,m_oIndex,dist,subdelta,true);
			
			distBackup[i2] = -1.0;

			indexMax = findMax(dist,m_nfiltered);
			if (indexMax == -1)
				break;

			i2 = findMax(distBackup,m_nMeasured);
		}
	} else 
		with_Constraint(m_nBase,m_bIndex,m_bData,
			m_nMeasured,m_mIndex,m_mData,m_nConstraints,m_cIndex,m_cLevel,
			&m_nfiltered,m_oTransformed,m_oIndex,dist,m_delta,m_bFiltered);
		
	// -- Write Result
	


	free(dist);
	
	free(distBackup);	
}


void KyAutoMatch::with_Filter_doit(bool bRotateZ)
{
	int indexMax;
	double subdelta;
	//�ڵ����� with filtering��� �Լ� ȣ��
	// nBase: ���� ����Ÿ�� ����
	// bIndex: ���� ����Ÿ�� �ε���
	// bData: ���� ����Ÿ ����Ʈ
	// nMeasured: ���� ����Ÿ�� ����
	// mIndex: ���� ����Ÿ�� �ε���
	// mData: ���� ����Ÿ ����Ʈ
	// nfiltered: ���յ� ������ ���� (���͸��� �� ������ ���� ����)
	// oTransformed: ���յ� ���� ��
	// filteringIndex: �� �������� ��Ī�Ǵ� ���԰��� �ε���
	// dist: �� ���� ���鰣�� �ּ� �Ÿ�
	// filteredYN: ���͸� ����

	double *dist = (double *)malloc(sizeof(double)*(m_nMeasured));
	double *distBackup = (double *)malloc(sizeof(double)*(m_nMeasured));
	/*with_Filtering(m_nBase,m_bIndex,m_bData,m_nMeasured,m_mIndex,m_mData,
				  &m_nfiltered,m_oTransformed,m_oIndex,dist,m_delta,m_bFiltered);*/
	int i,i2;

	if(m_bFiltered) 
	{

		with_Filtering(m_nBase,m_bIndex,m_bData,m_nMeasured,m_mIndex,m_mData,
					  &m_nfiltered,m_oTransformed,m_oIndex,dist,m_delta,false, bRotateZ);

		for(i=0;i<m_nMeasured;i++)
			distBackup[i] = dist[i];

		/* find the max dist */

		indexMax = findMax(dist,m_nfiltered);
		i2 = indexMax;

		while(dist[indexMax] > m_delta) 
		{

			for(i=0;i<m_nMeasured;i++)
					m_mIndex[i] = i;	

			subdelta = distBackup[i2]-1.0e-6;

			with_Filtering(m_nBase,m_bIndex,m_bData,m_nMeasured,m_mIndex,m_mData,
				  &m_nfiltered,m_oTransformed,m_oIndex,dist,subdelta,true, bRotateZ);
			distBackup[i2] = -1.0;

			indexMax = findMax(dist,m_nfiltered);
			i2 = findMax(distBackup,m_nMeasured);
		}
	}
	else 
		with_Filtering(m_nBase,m_bIndex,m_bData,m_nMeasured,m_mIndex,m_mData,
				  &m_nfiltered,m_oTransformed,m_oIndex,dist,m_delta,m_bFiltered, bRotateZ);

	// -- Write Result
	free(dist);
	free(distBackup);
}





void	with_Constraint(int nBase,int *normalIndex,double **bData,
						 int nMeasured,int *measuredIndex,double **mData,
						 int nconstraint,int *constraintIndex,int *constraintLevel,
						 int *nfiltered,double **transformedXYZ,int *filteringindex,
						 double *distData,double delta, bool filteredYN)
{
	int i;
	double **mRotatedData = (double **) malloc(sizeof(double*)*(nMeasured));
	for(int i=0;i<nMeasured;i++)
		mRotatedData[i]=(double*) malloc(sizeof(double)*3);

	// Compute the centroids
	double cmx,cmy,cmz,cbx,cby,cbz;
	cmx = 0.;
	cmy = 0.;
	cmz = 0.;
	for(i=0;i<nMeasured;i++) {
	cmx += mData[i][0];
	cmy += mData[i][1];
	cmz += mData[i][2];
	}
	cmx /= (double)nMeasured;
	cmy /= (double)nMeasured;
	cmz /= (double)nMeasured;

	cbx = 0.;
	cby = 0.;
	cbz = 0.;
	for(i=0;i<nBase;i++) {
	cbx += bData[i][0];
	cby += bData[i][0];
	cbz += bData[i][0];
	}
	cbx /= (double)nBase;
	cby /= (double)nBase;
	cbz /= (double)nBase;

	// Move both data sets with respect to the origin
	for(i=0;i<nMeasured;i++) {
		mData[i][0] -= cmx;
		mData[i][1] -= cmy;
		mData[i][2] -= cmz;
	}

	for(i=0;i<nBase;i++) {
		bData[i][0] -= cbx;
		bData[i][1] -= cby;
		bData[i][2] -= cbz;
	}


	// Compute the various initial positions for optimization
	double px,py,pz;
	double stx,sty,stz;
	double End_Loop,MinValue;
	double min_px,min_py,min_pz;

	MinValue = 1.0e+16;
	stx = sty = stz = 30.0;
	px = 0.0;
	while(px < 360.0) {
		py = 0.0;
		while(py < 360.0) {
			pz = 0.0;
			while(pz < 360.0) {
				change_InitialPosition(px,py,pz,nMeasured,mData,mRotatedData);
				End_Loop = find_motion(bData,nBase,mRotatedData,nMeasured); 
				if(End_Loop < MinValue) {
					min_px = px;
					min_py = py;
					min_pz = pz;
					MinValue = End_Loop;
				}
				pz += stz;
			}
			py += sty;
		}
		px += stx;
	}


	// Orientation angles for optimum initial position
	change_InitialPosition(min_px,min_py,min_pz,nMeasured,mData,mRotatedData);
	End_Loop = find_motion_final(bData,nBase,mRotatedData,nMeasured,cbx,cby,cbz);

	moveBackBaseData(bData,nBase,cbx,cby,cbz);
	createTable(bData,nBase,mRotatedData,nMeasured,distData,filteringindex,delta,filteredYN);
	find_motion_without_link(bData,nBase,mRotatedData,&nMeasured,distData,filteringindex);

	*nfiltered = nMeasured;

	for(i=0;i<nMeasured;i++) {
		transformedXYZ[i][0] = mRotatedData[i][0];
		transformedXYZ[i][1] = mRotatedData[i][1];
		transformedXYZ[i][2] = mRotatedData[i][2];
	}

	// Adjust based on the weights for each data point
	double newMx,newMy,newMz,newTx,newTy,newTz;
	int totalWeight,weight,indexBase;

	newMx = newMy = newMz = 0.0;
	newTx = newTy = newTz = 0.0;
	totalWeight = 0;
	for(i=0;i<nconstraint;i++) {
		weight = constraintLevel[i];
		indexBase = filteringindex[constraintIndex[i]]-1;

		newMx += weight*bData[indexBase][0];
		newMy += weight*bData[indexBase][1];
		newMz += weight*bData[indexBase][2];
		totalWeight += weight;

		newTx += weight*transformedXYZ[constraintIndex[i]][0];
		newTy += weight*transformedXYZ[constraintIndex[i]][1];
		newTz += weight*transformedXYZ[constraintIndex[i]][2];
	}

	double dtx,dty,dtz;
	dtx = dty = dtz = 0.0;
	if(totalWeight != 0) {
		// Compute the translation vector
		dtx = (newMx - newTx) / (double)totalWeight;
		dty = (newMy - newTy) / (double)totalWeight;
		dtz = (newMz - newTz) / (double)totalWeight;
	}

	// Final points
	for(i=0;i<nMeasured;i++) {
		transformedXYZ[i][0] += dtx;
		transformedXYZ[i][1] += dty;
		transformedXYZ[i][2] += dtz;
	}
	createTable(bData,nBase,transformedXYZ,nMeasured,distData,filteringindex,delta,filteredYN);

	for(i=0;i<nMeasured;i++)
		free(mRotatedData[i]);
	free(mRotatedData);
}

void moveBackBaseData(double **bD,int nB,double x,double y,double z)
{
  int i;
  for(i=0;i<nB;i++) {
    bD[i][0] += x;
    bD[i][1] += y;
    bD[i][2] += z;
  }
}



void createTable(double **bD,int nB,double **mR,int nM,double *mD,int *mC, double delta, bool useFileter)
{
  int i,j;
  double dist;
  double minDist;


  for(i=0;i<nM;i++) {
    minDist = 1.0e+16;
    for(j=0;j<nB;j++) {
      dist = (bD[j][0]-mR[i][0])*(bD[j][0]-mR[i][0])+(bD[j][1]-mR[i][1])*(bD[j][1]-mR[i][1])+(bD[j][2]-mR[i][2])*(bD[j][2]-mR[i][2]);
      dist = sqrt(dist);
      if(dist < minDist) {
	mC[i] = j+1;
	mD[i] = dist;
	minDist = dist;
      }
    }
  }

	//Marking LINK points
	//Fitering flag is set. Filtering is performed.
	if ( !useFileter ) return ;
   
    for(i=0;i<nM;i++) {
      if(mD[i] > delta) mC[i] = -9999;
    }
}


void change_InitialPosition(double x,double y,double z,int n,double **mData,double **mR)
{
  double R[3][3];

  double xc,xs,yc,ys,zc,zs;

  int i;

  //  [6/7/2011 KKY]
  x = x * 3.14159265 / 180.;
  y = y * 3.14159265 / 180.;
  z = z * 3.14159265 / 180.;


  xc = cos(x);xs = sin(x);
  yc = cos(y);ys = sin(y);
  zc = cos(z);zs = sin(z);

  R[0][0] = yc*zc;
  R[0][1] = xs*ys*zc - xc*zs;
  R[0][2] = xc*ys*zc + xs*zs;
  
  R[1][0] = yc * zs;
  R[1][1] = xs*ys*zs + xc*zc;
  R[1][2] = xc*ys*zs - xs * zc;
  
  R[2][0] = -ys;
  R[2][1] = xs*yc;
  R[2][2] = xc*yc;

  for(i=0;i<n;i++) {
    mR[i][0] = R[0][0]*mData[i][0] + R[0][1]*mData[i][1] + R[0][2]*mData[i][2];
    mR[i][1] = R[1][0]*mData[i][0] + R[1][1]*mData[i][1] + R[1][2]*mData[i][2];
    mR[i][2] = R[2][0]*mData[i][0] + R[2][1]*mData[i][1] + R[2][2]*mData[i][2];
  }
}





int findMax(double *dist,int nf)
{
	double max;
	int i, maxIndex;
	max = -9.999E9;
	maxIndex = -1;
	for(i=0;i<nf;i++) {
		if (max < dist[i]) {
			max = dist[i];
			maxIndex = i;
		}
	}
		return maxIndex;
}


void with_Filtering(int nNormal,int *normalIndex,double **normalXYZ,
					int nmeasured,int *measuredIndex,double **measuredXYZ,
					int *nfiltered,double **transformedXYZ,int *filteringindex,
					double *distData,double delta, bool filteredYN, bool bOnlyRotateZ)
{
  int nBase,nMeasured;

  double **mData,**bData,**mRotatedData;

  double px,py,pz;
  double stx,sty,stz;
  double End_Loop,MinValue;
  double min_px,min_py,min_pz;

  int i;

  double sum,cmx,cmy,cmz,cbx,cby,cbz;


  /* Input measured data points */
  nMeasured = nmeasured;

  mData = (double**) malloc(sizeof(double*)*(nMeasured+1));
  mRotatedData = (double **) malloc(sizeof(double*)*(nMeasured+1));
  for(i=0;i<nMeasured;i++) {

    mData[i]=(double*) malloc(sizeof(double)*3);
    mRotatedData[i]=(double*) malloc(sizeof(double)*3);

    mData[i][0]=measuredXYZ[i][0];
    mData[i][1]=measuredXYZ[i][1];
    mData[i][2]=measuredXYZ[i][2];
  }

  /* Input base data points */
  nBase = nNormal;
  bData = (double**) malloc(sizeof(double*)*(nBase+1));
  for(i=0;i<nBase;i++) {
    bData[i]=(double*) malloc(sizeof(double)*3);
    
    bData[i][0]=normalXYZ[i][0];
    bData[i][1]=normalXYZ[i][1];
    bData[i][2]=normalXYZ[i][2];
  }

  /* Compute the centroids */

  /* measured data */
  /* for x */
  sum = 0.0;
  for(i=0;i<nMeasured;i++) 
    sum += mData[i][0];
  cmx = sum / (double)nMeasured;

  /* for y */
  sum = 0.0;
  for(i=0;i<nMeasured;i++) 
    sum += mData[i][1];
  cmy = sum / (double)nMeasured;

  /* for z */
  sum = 0.0;
  for(i=0;i<nMeasured;i++) 
    sum += mData[i][2];
  cmz = sum / (double)nMeasured;


  /* base data */
  /* for x */
  sum = 0.0;
  for(i=0;i<nBase;i++) 
    sum += bData[i][0];
  cbx = sum / (double)nBase;

  /* for y */
  sum = 0.0;
  for(i=0;i<nBase;i++) 
    sum += bData[i][1];
  cby = sum / (double)nBase;

  /* for z */
  sum = 0.0;
  for(i=0;i<nBase;i++) 
    sum += bData[i][2];
  cbz = sum / (double)nBase;


  /* Move both data sets with respect to the origin */

  for(i=0;i<nMeasured;i++) {
    mData[i][0] -= cmx;
    mData[i][1] -= cmy;
    mData[i][2] -= cmz;
  }

  
  for(i=0;i<nBase;i++) {
    bData[i][0] -= cbx;
    bData[i][1] -= cby;
    bData[i][2] -= cbz;
  }

  /* Compute the various initial positions for optimization */


  MinValue = 1.0e+16;

  stx = sty = stz = 45.0;

  if (bOnlyRotateZ)
  {
	  stx = sty = 360;
	  stz = 45;
  }

  px = py = pz = 0.0;

  while(px < 360.0) {
    py = 0.0;
    while(py < 360.0) {
      pz = 0.0;
      while(pz < 360.0) {
		change_InitialPosition(px,py,pz,nMeasured,mData,mRotatedData);
		End_Loop = find_motion(bData,nBase,mRotatedData,nMeasured); 
		if(End_Loop < MinValue) {
		  min_px = px;
		  min_py = py;
		  min_pz = pz;
		  MinValue = End_Loop;
		  /*	  printf("%lf %lf %lf %lf\n",px, py, pz, End_Loop);*/
		}
		pz += stz;
      }
      py += sty;
    }
    px += stx;
  }


  //  [12/13/2011 KKY] - �߽��̵�
//   double** currData = (double**) malloc(sizeof(double*)*(nMeasured+1));
//   for(i=0;i<nMeasured;i++) 
// 	  currData[i]=(double*) malloc(sizeof(double)*3);
// 
// 	pz = 0.0;
// 	int minPos;
// 	while(pz < 360.0)
// 	{
// 		change_InitialPosition(px,py,pz,nMeasured,mData,mRotatedData);
// 		for(i=0;i<nMeasured;i++)
// 		{
// 			currData[i][0]=mRotatedData[i][0];
// 			currData[i][1]=mRotatedData[i][1];
// 			currData[i][2]=mRotatedData[i][2];
// 		}
// 
// 		for(int k=0; k<nBase+1; k++)
// 		{
// 			End_Loop = find_motion(bData,nBase,mRotatedData,nMeasured); 
// 			if(End_Loop < MinValue) {
// 			  min_px = px;
// 			  min_py = py;
// 			  min_pz = pz;
// 			  MinValue = End_Loop;
// 			  minPos = k;
// 			}
// 
// 			//���� ����Ʈ�� ���� ������ �̵�
// 			if(k < nBase)
// 			{
// 				for(i=0;i<nMeasured;i++) 
// 				{
// 					mRotatedData[i][0] = currData[i][0] + bData[k][0];
// 					mRotatedData[i][1] = currData[i][1] + bData[k][1];
// 					mRotatedData[i][2] = currData[i][2] + bData[k][2];
// 				}
// 			}
// 
// 		}
// 		pz += stz;
// 	}
// 	/* Free memory */
// 	for(i=0;i<nMeasured;i++) {
// 		free(currData[i]);
// 	}
// 	free(currData);



  /* Orientation angles for optimum initial position */

  change_InitialPosition(min_px,min_py,min_pz,nMeasured,mData,mRotatedData);

//   //  [12/13/2011 KKY] - �߽��̵�
//   if(minPos != 0)
//   {
// 	  for(i=0;i<nMeasured;i++) 
// 	  {
// 		  mRotatedData[i][0] = mRotatedData[i][0] + bData[minPos-1][0];
// 		  mRotatedData[i][1] = mRotatedData[i][1] + bData[minPos-1][1];
// 		  mRotatedData[i][2] = mRotatedData[i][2] + bData[minPos-1][2];
// 	  }
//   }

  End_Loop = find_motion_final(bData,nBase,mRotatedData,nMeasured,cbx,cby,cbz);


  moveBackBaseData(bData,nBase,cbx,cby,cbz);
  createTable(bData,nBase,mRotatedData,nMeasured,distData,filteringindex,delta,filteredYN);
  find_motion_without_link(bData,nBase,mRotatedData,&nMeasured,distData,filteringindex);
  *nfiltered = nMeasured;
  for(i=0;i<nMeasured;i++) {
    transformedXYZ[i][0] = mRotatedData[i][0];
    transformedXYZ[i][1] = mRotatedData[i][1];
    transformedXYZ[i][2] = mRotatedData[i][2];
  }
  /* ���͸� �Ǵ� ����Ÿ�� ����*/
  deleteFilteredData(measuredIndex,nfiltered,filteringindex,transformedXYZ,distData,filteredYN);
  updateCorres(bData,nBase,mRotatedData,nMeasured,distData,filteringindex,delta,filteredYN);	/* �߰��κ� 10/2/2007 */
  
  /* Free memory */
  for(i=0;i<nMeasured;i++) {
    free(mData[i]);
    free(mRotatedData[i]);
  }

  free(mData);
  free(mRotatedData);
  for(i=0;i<nBase;i++) {
    free(bData[i]);
  }
  free(bData);

  return;
}


/* �߰��� �Լ� 10/2/2007 */
void updateCorres(double **bD,int nB,double **mR,int nM,double *mD,int *mC, double delta, bool useFileter)
{
  int i,j;
  double dist;
  double minDist;


  for(i=0;i<nM;i++) {
    minDist = 1.0e+16;
    for(j=0;j<nB;j++) {
      dist = (bD[j][0]-mR[i][0])*(bD[j][0]-mR[i][0])+(bD[j][1]-mR[i][1])*(bD[j][1]-mR[i][1])+(bD[j][2]-mR[i][2])*(bD[j][2]-mR[i][2]);
      dist = sqrt(dist);
      if(dist < minDist) {
	     mC[i] = j+1;
	     mD[i] = dist;
	     minDist = dist;
      }
    }
  }
}



void deleteFilteredData(int *measuredIndex,int *nfiltered,int *filteringIndex,double **oTransformed,double *dist,bool filteredYN)
{
  /* When filteredYN == 0, do nothing */

  int finalDataIndex;
  int i;
  if(filteredYN) {
    finalDataIndex = 0;
    for(i=0;i<*nfiltered;i++) {
      if(filteringIndex[i] != -9999) {
	measuredIndex[finalDataIndex] = measuredIndex[i]+1; 
	filteringIndex[finalDataIndex] = filteringIndex[i];
/*	oTransformed[finalDataIndex][0] = oTransformed[i][0];
	oTransformed[finalDataIndex][1] = oTransformed[i][1];
	oTransformed[finalDataIndex][2] = oTransformed[i][2];
	dist[finalDataIndex] = dist[i];*/
	++finalDataIndex;
      }
    }
    *nfiltered = finalDataIndex;
  } 
}
//
//void createTable(double **bD,int nB,double **mR,int nM,double *mD,int *mC,bool filteredYN)
//{
//  int i,j;
//  double dist;
//  double minDist;
//
//
//  for(i=0;i<nM;i++) {
//    minDist = 1.0e+16;
//    for(j=0;j<nB;j++) {
//      dist = (bD[j][0]-mR[i][0])*(bD[j][0]-mR[i][0])+(bD[j][1]-mR[i][1])*(bD[j][1]-mR[i][1])+(bD[j][2]-mR[i][2])*(bD[j][2]-mR[i][2]);
//      dist = sqrt(dist);
//      if(dist < minDist) {
//	mC[i] = j+1;
//	mD[i] = dist;
//	minDist = dist;
//      }
//    }
//  }
//
// /* Marking LINK points */
//  
//  if(filteredYN) { /* Filtering flag is set. Filtering is performed.*/
//
//    //printf("\n\nEnter the delta value: ");
//    //scanf("%lf",&delta);
//
//    for(i=0;i<nM;i++) {
//      if(mD[i] > delta) 
//	mC[i] = -9999;
//    }
//  }
//
//  /* For a debugging purpose 
//  printf("+++++++++++++++++++++++ Correspondence Number to Base List\n");
//  for(i=0;i<nM;i++) 
//  printf("dist = %10.5f \t   No. = %d \n",mD[i],mC[i]);*/
//}
//
//void change_InitialPosition(double x,double y,double z,int n,double **mData,double **mR)
//{
//  double R[3][3];
//
//  double xc,xs,yc,ys,zc,zs;
//
//  int i;
//
//  xc = cos(x);xs = sin(x);
//  yc = cos(y);ys = sin(y);
//  zc = cos(z);zs = sin(z);
//
//  R[0][0] = yc*zc;
//  R[0][1] = xs*ys*zc - xc*zs;
//  R[0][2] = xc*ys*zc + xs*zs;
//  
//  R[1][0] = yc * zs;
//  R[1][1] = xs*ys*zs + xc*zc;
//  R[1][2] = xc*ys*zs - xs * zc;
//  
//  R[2][0] = -ys;
//  R[2][1] = xs*yc;
//  R[2][2] = xc*yc;
//
//  for(i=0;i<n;i++) {
//    mR[i][0] = R[0][0]*mData[i][0] + R[0][1]*mData[i][1] + R[0][2]*mData[i][2];
//    mR[i][1] = R[1][0]*mData[i][0] + R[1][1]*mData[i][1] + R[1][2]*mData[i][2];
//    mR[i][2] = R[2][0]*mData[i][0] + R[2][1]*mData[i][1] + R[2][2]*mData[i][2];
//  }
//}

double find_motion(double **bData,int nB,double **mRotatedData,int nM)
{
  double set1[NUM_OF_DATA][3],base[NUM_OF_DATA][3];
  double P0[NUM_OF_DATA][3];
  int n_set1;
  int n_base;

  int i,j,index;

  double MuP[3],MuX[3];
  double Covariance[3][3];
  double Rot[3][3];
  double qt[3];
  double Q_Matrix[4][4];

  double InputA[5][5];
  double d[5];
  double v[5][5];
  double quaternion[4];
  int nrot;
 
  double Prev_Loop,End_Loop;
  int loop;

  struct Closest_Pair *Head;
  struct Closest_Pair *Current;

  /* Input data from the disk */  
  
  loop = 0;

  Prev_Loop = 1.0e+16;
  End_Loop = 0.0;

  /* Input the measured data */
  n_set1 = nM;
  for(i=0;i<nM;i++) {
    P0[i][0]=set1[i][0]=mRotatedData[i][0];
    P0[i][1]=set1[i][1]=mRotatedData[i][1];
    P0[i][2]=set1[i][2]=mRotatedData[i][2];
  }

  /* Input the base data */
  n_base = nB;

  for(i=0;i<nB;i++) {
    base[i][0]=bData[i][0];
    base[i][1]=bData[i][1];
    base[i][2]=bData[i][2];
  }


  Head = Make_List_Of_Pair(n_set1);

  do {

    Prev_Loop = End_Loop;
    index=0;
    Find_Closest_Pair_List(set1,n_set1,base,n_base,MuP,MuX,Head);
   
    Current=Head;

    while(Current != NULL) {
      Current->First[0] = P0[index][0];
      Current->First[1] = P0[index][1];
      Current->First[2] = P0[index][2];
      MuP[0] = MuP[0] + P0[index][0];
      MuP[1] = MuP[1] + P0[index][1];
      MuP[2] = MuP[2] + P0[index++][2];
      Current=Current->next;
    }

    for(i=0;i<3;i++) 
      MuP[i] = MuP[i]/index;
    
    Calc_Covariance(Covariance,Head,MuP,MuX);
    
    Calc_Q(Q_Matrix,Covariance);
    
    for(i=0;i<4;i++) 
      for(j=0;j<4;j++) 
	InputA[i+1][j+1] = Q_Matrix[i][j];
    
    jacobi(InputA,4,d,v,&nrot);
    eigsrt(d,v,4);
    
    for(i=1;i<=4;i++) 
      quaternion[i-1]=v[i][1];
    
    Rotation_Matrix(Rot,quaternion);
    
    Make_Translation_Vector(qt,Rot,MuP,MuX);
    
    Update(Rot,qt,Head);
    
    End_Loop = Assign_Check_Final(set1,Head);
    /*    printf("The average distance : %10.5f\n",End_Loop);*/

  } while(fabs(Prev_Loop-End_Loop)>0.00001);

  Free_List(Head);

  return End_Loop;
}



double find_motion_final(double **bData,int nB,double **mRotatedData,int nM,double cbx,double cby,double cbz)
{
  double set1[NUM_OF_DATA][3],base[NUM_OF_DATA][3];
  double P0[NUM_OF_DATA][3];
 
  int n_set1;
  int n_base;

  int i,j,index;
 

 
  double MuP[3],MuX[3];
  double Covariance[3][3];
  double Rot[3][3];
  double qt[3];
  double Q_Matrix[4][4];

  double InputA[5][5];
  double d[5];
  double v[5][5];
  double quaternion[4];
  int nrot;


  double Prev_Loop,End_Loop;
  int loop;

  struct Closest_Pair *Head;
  struct Closest_Pair *Current;

  /* Input data from the disk */  
  
  loop = 0;

  Prev_Loop = 1.0e+16;
  End_Loop = 0.0;

  /* Input the measured data */
  n_set1 = nM;
  for(i=0;i<nM;i++) {
    P0[i][0]=set1[i][0]=mRotatedData[i][0];
    P0[i][1]=set1[i][1]=mRotatedData[i][1];
    P0[i][2]=set1[i][2]=mRotatedData[i][2];
  }

  /* Input the base data */
  n_base = nB;
  for(i=0;i<nB;i++) {
    base[i][0]=bData[i][0];
    base[i][1]=bData[i][1];
    base[i][2]=bData[i][2];
  }

  Head = Make_List_Of_Pair(n_set1);

  do {

    Prev_Loop = End_Loop;
    index=0;
    Find_Closest_Pair_List(set1,n_set1,base,n_base,MuP,MuX,Head);
   
    Current=Head;

    while(Current != NULL) {
      Current->First[0] = P0[index][0];
      Current->First[1] = P0[index][1];
      Current->First[2] = P0[index][2];
      MuP[0] = MuP[0] + P0[index][0];
      MuP[1] = MuP[1] + P0[index][1];
      MuP[2] = MuP[2] + P0[index++][2];
      Current=Current->next;
    }

    for(i=0;i<3;i++) 
      MuP[i] = MuP[i]/index;
    
    Calc_Covariance(Covariance,Head,MuP,MuX);
    
    Calc_Q(Q_Matrix,Covariance);
    
    for(i=0;i<4;i++) 
      for(j=0;j<4;j++) 
	InputA[i+1][j+1] = Q_Matrix[i][j];
    
    jacobi(InputA,4,d,v,&nrot);
    eigsrt(d,v,4);
    
    for(i=1;i<=4;i++) 
      quaternion[i-1]=v[i][1];
    
    Rotation_Matrix(Rot,quaternion);
    
    Make_Translation_Vector(qt,Rot,MuP,MuX);
    
    Update(Rot,qt,Head);
    
    End_Loop = Assign_Check_Final(set1,Head);
    /*    printf("The average distance : %10.5f\n",End_Loop);*/

  } while(fabs(Prev_Loop-End_Loop)>0.001);

  /*  Out_Trans_Rot(Rot,qt);  */
  Out_Result(Head,n_set1,cbx,cby,cbz,mRotatedData);

  Free_List(Head);

  return End_Loop;
}


double find_motion_without_link(double **bData,int nB,double **mRotatedData,int *nM,double *mDist, int *mCorres)
{
  double set1[NUM_OF_DATA][3],base[NUM_OF_DATA][3];
  double P0[NUM_OF_DATA][3];
 
  int n_set1;
  int n_base;

  int i,j,index;
  

  double MuP[3],MuX[3];
  double Covariance[3][3];
  double Rot[3][3];
  double qt[3];
  double Q_Matrix[4][4];

  double InputA[5][5];
  double d[5];
  double v[5][5];
  double quaternion[4];
  int nrot;
 

  double Prev_Loop,End_Loop;
  int loop;

  struct Closest_Pair *Head;
  struct Closest_Pair *Current;



  /* Input data from the disk */  
  
  loop = 0;

  Prev_Loop = 1.0e+16;
  End_Loop = 0.0;

  /* Input the measured data */
  int netNumber = 0;
  for(i=0;i<*nM;i++) {
    if(mCorres[i] != -9999) { 
      P0[netNumber][0]=set1[netNumber][0]=mRotatedData[i][0];
      P0[netNumber][1]=set1[netNumber][1]=mRotatedData[i][1];
      P0[netNumber][2]=set1[netNumber][2]=mRotatedData[i][2];
      ++netNumber;
    } 
  }

  if ( netNumber == 0 ) {
	  return End_Loop;	// Added by fish
  }

  n_set1 = netNumber;

  /* Input the base data */
  n_base = nB;
  for(i=0;i<nB;i++) {
    base[i][0]=bData[i][0];
    base[i][1]=bData[i][1];
    base[i][2]=bData[i][2];
  }

  Head = Make_List_Of_Pair(n_set1);

  do {

    Prev_Loop = End_Loop;
    index=0;
    Find_Closest_Pair_List(set1,n_set1,base,n_base,MuP,MuX,Head);
   
    Current=Head;

    while(Current != NULL) {
      Current->First[0] = P0[index][0];
      Current->First[1] = P0[index][1];
      Current->First[2] = P0[index][2];
      MuP[0] = MuP[0] + P0[index][0];
      MuP[1] = MuP[1] + P0[index][1];
      MuP[2] = MuP[2] + P0[index++][2];
      Current=Current->next;
    }

    for(i=0;i<3;i++) 
      MuP[i] = MuP[i]/index;


    Calc_Covariance(Covariance,Head,MuP,MuX);
    
    Calc_Q(Q_Matrix,Covariance);
    
    for(i=0;i<4;i++) 
      for(j=0;j<4;j++) 
	InputA[i+1][j+1] = Q_Matrix[i][j];
    
    jacobi(InputA,4,d,v,&nrot);
    eigsrt(d,v,4);
    
    for(i=1;i<=4;i++) 
      quaternion[i-1]=v[i][1];
    
    Rotation_Matrix(Rot,quaternion);
    
    Make_Translation_Vector(qt,Rot,MuP,MuX);
    
    Update(Rot,qt,Head);
    
    End_Loop = Assign_Check_Final(set1,Head);
    /*    printf("The average distance without link data : %10.5f\n",End_Loop);*/

  } while(fabs(Prev_Loop-End_Loop)>0.0001);

  Out_Result_without_link(Head,n_set1,mRotatedData,mDist);

  Free_List(Head);

  return End_Loop;
}

void Out_Result_without_link(struct Closest_Pair *Head,int n_set1,double **mR,double *mD)
{
// 
//   struct Closest_Pair *Current;
//   int index;
// 
//   double tx,ty,tz;
// 
//   /*  printf("=====================================\n Optimized distances \n===============================\n");*/
//   Current=Head;
//   index =0;
//   /*  fp = fopen("result_without_link.dat","w");*/
//   
//   /*  fprintf(fp,"%d\n",n_set1);*/
// 
//   FILE *fp2 = NULL;
// 
//   _tfopen_s(&fp2, "dist.txt", _T("w"));
//   while(Current != NULL) {
// 
//     tx = Current->First[0];
//     ty = Current->First[1];
//     tz = Current->First[2];
// 
//     mR[index][0] = tx;
//     mR[index][1] = ty;
//     mR[index][2] = tz;
// 
// 
// 
// 	mD[index] = Current->dist;
// 
// 
//     index += 1;
// 
// 	
// 
//     _ftprintf_s(fp2,_T("%10.5f %10.5f %10.5f "),tx,ty,tz);
// 
//     _ftprintf_s(fp2,_T("dist = %10.5f \n"),Current->dist);
// 
//     Current=Current->next;
//   }
//   fclose(fp2);
//   /*  printf("=================================\n");
//   printf("The total number of points: %i\n",index);
//   fclose(fp);*/
}


void Out_Result(struct Closest_Pair *Head,int n_set1,double cbx,double cby,double cbz,double **mR)
{

  struct Closest_Pair *Current;
  int index;
 
  double tx,ty,tz;

  Current=Head;
  index =0;

  /*  fp = fopen("result.dat","w");
  
  fprintf(fp,"%d\n",n_set1);*/


  while(Current != NULL) {

    tx = Current->First[0]+cbx;
    ty = Current->First[1]+cby;
    tz = Current->First[2]+cbz;

    mR[index][0] = tx;
    mR[index][1] = ty;
    mR[index][2] = tz;
    index += 1;

    /*    fprintf(fp,"%10.5f %10.5f %10.5f\n",tx,ty,tz);*/
    Current=Current->next;
  }

  /*  fclose(fp); */
}

void Out_Trans_Rot(double Rot[3][3],double qt[3])
{
//   int i;
//   FILE *fp = NULL;
//   _tfopen_s(&fp, _T("RT_Vector.dat"),_T("w"));
//   for(i=0;i<3;i++) 
//     _ftprintf_s(fp, _T("%10.5f %10.5f %10.5f\n"), Rot[i][0],Rot[i][1],Rot[i][2]);
//   _ftprintf_s(fp, _T("%10.5f %10.5f %10.5f\n"), qt[0],qt[1],qt[2]);
// 
//   fclose(fp);
}

      
  
double Assign_Check_Final(double set1[NUM_OF_DATA][3],struct Closest_Pair *Head)
{
  struct Closest_Pair* Current;
  int index;
  double distance;

  index = 0;
  distance = 0.0;
  
  Current = Head;
  while(Current != NULL) {
    set1[index][0] = Current->First[0];
    set1[index][1] = Current->First[1];
    set1[index++][2] = Current->First[2];
    distance = distance + Current->dist;
    Current=Current->next;
  }
  distance = distance / index;
  return distance;
}


void Make_Translation_Vector(double qt[3],double Rot[3][3],double MuP[3],double MuX[3])
{
  int i,j;
  qt[0] = qt[1] = qt[2] = 0.0;
  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) 
      qt[i] = qt[i] + Rot[i][j] * MuP[j]; 
  
  for(i=0;i<3;i++) 
    qt[i] = MuX[i] - qt[i];
}

void Update(double Rot[3][3],double p[3],struct Closest_Pair *Head)
{
  double temp[3];
  int i,j;
  double t_res;
  struct Closest_Pair *Current;
  
  Current = Head;

  while(Current != NULL) {
    for(i=0;i<3;i++) {
      t_res = 0.0;
      for(j=0;j<3;j++) 
	t_res = t_res + Rot[i][j] * Current->First[j];
      temp[i] = t_res;
    }
    for(i=0;i<3;i++) 
      Current->First[i] = temp[i] + p[i];
    Current=Current->next;
  }
}


void Calc_Q(double Q[4][4],double Covariance[3][3])
{
  double tr;
  double A[3][3];

  double Identity[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
  int i,j;

  tr = 0.0;

  tr = Covariance[0][0] + Covariance[1][1] + Covariance[2][2];
  Q[0][0] = tr;


  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) 
      A[i][j] = Covariance[i][j]-Covariance[j][i];
  
  Q[0][1] = A[1][2];
  Q[0][2] = A[2][0];
  Q[0][3] = A[0][1];

  Q[1][0] = A[1][2];
  Q[2][0] = A[2][0];
  Q[3][0] = A[0][1];

  for(i=0;i<3;i++) 
    for(j=0;j<3;j++) 
      Q[i+1][j+1] = Covariance[i][j]+Covariance[j][i]-tr*Identity[i][j];

}

void Calc_Covariance(double Covariance[3][3],struct Closest_Pair *Head,double MuP[3],double MuX[3])
{
  int i,j;
  double p1,p2,p3,x1,x2,x3;
  struct Closest_Pair *Current;
  int index;

  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      Covariance[i][j] = 0.0;

  Current=Head;
  index = 0;

  while(Current != NULL) {

    p1 = Current->First[0];
    p2 = Current->First[1];
    p3 = Current->First[2];
    
    x1 = Current->Base[0];
    x2 = Current->Base[1];
    x3 = Current->Base[2];
    
    Covariance[0][0] = Covariance[0][0] + p1*x1;
    Covariance[0][1] = Covariance[0][1] + p1*x2;
    Covariance[0][2] = Covariance[0][2] + p1*x3;
    
    Covariance[1][0] = Covariance[1][0] + p2*x1;
    Covariance[1][1] = Covariance[1][1] + p2*x2;
    Covariance[1][2] = Covariance[1][2] + p2*x3;
    
    Covariance[2][0] = Covariance[2][0] + p3*x1;
    Covariance[2][1] = Covariance[2][1] + p3*x2;
    Covariance[2][2] = Covariance[2][2] + p3*x3;
    index++;
    Current=Current->next;
  }
  
  Covariance[0][0] = Covariance[0][0]/index - MuP[0]*MuX[0];
  Covariance[0][1] = Covariance[0][1]/index - MuP[0]*MuX[1];
  Covariance[0][2] = Covariance[0][2]/index - MuP[0]*MuX[2];
  
  Covariance[1][0] = Covariance[1][0]/index - MuP[1]*MuX[0];
  Covariance[1][1] = Covariance[1][1]/index - MuP[1]*MuX[1];
  Covariance[1][2] = Covariance[1][2]/index - MuP[1]*MuX[2];
  
  Covariance[2][0] = Covariance[2][0]/index - MuP[2]*MuX[0];
  Covariance[2][1] = Covariance[2][1]/index - MuP[2]*MuX[1];
  Covariance[2][2] = Covariance[2][2]/index - MuP[2]*MuX[2];
  
}

void Rotation_Matrix(double Rot[3][3], double quaternion[4])
{
  double q1,q2,q3,q0;

  q0=quaternion[0];
  q1=quaternion[1];
  q2=quaternion[2];
  q3=quaternion[3];

  Rot[0][0] = q0*q0+q1*q1-q2*q2-q3*q3;
  Rot[1][0] = 2*(q1*q2+q3*q0);
  Rot[2][0] = 2*(q1*q3-q2*q0);
  
  Rot[0][1] = 2*(q1*q2-q3*q0);
  Rot[1][1] = q0*q0-q1*q1+q2*q2-q3*q3;
  Rot[2][1] = 2*(q2*q3+q1*q0);

  Rot[0][2] = 2*(q1*q3+q2*q0);
  Rot[1][2] = 2*(q2*q3-q1*q0);
  Rot[2][2] = q0*q0-q1*q1-q2*q2+q3*q3;
}

struct Closest_Pair *Make_Node(double First[3],double Base[3],double dist)
{
  struct Closest_Pair *Dummy;
  int i;
  Dummy = (struct Closest_Pair*) malloc(sizeof(struct Closest_Pair));
  for(i=0;i<3;i++) {
    Dummy->First[i] = First[i];
    Dummy->Base[i] = Base[i];
  }
  Dummy->dist = dist;
  Dummy->next = NULL;

  return Dummy;
}

void Find_Closest_Pair_List(double set1[NUM_OF_DATA][3],int n_set1,double base[NUM_OF_DATA][3],int n_base,double MuP[3],double MuX[3],struct Closest_Pair* Head)
{
  int i;
  double Points[3],Min_Pair[3];
  double Min_Dist;
  int n_Points;
  struct Closest_Pair *Current;
 
  double First[3]={0.0,0.0,0.0};
  double Base[3]={0.0,0.0,0.0};

  n_Points = 0;
  for(i=0;i<3;i++) {
    MuP[i] = 0.0;
    MuX[i] = 0.0;
  }

  Current = Head;

  for(i=0;i<n_set1;i++) {
    n_Points++;
    Points[0] = set1[i][0];
    Points[1] = set1[i][1];
    Points[2] = set1[i][2];

    Min_Dist = Find_Min_Dist_Pair(Points,Min_Pair,n_base,base);

    MuX[0] = MuX[0] + Min_Pair[0];
    MuX[1] = MuX[1] + Min_Pair[1];
    MuX[2] = MuX[2] + Min_Pair[2];

    Assign_Points(Current,Points,Min_Pair,Min_Dist);
    Current = Current->next;
  }

  for(i=0;i<3;i++) {
    MuP[i] = MuP[i]/n_set1;
    MuX[i] = MuX[i]/n_set1;
  }
}

void Assign_Points(struct Closest_Pair* Current,double Points[3],double Min_Pair[3],double Min_Dist)
{
  int i;

  for(i=0;i<3;i++) {
    Current->First[i]=Points[i];
    Current->Base[i]=Min_Pair[i];
  }
  Current->dist = Min_Dist;
}

struct Closest_Pair *Make_List_Of_Pair(int n_set)
{
  double Zero[3] = {0.0,0.0,0.0};
  int i;
  struct Closest_Pair *Head, *Current;
  Head = Make_Node(Zero,Zero,0.0);

  Current = Head;
  for(i=0;i<n_set-1;i++) {
    Current->next = Make_Node(Zero,Zero,0.0);
    Current=Current->next;
  }
  return Head;
}

void Free_List(struct Closest_Pair *Head)
{
  struct Closest_Pair *Current;
  struct Closest_Pair *temp;
  Current = Head;

  while(Current != NULL) {
    temp = Current;
    Current=Current->next;
    free(temp);
  }
}
  
    
double Find_Min_Dist_Pair(double Points[3],double Min_Pair[3],int n_base,double base[NUM_OF_DATA][3])
{
  int j;
  double Min_Dist,dist;

  Min_Dist = 1.0e+12;
  for(j=0;j<n_base;j++) {
    dist = (Points[0]-base[j][0])*(Points[0]-base[j][0]) +
      (Points[1]-base[j][1])*(Points[1]-base[j][1]) +
      (Points[2]-base[j][2])*(Points[2]-base[j][2]);
    dist = sqrt(dist);
    if(dist < Min_Dist) {
      Min_Dist = dist;
      Min_Pair[0] = base[j][0];
      Min_Pair[1] = base[j][1];
      Min_Pair[2] = base[j][2];
    }
  }
  return Min_Dist;
}

void eigsrt(double d[5], double v[5][5], int n)
{
  int k,j,i;
  double p;
  
  for (i=1;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<=n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=1;j<=n;j++) {
	p=v[j][i];
	v[j][i]=v[j][k];
	v[j][k]=p;
      }
    }
  }
}

#define NRANSI
//#include <KyBase/KyNrutil.h>
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);

void jacobi(double a[5][5], int n, double d[5], double v[5][5], int *nrot)
{
  int j,iq,ip,i;
  double tresh, theta, tau, t, sm, s, h, g, c;
  //double * b, * z;
  //b=vector(1,n);
  //z=vector(1,n);
  Eigen::VectorXd b(1, n);
  Eigen::VectorXd z(1, n);
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++)
	sm += fabs(a[ip][iq]);
    }
    if (fabs(sm) < 0.01) {
      //free_vector(z,1,n);
      //free_vector(b,1,n);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
	g=100.0*fabs(a[ip][iq]);
	if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
	    && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
	  a[ip][iq]=0.0;
	else if (fabs(a[ip][iq]) > tresh) {
	  h=d[iq]-d[ip];
	  if ((double)(fabs(h)+g) == (double)fabs(h))
	    t=(a[ip][iq])/h;
	  else {
	    theta=0.5*h/(a[ip][iq]);
	    t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
	    if (theta < 0.0) t = -t;
	  }
	  c=1.0/sqrt(1+t*t);
	  s=t*c;
	  tau=s/(1.0+c);
	  h=t*a[ip][iq];
	  z[ip] -= h;
	  z[iq] += h;
	  d[ip] -= h;
	  d[iq] += h;
	  a[ip][iq]=0.0;
	  for (j=1;j<=ip-1;j++) {
	    ROTATE(a,j,ip,j,iq)
	      }
	  for (j=ip+1;j<=iq-1;j++) {
	    ROTATE(a,ip,j,j,iq)
	      }
	  for (j=iq+1;j<=n;j++) {
	    ROTATE(a,ip,j,iq,j)
	      }
	  for (j=1;j<=n;j++) {
	    ROTATE(v,j,ip,j,iq)
	      }
	  ++(*nrot);
	}
      }
    }
    for (ip=1;ip<=n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  //nrerror("Too many iterations in routine jacobi");
}
#undef ROTATE
#undef NRANSI