#pragma once
#include "V3D.h"

class CSGNode;

class CSGObjectOwner
{
public:
	virtual void OnSgObjMoveDone(CSGNode* pDispObj, double x, double y) {}
};

class CSGObjectAttrib
{
	friend class CSGNode;

public:
	CSGObjectAttrib();
	virtual ~CSGObjectAttrib();
	
protected:
	CSGObjectOwner* m_pOwner;
	V3D m_BV;
	BOOL m_bValidBV;	
	std::vector<CSGNode*> m_SceneObjs;
	std::vector<KyVertex*> m_Vertexes;
	std::vector<KyEdge*> m_Edges;
	BOOL m_bValidFeatureLine;
	BOOL m_bEnableData;
public:
	UINT m_DisappearStep;

public:
	void SetDataEnable(BOOL bFlag) { m_bEnableData = bFlag; }
	BOOL IsDataEnable() { return m_bEnableData; }
	void AddSceneObject(CSGNode* pObj);
	void RemoveSceneObject(CSGNode* pObj);
	BOOL GetBV(KyMB3D &bv);		
	void InvalidBV() { m_bValidBV=FALSE; }
	CString GetName() { return m_Name; }
	void SetName(const TCHAR* fmt, ...);
	void SetAttribOwner(CSGObjectOwner* pOwner) { m_pOwner = pOwner; }
	CSGObjectOwner* GetAttribOwner() { return m_pOwner; }
	virtual void GenerateHiddenObject() {}
	virtual void GenerateHiddenObjectWhenReadFile() { GenerateHiddenObject(); }
	void ClearHiddenObject();
	virtual void TransformHiddenObject(KyTMatrix &tMat);
	void NotifyHiddenObjectUpdate();
	std::vector<KyVertex*>& GetHiddenVertexes() { return m_Vertexes; }
	std::vector<KyEdge*>& GetHiddenEdges() { return m_Edges; }	
	virtual void GenerateFeatureObject();
	virtual void ClearFeatureObject();
	virtual void GenerateFeatureObjectWhenReadFile();
	BOOL IsHasValidFeatureLine() { return m_bValidFeatureLine; }

public:
	virtual void GetSelectedPoints(std::vector<KyP3D> &pts) {}
	virtual BOOL CalcBV() { return FALSE; }
	virtual void Transform(KyTMatrix &tMat) {}	
	void WriteBaseAttribToFile(CFileIFData* pFileData);
	virtual CFileIFData* WriteToFile();
	virtual void ReadFromFile(CFileIFData* pFileData);	
};