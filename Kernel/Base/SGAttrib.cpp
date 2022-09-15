#include "stdafx.h"
#include <KyView/Entity/SGAttrib.h>
#include <KyView/View/RenderConfig.h>
#include <KyView/Entity/SGNode.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

CSGObjectAttrib::CSGObjectAttrib():m_pOwner(0),m_bValidBV(FALSE),m_bValidFeatureLine(FALSE),m_bEnableData(TRUE)
{	
	m_DisappearStep = 0;
}

CSGObjectAttrib::~CSGObjectAttrib() 
{ 
	 ClearHiddenObject();
}

void CSGObjectAttrib::SetName(const TCHAR* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	m_Name.FormatV(fmt, vargs);
	va_end(vargs);	
}

BOOL CSGObjectAttrib::GetBV(KyMB3D &bv)
{
	if(!m_bEnableData)
		return FALSE;

	if(!m_bValidBV)
		CalcBV();
	if(m_bValidBV)
		bv= m_BV;
	return m_bValidBV;	
}

void CSGObjectAttrib::AddSceneObject(CSGNode* pObj)
{
	std::vector<CSGNode*>::iterator it = std::find(m_SceneObjs.begin(), m_SceneObjs.end(), pObj);
	if(it==m_SceneObjs.end())
	{
		pObj->GenerateHiddenSGObject();
		m_SceneObjs.push_back(pObj);
	}
}

void CSGObjectAttrib::RemoveSceneObject(CSGNode* pObj)
{
	if (m_SceneObjs.empty())
		return;

	std::vector<CSGNode*>::iterator it = std::find(m_SceneObjs.begin(), m_SceneObjs.end(), pObj);
	if(it!=m_SceneObjs.end())
		m_SceneObjs.erase(it);
}

void CSGObjectAttrib::ClearHiddenObject()
{
	for(size_t i=0; i<m_Vertexes.size(); i++)
		delete m_Vertexes[i];
	m_Vertexes.clear();

	for(size_t i=0; i<m_Edges.size(); i++)
		delete m_Edges[i];
	m_Edges.clear();
}

void CSGObjectAttrib::TransformHiddenObject(KyTMatrix &tMat)
{
	for(size_t i=0; i<m_Vertexes.size(); i++)	
	{
		m_Vertexes[i]->Transform(tMat);
	}

	for(size_t i=0; i<m_Edges.size(); i++)
	{
		m_Edges[i]->Transform(tMat);
	}
}

void CSGObjectAttrib::NotifyHiddenObjectUpdate()
{
	for(size_t i=0; i<m_SceneObjs.size(); i++)
		m_SceneObjs[i]->GenerateHiddenSGObject();
}

void CSGObjectAttrib::GenerateFeatureObject()
{
	if(!m_bValidFeatureLine)
	{
		GenerateHiddenObject();		
		NotifyHiddenObjectUpdate();
		m_bValidFeatureLine = TRUE;
	}
}

void CSGObjectAttrib::ClearFeatureObject()
{
	ClearHiddenObject();
	NotifyHiddenObjectUpdate();
	m_bValidFeatureLine = FALSE;
}

void CSGObjectAttrib::GenerateFeatureObjectWhenReadFile()
{
	if(!m_bValidFeatureLine)
	{
		GenerateHiddenObjectWhenReadFile();		
		NotifyHiddenObjectUpdate();
		m_bValidFeatureLine = TRUE;
	}
}

CFileIFData* CSGObjectAttrib::WriteToFile()
{
	CFileIFData* pFileData = new CFileIFData;
	WriteBaseAttribToFile(pFileData);
	return pFileData;
}

void CSGObjectAttrib::WriteBaseAttribToFile(CFileIFData* pFileData)
{
	pFileData->AddValueString(_T("Base.Name"), GetName());
	pFileData->AddValueBox3d(_T("m_BV"), m_BV);
	pFileData->AddValueBoolean(_T("m_bValidBV"), m_bValidBV);
	pFileData->AddValueBoolean(_T("m_bValidFeatureLine"), m_bValidFeatureLine);	
	pFileData->AddValueBoolean(_T("m_bEnableData"), m_bEnableData);	
	pFileData->AddValueKyVertexIndexList(_T("KyVertexIntIndexes"), GetHiddenVertexes());
	pFileData->AddValueKyEdgeIndexList(_T("KyEdgeIndexes"), GetHiddenEdges());
}

void CSGObjectAttrib::ReadFromFile(CFileIFData* pFileData)
{	
	pFileData->GetValueString(_T("Base.Name"), m_Name);
	pFileData->GetValueBox3d(_T("m_BV"), m_BV);	
	pFileData->GetValueBoolean(_T("m_bValidBV"), m_bValidBV);
	pFileData->GetValueKyVertexIndexList(_T("KyVertexIntIndexes"), m_Vertexes);
	pFileData->GetValueBoolean(_T("m_bValidFeatureLine"), m_bValidFeatureLine);
 	if(!pFileData->GetValueKyEdgeIndexList(_T("KyEdgeIndexes"), m_Edges))
	{
		if(m_bValidFeatureLine)
		{
			m_bValidFeatureLine = FALSE;
			GenerateFeatureObjectWhenReadFile();
		}
	}

	pFileData->GetValueBoolean(_T("m_bEnableData"), m_bEnableData);
}
