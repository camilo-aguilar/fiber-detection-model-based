/*************************************************************************/
/*                                                                       */
/*                      Copyright © 2018                                 */
/* The Board of Trustees of Purdue University,  All Rights Reserved      */
/*                       Author: Tianyu Li                               */
/*                                                                       */
/*************************************************************************/
//  CylinderModel3D.hpp
//  HybridMPP_MFR
//
//  Created by Tianyu Li on 5/17/18.
//

#ifndef CylinderModel3D_h
#define CylinderModel3D_h

#include <stdio.h>
#include "opencv2/opencv.hpp"
#include "model3D_types.h"
#include <fstream>
#endif /* CylinderModel3D_hpp */

#define MAXSequence 1500
#define SIZEQue2D 16384//128*128
#define AddOne3DP(x,y,z,W,pnadd) pnadd[(z)][(y)*(W)+(x)]=(pnadd[(z)][(y)*(W)+(x)]+1)
#define DedOne3DP(x,y,z,W,pnadd) pnadd[(z)][(y)*(W)+(x)]=(pnadd[(z)][(y)*(W)+(x)]-1)
#define Get3DPixel(x,y,z,W,pchadd) pchadd[(z)][(y)*(W)+(x)]
#define Set3DPixel(x,y,z,W,pchadd,Value) pchadd[(z)][(y)*(W)+(x)]=(Value)
#define Set3DPixelR(x,y,z,W,pchadd,R) pchadd[(z)][(y)*(W)*3+(x)*3+2] = (R)
#define Set3DPixelG(x,y,z,W,pchadd,G) pchadd[(z)][(y)*(W)*3+(x)*3+1] = (G)
#define Set3DPixelB(x,y,z,W,pchadd,B) pchadd[(z)][(y)*(W)*3+(x)*3+0] = (B)
#define SqDist2Points3D(x1,y1,z1,x2,y2,z2) (((x1)-(x2))*((x1)-(x2))+((y1)-(y2))*((y1)-(y2))+((z1)-(z2))*((z1)-(z2)))
#define Dist2Points3D(x1,y1,z1,x2,y2,z2) sqrt(((x1)-(x2))*((x1)-(x2))+((y1)-(y2))*((y1)-(y2))+((z1)-(z2))*((z1)-(z2)))
#define IsInnerPoint(x,y,z,W,H,Layers) (((x)>=0&&(x)<(W))&& ((y)>=0&&(y)<(H))&& ((z)>=0&&(z)<(Layers)))


#define RandomX(PT,Range)  { (PT).m_dX=(PT).m_dX+random2()*2*(Range)-(Range); (PT).m_dY=(PT).m_dY+random2()*2*(Range)-(Range); (PT).m_dZ=(PT).m_dZ+random2()*2*(Range)-(Range);}

void Test3DFramework();
void TestFiberDetection(int SubVn);

void Load3DGrayImages(std::string strPre,int nLayers, cv::Mat* matSequence, unsigned char** ppuch3Dadd, bool bIsTif=true);      // name format strPre+[1,nLayers]   1~nLayers;

void CreateOneCylinder(CMark markC,std::map<CMark,std::vector<TPoint> >& CylinderPool, std::map<CMark,PreCMark>* pPredPool);
void SaveCylinderPool(std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string fileName="CylinderPool.data");
void ReadCylinderPool(std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string fileName="CylinderPool.data");


void CreateVoid3DColorImages(int nLayers, int nWidth, int nHeight, cv::Mat* matSequence, unsigned char** ppuch3Dadd);
void CreateVoid3DGrayImages(int nLayers, int nWidth, int nHeight, cv::Mat* matSequence, unsigned char** ppuch3Dadd);

void Clean3DGrayImages(int nLayers, int nWidth, int nHeight, unsigned char** ppuch3Dadd);
void Clean3DIntImages(int nLayers, int nWidth, int nHeight, int** ppn3Dadd);


void CreateVoid3DIntImages(int nLayers, int nWidth, int nHeight, cv::Mat* matSequence, int** ppn3Dadd);
void DrawCylinder(int nCenX, int nCenY,int nCenZ, CMark& Cmark, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);
void DrawCylinderColor(int nCenX, int nCenY,int nCenZ, CMark& Cmark, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb,int ntype=1);
void DrawCylHDirectColor(Que3D* pStartQue, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb,int nBreak=-1);
// Camilo
int DrawCylHDirectColor_Modified(double *R,double *H,double *Ty,double *Tz,double *fx,double *fy,double *fz, Que3D* pStartQue, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,int id, int nBreak=-1);
void DrawFiberResults3D(std::vector<Que3D*>& vFiber,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);



void Save3DGrayImages(std::string strPre,int nLayers, cv::Mat* matSequence);
void Save3DImages(std::string strPre,int nLayers, cv::Mat* matSequence);

void Convert3DGrayToColor(int nLayers, cv::Mat* matSequence, cv::Mat* matColorSeq, unsigned char** puchColor3D);



Que3D* GrowOneSide(Que3D* pStartQue, unsigned char** ppuchGray,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bIsHead)
;


double CalculateEnergy(Que3D* pObject, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool =NULL);

double OverlapRatio3D(Que3D* pObject, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);

double OverlapRatio3D(Que3D* pObject, int** ppn3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool=NULL);


Que3D* CreateStraighTubes(Que3D* pStartQue, int nNumTubes, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);

bool SelectDirect(TDPoint& OneSideCent, PreCMark& newPreMark); // determine the slide direction of new tube object, if true, slide to the PreMark.nfront direction, else slide to the PreMark.nback direction, based on the distance between new object and current object, the farther one is selected
void AdaptOneLocal3D(Que3D* qOneObj,unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, int nRange=2);

Que3D* GrowOneSideT(Que3D* pStartQue,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bIsHead);


/**********************************Synthsize******************************/

void SynFibers3DSeq();
void SynCylinder(int nLen, Que3D* qStartObj,CMark cmark,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsBoth=false,int nBreak=-1);

void SynCylinderData1(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides);

void SynCylinderData3(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides);

void SynCylinderData4(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides);

void SynCylinderData6(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides);    //sigle fiber

void SynCylinderData7(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides);

/**************************************************************************/


void GrayAndNoiseBG_GRAY(unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayer);

void CylinderMPPBirthDeath(std::vector<Que3D*>& vFiber,unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, int SubVn);


bool cmp3D(Que3D* m, Que3D* n);

void GrowStage3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);



void BirthStage3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);

void DeathStage3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP,bool bIsOverOnly=true);

void AdaptLocal3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);


void CalculateEnergy(Que3D* pObject, unsigned char** ppuch3Dadd,int** ppnConn,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,double& dDataE,double& dConnFRatio,double& dConnBRatio, double& dOvelapE);

void CalculateEnergy(Que3D* pObject, unsigned char** ppuch3Dadd,int** ppnConn,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,double& dDataE,double& dConnFRatio,double& dConnBRatio);
void AddObjectOverLap_CONN(Que3D* pObject,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);



void EraseObjectOverLap_CONN(Que3D* pObject,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);

void UpdatePriorConnect3D(std::vector<Que3D*>& vObjects, int** ppn3DConnect,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool);

void AdaptOneLocal3D(Que3D* qOneObj,unsigned char** ppuch3DGray, int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, int nRange=2);

void DrawResult3D(std::vector<Que3D*>& vFiber,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, bool bSingleColor=false);

void SaveResult3D(cv::Mat* pmatColorSeq,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, int nFrame);



void EvaluatingResults(double& dMissingRate, double& dFalseDRate, cv::Mat* pMat3DSeg, unsigned char** ppuch3DSeg, unsigned char** ppuch2DSeg,int nWidth, int nHeight, int nLayers);

void EvaluatingSegResultsFiles(double& dMissingRate, double& dFalseDRate, int nWidth, int nHeight, int nLayers, std::string str3DSeg, std::string str2DSeg);


double dPriorCaculation(double dConnFRatio, double dConnBRatio);


void GetCoarseSeg3DImages(int nLayers, cv::Mat* matOrigin,cv::Mat* matSeg, double dThresh);

/**********************Merge And Split Kernel******************************/

void MergeSplit3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,unsigned char** ppuch3DConLocal,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);

void AddObjectLocalCONN(Que3D* pObject,unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,int ntype=3);

void EraseObjectLocalCONN(Que3D* pObject,unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);

int CheckConnObjIJ(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bCheckF=true,int nDistance=101);      //0: No Connection, 1: front-front, 2: front-back, 3:back-front, 4,back-back


void ConvertFromXYZ2DYZ(double& dTY, double& dTZ,double dx,double dy, double dz);   // Convert from direction (x,y,z) to theta YZ

void ConvertFromXYZ2DYZ_N(int& nTY, int& nTZ,double dx,double dy, double dz);

void DrawObjByGiven2Points(int x1,int y1, int z1,
                           int x2,int y2, int z2,
                           unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);


bool ReplaceObjects(int i,TDPoint ptMid, TDPoint ptOrigin,std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,unsigned char** ppuch3DConLocal,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);


void SplitByPoints(std::vector<Que3D>& vSplitObjects,int nR, TDPoint ptMid, TDPoint ptOrigin,unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);


void SaveLooseConn(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string strSaveAdd,int nConnValue);

void SaveDrawObj(std::vector<Que3D>& vSplitObjects,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string strSaveAdd);
/**************************************************************************/



/**********Post Process, Check connection and create new (long) object queues*******/

void MergeShortToLongNew(std::vector<Que3D*>& vp3DStartLong, std::vector<Que3D*>& vp3DQues,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,unsigned char** ppuch3DConLocal=NULL);

void MergeShortToLong(std::vector<Que3D*>& vp3DStartLong, std::vector<Que3D*>& vp3DQues,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,unsigned char** ppuch3DConLocal=NULL);

void Convert2Long(std::vector<Que3D*>& vp3DStartLong, std::vector<Que3D*>& vp3DLongF,unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);

Que3D ConvertOne2Long(Que3D* pStartQue,TDPoint& TDBegin, TDPoint& TDEnd,int& nIsSingle);


double DistP2P(Que3D* pQobject1, Que3D* pQobject2, int nConCode,std::map<CMark,PreCMark>* pPredPool);


void DrawMergeShortResult(std::vector<Que3D*>& vp3DStartLong, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayeris, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);

void DrawMergeShortID(std::vector<Que3D*>& vp3DStartLong, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayeris, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);
void DrawMergeShortID_modified(std::vector<Que3D*>& vp3DStartLong, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayeris, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, int SubVn);


int MCheckConnObjIJ(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bCheckF=true,int nDistance=101);

int CMCheckConnObjIJ(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bCheckF=true,unsigned char** ppuch3DConLocal=NULL);

bool CheckSpecialAngle(int nTY1,int nTZ1, int nTY2, int nTZ2);

bool CheckSpecialAngleBigD(int nTY1,int nTZ1, int nTY2, int nTZ2);


void DrawLongFiberB2PT(unsigned char** ppuch3Dimgs, int nWidth,int nHeight, int nLayer, TDPoint TDBegin, TDPoint TDEnd,int nR,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);
void DrawMergeFiber(Que3D* pStartQue, unsigned char** ppuch3Dimgs, int nWidth,int nHeight, int nLayer,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);

double CalMatchScore(unsigned char** ppuch3DimgsA,unsigned char** ppuch3DimgsB,int nWidth,int nHeight, int nLayer);

void GrowStage3DFinal(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);

void SplitLongFiber(std::vector<Que3D*>& vp3DStartLong,std::vector<Que3D*>& vp3DFinal, unsigned char** ppuch3DGray,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP);

void CalculateNeigborID(std::vector<Que3D*>& vObjects, int nSqNbDist);


/************************************************************************************/
//**stochastic the fiber information
void WriteFiberResults(std::vector<Que3D*>& vp3DFinal, int SubVn);
void ReadFiberResults(std::vector<Que3D*>& vp3DFinal, std::string fileName);
void ShowSavedFibers(int nWidth, int nHeight, int nLayer, std::string fileName, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool);




