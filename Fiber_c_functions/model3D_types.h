/*************************************************************************/
/*                                                                       */
/*                      Copyright © 2018                                 */
/* The Board of Trustees of Purdue University,  All Rights Reserved      */
/*                       Author: Tianyu Li                               */
/*                                                                       */
/*************************************************************************/
//  CylinderModel3D.hpp
//  HybridMPP_MFR
//  Created by Tianyu Li on 5/17/18.

#ifndef model3D_types_h
#define model3D_types_h

#include <stdio.h>
#include "opencv2/opencv.hpp"
#include <fstream>
#endif /* model3D_types_hpp */


#define  TESTMERGE 0

#define INNERU_T  127.5//120.5//128
#define OUTER_T 180  //128
#define DIFU_T 17.5//15.5 //17.5
#define DTHRED 50
#define INNCOV 0.7

//8

#define PI_L 3.14159265358

#define FIX_JUNC_R 4

#define NMINR 2
#define NDIFR 3
#define NMAXR 5    ///try 6 !!!!

#define NMINH 2
#define NDIFH 6
#define NMAXH 8

#define NDEGREE_Y 89
#define NDEGREE_Z 89

#define DEGREE_Y 90
#define DEGREE_Z 90



#define CYLINDER_ENG_FT   0.68//0.24//0.46   //0.46//0.22


#define CYLINDER_ENG_T   0.46//0.24//0.46   //0.46//0.22
#define CYLINDER_ENG_BT  0.385//0.18//0.36//0.38//0.2
#define BirthCylinder_ENG 0.48//0.48//0.24//0.36
#define OVERLAP_T 0.26
#define OVERLAP_T1 0.13
#define OVERLAP_GT 0.1

#define ITER_NUM 41

#define CONNECTPRIOR3D 0.75//72
#define CONNECTFINEP3D 0.15
#define LENPRIOR3D 0.12

#define GrowThresh 0.9
#define DoubleConnBonus -1.8
#define DoubleConnThresh 0.85


#define T_LOSSCON_MIN 0
#define T_LOSSCON_MAX 0.25//0.12//0.16//0.15
struct PreCMark
{
    double m_dFCenX;
    double m_dFCenY;
    double m_dFCenZ;
    double m_dBCenX;
    double m_dBCenY;
    double m_dBCenZ;

};

class CMark    //Mask
{
public:
    CMark(int nR,int nH, int nthetaY, int nthetaZ):m_nR(nR),m_nH(nH),m_nthetaY(nthetaY),m_nthetaZ(nthetaZ){}
    CMark(){};
    bool operator> (const CMark& markB) const;
    bool operator< (const CMark& markB) const;
    
    int m_nR;
    int m_nH;
    int m_nthetaY; //
    int m_nthetaZ;
    
};

class Que3D
{
public:
    Que3D(){m_pQueHead= NULL; m_pQueTail=NULL; m_bFrontFil =false; m_bBackFil= false; m_bIsMerged = false;}
    Que3D(double dX, double dY, double dZ, CMark mark, PreCMark pmark):m_dX(dX),m_dY(dY),m_dZ(dZ), m_mark(mark), m_Pmark(pmark){m_pQueHead= NULL; m_pQueTail=NULL; m_bFrontFil =false; m_bBackFil= false;m_bFrontM=false; m_bBackM=false; m_bIsMerged=false;}
 /*   Que3D& operator = (Que3D& qobject)
    {
        this->m_dX = qobject.m_dX;
        this->m_dY = qobject.m_dY;
        this->m_dZ = qobject.m_dZ;
        
        this->m_mark = qobject.m_mark;
        this->m_Pmark = qobject.m_Pmark;
        
        this->m_bFrontFil = qobject.m_bFrontFil;
        this->m_bBackFil = qobject.m_bBackFil;
        this->m_bFrontM = qobject.m_bFrontM;
        this->m_bBackM =qobject.m_bBackM;
        
        this->m_dOverlapR = qobject.m_dOverlapR;
        this->m_dDateEng = qobject.m_dDateEng;
        this->m_dPriorCon = qobject.m_dPriorCon;
        
        this->m_dFConnR = qobject.m_dFConnR;
        this->m_dBConnR = qobject.m_dBConnR;
        
        this->m_pQueHead = qobject.m_pQueHead;
        this->m_pQueTail =qobject.m_pQueTail;
        
        
        return *this;
    }*/
    
    double m_dX;
    double m_dY;
    double m_dZ;
    CMark m_mark;
    PreCMark m_Pmark;
    
    bool m_bFrontFil;
    bool m_bBackFil;
    
    
    bool m_bFrontM;   //Do need merge?
    bool m_bBackM;

    
    double m_dOverlapR;
    double m_dDateEng;
    double m_dPriorCon;
    
    double m_dFConnR;
    double m_dBConnR;
    

    std::vector<int> m_vnNeighborIDs;
    bool m_bIsMerged;
    Que3D* m_pQueHead;
    Que3D* m_pQueTail;
    
};
class TPoint
{
public:
    TPoint(int nX, int nY, int nZ, int nType):m_nX(nX),m_nY(nY),m_nZ(nZ),m_nType(nType){}
    TPoint(){};
    
    int m_nX;
    int m_nY;
    int m_nZ;
    int m_nType;
    
};

class TDPoint
{
public:
    TDPoint(double dX, double dY, double dZ):m_dX(dX),m_dY(dY),m_dZ(dZ){}
    TDPoint(){};
    
    double m_dX;
    double m_dY;
    double m_dZ;

};

struct SetParam3D
{
    int nMaxR;      //for ellipse nMaxA is the max value of x axis
    int nMinR;
   
    
    int nMaxH;      // ......... max value of y axis
    int nMinH;
    int nMinTubeH;
    int nDegree;
    int nMDegree;
    int nDegreeY;
    int nDegreeZ;
    SetParam3D()
    {
        nMaxR = NMAXR;
        nMinR = NMINR;
        nMaxH = NMAXH;
        nMinH = NMINH;
        nMinTubeH = 2;
        nDegree = 90;    //we need 180 degree, but consider the memory limit, we divide 180 to 45
        nMDegree=0;
        nDegreeY =DEGREE_Y;
        nDegreeZ =DEGREE_Z;
    }
};


struct sRGB
{
    sRGB(unsigned char R, unsigned char G, unsigned char B):m_R(R),m_G(G),m_B(B){}
    sRGB()
    {
        m_R=0;
        m_B=0;
        m_G=0;
    }
    unsigned char m_R;
    unsigned char m_G;
    unsigned char m_B;
};

struct MppControPara
{
    double dBirthTemp;
    double dBirthERatio;
    MppControPara(){dBirthTemp=0.88; dBirthERatio=0.999;}
};

struct FiberInfo
{
    
    int m_nTotalFiber;
    double m_dMeanDia;  //mean diameter
    double m_dMaxDia; //max diameter
    double m_dMinDia; //min diameter
    
    double m_dMeanLen; //mean length
    double m_dMaxLen; // max length
    double m_dMinLen; // min length
    
    int m_nHistLen[500];// Hist Len
    int m_nHistDia[17];// Hist Dia
    int m_nHistAngle[73]; //Hist Y

    double m_dHistLen[500];// Hist Len
    double m_dHistDia[9];// Hist Dia
    double m_dHistAngle[73]; //Hist Y
    
    
    double m_dRM2F;  //matrix to fiber
    
    
    FiberInfo()
    {
        
        m_nTotalFiber=0;
        
        m_dMeanDia=0;
        m_dMaxDia=0;
        m_dMinDia=0;
        
        
        m_dMeanLen=0;
        m_dMaxLen=0;
        m_dMinLen=0;
        
        memset((char*)m_nHistLen, 0, sizeof(int)*500);
        memset((char*)m_nHistDia, 0, sizeof(int)*17);
        memset((char*)m_nHistAngle, 0, sizeof(int)*73);
    }
    
};

void OutputObjInfo(Que3D* qObj, std::string strPreInfo);
int CalAngleFromYZ(int nTy, int nTz);
void NormalizeHist(int* pnHist, double* pdHist, int nChan);
void SaveObjStatInfo(FiberInfo& FI,std::string strInfor);


