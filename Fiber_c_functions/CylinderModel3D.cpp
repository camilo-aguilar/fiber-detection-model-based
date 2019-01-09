/*************************************************************************/
/*                                                                       */
/*                      Copyright © 2018                                 */
/* The Board of Trustees of Purdue University,  All Rights Reserved      */
/*                       Author: Tianyu Li                               */
/*                                                                       */
/*************************************************************************/
//  CylinderModel3D.cpp
//  HybridMPP_MFR
//  Created by Tianyu Li on 5/17/18.
//

#include "CylinderModel3D.h"
#include  <ctime>
#include<algorithm>
#include "randlib.h"
#define ABS(A) ((A)>0?(A):-(A))

int nGlobalCnt =49;

bool cmp3D(Que3D*m, Que3D*n)
{
    return (m->m_dDateEng+m->m_dPriorCon> n->m_dDateEng+n->m_dPriorCon);
    // return (m.dlikelyhood > n.dlikelyhood);
    
}
void Load3DGrayImages(std::string strPre,int nLayers, cv::Mat* matSequence, unsigned char** ppuch3Dadd,bool bIsTif)      // name format strPre+[1,nLayers]   1~nLayers
{
    char chLayer[25];
    std::string strFull;
    for(int i=1; i<=nLayers; i++)
    {
        
        if(bIsTif)
             if(i<10)
               sprintf(chLayer, "000%d.tif", i);
             else if(i<100)
               sprintf(chLayer,"00%d.tif",i);
             else if(i<1000)
               sprintf(chLayer,"0%d.tif",i);
             else
               sprintf(chLayer,"%d.tif",i);
        else
             sprintf(chLayer, "%d.jpg", i-1);
        
         
        strFull=strPre + chLayer;
        cv::Mat matRead = cv::imread(strFull);
        
        std::cout<<strFull<<'\n';
        
        cv::cvtColor(matRead, matSequence[i-1], CV_BGR2GRAY);
        
        ppuch3Dadd[i-1]=matSequence[i-1].data;
        
       /* cv::imshow("mat", matSequence[i-1]);
        cv::waitKey(0);*/
        
    }
    
}

void CreateOneCylinder(CMark markC,std::map<CMark,std::vector<TPoint> >& CylinderPool, std::map<CMark,PreCMark>* pPredPool)
{
    int nH = markC.m_nH;
    int nR = markC.m_nR;
    int nThetaY = markC.m_nthetaY;
    int nThetaZ = markC.m_nthetaZ;
    
    double dRs = (nR+1.2)*(nR+1.2);
    double dRIs = (nR)*(nR);
    int ndistance = 0;
    std::vector<TPoint> vTPCylinder;
    
    for(int nz = -nH; nz<=nH; nz++)
        for(int ny = -nR-2; ny<= nR+2; ny++)
            for(int nx= -nR-2; nx <= nR+2; nx++)
            {
                ndistance = nx*nx+ny*ny;
                if(ndistance<= dRIs)    //change from <= to <  on July 2
                {
                    vTPCylinder.push_back(TPoint(nx,ny,nz,1));
                }
                else if(ndistance<= dRs)
                {
                    vTPCylinder.push_back(TPoint(nx,ny,nz,2));
                }
            }
    
    double dThetaY = 2*nThetaY/180.0*PI_L;
    double dThetaZ = 2*nThetaZ/180.0*PI_L;
    //Rotation
    
    int nSize = vTPCylinder.size();
    for (int ni=0; ni< nSize; ni++)
    {
        int nx = vTPCylinder[ni].m_nX;
        int ny = vTPCylinder[ni].m_nY;
        int nz = vTPCylinder[ni].m_nZ;
        /*   double dx = cos(dThetaZ)*(double(nx)*cos(dThetaY)+double(nz)*sin(dThetaY))+double(ny)*sin(dThetaZ);
         double dy = -sin(dThetaZ)*(double(nx)*cos(dThetaY)+double(nz)*sin(dThetaY))+double(ny)*cos(dThetaZ);
         double dz = -sin(dThetaY)*double(nx) +double(nz)*cos(dThetaY);
         */ //recorrect on June 24
        
        double dx = cos(dThetaZ)*(double(nx)*cos(dThetaY)+double(nz)*sin(dThetaY))-double(ny)*sin(dThetaZ);
        double dy = sin(dThetaZ)*(double(nx)*cos(dThetaY)+double(nz)*sin(dThetaY))+double(ny)*cos(dThetaZ);
        double dz = -sin(dThetaY)*double(nx) +double(nz)*cos(dThetaY);
        
        
        
        
        nx= round(dx); ny=round(dy); nz=round(dz);
        vTPCylinder[ni].m_nX = nx;
        vTPCylinder[ni].m_nY = ny;
        vTPCylinder[ni].m_nZ = nz;
        
        for(int nnx=0; nnx<=1; nnx++)
        {
            if(nnx==1&&nx==dx)
                continue;
            for(int nny=0; nny<=1; nny++)
            {
                if(nny==1&&ny==dy)
                    continue;
                for(int nnz=0; nnz<=1;nnz++)
                {
                    if(nnz==1&&nz==dz)
                        continue;
                    if(nnx!=0||nny!=0||nnz!=0)
                    {
                        vTPCylinder.push_back(TPoint(nx+nnx,ny+nny,nz+nnz,vTPCylinder[ni].m_nType));
                    }
                }
            }
        }
        
    }
    
    std::vector<TPoint> vTPCylinderF;
    
    for(int ni=0; ni<vTPCylinder.size(); ni++)
    {
        if(vTPCylinder[ni].m_nType==0)
            continue;
        vTPCylinderF.push_back(vTPCylinder[ni]);
        for(int nj=0; nj<vTPCylinder.size(); nj++)
        {
            
            if(vTPCylinder[ni].m_nX==vTPCylinder[nj].m_nX&&
               vTPCylinder[ni].m_nY==vTPCylinder[nj].m_nY&&
               vTPCylinder[ni].m_nZ==vTPCylinder[nj].m_nZ)
            {
                vTPCylinder[nj].m_nType=0;
            }
        }
    }
    
    PreCMark preMark;
    
    
    
    
    
    
    preMark.m_dFCenX = cos(dThetaZ)*(double(0)*cos(dThetaY)+double(nH+0.8)*sin(dThetaY))-double(0)*sin(dThetaZ);
    
    preMark.m_dFCenY = sin(dThetaZ)*(double(0)*cos(dThetaY)+double(nH+0.8)*sin(dThetaY))+double(0)*cos(dThetaZ);
    
    preMark.m_dFCenZ = -sin(dThetaY)*double(0) +double(nH+0.8)*cos(dThetaY);
    

    
    preMark.m_dBCenX = cos(dThetaZ)*(double(0)*cos(dThetaY)+double(-nH-0.8)*sin(dThetaY))-double(0)*sin(dThetaZ);
    
    preMark.m_dBCenY = sin(dThetaZ)*(double(0)*cos(dThetaY)+double(-nH-0.8)*sin(dThetaY))+double(0)*cos(dThetaZ);
    
    preMark.m_dBCenZ = -sin(dThetaY)*double(0) +double(-nH-0.8)*cos(dThetaY);
    
    
    
    
    
    (*pPredPool)[CMark(nR,nH,nThetaY,nThetaZ)] = preMark;
    
    int nFCenX = int(0.5+preMark.m_dFCenX);
    int nFCenY = int(0.5+preMark.m_dFCenY);
    int nFCenZ = int(0.5+preMark.m_dFCenZ);
    int nBCenX = int(0.5+preMark.m_dBCenX);
    int nBCenY = int(0.5+preMark.m_dBCenY);
    int nBCenZ = int(0.5+preMark.m_dBCenZ);
    int nRP2 =  nR*nR; // Let's Fix the junction area for different R
    /////Save the junction of front end
    for(int nni=-nR; nni<=nR; nni++)
        for(int nnj=-nR; nnj<=nR; nnj++)
            for(int nnz=-nR;nnz<=nR; nnz++)
            {
                if(nni*nni+nnj*nnj+nnz*nnz<=nRP2)
                {
                    vTPCylinderF.push_back(TPoint(nFCenX+nni,nFCenY+nnj,nFCenZ+nnz,3));
                    vTPCylinderF.push_back(TPoint(nBCenX+nni,nBCenY+nnj,nBCenZ+nnz,4));
                }
                
            }
    ///////////

    CylinderPool[CMark(nR,nH,nThetaY,nThetaZ)] = vTPCylinderF;
    
    
}


void SaveCylinderPool(std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string fileName)

{
    std::ofstream outfile;
    outfile.open(fileName.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    int nCnt=0;
    for(std::map<CMark,std::vector<TPoint> >::iterator iter = CylinderPool.begin(); iter!=  CylinderPool.end();++iter)
    {
        CMark tempMark = iter->first;
        PreCMark cpremark = (*pPredPool)[tempMark];
        
        std::vector<TPoint> temptVector = iter->second;
        int nSize = temptVector.size();
       // std::cout<< "The "<< nCnt++<<" size is"<<nSize<<'\n';

        outfile.write(((char*)&tempMark), sizeof(CMark));
        
        outfile.write(((char*)&cpremark.m_dFCenX),sizeof(double));
        outfile.write(((char*)&cpremark.m_dFCenY),sizeof(double));
        outfile.write(((char*)&cpremark.m_dFCenZ),sizeof(double));
        
        outfile.write(((char*)&cpremark.m_dBCenX),sizeof(double));
        outfile.write(((char*)&cpremark.m_dBCenY),sizeof(double));
        outfile.write(((char*)&cpremark.m_dBCenZ),sizeof(double));
        
        
        
        outfile.write(((char*)&nSize),sizeof(int));
        for(int i = 0; i<nSize; i++)
        {
            TPoint tp = temptVector[i];
            outfile.write(((char*)&tp),sizeof(TPoint));
        }
    }
    
    outfile.close();
    
}

void ReadCylinderPool(std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string fileName)

{
    std::ifstream infile;
    infile.open(fileName.c_str(), std::ios::in|std::ios::binary);
    int nCnt2=0;
   // int Test[10000];
    
    CylinderPool.clear();
    while(!infile.eof())
    {
        CMark tempMark;
        PreCMark cpremark;
        std::vector<TPoint> temptVector;
        temptVector.clear();
        infile.read(((char*)&tempMark), sizeof(CMark));
        
        
        infile.read(((char*)&cpremark.m_dFCenX),sizeof(double));
        infile.read(((char*)&cpremark.m_dFCenY),sizeof(double));
        infile.read(((char*)&cpremark.m_dFCenZ),sizeof(double));
        
        infile.read(((char*)&cpremark.m_dBCenX),sizeof(double));
        infile.read(((char*)&cpremark.m_dBCenY),sizeof(double));
        infile.read(((char*)&cpremark.m_dBCenZ),sizeof(double));
        
        
        
        
        
        
        int nSize=0;
        infile.read(((char*)&nSize),sizeof(int));
       // std::cout<< "The "<< nCnt2<<" size is"<<nSize<<'\n';
        
        if(nSize==0)
        {
            break;
        }
        //Test[nCnt2++]=nSize;
        for(int i=0; i<nSize; i++)
        {
            TPoint tp;
            infile.read(((char*)&tp),sizeof(TPoint));
            temptVector.push_back(tp);
        }
        CylinderPool[tempMark] = temptVector;
        (*pPredPool)[tempMark] = cpremark;
    }
    infile.close();
    int nCnt=0;
   // for(std::map<CMark,std::vector<TPoint> >::iterator iter = CylinderPool.begin(); iter!=  CylinderPool.end();++iter)
    {
        
     //   if(iter->second.size()!=Test[nCnt++])
       //     std::cout<<"There is inequlity!";
       // std::cout<<"The "<<nCnt<<" Mark is "<< iter->second.size()<<'\n';
       // CMark tempMark = iter->first;
       // std::cout<<"nR="<<tempMark.m_nR<<" nH="<<tempMark.m_nH<<"\n";
       // std::cout<<"nThetaY="<<tempMark.m_nthetaY<<" nThetaZ="<<tempMark.m_nthetaZ<< "NumofPoints="<< iter->second.size()<<'\n';
        
    }
    
}

void CreateVoid3DGrayImages(int nLayers, int nWidth, int nHeight, cv::Mat* matSequence, unsigned char** ppuch3Dadd)
{
    for(int i=0; i<nLayers; i++)
    {
        matSequence[i] = cv::Mat::zeros(nHeight,nWidth,CV_8UC1);
        ppuch3Dadd[i]= matSequence[i].data;
        
    }
}
void CreateVoid3DColorImages(int nLayers, int nWidth, int nHeight, cv::Mat* matSequence, unsigned char** ppuch3Dadd)
{
    for(int i=0; i<nLayers; i++)
    {
        matSequence[i] = cv::Mat::zeros(nHeight,nWidth,CV_8UC3);
        ppuch3Dadd[i]= (unsigned char*)matSequence[i].data;
    }
}

void CreateVoid3DIntImages(int nLayers, int nWidth, int nHeight, cv::Mat* matSequence, int** ppn3Dadd)
{
    for(int i=0; i<nLayers; i++)
    {
        matSequence[i] = cv::Mat::zeros(nHeight,nWidth,CV_32SC1);
        ppn3Dadd[i]= (int*)matSequence[i].data;
        
    }
}

void DrawCylinder(int nCenX, int nCenY, int nCenZ,CMark& Cmark, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
   // std::cout<<"Draw x="<<nCenX<<"y="<<nCenY<<"z="<<nCenZ<<'\n';
    std::vector<TPoint> CylinderPoints = CylinderPool[Cmark];
    for(int i=0; i<CylinderPoints.size(); i++)
    {
        int nx = CylinderPoints[i].m_nX+nCenX;
        int ny = CylinderPoints[i].m_nY+nCenY;
        int nz = CylinderPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        if(CylinderPoints[i].m_nType==2)
        {
           Set3DPixel(nx, ny, nz, nWidth, ppuch3Dadd, 255);
        }
        
    }
}

void DrawCylinderColor(int nCenX, int nCenY,int nCenZ, CMark& Cmark, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb,int ntype)
{
    
    
    
    
    if(CylinderPool.find(Cmark)==CylinderPool.end())
    {
        
        //         OutputObjInfo(qLongNew, "Long Fiber");
        CreateOneCylinder(Cmark, CylinderPool, pPredPool);
    }
    
    
    std::vector<TPoint> CylinderPoints = CylinderPool[Cmark];

    for(int i=0; i<CylinderPoints.size(); i++)
    {
        int nx = CylinderPoints[i].m_nX+nCenX;
        int ny = CylinderPoints[i].m_nY+nCenY;
        int nz = CylinderPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        if(CylinderPoints[i].m_nType==ntype)//||CylinderPoints[i].m_nType==2)
        {
            
        
            Set3DPixelR(nx, ny, nz, nWidth, ppuch3Dadd, rgb.m_R);
            Set3DPixelB(nx, ny, nz, nWidth, ppuch3Dadd, rgb.m_B);
            Set3DPixelG(nx, ny, nz, nWidth, ppuch3Dadd, rgb.m_G);
        }
        
    }

}

void Save3DImages(std::string strPre,int nLayers, cv::Mat* matSequence)
{
    char chLayer[25];
    std::string strFull;
    for(int i=1; i<=nLayers; i++)
    {
        sprintf(chLayer, "%d.tif", i);
        strFull=strPre + chLayer;
        cv::imwrite(strFull, matSequence[i-1]);
    }
}

void Convert3DGrayToColor(int nLayers, cv::Mat* matSequence, cv::Mat* matColorSeq,unsigned char** puchColor3D)
{
    for(int i=1; i<=nLayers; i++)
    {

        cv::Mat matTempt;
        cv::cvtColor(matSequence[i-1], matTempt, CV_GRAY2BGR);
        matColorSeq[i-1] = matTempt.clone();
        puchColor3D[i-1] = matColorSeq[i-1].data;
    }
}

void Save3DGrayImages(std::string strPre,int nLayers, cv::Mat* matSequence)
{
    char chLayer[25];
    std::string strFull;
    for(int i=1; i<=nLayers; i++)
    {
        sprintf(chLayer, "%d.tif", i);
        strFull=strPre + chLayer;
        cv::Mat matTempt;
        cv::cvtColor(matSequence[i-1], matTempt, CV_GRAY2BGR);
        cv::imwrite(strFull, matTempt);
    }
}

void Test3DFramework()
{
    std::cout<<"syn!\n";
    std::map<CMark,std::vector<TPoint> > CylinderPool;
    std::map<CMark,PreCMark> PredictPool;
    SetParam3D param3D;
    time_t t;
    srandom2(((unsigned)time(&t)));
    /*
    param3D.nMinH=2;
    param3D.nMaxH=8;
    param3D.nMinR=4;
    param3D.nMaxR=5;
    param3D.nDegree=11;
    param3D.nMDegree=10;*/

    ReadCylinderPool(CylinderPool,&PredictPool);
    
    
    std::cout<<"load complete!\n";
    cv::Mat matSequence[MAXSequence];
    cv::Mat matColorSeq[MAXSequence];
    unsigned char* ppuch3Dadd[MAXSequence];
    unsigned char* ppuch3DRGB[MAXSequence];
    std::string strPre = "./CuttedSlices/sliced";
    
    
    CreateVoid3DGrayImages(99, 128, 128,matSequence , ppuch3Dadd);
    Convert3DGrayToColor(99, matSequence, matColorSeq, ppuch3DRGB);
    
    
    int nWidth=128;
    int nHeight =128;
    int nLayers =99;
   // return;
    
    for(int nty=0;nty<90;nty+=3)
        for(int ntz=0;ntz<90;ntz+=3)
        {
    CMark cmark(3,4,nty,ntz);
    cmark.m_nH=4;
    Que3D qStart(50,50,50, cmark, PredictPool[cmark]);
    Que3D* pHeader,*pTail;
    pHeader = &qStart;
    pTail = &qStart;
    Que3D* pPrevious=pHeader;
    
    do
    {
        pPrevious=pHeader;
        pHeader = GrowOneSideT(pHeader, ppuch3Dadd, nWidth , nHeight, nLayers, CylinderPool, &PredictPool, true);
    }
    while(pPrevious!=pHeader);
    
    do
    {
        pPrevious=pTail;
        pTail = GrowOneSideT(pTail, ppuch3Dadd, nWidth , nHeight, nLayers, CylinderPool, &PredictPool, false);
    }
    while(pPrevious!=pTail);
    
    DrawCylHDirectColor(pTail, ppuch3DRGB , 128, 128, 99,CylinderPool,&PredictPool,sRGB(55+200*random2(),55+200*random2(),55+200*random2()));
        }
    Save3DImages("./ColorSeqence/cseq", 99, matColorSeq);
    std::cout<<"Hello world!\n";
   // Save3DGrayImages("./CreateSequence/seq", 99, matSequence);
 //   Load3DGrayImages(strPre,99, matSequence, ppuch3Dadd);
    
}
Que3D* CreateStraighTubes(Que3D* pStartQue, int nNumTubes, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{

    
    Que3D* pHeader = pStartQue;
    int nCount =0;
    int i;
    for(i=0; i<nNumTubes;i++)
    {
        Que3D* pNewQue = new Que3D;
        int nNewTubeH = pHeader->m_mark.m_nH;
        int nNewTubeR = pHeader->m_mark.m_nR;
        int nNewThetaY = pHeader->m_mark.m_nthetaY;
        int nNewThetaZ = pHeader->m_mark.m_nthetaZ;
        
        PreCMark temptPre = pHeader->m_Pmark;
        CMark newCmark = CMark(nNewTubeR, nNewTubeH,nNewThetaY,nNewThetaZ);
        PreCMark premarkNew = (*pPredPool)[newCmark];
        TDPoint pNewObjeCent;
        
        if(pHeader->m_bBackFil&&pHeader->m_bFrontFil)
        {
            std::cout<<"Exit connected\n";
            break;
        }
  
        if(pHeader->m_bFrontFil==false)
        {
           
            TDPoint pFrontP;
            pFrontP.m_dX = temptPre.m_dFCenX;
            pFrontP.m_dY = temptPre.m_dFCenY;
            pFrontP.m_dZ = temptPre.m_dFCenZ;
            if(SelectDirect(pFrontP, premarkNew))
            {
                pNewObjeCent.m_dX = pFrontP.m_dX+pHeader->m_dX+ premarkNew.m_dFCenX;
                pNewObjeCent.m_dY = pFrontP.m_dY+pHeader->m_dY+ premarkNew.m_dFCenY;
                pNewObjeCent.m_dZ = pFrontP.m_dZ+pHeader->m_dZ+ premarkNew.m_dFCenZ;
                pNewQue->m_bBackFil = true;
            }
            else
            {
                pNewObjeCent.m_dX = pFrontP.m_dX+pHeader->m_dX+ premarkNew.m_dBCenX;
                pNewObjeCent.m_dY = pFrontP.m_dY+pHeader->m_dY+ premarkNew.m_dBCenY;
                pNewObjeCent.m_dZ = pFrontP.m_dZ+pHeader->m_dZ+ premarkNew.m_dBCenZ;
                pNewQue->m_bFrontFil = true;
                
            }
            
            pHeader->m_bFrontFil= true;
        }
        else
        {
            TDPoint pBackP;
            pBackP.m_dX = temptPre.m_dBCenX;
            pBackP.m_dY = temptPre.m_dBCenY;
            pBackP.m_dZ = temptPre.m_dBCenZ;
            if(SelectDirect(pBackP, premarkNew))
            {
                pNewObjeCent.m_dX = pBackP.m_dX+pHeader->m_dX+ premarkNew.m_dFCenX;
                pNewObjeCent.m_dY = pBackP.m_dY+pHeader->m_dY+ premarkNew.m_dFCenY;
                pNewObjeCent.m_dZ = pBackP.m_dZ+pHeader->m_dZ+ premarkNew.m_dFCenZ;
                pNewQue->m_bBackFil = true;
            }
            else
            {
                pNewObjeCent.m_dX = pBackP.m_dX+pHeader->m_dX+ premarkNew.m_dBCenX;
                pNewObjeCent.m_dY = pBackP.m_dY+pHeader->m_dY+ premarkNew.m_dBCenY;
                pNewObjeCent.m_dZ = pBackP.m_dZ+pHeader->m_dZ+ premarkNew.m_dBCenZ;
                pNewQue->m_bFrontFil = true;
                
            }
            
            pHeader->m_bBackFil= true;
            
        }
        
        if(!IsInnerPoint(pNewObjeCent.m_dX, pNewObjeCent.m_dY,pNewObjeCent.m_dZ, nWidth, nHeight, nLayers))
        {
            std::cout<<"\n exit out border\n";
            break;
        }
        
        pNewQue->m_mark = newCmark;
        pNewQue->m_Pmark = premarkNew;
        pNewQue->m_dX = pNewObjeCent.m_dX;
        pNewQue->m_dY = pNewObjeCent.m_dY;
        pNewQue->m_dZ = pNewObjeCent.m_dZ;
        pNewQue->m_pQueTail =pHeader;
        pNewQue->m_pQueHead = NULL;
        
        pHeader->m_pQueHead = pNewQue;
        pHeader = pNewQue;
        nCount++;
   
        
    }
    
    
    std::cout<<"created "<<nCount<<" objects in que!\n";
    return pHeader;
}



bool SelectDirect(TDPoint& OneSideCent, PreCMark& newPreMark)
{
    TDPoint TPFront, TPBack;
    
    TPFront.m_dX = newPreMark.m_dFCenX + OneSideCent.m_dX;
    TPFront.m_dY = newPreMark.m_dFCenY + OneSideCent.m_dY;
    TPFront.m_dZ = newPreMark.m_dFCenZ + OneSideCent.m_dZ;
    
    TPBack.m_dX = newPreMark.m_dBCenX + OneSideCent.m_dX;
    TPBack.m_dY = newPreMark.m_dBCenY + OneSideCent.m_dY;
    TPBack.m_dZ = newPreMark.m_dBCenZ + OneSideCent.m_dZ;
    
    if(SqDist2Points3D(TPFront.m_dX,TPFront.m_dY, TPFront.m_dZ, 0, 0, 0)>=SqDist2Points3D(TPBack.m_dX,TPBack.m_dY, TPBack.m_dZ, 0, 0, 0))
    {
        return true;
    }
    else
        return false;
}


int DrawCylHDirectColor_Modified(double *R,double *H,double *Ty,double *Tz,double *fx,double *fy,double *fz, Que3D* pStartQue, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,int id,int nBreak)
{
    Que3D* pHeader = pStartQue;
    int nCount =0;
    int nDegreeY =0;
    int nDegreeZ =0;
    int nR=0;
    int nH=0;
    double dSx = pHeader->m_dX;
    double dSy = pHeader->m_dY;
    double dSz = pHeader->m_dZ;
    double dEx=0,dEy=0,dEz=0;
    double dETy = 0, dETz =0;
    nDegreeY=pHeader->m_mark.m_nthetaY;
    nDegreeZ=pHeader->m_mark.m_nthetaZ;
    int iBreakCnt=0;

    char r=0,g=0,b=0;

    r = 0;
    g = floor(id / 256);
    b = id%256;
    sRGB rgb = sRGB(r,g,b);
    id += 1;
    while(pHeader!=NULL)
    {
      //  std::cout<<"The "<< nCount<<"th "<<"X="<<pHeader->m_dX-50<<"  ";
      //  std::cout<<"Y="<<pHeader->m_dY-50<<" Z="<<pHeader->m_dZ-50<<'\n';
       iBreakCnt++;
       if(nBreak!=iBreakCnt)
            DrawCylinderColor(round(pHeader->m_dX), round(pHeader->m_dY),round(pHeader->m_dZ), pHeader->m_mark,ppuch3Dadd,nWidth, nHeight, nLayers,CylinderPool,pPredPool,rgb);
        nCount++;
        nH+=pHeader->m_mark.m_nH;
        nR+=pHeader->m_mark.m_nR;
        
        dEx = pHeader->m_dX;
        dEy = pHeader->m_dY;
        dEz = pHeader->m_dZ;

        dETy += pHeader->m_mark.m_nthetaY;
        dETz += pHeader->m_mark.m_nthetaZ;
        
        
        if(pHeader->m_pQueHead!=NULL)
        {
            
           pHeader = pHeader->m_pQueHead;

        }
        else{
            break;
        }
    }
    
    if(nCount==1)
    {
        DrawCylinderColor(round(pHeader->m_dX), round(pHeader->m_dY),round(pHeader->m_dZ), pHeader->m_mark,ppuch3Dadd,nWidth, nHeight, nLayers,CylinderPool,pPredPool,rgb);
    }
    double dMx,dMy,dMz;
    double dMTy, dMTz;
    double dMnH;
    //length
    dMnH = nH;
    // mean radius
    nR = nR/nCount;

    // mean angles
    dMTy = dETy/nCount;
    dMTz = dETz/nCount;

    // mean x,y,z coordinates
    dMx= (dSx+dEx)/2;
    dMy= (dSy+dEy)/2;
    dMz= (dSz+dEz)/2;

    *R = nR;
    *H = dMnH;
    *Ty = dMTy;
    *Tz = dMTz;
    *fx = dMx;
    *fy = dMy;
    *fz = dMz;

    return id;
    
}

void DrawCylHDirectColor(Que3D* pStartQue, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb,int nBreak)
{
    Que3D* pHeader = pStartQue;
    int nCount =0;
    int nDegreeY =0;
    int nDegreeZ =0;
    int nR=0;
    int nH=0;
    double dSx = pHeader->m_dX;
    double dSy = pHeader->m_dY;
    double dSz = pHeader->m_dZ;
    double dEx=0,dEy=0,dEz=0;
    double dETy = 0, dETz =0;
    nDegreeY=pHeader->m_mark.m_nthetaY;
    nDegreeZ=pHeader->m_mark.m_nthetaZ;
    int iBreakCnt=0;
    while(pHeader!=NULL)
    {
      //  std::cout<<"The "<< nCount<<"th "<<"X="<<pHeader->m_dX-50<<"  ";
      //  std::cout<<"Y="<<pHeader->m_dY-50<<" Z="<<pHeader->m_dZ-50<<'\n';
       iBreakCnt++;
       if(nBreak!=iBreakCnt)
            DrawCylinderColor(round(pHeader->m_dX), round(pHeader->m_dY),round(pHeader->m_dZ), pHeader->m_mark,ppuch3Dadd,nWidth, nHeight, nLayers,CylinderPool,pPredPool,rgb);
        nCount++;
        nH+=pHeader->m_mark.m_nH;
        nR+=pHeader->m_mark.m_nR;
        
        dEx = pHeader->m_dX;
        dEy = pHeader->m_dY;
        dEz = pHeader->m_dZ;

        dETy += pHeader->m_mark.m_nthetaY;
        dETz += pHeader->m_mark.m_nthetaZ;
        
        
        if(pHeader->m_pQueHead!=NULL)
        {
            
           pHeader = pHeader->m_pQueHead;

        }
        else{
            break;
        }
    }
    
    if(nCount==1)
    {
        DrawCylinderColor(round(pHeader->m_dX), round(pHeader->m_dY),round(pHeader->m_dZ), pHeader->m_mark,ppuch3Dadd,nWidth, nHeight, nLayers,CylinderPool,pPredPool,sRGB(nDegreeY,nDegreeZ,0));
    }
    double dMx,dMy,dMz;
    double dMTy, dMTz;
    double dMnH;
    //length
    dMnH = nH;
    // mean radius
    nR = nR/nCount;

    // mean angles
    dMTy = dETy/nCount;
    dMTz = dETz/nCount;

    // mean x,y,z coordinates
    dMx= (dSx+dEx)/2;
    dMy= (dSy+dEy)/2;
    dMz= (dSz+dEz)/2;    
    
}
Que3D* GrowOneSide(Que3D* pStartQue, unsigned char** ppuchGray,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bIsHead)
{
   // std::cout<<"begin Grow!"<<'\n';
    Que3D* pPointer = pStartQue;
    
    if(pStartQue->m_bBackFil&&pStartQue->m_bFrontFil)
    {
     //   std::cout<<"tail error!";
        return pPointer;
    }
    pPointer->m_Pmark = (*pPredPool)[pPointer->m_mark];
    
    
    /**Grow new pointer this place should be replaced by fitting & searching for best new tube**/

    Que3D* pNewQue = new Que3D;
    Que3D best3D;
    best3D.m_dDateEng=100;
    for(int nTY=-4; nTY<=4; nTY++)
        for(int nTZ=-4;nTZ<=4; nTZ++)
            for(int nH =-1; nH<=2; nH++)
    {
    
    int nNewTubeH  = pPointer->m_mark.m_nH+nH;
    int nNewTubeR  = pPointer->m_mark.m_nR;

    
    if(nNewTubeH<NMINH||nNewTubeH>NMAXH||
       nNewTubeR<NMINR||nNewTubeR>NMAXR)
        continue;
        
        

        int nNewThetaY = pPointer->m_mark.m_nthetaY;
        int nNewThetaZ = pPointer->m_mark.m_nthetaZ+nTZ;
        
        if(nNewThetaZ<0)
        {
            nNewThetaZ+=DEGREE_Z;
            nNewThetaY= DEGREE_Y - nNewThetaY;
        }
        else if(nNewThetaZ>=DEGREE_Z)
        {
            nNewThetaZ-=DEGREE_Z;
            nNewThetaY= DEGREE_Y - nNewThetaY;
        }
        
        nNewThetaY+=nTY;
        
        if(nNewThetaY<0)
            nNewThetaY+=DEGREE_Y;
        else if(nNewThetaY>=DEGREE_Y)
            nNewThetaY-=DEGREE_Y;
        
        
    PreCMark temptPre = pPointer->m_Pmark;
    CMark newCmark = CMark(nNewTubeR, nNewTubeH,nNewThetaY,nNewThetaZ);
    PreCMark premarkNew = (*pPredPool)[newCmark];
    TDPoint pNewObjeCent;
   
    if(pPointer->m_bFrontFil==false)
    {
        
        TDPoint pFrontP;
        pFrontP.m_dX = temptPre.m_dFCenX;
        pFrontP.m_dY = temptPre.m_dFCenY;
        pFrontP.m_dZ = temptPre.m_dFCenZ;
        
        if(SelectDirect(pFrontP, premarkNew))
        {
            pNewObjeCent.m_dX = pFrontP.m_dX+pPointer->m_dX+ premarkNew.m_dFCenX;
            pNewObjeCent.m_dY = pFrontP.m_dY+pPointer->m_dY+ premarkNew.m_dFCenY;
            pNewObjeCent.m_dZ = pFrontP.m_dZ+pPointer->m_dZ+ premarkNew.m_dFCenZ;
            pNewQue->m_bBackFil = true;
            pNewQue->m_bFrontFil = false;
        }
        else
        {
            pNewObjeCent.m_dX = pFrontP.m_dX+pPointer->m_dX+ premarkNew.m_dBCenX;
            pNewObjeCent.m_dY = pFrontP.m_dY+pPointer->m_dY+ premarkNew.m_dBCenY;
            pNewObjeCent.m_dZ = pFrontP.m_dZ+pPointer->m_dZ+ premarkNew.m_dBCenZ;
            pNewQue->m_bFrontFil = true;
            pNewQue->m_bBackFil = false;
            
        }
        
      //  pPointer->m_bFrontFil= true;
    }
    else
    {
        TDPoint pBackP;
        pBackP.m_dX = temptPre.m_dBCenX;
        pBackP.m_dY = temptPre.m_dBCenY;
        pBackP.m_dZ = temptPre.m_dBCenZ;
        if(SelectDirect(pBackP, premarkNew))
        {
            pNewObjeCent.m_dX = pBackP.m_dX+pPointer->m_dX+ premarkNew.m_dFCenX;
            pNewObjeCent.m_dY = pBackP.m_dY+pPointer->m_dY+ premarkNew.m_dFCenY;
            pNewObjeCent.m_dZ = pBackP.m_dZ+pPointer->m_dZ+ premarkNew.m_dFCenZ;
            pNewQue->m_bBackFil = true;
            pNewQue->m_bFrontFil = false;
            
        }
        else
        {
            pNewObjeCent.m_dX = pBackP.m_dX+pPointer->m_dX+ premarkNew.m_dBCenX;
            pNewObjeCent.m_dY = pBackP.m_dY+pPointer->m_dY+ premarkNew.m_dBCenY;
            pNewObjeCent.m_dZ = pBackP.m_dZ+pPointer->m_dZ+ premarkNew.m_dBCenZ;
            pNewQue->m_bFrontFil = true;
            pNewQue->m_bBackFil = false;
        }
        
      //  pPointer->m_bBackFil= true;
        
    }
    if(!IsInnerPoint(pNewObjeCent.m_dX+0.5, pNewObjeCent.m_dY+0.5,pNewObjeCent.m_dZ+0.5, nWidth, nHeight, nLayers))
    {
        std::cout<<temptPre.m_dBCenX<<"  "<<temptPre.m_dBCenY<<"  "<<temptPre.m_dBCenZ<<'\n';
        std::cout<<temptPre.m_dFCenX<<"  "<<temptPre.m_dFCenY<<"  "<<temptPre.m_dFCenZ<<'\n';
        std::cout<<pNewObjeCent.m_dX<<"  "<<pNewObjeCent.m_dY<<"  "<<pNewObjeCent.m_dZ<<'\n';
        std::cout<<"board erro!\n";
        continue;
    }
    
    pNewQue->m_mark = newCmark;
    pNewQue->m_Pmark = premarkNew;
    pNewQue->m_dX = pNewObjeCent.m_dX;
    pNewQue->m_dY = pNewObjeCent.m_dY;
    pNewQue->m_dZ = pNewObjeCent.m_dZ;
    pNewQue->m_dDateEng = CalculateEnergy(pNewQue, ppuchGray, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    pNewQue->m_dOverlapR = OverlapRatio3D(pNewQue, ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    if(pNewQue->m_dOverlapR<OVERLAP_T&&pNewQue->m_dDateEng<best3D.m_dDateEng)
    {
        best3D =*pNewQue;
    }
        
    }  //for for nH,nR,nTY,nTZ
    *pNewQue = best3D;
    
   // AdaptOneLocal3D(pNewQue, ppuchGray,ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool,2);
    
    if(pPointer->m_bFrontFil==false)
        pPointer->m_bFrontFil = true;
    else
        pPointer->m_bBackFil = true;
    
    
    if(!IsInnerPoint(pNewQue->m_dX+0.5, pNewQue->m_dY+0.5,pNewQue->m_dZ+0.5, nWidth, nHeight, nLayers))
    {
        delete pNewQue;
        return pPointer;
    }
    
    
    
    if(pNewQue->m_dOverlapR>OVERLAP_T||pNewQue->m_dDateEng>CYLINDER_ENG_T)
    {
        /*
        if(bIsHead==false)
            std::cout<<"back ";
        else
            std::cout<<"front";
        std::cout<<" ty= "<<nNewThetaY<<" tz= "<<nNewThetaZ<<"overlap error!";*/
        
      if(pNewQue->m_dOverlapR>OVERLAP_T)
         std::cout<<"overlap error!";
        else
            std::cout<<"en error!";
        delete pNewQue;
        return pPointer;
    }
    
  //  std::cout<<"overlap="<<pNewQue->m_dOverlapR<<'\n';
    
    DrawCylinder(pNewQue->m_dX+0.5, pNewQue->m_dY+0.5, pNewQue->m_dZ+0.5, pNewQue->m_mark, ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    
    
        if(bIsHead)
        {
            pNewQue->m_pQueTail =pPointer;
            pNewQue->m_pQueHead = NULL;
            
            pPointer->m_pQueHead = pNewQue;
            pPointer = pNewQue;
        }
        else
        {
            pNewQue->m_pQueHead = pPointer;
            pNewQue->m_pQueTail = NULL;
            
            pPointer->m_pQueTail = pNewQue;
            pPointer = pNewQue;
        }
    
    return pPointer;
}
double OverlapRatio3D(Que3D* pObject, int** ppn3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    double dRatio =0;
    int nCntInnerPts=0;
    int nCntOverlapPts=0;
    std::vector<TPoint> CylinderPoints = CylinderPool[pObject->m_mark];
    int nCenX = pObject->m_dX+0.5;
    int nCenY = pObject->m_dY+0.5;
    int nCenZ = pObject->m_dZ+0.5;
    
    //  std::cout<<"x="<<nCenX<<"y="<<nCenY<<"z="<<nCenZ<<'\n';
    for(int i=0; i<CylinderPoints.size(); i++)
    {
        int nx = CylinderPoints[i].m_nX+nCenX;
        int ny = CylinderPoints[i].m_nY+nCenY;
        int nz = CylinderPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        if(CylinderPoints[i].m_nType==1)
        {
            
            nCntInnerPts++;
            
            if(Get3DPixel(nx,ny,nz,nWidth,ppn3Dadd)>1)
            {
                nCntOverlapPts++;
            }
        }
        
    }
    // std::cout<<"Inp="<<nCntInnerPts<<"overl="<<nCntOverlapPts<<'\n';
    dRatio = double(nCntOverlapPts)/nCntInnerPts;
    return dRatio;
    
}
double OverlapRatio3D(Que3D* pObject, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
    double dRatio =0;
    int nCntInnerPts=0;
    int nCntOverlapPts=0;
    std::vector<TPoint> CylinderPoints = CylinderPool[pObject->m_mark];
    int nCenX = pObject->m_dX+0.5;
    int nCenY = pObject->m_dY+0.5;
    int nCenZ = pObject->m_dZ+0.5;
    
  //  std::cout<<"x="<<nCenX<<"y="<<nCenY<<"z="<<nCenZ<<'\n';
    for(int i=0; i<CylinderPoints.size(); i++)
    {
        int nx = CylinderPoints[i].m_nX+nCenX;
        int ny = CylinderPoints[i].m_nY+nCenY;
        int nz = CylinderPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        if(CylinderPoints[i].m_nType==1)
        {
          
            nCntInnerPts++;
            
           if(Get3DPixel(nx,ny,nz,nWidth,ppuch3Dadd)!=0)
            {
                nCntOverlapPts++;
            }
        }
        
    }
   // std::cout<<"Inp="<<nCntInnerPts<<"overl="<<nCntOverlapPts<<'\n';
    dRatio = double(nCntOverlapPts)/nCntInnerPts;
    return dRatio;
    
}

double CalculateEnergy(Que3D* pObject, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    double dLikelihood=0;
    
    
    int nCntInnerB=0, nCntOuterB =0 ;
  
    double dSumPIn = 0, dSumPInP2=0;
    
    double dSumPOut=0, dSumPOutP2=0;
    
    double dInnerU =0,dOuterU=0;
    int nTempt;
    
    std::vector<TPoint> vTPoints = CylinderPool[pObject->m_mark];
    
    int nCenX = pObject->m_dX+0.5;
    int nCenY = pObject->m_dY+0.5;
    int nCenZ = pObject->m_dZ+0.5;
    
    
    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        
        if(vTPoints[i].m_nType==1)
        {
            nCntInnerB++;
            
            nTempt = Get3DPixel(nx, ny, nz, nWidth, ppuch3Dadd);
            dSumPIn+=nTempt;
            dSumPInP2 +=nTempt*nTempt;

            
        }
        else if(vTPoints[i].m_nType==2)
        {
            nCntOuterB++;
            
            nTempt = Get3DPixel(nx, ny, nz, nWidth, ppuch3Dadd);
            dSumPOut +=nTempt;
            dSumPOutP2 += nTempt*nTempt;
        }
        
    }
    
    dInnerU = (dSumPIn)/(double)nCntInnerB;
    dOuterU =(dSumPOut)/(double)nCntOuterB;
    double dDifInOut = ABS(dInnerU-dOuterU);
    
    
    
    if(dInnerU<INNERU_T||dDifInOut<DIFU_T||dInnerU<dOuterU||dOuterU>OUTER_T)
    {
        dLikelihood =1;
        return dLikelihood;
    }
    
    double dInnerCov = dSumPInP2/(double)nCntInnerB-dInnerU*dInnerU;
    double dOuterCov = dSumPOutP2/(double)nCntOuterB-dOuterU*dOuterU;
    
    double dBhattacharya = sqrt(dDifInOut)+dDifInOut*dDifInOut/(sqrt(dInnerCov+dOuterCov))-INNCOV*sqrt(dInnerCov);
    double dthredh = DTHRED;
    if (dBhattacharya<dthredh)
        dLikelihood = 1-dBhattacharya/dthredh;
    else
        dLikelihood = exp(-(dBhattacharya-dthredh)/(3*dBhattacharya))-1;
    
    return dLikelihood;
    
}




void DrawFiberResults3D(std::vector<Que3D*>& vFiber,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
   // int nTotalFibers = vFiber.size();
    
    std::cout<<"Total Fiber is "<<vFiber.size()<<'\n';
    for(int i=0; i<vFiber.size(); i++)
    {
        sRGB rgb(random2()*255+80,random2()*255+80,random2()*255+80);
        DrawCylHDirectColor(vFiber[i], ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool, rgb);
        
    }
    
    
    
}

void TestFiberDetection(int SubVn)
{
    // Camilo
    //Define the 3 length for 3 dimensional
    int nWidth=450;//301;
    int nHeight =450;//301;
    int nLayers =450;//512;//301;
    
    std::map<CMark,std::vector<TPoint> > CylinderPool;     //The pool store the cylinder information
    std::map<CMark,PreCMark> PredictPool;                  //The pool for predicting
    SetParam3D param3D;                                    //Parameter for the MPP process




    ReadCylinderPool(CylinderPool,&PredictPool);          //Load the CylinderPool and Predict Pool from files "CylinderPool.data"
      
    std::cout<<"Finish reading Pool\n";


    cv::Mat matSequence[MAXSequence];                     //Define the original 3D sequence space 
    cv::Mat matColorSeq[MAXSequence];                     //Define 3D sequence as same size as the original data, but with RGB value, for storage use
    cv::Mat matSeg[MAXSequence];                          //Define gray 3D sequence for given birth, it's a birth map 

    //Define the unsigned char pointer for locate the images sequences define above 
    unsigned char* ppuch3DGray[MAXSequence];
    unsigned char* ppuch3DRGB[MAXSequence];
    unsigned char* ppuch3DSeg[MAXSequence];
    std::string strPre = "./SUBVOLUMES/sV" + std::to_string(SubVn) + "/data/subV_"; //./originals/im_num";//"./JPG/sliced_";//"./SynSeq/seq";"//./CuttedSlices/sliced";
    std::string strPreSeg = "./SUBVOLUMES/sV" + std::to_string(SubVn) + "/seg/subV_seg_";//./originals/im_num";//"./JPG/sliced_";//"./SynSeq/seq";"//./CuttedSlices/sliced";
    //std::string strPreSeg2 =  "./CuttedSeg/seg";
    Load3DGrayImages(strPre,nLayers, matSequence, ppuch3DGray);   //Load the original 3D data from file strPre
    
    CreateVoid3DGrayImages(nLayers, nWidth, nHeight, matSeg,ppuch3DSeg);
    
    Load3DGrayImages(strPreSeg,nLayers, matSeg, ppuch3DSeg);   //Load the preseg 3D data from file strPre
    int threshold = (int)(255.0 / 6.0 * 5.0); //PRE-SEGMENTATION EM/MPM with 6 casses. Assume class 5 and class 6 represent fiber
    GetCoarseSeg3DImages(nLayers, matSeg, matSeg, threshold);
    //Save3DImages(strPreSeg2, nLayers, matSeg);
    

    std::cout<<"Load Image and Pools finished! \n";
    std::cout<<"Size CylinderPool="<<CylinderPool.size()<<'\n';  
    //MPP ing
    std::vector<Que3D*> vFiber;

    std::cout<<"Fiber Detecting ... \n";    

    //To Camilo: the parallelization should happen in the function CylinderMPPBirthDeath, which is the core function for MPP process
    CylinderMPPBirthDeath(vFiber, ppuch3DGray, ppuch3DSeg, nWidth, nHeight, nLayers, CylinderPool, &PredictPool, SubVn);
    
    
    return;
    


}

void AdaptOneLocal3D(Que3D* qOneObj,unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,int nRange)
{
    Que3D qTempt = *qOneObj;
    Que3D qBest =  *qOneObj;
    for(int nx=-nRange; nx<=nRange; nx++)
        for(int ny=-nRange; ny<=nRange; ny++)
            for(int nz=-nRange; nz<=nRange; nz++)
                for(int nTy=-3; nTy<=3; nTy++)
                    for(int nTz=-3; nTz<=3; nTz++)
                        for(int nR=-2; nR<=2; nR++)
                            for(int nH=-1;nH<=2;nH++)
                        {
                            
                          qTempt.m_mark.m_nH = qOneObj->m_mark.m_nH+nH;
                            
                          if(qTempt.m_mark.m_nH<NMINH||qTempt.m_mark.m_nH>NMAXH)
                              continue;
                            
                          qTempt.m_dX = nx+qOneObj->m_dX;
                          qTempt.m_dY = ny+qOneObj->m_dY;
                          qTempt.m_dZ = nz+qOneObj->m_dZ;
                          
                          if(!IsInnerPoint(int((qTempt).m_dX+0.5), int((qTempt).m_dY+0.5), int((qTempt).m_dZ+0.5), nWidth, nHeight, nLayers))
                          {
                              continue;
                          }
                          
                          qTempt.m_mark.m_nthetaY = qOneObj->m_mark.m_nthetaY+nTy;
                          qTempt.m_mark.m_nthetaZ = qOneObj->m_mark.m_nthetaZ+nTz;
                            
                        
                         if(qTempt.m_mark.m_nthetaY<0)
                              qTempt.m_mark.m_nthetaY+=DEGREE_Y;
                          else if(qTempt.m_mark.m_nthetaY>=DEGREE_Y)
                              qTempt.m_mark.m_nthetaY-=DEGREE_Y;
                            
                         if(qTempt.m_mark.m_nthetaZ<0)
                             qTempt.m_mark.m_nthetaZ+=DEGREE_Z;
                            else if(qTempt.m_mark.m_nthetaZ>=DEGREE_Z)
                                qTempt.m_mark.m_nthetaZ-=DEGREE_Z;
                            
                            
                        
                          qTempt.m_mark.m_nR = qOneObj->m_mark.m_nR+nR;
                         
                          if(qTempt.m_mark.m_nR<NMINR||NMINR>NMAXR)
                          {
                                 continue;
                          }
                            
                          qTempt.m_dOverlapR = OverlapRatio3D(&qTempt, ppuchSeg3D, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
                          if(qTempt.m_dOverlapR>OVERLAP_T)
                              continue;
                        
                          qTempt.m_dDateEng = CalculateEnergy(&qTempt, ppuch3DGray, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
                            
                          if(qTempt.m_dDateEng<qBest.m_dDateEng)
                          {
                              qBest = qTempt;
                          }
                        
                        
                       }
    
    
    *qOneObj = qBest;
    
    
}

Que3D* GrowOneSideT(Que3D* pStartQue,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bIsHead)
{
    Que3D* pPointer = pStartQue;
    
    if(pStartQue->m_bBackFil&&pStartQue->m_bFrontFil)
    {
       // std::cout<<"tail error!";
        return pPointer;
    }
    
    
    /**Grow new pointer this place should be replaced by fitting & searching for best new tube**/
    
    Que3D* pNewQue = new Que3D;
    
    
   
    int nNewTubeH  = pPointer->m_mark.m_nH;
    int nNewTubeR  = pPointer->m_mark.m_nR;
    int nNewThetaY = pPointer->m_mark.m_nthetaY;
    int nNewThetaZ = pPointer->m_mark.m_nthetaZ;
                    

    PreCMark temptPre = pPointer->m_Pmark;
    CMark newCmark = CMark(nNewTubeR, nNewTubeH,nNewThetaY,nNewThetaZ);
    PreCMark premarkNew = (*pPredPool)[newCmark];
    TDPoint pNewObjeCent;
                    
    if(pPointer->m_bFrontFil==false)
    {
                        
        TDPoint pFrontP;
        pFrontP.m_dX = temptPre.m_dFCenX;
        pFrontP.m_dY = temptPre.m_dFCenY;
        pFrontP.m_dZ = temptPre.m_dFCenZ;
        if(SelectDirect(pFrontP, premarkNew))
        {
            pNewObjeCent.m_dX = pFrontP.m_dX+pPointer->m_dX+ premarkNew.m_dFCenX;
            pNewObjeCent.m_dY = pFrontP.m_dY+pPointer->m_dY+ premarkNew.m_dFCenY;
            pNewObjeCent.m_dZ = pFrontP.m_dZ+pPointer->m_dZ+ premarkNew.m_dFCenZ;
            pNewQue->m_bBackFil = true;
            pNewQue->m_bFrontFil = false;
        }
        else
        {
            pNewObjeCent.m_dX = pFrontP.m_dX+pPointer->m_dX+ premarkNew.m_dBCenX;
            pNewObjeCent.m_dY = pFrontP.m_dY+pPointer->m_dY+ premarkNew.m_dBCenY;
            pNewObjeCent.m_dZ = pFrontP.m_dZ+pPointer->m_dZ+ premarkNew.m_dBCenZ;
            pNewQue->m_bFrontFil = true;
            pNewQue->m_bBackFil = false;
                            
        }
                        
        pPointer->m_bFrontFil= true;
    }
    else
    {
        TDPoint pBackP;
        pBackP.m_dX = temptPre.m_dBCenX;
        pBackP.m_dY = temptPre.m_dBCenY;
        pBackP.m_dZ = temptPre.m_dBCenZ;
        if(SelectDirect(pBackP, premarkNew))
        {
            pNewObjeCent.m_dX = pBackP.m_dX+pPointer->m_dX+ premarkNew.m_dFCenX;
            pNewObjeCent.m_dY = pBackP.m_dY+pPointer->m_dY+ premarkNew.m_dFCenY;
            pNewObjeCent.m_dZ = pBackP.m_dZ+pPointer->m_dZ+ premarkNew.m_dFCenZ;
            pNewQue->m_bBackFil = true;
            pNewQue->m_bFrontFil = false;
                            
        }
        else
        {
            pNewObjeCent.m_dX = pBackP.m_dX+pPointer->m_dX+ premarkNew.m_dBCenX;
            pNewObjeCent.m_dY = pBackP.m_dY+pPointer->m_dY+ premarkNew.m_dBCenY;
            pNewObjeCent.m_dZ = pBackP.m_dZ+pPointer->m_dZ+ premarkNew.m_dBCenZ;
            pNewQue->m_bFrontFil = true;
            pNewQue->m_bBackFil = false;
        }
                        
          pPointer->m_bBackFil= true;
                        
    }
    if(!IsInnerPoint(pNewObjeCent.m_dX+0.5, pNewObjeCent.m_dY+0.5,pNewObjeCent.m_dZ+0.5, nWidth, nHeight, nLayers))
    {
        delete pNewQue;
        return pPointer;
    }
                    
    pNewQue->m_mark = newCmark;
    pNewQue->m_Pmark = premarkNew;
    pNewQue->m_dX = pNewObjeCent.m_dX;
    pNewQue->m_dY = pNewObjeCent.m_dY;
    pNewQue->m_dZ = pNewObjeCent.m_dZ;

    
    DrawCylinder(pNewQue->m_dX+0.5, pNewQue->m_dY+0.5, pNewQue->m_dZ+0.5, pNewQue->m_mark, ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    
    
    
    if(bIsHead)
    {
        pNewQue->m_pQueTail =pPointer;
        pNewQue->m_pQueHead = NULL;
        
        pPointer->m_pQueHead = pNewQue;
        pPointer = pNewQue;
    }
    else
    {
        pNewQue->m_pQueHead = pPointer;
        pNewQue->m_pQueTail = NULL;
        
        pPointer->m_pQueTail = pNewQue;
        pPointer = pNewQue;
    }
    
    return pPointer;
}

void SynCylinderData7(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides)
{
    std::map<CMark,PreCMark> PredictPool = *pPredPool;
    
    CMark cmark(5,4,16,38);//
    cmark.m_nH=2;
    Que3D qStart(62,48,39, cmark, PredictPool[cmark]);
    
    SynCylinder(12, &qStart, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(190,199,212),true,4);
    
    
    cmark = CMark(4,3,89,20);
    Que3D qStart1(30,66,69, cmark, PredictPool[cmark]);
    SynCylinder(15, &qStart1, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true,5);
    
    
    
    cmark = CMark(4,3,62,40);
    Que3D qStart2(45,45,45, cmark, PredictPool[cmark]);
    SynCylinder(15, &qStart2, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(222,229,252),true,6);
    
    
    
    cmark = CMark(4,3,20,70);
    Que3D qStart3(90,8,4, cmark, PredictPool[cmark]);
    SynCylinder(15, &qStart3, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(222,229,252),true);
    
}


void SynCylinderData6(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides)
{
    std::map<CMark,PreCMark> PredictPool = *pPredPool;
    
    CMark cmark(5,4,86,51);//
    cmark.m_nH=4;
    Que3D qStart(50,50,50, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(190,199,212),true);
    
    
    
    
}

void SynCylinderData1(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides)
{
    
    std::map<CMark,PreCMark> PredictPool = *pPredPool;
    
    CMark cmark(5,4,16,38);//
    cmark.m_nH=4;
    Que3D qStart(50,50,50, cmark, PredictPool[cmark]);
    
    SynCylinder(5, &qStart, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(190,199,212),true);
    
    cmark = CMark(4,4,60,40);
    Que3D qStart1(10,20,50, cmark, PredictPool[cmark]);
    SynCylinder(5, &qStart1, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(4,5,3,40);
    Que3D qStart2(30,40,60, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart2, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(5,6,22,23);
    Que3D qStart3(77,20,30, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart3, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
    
    cmark = CMark(5,6,21,23);
    Que3D qStart4(88,20,30, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart4, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
    
    
    
    cmark = CMark(5,6,33,17);
    Que3D qStart5(90,90,90, cmark, PredictPool[cmark]);
    SynCylinder(8, &qStart5, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);

}
void SynCylinderData3(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides)
{
    
    std::map<CMark,PreCMark> PredictPool = *pPredPool;
    
    CMark cmark(4,6,26,80);
    cmark.m_nH=4;
    Que3D qStart(10,20,30, cmark, PredictPool[cmark]);
    
    SynCylinder(8, &qStart, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(200,199,212),true);
    
    cmark = CMark(4,6,26,80);
    Que3D qStart1(10,20,50, cmark, PredictPool[cmark]);
    SynCylinder(5, &qStart1, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,179,232),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart2(10,30,50, cmark, PredictPool[cmark]);
    SynCylinder(8, &qStart2, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart3(10,40,60, cmark, PredictPool[cmark]);
    SynCylinder(8, &qStart3, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(199,179,202),true);
    
    cmark = CMark(4,6,26,80);
    Que3D qStart4(10,50,70, cmark, PredictPool[cmark]);
    SynCylinder(8, &qStart4, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(229,219,202),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart5(10,70,80, cmark, PredictPool[cmark]);
    SynCylinder(8, &qStart5, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,219,202),true);
    
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart11(30,20,50, cmark, PredictPool[cmark]);
    SynCylinder(5, &qStart11, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,179,232),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart12(30,30,50, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart12, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart13(30,40,60, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart13, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(199,179,202),true);
    
    cmark = CMark(4,6,26,80);
    Que3D qStart14(30,50,70, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart14, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(229,219,202),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart15(30,70,80, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart15, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,219,202),true);
    
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart21(50,20,50, cmark, PredictPool[cmark]);
    SynCylinder(5, &qStart21, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,179,232),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart22(50,30,50, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart22, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart23(50,40,60, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart23, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(199,179,202),true);
    
    cmark = CMark(4,6,26,80);
    Que3D qStart24(50,50,70, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart24, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(229,219,202),true);
    
    
    cmark = CMark(4,6,26,80);
    Que3D qStart25(50,70,80, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart25, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,219,202),true);
    
    
    
    cmark = CMark(4,6,70,26);
    Que3D qStart31(50,70,10, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart31, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,242,202),true);
 
    
    
    cmark = CMark(4,6,70,26);
    Que3D qStart32(50,90,10, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart32, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,242,202),true);
    
    
    cmark = CMark(4,6,65,28);
    Que3D qStart33(50,110,15, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart33, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,242,202),true);
    /*
    
    cmark = CMark(4,6,70,26);
    Que3D qStart33(50,70,40, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart33, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,242,202),true);
    
    cmark = CMark(4,6,70,26);
    Que3D qStart34(50,70,55, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart34, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,242,202),true);
    

    cmark = CMark(4,5,3,40);
    Que3D qStart2(30,40,60, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart2, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(5,6,22,23);
    Que3D qStart3(77,20,30, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart3, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
    
    cmark = CMark(5,6,21,23);
    Que3D qStart4(88,20,30, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart4, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
    
    
    
    cmark = CMark(5,6,33,17);
    Que3D qStart5(90,90,90, cmark, PredictPool[cmark]);
    SynCylinder(8, &qStart5, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
    */
}



void SynCylinderData4(int nLen,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides)
{
    
    
    std::map<CMark,PreCMark> PredictPool = *pPredPool;
    
    CMark cmark(5,4,0,12);
    cmark.m_nH=4;
    Que3D qStart(5,5,10, cmark, PredictPool[cmark]);
    
    SynCylinder(8, &qStart, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(245,229,212),true);
    
    cmark = CMark(4,4,2,13);
    Que3D qStart1(6,15,20, cmark, PredictPool[cmark]);
    SynCylinder(5, &qStart1, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(4,5,4,16);
    Que3D qStart2(20,25,30, cmark, PredictPool[cmark]);
    SynCylinder(10, &qStart2, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
    
    
    cmark = CMark(5,6,7,20);
    Que3D qStart3(77,20,30, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart3, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
    
    cmark = CMark(5,6,45,45);
    Que3D qStart4(88,20,30, cmark, PredictPool[cmark]);
    SynCylinder(4, &qStart4, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
     
     cmark = CMark(4,6,40,50);
     Que3D qStart33(50,70,40, cmark, PredictPool[cmark]);
     SynCylinder(10, &qStart33, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,242,202),true);
 
     cmark = CMark(4,6,20,70);
     Que3D qStart34(50,65,55, cmark, PredictPool[cmark]);
     SynCylinder(6, &qStart34, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(209,242,202),true);
     
     
     
    
     
     cmark = CMark(4,5,60,60);
     Que3D qStart22(30,90,90, cmark, PredictPool[cmark]);
     SynCylinder(4, &qStart22, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(210,209,252),true);
     
     /*
     cmark = CMark(5,6,22,23);
     Que3D qStart3(77,20,30, cmark, PredictPool[cmark]);
     SynCylinder(4, &qStart3, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
     
     cmark = CMark(5,6,21,23);
     Que3D qStart4(88,20,30, cmark, PredictPool[cmark]);
     SynCylinder(4, &qStart4, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
     
     
     
     cmark = CMark(5,6,33,17);
     Que3D qStart5(90,90,90, cmark, PredictPool[cmark]);
     SynCylinder(8, &qStart5, cmark,ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(220,249,252),true);
     */
}


void SynCylinder(int nLen, Que3D* qStartObj,CMark cmark,unsigned char** ppuch3DSeg,unsigned char** ppuch3DRGB,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,sRGB rgb, bool bIsbothSides,int nBreak)
{
    Que3D* pHeader,*pTail;
    pHeader = qStartObj;
    pTail = qStartObj;
    Que3D* pPrevious=pHeader;
    
    while(nLen-->0)
    {
        pHeader = GrowOneSideT(pHeader, ppuch3DSeg, nWidth , nHeight, nLayers, CylinderPool, pPredPool, true);
        if(bIsbothSides)
        {
            pTail = GrowOneSideT(pTail, ppuch3DSeg, nWidth , nHeight, nLayers, CylinderPool, pPredPool,false);
        }
    }
    

    DrawCylHDirectColor(pTail, ppuch3DRGB, nWidth , nHeight, nLayers, CylinderPool, pPredPool,rgb,nBreak);
    
}
void SynFibers3DSeq()
{
    
    
    std::map<CMark,std::vector<TPoint> > CylinderPool;
    std::map<CMark,PreCMark> PredictPool;
    SetParam3D param3D;
    
    // CreateCylinderPool(CylinderPool, param3D,&PredictPool);
    // SaveCylinderPool(CylinderPool,&PredictPool);
    
    time_t t;
    srandom2(((unsigned)time(&t)));
    
    
    // return;
    
    std::cout<<"Synthesizing ... \n";
    
    int nWidth = 129;
    int nHeight =129;
    int nLayers = 99;
    
    
    ReadCylinderPool(CylinderPool,&PredictPool);
    cv::Mat matSequence[MAXSequence];
    cv::Mat matColorSeq[MAXSequence];
    cv::Mat matSeg[MAXSequence];
    unsigned char* ppuch3DGray[MAXSequence];
    unsigned char* ppuch3DRGB[MAXSequence];
    unsigned char* ppuch3DSeg[MAXSequence];
    std::string strPre = "./CuttedSlices/sliced";
    
    
    
    CreateVoid3DGrayImages(nLayers, nWidth, nHeight, matSeg , ppuch3DSeg);
    
   
    
    Convert3DGrayToColor(nLayers, matSeg, matColorSeq, ppuch3DRGB);
    GrayAndNoiseBG_GRAY(ppuch3DRGB, nWidth, nHeight, nLayers);
   
    SynCylinderData7(5, ppuch3DSeg,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,&PredictPool,sRGB(190,199,212),true);
    
   
    Save3DImages("./SynSeq/seq", nLayers, matColorSeq);
    std::cout<<"Synthesis Finished!\n";
    
    
}


void GrayAndNoiseBG_GRAY(unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers)
{
  
    for(int z=0; z<nLayers; z++)
        for(int x=0; x<nWidth; x++)
            for(int y=0; y<nHeight; y++)
            {
                unsigned char uchBg =88;
                if(random2()<0.01)
                {
                    uchBg+= int(162*random2());
                }
                else if(random2()<0.02)
                {
                    uchBg-=int(40*random2());
                }
                Set3DPixelR(x, y, z, nWidth, ppuch3Dadd, uchBg);
                Set3DPixelB(x, y, z, nWidth, ppuch3Dadd, uchBg);

                Set3DPixelG(x, y, z, nWidth, ppuch3Dadd, uchBg);

                
            }
    
}
void CalculateEnergy(Que3D* pObject, unsigned char** ppuch3Dadd,int** ppnConn,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,double& dDataE,double& dConnFRatio,double& dConnBRatio, double& dOvelapE)
{
    double dLikelihood=0;
    
    int nCntFConnSum =0, nCntBConnSum=0;
    int nCntFConn =0, nCntBConn=0;
    
    int nCntInnerB=0, nCntOuterB =0 ;
    
    double dSumPIn = 0, dSumPInP2=0;
    
    double dSumPOut=0, dSumPOutP2=0;
    
    double dInnerU =0,dOuterU=0;
    int nTempt;
    int nOverlapCnt =0;
    
    std::vector<TPoint> vTPoints = CylinderPool[pObject->m_mark];
    
    int nCenX = pObject->m_dX+0.5;
    int nCenY = pObject->m_dY+0.5;
    int nCenZ = pObject->m_dZ+0.5;
    
    
    
    
    
    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        
        if(vTPoints[i].m_nType==1)
        {
            nCntInnerB++;
            
            nTempt = Get3DPixel(nx, ny, nz, nWidth, ppuch3Dadd);
            dSumPIn+=nTempt;
            dSumPInP2 +=nTempt*nTempt;
            
            if(Get3DPixel(nx, ny, nz, nWidth, ppnOverlap)>0)
            {
                nOverlapCnt++;
            }
            
        }
        else if(vTPoints[i].m_nType==2)
        {
            nCntOuterB++;
            
            nTempt = Get3DPixel(nx, ny, nz, nWidth, ppuch3Dadd);
            dSumPOut +=nTempt;
            dSumPOutP2 += nTempt*nTempt;
        }
        else if(vTPoints[i].m_nType==3)
        {
            nCntFConnSum++;
            if(Get3DPixel(nx, ny, nz, nWidth, ppnConn)!=0)
            {
                nCntFConn++;
            }
        }
        else if(vTPoints[i].m_nType==4)
        {
            nCntBConnSum++;
            if(Get3DPixel(nx, ny, nz, nWidth, ppnConn)!=0)
            {
                nCntBConn++;
            }
        }
        
    }
    
    dConnFRatio = double(nCntFConn)/nCntFConnSum;   //Connection Ratio
    dConnBRatio = double(nCntBConn)/nCntBConnSum;
    dOvelapE = double(nOverlapCnt)/nCntInnerB;
    
    
    dInnerU = (dSumPIn)/(double)nCntInnerB;
    dOuterU =(dSumPOut)/(double)nCntOuterB;
    double dDifInOut = ABS(dInnerU-dOuterU);
    
    
    
    if(dInnerU<INNERU_T||dDifInOut<DIFU_T||dInnerU<dOuterU)
    {
        dDataE =1;
        return;
    }
    
    double dInnerCov = dSumPInP2/(double)nCntInnerB-dInnerU*dInnerU;
    double dOuterCov = dSumPOutP2/(double)nCntOuterB-dOuterU*dOuterU;
    
    double dBhattacharya = sqrt(dDifInOut)+dDifInOut*dDifInOut/(sqrt(dInnerCov+dOuterCov))-INNCOV*sqrt(dInnerCov);
    
    double dthredh = DTHRED;
    if (dBhattacharya<dthredh)
        dLikelihood = 1-dBhattacharya/dthredh;
    else
        dLikelihood = exp(-(dBhattacharya-dthredh)/(3*dBhattacharya))-1;
    
    dDataE=dLikelihood;

}


void CalculateEnergy(Que3D* pObject, unsigned char** ppuch3Dadd,int** ppnConn,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,double& dDataE,double& dConnFRatio,double& dConnBRatio)
{
    double dLikelihood=0;
    
    int nCntFConnSum =0, nCntBConnSum=0;
    int nCntFConn =0, nCntBConn=0;
    
    int nCntInnerB=0, nCntOuterB =0 ;
    
    double dSumPIn = 0, dSumPInP2=0;
    
    double dSumPOut=0, dSumPOutP2=0;
    
    double dInnerU =0,dOuterU=0;
    int nTempt;
    
    std::vector<TPoint> vTPoints = CylinderPool[pObject->m_mark];
    
    int nCenX = pObject->m_dX+0.5;
    int nCenY = pObject->m_dY+0.5;
    int nCenZ = pObject->m_dZ+0.5;
    
    
    
    
    
    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        
        if(vTPoints[i].m_nType==1)
        {
            nCntInnerB++;
            
            nTempt = Get3DPixel(nx, ny, nz, nWidth, ppuch3Dadd);
            dSumPIn+=nTempt;
            dSumPInP2 +=nTempt*nTempt;
            
            
        }
        else if(vTPoints[i].m_nType==2)
        {
            nCntOuterB++;
            
            nTempt = Get3DPixel(nx, ny, nz, nWidth, ppuch3Dadd);
            dSumPOut +=nTempt;
            dSumPOutP2 += nTempt*nTempt;
        }
        else if(vTPoints[i].m_nType==3)
        {
            nCntFConnSum++;
            if(Get3DPixel(nx, ny, nz, nWidth, ppnConn)>0)
            {
                nCntFConn++;
            }
        }
        else if(vTPoints[i].m_nType==4)
        {
            nCntBConnSum++;
            if(Get3DPixel(nx, ny, nz, nWidth, ppnConn)>0)
            {
                nCntBConn++;
            }
        }
        
    }
    
    dConnFRatio = double(nCntFConn)/nCntFConnSum;   //Connection Ratio
    dConnBRatio = double(nCntBConn)/nCntBConnSum;
    
    
    dInnerU = (dSumPIn)/(double)nCntInnerB;
    dOuterU =(dSumPOut)/(double)nCntOuterB;
    double dDifInOut = ABS(dInnerU-dOuterU);
    
    
    if(dInnerU<INNERU_T||dDifInOut<DIFU_T||dInnerU<dOuterU||dOuterU>OUTER_T)
    {
        dDataE =1;
        return;
    }
    
    double dInnerCov = dSumPInP2/(double)nCntInnerB-dInnerU*dInnerU;
    double dOuterCov = dSumPOutP2/(double)nCntOuterB-dOuterU*dOuterU;
    
    double dBhattacharya = sqrt(dDifInOut)+dDifInOut*dDifInOut/(sqrt(dInnerCov+dOuterCov))-INNCOV*sqrt(dInnerCov);
    
    double dthredh = DTHRED;
    if (dBhattacharya<dthredh)
        dLikelihood = 1-dBhattacharya/dthredh;
    else
        dLikelihood = exp(-(dBhattacharya-dthredh)/(3*dBhattacharya))-1;
    
    dDataE=dLikelihood;
    
}

void AddObjectOverLap_CONN(Que3D* pObject,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
    std::vector<TPoint> vTPoints = CylinderPool[pObject->m_mark];
    
    int nCenX = round(pObject->m_dX);
    int nCenY = round(pObject->m_dY);
    int nCenZ = round(pObject->m_dZ);

    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        if(vTPoints[i].m_nType==1)
        {
       
            AddOne3DP(nx, ny, nz, nWidth,ppnOverlap);
            
        }
        else if(vTPoints[i].m_nType>2)
        {
            AddOne3DP(nx, ny, nz, nWidth,ppn3DConnect);
        }

        
    }
    
    
}

void EraseObjectOverLap_CONN(Que3D* pObject,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
    std::vector<TPoint> vTPoints = CylinderPool[pObject->m_mark];
    
    int nCenX = round(pObject->m_dX);
    int nCenY = round(pObject->m_dY);
    int nCenZ = round(pObject->m_dZ);
    
    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        if(vTPoints[i].m_nType==1)
        {
            
            DedOne3DP(nx, ny, nz, nWidth,ppnOverlap);
            
        }
        else if(vTPoints[i].m_nType>2)
        {
            DedOne3DP(nx, ny, nz, nWidth,ppn3DConnect);
        }
        
        
    }
    
    
}

void UpdatePriorConnect3D(std::vector<Que3D*>& vObjects, int** ppn3DConnect,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool)
{
    
    for(int i=0; i<vObjects.size(); i++)
    {
        std::vector<TPoint> vTPoints = CylinderPool[vObjects[i]->m_mark];
        
        int nCenX = vObjects[i]->m_dX+0.5;
        int nCenY = vObjects[i]->m_dY+0.5;
        int nCenZ = vObjects[i]->m_dZ+0.5;
        int nCntFConnSum =0, nCntBConnSum=0;
        int nCntFConn =0, nCntBConn=0;
        
        for(int j=0; j<vTPoints.size(); j++)
        {
            int nx = vTPoints[j].m_nX+nCenX;
            int ny = vTPoints[j].m_nY+nCenY;
            int nz = vTPoints[j].m_nZ+nCenZ;
            
            if(nx<0||nx>=nWidth||
               ny<0||ny>=nHeight||
               nz<0||nz>=nLayers)
            {
                continue;
            }
            

            
            if(vTPoints[j].m_nType==3)
            {
                nCntFConnSum++;
                if(Get3DPixel(nx, ny, nz, nWidth, ppn3DConnect)>1)
                {
                    nCntFConn++;
                }
            }
            else if(vTPoints[j].m_nType==4)
            {
                nCntBConnSum++;
                if(Get3DPixel(nx, ny, nz, nWidth, ppn3DConnect)>1)
                {
                    nCntBConn++;
                }
            }
        }
        
        double dConnFRatio = double(nCntFConn)/nCntFConnSum;   //Connection Ratio
        double dConnBRatio = double(nCntBConn)/nCntBConnSum;
            
        double dPriorC = (dConnFRatio+dConnBRatio)/2;
        
        if (dPriorC<CONNECTFINEP3D)
                //    dPriorC = -finePos/CONNECTFINEP*dPriorC+finePos;
            dPriorC = CONNECTPRIOR3D*(1-dPriorC/CONNECTFINEP3D);
        else
            //    dPriorC =1.5*(-dPriorC*dPriorC+CONNECTFINEP*CONNECTFINEP)/(1-CONNECTFINEP*CONNECTFINEP);
            dPriorC =-CONNECTPRIOR3D*1.8*(sqrt(dPriorC)-sqrt(CONNECTFINEP3D))/(1-sqrt(CONNECTFINEP3D));
        
        
        if(dConnFRatio>DoubleConnThresh&&dConnBRatio>DoubleConnThresh)
            dPriorC+=DoubleConnBonus;
        
        if(dConnBRatio>=T_LOSSCON_MIN&&dConnBRatio<T_LOSSCON_MAX)
        {
            vObjects[i]->m_bBackM=true;
        }
        else
        {
            vObjects[i]->m_bBackM=false;
        }
        
        if(dConnFRatio>=T_LOSSCON_MIN&&dConnFRatio<T_LOSSCON_MAX)
        {
            vObjects[i]->m_bFrontM=true;
        }
        else
        {
            vObjects[i]->m_bFrontM=false;
        }
        
        vObjects[i]->m_dPriorCon = dPriorC;
        vObjects[i]->m_dFConnR = dConnFRatio;
        vObjects[i]->m_dBConnR = dConnBRatio;
        
    }
}



void CylinderMPPBirthDeath(std::vector<Que3D*>& vFiber,unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, int SubVn)
{
    //ShowSavedFibers(nWidth, nHeight, nLayers,"./SUBVOLUMES/SV7/fibers_info/fiber1.data", CylinderPool, pPredPool);
    //return;
    time_t t;
    srandom2(((unsigned)time(&t)));   //t is defined for random seed

    MppControPara conP;     //this is the birth and death temperature control

    cv::Mat matOverlapS[MAXSequence];                 //The 3D image (same size as the original images) for save the overlapping information
    cv::Mat matConnectS[MAXSequence];                 //The 3D image (same size as the original images) for save the overlapping information
    cv::Mat matColorSeq[MAXSequence];     //Here it is defined for save result as color image
   
    cv::Mat matConnectLoose[MAXSequence];
////////////////////////////Draw fibers with ID///////////////////////


    cv::Mat matFiberID[MAXSequence];
    unsigned char* ppuchFiberID[MAXSequence];
    CreateVoid3DColorImages(nLayers, nWidth,nHeight,matFiberID,ppuchFiberID);


///////////////////////////////////////////////////////////////////////



    unsigned char* ppuch3DLos[MAXSequence];
    CreateVoid3DGrayImages(nLayers, nWidth,nHeight, matConnectLoose, ppuch3DLos);
    
    unsigned char* ppuch3DRGB[MAXSequence];
    int* ppnOverlap[MAXSequence];
    int* ppnConn[MAXSequence];
 

    CreateVoid3DColorImages(nLayers, nWidth,nHeight, matColorSeq, ppuch3DRGB);
    CreateVoid3DIntImages(nLayers, nWidth,nHeight,matOverlapS,ppnOverlap);
    CreateVoid3DIntImages(nLayers, nWidth,nHeight,matConnectS,ppnConn);

    std::vector<Que3D*> vpCylinderObjects;     //The que for store the detected short fibers
    
    time_t mpp_begin;
    time(&mpp_begin); 

    for(int nT=0; nT<ITER_NUM; nT++)
    //for(int nT=0; nT<1; nT++)
    // Camilo iter
    {
        time_t iter_begin;
        time(&iter_begin);

        nGlobalCnt++;
        std::cout<<"The "<<nT<<"th iteration ... \n";
        
       int nBefBir = vpCylinderObjects.size();

       GrowStage3D(vpCylinderObjects, ppuch3DGray, ppuchSeg3D, ppnConn, ppnOverlap,nWidth,nHeight, nLayers, CylinderPool,pPredPool,conP);
       
       int nGrowBir = vpCylinderObjects.size();
        
       time_t iter_grow;      time(&iter_grow);
       std::cout<<"The growing time is "<<difftime(iter_grow,iter_begin)<<"s"<<'\n';    
   
       BirthStage3D(vpCylinderObjects, ppuch3DGray, ppuchSeg3D, ppnConn, ppnOverlap,nWidth,nHeight, nLayers, CylinderPool,pPredPool,conP);
       
    
       time_t iter_birth;      time(&iter_birth);
       std::cout<<"The birth time is "<<difftime(iter_birth,iter_grow)<<"s"<<'\n';

       int nAftBir = vpCylinderObjects.size();
       UpdatePriorConnect3D(vpCylinderObjects, ppnConn,nWidth,nHeight, nLayers, CylinderPool);
    
       DeathStage3D(vpCylinderObjects, ppuch3DGray, ppuchSeg3D, ppnConn, ppnOverlap,nWidth,nHeight, nLayers, CylinderPool,pPredPool,conP);
        
       time_t iter_death;      time(&iter_death);
       std::cout<<"The death time is "<<difftime(iter_death,iter_birth)<<"s"<<'\n';       
        
       int nAftDeath = vpCylinderObjects.size();
        
       UpdatePriorConnect3D(vpCylinderObjects, ppnConn,nWidth,nHeight, nLayers, CylinderPool);

       AdaptLocal3D(vpCylinderObjects, ppuch3DGray, ppuchSeg3D, ppnConn, ppnOverlap,nWidth,nHeight, nLayers, CylinderPool,pPredPool,conP);
      
       time_t iter_adapt;      time(&iter_adapt);
       std::cout<<"The adapt time is "<<difftime(iter_adapt,iter_death)<<"s"<<'\n';      


       DeathStage3D(vpCylinderObjects, ppuch3DGray, ppuchSeg3D, ppnConn, ppnOverlap,nWidth,nHeight, nLayers, CylinderPool,pPredPool,conP,false);
        
        
       UpdatePriorConnect3D(vpCylinderObjects, ppnConn,nWidth,nHeight, nLayers, CylinderPool);
    
       MergeSplit3D(vpCylinderObjects, ppuch3DGray, ppuchSeg3D, ppnConn,ppuch3DLos, ppnOverlap,nWidth,nHeight, nLayers, CylinderPool,pPredPool,conP);

       UpdatePriorConnect3D(vpCylinderObjects, ppnConn,nWidth,nHeight, nLayers, CylinderPool);

       time_t iter_merge;      time(&iter_merge);
       std::cout<<"The merge time is "<<difftime(iter_merge,iter_adapt)<<"s"<<'\n';      
    
       

        
        std::cout<<"There are "<<vpCylinderObjects.size()<<" objects in "<< nT<<"th iteration \n";
        std::cout<<"Grow Birth to "<< nGrowBir-nBefBir<< "objects \n";
        std::cout<<"Map Birth to "<< nAftBir-nGrowBir<< "objects,Then Kill "<<nAftBir-nAftDeath<<" objects \n";

        time_t iter_end;
        time(&iter_end);

        double seconds = difftime(iter_end,iter_begin);

        std::cout<<"The time spend on "<<nT<<"th iteration is "<<seconds<<"s"<<'\n';
       
        
    //Save Results in every 10 iteration
    #if 0   
       if(nT
          %10==0)
       {
           
          std::cout<<"The start image number is "<<nT/10*99+1<<'\n';
           for(int i=0; i<nLayers;i ++)
               matColorSeq[i].setTo(0);
           
           int nMinR=10;
           for(int i=0; i<int(vpCylinderObjects.size());i++)
           {
               
               char infor[64];
               sprintf(infor, "%d th Ojbect Info:/n", i);
            
               
               
               if(vpCylinderObjects[i]->m_mark.m_nR<nMinR)
               {
                   nMinR =vpCylinderObjects[i]->m_mark.m_nR;
               }
           }
           std::cout<<"The min R = "<<nMinR<<'\n';
           DrawResult3D(vpCylinderObjects, ppuch3DRGB, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
           SaveResult3D(matColorSeq,ppnConn, ppnOverlap,nWidth,nHeight, nLayers, nT);
           
       }
    #endif
       std::cout<<'\n';
    }
    
    time_t mpp_end;
    time(&mpp_end);
    double dTotalTime =  difftime(mpp_end,mpp_begin);
    std::cout<<"The total time for Cylinder MPP process is "<<dTotalTime<<"s \n";


    std::cout<<"Post Process: Merging......\n";
    
    std::vector<Que3D*> vp3DStartLong;
    
    GrowStage3DFinal(vpCylinderObjects, ppuch3DGray, ppuchSeg3D, ppnConn, ppnOverlap,nWidth,nHeight, nLayers, CylinderPool,pPredPool,conP);
    MergeShortToLongNew(vp3DStartLong, vpCylinderObjects,nWidth,  nHeight, nLayers,CylinderPool, pPredPool,ppuch3DLos);
    
    for(int i=0; i<nLayers;i ++)
    matColorSeq[i].setTo(0);
    
    
    std::cout<<"Post Process: Merging finished\n";
    time_t merge_end;
    time(&merge_end);
    double dMergeTime =  difftime(merge_end,mpp_end);
    std::cout<<"Merge process time is "<< dMergeTime<<"s \n";
    
    //DrawMergeShortID(vp3DStartLong, ppuch3DRGB, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    DrawMergeShortID_modified(vp3DStartLong, ppuch3DRGB, nWidth, nHeight, nLayers, CylinderPool, pPredPool, SubVn);
    Save3DImages("./SUBVOLUMES/sV" + std::to_string(SubVn) + "/fibers/mer_",nLayers, matColorSeq);
    WriteFiberResults(vp3DStartLong, SubVn);
 
    return;   //don't need converting for draw fibers.    
    
    std::cout<<"Converting... \n";
    std::vector<Que3D*> vp3DLongF;
    
    Clean3DIntImages(nLayers, nWidth,  nHeight,  ppnConn);
    Clean3DIntImages(nLayers, nWidth,  nHeight,  ppnOverlap);
    Clean3DGrayImages(nLayers, nWidth, nHeight, ppuchSeg3D);
    
    Convert2Long(vp3DStartLong, vp3DLongF,ppuch3DGray, ppuchSeg3D,ppnConn,ppnOverlap, nWidth, nHeight, nLayers,CylinderPool,pPredPool, conP);
    time_t cvt_end;
    time(&cvt_end);
    double dCvtTime =  difftime(cvt_end,merge_end);
    std::cout<<"Converting from short to long time is "<<dCvtTime<<"s \n"; 
    

    std::cout<<"The whole process spend time is "<<difftime(cvt_end,mpp_begin);



    for(int i=0; i<nLayers;i ++)
        matColorSeq[i].setTo(0);
    
    DrawResult3D(vp3DLongF, ppuch3DRGB, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    

    
    Save3DImages("./LongFiber/long_",nLayers, matColorSeq);
    
    for(int i=0; i<nLayers;i ++)
        matColorSeq[i].setTo(0);
    
    DrawResult3D(vp3DLongF, ppuch3DRGB, nWidth, nHeight, nLayers, CylinderPool, pPredPool,true);
    Save3DImages("./LongSeg/longSeg_",nLayers, matColorSeq);
   
    std::cout<<"finish converting...\n";
    
    
    std::cout<<"finish save fibers\n";
    
   // ShowSavedFibers(nWidth, nHeight, nLayers,"./Output/fiber.dat", CylinderPool, pPredPool);
    
}



void BirthStage3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,MppControPara& conP)
{
   
    
    for(int nz = 0; nz<nLayers; nz++)
        for(int ny = 0; ny<nHeight; ny++)
            for(int nx = 0; nx<nWidth; nx++)
            {
                
                if((Get3DPixel(nx, ny, nz, nWidth, ppuchSeg3D)==0)||//if 2D segmentation in (x,y,z) not to 0, then give birth by random
                    (Get3DPixel(nx, ny, nz, nWidth, ppnOverlap)!=0)||//no birth on the pixel already been oqqupied||
                   Get3DPixel(nx, ny, nz, nWidth, ppn3DConnect)!=0||
                   random2()>conP.dBirthTemp)//birth Control by temperature
                {
                    continue;
                }
                
                Que3D* qNewQue3D= new Que3D;

                (*qNewQue3D).m_dX = nx;
                (*qNewQue3D).m_dY = ny;
                (*qNewQue3D).m_dZ = nz;

                CMark markNew;
                markNew.m_nH = NMINH+ int(random2()*NDIFH+0.5); // fix H
                //random thetaY thetaZ
                markNew.m_nR = floor((random2()*(NDIFR)+NMINR)+0.5);

                while(markNew.m_nH<markNew.m_nR)
                {
                    markNew.m_nH = NMINH+ int(random2()*NDIFH+0.5); // fix H
                    //random thetaY thetaZ
                    markNew.m_nR = floor((random2()*(NDIFR)+NMINR)+0.5);
                }
                
                markNew.m_nthetaY = floor(random2()*NDEGREE_Y+0.5);
                
                markNew.m_nthetaZ = floor(random2()*NDEGREE_Z+0.5);
                //random nR
                qNewQue3D->m_mark = markNew;
                
                double dFrontConRatio=0, dBackConRatio=0;
                double dDataEn =1;
                double dOverlapE=0;;
                CalculateEnergy(qNewQue3D, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth,nHeight,nLayers,CylinderPool, pPredPool,dDataEn,dFrontConRatio,dBackConRatio,dOverlapE);
                
               
                if(dOverlapE>OVERLAP_T)
                    continue;
                
                double dPriorC = (dFrontConRatio+dBackConRatio)/2;
                if (dPriorC<CONNECTFINEP3D)
                    dPriorC = CONNECTPRIOR3D*(1-dPriorC/CONNECTFINEP3D);
                else
                    dPriorC =-CONNECTPRIOR3D*1.8*(sqrt(dPriorC)-sqrt(CONNECTFINEP3D))/(1-sqrt(CONNECTFINEP3D));
                
                //std::cout<<"Prior="<<dPriorC<<"Eng="<<dDataEn<<'\n';

                if((dDataEn+0.6*dPriorC<CYLINDER_ENG_BT||dDataEn<CYLINDER_ENG_BT)&& (dDataEn<BirthCylinder_ENG)) //give birth
                {
                    
                    qNewQue3D->m_dDateEng =LENPRIOR3D*exp(double(NMAXH-markNew.m_nH)/NMAXH)+dDataEn;
                    qNewQue3D->m_dPriorCon = dPriorC;
                    qNewQue3D->m_dFConnR = dFrontConRatio;
                    qNewQue3D->m_dBConnR = dBackConRatio;
                    
                    AddObjectOverLap_CONN(qNewQue3D,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers,CylinderPool, pPredPool);
                    
          //        AdaptOneLocal3D(qNewQue3D, ppuch3DGray, ppn3DConnect, ppnOverlap,nWidth, nHeight,nLayers, CylinderPool, pPredPool, 3);
                    
                    vObjects.push_back(qNewQue3D);
                    
                }
                else
                {
                    //printf("%f\n", dDataEn);
                    //printf("%d\n", dDataEn+0.6*dPriorC<CYLINDER_ENG_BT||dDataEn<CYLINDER_ENG_BT);
                    //printf("%d\n", dDataEn<BirthCylinder_ENG);
                    delete qNewQue3D;
                }
                
                
            }
    
    
    
    conP.dBirthTemp=conP.dBirthTemp*conP.dBirthERatio;
    
                
}

void DeathStage3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP,bool bIsOverOnly)
{
    std::sort(vObjects.begin(),vObjects.end(),cmp3D);
    
    for(int i=0; i<vObjects.size(); i++)
    {
        double dOverlap = OverlapRatio3D(vObjects[i], ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
        
        double dConnC =0;
        if(bIsOverOnly)
        {
            dConnC=0;
        }
        else if(vObjects[i]->m_dFConnR<0.2&&vObjects[i]->m_dBConnR<0.2)
        {
            dConnC= 0.25*conP.dBirthTemp;
            
            if(vObjects[i]->m_dDateEng<=1.8)
                dConnC=0.25*dConnC;
        }

        if(dOverlap>OVERLAP_T||random2()<dConnC)
        {
            
            EraseObjectOverLap_CONN(vObjects[i], ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
            Que3D* pqD = vObjects[i];
        
            
            vObjects.erase(vObjects.begin()+i);
            delete pqD;
            i--;
        }
        
    }
 //   conP.dBirthTemp*=0.98;
    
}

void AdaptOneLocal3D(Que3D* qOneObj,unsigned char** ppuch3DGray, int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, int nRange)
{
    
    Que3D qTempt = *qOneObj;
    Que3D qBest =  *qOneObj;
    
    EraseObjectOverLap_CONN(qOneObj, ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    

     for(int nTy=-12; nTy<=12; nTy++)
           for(int nTz=-12; nTz<=12; nTz++)
            {
                           
                 int nx=-2+random2()*4+0.5;
                 int ny=-2+random2()*4+0.5;  
                 int nz=-2+random2()*4+0.5;
                 int nH=-3+random2()*6+0.5;
                 int nR=-1+random2()*2+0.5;

                 qTempt.m_dX = nx+qOneObj->m_dX;
                 qTempt.m_dY = ny+qOneObj->m_dY;
                 qTempt.m_dZ = nz+qOneObj->m_dZ;
                                
                 qTempt.m_mark.m_nR = qOneObj->m_mark.m_nR+nR;
                                
                 if(qTempt.m_mark.m_nR<NMINR||NMINR>NMAXR)
                 {
                     continue;
                  }
                                
                  qTempt.m_mark.m_nH = qOneObj->m_mark.m_nH+nH;
                                
                 if(qTempt.m_mark.m_nH<NMINH||qTempt.m_mark.m_nH>NMAXH||qTempt.m_mark.m_nH<qTempt.m_mark.m_nR)
                    continue;
                                
                                
                 if(!IsInnerPoint(int((qTempt).m_dX+0.5), int((qTempt).m_dY+0.5), int((qTempt).m_dZ+0.5), nWidth, nHeight, nLayers))
                 {
                    continue;
                 }
                                
                 qTempt.m_mark.m_nthetaY = qOneObj->m_mark.m_nthetaY;
                 qTempt.m_mark.m_nthetaZ = qOneObj->m_mark.m_nthetaZ+nTz;
                                
                                
                                
                  if(qTempt.m_mark.m_nthetaZ<0)
                  {
                      qTempt.m_mark.m_nthetaZ+=DEGREE_Z;
                      qTempt.m_mark.m_nthetaY= DEGREE_Y - qTempt.m_mark.m_nthetaY;
                  }
                  else if(qTempt.m_mark.m_nthetaZ>=DEGREE_Z)
                  {
                       qTempt.m_mark.m_nthetaZ-=DEGREE_Z;
                       qTempt.m_mark.m_nthetaY= DEGREE_Y - qTempt.m_mark.m_nthetaY;
                  }
                                
                  qTempt.m_mark.m_nthetaY+=nTy;
                                
                  if(qTempt.m_mark.m_nthetaY<0)
                       qTempt.m_mark.m_nthetaY+=DEGREE_Y;
                   else if(qTempt.m_mark.m_nthetaY>=DEGREE_Y)
                       qTempt.m_mark.m_nthetaY-=DEGREE_Y;
                                
                                    
                   if(qTempt.m_mark.m_nR<NMINR||NMINR>NMAXR)
                    {
                         continue;
                     }
                                
                    double dDataE,dConnFRatio,dConnBRatio, dOverlapE;
                                
                    CalculateEnergy(&qTempt, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers, CylinderPool,pPredPool,dDataE, dConnFRatio,dConnBRatio, dOverlapE);
                                
                     qTempt.m_dOverlapR = dOverlapE;
                     if(qTempt.m_dOverlapR>OVERLAP_GT)
                     {
                         continue;
                     }
                                
                                
                     double dPriorC =dPriorCaculation(dConnFRatio, dConnBRatio);
                                
                                
                                
                     qTempt.m_dPriorCon = dPriorC;
                     qTempt.m_dFConnR = dConnFRatio;
                     qTempt.m_dBConnR = dConnBRatio;
                     qTempt.m_dDateEng = dDataE+LENPRIOR3D*exp(double(NMAXH-qTempt.m_mark.m_nH)/NMAXH);
                                
                         
                     if(qTempt.m_dDateEng<qBest.m_dDateEng&&qTempt.m_dDateEng<CYLINDER_ENG_T)
                     {
                        qBest = qTempt;
                            
                     }
                                
              }
    qTempt = qBest;
    
    
    
//Adapt Height

        for(int nR=-1; nR<=1;nR++)
        {
            
            qTempt.m_mark.m_nR = qOneObj->m_mark.m_nR+nR;
            
            if(qTempt.m_mark.m_nR<NMINR||NMINR>NMAXR)
            {
                continue;
            }
             
            int nTy = -1+random2()*2+0.5;
            int nTz = -1+random2()*2+0.5;
                 
            qTempt.m_mark.m_nthetaY = qOneObj->m_mark.m_nthetaY;
            qTempt.m_mark.m_nthetaZ = qOneObj->m_mark.m_nthetaZ+nTz;     
            
            if(qTempt.m_mark.m_nthetaZ<0)
            {
                qTempt.m_mark.m_nthetaZ+=DEGREE_Z;
                qTempt.m_mark.m_nthetaY= DEGREE_Y - qTempt.m_mark.m_nthetaY;
            }
            else if(qTempt.m_mark.m_nthetaZ>=DEGREE_Z)
            {
                qTempt.m_mark.m_nthetaZ-=DEGREE_Z;
                qTempt.m_mark.m_nthetaY= DEGREE_Y - qTempt.m_mark.m_nthetaY;
            }
            
            qTempt.m_mark.m_nthetaY+=nTy;
            
            if(qTempt.m_mark.m_nthetaY<0)
                qTempt.m_mark.m_nthetaY+=DEGREE_Y;
            else if(qTempt.m_mark.m_nthetaY>=DEGREE_Y)
                qTempt.m_mark.m_nthetaY-=DEGREE_Y;
            
            int nx=-1+random2()*2+0.5;
            int ny=-1+random2()*2+0.5;
            int nz=-1+random2()*2+0.5;
            int nH=-1+random2()*2+0.5;

            qTempt.m_dX = nx+qBest.m_dX;
            qTempt.m_dY = ny+qBest.m_dY;
            qTempt.m_dZ = nz+qBest.m_dZ;
                    
            qTempt.m_mark.m_nH = qBest.m_mark.m_nH+nH;

            if(qTempt.m_mark.m_nH<NMINH||qTempt.m_mark.m_nH>NMAXH||qTempt.m_mark.m_nH<qTempt.m_mark.m_nR)
                continue;
            
            
            if(!IsInnerPoint(int((qTempt).m_dX+0.5), int((qTempt).m_dY+0.5), int((qTempt).m_dZ+0.5), nWidth, nHeight, nLayers))
            {
                continue;
            }
          
          if(qTempt.m_mark.m_nH<NMINH||qTempt.m_mark.m_nH>NMAXH||qTempt.m_mark.m_nH<qTempt.m_mark.m_nR)
          continue;

            double dDataE,dConnFRatio,dConnBRatio, dOverlapE;
            
            CalculateEnergy(&qTempt, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers, CylinderPool,pPredPool,dDataE, dConnFRatio,dConnBRatio, dOverlapE);
            
            qTempt.m_dOverlapR = dOverlapE;
            if(qTempt.m_dOverlapR>OVERLAP_GT)
            {
                continue;
            }
            
            
            double dPriorC =dPriorCaculation(dConnFRatio, dConnBRatio);
            
            
            
            qTempt.m_dPriorCon = dPriorC;
            qTempt.m_dFConnR = dConnFRatio;
            qTempt.m_dBConnR = dConnBRatio;
            qTempt.m_dDateEng = dDataE+LENPRIOR3D*exp(double(NMAXH-qTempt.m_mark.m_nH)/NMAXH);
            
            if(qTempt.m_dDateEng<qBest.m_dDateEng&&qTempt.m_dDateEng<CYLINDER_ENG_T)
            {
                qBest = qTempt;
            }
            
        }
    
    qTempt = qBest;
    
    
    int nbH = qBest.m_mark.m_nH;
    int nbR = qBest.m_mark.m_nR;
//Adapt x,y,z
    for(int nx=-nRange; nx<=nRange; nx++)
     for(int ny=-nRange; ny<=nRange; ny++)
     for(int nz=-nRange; nz<=nRange; nz++)
     {
         
         int nTy = -3+random2()*6+0.5;
         int nTz = -3+random2()*6+0.5;
         
         
         qTempt.m_mark.m_nthetaY = qOneObj->m_mark.m_nthetaY;
         qTempt.m_mark.m_nthetaZ = qOneObj->m_mark.m_nthetaZ+nTz;
         
         
         
         if(qTempt.m_mark.m_nthetaZ<0)
         {
             qTempt.m_mark.m_nthetaZ+=DEGREE_Z;
             qTempt.m_mark.m_nthetaY= DEGREE_Y - qTempt.m_mark.m_nthetaY;
         }
         else if(qTempt.m_mark.m_nthetaZ>=DEGREE_Z)
         {
             qTempt.m_mark.m_nthetaZ-=DEGREE_Z;
             qTempt.m_mark.m_nthetaY= DEGREE_Y - qTempt.m_mark.m_nthetaY;
         }
         
         qTempt.m_mark.m_nthetaY+=nTy;
         
         if(qTempt.m_mark.m_nthetaY<0)
             qTempt.m_mark.m_nthetaY+=DEGREE_Y;
         else if(qTempt.m_mark.m_nthetaY>=DEGREE_Y)
             qTempt.m_mark.m_nthetaY-=DEGREE_Y;
         
         
         
          qTempt.m_dX = nx+qOneObj->m_dX;
          qTempt.m_dY = ny+qOneObj->m_dY;
          qTempt.m_dZ = nz+qOneObj->m_dZ;
         
         int nH=-nRange+random2()*2*nRange+0.5;
         int nR= -1+random2()*3+0.5;
         qTempt.m_mark.m_nH =nbH+nH;
         qTempt.m_mark.m_nR =MAX(NMINR,MIN(NMAXR,nbR+nR));

         if(qTempt.m_mark.m_nH<NMINH||qTempt.m_mark.m_nH>NMAXH||qTempt.m_mark.m_nH<qTempt.m_mark.m_nR)
             continue;
         
          if(!IsInnerPoint(int((qTempt).m_dX+0.5), int((qTempt).m_dY+0.5), int((qTempt).m_dZ+0.5), nWidth, nHeight, nLayers))
          {
          continue;
          }

         double dDataE,dConnFRatio,dConnBRatio, dOverlapE;
         
         CalculateEnergy(&qTempt, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers, CylinderPool,pPredPool,dDataE, dConnFRatio,dConnBRatio, dOverlapE);
         
         qTempt.m_dOverlapR = dOverlapE;
         if(qTempt.m_dOverlapR>OVERLAP_GT)
         {
             continue;
         }
         
         
         double dPriorC =dPriorCaculation(dConnFRatio, dConnBRatio);
         
            
            
            qTempt.m_dPriorCon = dPriorC;
            qTempt.m_dFConnR = dConnFRatio;
            qTempt.m_dBConnR = dConnBRatio;
            qTempt.m_dDateEng = dDataE+LENPRIOR3D*exp(double(NMAXH-qTempt.m_mark.m_nH)/NMAXH);
            
            if(qTempt.m_dDateEng+ qTempt.m_dPriorCon<qBest.m_dDateEng+ qBest.m_dPriorCon&&qTempt.m_dDateEng<CYLINDER_ENG_T)
            {
                qBest = qTempt;
            //    std::cout<<"adaptXYZ+++++++++++++++++++";
            }
            
        }
    qTempt = qBest;
  
    *qOneObj = qBest;
    
    AddObjectOverLap_CONN(qOneObj, ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    
}

void AdaptLocal3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP)
{
    
    for(int i=0; i<vObjects.size(); i++)
    {
        if(vObjects[i]->m_dFConnR>0.85&&vObjects[i]->m_dBConnR>0.85&&random2()>0.5)   //If both of the two ends of the fiber are connected well then, have 0.5 probability to adapt
        {
             continue;
        }

        AdaptOneLocal3D(vObjects[i], ppuch3DGray, ppn3DConnect, ppnOverlap,nWidth, nHeight,nLayers, CylinderPool, pPredPool, 3);
        
    }
}


void DrawResult3D(std::vector<Que3D*>& vFiber,unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, bool bSingleColor)
{
    std::cout<<"Total draw fiber is "<< vFiber.size()<<'\n';
    sRGB rgb[10]={ sRGB(0,255,0),sRGB(0,0,255),sRGB(0,255,255),sRGB(255,255,0),sRGB(255,0,0),sRGB(255,255,255),sRGB(255,0,255),sRGB(128,180,30),sRGB(20,60,180),sRGB(180,40,110)};
    
    
    
    for(int i=0; i<vFiber.size(); i++)
    {
        // Camilo Draw
        //sRGB stempt =rgb[i%10];
        sRGB stempt= sRGB(0,floor(i / 256), i%256);
        if(bSingleColor)
        {
            stempt = sRGB(255,255,255);
        }
        Que3D* pHeader = vFiber[i];
        DrawCylinderColor(round(pHeader->m_dX), round(pHeader->m_dY),round(pHeader->m_dZ), pHeader->m_mark,ppuch3Dadd,nWidth, nHeight, nLayers,CylinderPool,pPredPool,stempt);

    }
}

void SaveResult3D(cv::Mat* pmatColorSeq,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, int nFrame)
{
    

    char chSaveName[60];
    
    cv::Mat matTmpt(nHeight, nWidth, CV_8UC1);
    cv::Mat matTmptCon(nHeight, nWidth, CV_8UC1);
    unsigned char* puchMatTmp = (unsigned char*)matTmpt.data;
    unsigned char* puchMatCon = (unsigned char*)matTmptCon.data;

    for(int nz=0; nz<nLayers; nz++)
    {
       
        sprintf(chSaveName, ".././3DFiberR/detect%d_%d.jpg",nFrame,nz);
        cv::imwrite(chSaveName, pmatColorSeq[nz]);
        
        matTmpt.setTo(0);
        
        for(int ny=0; ny<nHeight; ny++)
            for(int nx=0; nx<nWidth; nx++)
            {
                if(Get3DPixel(nx, ny, nz, nWidth, ppnOverlap)>0)
                {
                    puchMatTmp[ny*nWidth+nx]=255;
                    
                }
                else
                {
                    puchMatTmp[ny*nWidth+nx]=0;
                }
         
                if(Get3DPixel(nx, ny, nz, nWidth, ppn3DConnect)>0)
                {
                    puchMatCon[ny*nWidth+nx]=255;
                    
                }
                else
                {
                    puchMatCon[ny*nWidth+nx]=0;
                }
                
            }
        
        cv::Mat matCTmpt, matCTmptCon;
        cv::cvtColor(matTmpt, matCTmpt, CV_GRAY2BGR);
        cv::cvtColor(matTmptCon, matCTmptCon, CV_GRAY2BGR);
        sprintf(chSaveName, "./3DSegFiber/seg%d_%d.jpg",nFrame,nz);
        cv::imwrite(chSaveName, matCTmpt);
        
        
        //sprintf(chSaveName, "./3DConFiber/con%d_%d.jpg",nFrame,nz);
        //cv::imwrite(chSaveName, matCTmptCon);
        
    }
}


void EvaluatingResults(double& dMissingRate, double& dFalseDRate, cv::Mat* pMat3DSeg, unsigned char** ppuch3DSeg, unsigned char** ppuch2DSeg,int nWidth, int nHeight, int nLayers)
{
    
    
    cv:: Mat element = cv::getStructuringElement(cv::MORPH_RECT, cv::Size(3,3),cv::Point(1,1));
    
    int nCnt2DTotalPts =0;
    int nCntMDetect =0;   //correct detect
    
    int nCntFDetect =0;   //false detect
    
    for(int nz = 0; nz<nLayers; nz++)
    {
        cv::Mat matTemptD, matTemptE;
        cv::dilate(pMat3DSeg[nz], matTemptD, element);
        cv::erode(pMat3DSeg[nz], matTemptE, element);
        unsigned char* puchDataD = matTemptD.data;//pMat3DSeg[nz].data;//matTemptD.data;
        unsigned char* puchDataE =matTemptE.data;//pMat3DSeg[nz].data;// matTemptE.data;
        
      //  unsigned char* puchDataD = pMat3DSeg[nz].data;
      //  unsigned char* puchDataE = pMat3DSeg[nz].data;
        for(int ny= 0; ny<nHeight; ny++)
            for(int nx=0; nx<nWidth; nx++)
            {
                
                if(Get3DPixel(nx, ny, nz, nWidth, ppuch2DSeg)!=0)
                {
                    nCnt2DTotalPts++;
                    if(puchDataD[ny*nWidth+nx]==0)
                    {
                        nCntMDetect++;
                    }
                }
                
                if(puchDataE[ny*nWidth+nx]!=0)
                {
                    if(Get3DPixel(nx, ny, nz, nWidth, ppuch2DSeg)==0)
                    {
                        nCntFDetect++;
                    }
                }
                
            }
        
    }
    
    dMissingRate = double(nCntMDetect)/nCnt2DTotalPts;
    dFalseDRate = double(nCntFDetect)/nCnt2DTotalPts;
    

}

void EvaluatingSegResultsFiles(double& dMissingRate, double& dFalseDRate, int nWidth, int nHeight, int nLayers, std::string str3DSeg, std::string str2DSeg)
{
    
    cv::Mat mat2DSeg[MAXSequence];
    cv::Mat mat3DSeg[MAXSequence];

    unsigned char* ppuch2DSeg[MAXSequence];
    unsigned char* ppuch3DSeg[MAXSequence];

    Load3DGrayImages(str2DSeg,nLayers, mat2DSeg, ppuch2DSeg);

    Load3DGrayImages(str3DSeg,nLayers, mat3DSeg, ppuch3DSeg, false);
    
    EvaluatingResults(dMissingRate, dFalseDRate, mat3DSeg, ppuch3DSeg, ppuch2DSeg, nWidth, nHeight, nLayers);
}

void GrowStage3DFinal(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP)
{
    Que3D* pTemptObj =NULL;
    int nNum = vObjects.size();
    for(int ni = 0; ni<nNum; ni++)
    {
        
        pTemptObj = vObjects[ni];
        
        //Grow from front end
        
        if(pTemptObj->m_dFConnR<GrowThresh)
        {
            //std::cout<<"Growing........ Front.................\n";
            Que3D* pNewQue = new Que3D;
            Que3D best3D;
            best3D.m_dDateEng=100;
            best3D.m_dPriorCon = 100;
            for(int nTY=-4; nTY<=4; nTY++)
                for(int nTZ=-4;nTZ<=4; nTZ++)
                    for(int nH =-8; nH<=8; nH++)
                    {
                        
                        int nNewTubeH  = pTemptObj->m_mark.m_nH+nH;
                        int nNewTubeR  = pTemptObj->m_mark.m_nR;//-1+random2()*2+0.5;
                        int nNewThetaY = pTemptObj->m_mark.m_nthetaY;
                        int nNewThetaZ = pTemptObj->m_mark.m_nthetaZ+nTZ;
                        
                        
                        if(nNewThetaZ<0)
                        {
                            nNewThetaZ+=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        else if(nNewThetaZ>=DEGREE_Z)
                        {
                            nNewThetaZ-=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        
                        nNewThetaY+=nTY;
                        
                        if(nNewThetaY<0)
                            nNewThetaY+=DEGREE_Y;
                        else if(nNewThetaY>=DEGREE_Y)
                            nNewThetaY-=DEGREE_Y;
                        
                        
                        if(nNewTubeH<nNewTubeR)
                            continue;
                        
                        if(nNewTubeH<NMINH||nNewTubeH>NMAXH||
                           nNewTubeR<NMINR||nNewTubeR>NMAXR)
                            continue;
                        
                        
                        
                        PreCMark temptPre = (*pPredPool)[pTemptObj->m_mark];
                        CMark newCmark = CMark(nNewTubeR, nNewTubeH,nNewThetaY,nNewThetaZ);
                        PreCMark premarkNew = (*pPredPool)[newCmark];
                        TDPoint pNewObjeCent;
                        
                        TDPoint pFrontP;
                        pFrontP.m_dX = temptPre.m_dFCenX;
                        pFrontP.m_dY = temptPre.m_dFCenY;
                        pFrontP.m_dZ = temptPre.m_dFCenZ;
                        
                        if(SelectDirect(pFrontP, premarkNew))
                        {
                            pNewObjeCent.m_dX = pFrontP.m_dX+pTemptObj->m_dX+ premarkNew.m_dFCenX;
                            pNewObjeCent.m_dY = pFrontP.m_dY+pTemptObj->m_dY+ premarkNew.m_dFCenY;
                            pNewObjeCent.m_dZ = pFrontP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dFCenZ;
                            
                        }
                        else
                        {
                            pNewObjeCent.m_dX = pFrontP.m_dX+pTemptObj->m_dX+ premarkNew.m_dBCenX;
                            pNewObjeCent.m_dY = pFrontP.m_dY+pTemptObj->m_dY+ premarkNew.m_dBCenY;
                            pNewObjeCent.m_dZ = pFrontP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dBCenZ;
                            
                        }
                        
                        
                        
                        
                        
                        if(!IsInnerPoint(pNewObjeCent.m_dX+0.5, pNewObjeCent.m_dY+0.5,pNewObjeCent.m_dZ+0.5, nWidth, nHeight, nLayers))
                        {
                            /* std::cout<<"-------------------------\n";
                             std::cout<<"tempt.x="<<pTemptObj->m_dX<<'\n';
                             std::cout<<"tempt.y="<<pTemptObj->m_dY<<'\n';
                             std::cout<<"tempt.z="<<pTemptObj->m_dZ<<'\n';
                             
                             std::cout<<"New.x="<<pNewObjeCent.m_dX<<'\n';
                             std::cout<<"New.y="<<pNewObjeCent.m_dY<<'\n';
                             std::cout<<"New.z="<<pNewObjeCent.m_dZ<<'\n';
                             std::cout<<"!!!!!!!!!!!!!!!!!!!\n";*/
                            continue;
                        }
                        
                        
                        
                        pNewQue->m_mark = newCmark;
                        pNewQue->m_Pmark = premarkNew;
                        pNewQue->m_dX = pNewObjeCent.m_dX;
                        pNewQue->m_dY = pNewObjeCent.m_dY;
                        pNewQue->m_dZ = pNewObjeCent.m_dZ;
                        
                        
                        double dDataE,dConnFRatio,dConnBRatio, dOverlapE;
                        
                        CalculateEnergy(pNewQue, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers, CylinderPool,pPredPool,dDataE, dConnFRatio,dConnBRatio, dOverlapE);
                        
                        pNewQue->m_dOverlapR = dOverlapE;
                        if(pNewQue->m_dOverlapR>OVERLAP_T1)
                        {
                            continue;
                        }
                        //      std::cout<<"-------------------------\n";
                        
                        double dPriorC =dPriorCaculation(dConnFRatio, dConnBRatio);
                        
                        
                        pNewQue->m_dPriorCon = dPriorC;
                        pNewQue->m_dFConnR = dConnFRatio;
                        pNewQue->m_dBConnR = dConnBRatio;
                        pNewQue->m_dDateEng = dDataE+LENPRIOR3D*exp(double(NMAXH-pNewQue->m_mark.m_nH)/NMAXH);
                        
                        //     std::cout<<"------bestEng="<<pNewQue->m_dDateEng<<'\n';
                        
                        
                        if(pNewQue->m_dPriorCon< best3D.m_dPriorCon&&(pNewQue->m_dDateEng<CYLINDER_ENG_FT))
                        {
                            best3D = *pNewQue;
                            
                        }
                        
                    }
            
            if(best3D.m_dDateEng<CYLINDER_ENG_FT&&pNewQue->m_dFConnR>0.5&&pNewQue->m_dBConnR>0.5)
            {
                //  std::cout<<"###########################\n";
                
                *pNewQue = best3D;
                
                AddObjectOverLap_CONN(pNewQue, ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
                
                AdaptOneLocal3D(pNewQue, ppuch3DGray, ppn3DConnect, ppnOverlap,nWidth, nHeight,nLayers, CylinderPool, pPredPool, 4);
                
                
                
                vObjects.push_back(pNewQue);
                
                //   std::cout<<"Grow FrontEnd..\n";
            }
            else
                delete pNewQue;
            
            
        }
        
        // Grow from back end
        
        if(pTemptObj->m_dBConnR<GrowThresh)
        {
            Que3D* pNewQue = new Que3D;
            Que3D best3D;
            best3D.m_dDateEng=100;
            best3D.m_dPriorCon=100;
            for(int nTY=-4; nTY<=4; nTY++)
                for(int nTZ=-4;nTZ<=4; nTZ++)
                    for(int nH =-8; nH<=8; nH++)
                    {
                        
                        int nNewTubeH  = pTemptObj->m_mark.m_nH+nH;
                        int nNewTubeR  = pTemptObj->m_mark.m_nR;//-1+random2()*2+0.5;
                        int nNewThetaY = pTemptObj->m_mark.m_nthetaY;
                        int nNewThetaZ = pTemptObj->m_mark.m_nthetaZ+nTZ;
                        
                        
                        if(nNewThetaZ<0)
                        {
                            nNewThetaZ+=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        else if(nNewThetaZ>=DEGREE_Z)
                        {
                            nNewThetaZ-=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        
                        nNewThetaY+=nTY;
                        
                        if(nNewThetaY<0)
                            nNewThetaY+=DEGREE_Y;
                        else if(nNewThetaY>=DEGREE_Y)
                            nNewThetaY-=DEGREE_Y;
                        
                        if(nNewTubeH<nNewTubeR)
                            continue;
                        
                        
                        if(nNewTubeH<NMINH||nNewTubeH>NMAXH||
                           nNewTubeR<NMINR||nNewTubeR>NMAXR)
                            continue;
                        
                        
                        PreCMark temptPre = (*pPredPool)[pTemptObj->m_mark];
                        CMark newCmark = CMark(nNewTubeR, nNewTubeH,nNewThetaY,nNewThetaZ);
                        PreCMark premarkNew = (*pPredPool)[newCmark];
                        TDPoint pNewObjeCent;
                        
                        TDPoint pBackP;
                        pBackP.m_dX = temptPre.m_dBCenX;
                        pBackP.m_dY = temptPre.m_dBCenY;
                        pBackP.m_dZ = temptPre.m_dBCenZ;
                        
                        if(SelectDirect(pBackP, premarkNew))
                        {
                            pNewObjeCent.m_dX = pBackP.m_dX+pTemptObj->m_dX+ premarkNew.m_dFCenX;
                            pNewObjeCent.m_dY = pBackP.m_dY+pTemptObj->m_dY+ premarkNew.m_dFCenY;
                            pNewObjeCent.m_dZ = pBackP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dFCenZ;
                            
                        }
                        else
                        {
                            pNewObjeCent.m_dX = pBackP.m_dX+pTemptObj->m_dX+ premarkNew.m_dBCenX;
                            pNewObjeCent.m_dY = pBackP.m_dY+pTemptObj->m_dY+ premarkNew.m_dBCenY;
                            pNewObjeCent.m_dZ = pBackP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dBCenZ;
                            
                        }
                        
                        /*
                         std::cout<<"-------------------------\n";
                         
                         
                         std::cout<<"backp.x="<<pBackP.m_dX<<'\n';
                         std::cout<<"backp.y="<<pBackP.m_dY<<'\n';
                         std::cout<<"backp.z="<<pBackP.m_dZ<<'\n';
                         
                         std::cout<<"premarkNew.m_dFCenX="<<premarkNew.m_dFCenX<<'\n';
                         std::cout<<"premarkNew.m_dFCenY="<<premarkNew.m_dFCenY<<'\n';
                         std::cout<<"premarkNew.m_dFCenZ="<<premarkNew.m_dFCenZ<<'\n';
                         
                         std::cout<<"premarkNew.m_dBCenX="<<premarkNew.m_dBCenX<<'\n';
                         std::cout<<"premarkNew.m_dBCenY="<<premarkNew.m_dBCenY<<'\n';
                         std::cout<<"premarkNew.m_dBCenZ="<<premarkNew.m_dBCenZ<<'\n';
                         
                         
                         std::cout<<"tempt.x="<<pTemptObj->m_dX<<'\n';
                         std::cout<<"tempt.y="<<pTemptObj->m_dY<<'\n';
                         std::cout<<"tempt.z="<<pTemptObj->m_dZ<<'\n';
                         
                         std::cout<<"New.x="<<pNewObjeCent.m_dX<<'\n';
                         std::cout<<"New.y="<<pNewObjeCent.m_dY<<'\n';
                         std::cout<<"New.z="<<pNewObjeCent.m_dZ<<'\n';
                         std::cout<<"!!!!!!!!!!!!!!!!!!!\n";*/
                        //   cv::namedWindow("sfsdfs");
                        //   cv::waitKey(0);
                        
                        if(!IsInnerPoint(pNewObjeCent.m_dX+0.5, pNewObjeCent.m_dY+0.5,pNewObjeCent.m_dZ+0.5, nWidth, nHeight, nLayers))
                        {
                            continue;
                        }
                        
                        pNewQue->m_mark = newCmark;
                        pNewQue->m_Pmark = premarkNew;
                        pNewQue->m_dX = pNewObjeCent.m_dX;
                        pNewQue->m_dY = pNewObjeCent.m_dY;
                        pNewQue->m_dZ = pNewObjeCent.m_dZ;
                        
                        
                        double dDataE,dConnFRatio,dConnBRatio, dOverlapE;
                        
                        CalculateEnergy(pNewQue, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers, CylinderPool,pPredPool,dDataE, dConnFRatio,dConnBRatio, dOverlapE);
                        
                        pNewQue->m_dOverlapR = dOverlapE;
                        if(pNewQue->m_dOverlapR>OVERLAP_T1)
                        {
                            continue;
                        }
                        
                        
                        double dPriorC =dPriorCaculation(dConnFRatio, dConnBRatio);
                        
                        
                        pNewQue->m_dPriorCon = dPriorC;
                        pNewQue->m_dFConnR = dConnFRatio;
                        pNewQue->m_dBConnR = dConnBRatio;
                        pNewQue->m_dDateEng = dDataE+LENPRIOR3D*exp(double(NMAXH-pNewQue->m_mark.m_nH)/NMAXH);
                        
                        if(pNewQue->m_dPriorCon< best3D.m_dPriorCon&&(pNewQue->m_dDateEng<CYLINDER_ENG_FT))
                        {
                            best3D = *pNewQue;
                        }
                        
                    }
            
            if(best3D.m_dDateEng<CYLINDER_ENG_FT&&pNewQue->m_dFConnR>0.5&&pNewQue->m_dBConnR>0.5)
            {
                
                
                *pNewQue = best3D;
                
                AddObjectOverLap_CONN(pNewQue, ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
                
                AdaptOneLocal3D(pNewQue, ppuch3DGray, ppn3DConnect, ppnOverlap,nWidth, nHeight,nLayers, CylinderPool, pPredPool, 4);
                
                
                
                vObjects.push_back(pNewQue);
            }
            else
                delete pNewQue;
            
        }
        
    }
    
    
    
}



void GrowStage3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP)
{
   
    Que3D* pTemptObj =NULL;
    int nNum = vObjects.size();
    for(int ni = 0; ni<nNum; ni++)
    {
        
        pTemptObj = vObjects[ni];

        //Grow from front end
        
        if(pTemptObj->m_dFConnR<GrowThresh&&random2()>pTemptObj->m_dFConnR)
        {
           //std::cout<<"Growing........ Front.................\n";
            Que3D* pNewQue = new Que3D;
            Que3D best3D;
            best3D.m_dDateEng=100;
            best3D.m_dPriorCon = 100;
            for(int nTY=-3; nTY<=3; nTY++)
                for(int nTZ=-3;nTZ<=3; nTZ++)
                    for(int nH =-4; nH<=3; nH++)
                    {
                        
                        int nNewTubeH  = pTemptObj->m_mark.m_nH+nH;
                        int nNewTubeR  = pTemptObj->m_mark.m_nR;//-1+random2()*2+0.5;
                        int nNewThetaY = pTemptObj->m_mark.m_nthetaY;
                        int nNewThetaZ = pTemptObj->m_mark.m_nthetaZ+nTZ;
                        
                        
                        if(nNewThetaZ<0)
                        {
                            nNewThetaZ+=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        else if(nNewThetaZ>=DEGREE_Z)
                        {
                            nNewThetaZ-=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        
                        nNewThetaY+=nTY;
                        
                        if(nNewThetaY<0)
                            nNewThetaY+=DEGREE_Y;
                        else if(nNewThetaY>=DEGREE_Y)
                            nNewThetaY-=DEGREE_Y;
                        
                        
                        if(nNewTubeH<nNewTubeR)
                            continue;
                        
                        if(nNewTubeH<NMINH||nNewTubeH>NMAXH||
                           nNewTubeR<NMINR||nNewTubeR>NMAXR)
                            continue;
                        

                        
                        PreCMark temptPre = (*pPredPool)[pTemptObj->m_mark];
                        CMark newCmark = CMark(nNewTubeR, nNewTubeH,nNewThetaY,nNewThetaZ);
                        PreCMark premarkNew = (*pPredPool)[newCmark];
                        TDPoint pNewObjeCent;
                        
                        TDPoint pFrontP;
                        pFrontP.m_dX = temptPre.m_dFCenX;
                        pFrontP.m_dY = temptPre.m_dFCenY;
                        pFrontP.m_dZ = temptPre.m_dFCenZ;
                       
                        if(SelectDirect(pFrontP, premarkNew))
                        {
                            pNewObjeCent.m_dX = pFrontP.m_dX+pTemptObj->m_dX+ premarkNew.m_dFCenX;
                            pNewObjeCent.m_dY = pFrontP.m_dY+pTemptObj->m_dY+ premarkNew.m_dFCenY;
                            pNewObjeCent.m_dZ = pFrontP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dFCenZ;
                   
                        }
                        else
                        {
                            pNewObjeCent.m_dX = pFrontP.m_dX+pTemptObj->m_dX+ premarkNew.m_dBCenX;
                            pNewObjeCent.m_dY = pFrontP.m_dY+pTemptObj->m_dY+ premarkNew.m_dBCenY;
                            pNewObjeCent.m_dZ = pFrontP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dBCenZ;
                     
                        }
                        
      
                        
                        
                        
                        if(!IsInnerPoint(pNewObjeCent.m_dX+0.5, pNewObjeCent.m_dY+0.5,pNewObjeCent.m_dZ+0.5, nWidth, nHeight, nLayers))
                        {
                           /* std::cout<<"-------------------------\n";
                            std::cout<<"tempt.x="<<pTemptObj->m_dX<<'\n';
                            std::cout<<"tempt.y="<<pTemptObj->m_dY<<'\n';
                            std::cout<<"tempt.z="<<pTemptObj->m_dZ<<'\n';
                            
                            std::cout<<"New.x="<<pNewObjeCent.m_dX<<'\n';
                            std::cout<<"New.y="<<pNewObjeCent.m_dY<<'\n';
                            std::cout<<"New.z="<<pNewObjeCent.m_dZ<<'\n';
                            std::cout<<"!!!!!!!!!!!!!!!!!!!\n";*/
                            continue;
                        }
                        
   
                        
                        pNewQue->m_mark = newCmark;
                        pNewQue->m_Pmark = premarkNew;
                        pNewQue->m_dX = pNewObjeCent.m_dX;
                        pNewQue->m_dY = pNewObjeCent.m_dY;
                        pNewQue->m_dZ = pNewObjeCent.m_dZ;
                        
                        
                        double dDataE,dConnFRatio,dConnBRatio, dOverlapE;
                        
                        CalculateEnergy(pNewQue, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers, CylinderPool,pPredPool,dDataE, dConnFRatio,dConnBRatio, dOverlapE);
                        
                        pNewQue->m_dOverlapR = dOverlapE;
                        if(pNewQue->m_dOverlapR>OVERLAP_T1)
                        {
                            continue;
                        }
                  //      std::cout<<"-------------------------\n";
                    
                        double dPriorC =dPriorCaculation(dConnFRatio, dConnBRatio);
                        
                        
                        pNewQue->m_dPriorCon = dPriorC;
                        pNewQue->m_dFConnR = dConnFRatio;
                        pNewQue->m_dBConnR = dConnBRatio;
                        pNewQue->m_dDateEng = dDataE+LENPRIOR3D*exp(double(NMAXH-pNewQue->m_mark.m_nH)/NMAXH);
                        
                   //     std::cout<<"------bestEng="<<pNewQue->m_dDateEng<<'\n';

                        
                        if(pNewQue->m_dDateEng+ pNewQue->m_dPriorCon< best3D.m_dDateEng+ best3D.m_dPriorCon&&pNewQue->m_dDateEng<CYLINDER_ENG_T)
                        {
                            best3D = *pNewQue;
                           
                        }
                        
                    }
            
            if(best3D.m_dDateEng<CYLINDER_ENG_T)
            {
               //  std::cout<<"###########################\n";
                
                *pNewQue = best3D;
                
                AddObjectOverLap_CONN(pNewQue, ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
                
//              AdaptOneLocal3D(pNewQue, ppuch3DGray, ppn3DConnect, ppnOverlap,nWidth, nHeight,nLayers, CylinderPool, pPredPool, 4);
                
                
                vObjects.push_back(pNewQue);

             //   std::cout<<"Grow FrontEnd..\n";
            }
            else
                delete pNewQue;
        
            
        }
        
// Grow from back end
        
        if(pTemptObj->m_dBConnR<GrowThresh&&random2()>pTemptObj->m_dBConnR)
        {
            Que3D* pNewQue = new Que3D;
            Que3D best3D;
            best3D.m_dDateEng=100;
            best3D.m_dPriorCon=100;
            for(int nTY=-3; nTY<=3; nTY++)
                for(int nTZ=-3;nTZ<=3; nTZ++)
                    for(int nH =-4; nH<=3; nH++)
                    {
                        
                        int nNewTubeH  = pTemptObj->m_mark.m_nH+nH;
                        int nNewTubeR  = pTemptObj->m_mark.m_nR;//-1+random2()*2+0.5;
                        int nNewThetaY = pTemptObj->m_mark.m_nthetaY;
                        int nNewThetaZ = pTemptObj->m_mark.m_nthetaZ+nTZ;
                        
                        
                        if(nNewThetaZ<0)
                        {
                            nNewThetaZ+=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        else if(nNewThetaZ>=DEGREE_Z)
                        {
                            nNewThetaZ-=DEGREE_Z;
                            nNewThetaY= DEGREE_Y - nNewThetaY;
                        }
                        
                        nNewThetaY+=nTY;
                        
                        if(nNewThetaY<0)
                            nNewThetaY+=DEGREE_Y;
                        else if(nNewThetaY>=DEGREE_Y)
                            nNewThetaY-=DEGREE_Y;
                        
                        if(nNewTubeH<nNewTubeR)
                            continue;
                        
                        
                        if(nNewTubeH<NMINH||nNewTubeH>NMAXH||
                           nNewTubeR<NMINR||nNewTubeR>NMAXR)
                            continue;
  
                        
                        PreCMark temptPre = (*pPredPool)[pTemptObj->m_mark];
                        CMark newCmark = CMark(nNewTubeR, nNewTubeH,nNewThetaY,nNewThetaZ);
                        PreCMark premarkNew = (*pPredPool)[newCmark];
                        TDPoint pNewObjeCent;
                        
                        TDPoint pBackP;
                        pBackP.m_dX = temptPre.m_dBCenX;
                        pBackP.m_dY = temptPre.m_dBCenY;
                        pBackP.m_dZ = temptPre.m_dBCenZ;
                        
                        if(SelectDirect(pBackP, premarkNew))
                        {
                            pNewObjeCent.m_dX = pBackP.m_dX+pTemptObj->m_dX+ premarkNew.m_dFCenX;
                            pNewObjeCent.m_dY = pBackP.m_dY+pTemptObj->m_dY+ premarkNew.m_dFCenY;
                            pNewObjeCent.m_dZ = pBackP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dFCenZ;
                            
                        }
                        else
                        {
                            pNewObjeCent.m_dX = pBackP.m_dX+pTemptObj->m_dX+ premarkNew.m_dBCenX;
                            pNewObjeCent.m_dY = pBackP.m_dY+pTemptObj->m_dY+ premarkNew.m_dBCenY;
                            pNewObjeCent.m_dZ = pBackP.m_dZ+pTemptObj->m_dZ+ premarkNew.m_dBCenZ;
                            
                        }
                        
                      /*
                        std::cout<<"-------------------------\n";
                        
                        
                        std::cout<<"backp.x="<<pBackP.m_dX<<'\n';
                        std::cout<<"backp.y="<<pBackP.m_dY<<'\n';
                        std::cout<<"backp.z="<<pBackP.m_dZ<<'\n';
                        
                        std::cout<<"premarkNew.m_dFCenX="<<premarkNew.m_dFCenX<<'\n';
                        std::cout<<"premarkNew.m_dFCenY="<<premarkNew.m_dFCenY<<'\n';
                        std::cout<<"premarkNew.m_dFCenZ="<<premarkNew.m_dFCenZ<<'\n';
                        
                        std::cout<<"premarkNew.m_dBCenX="<<premarkNew.m_dBCenX<<'\n';
                        std::cout<<"premarkNew.m_dBCenY="<<premarkNew.m_dBCenY<<'\n';
                        std::cout<<"premarkNew.m_dBCenZ="<<premarkNew.m_dBCenZ<<'\n';
                        
                        
                        std::cout<<"tempt.x="<<pTemptObj->m_dX<<'\n';
                        std::cout<<"tempt.y="<<pTemptObj->m_dY<<'\n';
                        std::cout<<"tempt.z="<<pTemptObj->m_dZ<<'\n';
                        
                        std::cout<<"New.x="<<pNewObjeCent.m_dX<<'\n';
                        std::cout<<"New.y="<<pNewObjeCent.m_dY<<'\n';
                        std::cout<<"New.z="<<pNewObjeCent.m_dZ<<'\n';
                        std::cout<<"!!!!!!!!!!!!!!!!!!!\n";*/
                     //   cv::namedWindow("sfsdfs");
                     //   cv::waitKey(0);
                       
                        if(!IsInnerPoint(pNewObjeCent.m_dX+0.5, pNewObjeCent.m_dY+0.5,pNewObjeCent.m_dZ+0.5, nWidth, nHeight, nLayers))
                        {
                            continue;
                        }
                        
                        pNewQue->m_mark = newCmark;
                        pNewQue->m_Pmark = premarkNew;
                        pNewQue->m_dX = pNewObjeCent.m_dX;
                        pNewQue->m_dY = pNewObjeCent.m_dY;
                        pNewQue->m_dZ = pNewObjeCent.m_dZ;
                        
                        
                        double dDataE,dConnFRatio,dConnBRatio, dOverlapE;
                        
                        CalculateEnergy(pNewQue, ppuch3DGray,ppn3DConnect,ppnOverlap,nWidth, nHeight,nLayers, CylinderPool,pPredPool,dDataE, dConnFRatio,dConnBRatio, dOverlapE);
                        
                        pNewQue->m_dOverlapR = dOverlapE;
                        if(pNewQue->m_dOverlapR>OVERLAP_T1)
                        {
                            continue;
                        }
                        
                        
                        double dPriorC =dPriorCaculation(dConnFRatio, dConnBRatio);
                        
                        
                        pNewQue->m_dPriorCon = dPriorC;
                        pNewQue->m_dFConnR = dConnFRatio;
                        pNewQue->m_dBConnR = dConnBRatio;
                        pNewQue->m_dDateEng = dDataE+LENPRIOR3D*exp(double(NMAXH-pNewQue->m_mark.m_nH)/NMAXH);
                        
                        if(pNewQue->m_dDateEng+ pNewQue->m_dPriorCon< best3D.m_dDateEng+ best3D.m_dPriorCon&&pNewQue->m_dDateEng<CYLINDER_ENG_T)
                        {
                            best3D = *pNewQue;
                        }
                        
                    }
            
            if(best3D.m_dDateEng<CYLINDER_ENG_T)
            {
                
               
                *pNewQue = best3D;
                
                AddObjectOverLap_CONN(pNewQue, ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
                
             //   AdaptOneLocal3D(pNewQue, ppuch3DGray, ppn3DConnect, ppnOverlap,nWidth, nHeight,nLayers, CylinderPool, pPredPool, 4);
                           
                vObjects.push_back(pNewQue);
            }
            else
                delete pNewQue;
            
        }
        
    }
    
    
    
}


double dPriorCaculation(double dConnFRatio, double dConnBRatio)
{
    double dPriorC = (dConnFRatio+dConnBRatio)/2;
    
    if (dPriorC<CONNECTFINEP3D)
        //    dPriorC = -finePos/CONNECTFINEP*dPriorC+finePos;
        dPriorC = CONNECTPRIOR3D*(1-dPriorC/CONNECTFINEP3D);
    else
        //    dPriorC =1.5*(-dPriorC*dPriorC+CONNECTFINEP*CONNECTFINEP)/(1-CONNECTFINEP*CONNECTFINEP);
        dPriorC =-CONNECTPRIOR3D*1.8*(sqrt(dPriorC)-sqrt(CONNECTFINEP3D))/(1-sqrt(CONNECTFINEP3D));
    
    if(dConnBRatio>DoubleConnThresh&&dConnFRatio>DoubleConnThresh)
    {
        dPriorC+=DoubleConnBonus;
    }
    
        
    return dPriorC;
}



void GetCoarseSeg3DImages(int nLayers, cv::Mat* matOrigin,cv::Mat* matSeg, double dThresh)
{

    for(int i=0; i<nLayers; i++)
    {
        
        cv::threshold(matOrigin[i], matSeg[i], dThresh, 255, CV_THRESH_BINARY);
        
    }
    
}
void Clean3DGrayImages(int nLayers, int nWidth, int nHeight, unsigned char** ppuch3Dadd)
{
    int nLenImage = nWidth*nHeight;
    
    for(int z=0; z<nLayers;z++)
    {
        memset(ppuch3Dadd[z], 0, nLenImage);
    }

}

void Clean3DIntImages(int nLayers, int nWidth, int nHeight, int** ppn3Dadd)
{
    int nLenImage = nWidth*nHeight*sizeof(int);
    
    for(int z=0; z<nLayers;z++)
    {
        memset(ppn3Dadd[z], 0, nLenImage);
    }
}
bool ReplaceObjects(int i,TDPoint ptMid, TDPoint ptI,std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,unsigned char** ppuch3DConLocal,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP)
{
    
    std::vector<Que3D> vSplitObj;
    vSplitObj.clear();
    
    SplitByPoints(vSplitObj,vObjects[i]->m_mark.m_nR, ptMid,ptI,NULL,nWidth, nHeight,  nLayers, CylinderPool, pPredPool);
    
    
   
    
    bool bShallR =true;
    
    for(int nS=0; nS<vSplitObj.size();nS++)
    {
        double dE = CalculateEnergy(&vSplitObj[nS], ppuch3DGray, nWidth, nHeight,nLayers,  CylinderPool,pPredPool);
        
        vSplitObj[nS].m_dDateEng =dE+LENPRIOR3D*exp(double(NMAXH-vSplitObj[nS].m_mark.m_nH)/NMAXH);
        
        
        if(vSplitObj[nS].m_mark.m_nH<vSplitObj[nS].m_mark.m_nR)
        {
         //   std::cout<<"H<R Error! Not replaced!\n";
            bShallR=false;
            break;
        }
        
        if(dE>CYLINDER_ENG_T)
        {
            
          //  std::cout<<"Energy error! Not replaced!\n";
            bShallR=false;
            break;
        }
    }
    
    if(bShallR&&vSplitObj.size()>0)
    {
        
      //  std::cout<<"Replace!\n";
        EraseObjectOverLap_CONN(vObjects[i], ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
        
     //   std::vector<Que3D> vSplitObj1;
        
    //    vSplitObj1.push_back(*vObjects[i]);
    
        *vObjects[i]=vSplitObj[0];
        
     //   OutputObjInfo(vObjects[i], "BBBBBBBBBBBBBBBBBB/n");
        
        
        AddObjectOverLap_CONN(vObjects[i], ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
        
        
    //    SaveDrawObj(vSplitObj, nWidth, nHeight, nLayers, CylinderPool, pPredPool, "./looseCon/los_");
        
        
        
    //    SaveDrawObj(vSplitObj1, nWidth, nHeight, nLayers, CylinderPool, pPredPool, "./looseCon1/los_");
        
    //    cv::Mat mat=cv::Mat::zeros(100, 100, CV_8UC1);
    //    cv::imshow("df", mat);
        
    //    cv::waitKey(0);
        
        if(vSplitObj.size()>1)
        {
            Que3D* pnewQue=new Que3D;
            *pnewQue = vSplitObj[1];
            
            
            
  //          OutputObjInfo(pnewQue, "AAAAAAAAAAAAAAAAAAAAA/n");


            AddObjectOverLap_CONN(pnewQue, ppn3DConnect, ppnOverlap, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
            
            vObjects.push_back(pnewQue);
        }
    }
    else
    {
      //   SaveDrawObj(vSplitObj, nWidth, nHeight, nLayers, CylinderPool, pPredPool, "./looseCon/los_");
    }
    
    
    return bShallR;
    
}
void MergeSplit3D(std::vector<Que3D*>& vObjects, unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,unsigned char** ppuch3DConLocal, int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP)
{
    
    int nRTimes=5;
    Clean3DGrayImages(nLayers, nWidth, nHeight, ppuch3DConLocal);

    
    
    CalculateNeigborID(vObjects, pow(NMAXH*6,2));
    
    int nNumObj = vObjects.size();
    
    for(int i=0; i<nNumObj; i++)
    {
        
        if(vObjects[i]->m_bFrontM==true)
        {
            
       //   std::cout<<"******************************************";
          AddObjectLocalCONN(vObjects[i],ppuch3DConLocal,nWidth, nHeight,nLayers,CylinderPool,pPredPool,3);
          for(int j=i+1; j<nNumObj; j++)
          {
            
              
             if(std::find(vObjects[i]->m_vnNeighborIDs.begin(), vObjects[i]->m_vnNeighborIDs.end(), j) == vObjects[i]->m_vnNeighborIDs.end())
             {
                 continue;
             }
              
              
            
              
             if(vObjects[j]->m_bBackM==false&&vObjects[j]->m_bFrontM==false)
             {
                 continue;
             }
             int nConInf = CheckConnObjIJ(vObjects[i], vObjects[j], nWidth, nHeight, nLayers,CylinderPool,pPredPool,true);
            
            
              
              
             if(nConInf==0)
                 continue;
            
        //      std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&";
              
          //  vObjects[i]->m_bFrontM=false;
              
#if TESTMERGE==1
              /**************/
            int nTemptNum = vObjects.size();
            Que3D objectI,objectJ;
            objectI = *vObjects[i];
            objectJ = *vObjects[j];
      
              /************/
#endif
            bool bIsRF=true,bIsRB=true;
              
            if(nConInf==1)
            {
                if(vObjects[j]->m_bFrontM==false)
                    break;
                
           //     std::cout<<"I F rate is "<<vObjects[i]->m_dFConnR<<"\n";
           //     std::cout<<"J F rate is "<<vObjects[j]->m_dFConnR<<"\n";
                
                
                TDPoint ptI,ptJ,ptMid;
                
                ptMid.m_dX= (vObjects[i]->m_dX+(*pPredPool)[vObjects[i]->m_mark].m_dFCenX+vObjects[j]->m_dX+(*pPredPool)[vObjects[j]->m_mark].m_dFCenX)/2.0;
                
                ptMid.m_dY= (vObjects[i]->m_dY+(*pPredPool)[vObjects[i]->m_mark].m_dFCenY+vObjects[j]->m_dY+(*pPredPool)[vObjects[j]->m_mark].m_dFCenY)/2.0;
                
                ptMid.m_dZ= (vObjects[i]->m_dZ+(*pPredPool)[vObjects[i]->m_mark].m_dFCenZ+vObjects[j]->m_dZ+(*pPredPool)[vObjects[j]->m_mark].m_dFCenZ)/2.0;
                
                ptI.m_dX = (*pPredPool)[vObjects[i]->m_mark].m_dBCenX+vObjects[i]->m_dX;
                ptI.m_dY = (*pPredPool)[vObjects[i]->m_mark].m_dBCenY+vObjects[i]->m_dY;
                ptI.m_dZ = (*pPredPool)[vObjects[i]->m_mark].m_dBCenZ+vObjects[i]->m_dZ;
                
                
                ptJ.m_dX = (*pPredPool)[vObjects[j]->m_mark].m_dBCenX+vObjects[j]->m_dX;
                ptJ.m_dY = (*pPredPool)[vObjects[j]->m_mark].m_dBCenY+vObjects[j]->m_dY;
                ptJ.m_dZ = (*pPredPool)[vObjects[j]->m_mark].m_dBCenZ+vObjects[j]->m_dZ;
               
                TDPoint ptIC=ptI;
                TDPoint ptJC=ptJ;
                TDPoint ptMidC=ptMid;
                int nRT=nRTimes;
                
             //   OutputObjInfo(vObjects[i], "IIIIIIIIIIIIIIIIIII/n");
             //   OutputObjInfo(vObjects[j], "JJJJJJJJJJJJJJJJJJJ/n");
                
                while(nRT-->0)
                {
                    
                    ptI = ptIC; ptJ=ptJC; ptMid=ptMidC;
                    RandomX(ptI,1);
                    RandomX(ptJ,1);
                    RandomX(ptMid,1);
                
                    
                    
                    
                    bIsRF=ReplaceObjects(i, ptMid, ptI, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                
                    bIsRB=ReplaceObjects(j, ptMid, ptJ, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                    if(bIsRF==true||bIsRB==true)
                        break;
                }
                
                vObjects[j]->m_bFrontM=false;
                vObjects[i]->m_bFrontM=false;
                
                
            }
          
              
            if(nConInf==2)
            {
                if(vObjects[j]->m_bBackM==false)
                    break;
                
          //      std::cout<<"I F rate is "<<vObjects[i]->m_dFConnR<<"\n";
          //      std::cout<<"J B rate is "<<vObjects[j]->m_dBConnR<<"\n";
                TDPoint ptI,ptJ,ptMid;
                
                ptMid.m_dX= (vObjects[i]->m_dX+(*pPredPool)[vObjects[i]->m_mark].m_dFCenX+vObjects[j]->m_dX+(*pPredPool)[vObjects[j]->m_mark].m_dBCenX)/2.0;
                
                ptMid.m_dY= (vObjects[i]->m_dY+(*pPredPool)[vObjects[i]->m_mark].m_dFCenY+vObjects[j]->m_dY+(*pPredPool)[vObjects[j]->m_mark].m_dBCenY)/2.0;
                
                ptMid.m_dZ= (vObjects[i]->m_dZ+(*pPredPool)[vObjects[i]->m_mark].m_dFCenZ+vObjects[j]->m_dZ+(*pPredPool)[vObjects[j]->m_mark].m_dBCenZ)/2.0;
                
                ptI.m_dX = (*pPredPool)[vObjects[i]->m_mark].m_dBCenX+vObjects[i]->m_dX;
                ptI.m_dY = (*pPredPool)[vObjects[i]->m_mark].m_dBCenY+vObjects[i]->m_dY;
                ptI.m_dZ = (*pPredPool)[vObjects[i]->m_mark].m_dBCenZ+vObjects[i]->m_dZ;
                
                
                ptJ.m_dX = (*pPredPool)[vObjects[j]->m_mark].m_dFCenX+vObjects[j]->m_dX;
                ptJ.m_dY = (*pPredPool)[vObjects[j]->m_mark].m_dFCenY+vObjects[j]->m_dY;
                ptJ.m_dZ = (*pPredPool)[vObjects[j]->m_mark].m_dFCenZ+vObjects[j]->m_dZ;
                
                TDPoint ptIC=ptI;
                TDPoint ptJC=ptJ;
                TDPoint ptMidC=ptMid;
                int nRT=nRTimes;
                while(nRT-->0)
                {
                    
                    ptI = ptIC; ptJ=ptJC; ptMid=ptMidC;
                    RandomX(ptI,1);
                    RandomX(ptJ,1);
                    RandomX(ptMid,2);
                    
                    bIsRF=ReplaceObjects(i, ptMid, ptI, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                    
                    bIsRB=ReplaceObjects(j, ptMid, ptJ, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                    if(bIsRF==true||bIsRB==true)
                        break;
                }
                
                vObjects[i]->m_bFrontM=false;
                vObjects[j]->m_bBackM=false;
                
            }
           
 
#if TESTMERGE==1
           if((bIsRF==false||bIsRB==false)&&nGlobalCnt>100)
           {
               
            
            
            OutputObjInfo(vObjects[i], "aftermergeI");
            OutputObjInfo(&objectI, "beforemergeI");
            std::vector<Que3D> vSplitObjects;
            vSplitObjects.push_back(*vObjects[i]);
            vSplitObjects.push_back(*vObjects[j]);
            for(int nnni=nTemptNum;nnni<vObjects.size();nnni++)
            {
                vSplitObjects.push_back(*vObjects[nnni]);
            }
              
              SaveDrawObj(vSplitObjects,nWidth, nHeight, nLayers,CylinderPool,pPredPool,"./MergeSplit/meg_");


            
          //   SaveLooseConn(&objectI, &objectJ,nWidth, nHeight, nLayers,CylinderPool,pPredPool,"./looseCon/los_", nConInf);
             std::cout<<"ConnCode="<<nConInf<<'\n';
             cv::Mat matz=cv::Mat::zeros(nHeight,nWidth,CV_8UC1);
             cv::imshow("loose_wait",matz );
             cv::waitKey(0);
             cv::destroyWindow("loose_wait");
           }
#endif

  
              if(vObjects[i]->m_bFrontM==false)
                  break;
        
         }
            EraseObjectLocalCONN(vObjects[i],ppuch3DConLocal,nWidth, nHeight,nLayers,CylinderPool,pPredPool);
        }
        
        
        if(vObjects[i]->m_bBackM==true)
        {
            
        //    std::cout<<"-----------------------------------------";
            
            AddObjectLocalCONN(vObjects[i],ppuch3DConLocal,nWidth, nHeight,nLayers,CylinderPool,pPredPool,4);
            for(int j=i+1; j<nNumObj; j++)
            {
                
                
                if(std::find(vObjects[i]->m_vnNeighborIDs.begin(), vObjects[i]->m_vnNeighborIDs.end(), j) == vObjects[i]->m_vnNeighborIDs.end())
                {
                    continue;
                }
                
                
                if(vObjects[j]->m_bBackM==false&&vObjects[j]->m_bFrontM==false)
                {
                    continue;
                }
                int nConInf = CheckConnObjIJ(vObjects[i], vObjects[j], nWidth, nHeight, nLayers,CylinderPool,pPredPool,false);
                
                
                
                if(nConInf==0)
                    continue;
                
      //          std::cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&";
                vObjects[i]->m_bBackM=false;
                
  
#if TESTMERGE==1
                int nTemptNum = vObjects.size();
                Que3D objectI,objectJ;
                objectI = *vObjects[i];
                objectJ = *vObjects[j];
#endif
                bool bIsRF=true,bIsRB=true;
                
  
                if(nConInf==3)
                {
                    if(vObjects[j]->m_bFrontM==false)
                        break;
               //     std::cout<<"I B rate is "<<vObjects[i]->m_dBConnR<<"\n";
               //     std::cout<<"J F rate is "<<vObjects[j]->m_dFConnR<<"\n";
                    
                    TDPoint ptI,ptJ,ptMid;
                    
                    ptMid.m_dX= (vObjects[i]->m_dX+(*pPredPool)[vObjects[i]->m_mark].m_dBCenX+vObjects[j]->m_dX+(*pPredPool)[vObjects[j]->m_mark].m_dFCenX)/2.0;
                    
                    ptMid.m_dY= (vObjects[i]->m_dY+(*pPredPool)[vObjects[i]->m_mark].m_dBCenY+vObjects[j]->m_dY+(*pPredPool)[vObjects[j]->m_mark].m_dFCenY)/2.0;
                    
                    ptMid.m_dZ= (vObjects[i]->m_dZ+(*pPredPool)[vObjects[i]->m_mark].m_dBCenZ+vObjects[j]->m_dZ+(*pPredPool)[vObjects[j]->m_mark].m_dFCenZ)/2.0;
                    
                    ptI.m_dX = (*pPredPool)[vObjects[i]->m_mark].m_dFCenX+vObjects[i]->m_dX;
                    ptI.m_dY = (*pPredPool)[vObjects[i]->m_mark].m_dFCenY+vObjects[i]->m_dY;
                    ptI.m_dZ = (*pPredPool)[vObjects[i]->m_mark].m_dFCenZ+vObjects[i]->m_dZ;
                    
                    
                    ptJ.m_dX = (*pPredPool)[vObjects[j]->m_mark].m_dBCenX+vObjects[j]->m_dX;
                    ptJ.m_dY = (*pPredPool)[vObjects[j]->m_mark].m_dBCenY+vObjects[j]->m_dY;
                    ptJ.m_dZ = (*pPredPool)[vObjects[j]->m_mark].m_dBCenZ+vObjects[j]->m_dZ;
                    
                    TDPoint ptIC=ptI;
                    TDPoint ptJC=ptJ;
                    TDPoint ptMidC=ptMid;
                    int nRT=nRTimes;
                    while(nRT-->0)
                    {
                        
                        ptI = ptIC; ptJ=ptJC; ptMid=ptMidC;
                        RandomX(ptI,1);
                        RandomX(ptJ,1);
                        RandomX(ptMid,2);
                        
                        bIsRF=ReplaceObjects(i, ptMid, ptI, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                        
                        bIsRB=ReplaceObjects(j, ptMid, ptJ, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                        if(bIsRF==true||bIsRB==true)
                            break;
                    }
                    
                    
                    vObjects[i]->m_bBackM=false;
                    vObjects[j]->m_bFrontM=false;
                    
                
                }
                
                if(nConInf==4)
                {
                    if(vObjects[j]->m_bBackM==false)
                        break;
                    
               //     std::cout<<"I B rate is "<<vObjects[i]->m_dBConnR<<"\n";
               //     std::cout<<"J B rate is "<<vObjects[j]->m_dBConnR<<"\n";
                    
                    TDPoint ptI,ptJ,ptMid;
                    
                    ptMid.m_dX= (vObjects[i]->m_dX+(*pPredPool)[vObjects[i]->m_mark].m_dBCenX+vObjects[j]->m_dX+(*pPredPool)[vObjects[j]->m_mark].m_dBCenX)/2.0;
                    
                    ptMid.m_dY= (vObjects[i]->m_dY+(*pPredPool)[vObjects[i]->m_mark].m_dBCenY+vObjects[j]->m_dY+(*pPredPool)[vObjects[j]->m_mark].m_dBCenY)/2.0;
                    
                    ptMid.m_dZ= (vObjects[i]->m_dZ+(*pPredPool)[vObjects[i]->m_mark].m_dBCenZ+vObjects[j]->m_dZ+(*pPredPool)[vObjects[j]->m_mark].m_dBCenZ)/2.0;
                    
                    ptI.m_dX = (*pPredPool)[vObjects[i]->m_mark].m_dFCenX+vObjects[i]->m_dX;
                    ptI.m_dY = (*pPredPool)[vObjects[i]->m_mark].m_dFCenY+vObjects[i]->m_dY;
                    ptI.m_dZ = (*pPredPool)[vObjects[i]->m_mark].m_dFCenZ+vObjects[i]->m_dZ;
                    
                    
                    ptJ.m_dX = (*pPredPool)[vObjects[j]->m_mark].m_dFCenX+vObjects[j]->m_dX;
                    ptJ.m_dY = (*pPredPool)[vObjects[j]->m_mark].m_dFCenY+vObjects[j]->m_dY;
                    ptJ.m_dZ = (*pPredPool)[vObjects[j]->m_mark].m_dFCenZ+vObjects[j]->m_dZ;
                    
                    TDPoint ptIC=ptI;
                    TDPoint ptJC=ptJ;
                    TDPoint ptMidC=ptMid;
                    int nRT=nRTimes;
                    while(nRT-->0)
                    {
                        
                        ptI = ptIC; ptJ=ptJC; ptMid=ptMidC;
                        RandomX(ptI,1);
                        RandomX(ptJ,1);
                        RandomX(ptMid,2);
                        
                        bIsRF=ReplaceObjects(i, ptMid, ptI, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                        
                        bIsRB=ReplaceObjects(j, ptMid, ptJ, vObjects, ppuch3DGray,  ppuchSeg3D, ppn3DConnect,ppuch3DConLocal,ppnOverlap,nWidth, nHeight, nLayers,  CylinderPool, pPredPool, conP);
                        if(bIsRF==true||bIsRB==true)
                            break;
                    }
                    
                    
                    vObjects[i]->m_bBackM=false;
                    vObjects[j]->m_bBackM=false;
                    
                }
                
  
#if TESTMERGE==1
           if((bIsRF==false||bIsRB==false)&&nGlobalCnt>100)
            {
                OutputObjInfo(vObjects[i], "aftermergeI");
                OutputObjInfo(&objectI, "beforemergeI");
                std::vector<Que3D> vSplitObjects;
                vSplitObjects.push_back(*vObjects[i]);
                vSplitObjects.push_back(*vObjects[j]);
                for(int nnni=nTemptNum;nnni<vObjects.size();nnni++)
                {
                    vSplitObjects.push_back(*vObjects[nnni]);
                }
                
                SaveDrawObj(vSplitObjects,nWidth, nHeight, nLayers,CylinderPool,pPredPool,"./MergeSplit/meg_");
                
                
                
           //     SaveLooseConn(&objectI, &objectJ,nWidth, nHeight, nLayers,CylinderPool,pPredPool,"./looseCon/los_", nConInf);
                std::cout<<"ConnCode="<<nConInf<<'\n';
                cv::Mat matz=cv::Mat::zeros(nHeight,nWidth,CV_8UC1);
                cv::imshow("loose_wait",matz );
                cv::waitKey(0);
                cv::destroyWindow("loose_wait");
            }
#endif
  
                
                
                if(vObjects[i]->m_bBackM==false)
                    break;
                
            }
  
            EraseObjectLocalCONN(vObjects[i],ppuch3DConLocal,nWidth, nHeight,nLayers,CylinderPool,pPredPool);
        }

    }
        

    
    
    
}




void AddObjectLocalCONN(Que3D* pObject,unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,int ntype)
{
    
    std::vector<TPoint> vTPoints = CylinderPool[pObject->m_mark];
    
    int nCenX = round(pObject->m_dX);
    int nCenY = round(pObject->m_dY);
    int nCenZ = round(pObject->m_dZ);
    
    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        if(vTPoints[i].m_nType==ntype)
        {
            
            ppuch3DConnect[nz][nx+ny*nWidth]=ntype;
            
        }
        
        
    }
    
}

void EraseObjectLocalCONN(Que3D* pObject,unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    std::vector<TPoint> vTPoints = CylinderPool[pObject->m_mark];
    
    int nCenX = round(pObject->m_dX);
    int nCenY = round(pObject->m_dY);
    int nCenZ = round(pObject->m_dZ);
    
    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        if(vTPoints[i].m_nType==3)
        {
            
            ppuch3DConnect[nz][nx+ny*nWidth]=0;
            
        }
        else if(vTPoints[i].m_nType==4)
        {
            ppuch3DConnect[nz][nx+ny*nWidth]=0;
            
        }
        
    }
}
int MCheckConnObjIJ(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bCheckF,int nDistance)
{
    nDistance=101;
    //Check radius
    if(ABS(pObjectI->m_mark.m_nR-pObjectJ->m_mark.m_nR)>2)
    return 0;
    
    
    int nAngleDY = std::min(abs(pObjectI->m_mark.m_nthetaY-pObjectJ->m_mark.m_nthetaY),DEGREE_Y-abs((pObjectI->m_mark.m_nthetaY-pObjectJ->m_mark.m_nthetaY)));
    
    int nAngleDZ = std::min(abs(pObjectI->m_mark.m_nthetaZ-pObjectJ->m_mark.m_nthetaZ),DEGREE_Z-abs((pObjectI->m_mark.m_nthetaZ-pObjectJ->m_mark.m_nthetaZ)));
    //Check obj angle
    
  //  bool bSpecial = false;
    
    //Special case
    
 //   bSpecial= CheckSpecialAngle(pObjectI->m_mark.m_nthetaY,pObjectI->m_mark.m_nthetaZ,pObjectJ->m_mark.m_nthetaY,pObjectJ->m_mark.m_nthetaZ);
    
    
    if(nAngleDY>20||nAngleDZ>20)//&&bSpecial==false)
    return 0;
    
    
    //Check connect line angle
    
    double dx = pObjectJ->m_dX - pObjectI->m_dX;
    double dy = pObjectJ->m_dY - pObjectI->m_dY;
    double dz = pObjectJ->m_dZ - pObjectI->m_dZ;
    
    
    
    
    int nAngleCLY;
    int nAngleCLZ ;
    
    ConvertFromXYZ2DYZ_N(nAngleCLY, nAngleCLZ,dx,dy,dz);
    
    int nAngleDY1= std::min(abs(pObjectI->m_mark.m_nthetaY-nAngleCLY),DEGREE_Y-abs((pObjectI->m_mark.m_nthetaY-nAngleCLY)));
    
    int nAngleDZ1= std::min(abs(pObjectI->m_mark.m_nthetaZ-nAngleCLZ),DEGREE_Z-abs((pObjectI->m_mark.m_nthetaZ-nAngleCLZ)));
    
    
    int nAngleDY2= std::min(abs(pObjectJ->m_mark.m_nthetaY-nAngleCLY),DEGREE_Y-abs((pObjectJ->m_mark.m_nthetaY-nAngleCLY)));
    
    int nAngleDZ2= std::min(abs(pObjectJ->m_mark.m_nthetaZ-nAngleCLZ),DEGREE_Z-abs((pObjectJ->m_mark.m_nthetaZ-nAngleCLZ)));
    
    
 //   bool bSpecialI= CheckSpecialAngle(pObjectI->m_mark.m_nthetaY,pObjectI->m_mark.m_nthetaZ,nAngleCLY,nAngleCLZ);
    
 //   bool bSpecialJ= CheckSpecialAngle(pObjectJ->m_mark.m_nthetaY,pObjectJ->m_mark.m_nthetaZ,nAngleCLY,nAngleCLZ);
    
 //   bSpecial = (bSpecialI||bSpecialJ);
    
    
    if(MIN(nAngleDY1,nAngleDY2)>20||MIN(nAngleDZ1,nAngleDZ2)>20)//&&bSpecial==false)
    return 0;
    
    double dIx,dIy,dIz;
    
    if(bCheckF)
    {
        dIx = (*pPredPool)[pObjectI->m_mark].m_dFCenX+pObjectI->m_dX;
        dIy = (*pPredPool)[pObjectI->m_mark].m_dFCenY+pObjectI->m_dY;
        dIz = (*pPredPool)[pObjectI->m_mark].m_dFCenZ+pObjectI->m_dZ;
    }
    else
    {
        dIx = (*pPredPool)[pObjectI->m_mark].m_dBCenX+pObjectI->m_dX;
        dIy = (*pPredPool)[pObjectI->m_mark].m_dBCenY+pObjectI->m_dY;
        dIz = (*pPredPool)[pObjectI->m_mark].m_dBCenZ+pObjectI->m_dZ;
    }
    
    
    double dJFx,dJFy,dJFz,dJBx,dJBy,dJBz;
    
    dJFx = (*pPredPool)[pObjectJ->m_mark].m_dFCenX+pObjectJ->m_dX;
    dJFy = (*pPredPool)[pObjectJ->m_mark].m_dFCenY+pObjectJ->m_dY;
    dJFz = (*pPredPool)[pObjectJ->m_mark].m_dFCenZ+pObjectJ->m_dZ;
    
    dJBx = (*pPredPool)[pObjectJ->m_mark].m_dBCenX+pObjectJ->m_dX;
    dJBy = (*pPredPool)[pObjectJ->m_mark].m_dBCenY+pObjectJ->m_dY;
    dJBz = (*pPredPool)[pObjectJ->m_mark].m_dBCenZ+pObjectJ->m_dZ;
    
    
    double dIF = SqDist2Points3D(dIx,dIy,dIz,dJFx,dJFy,dJFz);
    double dIB = SqDist2Points3D(dIx,dIy,dIz,dJBx,dJBy,dJBz);
    
    
    double dDd = pObjectI->m_mark.m_nR;
    // dDd=MIN(dDd*dDd*4,100);
    
    if(dIF<nDistance&&dIF<=dIB)
    {
        if(bCheckF)
        {
            return 1;
        }
        else
        {
            return 3;
        }
    }
    
    
    if(dIB<nDistance&&dIB<=dIF)
    {
        if(bCheckF)
        {
            return 2;
        }
        else
        {
            return 4;
        }
    }
    
    return 0;
    
    
}
int CheckConnObjIJ(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bCheckF,int nDistance)
{
    
    
    if(ABS(pObjectI->m_mark.m_nR-pObjectJ->m_mark.m_nR)>2)
        return 0;
    
    int nAngleDZ = abs(pObjectI->m_mark.m_nthetaZ-pObjectJ->m_mark.m_nthetaZ);
    
    int nAngleDY = abs(pObjectI->m_mark.m_nthetaY-pObjectJ->m_mark.m_nthetaY);
    
    
    
    //Check obj angle
    
    bool bSpecial = false;
    
    //Special case
    
    bSpecial= CheckSpecialAngleBigD(pObjectI->m_mark.m_nthetaY,pObjectI->m_mark.m_nthetaZ,pObjectJ->m_mark.m_nthetaY,pObjectJ->m_mark.m_nthetaZ);
    
    
    if((nAngleDY>18||nAngleDZ>18)&&bSpecial==false)
        return 0;
    
    
    //Check connect line angle
    
    double dx = pObjectJ->m_dX - pObjectI->m_dX;
    double dy = pObjectJ->m_dY - pObjectI->m_dY;
    double dz = pObjectJ->m_dZ - pObjectI->m_dZ;
    
    
    
    
    int nAngleCLY;
    int nAngleCLZ ;
    
    ConvertFromXYZ2DYZ_N(nAngleCLY, nAngleCLZ,dx,dy,dz);
    
    
    int nAngleDZ1= abs(pObjectI->m_mark.m_nthetaZ-nAngleCLZ);
    int nAngleDY1= abs(pObjectI->m_mark.m_nthetaY-nAngleCLY);
    
    
    
    int nAngleDY2= abs(pObjectJ->m_mark.m_nthetaY-nAngleCLY);
    
    int nAngleDZ2= abs(pObjectJ->m_mark.m_nthetaZ-nAngleCLZ);
    
    
    
    bool bSpecialI= CheckSpecialAngleBigD(pObjectI->m_mark.m_nthetaY,pObjectI->m_mark.m_nthetaZ,nAngleCLY,nAngleCLZ);
    
    bool bSpecialJ= CheckSpecialAngleBigD(pObjectJ->m_mark.m_nthetaY,pObjectJ->m_mark.m_nthetaZ,nAngleCLY,nAngleCLZ);
    
    bSpecial = (bSpecialI||bSpecialJ);
    
    
    if((MIN(nAngleDY1,nAngleDY2)>15||MIN(nAngleDZ1,nAngleDZ2)>15)&&bSpecial==false)
        return 0;
    
    

    


    

    
    double dIx,dIy,dIz;
    
    if(bCheckF)
    {
        dIx = (*pPredPool)[pObjectI->m_mark].m_dFCenX+pObjectI->m_dX;
        dIy = (*pPredPool)[pObjectI->m_mark].m_dFCenY+pObjectI->m_dY;
        dIz = (*pPredPool)[pObjectI->m_mark].m_dFCenZ+pObjectI->m_dZ;
    }
    else
    {
        dIx = (*pPredPool)[pObjectI->m_mark].m_dBCenX+pObjectI->m_dX;
        dIy = (*pPredPool)[pObjectI->m_mark].m_dBCenY+pObjectI->m_dY;
        dIz = (*pPredPool)[pObjectI->m_mark].m_dBCenZ+pObjectI->m_dZ;
    }
    
    
    double dJFx,dJFy,dJFz,dJBx,dJBy,dJBz;
    
    dJFx = (*pPredPool)[pObjectJ->m_mark].m_dFCenX+pObjectJ->m_dX;
    dJFy = (*pPredPool)[pObjectJ->m_mark].m_dFCenY+pObjectJ->m_dY;
    dJFz = (*pPredPool)[pObjectJ->m_mark].m_dFCenZ+pObjectJ->m_dZ;
    
    dJBx = (*pPredPool)[pObjectJ->m_mark].m_dBCenX+pObjectJ->m_dX;
    dJBy = (*pPredPool)[pObjectJ->m_mark].m_dBCenY+pObjectJ->m_dY;
    dJBz = (*pPredPool)[pObjectJ->m_mark].m_dBCenZ+pObjectJ->m_dZ;
    
    
    double dIF = SqDist2Points3D(dIx,dIy,dIz,dJFx,dJFy,dJFz);
    double dIB = SqDist2Points3D(dIx,dIy,dIz,dJBx,dJBy,dJBz);
    
    
    double dDd = pObjectI->m_mark.m_nR;
   // dDd=MIN(dDd*dDd*4,100);
    
    if(pObjectJ->m_bFrontM==true&&dIF<nDistance&&dIF<=dIB)
    {
        if(bCheckF)
        {
            return 1;
        }
        else
        {
            return 3;
        }
    }
    

    if(pObjectJ->m_bBackM==true&&dIB<nDistance&&dIB<=dIF)
    {
        if(bCheckF)
        {
            return 2;
        }
        else
        {
            return 4;
        }
    }
    
    return 0;
    
    
    
    
    
  /*
    std::vector<TPoint> vTPoints = CylinderPool[pObjectJ->m_mark];
    
    int nCenX = round(pObjectJ->m_dX);
    int nCenY = round(pObjectJ->m_dY);
    int nCenZ = round(pObjectJ->m_dZ);
    
    
    int nFFCount =0, nFBCount=0, nBFCount=0, nBBCount=0;
    int nCountF =0; int nCountB=0;
    for(int i=0; i<vTPoints.size(); i++)
    {
        
        int nx = vTPoints[i].m_nX+nCenX;
        int ny = vTPoints[i].m_nY+nCenY;
        int nz = vTPoints[i].m_nZ+nCenZ;
        
        if(nx<0||nx>=nWidth||
           ny<0||ny>=nHeight||
           nz<0||nz>=nLayers)
        {
            continue;
        }
        
        if(vTPoints[i].m_nType==3)
        {
       
            if(ppuch3DConnect[nz][nx+ny*nWidth]==3)
            {
                nFFCount++;
            }
            if(ppuch3DConnect[nz][nx+ny*nWidth]==4)
            {
                nBFCount++;
            }
            
        }
        else if(vTPoints[i].m_nType==4)
        {
         
            if(ppuch3DConnect[nz][nx+ny*nWidth]==3)
            {
                nFBCount++;
            }
            if(ppuch3DConnect[nz][nx+ny*nWidth]==4)
            {
                nBBCount++;
            }
        }
        
    }
    
    if(nFFCount==0&&nFBCount==0&&nBFCount==0&&nBBCount==0)
    {
        return 0;
    }
    
    
    int nMaxCnt = MAX(MAX(nFFCount,nFBCount),MAX(nBFCount,nBBCount));
    
    if(nFFCount==nMaxCnt)
        return 1;
    if(nFBCount==nMaxCnt)
        return 2;
    if(nBFCount==nMaxCnt)
        return 3;
    if(nBBCount==nMaxCnt)
        return 4;
    
    return 0;*/
}


void ConvertFromXYZ2DYZ(double& dTY, double& dTZ,double dx,double dy, double dz)
{
    
    //normalizeing
    
   // dx =-1; dy=-1;dz=-1;
    double dnx, dny,dnz;
    
    double dLen = sqrt(dx*dx+dy*dy+dz*dz);
    
    dnx= dx/dLen;
    dny= dy/dLen;
    dnz= dz/dLen;
    

    
 //   dTZ = atan(dny/dnx);
    
    dTZ = cvFastArctan(dny, dnx);
    double dThetaZ = dTZ/180*PI_L;
    double drx =(cos(dThetaZ)*dnx+sin(dThetaZ)*dny);
    double dry = -sin(dThetaZ)*dnx+cos(dThetaZ)*dny;
    

  //  dTY = atan(drx/dnz);
    
    dTY = cvFastArctan(drx, dnz);
    
    double dThetaY =dTY/180*PI_L;
    
    
    double drz = dnz*cos(dThetaY)+sin(dThetaY)*drx;
 
    drx =-sin(dThetaY)*dnz+cos(dThetaY)*drx;
    
    
    std::cout<<"drx= "<<drx<<" dry="<< dry<<" drz="<<drz<<'\n';
    int nx=0,ny=0,nz=1;
    
    double dmx = double(nx)*cos(dThetaY)+double(nz)*sin(dThetaY);
    double dmy = ny;
    double dmz = -double(nx)*sin(dThetaY)+double(nz)*cos(dThetaY);
    
    
   // std::cout<<"dmx="<<dmx<<" dmy="<< dmy<<" dmz="<<dmz<<'\n';
  
    
    
    double dox = cos(dThetaZ)*dmx-sin(dThetaZ)*dmy;
    double doy = sin(dThetaZ)*dmx+cos(dThetaZ)*dmy;
    double doz = dmz;
    

    
    

    
    std::cout<<"Origin x= "<<dnx<<" y= "<<dny<<" z= "<<dnz<<'\n';
    
    std::cout<<"Recovery x = "<<dox<<"  y = "<<doy<<"  z = "<<doz<<'\n';
    
    
    
    
   // dTZ=dTZ*180/PI_L;
   /// dTY=dTY*180/PI_L;
    

    
}


void ConvertFromXYZ2DYZ_N(int& nTY, int& nTZ,double dx,double dy, double dz)
{
    
    //normalizeing
   if(dy<0)
   {
       dx=-dx;
       dz=-dz;
       dy=-dy;
   }
    
    
    
    
    double dTY,dTZ;
    double dnx, dny,dnz;
    
    double dLen = sqrt(dx*dx+dy*dy+dz*dz);
    
    dnx= dx/dLen;
    dny= dy/dLen;
    dnz= dz/dLen;
    
    
    
    //   dTZ = atan(dny/dnx);
    
    dTZ = cvFastArctan(dny, dnx);
    double dThetaZ = dTZ/180*PI_L;
    double drx =(cos(dThetaZ)*dnx+sin(dThetaZ)*dny);
    double dry = -sin(dThetaZ)*dnx+cos(dThetaZ)*dny;

    
    dTY = cvFastArctan(drx, dnz);
    
 
    nTY=  int(dTY)%180/2;
    nTZ = int(dTZ)%180/2;
    
    
}
void SplitByPoints(std::vector<Que3D>& vSplitObjects,int nR, TDPoint ptMid, TDPoint ptOrigin,unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
 /*
    for(int i=-1;i<=2;i++)
        for(int j=-1;j<=2;j++)
        {
            Set3DPixel( int(ptMid.m_dX+i+0.5),  int(ptMid.m_dY+j+0.5),  int(0.5+ptMid.m_dZ), nWidth, ppuch3DConnect, 255);
            Set3DPixel(int(0.5+ptOrigin.m_dX+i), int(0.5+ptOrigin.m_dY+j), int(0.5+ptOrigin.m_dZ), nWidth, ppuch3DConnect, 255);
            
        }
    
  */
    
    
    double dH = Dist2Points3D(ptMid.m_dX, ptMid.m_dY, ptMid.m_dZ, ptOrigin.m_dX, ptOrigin.m_dY, ptOrigin.m_dZ);
    dH=dH/2.0;
    
    
    if(dH>2*NMAXH)
    {
        return;
    }
    if(dH>NMAXH) //Then split
    {
        
   //     std::cout<<"Splitting!\n";
        double dx = (ptMid.m_dX+ptOrigin.m_dX)/2.0;
        double dy = (ptMid.m_dY+ptOrigin.m_dY)/2.0;
        double dz = (ptMid.m_dZ+ptOrigin.m_dZ)/2.0;
     //   int nRR = nR-1+random2()*2+0.5;
     //   nR= MAX(NMINR,MIN(NMAXR, nRR));
        
        dH = dH/2;
        
        int nThetaY,nThetaZ;
        
        ConvertFromXYZ2DYZ_N(nThetaY, nThetaZ, ptMid.m_dX-ptOrigin.m_dX, ptMid.m_dY-ptOrigin.m_dY, ptMid.m_dZ-ptOrigin.m_dZ);
        
        CMark mark(nR,dH+0.5,nThetaY,nThetaZ);
        
        double dx0 = (dx+ptOrigin.m_dX)/2.0;
        double dy0 = (dy+ptOrigin.m_dY)/2.0;
        double dz0 = (dz+ptOrigin.m_dZ)/2.0;
        
        
        
        double dx1 = (ptMid.m_dX+ dx)/2.0;
        double dy1 = (ptMid.m_dY+ dy)/2.0;
        double dz1 = (ptMid.m_dZ+ dz)/2.0;
    
        Que3D qObject(dx0,dy0,dz0,mark,(*pPredPool)[mark]);
        vSplitObjects.push_back(qObject);
        
        
        Que3D qObject1(dx1,dy1,dz1,mark,(*pPredPool)[mark]);
        vSplitObjects.push_back(qObject1);
     /*
        DrawCylinder(dx0+0.5,dy0+0.5, dz0+0.5, mark, ppuch3DConnect, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
        
        DrawCylinder(dx1+0.5,dy1+0.5, dz1+0.5, mark, ppuch3DConnect, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
        */
    //    std::cout<<"dH is "<<dH<<'\n';
        
    }
    else
    {
  //      std::cout<<"No splitting!\n";
        double dx = (ptMid.m_dX+ptOrigin.m_dX)/2.0;
        double dy = (ptMid.m_dY+ptOrigin.m_dY)/2.0;
        double dz = (ptMid.m_dZ+ptOrigin.m_dZ)/2.0;
       
     //   int nRR = nR-1+random2()*2+0.5;
     //   nR= MAX(NMINR,MIN(NMAXR, nRR));
       
        int nThetaY,nThetaZ;
        
        ConvertFromXYZ2DYZ_N(nThetaY, nThetaZ, ptMid.m_dX-ptOrigin.m_dX, ptMid.m_dY-ptOrigin.m_dY, ptMid.m_dZ-ptOrigin.m_dZ);
        
        
        CMark mark(nR,dH+0.5,nThetaY,nThetaZ);

        Que3D qObject(dx,dy,dz,mark,(*pPredPool)[mark]);
        vSplitObjects.push_back(qObject);
     /*
        
        DrawCylinder(dx+0.5,dy+0.5, dz+0.5, mark, ppuch3DConnect, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
      */
   //     std::cout<<"dH is "<<dH<<'\n';
        
    }
    
    
}
void DrawObjByGiven2Points(int x1,int y1, int z1,
                           int x2,int y2, int z2,
                           unsigned char** ppuch3DConnect,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
    if(IsInnerPoint(x1, y1, z1, nWidth, nHeight, nLayers)==false||
       IsInnerPoint(x2, y2, z2, nWidth, nHeight, nLayers)==false)
    {
        std::cout<<"Out of bound!\n";
    }
    std::cout<<"drawing!\n";
    double dx = (x2+x1)/2.0;
    double dy = (y2+y1)/2.0;
    double dz = (z2+z1)/2.0;
    
    double dH = Dist2Points3D(x1, y1, z1, dx, dy, dz);
    double dR = 5;
    Que3D object;
    object.m_dX=dx;
    object.m_dY=dy;
    
    int nThetaY,nThetaZ;
    
    ConvertFromXYZ2DYZ_N(nThetaY, nThetaZ, x2-x1, y2-y1, z2-z1);
    
    CMark mark(dR,dH,nThetaY,nThetaZ);
    
    object.m_mark =mark;
    
    OutputObjInfo(&object,"infor");
    
    for(int i=-1;i<=1;i++)
        for(int j=-1;j<=1;j++)
        {
            Set3DPixel(x1+i, y1+j, z1, nWidth, ppuch3DConnect, 255);
            Set3DPixel(x2+i, y2+j, z2, nWidth, ppuch3DConnect, 255);
        
        }
    
    DrawCylinder(dx+0.5,dy+0.5, dz+0.5, mark, ppuch3DConnect, nWidth, nHeight, nLayers, CylinderPool, pPredPool);
    
}
void SaveDrawObj(std::vector<Que3D>& vSplitObjects,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string strSaveAdd)
{

    sRGB rgb[10]={ sRGB(0,255,0),sRGB(0,0,255),sRGB(0,255,255),sRGB(255,255,0),sRGB(255,0,0),sRGB(255,255,255),sRGB(255,0,255),sRGB(128,180,30),sRGB(20,60,180),sRGB(180,40,110)};
    cv::Mat matColorSeq[MAXSequence];
    unsigned char* ppuch3DRGB[MAXSequence];
    CreateVoid3DColorImages(nLayers, nWidth,nHeight, matColorSeq, ppuch3DRGB);
    
    for(int i=0; i<vSplitObjects.size();i++)
    {
        DrawCylinderColor(round(vSplitObjects[i].m_dX), round(vSplitObjects[i].m_dY),round(vSplitObjects[i].m_dZ),vSplitObjects[i].m_mark,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,pPredPool,rgb[i%10]);
    }
    
    Save3DImages(strSaveAdd,nLayers, matColorSeq);
    
    
}
void SaveLooseConn(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,std::string strSaveAdd,int nConnValue)
{

    if(nConnValue==0)
        return;
    
    
    cv::Mat matColorSeq[MAXSequence];
    unsigned char* ppuch3DRGB[MAXSequence];
    CreateVoid3DColorImages(nLayers, nWidth,nHeight, matColorSeq, ppuch3DRGB);
    
    
     DrawCylinderColor(round(pObjectI->m_dX), round(pObjectI->m_dY),round(pObjectI->m_dZ),pObjectI->m_mark,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,pPredPool,sRGB(255,0,0));
    
     DrawCylinderColor(round(pObjectJ->m_dX), round(pObjectJ->m_dY),round(pObjectJ->m_dZ),pObjectJ->m_mark,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,pPredPool,sRGB(0,255,0));
    
    int ntype1=1,ntype2=1;
    
    
    switch (nConnValue) {
        case 1:
            ntype1=3;ntype2=3;break;
        case 2:
            ntype1=3;ntype2=4;break;
        case 3:
            ntype1=4;ntype2=3;break;
        case 4:
            ntype1=4;ntype2=4;break;
        default:
            return;
            break;
    }
    
    /*
    DrawCylinderColor(round(pObjectI->m_dX), round(pObjectI->m_dY),round(pObjectI->m_dZ),pObjectI->m_mark,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,pPredPool,sRGB(0,0,255),ntype1);
    
    DrawCylinderColor(round(pObjectJ->m_dX), round(pObjectJ->m_dY),round(pObjectJ->m_dZ),pObjectJ->m_mark,ppuch3DRGB,nWidth, nHeight, nLayers,CylinderPool,pPredPool,sRGB(255,255,255),ntype2);
    */
    
    
    Save3DImages(strSaveAdd,nLayers, matColorSeq);
    
}
void MergeShortToLongNew(std::vector<Que3D*>& vp3DStartLong, std::vector<Que3D*>& vp3DQues,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,unsigned char** ppuch3DConLocal)
{
    
        vp3DStartLong.clear();
        CalculateNeigborID(vp3DQues, pow(NMAXH*4+1,2));
    
        for(int i=0; i<vp3DQues.size(); i ++)
        {
            if(vp3DQues[i]->m_bIsMerged)
                continue;
            
            Que3D* pFront =   vp3DQues[i];
            Que3D* pBack  =   vp3DQues[i];
            
            vp3DQues[i]->m_bIsMerged = true;
            
            bool bMoreFront= true;
            bool bMoreBack = true;
            bool bFIsCheckFront =true;
            bool bBIsCheckFront =false;
            int ncnt=0;
            while(bMoreBack||bMoreFront)
            {
                int nF_NearJ=-1;
                int nB_NearJ=-1;
                int nFConCode =0;
                int nBConCode =0;
                int nTCodeF=0;
                int nTCodeB=0;
                double dNearestDistF=10000;
                double dNearestDistB=10000;
                double dDistF =0,dDistB=0;
                
                
                if(bMoreFront)
                {
                    bMoreFront=false;
                    int nNbT = pFront->m_vnNeighborIDs.size();
                    for(int j=0; j<nNbT; j++)
                    {
                         if(vp3DQues[pFront->m_vnNeighborIDs[j]]->m_bIsMerged)
                             continue;
                             
                         nTCodeF =CMCheckConnObjIJ(pFront, vp3DQues[pFront->m_vnNeighborIDs[j]], nWidth, nHeight, nLayers,CylinderPool,pPredPool,bFIsCheckFront,ppuch3DConLocal);
                         if(nTCodeF!=0)
                         {

                             dDistF = DistP2P(pFront, vp3DQues[pFront->m_vnNeighborIDs[j]], nTCodeF,pPredPool);
                             
                             if(dDistF<dNearestDistF)
                             {
                                 dNearestDistF = dDistF;
                                 nF_NearJ = pFront->m_vnNeighborIDs[j];
                                 nFConCode = nTCodeF;
                             }
                         }
                    }
                    if(nF_NearJ==-1)
                    {
                        bMoreFront=false;
                        pFront->m_pQueHead=NULL;
                    }
                    else
                    {
                        bMoreFront=true;
                        vp3DQues[nF_NearJ]->m_bIsMerged=true;
                        vp3DQues[nF_NearJ]->m_pQueTail = pFront;
                        pFront->m_pQueHead = vp3DQues[nF_NearJ];

                        pFront = vp3DQues[nF_NearJ];
                        
                        if(nFConCode == 1||nFConCode == 3)
                        {
                            bFIsCheckFront = false;
                        }
                        else
                        {
                            bFIsCheckFront = true;
                        }
                    }
                    
                }
                
                
                if(bMoreBack)
                {
                     bMoreBack=false;
                     int nNbT = pBack->m_vnNeighborIDs.size();
                     for(int j=0; j<nNbT; j++)
                     {
                         
                         if(vp3DQues[pBack->m_vnNeighborIDs[j]]->m_bIsMerged)
                             continue;
                         
                         nTCodeB =CMCheckConnObjIJ(pBack, vp3DQues[pBack->m_vnNeighborIDs[j]], nWidth, nHeight, nLayers,CylinderPool,pPredPool,bBIsCheckFront,ppuch3DConLocal);
                         if(nTCodeB!=0)
                         {
                          
                             dDistB = DistP2P(pBack, vp3DQues[pBack->m_vnNeighborIDs[j]], nTCodeB,pPredPool);
                            
                             if(dDistB<dNearestDistB)
                             {
                                 dNearestDistB = dDistB;
                                 nB_NearJ = pBack->m_vnNeighborIDs[j];
                                 nBConCode = nTCodeB;
                              }
                           }
                       }
                    
                    if(nB_NearJ==-1)
                    {
                        bMoreBack=false;
                        pBack->m_pQueTail=NULL;
                    }
                    else
                    {
                        
                        bMoreBack=true;
                        pBack->m_pQueTail = vp3DQues[nB_NearJ];
                        vp3DQues[nB_NearJ]->m_pQueHead = pBack;
                        vp3DQues[nB_NearJ]->m_bIsMerged=true;
                        pBack = vp3DQues[nB_NearJ];
                        
                        if(nBConCode == 1||nBConCode == 3)
                        {
                            bBIsCheckFront = false;
                        }
                        else
                        {
                            bBIsCheckFront = true;
                        }
                    }
                    
                    
                    
                    
                }

                
            }
            vp3DStartLong.push_back(pBack);
            
            
            
        }
    
    
    
}
void MergeShortToLong(std::vector<Que3D*>& vp3DStartLong, std::vector<Que3D*>& vp3DQues,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,unsigned char** ppuch3DConLocal)
{
  


    for(int i=0; i<vp3DQues.size(); i ++)
    {
        
    
        
        if(vp3DQues[i]->m_bIsMerged==false)
        {

            
            Que3D* pFront =   vp3DQues[i];
            Que3D* pBack  =   vp3DQues[i];
            
            vp3DQues[i]->m_bIsMerged = true;
            
            
            bool bMoreFront= true;
            bool bMoreBack = true;
            bool bFIsCheckFront =true;
            bool bBIsCheckFront =false;
            int ncnt=0;
            while(bMoreBack||bMoreFront)
            {
                /*ncnt++;
                if(ncnt>10000)
                {
                    std::cout<<"???????????????????????????????????????/n";
                }*/
  
                int nF_NearJ=-1;
                int nB_NearJ=-1;
                int nFConCode =0;
                int nBConCode =0;
                int nTCodeF=0;
                int nTCodeB=0;
                double dNearestDistF=10000;
                double dNearestDistB=10000;
                double dDistF =0,dDistB=0;
                
                
                
                
                
                
                for(int j=0; j<vp3DQues.size(); j++)
                {
                    
                    dDistF=1000; dDistB=1000;
                    if(vp3DQues[j]->m_bIsMerged)
                         continue;
                    
                    
                    if(bMoreFront)
                    {
                       nTCodeF =CMCheckConnObjIJ(pFront, vp3DQues[j], nWidth, nHeight, nLayers,CylinderPool,pPredPool,bFIsCheckFront,ppuch3DConLocal);
                    }
                    else
                    {
                        nTCodeF=0;
                    }
                    if(bMoreBack)
                    {
                       nTCodeB = CMCheckConnObjIJ(pBack, vp3DQues[j], nWidth, nHeight, nLayers,CylinderPool,pPredPool,bBIsCheckFront,ppuch3DConLocal);
                    }
                    else
                    {
                        nTCodeB=0;
                    }
                    
                    if(nTCodeF!=0)
                    {
                        
                        dDistF = DistP2P(pFront, vp3DQues[j], nTCodeF,pPredPool);
                        std::cout<<"nTCodeF="<<nTCodeF<<'\n'<<"dDistF="<<dDistF<<'\n';
                    }
                    
                    if(nTCodeB!=0)
                    {
                        dDistB = DistP2P(pBack, vp3DQues[j], nTCodeB,pPredPool);
                        std::cout<<"nTCodeB="<<nTCodeB<<'\n'<<"dDistB="<<dDistB<<'\n';
                    }
                    
                    if(nTCodeB!=0&&nTCodeF!=0)
                    {
                      if(dDistF<=dDistB)
                      {
                         if(dDistF<dNearestDistF)
                         {
                             dNearestDistF = dDistF;
                             nF_NearJ = j;
                             nFConCode = nTCodeF;
                          }
                       }
                      else
                      {
                         if(dDistB<dNearestDistB)
                         {
                            dNearestDistB = dDistB;
                            nB_NearJ = j;
                            nBConCode = nTCodeB;
                         }
                       }
                     }
                    else if(nTCodeB!=0)
                    {
                        if(dDistB<dNearestDistB)
                        {
                            dNearestDistB = dDistB;
                            nB_NearJ = j;
                            nBConCode = nTCodeB;
                        }
                    }
                    else if(nTCodeF!=0)
                    {
                        if(dDistF<dNearestDistF)
                        {
                            dNearestDistF = dDistF;
                            nF_NearJ = j;
                            nFConCode = nTCodeF;
                        }
                    }
                    
                    
                    
                }
                
                if(nF_NearJ==-1)
                {
                    bMoreFront=false;
                    pFront->m_pQueHead=NULL;
                }
                else
                {
                    std::cout<<"Front+1 /n";
                    bMoreFront=true;
                    vp3DQues[nF_NearJ]->m_bIsMerged=true;
                    vp3DQues[nF_NearJ]->m_pQueTail = pFront;
                    pFront->m_pQueHead = vp3DQues[nF_NearJ];

   
                    pFront = vp3DQues[nF_NearJ];
                    
                    if(nFConCode == 1||nFConCode == 3)
                    {
                        bFIsCheckFront = false;
                    }
                    else
                    {
                        bFIsCheckFront = true;
                    }
                }
                
                
                
                if(nB_NearJ==-1)
                {
                    bMoreBack=false;
                    pBack->m_pQueTail=NULL;
                }
                else
                {
                    std::cout<<"Back+1 /n";
                    
                    bMoreBack=true;
                    pBack->m_pQueTail = vp3DQues[nB_NearJ];
                    vp3DQues[nB_NearJ]->m_pQueHead = pBack;
                    vp3DQues[nB_NearJ]->m_bIsMerged=true;
                    pBack = vp3DQues[nB_NearJ];
                    
                    if(nBConCode == 1||nBConCode == 3)
                    {
                        bBIsCheckFront = false;
                    }
                    else
                    {
                        bBIsCheckFront = true;
                    }
                }
                
                
                
            }
            
            vp3DStartLong.push_back(pBack);
        }
        

    }
    
    
    
    
}

double DistP2P(Que3D* pQobject1, Que3D* pQobject2, int nConCode,std::map<CMark,PreCMark>* pPredPool)
{
    double dx1=pQobject1->m_dX,dy1=pQobject1->m_dY,dz1=pQobject1->m_dZ;
    double dx2=pQobject2->m_dX,dy2=pQobject2->m_dY,dz2=pQobject2->m_dZ;

    PreCMark pm1= (*pPredPool)[pQobject1->m_mark];
    PreCMark pm2= (*pPredPool)[pQobject2->m_mark];
    
    switch (nConCode) {
        case 1:
            dx1+= pm1.m_dFCenX; dy1+= pm1.m_dFCenY; dz1+= pm1.m_dFCenZ;
            dx2+= pm2.m_dFCenX; dy2+= pm2.m_dFCenY; dz2+= pm2.m_dFCenZ;
            break;
        case 2:
            dx1+= pm1.m_dFCenX; dy1+= pm1.m_dFCenY; dz1+= pm1.m_dFCenZ;
            dx2+= pm2.m_dBCenX; dy2+= pm2.m_dBCenY; dz2+= pm2.m_dBCenZ;
            break;
        case 3:
            dx1+= pm1.m_dBCenX; dy1+= pm1.m_dBCenY; dz1+= pm1.m_dBCenZ;
            dx2+= pm2.m_dFCenX; dy2+= pm2.m_dFCenY; dz2+= pm2.m_dFCenZ;
            break;
        case 4:
            dx1+= pm1.m_dBCenX; dy1+= pm1.m_dBCenY; dz1+= pm1.m_dBCenZ;
            dx2+= pm2.m_dBCenX; dy2+= pm2.m_dBCenY; dz2+= pm2.m_dBCenZ;
            break;
        default:
            break;
    }
    
    return Dist2Points3D(dx1,dy1,dz1,dx2,dy2,dz2);
    
}

void DrawMergeShortID_modified(std::vector<Que3D*>& vp3DStartLong, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, int SubVn)
{
    int fiber_counter = 0;
    
    std::string txt_fileName = "./SUBVOLUMES/sV" + std::to_string(SubVn) + "/fibers_info/merged_fiber.txt";
    std::ofstream outfile, outfile_txt;
    outfile_txt.open(txt_fileName.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    double R,H,Ty,Tz,x,y,z;

     for(int i=0; i<vp3DStartLong.size(); i++)
     {
         //if(vp3DStartLong[i]->m_pQueHead==NULL && vp3DStartLong[i]->m_mark.m_nH<10)
         //        continue;
         
         //Camilo Draw
         //sRGB temptID = sRGB(((i&0xf00)>>8)*10,((i&0xf0)>>4)*10,(i&0xf)*10);
         
         fiber_counter = DrawCylHDirectColor_Modified(&R,&H,&Ty,&Tz,&x,&y,&z,vp3DStartLong[i], ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool, fiber_counter);
         outfile_txt << fiber_counter << "," << R << "," << H << "," << Ty << "," << Tz << "," << x << "," << y << "," << z << "\n";
      }
      outfile_txt.close();
}

void DrawMergeShortID(std::vector<Que3D*>& vp3DStartLong, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    //int fiber_counter = 1;
     for(int i=0; i<vp3DStartLong.size(); i++)
     {
         //if(vp3DStartLong[i]->m_pQueHead==NULL && vp3DStartLong[i]->m_mark.m_nH<10)
         //        continue;
         
         //Camilo Draw
         //sRGB temptID = sRGB(((i&0xf00)>>8)*10,((i&0xf0)>>4)*10,(i&0xf)*10);
         
         //fiber_counter = DrawCylHDirectColor_Modified(vp3DStartLong[i], ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool, fiber_counter);         
      }
}

void DrawMergeShortResult(std::vector<Que3D*>& vp3DStartLong, unsigned char** ppuch3Dadd,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    sRGB rgbTable[20]={ sRGB(0,255,0),sRGB(0,0,255),sRGB(0,255,255),sRGB(255,255,0),sRGB(255,0,0),sRGB(255,255,255),sRGB(255,0,255),sRGB(128,180,30),sRGB(20,60,180),sRGB(180,40,110),sRGB(20,255,20),sRGB(0,50,255),sRGB(50,205,255),sRGB(205,225,20),sRGB(255,50,10),sRGB(160,160,160),sRGB(180,10,180),sRGB(100,160,50),sRGB(26,66,199),sRGB(200,70,80)
    };
    
    
    for(int i=0; i<vp3DStartLong.size(); i++)
    {
        
        
        DrawCylHDirectColor(vp3DStartLong[i], ppuch3Dadd, nWidth, nHeight, nLayers, CylinderPool, pPredPool, rgbTable[i%20]);

    }
    
    
}


bool CheckSpecialAngle(int nTY1,int nTZ1, int nTY2, int nTZ2)
{
    bool bSpecial = false;
    
    //Special case
    
    if((nTY1+nTY2<DEGREE_Y+10)&&(nTY1+nTY2>DEGREE_Y-10))
    {
        if(abs(nTZ1-nTZ2)>80)
            bSpecial=true;
        
        // std::cout<<"aaaaaaaaaaaaaa!/n";
    }
    
    
    
    if((nTY1<10||nTY1>80)&&(nTY2<10||nTY2>80))
    {
        bSpecial=true;
    }
    
    return bSpecial;
}

bool CheckSpecialAngleBigD(int nTY1,int nTZ1, int nTY2, int nTZ2)
{
    bool bSpecial = false;
    
    //Special case
    
    if((nTY1+nTY2<DEGREE_Y+15)&&(nTY1+nTY2>DEGREE_Y-15))
    {
        if(abs(nTZ1-nTZ2)>75)
            bSpecial=true;
        
       // std::cout<<"aaaaaaaaaaaaaa!/n";
    }
    
    
    
    if((nTY1<15||nTY1>75)&&(nTY2<15||nTY2>75))
    {
        bSpecial=true;
    }
    
    return bSpecial;
    
}

int CMCheckConnObjIJ(Que3D* pObjectI, Que3D* pObjectJ,int nWidth, int nHeight, int nLayers,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool,bool bCheckF,unsigned char** ppuch3DConLocal)
{
    
    //Check radius
    if(ABS(pObjectI->m_mark.m_nR-pObjectJ->m_mark.m_nR)>2)
        return 0;
    
    int nAngleDZ = abs(pObjectI->m_mark.m_nthetaZ-pObjectJ->m_mark.m_nthetaZ);
    
    int nAngleDY = abs(pObjectI->m_mark.m_nthetaY-pObjectJ->m_mark.m_nthetaY);

    

    //Check obj angle
    
    bool bSpecial = false;
    
    //Special case
    
    bSpecial= CheckSpecialAngleBigD(pObjectI->m_mark.m_nthetaY,pObjectI->m_mark.m_nthetaZ,pObjectJ->m_mark.m_nthetaY,pObjectJ->m_mark.m_nthetaZ);
    
    
    if((nAngleDY>15||nAngleDZ>15)&&bSpecial==false)
          return 0;
    
    
    //Check connect line angle
    
    double dx = pObjectJ->m_dX - pObjectI->m_dX;
    double dy = pObjectJ->m_dY - pObjectI->m_dY;
    double dz = pObjectJ->m_dZ - pObjectI->m_dZ;
    
    
    
    
    int nAngleCLY;
    int nAngleCLZ ;
    
    ConvertFromXYZ2DYZ_N(nAngleCLY, nAngleCLZ,dx,dy,dz);
    
   
    int nAngleDZ1= abs(pObjectI->m_mark.m_nthetaZ-nAngleCLZ);
    int nAngleDY1= abs(pObjectI->m_mark.m_nthetaY-nAngleCLY);
    

    
    int nAngleDY2= abs(pObjectJ->m_mark.m_nthetaY-nAngleCLY);
    
    int nAngleDZ2= abs(pObjectJ->m_mark.m_nthetaZ-nAngleCLZ);

    
    
    bool bSpecialI= CheckSpecialAngleBigD(pObjectI->m_mark.m_nthetaY,pObjectI->m_mark.m_nthetaZ,nAngleCLY,nAngleCLZ);
    
    bool bSpecialJ= CheckSpecialAngleBigD(pObjectJ->m_mark.m_nthetaY,pObjectJ->m_mark.m_nthetaZ,nAngleCLY,nAngleCLZ);
    
    bSpecial = (bSpecialI||bSpecialJ);
    
    
    if((MIN(nAngleDY1,nAngleDY2)>12||MIN(nAngleDZ1,nAngleDZ2)>12)&&bSpecial==false)
        return 0;
    
    
    Clean3DGrayImages(nLayers, nWidth, nHeight, ppuch3DConLocal);
    
    if(bCheckF)
    {
        AddObjectLocalCONN(pObjectI,ppuch3DConLocal,nWidth, nHeight,nLayers,CylinderPool,pPredPool,3);
    }
    else
    {
        AddObjectLocalCONN(pObjectI,ppuch3DConLocal,nWidth, nHeight,nLayers,CylinderPool,pPredPool,4);
    }
    
    
    
    
     std::vector<TPoint> vTPoints = CylinderPool[pObjectJ->m_mark];
     
     int nCenX = round(pObjectJ->m_dX);
     int nCenY = round(pObjectJ->m_dY);
     int nCenZ = round(pObjectJ->m_dZ);
     
     
     int nFFCount =0, nFBCount=0, nBFCount=0, nBBCount=0;
     int nCountF =0; int nCountB=0;
     for(int i=0; i<vTPoints.size(); i++)
     {
     
     int nx = vTPoints[i].m_nX+nCenX;
     int ny = vTPoints[i].m_nY+nCenY;
     int nz = vTPoints[i].m_nZ+nCenZ;
     
     if(nx<0||nx>=nWidth||
     ny<0||ny>=nHeight||
     nz<0||nz>=nLayers)
     {
     continue;
     }
     
     if(vTPoints[i].m_nType==3)
     {
     
     if(ppuch3DConLocal[nz][nx+ny*nWidth]==3)
     {
     nFFCount++;
     }
     if(ppuch3DConLocal[nz][nx+ny*nWidth]==4)
     {
     nBFCount++;
     }
     
     }
     else if(vTPoints[i].m_nType==4)
     {
     
     if(ppuch3DConLocal[nz][nx+ny*nWidth]==3)
     {
     nFBCount++;
     }
     if(ppuch3DConLocal[nz][nx+ny*nWidth]==4)
     {
     nBBCount++;
     }
     }
     
     }
     
     if(nFFCount==0&&nFBCount==0&&nBFCount==0&&nBBCount==0)
     {
     return 0;
     }
     
     
     int nMaxCnt = MAX(MAX(nFFCount,nFBCount),MAX(nBFCount,nBBCount));
     
     if(nFFCount==nMaxCnt)
     return 1;
     if(nFBCount==nMaxCnt)
     return 2;
     if(nBFCount==nMaxCnt)
     return 3;
     if(nBBCount==nMaxCnt)
     return 4;
     
     return 0;
    
}


void Convert2Long(std::vector<Que3D*>& vp3DStartLong, std::vector<Que3D*>& vp3DLongF,unsigned char** ppuch3DGray, unsigned char** ppuchSeg3D,int** ppn3DConnect,int** ppnOverlap,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP)
{
    
    vp3DLongF.clear();
    cv::Mat matImgA[MAXSequence];
    unsigned char* ppuch3DA[MAXSequence];
    CreateVoid3DGrayImages(nLayers, nWidth,nHeight, matImgA, ppuch3DA);
    
    cv::Mat matImgB[MAXSequence];
    unsigned char* ppuch3DB[MAXSequence];
    CreateVoid3DGrayImages(nLayers, nWidth,nHeight, matImgB, ppuch3DB);
    
    for(int i=0; i<vp3DStartLong.size(); i++)
    {
        Que3D* qLongNew = new Que3D;
        TDPoint tpBegin,tpEnd;
        vp3DStartLong[i]->m_Pmark = (*pPredPool)[vp3DStartLong[i]->m_mark];
        int nIsSingle;
        *qLongNew = ConvertOne2Long(vp3DStartLong[i],tpBegin, tpEnd,nIsSingle);
        
        if(nIsSingle==1)
        {
            vp3DLongF.push_back(qLongNew);
            continue;
        }
        
        DrawMergeFiber(vp3DStartLong[i], ppuch3DB, nWidth,nHeight, nLayers,CylinderPool, pPredPool);
        
        TDPoint tpBegin2,tpEnd2;
        
        Que3D qTempt = *qLongNew;

        int nRange =1;
        int nR = qLongNew->m_mark.m_nR;
        double dBestScore=-100;
        double dTScore=0;
        TDPoint tpBestBegin=tpBegin,tpBestEnd=tpEnd;
        
      //  std::cout<<"the length is "<< Dist2Points3D(tpBegin.m_dX, tpBegin.m_dY, tpBegin.m_dZ, tpEnd.m_dX, tpEnd.m_dY, tpEnd.m_dZ)<<'\n';      //  std::cout<<"Drawing....\n";
        int nCntIt=0;
        for(int nx =-nRange; nx<nRange; nx++)
            for(int ny = -nRange; ny<nRange; ny++)
                for(int nz = -nRange; nz<nRange; nz++)
                    for(int nx1 =-nRange; nx1<nRange; nx1++)
                        for(int ny1 = -nRange; ny1<nRange; ny1++)
                            for(int nz1 = -nRange; nz1<nRange; nz1++)
                            {
                       //         std::cout<<"Current "<<nCntIt++<<" th iter\n";
                                
                                
                                tpBegin2.m_dX=tpBegin.m_dX+nx;
                                tpBegin2.m_dY=tpBegin.m_dY+ny;
                                tpBegin2.m_dZ=tpBegin.m_dZ+nz;
                                
                                tpEnd2.m_dX = tpEnd.m_dX+nx1;
                                tpEnd2.m_dY = tpEnd.m_dY+ny1;
                                tpEnd2.m_dZ = tpEnd.m_dZ+nz1;
                                
                                DrawLongFiberB2PT(ppuch3DA, nWidth, nHeight,nLayers, tpBegin2, tpEnd2, nR,CylinderPool,pPredPool);
                                dTScore = CalMatchScore(ppuch3DA,ppuch3DB,nWidth,nHeight,  nLayers);
                                if(dTScore>dBestScore)
                                {
                                    dBestScore=dTScore;
                                    tpBestBegin = tpBegin2;
                                    tpBestEnd = tpEnd2;
                                }

                            
                            }
        
        
        qLongNew->m_dX = (tpBestBegin.m_dX+tpBestEnd.m_dX)/2;
        qLongNew->m_dY = (tpBestBegin.m_dY+tpBestEnd.m_dY)/2;
        qLongNew->m_dZ = (tpBestBegin.m_dZ+tpBestEnd.m_dZ)/2;
        qLongNew->m_mark.m_nH = Dist2Points3D(tpBestBegin.m_dX, tpBestBegin.m_dY, tpBestBegin.m_dZ, tpBestEnd.m_dX, tpBestEnd.m_dY, tpBestEnd.m_dZ)/2.0;
        
        
     //   std::cout<<"After search length = " <<Dist2Points3D(tpBestBegin.m_dX, tpBestBegin.m_dY, tpBestBegin.m_dZ, tpBestEnd.m_dX, tpBestEnd.m_dY, tpBestEnd.m_dZ)<<'\n';
        int nThetaY, nThetaZ;
        
        ConvertFromXYZ2DYZ_N(nThetaY, nThetaZ,tpBestBegin.m_dX-tpBestEnd.m_dX,tpBestBegin.m_dY-tpBestEnd.m_dY, tpBestBegin.m_dZ-tpBestEnd.m_dZ);
        
        qLongNew->m_mark.m_nthetaY = nThetaY;
        qLongNew->m_mark.m_nthetaZ = nThetaZ;
        qLongNew->m_mark.m_nR =nR;
        
        
    //    std::cout<<"Finsh search!\n";
        if(CylinderPool.find(qLongNew->m_mark)==CylinderPool.end())
        {
            
   //         OutputObjInfo(qLongNew, "Long Fiber");
            CreateOneCylinder(qLongNew->m_mark, CylinderPool, pPredPool);
        }
        
        vp3DLongF.push_back(qLongNew);
        
    }
    
    
    
}



Que3D ConvertOne2Long(Que3D* pStartQue,TDPoint& TDBegin, TDPoint& TDEnd,int& nIsSingle)
{
    Que3D qLongCylinder;
    
    Que3D* pHeader = pStartQue;
    nIsSingle =0;
    
    if(pHeader->m_pQueHead==NULL)
    {
        nIsSingle=1;
        qLongCylinder = *pStartQue;
        return qLongCylinder;
    }
    
    
    
    double dTx0,dTy0,dTz0,dHx0,dHy0,dHz0;
    if(pHeader->m_dFConnR<=pHeader->m_dBConnR)
    {
        dTx0 = pHeader->m_dX+ pHeader->m_Pmark.m_dFCenX;
        dTy0 = pHeader->m_dY+ pHeader->m_Pmark.m_dFCenY;
        dTz0 = pHeader->m_dZ+ pHeader->m_Pmark.m_dFCenZ;
    }
    else
    {
        dTx0 = pHeader->m_dX+ pHeader->m_Pmark.m_dBCenX;
        dTy0 = pHeader->m_dY+ pHeader->m_Pmark.m_dBCenY;
        dTz0 = pHeader->m_dZ+ pHeader->m_Pmark.m_dBCenZ;
    }
    
    
    
    
    int nCntR =0;
    int nCntH =0;
    double dCx=0;
    double dCy=0;
    double dCz=0 ;
    double dTx = pHeader->m_dX;
    double dTy = pHeader->m_dY;
    double dTz = pHeader->m_dZ;
    
    int nTR,nHR;
    nTR = pHeader->m_mark.m_nR;
    
    do
    {
        
        
        dCx += pHeader->m_dX*pHeader->m_mark.m_nH;
        dCy += pHeader->m_dY*pHeader->m_mark.m_nH;
        dCz += pHeader->m_dZ*pHeader->m_mark.m_nH;
        nCntR += pHeader->m_mark.m_nR*pHeader->m_mark.m_nH;
        nCntH += pHeader->m_mark.m_nH;
        pHeader = pHeader->m_pQueHead;
        
    }
    while(pHeader->m_pQueHead!=NULL);
    
    

    if(pHeader->m_dFConnR<=pHeader->m_dBConnR)
    {
        dHx0 = pHeader->m_dX+ pHeader->m_Pmark.m_dFCenX;
        dHy0 = pHeader->m_dY+ pHeader->m_Pmark.m_dFCenY;
        dHz0 = pHeader->m_dZ+ pHeader->m_Pmark.m_dFCenZ;
    }
    else
    {
        dHx0 = pHeader->m_dX+ pHeader->m_Pmark.m_dBCenX;
        dHy0 = pHeader->m_dY+ pHeader->m_Pmark.m_dBCenY;
        dHz0 = pHeader->m_dZ+ pHeader->m_Pmark.m_dBCenZ;
    }
    
    
    nHR = pHeader->m_mark.m_nR;
    
    double dMx = dCx/nCntH;
    double dMy = dCy/nCntH;
    double dMz = dCz/nCntH;
    double dMx0 = (dHx0+dTx0)/2.0;
    double dMy0 = (dHy0+dTy0)/2.0;
    double dMz0 = (dHz0+dTz0)/2.0;
    
    dMx = dMx;//(dMx0+dMx)/2.0;
    dMy = dMy;//(dMy0+dMy)/2.0;
    dMz = dMz;//(dMz0+dMz)/2.0;
    
    
    double dHx = pHeader->m_dX;
    double dHy = pHeader->m_dY;
    double dHz = pHeader->m_dZ;

    int nR = int(nCntR/double(nCntH)+0.5);
    int nH = (nCntH+int (Dist2Points3D(dHx, dHy, dHz, dTx, dTy, dTz)+nHR+nTR))/2;
    
    
    int nThetaY, nThetaZ;
    
    
    ConvertFromXYZ2DYZ_N(nThetaY, nThetaZ,dHx-dTx,dHy-dTy, dHz-dTz);
    
    CMark  mark(nR, nH, nThetaY, nThetaZ);
    qLongCylinder.m_dX = dMx;
    qLongCylinder.m_dY = dMy;
    qLongCylinder.m_dZ = dMz;
    qLongCylinder.m_mark = mark;
    
    TDBegin.m_dX = dTx0;
    TDBegin.m_dY = dTy0;
    TDBegin.m_dZ = dTz0;
    
    TDEnd.m_dX = dHx0;
    TDEnd.m_dY = dHy0;
    TDEnd.m_dZ = dHz0;
    
    return qLongCylinder;

}


void DrawLongFiberB2PT(unsigned char** ppuch3Dimgs, int nWidth,int nHeight, int nLayer, TDPoint TDBegin, TDPoint TDEnd,int nR,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    Clean3DGrayImages(nLayer, nWidth, nHeight, ppuch3Dimgs);
    
    int nThetaY, nThetaZ;
    
    
    ConvertFromXYZ2DYZ_N(nThetaY, nThetaZ,TDBegin.m_dX-TDEnd.m_dX,TDBegin.m_dY-TDEnd.m_dY, TDBegin.m_dZ-TDEnd.m_dZ);
    
    double dlen = Dist2Points3D(TDBegin.m_dX, TDBegin.m_dY, TDBegin.m_dZ, TDEnd.m_dX, TDEnd.m_dY, TDEnd.m_dZ);
    
    TDPoint TDtepmpt = TDBegin;
    
    while(dlen>2*nR)
    {
        int nk=NMAXH*2;
        while(nk>dlen)
        {
            nk-=2;
        }
        if(nk<NMINH*2)
        {
            break;
        }
        
        CMark mark(nR,nk/2,nThetaY,nThetaZ);
        PreCMark PMark = (*pPredPool)[mark];
        Que3D qobject;
        qobject.m_mark = mark;
        
        double dFx = PMark.m_dFCenX+TDtepmpt.m_dX;
        double dFy = PMark.m_dFCenY+TDtepmpt.m_dY;
        double dFz = PMark.m_dFCenZ+TDtepmpt.m_dZ;
        
        double dBx = PMark.m_dBCenX+TDtepmpt.m_dX;
        double dBy = PMark.m_dBCenY+TDtepmpt.m_dY;
        double dBz = PMark.m_dBCenZ+TDtepmpt.m_dZ;
        bool bDirectF = true;
        if(SqDist2Points3D(dFx,dFy,dFz,TDEnd.m_dX,TDEnd.m_dY,TDEnd.m_dZ)<=SqDist2Points3D(dBx,dBy,dBz,TDEnd.m_dX,TDEnd.m_dY,TDEnd.m_dZ))
        {
            bDirectF =true;
            qobject.m_dX = dFx;
            qobject.m_dY = dFy;
            qobject.m_dZ = dFz;
        }
        else
        {
            bDirectF=false;
            qobject.m_dX = dBx;
            qobject.m_dY = dBy;
            qobject.m_dZ = dBz;
            
        }
        
        DrawCylinder(qobject.m_dX+0.5, qobject.m_dY+0.5, qobject.m_dZ+0.5,mark, ppuch3Dimgs, nWidth, nHeight, nLayer,  CylinderPool,pPredPool);
        
        dlen-=nk*2;
        
        if(bDirectF)
        {
            TDtepmpt.m_dX = dFx+PMark.m_dFCenX;
            TDtepmpt.m_dY = dFy+PMark.m_dFCenY;
            TDtepmpt.m_dZ = dFz+PMark.m_dFCenZ;
        }
        else
        {
            TDtepmpt.m_dX = dBx+PMark.m_dBCenX;
            TDtepmpt.m_dY = dBy+PMark.m_dBCenY;
            TDtepmpt.m_dZ = dBz+PMark.m_dBCenZ;
        }
        
    }
    
}

void DrawMergeFiber(Que3D* pStartQue, unsigned char** ppuch3Dimgs, int nWidth,int nHeight, int nLayer,std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
    Clean3DGrayImages(nLayer, nWidth, nHeight, ppuch3Dimgs);
    Que3D* pHead = pStartQue;
    
    do
    {
        DrawCylinder(pHead->m_dX+0.5, pHead->m_dY+0.5, pHead->m_dZ+0.5,pHead->m_mark, ppuch3Dimgs, nWidth, nHeight, nLayer,  CylinderPool,pPredPool);
        
        pHead = pHead->m_pQueHead;
    }
    while(pHead->m_pQueHead!=NULL);
}


double CalMatchScore(unsigned char** ppuch3DimgsA,unsigned char** ppuch3DimgsB,int nWidth,int nHeight, int nLayer)
{
    double dScore=0;
    int nCntA=0, nCntB=0, nCntSame=0;
    double dTA,dTB;
    
    for(int nz=0; nz<nLayer; nz++)
        for(int ny=0; ny<nHeight; ny++)
            for(int nx=0; nx<nWidth; nx++)
            {
                if(Get3DPixel(nx,ny,nz,nWidth,ppuch3DimgsA)!=0)
                {
                    nCntA++;
                }
                
                if(Get3DPixel(nx,ny,nz,nWidth,ppuch3DimgsB)!=0)
                {
                    nCntB++;
                }
                if(Get3DPixel(nx,ny,nz,nWidth,ppuch3DimgsA)!=0&&Get3DPixel(nx,ny,nz,nWidth,ppuch3DimgsB)!=0)
                {
                    nCntSame++;
                }
                
            }
    
    dTA = double(nCntSame)/nCntA;
    dTB = double(nCntSame)/nCntB;
    
    dScore = dTB*0.65+dTA*0.35;
    
    return dScore;
}



void SplitLongFiber(std::vector<Que3D*>& vp3DStartLong,std::vector<Que3D*>& vp3DFinal, unsigned char** ppuch3DGray,int nWidth, int nHeight, int nLayers, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool, MppControPara& conP)
{
    
}


void CalculateNeigborID(std::vector<Que3D*>& vObjects, int nSqNbDist)
{
    int nTotal = vObjects.size();
    int nTmptDistIJ=0;
    for(int i=0; i<nTotal; i++)
    {
        vObjects[i]->m_vnNeighborIDs.clear();
    }
    
    for(int i=0; i<nTotal-1; i++)
        for(int j=i+1; j<nTotal; j++)
        {
            
             nTmptDistIJ = SqDist2Points3D(vObjects[i]->m_dX,vObjects[i]->m_dY,vObjects[i]->m_dZ, vObjects[j]->m_dX,vObjects[j]->m_dY,vObjects[j]->m_dZ);
            
             if(nTmptDistIJ<=nSqNbDist)
             {
                 vObjects[i]->m_vnNeighborIDs.push_back(j);
                 vObjects[j]->m_vnNeighborIDs.push_back(i);
             }
    
        }
    
}

void WriteFiberResults(std::vector<Que3D*>& vp3DFinal, int SubVn)
{
    int nLen = vp3DFinal.size();
    
    std::vector<Que3D> v3DObjects;
    
    std::cout<<"TotalWriteFiber is "<<nLen<<'\n';
    
    for(int i=0; i<nLen; i++)
    {
        Que3D tempt = *(vp3DFinal[i]);
        
        v3DObjects.push_back(tempt);
    }
    
    
    
    std::string fileName = "./SUBVOLUMES/sV" + std::to_string(SubVn) + "/fibers_info/fiber.data";
    std::ofstream outfile, outfile_txt;
    outfile.open(fileName.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    std::string txt_fileName = "./SUBVOLUMES/sV" + std::to_string(SubVn) + "/fibers_info/fiber.txt";
    outfile_txt.open(txt_fileName.c_str(), std::ios::out|std::ios::binary|std::ios::trunc);
    
    int nCnt=0;

    outfile_txt << nLen << "\n";
    outfile.write(((char*)&nLen), sizeof(int));
    
    int fnum, R, H, thetaY,thetaZ;
    double x,y,z;
    //double m_dX;
    //double m_dY;
    //double m_dZ;

    for(int i=0; i<nLen; i++)
    {
        outfile.write(((char*)&v3DObjects[i]), sizeof(Que3D));
        
        // Camilo save stuff to a text file
        fnum = i;
        R = v3DObjects[i].m_mark.m_nR;
        H = v3DObjects[i].m_mark.m_nR;
        thetaY = v3DObjects[i].m_mark.m_nthetaY;
        thetaZ = v3DObjects[i].m_mark.m_nthetaZ;
        x = v3DObjects[i].m_dX;
        y = v3DObjects[i].m_dY;
        z = v3DObjects[i].m_dZ;

        outfile_txt << fnum << "," << R << "," << H << "," << thetaY << "," << thetaZ << "," << x << "," << y << "," << z << "\n";

    }
    
    v3DObjects.clear();
    
    outfile.close();
    outfile_txt.close();
    
    
    
}

void ReadFiberResults(std::vector<Que3D*>& vp3DFinal, std::string fileName)
{
    std::ifstream infile;
    infile.open(fileName.c_str(), std::ios::in|std::ios::binary);
    int nCnt2=0;
    // int Test[10000];
    
    int nTotalFiber = 0;
    
    infile.read(((char*)& nTotalFiber), sizeof(int));
    vp3DFinal.clear();
    
    
    std::cout<<"read toatal fiber is "<<nTotalFiber<<'\n';
    
    for(int i=0; i<nTotalFiber; i++)
    {
        
        
        Que3D* pQue3d = new Que3D;
        infile.read((char*)pQue3d, sizeof(Que3D));
        
        
        
        
        vp3DFinal.push_back(pQue3d);
        
        
        //OutputObjInfo(pQue3d, "^_^");
        
    }

    
    
    infile.close();
    
    
}


void ShowSavedFibers(int nWidth, int nHeight, int nLayer, std::string fileName, std::map<CMark,std::vector<TPoint> >& CylinderPool,std::map<CMark,PreCMark>* pPredPool)
{
    
    cv::Mat matColorSeq[MAXSequence];
    cv::Mat matGraySeq[MAXSequence];
    unsigned char* ppuch3DRGB[MAXSequence];
    unsigned char* ppuch3DGray[MAXSequence];
    
    
    CreateVoid3DColorImages(nLayer, nWidth, nHeight, matColorSeq, ppuch3DRGB);
    CreateVoid3DGrayImages(nLayer, nWidth, nHeight, matGraySeq, ppuch3DGray);
    //Load3DGrayImages("./ReadFiberSeg/RSeg_", nLayer, matGraySeq, ppuch3DGray);
    std::vector<Que3D*> vp3DFinal;
    
    
    std::cout<<"beginning reading fiber results!\n";
    ReadFiberResults(vp3DFinal, fileName);
    
    std::cout<<"finished read fiber results!\n";
    
    
    
    
    
    
    //return;
    
    
    
    for(int i=0; i<nLayer;i ++)
        matColorSeq[i].setTo(0);
    
    DrawResult3D(vp3DFinal, ppuch3DRGB, nWidth, nHeight, nLayer, CylinderPool, pPredPool);
    
    Save3DImages("./SUBVOLUMES/temp/r_",nLayer, matColorSeq);
    
    for(int i=0; i<nLayer;i ++)
        matColorSeq[i].setTo(0);
    
    
    
    
    DrawResult3D(vp3DFinal, ppuch3DRGB, nWidth, nHeight, nLayer, CylinderPool, pPredPool,true);
    Save3DImages("./SUBVOLUMES/temp2/RSeg_",nLayer, matColorSeq);
    
    
    
    for(int i=0; i<nLayer; i++)
    {
        cv::cvtColor(matColorSeq[i], matGraySeq[i], CV_BGR2GRAY);
    }
        
    
}

