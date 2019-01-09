/*************************************************************************/
/*                                                                       */
/*                      Copyright © 2018                                 */
/* The Board of Trustees of Purdue University,  All Rights Reserved      */
/*                       Author: Tianyu Li                               */
/*                                                                       */
/*************************************************************************/
//  model3D_types.cpp
//  HybridMPP_MFR
//
//  Created by Tianyu Li on 5/17/18.

#include "model3D_types.h"

bool CMark::operator> (const CMark& markB) const
{
    if(m_nR>markB.m_nR)
    {
        return true;
    }
    else if(m_nR==markB.m_nR)
    {
        if(m_nH>markB.m_nH)
            return true;
        else if(m_nH==markB.m_nH)
        {
            if(m_nthetaY>markB.m_nthetaY)
            {
                return true;
            }
            else if(m_nthetaY==markB.m_nthetaY)
            {
                if(m_nthetaZ>markB.m_nthetaZ)
                {
                    return true;
                }
                else return false;
            }
            else
                return false;
        }
        else
            return false;
    }
    else return false;
}
bool CMark::operator< (const CMark& markB) const
{
    
        if(m_nR<markB.m_nR)
        {
            return true;
        }
        else if(m_nR==markB.m_nR)
        {
            if(m_nH<markB.m_nH)
                return true;
            else if(m_nH==markB.m_nH)
            {
                if(m_nthetaY<markB.m_nthetaY)
                {
                    return true;
                }
                else if(m_nthetaY==markB.m_nthetaY)
                {
                    if(m_nthetaZ<markB.m_nthetaZ)
                    {
                        return true;
                    }
                    else return false;
                }
                else
                    return false;
            }
            else
                return false;
        }
        else return false;
    }

void OutputObjInfo(Que3D* qObj, std::string strPreInfo)
{
    std::cout<<"The Object "<<strPreInfo<<" is:\n";
    std::cout<<"dx = "<<qObj->m_dX<<" dy= "<<qObj->m_dY<<"dz= "<<qObj->m_dZ<<'\n';
    std::cout<<"nR = "<<qObj->m_mark.m_nR<<" nH = "<<qObj->m_mark.m_nH<<" nTY= "<<qObj->m_mark.m_nthetaY<<" nTZ= "<<qObj->m_mark.m_nthetaZ<<'\n';
    std::cout<<"dEn = "<<qObj->m_dDateEng<<" dOlp= "<<qObj->m_dOverlapR<<'\n';
    std::cout<<"Infor listed above\n";
    
}


int CalAngleFromYZ(int nTy, int nTz)
{
    if(nTy<3||nTy>=88)
    {
        return 0;
    }
    
    int nCodeY,nCodeZ;
    
    
    if(nTz>87)
    {
        nTz=0;
        nTy=89-nTy;
    }
    
    
    if(nTy<15)
    {
        nCodeY=0;
    }
    else if(nTy<25)
    {
        nCodeY=1;
    }
    else if(nTy<35)
    {
        nCodeY=2;
    }
    else if(nTy<45)
    {
        nCodeY=3;
    }
    else if(nTy<55)
    {
        nCodeY=4;
    }
    else if(nTy<65)
    {
        nCodeY=5;
    }
    else if(nTy<75)
    {
        nCodeY=6;
    }
    else
    {
        nCodeY=7;
    }
    

    
    
    if(nTz<8)
    {
        nCodeZ=0;
    }
    else if(nTz<18)
    {
        nCodeZ=1;
    }
    else if(nTz<28)
    {
        nCodeZ=2;
    }
    else if(nTz<38)
    {
        nCodeZ=3;
    }
    else if(nTz<48)
    {
        nCodeZ=4;
    }
    else if(nTz<58)
    {
        nCodeZ=5;
    }
    else if(nTz<68)
    {
        nCodeZ=6;
    }
    else if(nTz<78)
    {
        nCodeZ=7;
    }
    else
    {
        nCodeZ=8;
    }
    
    int nCode = nCodeY*9+nCodeZ+1;
    return nCode;
    
}

void NormalizeHist(int* pnHist, double* pdHist, int nChan)
{
    int nTotal=0;
    
    for(int i=0; i<nChan; i++)
    {
        nTotal+= pnHist[i];
        
    }
    
    if(nTotal ==0)
    {
        return;
    }
    
    for(int i=0;i<nChan;i++)
    {
        pdHist[i]= double(pnHist[i])/double(nTotal);
    }
    
    
    
}

void SaveObjStatInfo(FiberInfo& FI, std::string strInfor)
{
    std::ofstream outfile;
    
    outfile.open(strInfor.c_str(),std::ios::trunc|std::ios::out);
    
    outfile<<"Total fiber is "<<FI.m_nTotalFiber<<'\n';
    outfile<<"FB ratio is "<<FI.m_dRM2F<<'\n';
    outfile<<"Max Fiber Length is "<< FI.m_dMaxLen<<'\n';
    outfile<<"Min Fiber Length is "<< FI.m_dMinLen<<'\n';
    outfile<<"Mean Fiber Length is "<< FI.m_dMeanLen<<'\n';
    
    outfile<<"Max Fiber Diameter is "<< FI.m_dMaxDia<<'\n';
    outfile<<"Min Fiber Diameter is "<< FI.m_dMinDia<<'\n';
    outfile<<"Mean Fiber Diameter is "<< FI.m_dMeanDia<<'\n';
    
    
     outfile<<"Length Histogram is "<<'\n';
    for(int i=0; i<500; i++)
    {
        outfile<<FI.m_nHistLen[i]<<", ";
    }
    outfile<<'\n';
    
    for(int i=0; i<500; i++)
    {
        outfile<<FI.m_dHistLen[i]<<", ";
    }
    outfile<<'\n';
    
    
    
    outfile<<"Diameter Histogram is "<<'\n';
    for(int i=0; i<18; i++)
    {
        outfile<<FI.m_nHistDia[i]<<", ";
    }
    outfile<<'\n';
    
    for(int i=0; i<18; i++)
    {
        outfile<<FI.m_dHistDia[i]<<", ";
    }
    outfile<<'\n';
    
    
    outfile<<"Angle Histogram is "<<'\n';
    for(int i=0; i<73; i++)
    {
        outfile<<FI.m_nHistAngle[i]<<", ";
    }
    outfile<<'\n';
    
    for(int i=0; i<73; i++)
    {
        outfile<<FI.m_dHistAngle[i]<<", ";
    }
    outfile<<'\n';
    
    
    outfile.close();
    
    
    
}
