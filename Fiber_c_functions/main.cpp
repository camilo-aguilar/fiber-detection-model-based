/*************************************************************************/
/*                                                                       */
/*                      Copyright 2018-2020                              */
/* The Board of Trustees of Purdue University,  All Rights Reserved      */
/*                       Author: Tianyu Li                               */
/*                                                                       */
/*************************************************************************/

//  main.cpp
//  HybridMPP_MFR
//
//  Created by Tianyu Li on 3/17/17.
//  Copyright © 2017 StudyMan. All rights reserved.
//

#include <iostream>
#include "opencv2/opencv.hpp"
#include <map>
#include "CylinderModel3D.h"



int main(int argc, const char * argv[]) {

	int SubVn = 7;
	if(argc > 1)
		{
			SubVn = atoi(argv[1]);
			if(!SubVn)
			{
				std::cout << "Enter a valid number\n";
				exit(-1);
			}


		}
    
  TestFiberDetection(SubVn);  

    return 0;
}
