# Short Fiber Reinforced Composite Algorithms
Scripts to extract fibers and voids from short fiber reinforced polymers. 

REQUIREMENTS:  
  
	-OpenCV (Tested with Version 3.4.3)  
	-cmake (Tested with version 2.8.12 and 3.12.4)  
	-Access to Rice:standby cluster  
	-Matlab (Tested with R2016a)  
    
COMPILE C FILES (this will generate "FiberMpp" executable):  
  
	-chmod +x compile  
	-./compile  

DOWNLOAD CylinderPool.data:  
  
  	https://purdue0-my.sharepoint.com/:f:/g/personal/aguilarh_purdue_edu/EuNfn57UU5ZLoWVXmMGtmS8Bo-ntwdkaFqEtVYs7uZhGRg?e=aDjur9     
  
  
SET UP:  
  
	The root directory must contain:  
		-Folder "INPUT_FILES" with 1395 uint16 2560x2560 tiff images  
		-FiberMpp C executable 
		-CylinderPool.data (Precomputed clynders for C function)  
		-Folder "include" with 12 .m functions  
		-8 MATLAB scripts  

TO RUN:  
  
	-run script1.m (7 hrs)  
	-When all the subvolumes have been processed:  
	-run script2.m (4 hrs)  

