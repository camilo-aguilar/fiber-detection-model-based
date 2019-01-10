%% File: split_subvolumes
%  Author: Camilo Aguilar
%  Function: Send 75 scripts to cluster to find fibers using C program.
%
%  Input:  
%          FiberMpp C program
%          fiber.dat binary file
%          Folder: SUBVOLUMES/sV# (where # = 1,2,3,...75)
%          Each sV# Folder must have:
%                  Data (for input files)
%                  Seg (presegmentation)
%                  fibers (for output files)
%
%  Output: 
%           Fibers detected at: SUBVOLUMES/sV#/fibers
%           Fiber information stored at: SUBVOLUMES/sV#/fibers_info/
%


if(~exist('FIBER_SCRIPTS', 'dir' ))
    mkdir('FIBER_SCRIPTS')
end

cd FIBER_SCRIPTS

for i=1:75
    fileID = fopen(['myqsub' num2str(i) '.sh'],'w');
    command = [...
    '#!/bin/sh -l \n' ...
    '# FILENAME:  myqsub' num2str(i) '.sh\n'...
    '#PBS -q standby\n'...
    '#PBS -l nodes=1:ppn=1 \n'...
    '#PBS -l walltime=4:00:00\n\n'...
    '# Change to the directory \n' ...
    'cd $PBS_O_WORKDIR\n'...
    'cd .. \n'...
    './FiberMpp ' num2str(i) '\n'...
    ];

    fprintf(fileID,command);
    fclose(fileID);
    
    
    system(['qsub myqsub' num2str(i) '.sh']);
end

cd ..