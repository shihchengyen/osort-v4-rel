function setPath

cd ./DPV; 
nptAddPath; 
cd ..; 
cd ./newNpt; 
nptAddPath; 
cd ..; 
cd ./Hippocampus; 
addpath(pwd); 
cd Compiler; 
addpath(genpath(pwd)); 
cd ../..; 
cd neuroshare; 
addpath(pwd); 
cd ..; 
cd hmmsort.py; 
addpath(pwd); 
cd ..; 
cd hmmsort.py/helper_functions; 
addpath(pwd); 
cd ../..; 
cd osort-v4-rel; 
nptAddPath; 
cd ..;