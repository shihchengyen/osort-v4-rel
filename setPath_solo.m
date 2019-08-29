function setPath_solo()

    cd osort-v4-rel;
    addpath(genpath(pwd));
    cd ..;
    
cd ./DPV; 
nptAddPath; 
cd ..; 
cd ./Hippocampus; 
addpath(pwd);
cd ..;    

end