function blnParallel = autoParallel

blnParallel = license('test','Distrib_Computing_Toolbox');

if blnParallel
    nCores = feature('numCores');
    if isempty(gcp('nocreate'))
        parpool('local',nCores-1);  % run on all, except one to stop machine from freezing
    end
end