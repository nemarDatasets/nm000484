function par = load_par(paraIdx)
  %% parameters for EMG detector  
  %% we used Case 2


  featureNames = {'lenStd','lenMu','lenValue'}; %% alternative methods to compute relative line-length (for detemining EMG artifect) in comparison with baseline
  acpMethods = {'std','mu'}; %% alternative methods compute relative line-length (for detemining physiological variation) in comparison with baseline
switch paraIdx 
  case 1
    par.fs = 600;
    par.featureName = featureNames{2};
    par.acpMethod = acpMethods{2};
    par.thresLR = 6;
    par.thresDet = 2.5;
    par.thresAcp = 2;
    par.agreeNum = 1;
    par.wdSec = 0.2;
    par.shSec = 0.1;
    par.blSec = 5;
    par.resetBaseline = false;

  case 2
    par.fs = 500;
    par.featureName = featureNames{2};
    par.acpMethod = acpMethods{2};
    par.thresLR = 6;
    par.thresDet = 2.5;
    par.thresAcp = 2;
    par.agreeNum = 2;
    par.wdSec = 0.2;
    par.shSec = 0.1;
    par.blSec = 5;
    par.resetBaseline = true;

  case 3
    par.fs = 600;
    par.featureName = featureNames{2};
    par.acpMethod = acpMethods{2};
    par.thresLR = 6;
    par.thresDet = 2.5;
    par.thresAcp = 2;
    par.agreeNum = 1;
    par.wdSec = 0.2;
    par.shSec = 0.1;
    par.blSec = 5;
    par.resetBaseline = true;

  case 4
    par.fs = 500;
    par.featureName = featureNames{2};
    par.acpMethod = acpMethods{2};
    par.thresLR = 6;
    par.thresDet = 2.5;
    par.thresAcp = 4;
    par.agreeNum = 1;
    par.wdSec = 0.2;
    par.shSec = 0.1;
    par.blSec = 5;
    par.resetBaseline = false;
end
