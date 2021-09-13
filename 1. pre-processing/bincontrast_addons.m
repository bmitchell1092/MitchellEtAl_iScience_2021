%% bincontrast_addons. 
%  Adds relevant variables required for detailed analysis. 
% clear; 
% USER = getenv('username');
% cd(strcat('C:/Users/',USER,'/OneDrive - Vanderbilt/Maier Lab/Analysis/MATLAB/workspaces/auto/'));
% load('bincontrast_276_raw_IOS.mat');

%% User selection 
occAdjustment = 1; % this will use a specified method of refining ocular assignment. 
conditions = upper({'de','nde','bin'});
colors     = colorstruct();

%% 1. CRF parameters

fprintf('Running CRF parameters... \n');

% Parameters
[PARAMS] = getPARAMS(UNIT, IDX, info, false);
disp('PARAMS created');
    
[PARAMSbsl] = getPARAMS(UNIT, IDX, info, true);
disp('PARAMSbsl created');

units = 1:242;
[PARAMS2] = getPARAMS_unique(UNIT, IDX, info, false,units);
disp('PARAMS created');

[PARAMS_SUA,SUA_errlog] = getPARAMS(UNIT_SUA,IDX_SUA,info_SUA,false);


%% 2. Ocularity adjustment
% Ocularity correction
if occAdjustment == 1 && info.occSwap == 0
    disp('Checking ocularity values...')
    occSwap;
    fprintf('OccSwap finished. \n');
    
    % Re-generate PARAMS
    disp('Re-generating PARAMS...')
    
    [PARAMS] = getPARAMS(UNIT, IDX, info, false);
    disp('PARAMS re-created');
    
    [PARAMSbsl] = getPARAMS(UNIT, IDX, info, true);
    disp('PARAMS_bsl re-created');
    
    
    if strcmp(method,'rmax')
        disp('Checking for error...')
        delta = PARAMS.DE.a(info.windows-1,:) - PARAMS.NDE.a(info.windows-1,:); % full response window
        idx = delta < tolerance_thresh;
        fprintf('Error = %d \n',sum(idx));
    end
    
end


%% 2. Structures containing relevant, organized data. 

fprintf('Initializing data structures for analysis... \n');

CRF        = getCRF(PARAMS,UNIT,conditions,true, info); fprintf('1) CRF \n');
CRF2       = getCRF(PARAMS2,UNIT,conditions,true,info); fprintf('2) CRF2 \n');
CRF_SUA    = getCRF(PARAMS_SUA, UNIT_SUA, conditions, true, info_SUA);
GAIN_SUA       = fitGain(UNIT_SUA, PARAMS_SUA, IDX_SUA, info_SUA);
OCC        = getOCC(PARAMS2,CRF2, GAIN, info,conditions, 4,'auc'); fprintf('2) OCC \n'); % 'rmax','resp','diUnitTuning'
OCC_SUA        = getOCC(PARAMS_SUA,CRF_SUA, [], info_SUA,conditions, 3,'auc'); fprintf('2) OCC \n'); % 'rmax','resp','diUnitTuning'
%LAY        = getLAY(IDX,PARAMS2,CRF2); fprintf('3) LAY \n');

fprintf('Addons initialized! \n');

%% Creating a modified workspace from the original.

USER = getenv('username');
cd(strcat('C:/Users/',USER,'/OneDrive - Vanderbilt/Maier Lab/Analysis/MATLAB/workspaces/',info.datatype,'/'));
info.filename = strcat(info.filename,'_auc_OA');
save(info.filename,'IDX','PARAMS','PARAMSbsl','UNIT','ERR','CRF','colors','conditions','info'); % add NDE and BIN back if needed
fprintf('Workspace Saved! \n');

%save('bincontrast_final_2','IDX*','PARAMS*','UNIT*','CRF*','info*','colors','conditions','ExpVar'); % add NDE and BIN back if needed




