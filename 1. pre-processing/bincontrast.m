%% initBincontrast.m
% Loads in ditask units. IDX is the info struct. UNIT and PEN contain data
% Author: Blake M
% Release: 6-28-2021

% Contents: 
% 1. User choice
% 2. Retrieve files
% 3. Pre-loop setup
% 4. Pull out data (by condition)
% 5. Data formatting 
% 6. Data organization (UNIT structure)
% 7. Unit properties (IDX structure)
% 8. General information (info structure)
% 9. Save workspace

clear

%% 1. User choices

% choose dataset
dataset = 'bincontrast';
datatype = 'kls';
animal = {'I','E'};


resptype = 3;               % 1 is normal; 2 is Michele's JOV paper; 3 is 10ms sliding window
dataform = 1;               % 1 = raw, 2 = baseline correct, 3 = tonic correct, 4 = % change, 5 = z-score, 6 = normalize 
normtype = 1;               % 1 = norm to median 95th percentile, 2 = norm to max, 3 = norm DE
cbins    = 2;               % 1 = four nearest levels, 2 = binning, 3 = all levels

% additional options
flag_addblank = 1;          % add spontaneous firing rate as 0 contrast (only to RESP, aids in CRF)?
flag_balanced = 1;          % balance the trial counts
flag_ctuned   = 0;          % get rid of units that are not tuned to contrast
flag_missing  = 0;          % require all contrast levels tested for each unit (necessary for CRF)
flag_spectral  = 0;         % adds multitaper spectral analysis from Chronux toolbox
flag_occCheck  = 1;         % Doesn't swap eye dominance, but does indicate there may be an issue. 
flag_save = 1;              % save the dataset?


%% 2. Retrieve Files

% Setup directory for files of interest
USER = getenv('username');
didir = strcat('C:\Users\',USER,'\OneDrive - Vanderbilt\Maier Lab\Analysis\MATLAB\data\',dataset,'\');


% create a list of data files for this analysis
list    = dir([didir '*_',upper(datatype),'.mat']);

for j = 1: length(list)
    idx(j) = contains(list(j).name, animal);
end

list = list(idx);

%% 3. Pre-loop setup

% Counts
N = 0;
uct = 0;

% constants for all files
sdfWin = -.150:.001:.250;
tw = 1:length(sdfWin);
bsl = [51,151];  % -100 to 0
evoked_w = bsl(2):tw(end); % for spectral analysis

if resptype == 1
    win_idx = [151 251; 301 401; 151 401; bsl];  % 0 - 100ms, 150-250ms, 0 - 250ms
elseif resptype == 2
    win_idx = [201 251; 301 401; 201 401; bsl];  % 50 - 100ms, 150-250ms, 50 - 250ms
elseif resptype == 3
    numWin = ((tw(end) - (bsl(2)+100)) / 10) + 1;
    win_idx = nan(numWin,2); win_idx(1,:) = [151 251]; % this is the first window
    for w = 2:numWin
        next = win_idx(w-1,:) + 10;
        win_idx(w,:) = next;
    end
    win_idx = [win_idx; bsl(2), tw(end); bsl];
elseif resptype == 4
    numWin = ((tw(end) - (bsl(2)+50)) / 20) + 5;
    win_idx = nan(numWin,2); 
    win_idx(1,:) = [201 251]; % this is the first window
    for w = 2:numWin
        next = win_idx(w-1,:) + 10; % 10 is the slide
        win_idx(w,:) = next;
    end
    win_idx = [win_idx; 201, tw(end); bsl];
    
end

respWin = sdfWin(win_idx); 
respWinNum = size(respWin,1);

%% 4. Pull data for each penetration, one electrode at a time

ERR = struct;

% Penetration Loop
for pen = 1:length(list)
    tic
    
    % Load penetration data
    clear penetration
    penetration = list(pen).name(1:11);
    
    load([didir penetration '_',upper(datatype),'.mat'],'STIM')
    matobj = matfile([didir penetration '_',upper(datatype),'.mat']);
    
    N = N+1; % will have a running count of penetrations
    
    switch datatype
        case 'kls'
            eLength = length(STIM.units);
        otherwise
            eLength = length(STIM.depths);
    end
    
    % Electrode Loop
    for e = 1:eLength
        
        uct = uct+1;
        ERR(uct).penetration = STIM.penetration;
        ERR(uct).depth       = STIM.depths(e);
        ERR(uct).message     = [];
        
        % Task selection
        switch datatype
            case 'kls'
                goodtasks = STIM.units(e).fileclust(:,1);
            otherwise
                goodtasks = unique(STIM.filen);
        end
        
        % resp used for diUnitTuning
        sdf = squeeze(matobj.SDF(e,:,:));
        resp = nan(size(win_idx,1),size(sdf,2));
        for w = 1:size(win_idx,1)
            resp(w,:) = nanmean(sdf(win_idx(w,1):win_idx(w,2),:),1);
        end
        resp = squeeze(bsxfun(@minus,resp(end-1,:), resp(end,:)))';% baseline corrects full response (end - 1) by baseline (end)
        

        % Unit Tuning. Decide DE, NDE, PS, and NS
        try
            X = diUnitTuning(resp,STIM,goodtasks); %get tuning info for the unit
        catch
            ERR(uct).message = {'diUnitTuning failed'};
            warning('diUnitTuning failed for unit %d',uct)
            continue
        end
        
        switch datatype
            case 'lfp'  % this entire thing is flipped for LFP, since negative voltage = stronger response.
                NDE = X.dipref(1);  % preferred eye
                DE = X.dinull(1); % non-preferred eye
                NS = X.dipref(2);  % preferred stimulus
                PS = X.dinull(2);  % null stimulus
            otherwise
                DE = X.dipref(1);  % preferred eye
                NDE = X.dinull(1); % non-preferred eye
                PS = X.dipref(2);  % preferred stimulus
                NS = X.dinull(2);  % null stimulus
        end
        
        if isnan(DE) || isnan(NDE) || isnan(PS) || isnan(NS) %%|| all(isnan(X.bio))
            warning('DiUnitTuning passed for unit %d, but produced NaNs',e)
            ERR(uct).message = {'diUnitTuning produced NaNs'};
            continue
        end
        
        if flag_ctuned == true
            if X.dianp(3) > .10
                ERR(uct).message = {'Not tuned to contrast'};
                disp('Unit not tuned to contrast. Next unit...');
                continue
            end
        end
        
        % sort data so that they are [prefeye nulleye]
        clear eyes sortidx contrasts tilts
        eyes      = STIM.eyes;
        contrasts = STIM.contrast;
        tilts     = STIM.tilt;
        if DE == 2
            [eyes,sortidx] = sort(eyes,2,'ascend');
        else
            [eyes,sortidx] = sort(eyes,2,'descend');
        end
        for w = 1:length(eyes)
            contrasts(w,:) = contrasts(w,sortidx(w,:)); % sort contrasts in dominant eye and non-dominant eye
            tilts(w,:)     = tilts(w,sortidx(w,:));
        end; clear w
        
        % establish constant conditions
        I = STIM.ditask ...
            & STIM.adapted == 0 ...           % is not adapted
            & STIM.rns == 0 ...               % not random noise stimulus
            & STIM.cued == 0 ...              % not cued or uncued
            & STIM.motion == 0 ...            % not moving
            & ismember(STIM.filen,goodtasks); % tasks that should be included.
        
        % pull out the data for single electrode
        clear sdf sdftm resp
        sdftm =  matobj.sdftm;
        psthtm =  matobj.psthtm;
        sdf   = squeeze(matobj.SDF(e,:,:));
        
        for w = 1:size(win_idx,1) % for each resp window
            resp(w,:) = nanmean(sdf(win_idx(w,1):win_idx(w,2),:),1); % take the average of each window and place in resp variable
        end
        
        if strcmp('kls',datatype)
            sua = squeeze(matobj.SUA(e,:,:));
        end
        
        
        % Contrast levels
        switch cbins
            case 1
                cLevels = [0, 0; 0.20 0.225; 0.40 0.45; 0.80 0.90]; % four levels: only one contrast per level
            case 2
                cLevels = [0, 0; 0.15, 0.30; 0.40 0.60; 0.80, 1.00];  % four levels: some levels get binned together
            case 3
                cLevels = [0 X.dicontrasts]; % variable levels: all contrast levels included
        end
        
        numC = length(cLevels);
        
        % spectral analysis parameters
        if flag_spectral == 1
            movingwin = [0.04 0.004];
            params.Fs = 1000;
            params.tapers = [0.5 2];
            params.fpass = [0 120];
            params.pad = 2;
            params.trialave = 1;
            params.err = 0;
        end
        
        %% Dioptic Conditions
        
        condition     = {'BIN','DE','NDE'};
        trl_index     = false(1,size(resp,2))';
        trl_count     = nan(numC,3);  % contrast x condition
        
        for cond = 1:size(condition,2) % for each condition
            for c = 1:length(cLevels) % for each contrast level
                switch condition{cond}
                    case 'BIN'
                        if c == 1
                            trls = STIM.blank;
                        else
                            switch cbins
                                case {1,2}
                                    trls = I & STIM.botheyes...
                                        & contrasts(:,1) >= cLevels(c,1) & contrasts(:,1) <= cLevels(c,2)... % contrast in de
                                        & contrasts(:,2) >= cLevels(c,1) & contrasts(:,2) <= cLevels(c,2)... % contrast in nde
                                        & tilts(:,1) == PS... % pref orientation in de
                                        & tilts(:,2) == PS; % pref orientation in nde 
                                case 3
                                    trls = I & STIM.botheyes...
                                        & contrasts(:,1) == cLevels(c)... % contrast in de
                                        & contrasts(:,2) == cLevels(c)... % contrast in nde
                                        & tilts(:,1) == PS... % pref orientation in de
                                        & tilts(:,2) == PS; % pref orientation in nde
                            end
                        end       
                        
                    case 'DE'
                        if c == 1
                            trls = STIM.blank; % zero contrast in both eyes
                        else
                            switch cbins
                                case {1,2}
                                    trls = I & STIM.monocular & DE... % is monocular and dominant eye
                                        & contrasts(:,1) >= cLevels(c,1) & contrasts(:,1) <= cLevels(c,2)... % contrast in dom eye
                                        & tilts(:,1) == PS; % pref orientation in dom eye
                                case 3
                                    trls = I & STIM.monocular & DE... % is monocular and dominant eye
                                        & contrasts(:,1) == cLevels(c)... % contrast in dom eye
                                        & tilts(:,1) == PS; % pref orientation in dom eye
                            end
                        end
                        
                        if flag_balanced == true
                            n = trl_count(c,1); % number of trials to keep
                            f = find(trls); % find the location of the logical 1's in monocular trials
                            f = f(randperm(numel(f))); % randomize the find results
                            trls(f(n+1:end)) = false; % get rid of the other trials beyond the random n
                        end
                        
                    case 'NDE'
                        if c == 1
                            trls = STIM.blank;
                        else
                            switch cbins
                                case {1,2}
                                    trls = I & STIM.monocular & NDE...  % is monocular and non-dominant eye
                                        & contrasts(:,2) >= cLevels(c,1) & contrasts(:,2) <= cLevels(c,2)... % contrast in ndom eye
                                        & tilts(:,2) == PS; % pref orientation in non-dom eye
                                case 3
                                    trls = I & STIM.monocular & NDE...  % is monocular and non-dominant eye
                                        & contrasts(:,2) == cLevels(c) ... % contrast in non-dom eye
                                        & tilts(:,2) == PS; % pref orientation in non-dom eye
                            end
                        end
                        
                        if flag_balanced == true 
                            n = trl_count(c,1); % number of trials to keep
                            f = find(trls); % find the location of the logical 1's in monocular trials
                            f = f(randperm(numel(f))); % randomize the find results
                            trls(f(n+1:end)) = false; % get rid of the other trials beyond the random n
                        end
                end
                
                
                % remove nans
                [~, sdftrls] = rmnantrls(trls, sdf, tw);
                [goodtrls, resptrls] = rmnantrls(trls, resp);
                
                % organize
                SDFtrls.(condition{cond}){c} = sdftrls;
                RESPtrls.(condition{cond}){c} = resptrls;
                
                % track
                trl_index(goodtrls) = true;
                trl_count(c,cond) = size(goodtrls,1);
                
                % spectral analysis
                if flag_spectral == 1
                    [~,data] = rmnantrls(trls,sdf, evoked_w);
                    
                    [SPEC_S{cond}(c,:,:), SPEC_t{cond}(c,:), SPEC_f{cond}(c,:)] = mtspecgramc(data,movingwin,params);
                    [FREQ_S{cond}(c,:,:), FREQ_f{cond}(c,:)] = mtspectrumc(data,params);
                end
                
                % single-unit activity (for KLS)
                if strcmp(datatype,'kls')
                    SUA.(condition{cond}){c}  = sua(tw,goodtrls);
                end
                
                clear trls sdftrls resptrls 
            end
        end
        
        % (optional) checking that all conditions are present
        if flag_missing == true
            switch dataform
                case 1
                    threshold = 0; %
                otherwise
                    threshold = 4;
            end
            if any(trl_count(2:end,:) <= threshold)
                ERR(uct).message = {'Missing trials/conditions'};
                disp('Missing trials/conditions. Next unit...');
                continue
            end      
        end
        
        % always check if there are atleast binocular AND monocular trials
        if ~any(trl_count(2:end,1)) 
            ERR(uct).message = {'Missing trials/conditions'};
            disp('Missing trials/conditions. Next unit...');
            continue
        end
        
        
        %% 5. all contrast combinations
        
        diLevels = [0.0, 0.0,   0.0, 0.0;
                    0.0, 0.0,   0.2, 0.225;
                    0.0, 0.0,   0.4, 0.45;
                    0.0, 0.0,   0.8, 0.9;
                    0.2, 0.225, 0.0, 0.0;
                    0.2, 0.225, 0.2, 0.225; 
                    0.2, 0.225, 0.4, 0.45;
                    0.2, 0.225, 0.8, 0.9;
                    0.4, 0.45,  0.0, 0.0;
                    0.4, 0.45,  0.2, 0.225;
                    0.4, 0.45,  0.4, 0.45;
                    0.4, 0.45,  0.8, 0.9;
                    0.8, 0.9,   0.0, 0.0;
                    0.8, 0.9,   0.2, 0.225;
                    0.8, 0.9,   0.4, 0.45;
                    0.8, 0.9,   0.8, 0.9]; 
                
                
        % Determine which levels were shown to this unit
        if any(X.dicontrasts == 0.22)
            diType = diLevels(:,[2,4]);
        else 
            diType = diLevels(:,[1,3]);
        end
        
        numC = length(diLevels);       
        di_trl_count     = nan(numC,1);  % total contrast conditions  
        
        for c = 1:length(diLevels) % for each contrast level
            if c == 1
                trls = STIM.blank;
            elseif c == 2 || c == 3 || c == 4
                trls = I & STIM.monocular & NDE...
                    & contrasts(:,2) >= diLevels(c,3) & contrasts(:,2) <= diLevels(c,4)... % contrast in ndom eye
                    & tilts(:,2) == PS; % pref orientation in null eye
                balance = true;
                
            elseif c == 5 || c == 9 || c == 13
                trls = I & STIM.monocular & DE...
                    & contrasts(:,1) >= diLevels(c,1) & contrasts(:,1) <= diLevels(c,2)... % contrast in ndom eye
                    & tilts(:,1) == PS; % pref orientation in null eye
                balance = true;
            else
                trls = I & STIM.botheyes...
                    & contrasts(:,1) >= diLevels(c,1) & contrasts(:,1) <= diLevels(c,2)... % contrast in dom eye
                    & contrasts(:,2) >= diLevels(c,3) & contrasts(:,2) <= diLevels(c,4)... % contrast in ndom eye
                    & tilts(:,1) == PS... % pref orientation in dom eye
                    & tilts(:,2) == PS; % pref orientation in null eye
                balance = false;
            end
            
%             if flag_balanced == true && balance == true
%                 n = trl_count(c,1); % number of trials to keep
%                 f = find(trls); % find the location of the logical 1's in monocular trials
%                 f = f(randperm(numel(f))); % randomize the find results
%                 trls(f(n+1:end)) = false; % get rid of the other trials beyond the random n
%             end


                % remove nans
                [~, sdftrls] = rmnantrls(trls, sdf, tw);
                [goodtrls, resptrls] = rmnantrls(trls, resp);
                
                % organize
                SDFtrls.DI{c} = sdftrls;
                RESPtrls.DI{c} = resptrls;
                
                % track
                trl_index(goodtrls) = true;
                di_trl_count(c) = size(goodtrls,1);
                
                % spectral analysis
                if flag_spectral == 1
                    [~,data] = rmnantrls(trls,sdf, evoked_w);
                    
                    [SPEC_S{DI}(c,:,:), SPEC_t{DI}(c,:), SPEC_f{DI}(c,:)] = mtspecgramc(data,movingwin,params);
                    [FREQ_S{DI}(c,:,:), FREQ_f{DI}(c,:)] = mtspectrumc(data,params);
                end
                
                % single-unit activity (for KLS)
                if strcmp(datatype,'kls')
                    SUA.DI{c}  = sua(tw,goodtrls);
                end
                
                clear trls sdftrls resptrls 
        end

        
        % (optional) checking that all conditions are present
        if flag_missing == true
            switch dataform 
                case 1
                    threshold = 4;
                otherwise
                    threshold = 4;
            end
            if any(di_trl_count(1:end,:) <= threshold)
                ERR(uct).message = {'Missing trials/conditions'};
                disp('Missing DI trials/conditions. Logged but not skipped');
%                 continue
            end
        end
        
        
        
        
        %% 5. Data form        
        switch dataform
            case 1 % raw
                for cond = 1:size(condition,2)
                    if flag_addblank == true
                        blank = resp(respWinNum,trl_index); % adds to last window: respWinNum
                        RESPtrls.(condition{cond}){1,1} = repmat(blank,respWinNum,1); 
                        RESPtrls.DI{1,1} = repmat(blank,respWinNum,1);
                    end
                end
            case 2 % baseline correct
                for cond = 1:size(condition,2)
                    
                    if flag_addblank == true
                        blank = resp(respWinNum,trl_index);
                        RESPtrls.(condition{cond}){1,1} = repmat(blank,length(cLevels),1);
                    end
                    
                    for c = 1:length(cLevels)
                        r = RESPtrls.(condition{cond}){c};
                        r = bsxfun(@minus,r, r(respWinNum,:)); % last column in resp matrix is baseline
                        
                        s  = SDFtrls.(condition{cond}){c}; 
                        s   = bsxfun(@minus,s, nanmean(s(bsl(1):bsl(2),:),1)); % bsl is defined earlier
                        
                        s(s < 0) = 0;
                        r(r < 0) = 0; % half wave rectify
                        
                        RESPtrls.(condition{cond}){c} = r;
                        SDFtrls.(condition{cond}){c} = s;
                    end
                end
                
            case 3 % tonic correcting
                for cond = 1:size(condition,2)
                    
                    if flag_addblank == true
                        blank = resp(respWinNum,trl_index);
                        RESPtrls.(condition{cond}){1,1} = repmat(blank,length(cLevels),1);
                    end
                    
                    tonic = nanmean(resp(respWinNum,trl_index));
                    if tonic == 0
                        tonic = 1; % half wave rectify
                    end
                    
                    for c = 1:length(cLevels)
                        r = RESPtrls.(condition{cond}){c};
                        r = bsxfun(@minus,r, tonic);
                        
                        s  = SDFtrls.(condition{cond}){c}; % last column in resp matrix is baseline
                        s   = bsxfun(@minus,s, tonic); % bsl is defined earlier
                        
                        s(s < 0) = 0;
                        r(r < 0) = 0;
                        
                        RESPtrls.(condition{cond}){c} = r;
                        SDFtrls.(condition{cond}){c} = s;
                    end
                end
                
            case 4 % percent change from tonic firing rate
                for cond = 1:size(condition,2)
                    
                    if flag_addblank == true
                        blank = resp(respWinNum,trl_index);
                        RESPtrls.(condition{cond}){1,1} = repmat(blank,length(cLevels),1);
                    end
                    
                    % calculate tonic firing rate
                    tonic = nanmean(resp(respWinNum,trl_index)); % takes average of all trial baselines for this unit
                    if tonic == 0
                        tonic = 1;
                    end
                    
                    for c = 1:length(cLevels)
                        % define percent change function
                        fun = @(a,b) ((a - b) / b) * 100;
                        r = RESPtrls.(condition{cond}){c};
                        r = bsxfun(fun, r, tonic);
                        
                        s = SDFtrls.(condition{cond}){c};
                        s = bsxfun(fun, s, tonic);
                        
                        s(s < 0) = 0;
                        r(r < 0) = 0;
                        
                        RESPtrls.(condition{cond}){c} = r;
                        SDFtrls.(condition{cond}){c} = s;
                    end
                end
                
            case 5 % Z-score (currently just normalized to sigma)
                
                for cond = 1:size(condition,2)
                    
                    if flag_addblank == true
                        blank = resp(respWinNum,trl_index);
                        RESPtrls.(condition{cond}){1,1} = repmat(blank,length(cLevels),1);
                    end
                    
                    for c = 1:length(cLevels)
                        
                        mu = nanmean(resp(respWinNum-1,trl_index)); % average window is respWinNum - 1
                        sigma = nanstd(resp(respWinNum-1,trl_index),0);
                        
                        if sigma == 0
                            sigma = 0.1;
                        end
                        
                        r = RESPtrls.(condition{cond}){c};
                        %r = bsxfun(@minus, r, mu);
                        r = bsxfun(@rdivide, r, sigma);
                        
                        s = SDFtrls.(condition{cond}){c};
                        %s = bsxfun(@minus, s, mu);
                        s = bsxfun(@rdivide, s, sigma);
                        
                        RESPtrls.(condition{cond}){c} = r;
                        SDFtrls.(condition{cond}){c} = s;
                    end
                    
                end
                
            case 6 % Normalize to some central tendency
                
                % 1. Define central tendency
                tonic = nanmean(resp(respWinNum,trl_index));
                if tonic == 0
                    tonic = 1;
                end
                
                switch normtype
                    case 1
                        ct = median(prctile(resp(respWinNum-1,trl_index),95)) - tonic;
                    case 2
                        ct = max(resp(respWinNum-1,trl_index)) - tonic;
                    case 3
                        ct = max(RESPtrls.DE{1,4}(3,:)) - tonic; % this needs to be checked
                end
                
                % 2. Loop through organized data
                for cond = 1:size(condition,2)
                    
                    if flag_addblank == true
                        blank = resp(respWinNum,trl_index);
                        RESPtrls.(condition{cond}){1,1} = repmat(blank,length(cLevels),1);
                    end
                    
                    for c = 1:length(cLevels)
                        
                        % Remove baseline
                        r = RESPtrls.(condition{cond}){c};
                        %r = bsxfun(@minus,r, r(4,:));
                        r = bsxfun(@minus,r, tonic);
                        
                        s  = SDFtrls.(condition{cond}){c};
                        %s   = bsxfun(@minus,s, nanmean(s(bsl(1):bsl(2),:),1)); % bsl is defined earlier
                        s   = bsxfun(@minus,s, tonic);
                        
                        % Normalize by central tendency (ct)
                        r = bsxfun(@rdivide, r, ct);
                        s = bsxfun(@rdivide, s, ct);
                        
                        RESPtrls.(condition{cond}){c} = r;
                        SDFtrls.(condition{cond}){c} = s;
                    end
                    
                end
        end
        
        % for clarity
        if dataform ~= 6
            normtype = 0;
        end
        
        %% 6. Data organization: Unit structure
        
        % Unit structure
        for cond = 1:size(condition,2)
            for c = 1:length(cLevels)
                UNIT(uct).(condition{cond}).SDF_trls{1,c}  = SDFtrls.(condition{cond}){1,c};
                UNIT(uct).(condition{cond}).SDF_avg(c,:) = nanmean(SDFtrls.(condition{cond}){1,c},2);
                UNIT(uct).(condition{cond}).SDF_sd(c,:)  = nanstd(SDFtrls.(condition{cond}){1,c},[],2);
                UNIT(uct).(condition{cond}).SDF_error(c,:) = nanstd(SDFtrls.(condition{cond}){1,c},[],2) ./ ...
                    sqrt(size(SDFtrls.(condition{cond}){1,c},2));
                
                UNIT(uct).(condition{cond}).RESP_trls{1,c}  = RESPtrls.(condition{cond}){1,c};
                UNIT(uct).(condition{cond}).RESP_avg(c,:) = nanmean(RESPtrls.(condition{cond}){1,c},2);
                UNIT(uct).(condition{cond}).RESP_sd(c,:)  = nanstd(RESPtrls.(condition{cond}){1,c},[],2);
                UNIT(uct).(condition{cond}).RESP_error(c,:) = nanstd(RESPtrls.(condition{cond}){1,c},[],2) ./ ...
                    sqrt(size(RESPtrls.(condition{cond}){1,c},2));
            end
        end
        
  
        for c = 1:length(diLevels)
            UNIT(uct).DI.SDF_trls{1,c}  = SDFtrls.DI{1,c};
            UNIT(uct).DI.SDF_avg(c,:) = nanmean(SDFtrls.DI{1,c},2);
            UNIT(uct).DI.SDF_sd(c,:)  = nanstd(SDFtrls.DI{1,c},[],2);
            UNIT(uct).DI.SDF_error(c,:) = nanstd(SDFtrls.DI{1,c},[],2) ./ ...
                sqrt(size(SDFtrls.DI{1,c},2));
            
            UNIT(uct).DI.RESP_trls{1,c}  = RESPtrls.DI{1,c};
            UNIT(uct).DI.RESP_avg(c,:) = nanmean(RESPtrls.DI{1,c},2);
            UNIT(uct).DI.RESP_sd(c,:)  = nanstd(RESPtrls.DI{1,c},[],2);
            UNIT(uct).DI.RESP_error(c,:) = nanstd(RESPtrls.DI{1,c},[],2) ./ ...
                sqrt(size(RESPtrls.DI{1,c},2));
        end
  
        
        % SUA (kls specific)
        if strcmp(datatype,'kls')
            for cond = 1:size(condition,2)
                for c = 1:length(cLevels)
                    UNIT(uct).(condition{cond}).SUA_trls{c}  = SUA.(condition{cond}){c};
                end
            end
        end
        
        
        %% 7. Unit properties: IDX structure
        
        IDX(uct).penetration = STIM.penetration;
        
        if contains('kls',datatype)
            IDX(uct).clusterID = STIM.kls.cluster(e);
            IDX(uct).depth = STIM.units(e).depth';  % This could use some sorting by depth
            IDX(uct).rate  = STIM.units(e).rate;
        else
            IDX(uct).depth = STIM.depths(e,:)';
        end
        
        % layer identification
        if STIM.depths(e,2) >= 5
            layer = 1;
        elseif STIM.depths(e,2) < 5 && STIM.depths(e,2) >= 0
            layer = 2;
        elseif STIM.depths(e,2) < 0
            layer = 3;
        end
        
        IDX(uct).layer      = layer;
        IDX(uct).v1lim      = STIM.v1lim;
        IDX(uct).X          = X;
        IDX(uct).effects    = X.dianp; % p for main effect of each 'eye' 'tilt' 'contrast'    
        
        IDX(uct).DE         = DE;
        IDX(uct).NDE        = NDE;
        IDX(uct).occ        = X.occ';    % how much it prefers one eye over the other
        IDX(uct).occana     = X.occana;
        
        IDX(uct).prefori    = PS;
        IDX(uct).nullori    = NS;
        IDX(uct).ori        = X.ori';    % how much it prefers one orientation over the other
        IDX(uct).oriana     = X.oriana;

        IDX(uct).binocularity   = X.bio';    % How much it prefers both eyes over one
        IDX(uct).diana          = X.diana;

        
        IDX(uct).SDFlength     = length(matobj.sdftm);
        IDX(uct).cLevels       = cLevels;
        IDX(uct).diLevels      = diType;
        IDX(uct).trial_index   = trl_index;
        IDX(uct).trial_count   = trl_count;
        IDX(uct).di_trial_count = di_trl_count;
        IDX(uct).occCheck      = false;
        
        toc
        
        
        
    end % end electrode loop
    
end % end penetration loop

%% 8. General: info structure

dataformstring = {'raw','bsl','norm','norm2'};

info = struct;
info.filename       = [];
info.dataset        = dataset;
info.animal         = animal;
info.datatype       = datatype;
info.N              = N;
info.uct            = uct;
info.sdfWin         = sdfWin;
info.respWin        = respWin;
info.windows        = respWinNum;
info.cbinning       = cbins;
info.cLevels        = cLevels;
info.diLevels       = diLevels;
info.resptype       = resptype;
info.windows        = respWinNum;
info.dataform       = char(dataformstring(dataform));
info.normtype       = normtype;
info.ctuned         = flag_ctuned;
info.allTrials      = flag_missing;
info.balanced       = flag_balanced;
info.spectral       = flag_spectral;


% get rid of units that are deemed "bad"
if exist('IDX','var') 
    IDX( all( cell2mat( arrayfun( @(x) structfun( @isempty, x ), IDX, 'UniformOutput', false ) ), 1 ) ) = [];
    UNIT( all( cell2mat( arrayfun( @(x) structfun( @isempty, x ), UNIT, 'UniformOutput', false ) ), 1 ) ) = []; 
end

% re-calculate N and uct after matrix cleanup
info.uct = length(IDX);
[uc, ~, idc] = unique({IDX.penetration});
info.N   = length(uc);
info.occSwap = 0;

if flag_occCheck == true
    
    % counting # of units that exceed a threshold of ocular dominance identification error
    count = 0;
    thresh = 0; 
    % identify problematic ocularity
    for i = 1:length(IDX)
        sumDE = UNIT(i).DE.RESP_avg(2,3) + UNIT(i).DE.RESP_avg(3,3) + UNIT(i).DE.RESP_avg(4,3);
        sumNDE = UNIT(i).NDE.RESP_avg(2,3) + UNIT(i).NDE.RESP_avg(3,3) + UNIT(i).NDE.RESP_avg(4,3) ;
 
        if (sumDE - sumNDE) < thresh
            count = count + 1;
            IDX(i).occCheck = true;
            fprintf('Unit %d flagged for ocularity check',i); fprintf('\n')
        else
            IDX(i).occCheck = false;
        end
    end
    
    info.flaggedOcc = count;
end

flag_code = [resptype,dataform,normtype,cbins,flag_addblank,flag_balanced,flag_ctuned, flag_missing];

clearvars -except IDX UNIT ERR info flag_save flag_code
     

%% 9. SAVE WORKSPACE

% Figure save directory

USER = getenv('username');
figDir = strcat('C:\Users\',USER,'\OneDrive - Vanderbilt\Maier Lab\Analysis\Figures\');


if flag_save == true
    USER = getenv('username');
    cd(strcat('C:/Users/',USER,'/OneDrive - Vanderbilt/Maier Lab/Analysis/MATLAB/workspaces/',info.datatype,'/'));
    info.filename = strcat(sprintf('%s_',string(info.uct)),'_slidingwindows_I_2'); % IOS for interocular suppression analysis
    save(info.filename,'IDX','ERR','UNIT','info','flag_code','figDir');
    fprintf('Workspace saved\n');
end

fprintf('bincontrast.m has finished. \n');

