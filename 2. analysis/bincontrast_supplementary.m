%% Figure_1.m 
% Author: Blake A. Mitchell
% Creation date: 11/3/2020
% Much of Figure 1 is created in illustrator

figDir = 'C:\Users\bmitc_000\OneDrive - Vanderbilt\Maier Lab\Manuscripts\bincontrast\figures\';

%% B - RFs for all units

global ALIGNDIR
if isempty(ALIGNDIR)
    setup('blakeTeba')
end

% directory for V1 Limits
try
    aligndir = ALIGNDIR;
catch
    error('ALIGNDIR not found');
end


dataset = 'bincontrast_12101_08-Oct-2020-250ms';
datatype = 'auto';

USER = getenv('username');

if ~exist('IDX','var')
    load(strcat('C:\Users\',USER,'\OneDrive - Vanderbilt\Maier Lab\MATLAB\workspaces\',datatype,'\',dataset))
end


for i = 1:length(IDX)
list{i} = IDX(i).penetration;
end

list = unique(list);


% Collect putative RFs across all penetrations

for i = 1:length(list)
    load([aligndir list{i} '.mat'],'fRF')
    fRF = rmmissing(fRF);
    
    uCentroid = nanmedian(fRF(:,1:2));
    uWidth    = nanmedian(fRF(:,3:4));
    ecc   = sqrt(sum(uCentroid .^2,2));
    diam   = mean(uWidth,2);
    area   = sqrt(pi .* uWidth(:,1) .* uWidth(:,2));
    
    CENT(:,i) = uCentroid;
    ECC(:,i) = ecc;
end

% Plot
fig = get(groot,'CurrentFigure');
if isempty(fig) 
h = figure('position',[411,455,894,382]); 
end

figure('position',[411,558,356,279]); 
xlim([-7.5 7.5])
ylim([-7.5 7.5]);
hl = hline(0);
vl = vline(0);

for i = 1:size(CENT,2)
    x = CENT(1,i);
    y = CENT(2,i);
    cp = circle(x, y, (.21*(sqrt(x^2 + y^2))) /2); hold on
    if strcmp(dataset,'bincontrast_I')
        set(cp,'Color','r');
    else
        set(cp,'Color','b');
    end
end

axis square

set(gca,'linewidth',0.7,'fontsize',16);
set(vl,'linewidth',0.5,'linestyle','-','color','k');
set(hl,'linewidth',0.5,'linestyle','-','color','k');
xlabel('Horizontal DVA')
ylabel('Vertical DVA')

% Save
cd(figDir);
saveas(gcf, strcat('1-b','.svg'));
disp("Figure saved");



%% C

%% G: Raster plot examples

% user choice
unit = 20; %8, %11 %14 % 15-monocular, 16-binocular
cont = 4; % specify contrast
eye = 1; % 1 = DE, 2 = NDE, 3 = BIN

tw = 1:length(info.sdfWin);clc
cond = {'de','nde','bin'};

clear *binarySpks*
for c = 1:length(IDX(unit).cLevels)
    binarySpks{c}.de(:,:)  = SUA(unit).MON.DE_PS{1,c}(tw,:);
    binarySpks{c}.nde(:,:) = SUA(unit).MON.NDE_PS{1,c}(tw,:);
    binarySpks{c}.bin(:,:) = SUA(unit).BIN.PS{1,c}(tw,:);
end

figure('Position',[560,423,358,146]);
clear spikeTimes spks xspikes yspikes
for iTrials = 1:size(binarySpks{cont}.(cond{eye}),2)
    spikeTimes{iTrials} = info.sdfWin(:,binarySpks{cont}.(cond{eye})(:,iTrials)==1);
end

% Raster plot
MarkerFormat.MarkerSize = 5;
plotSpikeRaster(spikeTimes,'PlotType','scatter','MarkerFormat',MarkerFormat,'XLimForCell',[-0.05 0.250]); hold on

% Onset vertical line
vl = vline(0);
set(vl, 'linewidth',0.4,'color',colors.black)

% beautify
set(gca,'linewidth',1.5,'fontsize',14);
xlabel('Time [s]');
ylabel('Trials');

% Save
cd(figDir);
saveas(gcf, strcat('1-g','.svg'));
disp("Figure saved");


%% H: Example CRF

clear -global Data 

% choices
tw = 3; % time window
error = 1;
unit = 19; % good: 1,2,4, 13, 14, 16, 19, 20, 21

% strings
curves = {'de','nde','bin'};

% contrast levels
if isfield(IDX,'diLevels')
    x = IDX(unit).cLevels*100; x(1) = x(1)+1;
else
    disp('Using average of binned contrast levels'); fprintf('\n')
    x = mean(IDX(unit).cLevels,2)*100; x(1) = x(1)+1;
end

% data
de_all.unit = squeeze(UNIT(unit).MON.DE_PS.RESP);
nde_all.unit = squeeze(UNIT(unit).MON.NDE_PS.RESP);
bin_all.unit = squeeze(UNIT(unit).BIN.PS.RESP);

de_all.err = squeeze(UNIT(unit).MON.DE_PS.RESP_error);
nde_all.err = squeeze(UNIT(unit).MON.NDE_PS.RESP_error);
bin_all.err = squeeze(UNIT(unit).BIN.PS.RESP_error);

global Data
Data(1,:) = x;
figure('position',[891,371,423,420]);
    
    DE = de_all.unit(:,tw);
    de_ref = de_all.unit(:,1);
    NDE = nde_all.unit(:,tw);
    BIN = bin_all.unit(:,tw);
    
    mn       = min(de_ref);
    mx       = max(de_ref);
    nDE      = (DE - mn)./(mx - mn);
    nNDE     = (NDE - mn)./(mx - mn);
    nBIN     = (BIN - mn)./(mx - mn);
    
    for curve = 1:size(curves,2) % for each curve (bin and mon)
        switch curves{curve}
            case 'de'
                Data(2,:) = DE;
            case 'nde'
                Data(2,:) = NDE;
            case 'bin'
                Data(2,:) = BIN;
        end
        
        [a,K,n,b] = BMfitCRdata;
        predictions = 1:100;
        for c = 1:length(predictions) % generate prediction for mon
            prd(curve,c) = a*[(c^n)/((c^n) + (K^n))+ b]; % mon prediction
            eG(curve,c)  = 100 / (K^n + c^n);
            eT(curve,c)  = 1*(1+(2*K)/(2*K)).*sqrt(eG(curve,c));
        end
        
        params(curve,:) = [a, K, n, b];
        
        [~, ~, r2(curve)] = BMNakaRushton(params(curve,:));
        
    end
    
    
    % Plot data points and CRF fit
    semilogx(predictions,prd(1,:),'color',colors.black,'linestyle','-','linewidth',1.5,'HandleVisibility','off');  hold on
    semilogx(Data(1,:),DE,'o','color',colors.black,'linewidth',2,'markersize',5); hold on
    
    % error bars
    ci1 = errorbar(Data(1,:),de_all.unit(:,tw),de_all.err(:,tw));
    set(ci1,'linestyle','none','linewidth',1.5,'color',colors.black,'handleVisibility','off');
    
    % beautify
    set(gca,'FontSize',16,'linewidth',1.5,'box','off',...
        'xlim',[0.5 100],...
        'XTick',x);
    
    ylabel(yL,'fontsize',16);
    xticks([1,5,10,20,50,100]);
    xlabel('contrast','fontsize',16);
    xticklabels({'0','5','10','20','50','100'});
    
    
    % Save
    cd(figDir);
    saveas(gcf, strcat('1-g','.svg'));
    disp("Figure saved");

