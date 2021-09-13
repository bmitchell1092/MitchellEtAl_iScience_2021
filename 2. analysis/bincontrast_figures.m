%% Bincontrast Figures
% Original code for Figures 1 - 5
% Author: Blake Mitchell 
% Date: 9/13/2021

set(0,'DefaultFigureWindowStyle','normal') 
flag_save = 1;

if flag_save == 1 
    USER = getenv('username');
    figDir = strcat('C:\Users\',USER,'\OneDrive - Vanderbilt\Maier Lab\Analysis\Figures\bincontrast\');
    cd(figDir);
    disp("Figure directory created");
end

clear SAMPLE*
% initialize workspace with sample to isolate

SAMPLE.UNIT = UNIT;
SAMPLE.IDX = IDX;
SAMPLE.PARAMS = PARAMS;
SAMPLE.CRF = CRF;
SAMPLE.GAIN = GAIN;
SAMPLE.info = info;
uct = length(SAMPLE.UNIT);

%% Figure 1. Experimental design and sublinear binocular responses
% A. Stereoscope
% B. Stimulus conditions
% C. Single unit example of V1 responses
% D. Sublinearity of binocular responses

%% 1C. Example Unit

clear *MarkerWidth*
unit = 52; % 1,2,4, 11, 13, 14, 16, 19, 20, 21
ct = 0;
tw = 1:length(info.sdfWin);
cond = {'de','nde','bin'};

clear *binarySpks*

for c = 2:length(IDX(unit).cLevels)
    binarySpks{c-1}.de(:,:)  = UNIT_SUA(unit).DE.SUA_trls{1,c}(tw,:);
    binarySpks{c-1}.nde(:,:)  = UNIT_SUA(unit).NDE.SUA_trls{1,c}(tw,:);
    binarySpks{c-1}.bin(:,:) = UNIT_SUA(unit).BIN.SUA_trls{1,c}(tw,:);
end

figure('Units','normalized','Position',[0.1875,0.465740740740741,0.296875,0.250462962962963]);

for eye = 1:3  % for DE, NDE, and BIN
    
    for c = 1:size(binarySpks,2)
        
       
        ct = ct+1;
        
        clear spikeTimes spks xspikes yspikes
        for iTrials = 1:size(binarySpks{c}.(cond{eye}),2)
            spikeTimes{iTrials} = info.sdfWin(:,binarySpks{c}.(cond{eye})(:,iTrials)==1);
        end
        
        % 2. Raster plot 
        ax = subplot(3,size(binarySpks,2),ct);
        MarkerFormat.MarkerSize = 3.5;
        plotSpikeRaster(spikeTimes,'PlotType','scatter','MarkerFormat',MarkerFormat,'XLimForCell',[-0.1 0.250]); hold on
        vl = vline(0);
        set(vl, 'linewidth',0.4,'color',colors.black)
        set(gca,'linewidth',1);
        

        ax.XLabel.String  	= 'Time (s)';
        ax.YLabel.String  	= 'Trials';

        xticklabels([]);
        xlabel([]);
        ylabel([]);
        yticklabels([]);
  
        
        set(gca,'tickDir','out')
        box on
        
    end  
end


cd(figDir);
saveas(gcf, strcat('1C_raster',IDX_SUA(unit).penetration(end-1:end),num2str(IDX_SUA(unit).clusterID),'.svg'));
disp("Figure saved");

%% 1C. SDF overlay

sua_ID = 52; % 1,2,4, 11, 13, 14, 16, 19, 20, 21
mua_ID = 222;
ct = 0;

conditions = {'BIN','DE','NDE'};

color = [colors.blue; colors.cyan; colors.dark_red];

for cond = 1:length(conditions)
    dataA.(conditions{cond}) = squeeze(UNIT_SUA(sua_ID).(conditions{cond}).SDF_avg(2:4,:));
    dataB.(conditions{cond}) = squeeze(UNIT(mua_ID).(conditions{cond}).SDF_avg(2:4,:));
end

dataA.de_err = squeeze(UNIT_SUA(sua_ID).DE.SDF_error(2:4,:));
dataA.nde_err = squeeze(UNIT_SUA(sua_ID).NDE.SDF_error(2:4,:));
dataA.bin_err = squeeze(UNIT_SUA(sua_ID).BIN.SDF_error(2:4,:));

dataB.de_err = squeeze(UNIT(mua_ID).DE.SDF_error(2:4,:));
dataB.nde_err = squeeze(UNIT(mua_ID).NDE.SDF_error(2:4,:));
dataB.bin_err = squeeze(UNIT(mua_ID).BIN.SDF_error(2:4,:));

figure('Units','normalized','Position',[0.1875,0.465740740740741,0.296875,0.250462962962963]);
cond = {'DE','NDE','BIN'; 'de_err','nde_err','bin_err'};
for eye = 1:3
    for c = 1:3
        
        ct = ct+1;
        ax = subplot(size(cond,2),length(info.cLevels)-1,ct);
        
        % plot SUA
        plot(info.sdfWin,dataA.(cond{1,eye})(c,:),'color',[color(eye,:),0.2],'linewidth',1.5); hold on    
        ci1 = ciplot(dataA.(cond{1,eye})(c,:)+dataA.(cond{2,eye})(c,:),dataA.(cond{1,eye})(c,:)-dataA.(cond{2,eye})(c,:),info.sdfWin,color(eye,:),0.3);
        
        % plot MUA
        plot(info.sdfWin,dataB.(cond{1,eye})(c,:),'color',color(eye,:),'linewidth',1.5); hold on    
        ci2 = ciplot(dataB.(cond{1,eye})(c,:)+dataB.(cond{2,eye})(c,:),dataB.(cond{1,eye})(c,:)-dataB.(cond{2,eye})(c,:),info.sdfWin,color(eye,:),0.3);
        
        %lin_sum = dataA.de(3,:)+dataA.nde(3,:);
        limit = max(dataB.BIN(3,:));
        
        % beautify
        set(ci1,'linestyle','none','handleVisibility','off');
        set(ci2,'linestyle','none','handleVisibility','off');
        set(gca,'box','off','linewidth',1);
        ylim([-1 limit*1.3]); 
        xlim([-.1 .250]);
        
        if ct == 1 || ct == 4 || ct == 7
            yticks('auto')
            xticklabels([]);
        else
            xticklabels([]); yticklabels([]);
            xlabel([]); ylabel([]);
        end
        
        if ct > 6 
            xticklabels('auto');
        end
        
        % onset line
        vl = vline(0);
        set(vl, 'linewidth',0.3,'color',colors.black)
        set(gca,'tickDir','out')
        box on
        
    end
end
    
cd(figDir);
saveas(gcf, strcat('1C_sdf',IDX_SUA(unit).penetration(end-1:end),num2str(IDX_SUA(unit).clusterID),'.svg'));
disp("Figure saved");

%% 1D. Depth of sublinearity
clear de* nde* bin* SUM AI

method = 'raw';
tw = info.windows-1;
% Observed responses
switch method
    case 'raw'
        for i = 1:uct
            for c = 1:4
                de_units(c,i)  = SAMPLE.UNIT(i).DE.RESP_avg(c,tw);
                nde_units(c,i) = SAMPLE.UNIT(i).NDE.RESP_avg(c,tw);
                bin_units(c,i) = SAMPLE.UNIT(i).BIN.RESP_avg(c,tw);
            end
        end
        
    case 'percent change'
        % percent increase
        for i = 1:uct
            for c = 1:4
                de_units(c,i)  = ((SAMPLE.UNIT(i).DE.RESP_avg(c,tw) - SAMPLE.UNIT(i).DE.RESP_avg(c,end)) / (SAMPLE.UNIT(i).DE.RESP_avg(c,tw))) * 100;
                nde_units(c,i) = ((SAMPLE.UNIT(i).NDE.RESP_avg(c,tw) - SAMPLE.UNIT(i).NDE.RESP_avg(c,end)) / (SAMPLE.UNIT(i).NDE.RESP_avg(c,tw))) * 100;
                bin_units(c,i) = ((SAMPLE.UNIT(i).BIN.RESP_avg(c,tw) - SAMPLE.UNIT(i).BIN.RESP_avg(c,end)) / (SAMPLE.UNIT(i).BIN.RESP_avg(c,tw))) * 100;
            end
        end
end

% Prediction based on linear sum
for i = 1:uct
    for c = 1:4
       SUM.units(c,i) = de_units(c,i)+nde_units(c,i);
    end
end

for i = 1:uct
    for c = 1:4
       AI.units(c,i) = bin_units(c,i) ./ SUM.units(c,i);
    end
end

% Plot
xlimit = [0 400];
ylimit = [0 250];
sz = 10;
x1 = linspace(0,250);
y1 = linspace(0,400);
labels = {'Low','Medium','High'};

figure('position',[359,599.6666666666666,711.3333333333333,205.3333333333334]);

for c = 1:3
    subplot(1,3,c)
    plot(x1,y1,'linewidth',1.5,'color','k'); hold on
    
    x = SUM.units(c+1,~isnan(SUM.units(c+1,:)))'; y = bin_units(c+1,~isnan(bin_units(c+1,:)))';
    X = [ones(length(x),1) x];
    b = X\y;
    b1 = x\y; slope(c) = b1;
    yCalc1 = b1*x;
    yCalc2 = X*b;
    
    s = scatter(x,y,sz,colors.black,'filled'); hold on
    set(gca,'linewidth',1.5,'fontsize',12,'xlim',xlimit,'ylim',ylimit,'box','off');
    plot(x,yCalc2,'linewidth',2,'color',colors.light_orange,'linestyle','-');
    xlabel([]);
    s.MarkerFaceAlpha = 0.5;
    title(labels{c});
    
    if c > 1
        yticklabels([]);
    end
    
    axis square
end

jamovi.f1D.additivity = AI.units';
jamovi.f1D.bin = bin_units';
jamovi.f1D.sum = SUM.units'; 

if exist('figDir','var') 
    cd(figDir);
    saveas(gcf, strcat('1D_sublinearity','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end


% Jamovi setup 

jamovi.fig1.sum = AI.units(2:4,:)';
additivity = jamovi.fig1.sum;


%% Figure 2. Main effect of binocular modulation
% B. Grand average binocular modulation
% C. Contrast-dependent binocular modulation: scatter with sword-tips

%% 2A: Grand average binocular modulation
tw = info.windows-1;

clear A B RFII
for u = 1:uct
    A(:,:,u) = SAMPLE.UNIT(u).DE.RESP_avg(2:end,1:tw); %-PARAMS.B(tw,u);
    B(:,:,u) = SAMPLE.UNIT(u).BIN.RESP_avg(2:end,1:tw); %-PARAMS.B(tw,u);
end

A = squeeze(nanmean(A,1));
B = squeeze(nanmean(B,1));

clear RFII*
for u = 1:uct
    
    for t = 1:size(A,1)
        RFII.units(t,u) = (B(t,u) - A(t,u)) ./ (A(t,u) + B(t,u));
    end
    
end

p = 95;
RFII.avg = nanmean(RFII.units,2);
RFII.sem = nanstd(RFII.units,[],2)./sqrt(size(RFII.units,2));
RFII.CI = CIFcn(RFII.units(tw,:),p);

tw = info.windows-1;
sz = 30;

% random index 
r = randperm(uct)'; r(:,2) = nan;
for i = 1:length(r)
    if r(i,1) <= count([IDX_SUA.penetration],'E','IgnoreCase',false) % # number of units from E
        r(i,2) = 1;
    else
        r(i,2) = 2;
    end
end

rRFII = RFII.units(tw,r(:,1));
rRFII(2,:) = r(:,2);
E = rRFII(2,:)==1; E = E';
I = rRFII(2,:)==2; I = I';
C = nan(length(rRFII),3);
for i = 1:length(C)
    if rRFII(2,i) == 1
        C(i,:) = colors.green;
    elseif rRFII(2,i) == 2
        C(i,:) = colors.orange;
    end
end


close all;
figure('position',[590.3333333333333,817.6666666666666,522,313.9999999999999]); hold on
set(gca,'linewidth',1.5,'fontsize',14,'ylim',[-0.3 0.3],'tickdir','out'); 
yticks([-0.3,0,0.3]);

% data
s=scatter(1:size(rRFII,2),rRFII(1,:),sz,C,'filled','linewidth',1.5); hold on
s.MarkerFaceAlpha = 0.5;

% midline
hl = hline(0); set(hl,'linestyle','-','linewidth',1,'color','k'); hold on

if flag_save == 1
    saveas(gcf, strcat('2A_main','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

data = RFII.units(tw,:);
up = sum(data > 0); down = sum(data < 0);
percent_up = (up / length(UNIT))*100;
percent_down = (down / length(UNIT))*100;

jamovi.f2A.BMI_grand = RFII.units(tw,:)';
%jamovi.f2A.BMI_grand(1:length(RFII.units),2) = 0;

%% 2A - Histogram

close all;
figure('position',[1000,918,198.3333333333333,420]); 

nbins = 25;
mean_mark = 75;
tw = info.windows-1;
     
plot(mean_mark,RFII.avg(tw,:),'d','linewidth',3,'color',colors.black+0.4); hold on

histogram(RFII.units(tw,:),nbins,'binWidth',0.025,'orientation','horizontal',...
    'faceAlpha',0.5,'linewidth',1.2,'faceColor',colors.black+0.4,'edgealpha',0.3); hold on

set(gca,'linewidth',1.5,'fontsize',14,'ylim',[-.25 0.25],'xlim',[0 mean_mark],'tickdir','out','box','off');
hl = hline(0); set(hl,'linestyle','-','linewidth',1,'color','k');
yticks([-0.25,0,0.25]); xticks([0, mean_mark]);
% Identify CI
L = BMI.CI(c,1);
U = BMI.CI(c,2);

% fill the 95% CI
x2 = get(gca,'xlim');
ypoints = [L, L, U, U]; xpoints = [x2, fliplr(x2)];
h1 = fill(xpoints, ypoints,colors.black,'linestyle','none'); % 'color','b','EdgeColor','none');
h1.FaceAlpha = 0.3;

jamovi.f2B.BMI_contrasts = BMI.units';

if flag_save == 1
    saveas(gcf, strcat('2A_histogram','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

%% 2B: Contrast-dependent binocular modulation: Scatter
clear de* nde* bin*

method = 'bslc';
tw = info.windows-1;
% Observed responses
switch method
    case 'raw'
        for i = 1:uct
            for c = 1:4
                de_units(c,i)  = SAMPLE.UNIT(i).DE.RESP_avg(c,tw);
                nde_units(c,i) = SAMPLE.UNIT(i).NDE.RESP_avg(c,tw);
                bin_units(c,i) = SAMPLE.UNIT(i).BIN.RESP_avg(c,tw);
            end
        end
    case 'bslc'
        for i = 1:uct
            for c = 1:4
                de_units(c,i)  = SAMPLE.UNIT(i).DE.RESP_avg(c,tw)-SAMPLE.UNIT(i).DE.RESP_avg(c,end);
                nde_units(c,i) = SAMPLE.UNIT(i).NDE.RESP_avg(c,tw)-SAMPLE.UNIT(i).NDE.RESP_avg(c,end);
                bin_units(c,i) = SAMPLE.UNIT(i).BIN.RESP_avg(c,tw)-SAMPLE.UNIT(i).BIN.RESP_avg(c,end);
            end
        end
    case 'percent change'
        for i = 1:uct
            for c = 1:4
                de_units(c,i)  = ((SAMPLE.UNIT(i).DE.RESP_avg(c,tw) - SAMPLE.UNIT(i).DE.RESP_avg(c,end)) / (SAMPLE.UNIT(i).DE.RESP_avg(c,tw))) * 100;
                nde_units(c,i) = ((SAMPLE.UNIT(i).NDE.RESP_avg(c,tw) - SAMPLE.UNIT(i).NDE.RESP_avg(c,end)) / (SAMPLE.UNIT(i).NDE.RESP_avg(c,tw))) * 100;
                bin_units(c,i) = ((SAMPLE.UNIT(i).BIN.RESP_avg(c,tw) - SAMPLE.UNIT(i).BIN.RESP_avg(c,end)) / (SAMPLE.UNIT(i).BIN.RESP_avg(c,tw))) * 100;
            end
        end
end
        
lim = 250;
xlimit = [0 lim];
ylimit = [0 lim];
sz = 10;
x1 = linspace(0,lim);
y1 = linspace(0,lim);

labels = {'Low','Med','High'};

figure('position',[312.3333333333333,767,749.3333333333333,205.3333333333333]);

for c = 1:3
    subplot(1,3,c)
    plot(x1,y1,'linewidth',1.5,'color','k'); hold on
    
    clear x X b b1 yCalc1 yCalc2 
    %plot(x2,y2,'linewidth',2,'color',colors.light_orange,'linestyle','-');
    x = de_units(c+1,~isnan(de_units(c+1,:)))'; y = bin_units(c+1,~isnan(bin_units(c+1,:)))';
    X = [ones(length(x),1) x];
    b = X\y;
    b1 = x\y; slope(c) = b1;
    yCalc1 = b1*x;
    yCalc2 = X*b;
    s = scatter(x,y,sz,colors.black,'filled'); hold on
    set(gca,'linewidth',1.5,'fontsize',12,'xlim',xlimit,'ylim',ylimit,'box','off','tickdir','out');
    yticks([0, 50, 100, 150, 200, 250]); xticks([0, 50, 100, 150, 200, 250]);
    plot(x,yCalc2,'linewidth',2,'color',colors.purple+0.3,'linestyle','-');
    if c > 1
        yticklabels([]);
    end
    xlabel([]);
    s.MarkerFaceAlpha = 0.5;
    title(labels{c});
end

if exist('figDir','var') 
    cd(figDir);
    saveas(gcf, strcat('2B_scatter','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end



%% Figure 2B: Contrast-dependent binocular modulation: Histogram

method = 'raw';
tw = info.windows-1;

% Observed responses
switch method
    case 'raw'
        for i = 1:uct
            for c = 2:4
                de_units(c-1,i)  = SAMPLE.UNIT(i).DE.RESP_avg(c,tw);
                nde_units(c-1,i) = SAMPLE.UNIT(i).NDE.RESP_avg(c,tw);
                bin_units(c-1,i) = SAMPLE.UNIT(i).BIN.RESP_avg(c,tw);
            end
        end
    case 'bslc'
        for i = 1:uct
            for c = 2:4
                de_units(c-1,i)  = SAMPLE.UNIT(i).DE.RESP_avg(c,tw)-SAMPLE.UNIT(i).DE.RESP_avg(c,end);
                nde_units(c-1,i) = SAMPLE.UNIT(i).NDE.RESP_avg(c,tw)-SAMPLE.UNIT(i).NDE.RESP_avg(c,end);
                bin_units(c-1,i) = SAMPLE.UNIT(i).BIN.RESP_avg(c,tw)-SAMPLE.UNIT(i).BIN.RESP_avg(c,end);
            end
        end
    case 'percent change'
        for i = 1:uct
            for c = 2:4
                de_units(c-1,i)  = ((SAMPLE.UNIT(i).DE.RESP_avg(c,tw) - SAMPLE.UNIT(i).DE.RESP_avg(c,end)) / (SAMPLE.UNIT(i).DE.RESP_avg(c,tw))) * 100;
                nde_units(c-1,i) = ((SAMPLE.UNIT(i).NDE.RESP_avg(c,tw) - SAMPLE.UNIT(i).NDE.RESP_avg(c,end)) / (SAMPLE.UNIT(i).NDE.RESP_avg(c,tw))) * 100;
                bin_units(c-1,i) = ((SAMPLE.UNIT(i).BIN.RESP_avg(c,tw) - SAMPLE.UNIT(i).BIN.RESP_avg(c,end)) / (SAMPLE.UNIT(i).BIN.RESP_avg(c,tw))) * 100;
            end
        end
end

close all
clear BMI
for i = 1:uct
    for c = 1:3
       BMI.units(c,i) = (bin_units(c,i)-de_units(c,i)) ./ (bin_units(c,i)+de_units(c,i));
    end
end

BMI.avg = nanmean(BMI.units,2);


p = 95;
for c = 1:3
    BMI.CI(c,:) = CIFcn(BMI.units(c,:),p);
end

figure; 

color = [colors.black+0.5;...
         colors.black+0.2;...
         colors.black-.2];
sz = 20;
nbins = 25;
mean_mark = 60;
     
for c = 1:3
    subplot(1,size(BMI.units,1),c)
    plot(mean_mark,BMI.avg(c,:),'d','linewidth',3,'color',color(c,:)); hold on

    histogram(BMI.units(c,:),nbins,'binWidth',0.025,'orientation','horizontal',...
        'faceAlpha',0.5,'linewidth',1.2,'faceColor',color(c,:),'edgealpha',0.3); hold on
    
    set(gca,'linewidth',1.5,'fontsize',14,'ylim',[-.35 0.35],'xlim',[0 mean_mark],'tickdir','out','box','off');
    hl = hline(0); set(hl,'linestyle','-','linewidth',1,'color','k');
    yticks([-0.35,0,0.35]); xticks([0, mean_mark]);
    % Identify CI
    L = BMI.CI(c,1);
    U = BMI.CI(c,2);
    
    % fill the 95% CI
    x2 = get(gca,'xlim');
    ypoints = [L, L, U, U]; xpoints = [x2, fliplr(x2)];
    h1 = fill(xpoints, ypoints,colors.black,'linestyle','none'); % 'color','b','EdgeColor','none');
    h1.FaceAlpha = 0.3;

    if c ~= 1
        yticklabels([]);
    end
end

jamovi.f2B.BMI_contrasts = BMI.units';

binocular_modulation = BMI.units;

% if flag_save == 1
%     saveas(gcf, strcat('2B_hist','.svg'));
%     disp("Figure saved");
% else
%     disp("Figure not saved");
% end



%% Figure 2B - SwordTips

clear SAMPLE
SAMPLE = UNIT(1:234);

tw = info.windows-1;
clear A B
for u = 1:length(SAMPLE)
    A(:,:,u) = SAMPLE(u).DE.RESP_avg(2:end,1:tw); %-PARAMS.B(tw,u);
    B(:,:,u) = SAMPLE(u).BIN.RESP_avg(2:end,1:tw); %-PARAMS.B(tw,u);
end

close all
clear RFII
for u = 1:size(A,3)
    for c = 1:size(A,1)
        for t = 1:size(A,2)
            RFII(c,t,u) = (B(c,t,u) - A(c,t,u)) ./ (A(c,t,u) + B(c,t,u));
        end
    end
end



RFII_avg = nanmean(RFII,3);
for c = 1:3
    RFII_CI(c,:) = CIFcn(RFII(c,tw,:),p);
end

color = [colors.black+0.5;...
         colors.black+0.2;...
         colors.black-.2];
sz = 20;
nbins = 25;
     
for c = 1:3
    figure('position',[178,419,202,151]);
    plot(RFII_avg(c,tw),50,'d','linewidth',1.5,'color',color(c,:)); hold on
    histogram(RFII(c,tw,:),nbins,'binWidth',0.025,'orientation','vertical',...
        'faceAlpha',0.5,'linewidth',1,'faceColor',color(c,:),'edgealpha',0.1); hold on
    
    set(gca,'linewidth',1.5,'fontsize',14,'xlim',[-.3 0.3],'ylim',[0 50],'tickdir','out','box','off');
    yticks([0,25,50]); xticks([-0.3,0, 0.3]);
    hl = vline(0); set(hl,'linestyle','-','linewidth',1,'color','k');
    % Identify CI
    L = RFII_CI(c,1);
    U = RFII_CI(c,2);
    
    % fill the 95% CI
    x2 = get(gca,'ylim');
    xpoints = [L, L, U, U]; ypoints = [x2, fliplr(x2)];
    h1 = fill(xpoints, ypoints,colors.black,'linestyle','none'); % 'color','b','EdgeColor','none');
    h1.FaceAlpha = 0.3;

    if c ~= 1
        yticklabels([]);
    end
    
    if flag_save == 1
        saveas(gcf, strcat('2B_sword_',string(c),'.svg'));
        disp("Figure saved");
    else
        disp("Figure not saved");
    end

end




%% Figure 3. Binocular facilitation dynamically evolves in a contrast-dependent manner
% A. Mean spike density functions (SDFs) 
% B. Normalized delta spiking
% C. Latency to Zero Crossing Analysis

%% Figure 3A. Spike-density functions

% data
clear de nde bin
for u = 1:uct
    de.units(:,:,u) = SAMPLE.UNIT(u).DE.SDF_avg(2:4,:)-SAMPLE.UNIT(u).DE.RESP_avg(2:end,end); %-PARAMS.B(tw,u);
    bin.units(:,:,u) = SAMPLE.UNIT(u).BIN.SDF_avg(2:4,:)-SAMPLE.UNIT(u).BIN.RESP_avg(2:end,end); %-PARAMS.B(tw,u);
end

%de.units(de.units < 0) = 0; bin.units(bin.units < 0) = 0; % rectify

de.avg = nanmean(de.units,3);
de.err = nanstd(de.units,[],3)./sqrt(size(de.units,3));
bin.avg = nanmean(bin.units,3);
bin.err = nanstd(bin.units,[],3)./sqrt(size(bin.units,3));

p = 95;
clear CI
for t = 1:size(de.units,2)
    for c = 1:3
        CI{c}.DE(:,t) = CIFcn(squeeze(de.units(c,t,:)),p);
        CI{c}.BIN(:,t) = CIFcn(bin.units(c,t,:),p);
    end
end

% Plot
figure('position',[265,717,826,178.6666666666666]);
ylimit = [-5 175];
for c = 1:3
    
    subplot(1,3,c)
    plot(info.sdfWin,de.avg(c,:),'color',colors.blue,'linewidth',2); hold on
    ci1 = ciplot(CI{c}.DE(1,:),CI{c}.DE(2,:),info.sdfWin,colors.blue,0.3); set(ci1,'linestyle','none','handleVisibility','off');
    
    plot(info.sdfWin,bin.avg(c,:),'color',colors.dark_red,'linewidth',2); hold on
    ci1 = ciplot(CI{c}.BIN(1,:),CI{c}.BIN(2,:),info.sdfWin,colors.dark_red,0.3); set(ci1,'linestyle','none','handleVisibility','off');
    
%     clear h p stats
%     [h,~, ~] = ttest_time(bin.units,de.units); tmh = find(h(c,:));
%     scatter(info.sdfWin(tmh), ones(1,numel(tmh)) * max(bin.avg(c,:)) * 1.1,0.1,'.k')
%     
    
    %lgnd = legend('Monocular','Binocular','location','northeast'); legend boxoff; set(lgnd,'fontsize',16);
    if c == 1
        ylabel('spiking activity (spikes / s)','fontsize',12);
    end
    
    % if c == 2
    %     xlabel('time (s) from stimulus onset','fontsize',14);
    % end
    
    if c > 1
        yticklabels([]);
            ylabel([]);  xlabel([]); xticklabels([]);
    end
    

    
    set(gca,'Box','off','TickDir','out','linewidth',1.5,'ylim',ylimit,'xlim',[-0.05 .250],...
        'Fontsize',12)
    
    vl = vline(0);
    set(vl,'linewidth',1,'linestyle',':');
end

if flag_save == 1
    saveas(gcf, strcat('5A_sdfs','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end


%% Residual: Pre-processing

% data

clear de nde bin nDE nBIN RES res
for u = 1:uct
    de.units(:,:,u) = SAMPLE.UNIT(u).DE.SDF_avg(2:4,:)-SAMPLE.UNIT(u).DE.RESP_avg(2:end,end); %-PARAMS.B(tw,u);
    bin.units(:,:,u) = SAMPLE.UNIT(u).BIN.SDF_avg(2:4,:)-SAMPLE.UNIT(u).BIN.RESP_avg(2:end,end); %-PARAMS.B(tw,u);
end

x = info.sdfWin;
de.units(de.units < 0) = 0; bin.units(bin.units < 0) = 0; % rectify

for u = 1:uct

    % determine maximum firing to binocular stimulation for each unit
    mx = max([max(bin.units(:,:,u),[],'all'),max(de.units(:,:,u),[],'all')]);
    %mx = max(bin.units(:,:,u),[],'all');

    % for each contrast, normalize binocular modulation 
    for c = 1:3
        BIN = squeeze(bin.units(c,:,u)); 
        DE = squeeze(de.units(c,:,u)); 

        y = (BIN - DE) ./ mx;
        y2 = y; y2(y2 < 0) = 0; % rectify only for area under the curve

        % area under the curve       
        % organize in matrix A, as in question
        A = [x' y2'];

        % find area between 0.05 and 0.250
        x1 = 0.05;
        x2 = 0.250;
        % find indices of these points in A
        idx1 = find(A(:,1) >= x1, 1);
        idx2 = find(A(:,1) >= x2, 1);

        % get area under curve from x1 to x2 using trapz
        a1 = trapz(A(idx1:idx2,1), A(idx1:idx2,2));

        % organize
        RES{c}.residual(:,u) = y; % un-rectified 
        RES{c}.residual_rectified(:,u) = y2; % rectified
        RES{c}.auc(:,u) = a1;
        RES{c}.magnitude(:,u) = max(y,[],'all');
        res.raw(c,:,u) = y;
        res.rectified(c,:,u) = y2;
    end
   
    
end

% descriptive statistics

for c = 1:3
    for t = 1:size(de.units,2)
        p = 95;
        RES{c}.CI(:,t) = CIFcn(RES{c}.residual(t,:),p);
    end
    
    RES{c}.avg = nanmean(RES{c}.residual,2); % swapped this to median
    RES{c}.upper = RES{c}.CI(2,:);
    RES{c}.lower = RES{c}.CI(1,:);
    res.units(c,:,:) = RES{c}.residual;
end


res.zeros = zeros(size(res.units,1),size(res.units,2),size(res.units,3));

%% 3B. Residual Spiking. 
ylimit2 = [-.1 .25];

figure('position',[265,739,686.6666666666666,156.6666666666666]);

for c = 1:3
    subplot(1,3,c)
    ylabel('Norm. Residual'); xlabel('Contrast');
    set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
        'xscale','linear','xlim',[-0.05,0.250],...
        'ylim',ylimit2,'linewidth',1.5,'fontsize',12);
    hl = hline(0); set(hl,'linewidth',1,'linestyle','-','color','k','HandleVisibility','off'); hold on
    pl = plot(info.sdfWin, RES{c}.avg,'linestyle','-','linewidth',1,'color',colors.purple); hold on
    ci = ciplot(RES{c}.lower,RES{c}.upper,info.sdfWin,colors.purple,0.4);
    set(ci,'linestyle',':','linewidth',0.25,'handleVisibility','off'); hold on
    area(RES{c}.lower,0);
    
    clear h p stats
    [h,~, ~] = ttest_time(res.units,res.zeros); tmh = find(h(c,:));
    scatter(info.sdfWin(tmh), ones(1,numel(tmh)) * max(RES{c}.avg) * 1.2,0.1,'.k')
    
    if c == 1
        ylabel('Delta spiking');
        yticklabels('auto'); xticklabels('auto'); xlabel('time (s)');
    else
        ylabel([]); xlabel([]); yticklabels([]);
    end
    grid on
        
end

h(:,1:150) = 0;
clear latency
for c =1:3
    onset = find(h(c,:),1,'first');
    offset = find(h(c,:),1,'last');
    duration = offset-onset;
    
    latency.onset(c) = info.sdfWin(onset);
    latency.offset(c) = info.sdfWin(offset);
    latency.duration(c) = duration;

end


if flag_save == 1
    saveas(gcf, strcat('5A_residual','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end


%% 3C. Top - Peak magnitude of facilitation
clear magnitude*
for c =1:3
    magnitudes(c,:) = RES{c}.magnitude;
    magnitude_avg(c) = nanmean(RES{c}.magnitude);
    magnitude_med(c) = nanmedian(RES{c}.magnitude);
    auc(c) = nanmean(RES{c}.auc);
end

values = reshape(magnitudes,[],1);
boxcolors = [colors.black;colors.black+0.3;colors.black+0.5];

clear G g
for i = 1:3
    g(i,:) = repmat(i,size(magnitudes,2),1);
end
G = reshape(g,[],1);

figure('position',[1000,1000.333333333333,336.3333333333333,317.3333333333333]);

% boxplot
sz = 2;
scatter(G,values,sz,colors.black+0.1); hold on

boxplot(values,G,'notch','on','whisker',1,'symbol','','plotstyle',...
    'traditional','BoxStyle','outline','MedianStyle','line','Widths',0.5,...
    'FactorDirection','list','colors',colors.black+0.5); hold on

set(gca,'fontsize',12','linewidth',1.5,'box','off','ylim',[0, 1]);
xticklabels({'Low','Med','High'}); xlabel('Contrast');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',0.4,'linewidth',0.1);
end
plot(1:3,magnitude_med,'ro-','linewidth',1);



jamovi.magnitudes = magnitudes';

if flag_save == 1
    saveas(gcf, strcat('5A_magnitude','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

%% (not currently reported). Area of facilitation
clear auc* values
for c =1:3
    auc(c,:) = RES{c}.auc;
    auc_mean(c) = nanmean(RES{c}.auc);
end

values = reshape(auc,[],1);
boxcolors = [colors.black;colors.black+0.3;colors.black+0.5];

clear G g
for i = 1:3
    g(i,:) = repmat(i,size(auc,2),1);
end
G = reshape(g,[],1);
close all
figure('position',[1000,1017,355,321]);

% boxplot
sz = 2;
scatter(G,values,sz,colors.black); hold on
boxplot(values,G,'notch','on','whisker',1,'symbol','','plotstyle',...
    'traditional','BoxStyle','outline','MedianStyle','target','Widths',0.5,...
    'FactorDirection','list'); hold on

xticklabels({'Low','Med','High'}); xlabel('Contrast');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',.6,'linewidth',1.2);
end

%% 3C. Bottom - Time to Zero Crossing, Pre-process

clear de* nde* bin* diff toZero zci ZEROX
sample = 314;
window = 1:401;

clear XLAT
for u = 1:sample
    for c = 1:3
        data = res.raw(c,:,u);
        if isempty(data)
            break
        else
            
            % prepare to fit the delta response with a smoothing function
            try
                [xData, yData] = prepareCurveData( window, data );
            catch
            end
            
            
            % sum of sine smoothing function
            ft = fittype( 'sin8' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Display = 'Off';
            opts.Lower = [-Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf -Inf 0 -Inf];
            opts.Normalize = 'on';
            try
                [fitresult, gof] = fit( xData, yData, ft, opts);
                XLAT.fitted(c,:,u) = feval(fitresult,window);
            catch
                warning('Could not fit'); % for missing observations
                fitresult = nan;
                XLAT.fitted(c,:,u) = nan;
            end

        end
    end
end

%% Identify Zero Crossing
% previous
sample = 314; 
clear ZEROX

for u = 1:sample
    for c = 1:3
        win = latency.onset(c); % timepoint of initial, sign. facilitation 
        window = find(info.sdfWin == win); % convert the timepoint
        
        % threshold to catch zero crossing
        thresh = nanmedian(RES{c}.CI(2,201:401)-RES{c}.CI(1,201:401)); % defined as median width of the error in our observation
    
        % data: the delta spiking (BIN - MON) that has been smoothed
        data = XLAT.fitted(c,:,u);
    
        if any(isnan(data)) % for missing observations
            ZEROX.offset(c,u) = nan;
            ZEROX.duration(c,u) = nan;
        else
            crossing_idx = find(data < -thresh); % find all timepoints where threshold is exceeded
            crossing_idx2 = crossing_idx(crossing_idx >= window); % narrow the search to just after sign. facil. 
            offset_idx = crossing_idx2(find(crossing_idx2,1,'first')); % identify first instance of zero crossing (offset of facil.)
            
            if ~isempty(offset_idx)
                offset = info.sdfWin(offset_idx);
            else
                offset = info.sdfWin(end); % if offset never occurs after facil. onset
            end
            
            ZEROX.offset(c,u) = offset;
            ZEROX.duration(c,u) = offset - (latency.onset(c));
            if ZEROX.duration(c,u) == 0 
                figure;
                plot(info.sdfWin, data); hold on
                plot(info.sdfWin, res.raw(c,:,u));
            end
        end
        

    end
end

% some descriptives
duration_mean = nanmean(ZEROX.duration,2); % report in text
duration_med = nanmedian(ZEROX.duration,2); % redline in 3C 

p = 95;
clear ZEROX_CI
for c = 1:3
    offset_CI(:,c) = CIFcn(squeeze(ZEROX.offset(c,:)),p);
    duration_CI(:,c) = CIFcn(squeeze(ZEROX.duration(c,:)),p);
end

facil_duration = ZEROX.duration;



%% 3C - Bottom. Duration of facilitation. 
values = reshape(ZEROX.duration,[],1);
boxcolors = [colors.black;colors.black+0.3;colors.black+0.5];

clear G g
for i = 1:3
    g(i,:) = repmat(i,size(ZEROX.duration,2),1);
end
G = reshape(g,[],1);
close all
figure('position',[1000,1000.333333333333,336.3333333333333,317.3333333333333]);

% boxplot
sz = 2;
scatter(G,values,sz,colors.black+0.1); hold on

boxplot(values,G,'notch','on','whisker',1,'symbol','','plotstyle',...
    'traditional','BoxStyle','outline','MedianStyle','target','Widths',0.5,...
    'FactorDirection','list','colors',colors.black+0.5); hold on
set(gca,'fontsize',12','linewidth',1.5,'box','off','ylim',[0, 0.25]);
xticklabels({'Low','Med','High'}); xlabel('Contrast');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
    patch(get(h(j),'XData'),get(h(j),'YData'),boxcolors(j,:),'FaceAlpha',0.4,'linewidth',0.1);
end
plot(1:3,duration_med,'ro-','linewidth',1);

jamovi.durations = ZEROX.duration';
facil_duration = ZEROX.duration;


if flag_save == 1
    saveas(gcf, strcat('5A_duration','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end


%% Figure 4. Temporal evolution of CRFs 
% A. Naka-Rushton function 
% B. 16 CRF pairs with 2D delta surface
% C. 3 CRF pairs isolated for closer inspection

%% Figure 4. Full Resolution Time X Contrast

addpath 'C:\Users\bmitc\OneDrive - Vanderbilt\Maier Lab\Analysis\Figures\dios'
clear y z x xMat yMat zMat
x = [0.01:0.01:1]';
windows = info.windows-2; % this is all but bsl and full
window_idx = 1:info.windows-2;
xMat = repmat(x, 1, windows);
y = 1:windows;
yMat = repmat(y, numel(x), 1); %//For plot3

zMat = squeeze(CRF2.DE.curve_avg(:,window_idx));
zMat2 = squeeze(CRF2.BIN.curve_avg(:,window_idx));
zMat3 = squeeze(CRF2.BIN.curve_avg(:,window_idx)-CRF2.DE.curve_avg(:,window_idx));

% zMat = dom. eye (contrast x time window)
% zMat2 =  binocular (contrast x time window)
% zMat3 = residual (contrast x time window)
% Plot

display = [5,10,15];
figure('position',[759,747,590,352.6666666666665]);
%figure('position',[582.3333333333333,766.3333333333333,766.6666666666667,333.3333333333333]);
plot3(xMat, yMat, zMat,'color',[colors.blue],'linestyle','-','linewidth',1.5,'HandleVisibility','off');  hold on
plot3(xMat, yMat, zMat2,'color',[colors.dark_red],'linestyle','-','linewidth',1.5,'HandleVisibility','off');  hold on
plot3(xMat(:,display), yMat(:,display), zMat(:,display),'color',colors.blue,'linestyle','-','linewidth',1.5,'HandleVisibility','off');  hold on
plot3(xMat(:,display), yMat(:,display), zMat2(:,display),'color',colors.dark_red,'linestyle','-','linewidth',1.5,'HandleVisibility','off');  hold on
ylabel('Time window'); xlabel('Contrast'); zlabel('Delta Spiking');
xticks([0.10,0.22,0.40,0.90]); xticklabels({'0.1','0.22','0.45','0.90'});
set(gca,'fontsize',10);
%view([60.7583711500765 15.8693628136953]); %// Adjust viewing angle so you can clearly see data
grid on

% if flag_save == 1
%     export_fig 7_3D_CRFs -pdf
%     disp("Figure saved");
% else
%     disp("Figure not saved");
% end


if ~exist('mycbar','var')
    load('colormap_rbw','mycbar');
end

% Normalized Residual
x = [0.01:0.01:1]';
windows = info.windows-2; % this is all but bsl and full
xMat = repmat(x, 1, windows);
y = 1:windows;
yMat = repmat(y, numel(x), 1); %//For plot3
window_idx = 1:16;

clear BIN DE RES

for w = 1:length(window_idx)
    tw = window_idx(w);
    for u = 1:length(PARAMS2.B)
        BIN(:,w,u) = feval(PARAMS2.BIN.fitobject(tw,u).f,x);
        DE(:,w,u) = feval(PARAMS2.DE.fitobject(tw,u).f,x);
    end
    
    clear nDE nBIN
    
    for u = 1:length(PARAMS2.B)
        nDE(:,w,u) = DE(:,w,u) ./ max(BIN(:,w,u));
        nBIN(:,w,u) = BIN(:,w,u) ./ max(BIN(:,w,u));
    end
    
    for u = 1:length(PARAMS2.B)
        RES.residual(:,w,u) = nBIN(:,w,u)' - nDE(:,w,u)';
    end
    
    clear CI
    for c = 1:length(x)
        p = 95;
        CI(:,w,c) = CIFcn(squeeze(RES.residual(c,w,:)),p);
    end
    
    RES.avg = nanmean(RES.residual,3); % swapped this to median
    RES.upper = CI(2,:,:);
    RES.lower = CI(1,:,:);
end

%figure; 
mesh(xMat, yMat, RES.avg(:,:)); hold on
colors_pgw = [0.5 0.5 0.5;0.533333361148834 0.533333361148834 0.533333361148834;0.566666662693024 0.566666662693024 0.566666662693024;0.600000023841858 0.600000023841858 0.600000023841858;0.633333325386047 0.633333325386047 0.633333325386047;0.666666686534882 0.666666686534882 0.666666686534882;0.699999988079071 0.699999988079071 0.699999988079071;0.733333349227905 0.733333349227905 0.733333349227905;0.766666650772095 0.766666650772095 0.766666650772095;0.800000011920929 0.800000011920929 0.800000011920929;0.833333313465118 0.833333313465118 0.833333313465118;0.866666674613953 0.866666674613953 0.866666674613953;0.899999976158142 0.899999976158142 0.899999976158142;0.933333337306976 0.933333337306976 0.933333337306976;0.966666638851166 0.966666638851166 0.966666638851166;1 1 1;0.986266672611237 0.96560001373291 0.990400016307831;0.972533345222473 0.93120002746582 0.980799973011017;0.95880001783371 0.896799981594086 0.971199989318848;0.945066690444946 0.862399995326996 0.961600005626678;0.931333363056183 0.828000009059906 0.952000021934509;0.917600035667419 0.793600022792816 0.942399978637695;0.903866708278656 0.759199976921082 0.932799994945526;0.890133321285248 0.724799990653992 0.923200011253357;0.876399993896484 0.690400004386902 0.913600027561188;0.862666666507721 0.656000018119812 0.903999984264374;0.848933339118958 0.621599972248077 0.894400000572205;0.835200011730194 0.587199985980988 0.884800016880035;0.821466684341431 0.552799999713898 0.875200033187866;0.807733356952667 0.518400013446808 0.865599989891052;0.794000029563904 0.483999997377396 0.856000006198883];
colormap(colors_pgw);
colorbar
caxis([-0.1 0.1]);
ylabel('Time window'); xlabel('Contrast'); zlabel('Spiking activity (imps. / sec)');
xticks([0.10,0.22,0.40,0.90]); xticklabels({'0.1','0.22','0.45','0.90'});
yticks([1.1,5,10,15]); yticklabels({'1','5','10','15'});
view([75.5224592198716 25.8215997485344]); %// Adjust viewing angle so you can clearly see data
set(gca,'XScale','log','xlim',[0.05 1.05],'ylim',[1.5,16])


if flag_save == 1
    export_fig 7_CRFs -transparent -pdf
    disp("Figure saved");
else
    disp("Figure not saved");
end

%% 4Ca. Early CRF, Window #5

figure('position',[336,526,414,334]);
x = 0.01:.01:1;
ylimit = [0, 140];
tw = 5;

% DE
plot(x, SAMPLE.CRF.DE.curve_avg(:,tw),'color',colors.blue,'linestyle','-','linewidth',2); hold on
ci1 = ciplot(CRF2.DE.curve_lower(:,tw),SAMPLE.CRF.DE.curve_upper(:,tw),x,colors.blue,0.15); 
set(ci1,'linestyle','none','edgecolor',colors.blue,'linewidth',1,'handleVisibility','on'); hold on

% BIN
plot(x, SAMPLE.CRF.BIN.curve_avg(:,tw),'color',colors.dark_red,'linestyle','-','linewidth',2); hold on
ci3 = ciplot(SAMPLE.CRF.BIN.curve_lower(:,tw),SAMPLE.CRF.BIN.curve_upper(:,tw),x,colors.dark_red,0.15); 
set(ci3,'linestyle','none','edgecolor',colors.dark_red,'linewidth',1,'handleVisibility','on'); hold on

ylabel('Spiking activity'); 
set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
  'xscale','log','xlim',[0.05,1],...
  'ylim',ylimit,'linewidth',1.5,'fontsize',12); 
xlabel([]); xticklabels([]); xticks([0.05,0.1,0.22,0.45,1]);
%lgn = legend('Monocular','MON 95% CI','Binocular','BIN 95% CI','location','northwest');
%set(lgn,'fontsize',10); legend boxoff


if flag_save == 1
    saveas(gcf, strcat('6_early_CRF','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

% Early delta
clear BIN DE


clear BIN DE
x = 0.01:.01:1;
for u = 1:length(SAMPLE.PARAMS.B)
    BIN(:,u) = feval(SAMPLE.PARAMS.BIN.fitobject(tw,u).f,x);
    DE(:,u) = feval(SAMPLE.PARAMS.DE.fitobject(tw,u).f,x);
end


clear nDE nBIN

for u = 1:length(SAMPLE.PARAMS.B)
    nDE(:,u) = DE(:,u) ./ max(BIN(:,u));
    nBIN(:,u) = BIN(:,u) ./ max(BIN(:,u));
end

clear RES*
for u = 1:length(SAMPLE.PARAMS.B)
    RES.residual(:,u) = nBIN(:,u)' - nDE(:,u)';
end

clear CI
for c = 1:length(x)
    p = 95;
    CI(:,c) = CIFcn(RES.residual(c,:),p);
end

RES.avg = mean(RES.residual,2); % swapped this to median
RES.upper = CI(2,:);
RES.lower = CI(1,:);


ylimit2 = [-.15 .15];

figure('position',[393.6666666666666,838.3333333333333,414,152]);

ylabel('Norm. Residual'); xlabel('Contrast');
set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
  'xscale','log','xlim',[.05,1],...
  'ylim',ylimit2,'linewidth',1.5,'fontsize',12); 
hl = hline(0); set(hl,'linewidth',1,'linestyle','-','color','k','HandleVisibility','off'); hold on
pl = plot(x, RES.avg,'linestyle','-','linewidth',1,'color',colors.purple); hold on
ci = ciplot(RES.lower,RES.upper,x,colors.purple,0.2);
set(ci,'linestyle',':','linewidth',0.25,'handleVisibility','off'); hold on
xticklabels([0.05,0.1,0.22,0.45,1]); xticks([0.05,0.1,0.22,0.45,1]);


if flag_save == 1
    saveas(gcf, strcat('6_earlydelta','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

%% 4Cb. Intermediate CRF, Window #10

figure('position',[336,526,414,334]);
x = 0.01:.01:1;
ylimit = [0, 150];
tw = 10;

% DE
plot(x, SAMPLE.CRF.DE.curve_avg(:,tw),'color',colors.blue,'linestyle','-','linewidth',2); hold on
ci1 = ciplot(CRF.DE.curve_lower(:,tw),SAMPLE.CRF.DE.curve_upper(:,tw),x,colors.blue,0.15); 
set(ci1,'linestyle','none','edgecolor',colors.blue,'linewidth',1,'handleVisibility','on'); hold on

% BIN
plot(x, SAMPLE.CRF.BIN.curve_avg(:,tw),'color',colors.dark_red,'linestyle','-','linewidth',2); hold on
ci3 = ciplot(SAMPLE.CRF.BIN.curve_lower(:,tw),SAMPLE.CRF.BIN.curve_upper(:,tw),x,colors.dark_red,0.15); 
set(ci3,'linestyle','none','edgecolor',colors.dark_red,'linewidth',1,'handleVisibility','on'); hold on

ylabel('Spiking activity'); 
set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
  'xscale','log','xlim',[0.05,1],...
  'ylim',ylimit,'linewidth',1.5,'fontsize',12); 
xlabel([]); xticklabels([]); xticks([0.05,0.1,0.22,0.45,1]);

if flag_save == 1
    saveas(gcf, strcat('6_intermediate_CRF','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

% Intermediate delta plot
clear BIN DE
tw = 11;


clear BIN DE
x = 0.01:.01:1;
for u = 1:length(SAMPLE.PARAMS.B)
    BIN(:,u) = feval(SAMPLE.PARAMS.BIN.fitobject(tw,u).f,x);
    DE(:,u) = feval(SAMPLE.PARAMS.DE.fitobject(tw,u).f,x);
end


clear nDE nBIN

for u = 1:length(SAMPLE.PARAMS.B)
    nDE(:,u) = DE(:,u) ./ max(BIN(:,u));
    nBIN(:,u) = BIN(:,u) ./ max(BIN(:,u));
end

clear RES*
for u = 1:length(SAMPLE.PARAMS.B)
    RES.residual(:,u) = nBIN(:,u)' - nDE(:,u)';
end

clear CI
for c = 1:length(x)
    p = 95;
    CI(:,c) = CIFcn(RES.residual(c,:),p);
end

RES.avg = mean(RES.residual,2); % swapped this to median
RES.upper = CI(2,:);
RES.lower = CI(1,:);


ylimit2 = [-.15 .15];

figure('position',[393.6666666666666,838.3333333333333,414,152]);

ylabel('Norm. Residual'); xlabel('Contrast');
set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
  'xscale','log','xlim',[.05,1],...
  'ylim',ylimit2,'linewidth',1.5,'fontsize',12); 
hl = hline(0); set(hl,'linewidth',1,'linestyle','-','color','k','HandleVisibility','off'); hold on
pl = plot(x, RES.avg,'linestyle','-','linewidth',1,'color',colors.purple); hold on
ci = ciplot(RES.lower,RES.upper,x,colors.purple,0.2);
set(ci,'linestyle',':','linewidth',0.25,'handleVisibility','off'); hold on
xticklabels([0.05,0.1,0.22,0.45,1]); xticks([0.05,0.1,0.22,0.45,1]);


if flag_save == 1
    saveas(gcf, strcat('6_inter_delta','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

%% 4Cc. Late CRF, Window #15

figure('position',[336,526,414,334]);
x = 0.01:.01:1;
ylimit = [0, 150];
tw = info.windows-2;

% DE
plot(x, SAMPLE.CRF.DE.curve_avg(:,tw),'color',colors.blue,'linestyle','-','linewidth',2); hold on
ci1 = ciplot(SAMPLE.CRF.DE.curve_lower(:,tw),SAMPLE.CRF.DE.curve_upper(:,tw),x,colors.blue,0.15); 
set(ci1,'linestyle','none','edgecolor',colors.blue,'linewidth',1,'handleVisibility','on'); hold on

% BIN
plot(x, SAMPLE.CRF.BIN.curve_avg(:,tw),'color',colors.dark_red,'linestyle','-','linewidth',2); hold on
ci3 = ciplot(SAMPLE.CRF.BIN.curve_lower(:,tw),SAMPLE.CRF.BIN.curve_upper(:,tw),x,colors.dark_red,0.15); 
set(ci3,'linestyle','none','edgecolor',colors.dark_red,'linewidth',1,'handleVisibility','on'); hold on

ylabel('Spiking activity'); 
set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
  'xscale','log','xlim',[0.05,1],...
  'ylim',ylimit,'linewidth',1.5,'fontsize',12); 
xlabel([]); xticklabels([]); xticks([0.05,0.1,0.22,0.45,1]);
lgn = legend('Monocular','MON 95% CI','Binocular','BIN 95% CI','location','northwest');
set(lgn,'fontsize',10); legend boxoff

if flag_save == 1
    saveas(gcf, strcat('6_sustained_CRF','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

clear BIN DE
x = 0.01:.01:1;
for u = 1:length(SAMPLE.PARAMS.B)
    BIN(:,u) = feval(SAMPLE.PARAMS.BIN.fitobject(tw,u).f,x);
    DE(:,u) = feval(SAMPLE.PARAMS.DE.fitobject(tw,u).f,x);
end

clear nDE nBIN

for u = 1:length(SAMPLE.PARAMS.B)
    nDE(:,u) = DE(:,u) ./ max(BIN(:,u));
    nBIN(:,u) = BIN(:,u) ./ max(BIN(:,u));
end

clear RES*
for u = 1:length(SAMPLE.PARAMS.B)
    RES.residual(:,u) = nBIN(:,u)' - nDE(:,u)';
end

clear CI
for c = 1:length(x)
    p = 95;
    CI(:,c) = CIFcn(RES.residual(c,:),p);
end

RES.avg = mean(RES.residual,2); % swapped this to median
RES.upper = CI(2,:);
RES.lower = CI(1,:);


ylimit2 = [-.15 .15];

figure('position',[393.6666666666666,838.3333333333333,414,152]);

ylabel('D Spiking (norm.)'); xlabel('Contrast');
set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
  'xscale','log','xlim',[.05,1],...
  'ylim',ylimit2,'linewidth',1.5,'fontsize',12); 
hl = hline(0); set(hl,'linewidth',1,'linestyle','-','color','k','HandleVisibility','off'); hold on
plot(x, RES.avg,'linestyle','-','linewidth',1,'color',colors.purple); hold on
ci = ciplot(RES.lower,RES.upper,x,colors.purple,0.2);
set(ci,'linestyle',':','linewidth',0.25,'handleVisibility','off'); hold on
xticklabels([0.05,0.1,0.22,0.45,1]); xticks([0.05,0.1,0.22,0.45,1]);

if flag_save == 1
    saveas(gcf, strcat('6_late_delta','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

%% Figure 5. Models of binocular 'gain' and 'gain-control' 
% A. Simulated models for illustration
% B. Top - Fitting response-gain and contrast-gain for 3 windows of time 
% B. Bottom - Model goodness of fit for all 16 windows of time. 
% C. Mean parameters (Rmax and C50) as a function of time. 

%% 5A. Effects of manipulating response (response-gain) and contrast (contrast-gain)

rMax = 1; %multiplicative response gain factor (=highest response amplitude)
K = 20;  %normalization pool (determines x position)
n = 2;   %exponent that determines rise and saturation
b = 0;    %baseline offset w/o stim

clear Rc
clear mRc
for c=1:100 
    %compute response for each conrast level c:  
    R(c) = rMax*(c^n/(c^(n) + K^(n)))+b; 
    
    % response-gain control prediction
    Rgc = 0.80;
    RGC(c) = Rgc*(c^n /(c^(n) + K^(n)))+b;
    
    % response-gain prediction
    Rg = 1.2;
    RG(c) = Rg*(c^n /(c^(n) + K^(n)))+b;
    
    % contrast-gain control prediction
    Cgc = 40;
    CGC(c) = rMax*(c^n /(c^(n) + Cgc^(n)))+b;
    
    % contrast-gain prediction (sensitivity)
    Cg = 10;
    CG(c) = rMax*(c^n /(c^(n) + Cg^(n)))+b;
    
    % additive gain
    ab = 1;
    aRc(c) = rMax*(c^n /(c^(n) + K^(n)))+ab;
    
    % n
    N = 1;
    nGC(c) = rMax*(c^N /(c^(N) + K^(N)))+b;
    
    % n
    N = 3.5;
    nG(c) = rMax*(c^N /(c^(N) + K^(N)))+b;
end

close all;
figure('position',[538.3333333333333,958.3333333333333,678,265.3333333333333]);
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
semilogx(1:100,R,'k','linewidth',2); hold on
semilogx([1:100],RGC,'k','linestyle','--')
semilogx([1:100],RG,'k','linestyle','--')
ylim([0 1.2])
xlabel([]);
xticks([1,10,100]);
set(gca,'box','off','tickdir','out','fontsize',14,'linewidth',1.5);
xticklabels([]);
ylabel('Response (a.u.)'); 


nexttile
semilogx(1:100,R,'k','linewidth',2); hold on
semilogx([1:100],CGC,'k','linestyle','--')
semilogx([1:100],CG,'k','linestyle','--')
set(gca,'ylim',[0 1.2]);
xlabel([])
xticks([1,10,100]);
set(gca,'box','off','tickdir','out','fontsize',14,'linewidth',1.5);
xticklabels([]);
ylabel([]); yticklabels([]);


% if flag_save == 1
%     saveas(gcf, strcat('5A_sim','.svg'));
%     disp("Figure saved");
% else
%     disp("Figure not saved");
% end

%% 5B - Top. Model fit: Response-gain set vs. Contrast-gain set

figure('position',[1065,959,550,155.3333333333333]); 
x = 0.01:.01:1;
ylimit = [0, 140];
windows = [5,10,15];


for t = 1:length(windows)
    % DE
    tw = windows(t);
    subplot(1,3,t)
    plot(x, CRF.DE.curve_avg(:,tw),'color',colors.black,'linestyle','-','linewidth',2); hold on
    % ci1 = ciplot(CRF.DE.curve_lower(:,tw),CRF.DE.curve_upper(:,tw),x,colors.blue,0.15);
    % set(ci1,'linestyle','none','edgecolor',colors.blue,'linewidth',1,'handleVisibility','on'); hold on
    
    % BIN
    %plot(x, CRF.BIN.curve_avg(:,tw),'color',colors.dark_red,'linestyle','-','linewidth',2); hold on
    ci3 = ciplot(CRF.BIN.curve_lower(:,tw),CRF.BIN.curve_upper(:,tw),x,colors.dark_red,0.5);
    set(ci3,'linestyle','none','edgecolor',colors.dark_red,'linewidth',1,'handleVisibility','on'); hold on
    
    % % Response-gain
    plot(x, GAIN.response.curve_avg(:,tw),'color',colors.black,'linestyle','--','linewidth',1); hold on
    %ci3 = ciplot(GAIN.response.curve_lower(:,tw),GAIN.response.curve_upper(:,tw),x,colors.black,0.2);
    % %set(ci3,'linestyle',':','edgecolor',colors.black,'linewidth',1,'handleVisibility','on'); hold on
    %
    % % Contrast-gain
    plot(x, GAIN.contrast.curve_avg(:,tw),'color',colors.black,'linestyle',':','linewidth',1); hold on
    % %ci3 = ciplot(GAIN.response.curve_lower(:,tw),GAIN.response.curve_upper(:,tw),x,colors.black,0.2);
    % %set(ci3,'linestyle',':','edgecolor',colors.black,'linewidth',1,'handleVisibility','on'); hold on
    
    plot(x, GAIN.additive.curve_avg(:,tw),'color',colors.black,'linestyle',':','linewidth',1); hold on
    
    
    set(gca,'box','off','tickdir','out','Xgrid','off','ygrid','off',...
        'xscale','log','xlim',[0.05,1],...
        'ylim',ylimit,'linewidth',1.5,'fontsize',12);
    xlabel([]); xticklabels([]); xticks([0.05,0.1,0.22,0.45,1]); yticklabels([]); 
    if t == 1
        ylabel('Spiking activity');
        yticklabels('auto');
    end
    %lgn = legend('Monocular','MON 95% CI','Binocular','BIN 95% CI','location','northwest');
    %set(lgn,'fontsize',10);
    box on
end


if flag_save == 1
    saveas(gcf, strcat('5a_models','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end

%% 5B - Bottom. Average Model-fit across time

x = 1:16;
y1 = mean(GAIN.response.rsquare,2);
y1_CI = GAIN.response.rsquare_CI;
y2 = mean(GAIN.contrast.rsquare,2);
y2_CI = GAIN.contrast.rsquare_CI;
y3 = mean(GAIN.additive.rsquare,2);
y3_CI = GAIN.additive.rsquare_CI;
figure('position',[1065,959,550,155.3333333333333]); 
plot(x,y1(x),'-o','linewidth',2,'color',colors.orange,'linewidth',2); hold on
plot(x,y2(x),'-o','linewidth',2,'color',colors.green,'linewidth',2);
legend('Response-gain set','Contrast-gain set','location','south');
xlabel([]);
ylabel([]);
set(gca,'linewidth',1.5,'fontsize',12,'xlim',[1,16],'ylim',[0.7, 0.92]);

if flag_save == 1
    saveas(gcf, strcat('7D_modelfit','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end


%% Figure 5C. Delta parameters (norm.) across time

win = 1:16; % only 1st and last window (trans. and sust.)
conditions = {'DE','BIN'};
clear a k n
for t = 1:length(win)
        a.units(t,:) = 100.*((PARAMS2.BIN.a(win(t),:) - PARAMS2.DE.a(win(t),:)) ./ (PARAMS2.DE.a(win(t),:)+PARAMS2.BIN.a(win(t),:))); 
        k.units(t,:) = 100.*((PARAMS2.BIN.k(win(t),:) - PARAMS2.DE.k(win(t),:)) ./ (PARAMS2.DE.k(win(t),:)+PARAMS2.BIN.k(win(t),:))); 
        n.units(t,:) = 100.*((PARAMS2.BIN.n(win(t),:) - PARAMS2.DE.n(win(t),:)) ./ (PARAMS2.DE.n(win(t),:)+PARAMS2.BIN.n(win(t),:)));  
end

% average
for t=1:length(win)
    a.avg(t) = nanmean(a.units(t,:),2);
    a.std(t) = nanstd(a.units(t,:),[],2);
    a.err(t) = nanstd(a.units(t,:),[],2) ./ sqrt(size(a.units,2));
    k.avg(t) = nanmean(k.units(t,:),2);
    k.std(t) = nanstd(k.units(t,:),[],2);
    k.err(t) = nanstd(k.units(t,:),[],2) ./ sqrt(size(k.units,2));
    n.avg(t) = nanmean(n.units(t,:),2);
    n.std(t) = nanstd(n.units(t,:),[],2);
    n.err(t) = nanstd(n.units(t,:),[],2) ./ sqrt(size(n.units,2));
    
end

figure('position',[1065,959,550,155.3333333333333]); 
plot(smooth(a.avg),'linewidth',1.5); hold on
plot(smooth(k.avg),'linewidth',1.5);
%plot(smooth(n.avg),'linewidth',1.5);
hl = hline(0); set(hl,'linewidth',1,'color','k')
grid off
set(gca,'linewidth',1,'fontsize',12,...
    'xlim',[1,16],'ylim',[-8 8]);
lgn = legend('R_m_a_x','C_5_0','location','northeast'); legend boxoff
set(lgn,'fontsize',10);

if flag_save == 1
    saveas(gcf, strcat('5C_deltaparams','.svg'));
    disp("Figure saved");
else
    disp("Figure not saved");
end



