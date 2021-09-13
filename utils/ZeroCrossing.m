
clear RES
% Residual and CI
for c = 2:4
    BIN = squeeze(bin.units(c,:,:));
    DE = squeeze(de.units(c,:,:));
    
    
    
    for u = 1:length(IDX)
        nDE(:,u) = DE(:,u) ./ max(DE(:,u));
        nBIN(:,u) = BIN(:,u) ./ max(DE(:,u));
        RES.residual(c-1,:,u) = nBIN(:,u)' - nDE(:,u)';
    end

end

%
newWin = info.sdfWin(251:401); % look for negatives at this point on-ward
for c = 1:3
    x = squeeze(RES.residual(c,251:401,:));
    nc=size(x,2);
    zeroX=zeros(1,nc);
    % number of columns
    
    for i=1:nc
        thresh = 0;
        idx=find(x(:,i)< thresh,1);  % the first negative location in ith column
        if ~isempty(idx)  % if none, find returns empty array
            zeroX(i)=newWin(idx);      % save the valid ones...
        end
        
        zeroX(zeroX == 0) = nan;
    end
    
    timepoints(c,:) = zeroX;
end



figure('position',[839,550,286,202]);
for c = 1:3
    subplot(3,1,c)
    nbins = 20;
    h = histogram(timepoints(c,:),nbins,'orientation','vertical',...
        'faceAlpha',0.5,'linewidth',1,'faceColor',colors.black+0.6,'edgealpha',0.5,'normalization','probability'); hold on
    %histogram(binmod(~mask),'orientation','horizontal','binwidth',0.05,'faceAlpha',0.5,'linewidth',1,'faceColor',colors.purple,'edgealpha',0.3);
    if c == 1
        set(h, 'facecolor', colors.black+0.6);
    elseif c == 2
        set(h, 'facecolor', colors.black+0.4);
    else
        set(h, 'facecolor', colors.black);
    end
    xlim([0.1 0.2]);
    ylim([0 0.30]);
    vl = vline(nanmean(timepoints(c,:),2));
    set(vl,'linewidth',2);
    
    vl = vline(nanmedian(timepoints(c,:),2));
    set(vl,'linewidth',2,'linestyle',':');
    %[h, p, ci, stats] = ttest(binmod,1);
end