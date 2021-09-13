function [OCC] = getOCC(PARAMS,CRF, GAIN, info,conditions, K, method)

clear OCC*

ftw = info.windows-1; %% this is the full time window of the response
x = 0.01:.01:1;

switch method
    case 'rmax'
        
        OCC.ocularity = (PARAMS.DE.a(ftw,:) - PARAMS.NDE.a(ftw,:)) ./ (PARAMS.DE.a(ftw,:));
    case 'resp'
        
    case 'auc'
        
        for u = 1:length(PARAMS.B) % unit
            de(u) = trapz(feval(PARAMS.DE.fitobject(ftw,u).f,x));
            nde(u) = trapz(feval(PARAMS.NDE.fitobject(ftw,u).f,x));
        end
        
        OCC.ocularity = (de-nde) ./ (nde + de);
        OCC.ocularity_abs = abs((de-nde) ./ (nde + de));
        
    case 'diUnitTuning'
        OCC.ocularity = IDX(:).occ(3);
end


[B,I] = sort(OCC.ocularity,'ComparisonMethod','abs');
OCC.sorted = [I; B];
varc = @(oldvar) mat2cell(oldvar(:), [fix(numel(oldvar)/K) *[ones(K-1,1)]', numel(oldvar)-(K-1)*fix(numel(oldvar)/K)], 1);     % Create New Matrix From Original Vector
OCC.groups = varc(I); % rows are: low to high
temp = varc(B);
OCC.groups(:,2) = temp;

for o = 1:K
    OCC.groups{o,3} = [abs(OCC.groups{o,2}(1)),abs(OCC.groups{o,2}(end))];
end
OCC.lengths = nan(K,1);
for s = 1:length(OCC.lengths)
    OCC.lengths(s) = numel(OCC.groups{s,1});
end

% remove units if there's unequal groups
if ~isequal(length(OCC.groups{end,1}),length(OCC.groups{end-1,1}))
    while length(OCC.groups{end,1}) - length(OCC.groups{end-1,1}) ~= 0
        OCC.groups{end,1}(end) = [];
    end
end

OCC.label = nan(1,length(OCC.ocularity));

for o = 1:K
    OCC.label(OCC.groups{o,1}) = o;
end



% pull out data into ocularity groups
for o = 1:K
    
    for cond = 1:length(conditions)
        
        
        OCC.(conditions{cond}){o}.curve = CRF.(conditions{cond}).curve(:,:,OCC.groups{o,1});
        OCC.(conditions{cond}){o}.data  = CRF.(conditions{cond}).data(:,:,OCC.groups{o,1});
        OCC.(conditions{cond}){o}.error  = CRF.(conditions{cond}).error(:,:,OCC.groups{o,1});
        
        % central tendency
        OCC.(conditions{cond}){o}.data_avg  = mean(OCC.(conditions{cond}){o}.data,3);
        OCC.(conditions{cond}){o}.data_sem  = std(OCC.(conditions{cond}){o}.data,[],3) ./ sqrt(OCC.lengths(o));
        
        
        for tw = 1:info.windows-1
            p = 95;
            a = mean(PARAMS.(conditions{cond}).a(tw,OCC.groups{o,1}));
            aC = CIFcn(PARAMS.(conditions{cond}).a(tw,OCC.groups{o,1}),p);
            
            k = mean(PARAMS.(conditions{cond}).k(tw,OCC.groups{o,1}));
            kC = CIFcn(PARAMS.(conditions{cond}).k(tw,OCC.groups{o,1}),p);
            
            n = mean(PARAMS.(conditions{cond}).n(tw,OCC.groups{o,1}));
            nC = CIFcn(PARAMS.(conditions{cond}).n(tw,OCC.groups{o,1}),p);
            
            b = 0; %mean(PARAMS2.B);
            
            OCC.(conditions{cond}){o}.curve_avg(:,tw) = a.* ((x.^n) ./ ((x.^n) + (k.^n)) ) + b;
            OCC.(conditions{cond}){o}.curve_upper(:,tw) = aC(2).* ((x.^nC(2)) ./ ((x.^nC(2)) + (kC(2).^nC(2))) ) + b;
            OCC.(conditions{cond}){o}.curve_lower(:,tw) = aC(1).* ((x.^nC(1)) ./ ((x.^nC(1)) + (kC(1).^nC(1))) ) + b;
        end
    end
    
%     % gain models
%     if ~isempty('GAIN')
%         model = {'response','contrast','exponent'};
%         
%         for tw = 1:info.windows-1
%             for m = 1:length(model)
%                 A = mean(PARAMS.DE.a(tw,OCC.groups{o,1}),2);
%                 Ac = CI95(PARAMS.DE.a(tw,OCC.groups{o,1}));
%                 
%                 K = mean(PARAMS.DE.k(tw,OCC.groups{o,1}),2);
%                 Kc = CI95(PARAMS.DE.k(tw,OCC.groups{o,1}));
%                 
%                 N = mean(PARAMS.DE.n(tw,OCC.groups{o,1}),2);
%                 Nc = CI95(PARAMS.DE.n(tw,OCC.groups{o,1}));
%                 
%                 G = mean(GAIN.(model{m}).value(tw,OCC.groups{o,1}),2);
%                 Gc = CI95(GAIN.(model{m}).value(tw,OCC.groups{o,1}));
%                 
%                 B = 0;
%                 
%                 switch model{m}
%                     case 'response'
%                         OCC.(model{m}){o}.curve_avg(:,tw) = G.*((A.*(x.^N)) ./ ((x.^N) + (K.^N))) + B;
%                         OCC.(model{m}){o}.curve_upper(:,tw) = (Gc(2).*Ac(2)).*((x.^Nc(2)) ./ ((x.*Nc(2)) + (Kc(2).^Nc(2)))) + B;
%                         OCC.(model{m}){o}.curve_lower(:,tw) = (Gc(1).*Ac(1)).*((x.^Nc(1)) ./ ((x.*Nc(1)) + (Kc(1).^Nc(1)))) + B;
%                     case 'contrast'
%                         
%                         OCC.(model{m}){o}.curve_avg(:,tw) = ((A.*(G.*x.^N)) ./ ((G.*x.^N) + (K.^N))) + B;
%                         OCC.(model{m}){o}.curve_upper(:,tw) = Ac(2).*((Gc(2).*(x.^Nc(2))) ./ ((Gc(2).*(x.^Nc(2))) + (Kc(2).^Nc(2)))) + B;
%                         OCC.(model{m}){o}.curve_lower(:,tw) = Ac(1).*((Gc(1).*(x.^Nc(1))) ./ ((Gc(1).*(x.^Nc(1))) + (Kc(1).^Nc(1)))) + B;
%                     case 'exponent'
%                         
%                         OCC.(model{m}){o}.curve_avg(:,tw) = ((A.*(x.^(G.*N))) ./ ((x.^(G.*N)) + (K.^(G.*N)))) + B;
%                         OCC.(model{m}){o}.curve_upper(:,tw) = ((Ac(2).*(x.^(Gc(2).*Nc(2)))) ./ ((x.^(Gc(2).*Nc(2))) + (Kc(2).^(Gc(2).*Nc(2))))) + B;
%                         OCC.(model{m}){o}.curve_lower(:,tw) = ((Ac(1).*(x.^(Gc(1).*Nc(1)))) ./ ((x.^(Gc(1).*Nc(1))) + (Kc(1).^(Gc(1).*Nc(1))))) + B;
%                 end
%             end
%         end
%     end
end
    
end