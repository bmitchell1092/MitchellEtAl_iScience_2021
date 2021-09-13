function [CRF] = getCRF(PARAMS, UNIT, conditions, bslc, info)

x = 0.01:.01:1;


for uct = 1:size(PARAMS.B,2)
    
    for cond = 1:length(conditions)
        
        for tw = 1:info.windows-1 % all windows but the bsl period
            
            if bslc == true
                CRF.(conditions{cond}).curve(:,tw,uct) = feval(PARAMS.(conditions{cond}).fitobject(tw,uct).f,x)-PARAMS.B(tw,uct)'; % for baseline correct, use tw =3
                CRF.(conditions{cond}).data(:,tw,uct)  = cellfun(@(x) median(x(tw,:)), UNIT(uct).(conditions{cond}).RESP_trls)-PARAMS.B(tw,uct)';
                
            else
                
                CRF.(conditions{cond}).curve(:,tw,uct) = feval(PARAMS.(conditions{cond}).fitobject(tw,uct).f,x)';
                CRF.(conditions{cond}).data(:,tw,uct)  = cellfun(@(x) median(x(tw,:)), UNIT(uct).(conditions{cond}).RESP_trls)';
            end
            
            CRF.(conditions{cond}).error(:,tw,uct) = UNIT(uct).(conditions{cond}).RESP_error(:,tw);
            CRF.(conditions{cond}).sd(:,tw,uct) = UNIT(uct).(conditions{cond}).RESP_sd(:,tw);

            % central tendency
            CRF.(conditions{cond}).curve_sem(:,tw)  = std(CRF.(conditions{cond}).curve(:,tw,:),[],3) ./ sqrt(size(CRF.(conditions{cond}).curve(:,tw,:),3));
            CRF.(conditions{cond}).data_median(:,tw)  = median(CRF.(conditions{cond}).data(:,tw,:),3);
            CRF.(conditions{cond}).data_mean(:,tw)  = mean(CRF.(conditions{cond}).data(:,tw,:),3);
        end
    end
    
end

for cond = 1:length(conditions)
    for tw = 1:info.windows-1
        p = 95;
        a = mean(PARAMS.(conditions{cond}).a(tw,:),2);
        %aC = CIFcn(PARAMS.(conditions{cond}).a(tw,:),p);
        aC = CI95(PARAMS.(conditions{cond}).a(tw,:));
        
        k = mean(PARAMS.(conditions{cond}).k(tw,:),2);
        %kC = CIFcn(PARAMS.(conditions{cond}).k(tw,:),p);
        kC = CI95(PARAMS.(conditions{cond}).k(tw,:));
        
        n = mean(PARAMS.(conditions{cond}).n(tw,:),2);
        %nC = CIFcn(PARAMS.(conditions{cond}).n(tw,:),p);
        nC = CI95(PARAMS.(conditions{cond}).n(tw,:));
        
        b = 0; 
        
        %CRF.(conditions{cond}).curve_avg(:,tw) = a.* ((x.^n) ./ ((x.^n) + (k.^n)) ) + b;
        CRF.(conditions{cond}).curve_avg(:,tw) = ((a.*(x.^n)) ./ ((x.^n) + (k.^n))) + b; 
        CRF.(conditions{cond}).curve_upper(:,tw) = aC(2).* ((x.^nC(2)) ./ ((x.^nC(2)) + (kC(1).^nC(2))) ) + b;
        CRF.(conditions{cond}).curve_lower(:,tw) = aC(1).* ((x.^nC(1)) ./ ((x.^nC(1)) + (kC(2).^nC(1))) ) + b;
    end
end

end
    