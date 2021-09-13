function [PARAMS,errlog] = getPARAMS(UNIT,IDX,info, bslc)
%CRFparams.m. Measure CRFs for UNIT structure.
%   UNIT is the structure with the data, flag_individualplots for visualization
condition = upper({'de','nde','bin'});
errlog = false(length(UNIT),1);
for unit = 1:length(IDX)
    flag = 0;
    contrast = IDX(unit).cLevels(:,2); contrast(1) = 0.01;
    for tw = 1:info.windows-1  % for all but the bsl window
        for cond = 1:length(condition)
            
            
            if flag == 1
                errlog(unit,:) = true; %'flag_break';
                break
            end
            
            % note, average response is mean
            try
            average_resp = cellfun(@(x) median(x(tw,:)), UNIT(unit).(condition{cond}).RESP_trls)';
            max_resp = cellfun(@(x) max(x(tw,:)), UNIT(unit).(condition{cond}).RESP_trls)';
            catch
                errlog(unit,:) = true;
                break
            end
            
            if cond == 1 %DE
                B = average_resp(1);
            end
            
            if bslc == true
                average_resp = average_resp - B;
                max_resp = max_resp - B;
            end
            
            fop = fitoptions('Method','NonlinearLeastSquares');
            fcn = @(a,k,n,x) a.* ((x.^n) ./ ((x.^n) + (k.^n)) ) + B; % B = DE baseline
            %%%% Fitting limits %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %                [ a (rmax)                K  n] %
            fop.Upper      = [         max_resp(end) .80  5];%
            fop.Lower      = [                    1  .05  1];%
            fop.StartPoint = [1.1*average_resp(end)  .40  2];%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % fit
            try
                [f, gof] = fit(contrast,average_resp,fcn,fop);
            catch
                warning('Unit %d was not fit',unit)
                errlog(unit,:) = true;
                flag = 1;
                break
            end
            
            PARAMS.(condition{cond}).fitobject(tw,unit).f = f;
            PARAMS.(condition{cond}).a(tw,unit) = f.a;
            PARAMS.(condition{cond}).k(tw,unit) = f.k;
            PARAMS.(condition{cond}).n(tw,unit) = f.n;
            PARAMS.(condition{cond}).rsquare(tw,unit) = gof.rsquare;
            PARAMS.(condition{cond}).gof(tw,unit) = gof;
            
        end
        
        PARAMS.B(tw,unit) = B;
    end
    
end
end

