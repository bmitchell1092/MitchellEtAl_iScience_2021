function [ExpVar] = fitBIN(UNIT, PARAMS, IDX, info)
% NDEfit.m. Refit NDE with some fixed parameters
%   UNIT is the structure with the data.

param = {'a','k','n'};

clear NDE



for uct = 1:size(PARAMS.B,2)
    contrast = IDX(uct).cLevels(:,2); contrast(1) = 0.01;
    
    for tw = 1:info.windows-1
        
        % Grab parameters from the DE
        A = PARAMS.DE.a(tw,uct); K = PARAMS.DE.k(tw,uct); N = PARAMS.DE.n(tw,uct);
        B = PARAMS.B(tw,uct);
        
        % Collect BIN responses
        average_resp = cellfun(@(x) median(x(tw,:)), UNIT(uct).BIN.RESP_trls)';
        max_resp = cellfun(@(x) max(x(tw,:)), UNIT(uct).BIN.RESP_trls)';
        
        for p = 1:length(param)
            switch param{p}
                case 'a'
                    fop = fitoptions('Method','NonlinearLeastSquares');
                    fcn = @(a,x) a.*((x.^N) ./ ((x.^N) + (K.^N)) ) + B; % a free to vary
                    %%%% Fitting limits %%%%%%%%%%%%%%%%%%%%%%
                    %                [ a                     %
                    fop.Upper      = [ max_resp(end)         ];%
                    fop.Lower      = [ 1                     ];%
                    fop.StartPoint = [ 1.1*average_resp(end) ];%
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    [f, gof] = fit(contrast,average_resp,fcn,fop);
                    value = f.a;
                    
                case 'k'
                    fop = fitoptions('Method','NonlinearLeastSquares');
                    fcn = @(k,x) A.*((x.^N) ./ ((x.^N) + (k.^N)) ) + B; % k free to vary
                    %%%% Fitting limits %%%%%%%%%%%%%%%%%%%%%%
                    %                [ k                     %
                    fop.Upper      = [ 0.90                  ];%
                    fop.Lower      = [ .10                   ];%
                    fop.StartPoint = [ 0.40                  ];%
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    [f, gof] = fit(contrast,average_resp,fcn,fop);
                    value = f.k;
                    
                case 'n'
                    fop = fitoptions('Method','NonlinearLeastSquares');
                    fcn = @(n,x) A.*((x.^n) ./ ((x.^n) + (K.^n)) ) + B; % n free to vary
                    %%%% Fitting limits %%%%%%%%%%%%%%%%%%%%%%
                    %                [ n                     %
                    fop.Upper      = [ 10                    ];%
                    fop.Lower      = [ 1                     ];%
                    fop.StartPoint = [ 2                     ];%
                    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    [f, gof] = fit(contrast,average_resp,fcn,fop);
                    value = f.n;

            end
            
            ExpVar.(param{p}).fitobject(tw,uct).f = f;
            ExpVar.(param{p}).curve(tw,:,uct) = feval(f,[0.01:0.01:1])-B';
            ExpVar.(param{p}).value(tw,uct) = value;
            
            gof.rsquare(gof.rsquare < 0) = 0;
            ExpVar.(param{p}).rsquare(tw,uct) = gof.rsquare;
            ExpVar.(param{p}).gof(tw,uct) = gof;
            
            
        end
    end
    
end

% Once all units collected, can calc CI
for p = 1:length(param)
    for tw = 1:info.windows-1
        P = 95;
        ExpVar.(param{p}).rsquare_CI(tw,:) = CIFcn(ExpVar.(param{p}).rsquare(tw,:),P);
    end
end