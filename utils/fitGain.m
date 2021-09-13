function [GAIN] = fitGain(UNIT, PARAMS, IDX, info)
% NDEfit.m. Refit NDE with some fixed parameters
%   UNIT is the structure with the data.

model = {'response','contrast','exponent','additive'};
x = 0.01:.01:1;

for uct = 1:size(PARAMS.B,2)
    contrast = IDX(uct).cLevels(:,2); contrast(1) = 0.01;
    
    for tw = 1:info.windows-1
        
        % Grab parameters from the DE
        A = PARAMS.DE.a(tw,uct); K = PARAMS.DE.k(tw,uct); N = PARAMS.DE.n(tw,uct);
        B = PARAMS.B(tw,uct);
        
        % Collect BIN responses
        average_resp = cellfun(@(x) median(x(tw,:)), UNIT(uct).BIN.RESP_trls)';
        
        for m = 1:length(model)
            switch model{m}
                
                case 'response'
                    fop = fitoptions('Method','NonlinearLeastSquares');
                    fcn = @(G,x) G.*((A.*(x.^N)) ./ ((x.^N) + (K.^N))) + B; % G free to vary
                    %%%% Fitting limits %%%
                    %                G    %
                    fop.Upper      = 10;  %
                    fop.Lower      = 0.25;% 
                    fop.StartPoint = 1;   %
                    % %%%%%%%%%%%%%%%%%%%%%
                    
                    [f, gof] = fit(contrast,average_resp,fcn,fop);
                    value = f.G;
                    
                case 'contrast'
                    fop = fitoptions('Method','NonlinearLeastSquares');
                    fcn = @(G,x) ((A.*(G.*x.^N)) ./ ((G.*x.^N) + (K.^N))) + B; % G free to vary
                    %%%% Fitting limits %
                    %                G 
                    fop.Upper      = 10;%
                    fop.Lower      = 0.25;
                    fop.StartPoint = 1;
                    % %%%%%%%%%%%%%%%%%%
                    
                    [f, gof] = fit(contrast,average_resp,fcn,fop);
                    value = f.G;
                    
                case 'exponent'
                    
                    fop = fitoptions('Method','NonlinearLeastSquares');
                    fcn = @(G,x) ((A.*(x.^(G.*N))) ./ ((x.^(G.*N)) + (K.^(G.*N)))) + B; % G free to vary
                    %%%% Fitting limits %
                    %                G
                    fop.Upper      = 10;%
                    fop.Lower      = 0.1;
                    fop.StartPoint = 0.1;
                    % %%%%%%%%%%%%%%%%%%
                    
                    [f, gof] = fit(contrast,average_resp,fcn,fop);
                    value = f.G;
                    
                case 'additive'
                    fop = fitoptions('Method','NonlinearLeastSquares');
                    fcn = @(G,x) G + ((A.*(x.^N)) ./ ((x.^N) + (K.^N))) + B; % G free to vary
                    %%%% Fitting limits %%%
                    %                G    %
                    fop.Upper      = 50;  %
                    fop.Lower      = 0.25;% 
                    fop.StartPoint = 1;   %
                    % %%%%%%%%%%%%%%%%%%%%%
                    
                    [f, gof] = fit(contrast,average_resp,fcn,fop);
                    value = f.G;
                    
            end
            
            GAIN.(model{m}).fitobject(tw,uct).f = f;
            GAIN.(model{m}).curve(tw,:,uct) = feval(f,x)-B';
            GAIN.(model{m}).value(tw,uct) = value;
            
            gof.rsquare(gof.rsquare < 0) = 0;
            GAIN.(model{m}).rsquare(tw,uct) = gof.rsquare;
            GAIN.(model{m}).gof(tw,uct) = gof;
            
            
        end
    end
    
end

% Once all units collected, can calc CI
for m = 1:length(model)
    for tw = 1:info.windows-1
        P = 95;
        GAIN.(model{m}).rsquare_CI(tw,:) = CIFcn(GAIN.(model{m}).rsquare(tw,:),P);
    end
end

% Calculate average curves and 95% CI


for tw = 1:info.windows-1
    for m = 1:length(model)
        A = mean(PARAMS.DE.a(tw,:),2);
        Ac = CI95(PARAMS.DE.a(tw,:));
        
        K = mean(PARAMS.DE.k(tw,:),2);
        Kc = CI95(PARAMS.DE.k(tw,:));
        
        N = mean(PARAMS.DE.n(tw,:),2);
        Nc = CI95(PARAMS.DE.n(tw,:));
        
        G = mean(GAIN.(model{m}).value(tw,:),2);
        Gc = CI95(GAIN.(model{m}).value(tw,:));
        
        B = 0;
        
        switch model{m}
            case 'response'
                GAIN.(model{m}).curve_avg(:,tw) = G.*((A.*(x.^N)) ./ ((x.^N) + (K.^N))) + B;
                GAIN.(model{m}).curve_upper(:,tw) = Gc(2).*((Ac(2).*(x.^Nc(2))) ./ ((x.^Nc(2)) + (Kc(2).^Nc(2)))) + B;
                GAIN.(model{m}).curve_lower(:,tw) = Gc(1).*((Ac(1).*(x.^Nc(1))) ./ ((x.^Nc(1)) + (Kc(1).^Nc(1)))) + B;
            case 'contrast'
                
          
                GAIN.(model{m}).curve_avg(:,tw) = ((A.*(G.*x.^N)) ./ ((G.*x.^N) + (K.^N))) + B;
                GAIN.(model{m}).curve_upper(:,tw) = Ac(2).*((Gc(2).*(x.^Nc(2))) ./ ((Gc(2).*(x.^Nc(2))) + (Kc(2).^Nc(2)))) + B;
                GAIN.(model{m}).curve_lower(:,tw) = Ac(1).*((Gc(1).*(x.^Nc(1))) ./ ((Gc(1).*(x.^Nc(1))) + (Kc(1).^Nc(1)))) + B;
            case 'exponent'
                
                GAIN.(model{m}).curve_avg(:,tw) = ((A.*(x.^(G.*N))) ./ ((x.^(G.*N)) + (K.^(G.*N)))) + B;
                GAIN.(model{m}).curve_upper(:,tw) = ((Ac(2).*(x.^(Gc(2).*Nc(2)))) ./ ((x.^(Gc(2).*Nc(2))) + (Kc(2).^(Gc(2).*Nc(2))))) + B;
                GAIN.(model{m}).curve_lower(:,tw) = ((Ac(1).*(x.^(Gc(1).*Nc(1)))) ./ ((x.^(Gc(1).*Nc(1))) + (Kc(1).^(Gc(1).*Nc(1))))) + B;
                
            case 'additive'
                
                GAIN.(model{m}).curve_avg(:,tw) = G + ((A.*(x.^N)) ./ ((x.^N) + (K.^N))) + B;
                GAIN.(model{m}).curve_upper(:,tw) = Gc(2) + ((Ac(2).*(x.^Nc(2))) ./ ((x.^Nc(2)) + (Kc(2).^Nc(2)))) + B;
                GAIN.(model{m}).curve_lower(:,tw) = Gc(1) + ((Ac(1).*(x.^Nc(1))) ./ ((x.^Nc(1)) + (Kc(1).^Nc(1)))) + B;
        end
    end
    
end