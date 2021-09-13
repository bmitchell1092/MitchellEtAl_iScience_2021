clear contrast

for u = 1:length(IDX)
    dicontrast = IDX(u).X.dicontrasts;
    
    if length(dicontrast) == 3
        contrast(:,u) = dicontrast;
    
    elseif dicontrast == [0.200000000000000,0.400000000000000,0.600000000000000,0.800000000000000,1]
        contrast(:,u) = dicontrast([1,2,4]);

    elseif dicontrast == [0.0500000000000000,0.150000000000000,0.300000000000000,0.500000000000000,1]
        contrast(:,u) = dicontrast([3,4,5]);
    end
end


contrasts = reshape(contrast',[],942);



dependent_measures.additivity = additivity';
dependent_measures.binocular_modulation = binocular_modulation;
dependent_measures.facil_magnitude = magnitudes;
dependent_measures.facil_duration = facil_duration;