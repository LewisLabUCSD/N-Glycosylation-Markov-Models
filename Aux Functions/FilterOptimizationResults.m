function     OptimizationResults = FilterOptimizationResults(Prof, OptimizationResults, Method, removeOutlierFlag)

%% load fval data
fval = OptimizationResults.(Prof).fval;
xval = OptimizationResults.(Prof).xval;

%% Apply filtering method

if strcmp(Method,'KernalDensity')
    
    % compute the number of points (resolution) for kernal smoothing
    % estimation
    Num = max(abs(diff(fval)))/min(abs(diff(fval)));
    if Num>1000
        Num = 1000;
    end

    [f,xi] = ksdensity(fval,'NumPoints',Num);
    TF = islocalmax([0,f,0]); % pad variable to reveal maxima on the ranges
    TF = xi(TF(2:end-1));
    
    % cluster fvals based on distance to maxima
    clusters = zeros(size(fval));
    for a = 1:length(fval)
        [~,clusters(a)] = min(abs(fval(a)-TF));
    end

    selflag = clusters == 1;

elseif strcmp(Method,'Outlier')

    selflag = ~isoutlier(fval);

end

if removeOutlierFlag
   OptimizationResults.(Prof).fval = fval(selflag);
   OptimizationResults.(Prof).xval = xval(selflag,:);
end

end