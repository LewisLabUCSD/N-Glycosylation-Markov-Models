function     OptimizationResults = FilterOptimizationResults(Prof, OptimizationResults, Method, removeOutlierFlag)

%% load fval data
fval = OptimizationResults.(Prof).fval;
xval = OptimizationResults.(Prof).xval;
selflag = true(size(xval,1),1);

%% Apply filtering method

if strcmp(Method,'KernalDensity') && size(xval,1)>1
    
    % compute the number of points (resolution) for kernal smoothing
    % estimation
    Num = max(abs(diff(fval)))/min(abs(diff(fval)))*10;
    if Num>7000
        Num = 7000;
    end

    [f,xi] = ksdensity(fval,'NumPoints',Num);
    TF = islocalmax([0,f,0]); % pad variable to reveal maxima in the range
    TF = xi(TF(2:end-1));

   % Find inflection points
    f2= diff(diff(f));
    infpoint = [];
    for a = 1:length(f2)-1
        if f2(a)<=0 && f2(a+1)>0
            infpoint = [infpoint,xi(a)];
        end
    end
    TF = sort([TF,infpoint]);

    % cluster fvals based on distance to maxima and inflection points
    clusters = zeros(size(fval));
    for a = 1:length(fval)
        [~,clusters(a)] = min(abs(fval(a)-TF));
    end

    selflag = clusters == 1;

elseif strcmp(Method,'Outlier') && size(xval,1)>1

    selflag = ~isoutlier(fval);

end


if removeOutlierFlag
   OptimizationResults.(Prof).fval = fval(selflag);
   OptimizationResults.(Prof).xval = xval(selflag,:);
   OptimizationResults.(Prof).LacNAcLenPenalty = OptimizationResults.(Prof).LacNAcLenPenalty(selflag);

   % further remove outliers
   selflag = ~isoutlier(OptimizationResults.(Prof).fval);
   OptimizationResults.(Prof).fval = OptimizationResults.(Prof).fval(selflag);
   OptimizationResults.(Prof).xval = OptimizationResults.(Prof).xval(selflag,:);
   OptimizationResults.(Prof).LacNAcLenPenalty = OptimizationResults.(Prof).LacNAcLenPenalty(selflag);
end

end