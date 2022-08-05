function OptimizationResults = runGlobalOptimization(Prof, OptimizationResults,OptimizationProblem)

method = OptimizationProblem.optimproblem.solver;
repeatflag = true; % repeat the fitting if algorithm failed to converge in rare cases
count = 0;

%%%%%%%%%%%%%%%%%%%%% Fitting %%%%%%%%%%%%%%%%%%%%%
if strcmp(method,'patternsearch')

    while repeatflag
        [xval,fval,~,output] = patternsearch(OptimizationProblem.optimproblem);
        count = count + 1;
        if output.iterations > 500 || count>2
            repeatflag = false;
        end
    end

elseif strcmp(method,'particleswarm')

    while repeatflag
        [xval,fval,~,output] = particleswarm(OptimizationProblem.optimproblem);
        count = count + 1;
        if output.iterations > 500 || count>2
            repeatflag = false;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%% Record Fitting Results %%%%%%%%%%%%%%%%%%%%%

% OptimizationResults = RecordOptimizationResults(ProfSel{a}, OptimizationResults, xval, fval,OptimizationProblem);

% Input
% 1. Prof: the string of the selected profile name for fitting
% 2. OptimizationResults: struct variable used to store all
% optimization results.
% 3. xval: optimized model parameters of the current run
% 4. fval: minimized objective function error of the current run
% 5. OptimizationProblem: struct storing the fitting problem set up
% created by the function SetUpFittingProblem.

% Output
% c. OptimizationResults: updated OptimizationResults variable. Has
% the following fields:
%   a.xval: double-type fitted transition probablities. Each row represent a fitted transition
%     probabilities vectors and each column represents a reaction type (in the order of variable RxnTypes).
%     Therefore, the number of rows is equal to the number of fitted
%     models for a specific profile, whereas the number of columns is
%     equal to the the length of variable RxnTypes.
%   b.fval: the errors computed by the objective function for each
%   fitted model.
%   c. OptimizationProblem: the struct variable storing the fitting
%   problem setups

OptimizationResults = RecordOptimizationResults(Prof, OptimizationResults, xval, fval,OptimizationProblem);

end