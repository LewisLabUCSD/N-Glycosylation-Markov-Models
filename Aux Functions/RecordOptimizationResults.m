function OptimizationResults = RecordOptimizationResults(Prof, OptimizationResults, xval, fval,OptimizationProblem)

% Load variables
Rxn_idx = OptimizationProblem.Rxn_idx;
AppliedGeneidx = OptimizationProblem.AppliedGeneidx;
stericRxns = OptimizationProblem.stericRxns;
WTSteric = OptimizationProblem.WTSteric;
stericRxnsIdx = length(AppliedGeneidx)-length(stericRxns)+1:length(AppliedGeneidx);
UseWTStericFlag = OptimizationProblem.UseWTStericFlag;


% determine positions of xval values
if OptimizationProblem.StericFlag
   xval_full = zeros(1, length(AppliedGeneidx));
else
    xval_full = zeros(1, length(OptimizationProblem.optimproblem.x0)+3);
end

if OptimizationProblem.StericFlag && OptimizationProblem.UseWTStericFlag
    xval_full(Rxn_idx) = xval;
    xval_full(stericRxnsIdx) = WTSteric;
elseif OptimizationProblem.StericFlag && ~OptimizationProblem.UseWTStericFlag
    xval_full(Rxn_idx) = xval(1:length(Rxn_idx));
    xval_full(stericRxnsIdx) = xval(length(Rxn_idx)+1:end);
elseif ~OptimizationProblem.StericFlag
    xval_full(Rxn_idx) = xval;
end

% create a sub-struct if the current Prof does not exist in OptimizationResults
if ~any(strcmp(Prof,fieldnames(OptimizationResults)))
    OptimizationResults.(Prof).xval = xval_full;
    OptimizationResults.(Prof).fval = fval;
else % otherwise concatenate
    OptimizationResults.(Prof).xval = [OptimizationResults.(Prof).xval;xval_full];
    OptimizationResults.(Prof).fval = [OptimizationResults.(Prof).fval;fval];
end

OptimizationResults.(Prof).OptimizationProblem = OptimizationProblem;

end