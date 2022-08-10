function VisualizePerformance(OptimizationResults)

% Initiate variables
profs = fieldnames(OptimizationResults);
err_mean = zeros(length(profs),2);
err_std = zeros(length(profs),2);
leakage_mean = zeros(length(profs),2);
leakage_std = zeros(length(profs),2);
NumofProf = length(profs);

% Compute variables
for a = 1:NumofProf

    err = OptimizationResults.(profs{a}).error;
    err_mean(a,1) = mean(err);
    err_std(a,1) = std(err);
    err = OptimizationResults.(profs{a}).RandomResults.error;
    err_mean(a,2) = mean(err);
    err_std(a,2) = std(err);

    leakage = 1-sum(OptimizationResults.(profs{a}).Predata_noRes,2);
    leakage_mean(a,1) = mean(leakage);
    leakage_std(a,1) = std(leakage);
    leakage = 1-sum(OptimizationResults.(profs{a}).RandomResults.Predata_noRes,2);
    leakage_mean(a,2) = mean(leakage);
    leakage_std(a,2) = std(leakage);

end

err_mean = log10(err_mean);
err_std = log10(err_std);
leakage_mean = log10(leakage_mean);
leakage_std = log10(leakage_std);


% Plot 
figure;
hold on
yyaxis left
errorbar(1:NumofProf, err_mean(:,1),err_std(:,1),'-o');
errorbar(1:NumofProf,err_mean(:,2),err_std(:,2),'-^');
ylabel('log10(RMSE)');
yyaxis right
errorbar(1:NumofProf,leakage_mean(:,1),leakage_std(:,1),'-o');
errorbar(1:NumofProf,leakage_mean(:,2),leakage_std(:,2),'-^');
ylabel('log10(Leakage)');
xticks(1:NumofProf);
xticklabels(profs);
axis([0 NumofProf+1 -inf inf])
legend('Fitted RMSE','Randome RMSE', 'Fitted Leakage', 'Random Leakge');
hold off

end