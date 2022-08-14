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

    err = log10(OptimizationResults.(profs{a}).error);
    err_mean(a,1) = mean(err);
    err_std(a,1) = std(err)./sqrt(length(err));
    err = log10(OptimizationResults.(profs{a}).RandomResults.error);
    err_mean(a,2) = mean(err);
    err_std(a,2) = std(err)./sqrt(length(err));

    leakage = log10(1-sum(OptimizationResults.(profs{a}).Predata_noRes,2));
    leakage_mean(a,1) = mean(leakage);
    leakage_std(a,1) = std(leakage)./sqrt(length(leakage));
    leakage = log10(1-sum(OptimizationResults.(profs{a}).RandomResults.Predata_noRes,2));
    leakage_mean(a,2) = mean(leakage);
    leakage_std(a,2) = std(leakage)./sqrt(length(leakage));

end


% Plot 
figure;
hold on
subplot(2,1,1);
hold on
errorbar(1:NumofProf, err_mean(:,1),err_std(:,1),'-o');
errorbar(1:NumofProf,err_mean(:,2),err_std(:,2),'-^');
ylabel('log10(RMSE)');
legend('Fitted RMSE','Random RMSE');
xticks(1:NumofProf);
xticklabels(profs);
axis([0 NumofProf+1 -inf inf])
hold off

subplot(2,1,2);
hold on
errorbar(1:NumofProf,leakage_mean(:,1),leakage_std(:,1),'-o');
errorbar(1:NumofProf,leakage_mean(:,2),leakage_std(:,2),'-^');
ylabel('log10(Leakage)');
xticks(1:NumofProf);
xticklabels(profs);
axis([0 NumofProf+1 -inf inf])
legend('Fitted Leakage', 'Random Leakge');
hold off

hold off

end