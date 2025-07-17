
% The predictive ability of baseline deviation values for treatment response in each StD subtype
clear;
load('D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\z_treatment.mat');
load('D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\clinical_base_follow.mat');
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_treatment.mat');

save_path='D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\brain_clinic\predict\';

Subjects_Data = z_base;
Subjects_Scores = clinical_base_follow(:,8);

Covariates = [];
Pre_Method = 'Normalize';
FoldQuantity=5;
C_Range = power(2,-5:10);
Weight_Flag=1;
Permutation_Flag=0;
ResultantFolder = [save_path,'5fold\'];
Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data, Subjects_Scores, FoldQuantity, Pre_Method,C_Range, Weight_Flag, Permutation_Flag, ResultantFolder)
score_predicted = cell2mat(Prediction.Score');
id = cell2mat(Prediction.Origin_ID');
corr_res = corr(score_predicted,Subjects_Scores(id,1));
mse_res = mean((score_predicted/100 - Subjects_Scores(id,1)/100).^2);
mae_res = mean(abs(score_predicted/100 - Subjects_Scores(id,1)/100));
SS_res = sum((score_predicted/100 - Subjects_Scores(id,1)/100).^2);
SS_tot = sum((Subjects_Scores(id,1)/100 - mean(Subjects_Scores(id,1)/100)).^2);
r2_res = 1 - SS_res/SS_tot;

clus_StD=clus_treatment;
ind1=find(clus_StD==1);
ind2=find(clus_StD==2);

Prediction1 = SVR_NFolds_Sort_CSelect(Subjects_Data(ind1,:), Subjects_Scores(ind1,:), FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag, [ResultantFolder,'1\']);
score_predicted1 = cell2mat(Prediction1.Score');
id1 = cell2mat(Prediction1.Origin_ID');
Subjects_Scores1=Subjects_Scores(ind1,:);
corr_res1 = corr(score_predicted1,Subjects_Scores1(id1,1)); %
mse_res1 = mean((score_predicted1/100 - Subjects_Scores1(id,1)/100).^2);
mae_res1 = mean(abs(score_predicted1/100 - Subjects_Scores1(id,1)/100));
SS_res = sum((score_predicted1/100 - Subjects_Scores1(id,1)/100).^2);
SS_tot = sum((Subjects_Scores1(id,1)/100 - mean(Subjects_Scores1(id,1)/100)).^2);
r2_res1 = 1 - SS_res/SS_tot;

Prediction2 = SVR_NFolds_Sort_CSelect(Subjects_Data(ind2,:), Subjects_Scores(ind2,:), FoldQuantity, Pre_Method, C_Range, Weight_Flag, Permutation_Flag, [ResultantFolder,'2\']);
score_predicted2 = cell2mat(Prediction2.Score');
id2 = cell2mat(Prediction2.Origin_ID');
Subjects_Scores2=Subjects_Scores(ind2,:);
corr_res2 = corr(score_predicted2,Subjects_Scores2(id2,1)); %
mse_res2 = mean((score_predicted2/100 - Subjects_Scores2(id,1)/100).^2);
mae_res2 = mean(abs(score_predicted2/100 - Subjects_Scores2(id,1)/100));
SS_res = sum((score_predicted2/100 - Subjects_Scores2(id,1)/100).^2);
SS_tot = sum((Subjects_Scores2(id,1)/100 - mean(Subjects_Scores2(id,1)/100)).^2);
r2_res2 = 1 - SS_res/SS_tot;

% permutation
perm_corr = zeros(1000,1);
Permutation_Flag = 1;
ResultantFolder = [save_path,'5fold\perm\'];
parfor i = 1:1000
    Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data, Subjects_Scores, FoldQuantity, Pre_Method,C_Range, Weight_Flag, Permutation_Flag, ResultantFolder);
    a = cell2mat(Prediction.Score');
    id = cell2mat(Prediction.Origin_ID');
    perm_corr(i) = corr(a,Subjects_Scores(id,1));
    perm_mse(i) = mean((a/100 - Subjects_Scores(id,1)/100).^2);
    perm_mae(i) = mean(abs(a/100 - Subjects_Scores(id,1)/100));
    SS_res = sum((a/100 - Subjects_Scores(id,1)/100).^2);
    SS_tot = sum((Subjects_Scores(id,1)/100 - mean(Subjects_Scores(id,1)/100)).^2);
    perm_r2(i) = 1 - SS_res/SS_tot;
end
p = length(find(perm_corr>corr_res))/1000;
p_mse = length(find(perm_mse<mse_res))/1000;
p_mae = length(find(perm_mae<mae_res))/1000;
p_r2 = length(find(perm_r2>r2_res))/1000;

perm_corr1 = zeros(1000,1);
parfor i = 1:1000
    Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data(ind1,:), Subjects_Scores(ind1,:), FoldQuantity, Pre_Method,C_Range, Weight_Flag, Permutation_Flag, ResultantFolder);
    a = cell2mat(Prediction.Score');
    id = cell2mat(Prediction.Origin_ID');
    Subjects_Scores1=Subjects_Scores(ind1,:);
    perm_corr1(i) = corr(a,Subjects_Scores1(id,1));
    perm_mse1(i) = mean((a/100 - Subjects_Scores1(id,1)/100).^2);
    perm_mae1(i) = mean(abs(a/100 - Subjects_Scores1(id,1)/100));
    SS_res = sum((a/100 - Subjects_Scores1(id,1)/100).^2);
    SS_tot = sum((Subjects_Scores1(id,1)/100 - mean(Subjects_Scores1(id,1)/100)).^2);
    perm_r21(i) = 1 - SS_res/SS_tot;
end
p1 = length(find(perm_corr1>corr_res1))/1000;
p_mse1 = length(find(perm_mse1<mse_res1))/1000;
p_mae1 = length(find(perm_mae1<mae_res1))/1000;
p_r21 = length(find(perm_r21>r2_res1))/1000;

perm_corr2 = zeros(1000,1);
parfor i = 1:1000
    Prediction = SVR_NFolds_Sort_CSelect(Subjects_Data(ind2,:), Subjects_Scores(ind2,:), FoldQuantity, Pre_Method,C_Range, Weight_Flag, Permutation_Flag, ResultantFolder);
    a = cell2mat(Prediction.Score');
    id = cell2mat(Prediction.Origin_ID');
    Subjects_Scores2=Subjects_Scores(ind2,:);
    perm_corr2(i) = corr(a,Subjects_Scores2(id,1));
    perm_mse2(i) = mean((a/100 - Subjects_Scores2(id,1)/100).^2);
    perm_mae2(i) = mean(abs(a/100 - Subjects_Scores2(id,1)/100));
    SS_res = sum((a/100 - Subjects_Scores2(id,1)/100).^2);
    SS_tot = sum((Subjects_Scores2(id,1)/100 - mean(Subjects_Scores2(id,1)/100)).^2);
    perm_r22(i) = 1 - SS_res/SS_tot;
end
p2 = length(find(perm_corr2>corr_res2))/1000;
p_mse2 = length(find(perm_mse2<mse_res2))/1000;
p_mae2 = length(find(perm_mae2<mae_res2))/1000;
p_r22 = length(find(perm_r22>r2_res2))/1000;
