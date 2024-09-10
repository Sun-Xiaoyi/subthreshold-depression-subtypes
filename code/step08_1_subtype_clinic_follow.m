
% Subtype differences in HDRS changes
clear;
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_treatment.mat');
load('D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\clinical_base_follow.mat');
title={'subtype','time1','time2'};
data=table(clus_treatment,clinical_base_follow(:,5),clinical_base_follow(:,6),...
    'VariableNames',title);
writetable(data,'D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\clinical\data_forR.csv');
%run step08_subtype_clinic_follow.R
