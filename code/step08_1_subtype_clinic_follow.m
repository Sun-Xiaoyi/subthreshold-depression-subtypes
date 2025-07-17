
% Subtype differences in HDRS changes
clear;
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_treatment.mat');
load('D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\clinical_base_follow.mat');
title={'subtype','time1','time2'};
data=table(clus_treatment,clinical_base_follow(:,5),clinical_base_follow(:,6),...
    'VariableNames',title);
writetable(data,'D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\clinical\data_forR.csv');
%run step08_2_subtype_clinic_follow.R

% Subtype differences in HDRS item changes
clear;
HDRS_items=readtable("D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\clinical\BLT_HDRS_item.xlsx");
HDRS_item_baseline=HDRS_items(:,3:26);
HDRS_item_follow=HDRS_items(:,28:51);

load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base_follow.mat');
load('D:\Data_Chen\With_DIDA_age_match_HC\subtype\treatment\clinical_base_follow.mat');
ind=find(clinical_base_follow(:,1)==1);

title={'subtype','time1','time2'};
for i=1:24
    data=table(clus_base_follow(ind,:),HDRS_item_baseline{:,i},HDRS_item_follow{:,i},'VariableNames',title);
    writetable(data,['D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\clinical\HDRS_item\data_forR',num2str(i),'.csv']);
end
%run step08_2_subtype_clinic_follow.R
