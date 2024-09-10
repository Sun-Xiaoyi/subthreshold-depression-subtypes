% Subtype differences in longitudinal changes in brain deviations
clear;
load('D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\z_treatment.mat');
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_treatment.mat');
title={'subtype','time1','time2'};
for i=1:220
    data=table(clus_treatment,z_base(:,i),z_follow(:,i),'VariableNames',title);
    writetable(data,['D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\brain\data_forR',num2str(i),'.csv']);
end
% run step08_subtype_brain_follow.R

% FDR corrected
clear;
F_P_value=[];
F_P_inter=[];
for i=1:220
    temp=importdata(['D:\Data_Chen\With_DIDA_all_HC\subtype\treatment\brain\result_rANOVA',num2str(i),'.csv']);
    data=zeros(3,2);
    data(1,1)=str2num(temp.textdata{3,6});
    data(1,2)=str2num(temp.textdata{3,7});
    data(2,1)=str2num(temp.textdata{4,6});
    data(2,2)=str2num(temp.textdata{4,7});
    data(3,1)=str2num(temp.textdata{5,6});
    data(3,2)=str2num(temp.textdata{5,7});
    F_P_value=[F_P_value;data];
    F_P_inter=[F_P_inter;data(3,:)];
end
p_corr=gretna_FDR(F_P_value(:,2),0.05);
ind=find(F_P_inter(:,2)<=p_corr);
F_sig=zeros(220,1);
F_sig(ind)=F_P_inter(ind);
show_z(F_sig,'brain_F_allFDR.nii');

