% To investigate the differences in deviation patterns between StD subtypes, 
% we computed the mean deviation map for each subtype and compared them at the network level. 

clear;
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat');
load('D:\Data_Chen\With_DIDA_all_HC\res_norm\res_norm\z_base.mat');
z_clus1=mean(z_base(ind_clus1,:),1);
show_z(z_clus1,'z_clus1.nii');
z_clus2=mean(z_base(ind_clus2,:),1);
show_z(z_clus2,'z_clus2.nii');
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','z_clus1.nii','option_deviation.mat','z_clus1.tif')
BrainNet_MapCfg('BrainMesh_ICBM152_smoothed.nv','z_clus2.nii','option_deviation.mat','z_clus2.tif')

load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat');
load('D:\Data_Chen\With_DIDA_all_HC\res_norm\res_norm\sys_level_dev.mat');
mean_z=sys_level_dev{9};
mean_z_clus1=mean_z(ind_clus1,:);
mean_z_clus2=mean_z(ind_clus2,:);
data{1}=mean_z_clus1;
data{2}=mean_z_clus2;
[T, P] = gretna_TTest2(data);
p_corr=gretna_FDR(P,0.05);

% To assess subtype differences in the overall extent of individual deviations, 
% we calculated the number of regions with extreme deviations, the sum of positive extreme deviations, 
% and the sum of negative extreme deviations for each participant with StD and compared them between subtypes.

clear;
z_base=importdata('D:\Data_Chen\With_DIDA_all_HC\res_norm\res_norm\z_base.mat');

thr=2.6;
z_sig=zeros(size(z_base));
z_sig(find(z_base<-thr))=-1;
z_sig(find(z_base>thr))=1;
num_sub_dev=sum(z_sig~=0,2);

z_pos=z_base;
z_pos(find(z_pos<thr))=0;
sum_z_pos=sum(z_pos,2);

z_neg=z_base;
z_neg(find(z_neg>-thr))=0;
sum_z_neg=sum(z_neg,2);

load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat');
data{1}=[num_sub_dev(ind_clus1,:), sum_z_pos(ind_clus1,:), sum_z_neg(ind_clus1,:)];
data{2}=[num_sub_dev(ind_clus2,:), sum_z_pos(ind_clus2,:), sum_z_neg(ind_clus2,:)];
[T(:,1), T(:,2)] = gretna_TTest2(data);
p_corr=gretna_FDR(T(:,2),0.05);
