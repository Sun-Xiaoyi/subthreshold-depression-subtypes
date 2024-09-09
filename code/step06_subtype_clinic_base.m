% To further validate the clinical differences between the StD subtypes,
% we conducted statistical comparisons of various demographic, clinical, and cognitive variables.

clear;
load('D:\Data_Chen\With_DIDA_all_HC\subtype\baseline\clinic\clinical_cog.mat');
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat');
clinical_cog_clus1=clinical_cog(ind_clus1,:);
clinical_cog_clus2=clinical_cog(ind_clus2,:);

t_p=zeros(23,2);
for i=2:24
    data1=clinical_cog_clus1(:,i);
    data2=clinical_cog_clus2(:,i);
    ind1=find(~isnan(data1));
    data1=data1(ind1);
    ind2=find(~isnan(data2));
    data2=data2(ind2);
    data{1}=data1;
    data{2}=data2;
    [t_p(i-1,1),t_p(i-1,2)]=gretna_TTest2(data);
end

% HAMD item scores
HAMD_item(:,1)=xlsread('D:\Data_Chen\With_DIDA_all_HC\subtype\baseline\clinic\HAMD_item_update_update.xlsx','sheet1','A2:A198');
HAMD_item(:,2:25)=xlsread('D:\Data_Chen\With_DIDA_all_HC\subtype\baseline\clinic\HAMD_item_update_update.xlsx','sheet1','H2:AE198');
HAMD_sum=sum(HAMD_item(:,2:25),2);
HAMD_item(find(isnan(HAMD_sum)),:)=[];
ind1=find(HAMD_item(:,1)==1);
ind2=find(HAMD_item(:,1)==2);
data1=HAMD_item(ind1,2:25);
data2=HAMD_item(ind2,2:25);
[t_p(:,1),t_p(:,2)]=gretna_TTest2({data1,data2});
