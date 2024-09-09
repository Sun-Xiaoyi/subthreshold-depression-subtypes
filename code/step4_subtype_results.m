% The optimal number of clusters was determined via a winner-take-all approach 
% across 21 effective indexes using the NbClust package
clear;
nc_all=zeros(21,1);
for i=1:21
    nc=importdata(['D:\Data_Chen\With_DIDA_all_HC\subtype\res_nbclus\nc_',num2str(i),'.csv']);
    nc_all(i,1)=nc.data(1,1);
end
for i=2:10
    num(i-1,1)=length(find(nc_all==i));
end
num=flip(num);
for i=2:10
    bar(i,num(i-1,1),'BarWidth',0.8,'LineWidth',0.5,'FaceColor',[173 203 222]/255,'FaceAlpha',num(i-1,1)/11, 'EdgeColor',[173 203 222]/255,'Horizontal','on');
    hold on
end
set(gca,'XLim',[0,12],'XTick',[0:4:12]);
set(gca,'YLim',[1,11],'YTick',[2:1:10], 'yticklabel',{'10','9','8','7','6','5','4','3','2'});
set(gca, 'Fontname', 'Arail','FontSize',6);
box off;
set(gca,'linewidth',0.5);
ylabel('Number of clusters')
xlabel('Frequency among all indices')
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 2 4]);
print(gcf,'res_Nbclust.tif','-dtiff','-r1000')

% subtype results
clus=importdata(['D:\Data_Chen\With_DIDA_all_HC\subtype\KMClus2.csv']);
clus=clus.data;
ind_clus1=find(clus==1);
ind_clus2=find(clus==2);
save('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat','clus','ind_clus1','ind_clus2');

% plot euclidean distance of patients
load('D:\Data_Chen\With_DIDA_all_HC\res_norm\res_norm\z_base.txt');
dif=squareform(pdist(z_base,'euclidean'));
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat')
[i1,j1]=sort(sum(dif(ind_clus1,ind_clus1),1));
[i2,j2]=sort(sum(dif(ind_clus2,ind_clus2),1));
dif_sort=dif([ind_clus1(j1);ind_clus2(j2)]',[ind_clus1(j1);ind_clus2(j2)]');
imagesc(dif_sort)
axis square
set(gca,'FontName','Arial','FontSize',25);