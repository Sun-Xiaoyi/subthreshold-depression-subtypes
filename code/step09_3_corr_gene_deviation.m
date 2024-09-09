% Associations between the transcriptional profiles and FCS deviations

% load parcellation, Z-maps, gene
clear;
path='D:\Data_Chen\With_DIDA_all_HC\subtype\gene_abagen\';

hdr_par = spm_vol([path,'code\L_shen268_group.nii']);
vol_par =spm_read_vols(hdr_par);

load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat');
load('D:\Data_Chen\With_DIDA_all_HC\res_norm\res_norm\z_base.mat');
z_clus1=mean(z_base(ind_clus1,110:220),1);

temp =  importdata([path, 'res\L_expression_filter.csv']);
expression = temp.data;
gene_name = temp.textdata;

% remove missing roi
missingdata_regions=find(isnan(expression(:,2)));
region_ind=setdiff(expression(:,1),missingdata_regions);

group_express=expression(region_ind,2:end);
gene_name = gene_name(2:end);

GENEdata=group_express;
MRIdata=z_clus1(region_ind)';

[I J]=find(isnan(GENEdata));
for i=1:length(I)
    GENEdata(I(i),J(i))=0;
end

% PLS_calculation
Y = zscore(MRIdata);
dim = 10;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(GENEdata,Y,dim,'CV',dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);

%Align PLS components with desired direction
R1 = corr([XS(:,1),XS(:,2),XS(:,3),XS(:,4),XS(:,5)],MRIdata);
if R1(1,1)<0
    XS(:,1)=-1*XS(:,1);
end
if R1(2,1)<0
    XS(:,2)=-1*XS(:,2);
end
if R1(3,1)<0
    XS(:,3)=-1*XS(:,3);
end
if R1(4,1)<0
    XS(:,4)=-1*XS(:,4);
end
if R1(5,1)<0
    XS(:,5)=-1*XS(:,5);
end

% permutation test
load('D:\Data_Chen\With_DIDA_all_HC\subtype\gene_abagen\res\surrogate_maps_z_clus1.mat');
vol_mask=spm_read_vols(spm_vol([path,'code\shen268_group.nii']));
ind = find(vol_mask);
PCTVARrand = zeros(10000,10);
Rsq = zeros(10000,1);
parfor j = 1:10000
    disp(j);
    z_sur = zeros(111,1);
    Y1 = zeros(hdr_par.dim);    
    Y1(ind) = surrogate_maps(j,:);
    for i=1:111
        z_sur(i,1)=mean(Y1(vol_par==i));
    end
    MRIdata_s = z_sur(region_ind);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(GENEdata,MRIdata_s,dim);
    PCTVARrand(j,:)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim);    
end
p_single = zeros(1,10);
for l=1:dim
    p_single(l)=length(find(PCTVARrand(:,l)>=PCTVAR(2,l)))/10000;
end
p_cum = length(find(Rsq>=Rsquared))/10000;
myStats=[PCTVAR; p_single];
csvwrite([path,'res\L_PLS_stats_z_clus1.csv'],myStats);

% Draw variance explanation
clear;
path='D:\Data_Chen\With_DIDA_all_HC\subtype\gene_abagen\';
load([path,'res\L_PLS_stats_z_clus1.csv']);
py = plot(sort(L_PLS_stats_z_clus1(2,:)','descend'),'.-','color', [181,43,33]./255,'LineWidth',2, ...
   'MarkerEdgeColor',[181,43,33]./255, 'MarkerFaceColor',[181,43,33]./255, 'MarkerSize',20);
hold on
plot(0:11,0.1*ones(1,12),'--','Color',[192,192,192]./255,'LineWidth',2);
hold off
xlabel('PLS Component');
ylabel('Explained variance');
set(gca,'XLim',[0,11]);
set(gca,'YLim',[0,0.3],'YTick',0:0.1:0.3);
t1 = text(1,0.28,'*','FontWeight','bold','FontName','Arial','HorizontalAlignment','Center','FontSize',20);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',12);
box off
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[1 1 9 5.6]);
print(gcf,[path,'res\L_PLS_Variance_z_clus1.tif'],'-dtiff','-r1000')

% calculate corrected weight
geneindex=1:size(GENEdata,2);
genes = gene_name;
bootnum=10000;
X=GENEdata;
Y=zscore(MRIdata);
dim=3;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
[R2,p2]=corr(XS(:,3),MRIdata);
if R2(1,1)<0
    stats.W(:,3)=-1*stats.W(:,3);
    XS(:,3)=-1*XS(:,3);
end
[PLS2w,x2] = sort(stats.W(:,3),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);
PLS2_ROIscores_200=XS(:,3);
save([path,'res\L_PLS_ROIscore_z_clus1.mat'],'PLS2_ROIscores_200');
csvwrite([path,'res\L_PLS_ROIscores_z_clus1.csv'],XS(:,3));

PLS2weights = zeros(12639,10000);
parfor i=1:bootnum
    i
    myresample = randsample(size(X,1),size(X,1));
    res(i,:)=myresample; 
    Xr=X(myresample,:); 
    Yr=Y(myresample,:); 
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); 
   
    temp=stats.W(:,3);
    newW=temp(x2); 
    if corr(PLS2w,newW)<0 
        newW=-1*newW;
    end
    PLS2weights(:,i) = newW; 
end
PLS2sw = std(PLS2weights');
temp2=PLS2w./PLS2sw';
[Z2,ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);
fid2 = fopen([path,'res\L_PLS_geneWeights_z_clus1.csv'],'w');
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f\n', PLS2{i},geneindex2(i), Z2(i));
end
fclose(fid2);

% generate weighted image for BrainNet Viewer
vol_new=zeros(size(vol_par));
for i=1:105
    ind=find(vol_par==region_ind(i));
    vol_new(ind)=PLS2_ROIscores_200(i);
end
hdr_new=hdr_par;
hdr_new.fname=[path,'res\L_PLS_weight_z_clus1.nii'];
spm_write_vol(hdr_new,vol_new);

% correlation between PLS score and deviation
V1 = spm_vol('D:\Data_Chen\With_DIDA_all_HC\subtype\baseline\brain\z_clus1.nii');
Y1 = spm_read_vols(V1);
z = zeros(111,1);
for i=1:111
    z(i,1)=mean(Y1(vol_par==i));
end
corr_real = corr(z(region_ind),PLS2_ROIscores_200);

load('D:\Data_Chen\With_DIDA_all_HC\subtype\gene_abagen\res\surrogate_maps_z_clus1.mat');
vol_mask=spm_read_vols(spm_vol([path,'code\shen268_group.nii']));
ind = find(vol_mask);
corr_surr = zeros(1,10000);
z_surr = zeros(111,1);
for j = 1:10000
    disp(j);
    Y_s = zeros(hdr_par.dim);
    Y_s(ind) = surrogate_maps(j,:);
    for i = 1:111
        z_surr(i) = mean(Y_s(vol_par==i));
    end
    corr_surr(1,j) = corr(z_surr(region_ind),PLS2_ROIscores_200);
end
p = length(find(corr_surr(1,:)>corr_real))/10000;
save([path,'res\L_corr_surr_z_clus1.mat'],'corr_surr','corr_real','p');

% Draw figures for correlations
close all
dotcolor = [0.71, 0.17, 0.13];
linecolor = [0.71, 0.17, 0.13];
ylable1 = 'PLS scores';
xlable1 = {'Deviation values'};
[xData, yData] = prepareCurveData(z(region_ind), PLS2_ROIscores_200);
ft = fittype( 'poly1' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf];
opts.Upper = [Inf Inf];
[fitresult, gof] = fit( xData, yData, ft, opts );
h=plot( fitresult, xData, yData);
set(h(1),'Marker','.','MarkerSize',6,'Color',dotcolor)
set(h(2),'LineWidth',0.5,'Color',linecolor)
hold on
xFit = linspace(min(xData),max(xData),100);
yPredict = predint(fitresult,xFit,0.95,'functional','off');
fy = cat(2,yPredict(:,2)',flip(yPredict(:,1),1)')';
fx = cat(2,xFit,flip(xFit',1)')';
fill(fx,fy,[0.71, 0.17, 0.13],'EdgeAlpha',0,'FaceAlpha',0.3);
hold off
legend off
ylabel(ylable1);
xlabel(xlable1);
set(gca,'LineWidth',0.5);
set(gca,'FontName','Arial','FontSize',6);
set(gca,'XTick',-1:0.5:1);
set(gca,'YLim',[-0.2,0.3]);
set(gca,'YTick',[-0.2:0.1:0.3]);
t1 = text(-0.9,0.25,{'{\itr} = 0.50';'{\itP} < 0.0001'},'FontName','Arial','FontSize',6);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf,'Paperposition',[0 0 3.95 3.8]);
grid off
box off
print(gcf,[path,'res\L_corr_PLS_z_clus1.tif'],'-dtiff','-r1000')
