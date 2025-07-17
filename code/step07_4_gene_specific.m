
clear;
load('D:\Data_Chen\With_DIDA_all_HC\subtype\clus_base.mat');
load('D:\Data_Chen\With_DIDA_all_HC\res_norm\res_norm\z_base_unique.mat');
img_data=mean(z_base(ind_clus1,110:220),1);

geneset = readtable('D:\Data_Chen\With_DIDA_all_HC\subtype\gene_abagen\res\L_PLS3_geneWeights_z_clus1.csv', 'TextType', 'string');
ind=find(geneset{:,4}<0.0132);
geneset_sig=cellstr(geneset{ind,1});

expression_all =  importdata('D:\Data_Chen\With_DIDA_all_HC\subtype\gene_abagen\res\L_expression_filter.csv');
expression = expression_all.data;
gene_name = expression_all.textdata;
missingdata_regions=find(isnan(expression(:,2)));
region_ind=setdiff(expression(:,1),missingdata_regions);
group_express=expression(region_ind,2:end);
gene_name = gene_name(2:end);

expressions=group_express;
img_data=img_data(region_ind)';

[I J]=find(isnan(expressions));
for i=1:length(I)
    expressions(I(i),J(i))=0;
end

res_nullcoexpGene = permutation_null_coexp(img_data, geneset_sig, expressions, gene_name)
res_nullbraingene = permutation_null_brain(img_data, geneset_sig, expressions, gene_name)