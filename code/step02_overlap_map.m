% To quantify the overall extent of individual deviations, we calculated the number of brain regions
% exhibiting extreme deviations for each participant.
% To assess the intersubject heterogeneity of deviations, we computed a spatial overlap map by 
% determining the percentage of participants with extreme deviation (Z>2.6 or Z<-2.6) in each brain region.

z_base=Z_estimate([1:197],:);
z_follow=Z_estimate([198:244],:);

thr=2.6;
z_sig=zeros(size(z_base));
z_sig(find(z_base<-thr))=-1;
z_sig(find(z_base>thr))=1;

num_sub_dev=sum(z_sig~=0,2);
pro1=sum(num_sub_dev~=0,1)/size(z_sig,1);
num_sub_dev_pos=sum(z_sig==1,2);
pro2=sum(num_sub_dev_pos~=0,1)/size(z_sig,1);
num_sub_dev_neg=sum(z_sig==-1,2);
pro3=sum(num_sub_dev_neg~=0,1)/size(z_sig,1);

num_dev_pos=sum(z_sig==1,1);
pro2=num_dev_pos/size(z_sig,1)*100;
num_vox_dev_pos=size(find(num_dev_pos),2);
pro22=num_vox_dev_pos/220;

num_dev_neg=sum(z_sig==-1,1);
pro3=num_dev_neg/size(z_sig,1)*100;
num_vox_dev_neg=size(find(num_dev_neg),2);
pro33=num_vox_dev_neg/220;

show_z(pro2,'pro_dev_pos.nii')
show_z(pro3,'pro_dev_neg.nii')
