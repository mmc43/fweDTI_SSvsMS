%% Load the data
close all 
clear all

metric = {"dti_SS_FA","dti_SS_MD","fw_MS_FA","fw_MS_FW","fw_MS_MD","fw_SS_FA","fw_SS_FW","fw_SS_MD"};

age=readmatrix('/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/age.txt');

data = zeros(size(age,1),48,size(metric,2));

for m=1:size(metric,2)
    
    char(metric{m})

    csv_filename = ['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data/data_thr07/ROI_analysis/mean_' char(metric{m}) '_JHU_thr07.csv'];
    
    data(:,:,m)=readmatrix(csv_filename);
    
end


%% ROI analyses combined hemipheres  

% Group hemisphere (leave first 6 as separate ROIs)
roi_names = textread('JHU-labels.txt','%s');

roig = {}; Nroi = 0; new_names = {}; 
for r = 1:6 % first 6 are unimodal ROIs
        Nroi = Nroi+1;
        roig{Nroi} = r;
        name = roi_names{r}(1:(end-2));
        name(find(name=='_')) = ' ';
        new_names{Nroi} = name;
end


for r = 1:((size(data,2)-6)/2) % assumes even number

    rind = (r-1)*2 + [0 1] + 6 + 1;

    Nroi = Nroi+1;
    roig{Nroi} = rind;
    name = roi_names{rind(1)}(1:(end-2));
    name(find(name=='_')) = ' ';
    new_names{Nroi} = name;
   
end
Nroi


% Re-average across hemispheres 
d = []; 
for r = 1:Nroi
    
    d=[d; mean(data(:,roig{r},:),2)];
        
end

d=squeeze(d);
d=reshape(d, [size(data,1),Nroi,size(data,3)]);

%%

close all

%load ROI data
%MD_MS_data = data(:,:,5);
%MD_SS_data = data(:,:,8);

FA_MS_sig = zeros(Nroi,1);
FA_SS_sig = zeros(Nroi,1);
FA_MS_SS_int_sig = zeros(Nroi,1);

for r=1:Nroi
    
    %multi-shell MD vs age
    y=d(:,r,5);  
    x=age;
    %x=age;
    
    mdl=fitlm(x,y);
    
    pval_MDm(r) = mdl.Coefficients.pValue(2);
    adjsR2_MDm(r) = mdl.Rsquared.Adjusted;
    coeff_age_MDm(r) = mdl.Coefficients.Estimate(2);
    
    %multi-shell FA vs age
    y=d(:,r,3);  
    x=age;
    
    mdl=fitlm(x,y);
    
    pval_FAm(r) = mdl.Coefficients.pValue(2);
    adjsR2_FAm(r) = mdl.Rsquared.Adjusted;
    coeff_age_FAm(r) = mdl.Coefficients.Estimate(2); 
    
    if mdl.Coefficients.pValue(2)<=0.05
        FA_MS_sig(r)=1;
    end
    
    
    %single-shell MD vs age
    y=d(:,r,8);  
    x=age;

    mdl=fitlm(x,y);
    
    pval_MD(r) = mdl.Coefficients.pValue(2);
    adjsR2_MD(r) = mdl.Rsquared.Adjusted;
    coeff_age_MD(r) = mdl.Coefficients.Estimate(2);
    
    %single-shell FA vs age
    y=d(:,r,6);  
    x=age;
    
    mdl=fitlm(x,y);
    
    pval_FA(r) = mdl.Coefficients.pValue(2);
    adjsR2_FA(r) = mdl.Rsquared.Adjusted;
    coeff_age_FA(r) = mdl.Coefficients.Estimate(2);  
    
    if mdl.Coefficients.pValue(2)<=0.05
        FA_SS_sig(r)=1;
    end
    
        
    
    %ANCOVA MD analysis
    
    yms=d(:,r,5);  
    xms=age;
    
    yss=d(:,r,8);  
    xss=age; 
    
    ancova_x = [xms; xss];
    ancova_y = [yms; yss];
    
    group_ancova = [ones(size(yms)); 2*ones(size(yss))];
    
    [h_out,a_out,c_out,s_out] = aoctool(ancova_x, ancova_y, group_ancova, 0.05,'age','MD','fitting','off','separate lines');
       
    inter_pval(r) = a_out{4,6};
    slope_ms(r) = c_out{5,2} + c_out{6,2};
    slope_ss(r) = c_out{5,2} + c_out{7,2};
    
    
    %ANCOVA FA analysis
    
    yms=d(:,r,3);      
    yss=d(:,r,6);  

    ancova_x = [age; age];
    ancova_y = [yms; yss];
    
    group_ancova = [ones(size(yms)); 2*ones(size(yss))];
    
    [h_out,a_out,c_out,s_out] = aoctool(ancova_x, ancova_y, group_ancova, 0.05,'age','FA','fitting','off','separate lines');
       
    interfa_pval(r) = a_out{4,6};
    slopefa_ms(r) = c_out{5,2} + c_out{6,2};
    slopefa_ss(r) = c_out{5,2} + c_out{7,2};
   
 
end

%correction for multiple comparisons

[h_FDR_MD, crit_p_MD]=fdr_bh(pval_MD, 0.05, 'pdep');
[h_FDR_MDm, crit_p_MDm]=fdr_bh(pval_MDm, 0.05, 'pdep');
[h_FDR_FA, crit_p_FA]=fdr_bh(pval_FA, 0.05, 'pdep');
[h_FDR_FAm, crit_p_FAm]=fdr_bh(pval_FAm, 0.05, 'pdep');


R2all=[adjsR2_MDm' adjsR2_MD' adjsR2_FAm' adjsR2_FA'];
[~,rind] = sort(mean(R2all,2),'descend');

lab={'fwe MS MD','fwe SS MD','fwe MS FA','fwe SS FA'};

figure(1); clf
imagesc(R2all(rind,:)),colorbar,colormap('jet')
set(gca,'Xtick', [1:size(lab,2)], 'XTickLabels',lab,'Ytick',[1:Nroi],'YTickLabels',{new_names{rind}})
title('Adjusted R2');


pvalall=[pval_MDm' pval_MD' pval_FAm' pval_FA'];
[~,rind] = sort(mean(pvalall,2),'ascend');

r2yfullcolourmap = [ones(255,1) linspace(1,0,255)' zeros(255,1); 0 0 0];

figure(2); clf
imagesc(pvalall(rind,:)),colorbar,colormap(r2yfullcolourmap)
set(gca,'Xtick', [1:size(lab,2)], 'XTickLabels',lab,'Ytick',[1:Nroi],'YTickLabels',{new_names{rind}})
title('p-value (age effects)');
caxis([0 0.05])

H_all=[h_FDR_MDm' h_FDR_MD' h_FDR_FAm' h_FDR_FA'];
[~,rind] = sort(mean(pvalall,2),'ascend');

figure(3); clf
imagesc(H_all(rind,:)),colorbar,colormap(r2yfullcolourmap)
set(gca,'Xtick', [1:size(lab,2)], 'XTickLabels',lab,'Ytick',[1:Nroi],'YTickLabels',{new_names{rind}})
title('FDR corrected significant tests');
caxis([0 1.1])

%define new colour map to match tbss results
r2ycolourmap = [ones(127,1) linspace(0,1,127)' zeros(127,1)];
b2lbcolourmap = [zeros(127,1) linspace(220/255,0,127)' ones(127,1)];
mycolourmap = [b2lbcolourmap; 0 0 0; 0 0 0; r2ycolourmap];

%slopemd=[slope_ms' slope_ss'];
%slopefa=[slopefa_ms' slopefa_ss'];

slopemd=[(slope_ss.*h_FDR_MD)' (slope_ms.*h_FDR_MDm)' ];
slopefa=[(slopefa_ss.*h_FDR_FA)' (slopefa_ms.*h_FDR_FAm)' ];

[~,rind] = sort(mean(slopemd,2),'descend');
%[~,rind] = sort(mean(slopefa,2),'descend');

labmd={'fwe MD SS','fwe MD MS'};

figure(4); clf
imagesc(slopemd(rind,:)),colorbar,colormap(mycolourmap)
set(gca,'Xtick', [1:size(labmd,2)], 'XTickLabels',labmd,'Ytick',[1:Nroi],'YTickLabels',{new_names{rind}})
title('rate of MD change with age (mm^2s^-^1y^-^1)       ');
caxis([-2*10^(-6) 2*10^(-6)])

labfa={'fwe FA SS','fwe FA MS'};

figure(5); clf
imagesc(slopefa(rind,:)),colorbar,colormap(mycolourmap)
set(gca,'Xtick', [1:size(labfa,2)], 'XTickLabels',labfa,'Ytick',[1:Nroi],'YTickLabels',{new_names{rind}})
title('rate of FA change with age (y^-^1)');
caxis([-min(abs(min(min(slopefa))), max(max(slopefa))) min(abs(min(min(slopefa))), max(max(slopefa)))]);

%plot cases with significant interaction betweeen FA MS and FA SS AND FA MS slope < FA SS slope  

compare_slopes=slopefa_ms < slopefa_ss;
[h_FDR_inter, crit_p_MD]=fdr_bh(inter_pval, 0.05, 'pdep');

greencolourmap= [ zeros(255,1) zeros(255,1) zeros(255,1); 0 1 0];


figure(6); clf
imagesc([h_FDR_inter(rind)' compare_slopes(rind)']),colormap(greencolourmap)
set(gca,'Xtick', [1:2], 'XTickLabels',{'significant interation', 'slope fwe FA MS < slope fwe FA SS'},'Ytick',[1:Nroi],'YTickLabels',{new_names{rind}})
title({'ANCOVA: Interaction between', 'FWE fitting method and Age'});
caxis([0 1])

figure(7); clf
scatter(slopemd(:,2),slopefa(:,1)-slopefa(:,2));
title({'Relationship between fwe MD MS and FA slopes', ''})
xlabel('slope fwe MD MS (mm^2s^-^1y^-^1)')
ylabel('slope fwe FA SS - slope fwe FA MS (y^-^1)')


lmfit_filename='fw_MS_MD_vs_age.csv';
writematrix([pval_MD; adjsR2_MD; coeff_age_MD], lmfit_filename);

lmfit_filename='fw_MS_FA_vs_age.csv';
writematrix([pval_FAm; adjsR2_FAm; coeff_age_FAm], lmfit_filename);

lmfit_filename='fw_SS_FA_vs_age.csv';
writematrix([pval_FA; adjsR2_FA; coeff_age_FA], lmfit_filename);

lmfit_filename='fw_MD_vs_age_ancova.csv';
writematrix([inter_pval; slope_ms; slope_ss], lmfit_filename);

lmfit_filename='fw_FA_vs_age_ancova.csv';
writematrix([interfa_pval; slopefa_ms; slopefa_ss], lmfit_filename);

%create significant ROI maps

JHU_rois=niftiread('/imaging/local/software/fsl/v5.0.11/x86_64/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz');
JHU_info=niftiinfo('/imaging/local/software/fsl/v5.0.11/x86_64/data/atlases/JHU/JHU-ICBM-labels-1mm.nii.gz');

JHU_rois_sig_FA_MS = cast(JHU_rois-JHU_rois, 'double');
JHU_rois_sig_FA_SS = cast(JHU_rois-JHU_rois, 'double');
%JHU_rois_sig_FA_MS_SS_int = JHU_rois-JHU_rois;
JHU_rois_negsig_FA_MS = cast(JHU_rois-JHU_rois, 'double');
JHU_rois_negsig_FA_SS = cast(JHU_rois-JHU_rois, 'double');

JHU_rois_sig_MD_MS = cast(JHU_rois-JHU_rois, 'double');
JHU_rois_sig_MD_SS = cast(JHU_rois-JHU_rois, 'double');
JHU_rois_negsig_MD_MS = cast(JHU_rois-JHU_rois, 'double');
JHU_rois_negsig_MD_SS = cast(JHU_rois-JHU_rois, 'double');

JHU_info.Datatype='double';


%generate nifti files showing significant ROIs for the single hemisphere
%ROIs first (ROIs 1 to 6)
for r=1:6
    
   if pval_FAm(r)<=0.05 && coeff_age_FAm(r)>0
       JHU_rois_sig_FA_MS(JHU_rois==r)=1-pval_FAm(r);
   end

   if pval_FA(r)<=0.05 && coeff_age_FA(r)>0
    
       JHU_rois_sig_FA_SS(JHU_rois==r)=1-pval_FA(r);
   end
       
   if pval_FAm(r)<=0.05 && coeff_age_FAm(r)<0
       JHU_rois_negsig_FA_MS(JHU_rois==r)=1-pval_FAm(r);
   end

   if pval_FA(r)<=0.05 && coeff_age_FA(r)<0
    
       JHU_rois_negsig_FA_SS(JHU_rois==r)=1-pval_FA(r);
   end
   
   %repeat for MD rois
   
   if pval_MDm(r)<=0.05 && coeff_age_MDm(r)>0
       JHU_rois_sig_MD_MS(JHU_rois==r)=1-pval_MDm(r);
   end

   if pval_MD(r)<=0.05 && coeff_age_MD(r)>0
    
       JHU_rois_sig_MD_SS(JHU_rois==r)=1-pval_MD(r);
   end
       
   if pval_MDm(r)<=0.05 && coeff_age_MDm(r)<0
       JHU_rois_negsig_MD_MS(JHU_rois==r)=1-pval_MDm(r);
   end

   if pval_MD(r)<=0.05 && coeff_age_MD(r)<0
    
       JHU_rois_negsig_MD_SS(JHU_rois==r)=1-pval_MD(r);
   end
   
end

%repeat for the merged (right + left) ROIs
for r=7:Nroi
    
   if pval_FAm(r)<=0.05 && coeff_age_FAm(r)>0
       JHU_rois_sig_FA_MS(JHU_rois==roig{r}(1,1)) = 1-pval_FAm(r); 
       JHU_rois_sig_FA_MS(JHU_rois==roig{r}(1,2)) = 1-pval_FAm(r);
   end

   if pval_FA(r)<=0.05 && coeff_age_FA(r)>0
       JHU_rois_sig_FA_SS(JHU_rois==roig{r}(1,1)) = 1-pval_FA(r);
       JHU_rois_sig_FA_SS(JHU_rois==roig{r}(1,2)) = 1-pval_FA(r); 
   end
       
   if pval_FAm(r)<=0.05 && coeff_age_FAm(r)<0
       JHU_rois_negsig_FA_MS(JHU_rois==roig{r}(1,1))= 1-pval_FAm(r);
       JHU_rois_negsig_FA_MS(JHU_rois==roig{r}(1,2))= 1-pval_FAm(r);
   end

   if pval_FA(r)<=0.05 && coeff_age_FA(r)<0
       JHU_rois_negsig_FA_SS(JHU_rois==roig{r}(1,1)) = 1-pval_FA(r);
       JHU_rois_negsig_FA_SS(JHU_rois==roig{r}(1,2)) = 1-pval_FA(r);
   end
   
   
   %repeat for MD rois
   
   if pval_MDm(r)<=0.05 && coeff_age_MDm(r)>0
       JHU_rois_sig_MD_MS(JHU_rois==roig{r}(1,1)) = 1-pval_MDm(r); 
       JHU_rois_sig_MD_MS(JHU_rois==roig{r}(1,2)) = 1-pval_MDm(r);
   end

   if pval_MD(r)<=0.05 && coeff_age_MD(r)>0
       JHU_rois_sig_MD_SS(JHU_rois==roig{r}(1,1)) = 1-pval_MD(r);
       JHU_rois_sig_MD_SS(JHU_rois==roig{r}(1,2)) = 1-pval_MD(r); 
   end
       
   if pval_MDm(r)<=0.05 && coeff_age_MDm(r)<0
       JHU_rois_negsig_MD_MS(JHU_rois==roig{r}(1,1))= 1-pval_MDm(r);
       JHU_rois_negsig_MD_MS(JHU_rois==roig{r}(1,2))= 1-pval_MDm(r);
   end

   if pval_MD(r)<=0.05 && coeff_age_MD(r)<0
       JHU_rois_negsig_MD_SS(JHU_rois==roig{r}(1,1)) = 1-pval_MD(r);
       JHU_rois_negsig_MD_SS(JHU_rois==roig{r}(1,2)) = 1-pval_MD(r);
   end
   
end

niftiwrite(JHU_rois_sig_FA_MS,'JHU_rois_sig_FA_MS_map_merged.nii', JHU_info);
niftiwrite(JHU_rois_sig_FA_SS,'JHU_rois_sig_FA_SS_map_merged.nii', JHU_info);
niftiwrite(JHU_rois_negsig_FA_MS,'JHU_rois_negsig_FA_MS_map_merged.nii', JHU_info);
niftiwrite(JHU_rois_negsig_FA_SS,'JHU_rois_negsig_FA_SS_map_merged.nii', JHU_info);

niftiwrite(JHU_rois_sig_MD_MS,'JHU_rois_sig_MD_MS_map_merged.nii', JHU_info);
niftiwrite(JHU_rois_sig_MD_SS,'JHU_rois_sig_MD_SS_map_merged.nii', JHU_info);
niftiwrite(JHU_rois_negsig_MD_MS,'JHU_rois_negsig_MD_MS_map_merged.nii', JHU_info);
niftiwrite(JHU_rois_negsig_MD_SS,'JHU_rois_negsig_MD_SS_map_merged.nii', JHU_info);

