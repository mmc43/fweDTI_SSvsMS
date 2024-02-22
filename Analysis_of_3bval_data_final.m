%% Load all data

subj = {"CBU220739","CBU220740","CBU220741"};

metric = {"FA","MD","FW"};

metric_filename = {"FA_fw_thresholded_0.7.nii.gz", "MD_fw_thresholded_0.7.nii.gz", "FW.nii"};

dti_FA_rois = zeros(3,2,48);

FA_MS_12_rois = zeros(3,2,48);
FA_MS_31_rois = zeros(3,2,48);
FA_SS_1_rois = zeros(3,2,48);

MD_MS_12_rois = zeros(3,2,48);
MD_MS_31_rois = zeros(3,2,48);
MD_SS_1_rois = zeros(3,2,48);

FW_MS_12_rois = zeros(3,2,48);
FW_MS_31_rois = zeros(3,2,48);
FW_SS_1_rois = zeros(3,2,48);

diff_FA_MS_12_GS_rois = zeros(3,2,48);
diff_FA_SS_1_GS_rois = zeros(3,2,48);

FA_MS_12_GS_rois_diff = zeros(3,2,48);
FA_SS_1_GS_rois_diff = zeros(3,2,48);


for s=1:3
    
    for r=1:2
    
        JHU_rois = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/atlases/JHU-ICBM-labels.nii.gz']);
        
        %load FA data
        m=1;
        FA_MS_12 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_2000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_multi/' char(metric_filename{m})]);
        FA_MS_31 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_300_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_multi/' char(metric_filename{m})]);
        FA_SS_1 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_single/' char(metric_filename{m})]);
    
        %load MD data
        m=2;
        MD_MS_12 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_2000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_multi/' char(metric_filename{m})]);
        MD_MS_31 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_300_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_multi/' char(metric_filename{m})]);
        MD_SS_1 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_single/' char(metric_filename{m})]);

        %load FW data
        m=3;
        FW_MS_12 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_2000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_multi/' char(metric_filename{m})]);
        FW_MS_31 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_300_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_multi/' char(metric_filename{m})]);
        FW_SS_1 = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/fw_single/FW.nii.gz']);
              
        %load dti FA
        dti_FA = niftiread(['/imaging/correia/users/mc04/MRI_METHODS/FreeWaterDiffusion/data_3bvals_new/data_subsets/' char(subj{s}) '/b_1000_run' num2str(r) '/proc/nipype/datasink/dwi_space/dti/dtifit_wls_FA.nii.gz']);

        
        
        % find mean for each ROI    
        
        for roi=1:48
            
            dti_FA_rois(s,r,roi) = mean(dti_FA(JHU_rois==roi));
            
            FA_MS_12_rois(s,r,roi)=mean(FA_MS_12(JHU_rois==roi));
            FA_MS_31_rois(s,r,roi)=mean(FA_MS_31(JHU_rois==roi));
            FA_SS_1_rois(s,r,roi)=mean(FA_SS_1(JHU_rois==roi));
            
            MD_MS_12_rois(s,r,roi)=mean(MD_MS_12(JHU_rois==roi & MD_MS_12<0.001));
            MD_MS_31_rois(s,r,roi)=mean(MD_MS_31(JHU_rois==roi & MD_MS_31<0.001));
            MD_SS_1_rois(s,r,roi)=mean(MD_SS_1(JHU_rois==roi & MD_SS_1<0.001));
            
            FW_MS_12_rois(s,r,roi)=mean(FW_MS_12(JHU_rois==roi));
            FW_MS_31_rois(s,r,roi)=mean(FW_MS_31(JHU_rois==roi));
            FW_SS_1_rois(s,r,roi)=mean(FW_SS_1(JHU_rois==roi));
        end 
        
        %calculate differences between FA estimates
        
        diff_FA_MS_12_GS = FA_MS_12-FA_MS_31;
        diff_FA_SS_1_GS = FA_SS_1-FA_MS_31;
        
       %Find mean difference per ROI
       
       for roi=1:48
       
           diff_FA_MS_12_GS_rois(s,r,roi) = mean(diff_FA_MS_12_GS(JHU_rois==roi));
           diff_FA_SS_1_GS_rois(s,r,roi) = mean(diff_FA_SS_1_GS(JHU_rois==roi));

           FA_MS_12_GS_rois_diff(s,r,roi) = FA_MS_12_rois(s,r,roi)-FA_MS_31_rois(s,r,roi);
           FA_SS_1_GS_rois_diff(s,r,roi) = FA_SS_1_rois(s,r,roi)-FA_MS_31_rois(s,r,roi);
       
       end
           
    end

    for roi=1:48
        
        coeff_var_FA_MS_12(s,roi)= std([FA_MS_12_rois(s,1,roi) FA_MS_12_rois(s,2,roi)])/mean([FA_MS_12_rois(s,1,roi) FA_MS_12_rois(s,2,roi)]);
        coeff_var_FA_SS_1(s,roi)= std([FA_SS_1_rois(s,1,roi) FA_SS_1_rois(s,2,roi)])/mean([FA_SS_1_rois(s,1,roi) FA_SS_1_rois(s,2,roi)]);
        coeff_var_MD_MS_12(s,roi)= std([MD_MS_12_rois(s,1,roi) MD_MS_12_rois(s,2,roi)])/mean([MD_MS_12_rois(s,1,roi) MD_MS_12_rois(s,2,roi)]);
        coeff_var_MD_SS_1(s,roi)= std([MD_SS_1_rois(s,1,roi) MD_SS_1_rois(s,2,roi)])/mean([MD_SS_1_rois(s,1,roi) MD_SS_1_rois(s,2,roi)]);
        %coeff_var_FW_MS_12(s,roi)= std([FW_MS_12_rois(s,1,roi) FW_MS_12_rois(s,2,roi)])/mean([FW_MS_12_rois(s,1,roi) FW_MS_12_rois(s,2,roi)]);
        %coeff_var_FW_SS_1(s,roi)= std([FW_SS_1_rois(s,1,roi) FW_SS_1_rois(s,2,roi)])/mean([FW_SS_1_rois(s,1,roi) FW_SS_1_rois(s,2,roi)]);
    
    end

end

%% Create plots for paper

close all

% Figure 1 FA MS vs GS, all subjects run 1
figure,scatter(FA_MS_31_rois(1,1,:),FA_MS_12_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FA_MS_31_rois(2,1,:),FA_MS_12_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FA_MS_31_rois(3,1,:),FA_MS_12_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
line(0:0.1:1,0:0.1:1,'Color', 'black');

title('fwe FA MS vs gold standard FA')
xlabel('gold standard FA')
ylabel('fwe FA MS (b=1000,2000)')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast') 

% Figure 2 FA MS vs GS, all subjects both runs
figure,scatter(FA_MS_31_rois(1,1,:),FA_MS_12_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FA_MS_31_rois(2,1,:),FA_MS_12_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FA_MS_31_rois(3,1,:),FA_MS_12_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 

%add run 2
scatter(FA_MS_31_rois(1,2,:),FA_MS_12_rois(1,2,:),'x','MarkerEdgeColor', [1 0.5 0.5]);
scatter(FA_MS_31_rois(2,2,:),FA_MS_12_rois(2,2,:),'x','MarkerEdgeColor', [0.5 1 0.5]);  
scatter(FA_MS_31_rois(3,2,:),FA_MS_12_rois(3,2,:),'x','MarkerEdgeColor', [0.5 0.5 1]); 

line(0:0.1:1,0:0.1:1,'Color', 'black');

title('fwe FA MS vs gold standard FA')
xlabel('gold standard FA')
ylabel('fwe FA MS (b=1000,2000)')
legend({'Subj 1, run 1','Subj 2, run 1','Subj 3, run 1','Subj 1, run 2','Subj 2, run 2','Subj 3, run 2','x=y'},'Location','Southeast')

% Figure 3 FA SS vs GS, all subjects run 1
figure,scatter(FA_MS_31_rois(1,1,:),FA_SS_1_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FA_MS_31_rois(2,1,:),FA_SS_1_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FA_MS_31_rois(3,1,:),FA_SS_1_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
line(0:0.1:1,0:0.1:1,'Color', 'black');

title('fwe FA SS vs gold standard FA')
xlabel('gold standard FA')
ylabel('fwe FA SS (b=1000)')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast')


% Figure 4 FA SS vs GS, all subjects both runs
figure,scatter(FA_MS_31_rois(1,1,:),FA_SS_1_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FA_MS_31_rois(2,1,:),FA_SS_1_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FA_MS_31_rois(3,1,:),FA_SS_1_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 

%add run 2
scatter(FA_MS_31_rois(1,2,:),FA_SS_1_rois(1,2,:),'x','MarkerEdgeColor', [1 0.5 0.5]);
scatter(FA_MS_31_rois(2,2,:),FA_SS_1_rois(2,2,:),'x','MarkerEdgeColor', [0.5 1 0.5]);  
scatter(FA_MS_31_rois(3,2,:),FA_SS_1_rois(3,2,:),'x','MarkerEdgeColor', [0.5 0.5 1]); 

line(0:0.1:1,0:0.1:1,'Color', 'black');

title('fwe FA SS vs gold standard FA')
xlabel('gold standard FA')
ylabel('fwe FA SS (b=1000)')
legend({'Subj 1, run 1','Subj 2, run 1','Subj 3, run 1','Subj 1, run 2','Subj 2, run 2','Subj 3, run 2','x=y'},'Location','Southeast')

 
% Figure 5 MD MS vs GS, all subjects run 1
figure,scatter(MD_MS_31_rois(1,1,:),MD_MS_12_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),MD_MS_12_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),MD_MS_12_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
line(0:0.00001:0.0007,0:0.00001:0.0007,'Color', 'black');

title('fwe MD MS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe MD MS (b=1000,2000)')
%legend({'Subj 1, run 1','Subj 2, run 1','Subj 3, run 1','Subj 1, run 2','Subj 2, run 2','Subj 3, run 2','x=y'},'Location','Southeast')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast')

% Figure 6 MD MS vs GS, all subjects both runs
figure,scatter(MD_MS_31_rois(1,1,:),MD_MS_12_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),MD_MS_12_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),MD_MS_12_rois(3,1,:),'MarkerEdgeColor', [0 0 1]);

%add run 2
scatter(MD_MS_31_rois(1,2,:),MD_MS_12_rois(1,2,:),'x', 'MarkerEdgeColor', [1 0.5 0.5]);
scatter(MD_MS_31_rois(2,2,:),MD_MS_12_rois(2,2,:),'x', 'MarkerEdgeColor', [0.5 1 0.5]);  
scatter(MD_MS_31_rois(3,2,:),MD_MS_12_rois(3,2,:),'x', 'MarkerEdgeColor', [0.5 0.5 1]);

line(0:0.00001:0.0007,0:0.00001:0.0007,'Color', 'black');

title('fwe MD MS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe MD MS (b=1000,2000)')
legend({'Subj 1, run 1','Subj 2, run 1','Subj 3, run 1','Subj 1, run 2','Subj 2, run 2','Subj 3, run 2','x=y'},'Location','Southeast')



% Figure 7 MD SS vs GS, all subjects run 1
figure,scatter(MD_MS_31_rois(1,1,:),MD_SS_1_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),MD_SS_1_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),MD_SS_1_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
line(0:0.00001:0.0007,0:0.00001:0.0007,'Color', 'black');

title('fwe MD SS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe MD SS (b=1000)')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast')

% Figure 8 MD SS vs GS, all subjects both runs
figure,scatter(MD_MS_31_rois(1,1,:),MD_SS_1_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),MD_SS_1_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),MD_SS_1_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 

%add run 2
scatter(MD_MS_31_rois(1,2,:),MD_SS_1_rois(1,2,:),'x','MarkerEdgeColor', [1 0.5 0.5]);
scatter(MD_MS_31_rois(2,2,:),MD_SS_1_rois(2,2,:),'x','MarkerEdgeColor', [0.5 1 0.5]);  
scatter(MD_MS_31_rois(3,2,:),MD_SS_1_rois(3,2,:),'x','MarkerEdgeColor', [0.5 0.5 1]);

line(0:0.00001:0.0007,0:0.00001:0.0007,'Color', 'black');

title('fwe MD SS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe MD SS (b=1000)')
legend({'Subj 1, run 1','Subj 2, run 1','Subj 3, run 1','Subj 1, run 2','Subj 2, run 2','Subj 3, run 2','x=y'},'Location','Southeast')


% Figure 9 FA MS - FA GS vs MD GS all subjects run 1
figure,scatter(MD_MS_31_rois(1,1,:),diff_FA_MS_12_GS_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),diff_FA_MS_12_GS_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),diff_FA_MS_12_GS_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
axis([2*10^-4 8*10^-4 -0.2 0.1]);

x=[squeeze(MD_MS_31_rois(1,1,:));squeeze(MD_MS_31_rois(2,1,:));squeeze((MD_MS_31_rois(3,1,:)))];
y=[squeeze(diff_FA_MS_12_GS_rois(1,1,:));squeeze(diff_FA_MS_12_GS_rois(3,1,:));squeeze(diff_FA_MS_12_GS_rois(3,1,:))];

coefficients = polyfit(x, y, 1);

% Create a new x axis with exactly 1000 points.
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

%plot(xFit, yFit, 'k-', 'LineWidth', 1); % Plot fitted line.

title('\DeltaFA MS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe FA MS (b=1000,2000) - gold standard FA')
legend({'Participant 1','Participant 2','Participant 3', 'line of best fit'},'Location','Southeast')

% Figure 10 FA MS - FA GS vs MD GS all subjects both runs
figure,scatter(MD_MS_31_rois(1,1,:),diff_FA_MS_12_GS_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),diff_FA_MS_12_GS_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),diff_FA_MS_12_GS_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 

%add run 2
scatter(MD_MS_31_rois(1,2,:),diff_FA_MS_12_GS_rois(1,2,:),'x','MarkerEdgeColor', [1 0.5 0.5]);
scatter(MD_MS_31_rois(2,2,:),diff_FA_MS_12_GS_rois(2,2,:),'x','MarkerEdgeColor', [0.5 1 0.5]);  
scatter(MD_MS_31_rois(3,2,:),diff_FA_MS_12_GS_rois(3,2,:),'x','MarkerEdgeColor', [0.5 0.5 1]); 

title('fwe FA MS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe FA MS (b=1000,2000) - gold standard FA')
legend({'Subj 1, run 1','Subj 2, run 1','Subj 3, run 1','Subj 1, run 2','Subj 2, run 2','Subj 3, run 2'},'Location','Southeast')


% Figure 11 FA SS - FA GS vs MD GS all subjects run 1
figure,scatter(MD_MS_31_rois(1,1,:),diff_FA_SS_1_GS_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),diff_FA_SS_1_GS_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),diff_FA_SS_1_GS_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
axis([2*10^-4 8*10^-4 -0.2 0.1]);

x=[squeeze(MD_MS_31_rois(1,1,:));squeeze(MD_MS_31_rois(2,1,:));squeeze((MD_MS_31_rois(3,1,:)))];
y=[squeeze(diff_FA_SS_1_GS_rois(1,1,:));squeeze(diff_FA_SS_1_GS_rois(3,1,:));squeeze(diff_FA_SS_1_GS_rois(3,1,:))];

coefficients = polyfit(x, y, 1);

% Create a new x axis with exactly 1000 points.
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

%plot(xFit, yFit, 'k-', 'LineWidth', 1); % Plot fitted line.

title('\DeltaFA SS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe FA SS (b=1000) - gold standard FA')
legend({'Participant 1','Participant 2','Participant 3'},'Location','Southeast')


% Figure 12 FA SS - FA GS vs MD GS all subjects both runs
figure,scatter(MD_MS_31_rois(1,1,:),diff_FA_SS_1_GS_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),diff_FA_SS_1_GS_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),diff_FA_SS_1_GS_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 

%add run 2
scatter(MD_MS_31_rois(1,2,:),diff_FA_SS_1_GS_rois(1,2,:),'x','MarkerEdgeColor', [1 0.5 0.5]);
scatter(MD_MS_31_rois(2,2,:),diff_FA_SS_1_GS_rois(2,2,:),'x','MarkerEdgeColor', [0.5 1 0.5]);  
scatter(MD_MS_31_rois(3,2,:),diff_FA_SS_1_GS_rois(3,2,:),'x','MarkerEdgeColor', [0.5 0.5 1]); 

title('fwe FA SS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe FA SS (b=1000) - gold standard FA')
legend({'Subj 1, run 1','Subj 2, run 1','Subj 3, run 1','Subj 1, run 2','Subj 2, run 2','Subj 3, run 2'},'Location','Southeast')

% Figure 13 FW MS vs GS, all subjects run 1
figure,scatter(FW_MS_31_rois(1,1,:),FW_MS_12_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FW_MS_31_rois(2,1,:),FW_MS_12_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FW_MS_31_rois(3,1,:),FW_MS_12_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
line(0:0.1:1,0:0.1:1,'Color', 'black');

title('fwe FW MS vs gold standard FW index')
xlabel('gold standard FW index')
ylabel('fwe FW MS (b=1000,2000)')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast') 

% Figure 14 FW SS vs GS, all subjects run 1
figure,scatter(FW_SS_1_rois(1,1,:),FW_MS_12_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FW_SS_1_rois(2,1,:),FW_MS_12_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FW_SS_1_rois(3,1,:),FW_MS_12_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
line(0:0.1:1,0:0.1:1,'Color', 'black');

title('fwe FW SS vs gold standard FW index')
xlabel('gold standard FW index')
ylabel('fwe FW SS (b=1000)')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast') 

%EXtra plot

% Figure 15 , all subjects run 1
figure,scatter(MD_MS_31_rois(1,1,:),MD_SS_1_rois(1,1,:)-MD_MS_31_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),MD_SS_1_rois(2,1,:)-MD_MS_31_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),MD_SS_1_rois(3,1,:)-MD_MS_31_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
axis([2*10^-4 8*10^-4  -2.2*10^-4 2*10^-4]);

x=[squeeze(MD_MS_31_rois(1,1,:));squeeze(MD_MS_31_rois(2,1,:));squeeze((MD_MS_31_rois(3,1,:)))];
y=[squeeze(MD_SS_1_rois(1,1,:)-MD_MS_31_rois(1,1,:));squeeze(MD_SS_1_rois(2,1,:)-MD_MS_31_rois(2,1,:));squeeze(MD_SS_1_rois(3,1,:)-MD_MS_31_rois(3,1,:))];

coefficients = polyfit(x, y, 1);

% Create a new x axis with exactly 1000 points.
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

%plot(xFit, yFit, 'k-', 'LineWidth', 1); % Plot fitted line.

title('\DeltaMD SS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe MD SS (b=1000) - gold standard MD')
legend({'Participant 1','Participant 2','Participant 3', 'line of best fit'},'Location','Northwest') 

% Figure 16 , all subjects run 1
figure,scatter(MD_MS_31_rois(1,1,:),MD_MS_12_rois(1,1,:)-MD_MS_31_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),MD_MS_12_rois(2,1,:)-MD_MS_31_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),MD_MS_12_rois(3,1,:)-MD_MS_31_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
axis([2*10^-4 8*10^-4  -2.2*10^-4 2*10^-4]);

x=[squeeze(MD_MS_31_rois(1,1,:));squeeze(MD_MS_31_rois(2,1,:));squeeze((MD_MS_31_rois(3,1,:)))];
y=[squeeze(MD_MS_12_rois(1,1,:)-MD_MS_31_rois(1,1,:));squeeze(MD_MS_12_rois(2,1,:)-MD_MS_31_rois(2,1,:));squeeze(MD_MS_12_rois(3,1,:)-MD_MS_31_rois(3,1,:))];

coefficients = polyfit(x, y, 1);

% Create a new x axis with exactly 1000 points.
xFit = linspace(min(x), max(x), 1000);
% Get the estimated yFit value for each of those 1000 new x locations.
yFit = polyval(coefficients , xFit);

%plot(xFit, yFit, 'k-', 'LineWidth', 1); % Plot fitted line.

title('\DeltaMD MS vs gold standard MD')
xlabel('gold standard MD')
ylabel('fwe MD MS (b=1000,2000) - gold standard MD')
legend({'Participant 1','Participant 2','Participant 3', 'line of best fit'},'Location','Northwest') 


% Figure 17 FA dti vs FA gold standard all subjects run 1
figure,scatter(FA_MS_31_rois(1,1,:),dti_FA_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FA_MS_31_rois(2,1,:),dti_FA_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FA_MS_31_rois(3,1,:),dti_FA_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 
line(0:0.1:1,0:0.1:1,'Color', 'black');

title('dti FA SS vs gold standard FA')
xlabel('gold standard FA')
ylabel('dti FA SS (b=1000)')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast') 

% Figure 18 Delta FA dti vs FA gold standard all subjects run 1
figure,scatter(FA_MS_31_rois(1,1,:),dti_FA_rois(1,1,:)-FA_MS_31_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(FA_MS_31_rois(2,1,:),dti_FA_rois(2,1,:)-FA_MS_31_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(FA_MS_31_rois(3,1,:),dti_FA_rois(3,1,:)-FA_MS_31_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 

title('DTI \Delta FA vs gold standard FA')
xlabel('gold standard FA')
ylabel('dti FA SS (b=1000) - gold standard FA')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast') 

% Figure 19 Delta FA dti vs MD gold standard all subjects run 1
figure,scatter(MD_MS_31_rois(1,1,:),dti_FA_rois(1,1,:)-FA_MS_31_rois(1,1,:),'MarkerEdgeColor', [1 0 0]);
hold on
scatter(MD_MS_31_rois(2,1,:),dti_FA_rois(2,1,:)-FA_MS_31_rois(2,1,:),'MarkerEdgeColor', [0 1 0]);  
scatter(MD_MS_31_rois(3,1,:),dti_FA_rois(3,1,:)-FA_MS_31_rois(3,1,:),'MarkerEdgeColor', [0 0 1]); 

title('\Delta dti FA SS vs gold standard MD')
xlabel('gold standard MD')
ylabel('dti FA SS (b=1000) - gold standard FA')
legend({'Participant 1','Participant 2','Participant 3','x=y'},'Location','Southeast') 




% %plot coeff of variation
% 
% figure,boxplot([coeff_var_FA_MS_12(1,:)' coeff_var_FA_SS_1(1,:)' coeff_var_FA_MS_12(2,:)' coeff_var_FA_SS_1(2,:)' coeff_var_FA_MS_12(3,:)' coeff_var_FA_SS_1(3,:)']);
% figure,boxplot([coeff_var_MD_MS_12(1,:)' coeff_var_MD_SS_1(1,:)' coeff_var_MD_MS_12(2,:)' coeff_var_MD_SS_1(2,:)' coeff_var_MD_MS_12(3,:)' coeff_var_MD_SS_1(3,:)']);


