% Copyright © 2016-2017 LREN CHUV, Jing Cui de Chambrier
% Licenced under the GNU Affero General Public License

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data from HBP workshop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir     = 'C:\DATA\ExternalDrive\Jing_work\LREN\Projects\DiseaseProgression\Staging\ADNI_Imaging_Data'
cd(dir)
RegionVolume    = readtable('ADNI_Data_Volumes_Neuromorphometric_withTIV.xls');
%[colnames,ID]  = xlsread(strcat(dir,'ADNI_Imaging_Data\baseline_imageID.csv'));

%bslID=ID(2:end,:);
bslid           = ismember(RegionVolume.VISCODE,{'bl'});
% bslid           = ~(RegionVolume.APOE==0);
% bslid           = (RegionVolume.APOE==0);
volume          = RegionVolume(:,34:end-1);
nmofsubjects    = sum(bslid);
nmofregions     = size(volume,2);
volumebsl       = table2array(volume(bslid,:));
agebsl          = RegionVolume.Exam_Age(bslid,:);
MRIseq          = RegionVolume.MRIScanner(bslid,:);
genderbsl       = RegionVolume.PTGENDER(bslid,:);
genderbsl       = double(~cellfun('isempty', strfind(genderbsl,'Male'))); %convert it into numeric
ptid            = RegionVolume.PTID(bslid,:);
tivbsl          = RegionVolume.TIV(bslid,:)*10^3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustment Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% =====================================================================
% First scale by total grey matter
% then adjust by age, gender, tiv, MR strength

% scale each region with total grey matter
grmtvlmbsl  = sum(volumebsl');
norvolbsl   = volumebsl.*repmat(1./grmtvlmbsl',1,nmofregions);

% Logit Transformation
log_vol     = log(norvolbsl./(1-norvolbsl));
% remove effects of age, gender
% COVARIATES  = [ones(length(agebsl),1),agebsl(:,1),genderbsl(:,1), MRIseq(:,1)];
COVARIATES  = [ones(length(agebsl),1),agebsl(:,1),tivbsl(:,1),genderbsl(:,1), MRIseq(:,1)];
B1          = inv(COVARIATES'*COVARIATES)*COVARIATES'*log_vol(:,:);
log_vol     = log_vol-COVARIATES(:,2:end)*B1(2:end,:);

% =====================================================================
% First scale by tiv
% then adjust by age, gender

% scale each region with total grey matter
norvolbsl   = volumebsl.*repmat(1./tivbsl,1,nmofregions);

% Logit Transformation
log_vol     = log(norvolbsl./(1-norvolbsl));
% remove effects of age, gender
COVARIATES  = [ones(length(agebsl),1),agebsl(:,1),genderbsl(:,1), MRIseq(:,1)];
B1          = inv(COVARIATES'*COVARIATES)*COVARIATES'*log_vol(:,:);
log_vol     = log_vol-COVARIATES(:,2:end)*B1(2:end,:);

% =====================================================================
% First adjust by age, gender
% then scale by total grey matter

% remove effects of age, gender
COVARIATES   = [ones(length(agebsl),1),agebsl(:,1),genderbsl(:,1)];
B1           = inv(COVARIATES'*COVARIATES)*COVARIATES'*volumebsl(:,:);
volbsl_res   = volumebsl-COVARIATES(:,2:end)*B1(2:end,:);

% scale each region with total grey matter
grmtvlmbsl  = sum(volbsl_res');
norvolbsl   = volbsl_res.*repmat(1./grmtvlmbsl',1,nmofregions);
log_vol     = log(norvolbsl./(1-norvolbsl));

% =====================================================================
% First adjust by age, gender and total grey matter

% calculate total grey matter
grmtvlmbsl  = sum(volumebsl')';

% remove effects of age, gender and total grey matter
COVARIATES      = [ones(length(agebsl),1),agebsl(:,1),genderbsl(:,1),grmtvlmbsl(:,1), MRIseq(:,1)];
B1              = inv(COVARIATES'*COVARIATES)*COVARIATES'*volumebsl(:,:);
norvolbsl_res   = volumebsl-COVARIATES(:,2:end)*B1(2:end,:);

% Logit Transformation
log_vol         = norvolbsl_res;

% =====================================================================
% First adjust by age, gender and TIV

% remove effects of age, gender and total grey matter
COVARIATES      = [ones(length(agebsl),1),agebsl(:,1),genderbsl(:,1),tivbsl(:,1), MRIseq(:,1)];
B1              = inv(COVARIATES'*COVARIATES)*COVARIATES'*volumebsl(:,:);
norvolbsl_res   = volumebsl-COVARIATES(:,2:end)*B1(2:end,:);

% Logit Transformation
log_vol         = norvolbsl_res;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% factoran function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% factor analysis with raw data
[lambda psi T stats F]=factoran(log_vol,1)
writetable(Regions_lambda,'adjust/Loadings_BSLvolume_region.txt','Delimiter',' ')


[F,h,psi,LogLikeli]=FA_EM(log_vol',3,100);
c              = mean(log_vol)';
c              = table(c);
writetable(c,'scale_adjust/Intercept_BSLvolume_region.txt','Delimiter',' ')
Regions        = RegionVolume.Properties.VariableNames(34:end-1)';
Loading1       = F(:,1);
Loading2       = F(:,2);
Loading3       = F(:,3);
Regions_lambda = table(Regions,Loading1,Loading2,Loading3);
writetable(Regions_lambda,'scale_adjust/LoadingMatrix_BSLvolume_region.txt','Delimiter',' ')

Factor1=h(1,:)';
Factor2=h(2,:)';
Factor3=h(3,:)';
T = table(ptid,Factor1,Factor2,Factor3);
writetable(T,'scale_adjust/HiddenFactor_BSLvolume_region.txt','Delimiter',' ')

Psi = diag(psi);
Psi = table(Psi);
writetable(Psi,'scale_adjust/Psi_BSLvolume_region.txt','Delimiter',' ')



% factor analysis with covariance data

[lambda_obs psi_obs T_obs stats_obs]=factoran(log_vol*log_vol',1,'xtype','cov');

[V, D]  = eig(log_vol*log_vol');
D       = D+diag(abs(min(D))); % or 0 are same
newMat  = V*D*V';
[lambda_obs psi_obs T_obs stats_obs]    = factoran(newMat,1,'xtype','cov');

% write the results into a file
ptid=ptid(~cellfun('isempty', bslid),:);
T=[ptid,num2cell(F),repmat({'bl'},nmofsubjects,1),num2cell(log_vol)];
xlswrite('C:\jcui\DiseaseProgression\CommonFactor_BSLvolume_region.csv',T);
Lambda=num2cell(lambda);
xlswrite('C:\jcui\DiseaseProgression\Lambda_BSLvolume_region.csv',Lambda);
Psi=num2cell(psi);
xlswrite('C:\jcui\DiseaseProgression\Psi_BSLvolume_region.csv',Psi);


% subject space

[F,h,LogLikeli] = FA_EM(log_vol,3,100);
log_vol         = log_vol-repmat(mean(log_vol')',1,size(log_vol,2));
weights_F       = log_vol'*F;
Regions_lambda  = table(RegionVolume.Properties.VariableNames(34:end)',weights_F(:,1));
