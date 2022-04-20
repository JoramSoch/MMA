function mma_mma(SPM_mat_DM, SPM_mat_IM, con_mat_IM, MMA_dir)
% _
% Multi-Modal Analyses using Voxel-Wise General Linear Models
% FORMAT mma_mma(SPM_mat_DM, SPM_mat_IM, con_mat_IM, MMA_dir)
%     SPM_mat_DM - a string indicating the filepath to the
%                  SPM.mat of the dependent modality
%     SPM_mat_IM - a cell array containing filepaths to the
%                  SPM.mats of the independent modalities
%     con_mat_IM - a cell array of contrast matrices used for
%                  reduction of the independent modalities
%     MMA_dir    - a string specifying the multi-modal analysis directory
% 
% FORMAT mma_mma(SPM_mat_DM, SPM_mat_IM, con_mat_IM, MMA_dir) loads SPM.mat
% files specifying general linear models (GLMs) used for analyzing a
% dependent modality (DM) and one or more independent modalities (IM), 
% builds a multi-modal analysis (MMA), using the design matrix from the DM
% analysis and adding the data/signals from the IM analyses as voxel-wise
% covariates, possibly reduced to sub-groups via contrast matrices, and
% writes the results into an MMA directory.
% 
% Note: Scans must be in the same order for each SPM.mat underlying the
% MMA model, i.e. the i-th row in each GLM must correspond to the same
% measurement unit (e.g. subject), acquired under different modalities
% (e.g. task-fMRI, rest-fMRI, structural MRI). Moreover, the IM models must
% only contain categorical regressors (e.g. one-sample t-test, two-sample
% t-test, two-way ANOVA etc., without covariates).
% 
% The purpose of the input variable "con_mat_IM" is to reduce the voxel-
% wise IM signals to sub-groups of the DM analysis, in case that they
% should not be distinguished in the same way that the DM signals are
% distinguished (for an example, see below). This is achieved by right-
% multiplying the design matrix from each independent modality with the
% respective entry of the cell array "con_mat_IM".
% 
% An exemplary call to this function could therefore be done as follows:
%     SPM_mat_DM = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MS_FADE_04_FADE_memory_as_anova2\SPM.mat';
%     SPM_mat_IM ={'C:\Joram\projects\DZNE\FADE\analyses_BS\resting_state_fMRI\mPerAF_PerAF_as_anova2\SPM.mat', ...
%                  'C:\Joram\projects\DZNE\FADE\analyses_BS\structural_analyses\s6_VBM_GMV_as_anova2\SPM.mat'};
%     con_mat_IM ={[1, 1, 0, 0; 0, 0, 1, 1]', [1, 1, 1, 1]'};
%     MMA_dir    = 'C:\Joram\projects\DZNE\FADE\analyses_BS\group_statistics\MMA_playground\FADE_memory_as_anova2_vs_Aprime_HCvol_PerAF_GMV\';
%     mma_mma(SPM_mat_DM, SPM_mat_IM, con_mat_IM, MMA_dir);
% 
% Example: Consider the following data set comprising three modalities
% (note that the values for DM, IM1 and IM2 are changing across voxels):
% Subj   DM A1 A2 B1 B2  C1  C2  C3   IM1 A1 A2 B1 B2   IM2 A1 A2 B1 B2
%    1   y1  1  0  0  0 c11 c21 c31   y11  1  0  0  0   y21  1  0  0  0
%    2   y2  1  0  0  0 c12 c22 c32   y12  1  0  0  0   y22  1  0  0  0
%    3   y3  0  1  0  0 c13 c23 c33   y13  0  1  0  0   y23  0  1  0  0
%    4   y4  0  1  0  0 c14 c24 c34   y14  0  1  0  0   y24  0  1  0  0
%    5   y5  0  0  1  0 c15 c25 c35   y15  0  0  1  0   y25  0  0  1  0
%    6   y6  0  0  1  0 c16 c26 c36   y16  0  0  1  0   y26  0  0  1  0
%    7   y7  0  0  0  1 c17 c27 c37   y17  0  0  0  1   y27  0  0  0  1
%    8   y8  0  0  0  1 c18 c28 c38   y18  0  0  0  1   y28  0  0  0  1
%       |------- SPM_mat_DM -------| |---------- SPM_mat_IM -----------|
% 
% Let the contrast matrices mapping the independent modalities to sub-
% groups of the dependent modality analysis be defined as follows:
%     con_mat_IM = {[1 1 0 0; 0 0 1 1]', [1 1 1 1]'}.
% This would specify that for IM1, the groups A1/A2 and B1/B2 are collapsed,
% whereas for IM2, all four groups are collapsed into a single voxel-wise
% regressor.
% 
% This data set and the contrast would result in the following MMA model
% (note that again, values for DM, IM1 and IM2 are voxel-dependent):
% Subj   DM A1 A2 B1 B2  C1  C2  C3 IM1-A IM1-B IM2-all
%    1   y1  1  0  0  0 c11 c21 c31   y11*    0     y21*
%    2   y2  1  0  0  0 c12 c22 c32   y12*    0     y22*
%    3   y3  0  1  0  0 c13 c23 c33   y13*    0     y23*
%    4   y4  0  1  0  0 c14 c24 c34   y14*    0     y24*
%    5   y5  0  0  1  0 c15 c25 c35     0   y15*    y25*
%    6   y6  0  0  1  0 c16 c26 c36     0   y16*    y26*
%    7   y7  0  0  0  1 c17 c27 c37     0   y17*    y27*
%    8   y8  0  0  0  1 c18 c28 c38     0   y18*    y28*
%       |------- SPM_mat_DM -------|---- SPM_mat_IM ----|
%    *  all IM values are mean-centered within each regressor
% 
% This design matrix would allow to test for a general effect of IM1 (using
% an F-contrast over regressors "IM1-A" and "IM1-B") and for an interaction
% of IM1 with A/B (by contrasting the regressors "IM1-A" and "IM1-B"), but
% only for a general effect of IM2 (using a contrast on "IM2-all") and not
% for interactions of IM2 with A/B or 1/2.
% 
% If the SPM.mat files specifying DM and IM models contain contrasts, these
% contrasts are automatically taken over to the MMA model and prefixed with
% "IM1", "IM2" etc. For independent modalities that are reduced to sub-
% groups, all estimable contrasts are executed. For example, in the above
% case where IM1 is reduced to A/B (collapsing 1/2), the main effect of A/B
% would still be estimable, but the main effect of 1/2 is discarded.
% 
% If no collapsing of independent modalities is desired, the contrast
% vector variable would consist of identity matrices:
%     con_mat_IM = {[eye(4)], [eye(4)]}
% This would distinguish the IM signals across all groups that are also
% used for modelling the DM, possibly allowing for the same contrasts.
% 
% The mask image used in the MMA model is the intersection of the masks
% from all SPM.mats, i.e. only voxels which have data on all modalities
% from all subjects are used for the analysis.
% 
% SPM's temporal non-sphericity estimation (cf. spm_est_non_sphericity) is
% performed before voxel-wise parameter estimation, using a mean design
% matrix, averaged across IM signals from all voxels. This is because SPM
% uses the same temporal covariance matrix for all voxels [1]. The average
% design matrix is also used for results display in SPM.
% 
% SPM's spatial smoothness estimation (cf. spm_est_smoothness) is
% performed after voxel-wise parameter estimation, using standardized
% residuals based on voxel-wise design matrices. This is important for
% family wise error (FWE) correction based on random field theory (RFT) 
% used for multiple comparison correction when interrogating results in SPM.
% 
% This code uses the following functions from the MACS toolbox [1,2]:
% o MA_load_mask
% o MA_load_data
% o MA_init_header
% o MS_create_mask
% o MF_int2str0
% 
% References:
% [1] Soch J & Allefeld C (2018). MACS - a new SPM toolbox for model
%     assessment, comparison and selection. Journal of Neuroscience
%     Methods, vol. 306, pp. 19-31; DOI: 10.1016/j.jneumeth.2018.05.017.
% [2] Soch J (2017). MACS – a new SPM toolbox for model assessment, 
%     comparison and selection. GitHub, first uploaded 2017-05-24;
%     URL: https://github.com/JoramSoch/MACS.
% 
% Author: Joram Soch, DZNE Göttingen
% E-Mail: Joram.Soch@DZNE.de
% 
% First edit: 13/01/2022, 00:18
%  Last edit: 19/01/2022, 21:51


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Load SPM.mats
%-------------------------------------------------------------------------%
orig_dir = pwd;
load(SPM_mat_DM);               % dependent modality SPM.mat
SPM1  = SPM;
nIMs  = numel(SPM_mat_IM);      % number of independent modalities
for h = 1:nIMs
    load(SPM_mat_IM{h});        % independent modality SPM.mat
    SPM2(h) = SPM;
end;                            % store MMA design information
MMA.design.SPM_mat_DM = SPM_mat_DM;
MMA.design.SPM_mat_IM = SPM_mat_IM;
MMA.design.con_mat_IM = con_mat_IM;
MMA.design.MMA_dir    = MMA_dir;
MMA.swd               = MMA_dir;

% Load mask images
%-------------------------------------------------------------------------%
cd(SPM1.swd);
[M1, m_dim, m_ind] = MA_load_mask(SPM1);
 V  = prod(m_dim);              % total number of voxels
 M2 = zeros(nIMs,V);
for h = 1:nIMs
    cd(SPM2(h).swd);
    [M2(h,:), m_dim, m_ind] = MA_load_mask(SPM2(h));
end;
clear m_dim m_ind
[M, m_hdr, m_ind]  = MS_create_mask([M1; M2], SPM1.VM);
 v  = numel(m_ind);             % number of in-mask voxels
 d  = floor(v/100);
MMA.xM.VM = m_hdr;

% Load data images
%-------------------------------------------------------------------------%
cd(SPM1.swd);
Y1 = MA_load_data(SPM1, m_ind);
Y2 = zeros(size(Y1,1),numel(m_ind),nIMs);
MMA.xY.DM.VY = SPM1.xY.VY;
for h = 1:nIMs
    cd(SPM2(h).swd);
    Y2(:,:,h) = MA_load_data(SPM2(h), m_ind);
    MMA.xY.IM(h).VY = SPM2(h).xY.VY;
end;

% Get design matrices
%-------------------------------------------------------------------------%
X1 = SPM1.xX.X;                 % basic design matrix
n1 = size(X1,1);                % number of scans (original)
p1 = size(X1,2);                % number of regressors (original)
i1 = [1:p1];                    % indices of regressors (original)
X2 = cell(1,nIMs);              % additional design matrices
p2 = zeros(1,nIMs);             % numbers of regressors (additional)
i2 = cell(1,nIMs);              % indices of regressors (additional)
for h = 1:nIMs
    X2{h} = SPM2(h).xX.X;       % design matrices for
    X2{h} = X2{h}*con_mat_IM{h};% independent modalities,
    p2(h) = size(X2{h},2);      % reduced by contrast matrix
    i2{h} = p1+sum(p2(1:h-1))+[1:p2(h)];
end;
n = n1;                         % number of scans (augmented)
p = p1 + sum(p2);               % number of regressors (augmented)


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   N E W   M O D E L                       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','mma_mma: estimate (1)');

% Assemble new design matrix
%-------------------------------------------------------------------------%
X = zeros(n,p);
X(:,i1) = X1;
for h = 1:nIMs
    Xh = X2{h};
    for j = 1:p2(h)
        % identify observations belonging to this category in this model
        ij = find(Xh(:,j)==1);
        % add average mean-centered modality values to average design matrix
        Xh(ij,j) = mean( Y2(ij,:,h)-repmat(mean(Y2(ij,:,h)),[numel(ij) 1]), 2 );
    end;
    X(:,i2{h}) = Xh;
    clear Xh ij
end;

% Assemble new SPM structure
%-------------------------------------------------------------------------%
SPM     = SPM1;
SPM.swd = MMA_dir;
SPM.VM  = m_hdr;
if isfield(SPM.xVi,'Vi')
    SPM.xVi = rmfield(SPM.xVi,{'h','V','Cy'});
end;
if isfield(SPM,'xCon')
    SPM     = rmfield(SPM,{'xCon'});
end;
SPM     = rmfield(SPM,{'Vbeta','VResMS'});

% Assemble new design matrix
%-------------------------------------------------------------------------%
SPM.xX.X    = X;                % average design matrix
SPM.xX.iH   = SPM1.xX.iH;
SPM.xX.iC   =[SPM1.xX.iC, [(p1+1):p]];
SPM.xX.name = cell(p,1);
SPM.xX.name(1:p1) = SPM1.xX.name;
for h = 1:nIMs                  % additional regressor names
    for j = 1:p2(h)
        jj = find(con_mat_IM{h}(:,j)==1);
        if numel(jj) == 0
            SPM.xX.name{i2{h}(j)} = sprintf('IM%d:unknown', h);
        elseif numel(jj) == 1
            SPM.xX.name{i2{h}(j)} = sprintf('IM%d:%s', h, SPM2(h).xX.name{jj});
        elseif numel(jj) > 1
            str = 'groups-';
            for k = 1:numel(jj), str = sprintf('%s%d,', str, jj(k)); end;
            SPM.xX.name{i2{h}(j)} = sprintf('IM%d:%s', h, str(1:end-1));
        end;
        clear jj str
    end;
end;
SPM.xX = rmfield(SPM.xX,{'V','W','pKX','xKXs','trRV','trRVRV','erdf','Bcov','nKX'});

% Write new mask image
%-------------------------------------------------------------------------%
if ~exist(SPM.swd,'dir'), mkdir(SPM.swd); end;
cd(SPM.swd);
spm_write_vol(SPM.VM, reshape(M, SPM.VM.dim));
[x,y,z] = ind2sub(SPM.VM.dim, m_ind');
XYZ     = [x, y, z]';
clear x y z

% Estimate new non-sphericity ("spm_spm.m", ll. 433, 443-444)
%-------------------------------------------------------------------------%
if isfield(SPM.xVi,'Vi')
    [xVi, msk] = spm_est_non_sphericity(SPM);
    W          = spm_sqrtm(spm_inv(xVi.V));
    W          = W.*(abs(W) > 1e-6);
    SPM.xX.W   = W;
    SPM.xVi    = xVi;
    clear xVi msk W
else
    SPM.xX.W  = eye(n);
    SPM.xVi.V = eye(n);
end;
fprintf('\n');

% Calculate new whitened design matrix ("spm_spm.m", ll. 448-453, 709-711)
%-------------------------------------------------------------------------%
xX = SPM.xX;
xX.xKXs   = spm_sp('Set',spm_filter(xX.K,xX.W*xX.X));
xX.xKXs.X = full(xX.xKXs.X);
xX.pKX    = spm_sp('x-',xX.xKXs);
xX.nKX    = spm_DesMtx('sca',xX.xKXs.X,xX.name);

% Calculate new residual degrees of freedom ("spm_spm.m", ll. 455-462)
%-------------------------------------------------------------------------%
xX.V           = spm_filter(xX.K,spm_filter(xX.K,xX.W*SPM.xVi.V*xX.W')');
[trRV, trRVRV] = spm_SpUtil('trRV',xX.xKXs,xX.V);
xX.trRV        = trRV;
xX.trRVRV      = trRVRV;
xX.erdf        = trRV^2/trRVRV;
xX.Bcov        = xX.pKX*xX.V*xX.pKX';
SPM.xX = xX;

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   P A R A M E T E R S                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','mma_mma: estimate (2)');
spm_progress_bar('Init', 100, 'Estimate voxel-wise parameters ...', '');

% Prepare voxel-wise estimation
%-------------------------------------------------------------------------%
Y = Y1; clear Y1;               % data matrix
X = SPM.xX.X;                   % design matrix
V = SPM.xVi.V;                  % covariance matrix
P = spm_inv(V);                 % precision matrix
B = zeros(p,v);                 % regression coefficients
R = zeros(n,v);                 % residual values
S2= zeros(1,v);                 % residual variance
covB = zeros(p,p,v);            % parameter covariance
erdf = SPM.xX.trRV;             % effective residual DOFs

% Estimate voxel-wise parameters
%-------------------------------------------------------------------------%
for i = 1:v
    
    % Create voxel-wise design matrix
    %---------------------------------------------------------------------%
    Xi = X;
    for h = 1:nIMs
        Xh = X2{h};
        for j = 1:p2(h)
            % identify observations belonging to this category
            ij = find(Xh(:,j)==1);
            % add mean-centered modality values to design matrix
            Xh(ij,j) = Y2(ij,i,h)-mean(Y2(ij,i,h));
        end;
        Xi(:,i2{h}) = Xh;
        clear Xh ij
    end;
    
    % Calculate voxel-wise parameters
    %---------------------------------------------------------------------%
    covB(:,:,i) = inv(Xi'*P*Xi);
    B(:,i)      = covB(:,:,i) * Xi'*P*Y(:,i);
    R(:,i)      = Y(:,i) - Xi*B(:,i);
    S2(i)       = (1/erdf) * R(:,i)'*P*R(:,i);
    
    % Update progress bar
    %---------------------------------------------------------------------%
    if mod(i,d) == 0, spm_progress_bar('Set',(i/v)*100); end;
    
end;
clear Xi Xh ij

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Store design information
%-------------------------------------------------------------------------%
MMA.xX.X    = X;
MMA.xX.V    = V;
MMA.xX.covB = covB;
MMA.xX.erdf = erdf;
MMA.xX.i1   = i1;
MMA.xX.i2   = i2;

% Save standardized residuals ("spm_spm.m", ll. 512, 524-528, 635-640)
%-------------------------------------------------------------------------%
nSRs = min([n, spm_get_defaults('stats.maxres')]);
Rstd = R(1:nSRs,:)./repmat(sqrt(S2),[nSRs 1]);
H    = MA_init_header(SPM, false);
I    = NaN(SPM.VM.dim);
for i = 1:nSRs
    I(m_ind)  = Rstd(i,:);
    H.fname   = strcat('ResI_',MF_int2str0(i,4),'.nii');
    H.descrip = sprintf('mma_mma: standardized residuals, data point %d', i);
    H = spm_data_hdr_write(H);
    H = spm_data_write(H,I);
    VResI(i) = H;
end;

% Estimate new smoothness ("spm_spm.m", ll. 675-676, 718-727)
%-------------------------------------------------------------------------%
[FWHM, VRpv, Rc] = spm_est_smoothness(VResI, SPM.VM, [n erdf]);
SPM.xVol.XYZ     = XYZ;         % in-mask voxel coordinates
SPM.xVol.FWHM    = FWHM;        % full width at half maximum
SPM.xVol.R       = Rc;          % resel count
SPM.xVol.S       = v;           % voxel volume
SPM.xVol.VRpv    = VRpv;        % resels per voxel
clear FHWM VRpv R
fprintf('\n');

% Delete standardized residuals ("spm_spm.m", ll. 697-702) 
%-------------------------------------------------------------------------%
for i = 1:nSRs
    delete(VResI(i).fname);
end;


%=========================================================================%
% E S T I M A T I O N   ( 3 ) :   C O N T R A S T S                       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','mma_mma: estimate (3)');
spm_progress_bar('Init', 100, 'Calculate voxel-wise contrasts ...', '');

% Assemble new contrasts
%-------------------------------------------------------------------------%
if isfield(SPM1,'xCon')         % dependent modality contrasts
    for j = 1:numel(SPM1.xCon)
        xCon(j).name = SPM1.xCon(j).name;
        xCon(j).STAT = SPM1.xCon(j).STAT;
        xCon(j).c    =[SPM1.xCon(j).c', zeros(size(SPM1.xCon(j).c',1),sum(p2))]';
    end;
    k = numel(SPM1.xCon);
else
    k = 0;
end;
k = k + 1;                      % effects of interest contrast (DM)
xCon(k).name = 'DM: EOI';
xCon(k).STAT = 'F';
xCon(k).c    = [eye(p1), zeros(p1,sum(p2))]';
for h = 1:nIMs
    if isfield(SPM2(h),'xCon')  % independent modality contrasts
        for j = 1:numel(SPM2(h).xCon)
            cj = SPM2(h).xCon(j).c' * con_mat_IM{h};
            cj = flipud(unique(cj,'rows'));
            if sum(abs(cj)) > 0 % if contrast vector is applicable
                k = k + 1;      % after reduction by modality contrast
                xCon(k).name = sprintf('IM %d: %s', h, SPM2(h).xCon(j).name);
                xCon(k).STAT = SPM2(h).xCon(j).STAT;
                c            = zeros(size(cj,1),p);
                c(:,i2{h})   = cj;
                xCon(k).c    = c';
            end;
        end;
    end;
    k = k + 1;                  % effects of interest contrast (IMs)
    xCon(k).name = sprintf('IM %d: EOI', h);
    xCon(k).STAT = 'F';
    c            = zeros(p2(h),p);
    c(:,i2{h})   = eye(p2(h));
    xCon(k).c    = c';
    clear c cj
end;
k = k + 1;                      % effects of interest contrast (all)
xCon(k).name = 'all: EOI';
xCon(k).STAT = 'F';
xCon(k).c    = eye(p);

% Finalize contrast information ("spm_contrasts.m", ll. 82-89)
%-------------------------------------------------------------------------%
for k = 1:numel(xCon)
    [X1o, X0]     = spm_SpUtil('c->TSp',SPM.xX.X,xCon(k).c);
    [trMV,trMVMV] = spm_SpUtil('trMV',X1o,SPM.xX.V);
    xCon(k).X0    = X0;
    xCon(k).X1o   = X1o;
    xCon(k).eidf  = trMV^2/trMVMV;
end

% Prepare voxel-wise inference
%-------------------------------------------------------------------------%
CONs = zeros(numel(xCon),v);    % contrast/ESS values
SPMs = zeros(numel(xCon),v);    % statistical parametric maps
    
% Calculate voxel-wise contrasts
%-------------------------------------------------------------------------%
for i = 1:v
    
    % Compute voxel-wise statistics
    %---------------------------------------------------------------------%
    for k = 1:numel(xCon)
        
        % Compute t-statistics ("ME_GLM_con.m", ll. 87-95)
        %-----------------------------------------------------------------%
        if strcmp(xCon(k).STAT,'T')
            c         = xCon(k).c;
            CONs(k,i) = c'*B(:,i);          % contrast values
            SPMs(k,i) = CONs(k,i)/sqrt(S2(i) * c'*covB(:,:,i)*c);
        end;
        
        % Compute F-statistics ("ME_GLM_con.m", ll. 104-115)
        %-----------------------------------------------------------------%
        if strcmp(xCon(k).STAT,'F')
            C         = xCon(k).c;
            q         = size(C,2);
            CON       = C'*B(:,i);          % extra sums of squares
            CONs(k,i) = CON' * inv(C'*covB(:,:,i)*C) * CON;
            SPMs(k,i) = (CONs(k,i)/q) / S2(i);
        end;
        
    end;
    
    % Update progress bar
    %---------------------------------------------------------------------%
    if mod(i,d) == 0, spm_progress_bar('Set',(i/v)*100); end;
    
end;
clear c C q CON

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','mma_mma: save');

% Initialise image files
%-------------------------------------------------------------------------%
cd(SPM.swd);
H = MA_init_header(SPM, false);
I = NaN(SPM.VM.dim);

% Save residual variance
%-------------------------------------------------------------------------%
H.fname   = 'ResMS.nii';
H.descrip = 'mma_mma: voxel-wise residual variance';
I(m_ind)  = S2;
H = spm_data_hdr_write(H);
H = spm_data_write(H,I);
SPM.VResMS = H;

% Save beta estimates
%-------------------------------------------------------------------------%
for j = 1:p
    % regression coefficients
    H.fname   = strcat('beta_',MF_int2str0(j,4),'.nii');
    H.descrip = sprintf('mma_mma: voxel-wise beta estimates, regression coefficient %d', j);
    I(m_ind)  = B(j,:);
    H = spm_data_hdr_write(H);
    H = spm_data_write(H,I);
    SPM.Vbeta(j) = H;
end;

% Save contrast maps
%-------------------------------------------------------------------------%
for k = 1:numel(xCon)
    % contrast value image
    if strcmp(xCon(k).STAT,'T')
        H.fname   = strcat('con_',MF_int2str0(k,4),'.nii');
        H.descrip = sprintf('mma_mma: contrast value image, contrast vector %d', k);
    else
        H.fname   = strcat('ess_',MF_int2str0(k,4),'.nii');
        H.descrip = sprintf('mma_mma: extra sums of squares, contrast matrix %d', k);
    end;
    I(m_ind)  = CONs(k,:);
    H = spm_data_hdr_write(H);
    H = spm_data_write(H,I);
    xCon(k).Vcon = H;
    % statistical parametric map
    H.fname   = strcat('spm',xCon(k).STAT,'_',MF_int2str0(k,4),'.nii');
    H.descrip = sprintf('mma_mma: statistical parametric map, contrast %d', k);
    I(m_ind)  = SPMs(k,:);
    H = spm_data_hdr_write(H);
    H = spm_data_write(H,I);
    xCon(k).Vspm = H;
end;
SPM.xCon = xCon;
clear xCon;

% Save MMA and SPM structures
%-------------------------------------------------------------------------%
save(strcat(SPM.swd,'/','MMA.mat'),'MMA');
save(strcat(SPM.swd,'/','SPM.mat'),'SPM');
cd(orig_dir);