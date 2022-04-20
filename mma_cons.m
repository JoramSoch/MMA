function mma_cons(MMA, xCon, dec)
% _
% Contrast-Based Inference for Multi-Modal Analyses
% FORMAT mma_cons(MMA, xCon, dec)
%     MMA    - a structure specifying an estimated MMA model
%     xCon   - an 1 x K structure variable with the following fields:
%     o name - a string specifying the name of the contrast
%     o STAT - a character indicating contrast type ('T' or 'F')
%     o c    - a p x 1 contrast vector or a p x q contrast matrix
%     dec    - a logical indicating whether to delete existing contrasts
% 
% FORMAT mma_cons(MMA, xCon, dec) takes a multi-modal analysis model
% specified by MMA and calculates t- or F-contrasts specified in xCon.
% 
% If the input variable "dec" is set to true, existing contrasts are
% deleted, otherwise not. The default value for this variable is false.
% 
% Author: Joram Soch, DZNE Göttingen
% E-Mail: Joram.Soch@DZNE.de
% 
% First edit: 10/02/2022, 09:44
%  Last edit: 10/02/2022, 10:35


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get MMA.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    MMA_mat = spm_select(1,'^MMA\.mat$','Select MMA.mat!');
    MMA_dir = fileparts(MMA_mat); load(MMA_mat);
    MMA.swd = MMA_dir;
    mma_cons(MMA);
    return
else
    orig_dir = pwd;
    SPM_mat  = strcat(MMA.swd,'/','SPM.mat');
    load(SPM_mat);
    cd(SPM.swd);
end;

% Set contrasts if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(xCon)
    i2 = MMA.xX.i2;
    k  = 0;
    for h = 1:numel(i2)
        for j = 1:numel(i2{h})
            k = k + 1;
            xCon(k).name = sprintf('IM %d: reg %d', h, j);
            xCon(k).STAT = 'F';
            c            = zeros(size(MMA.xX.X,2),1);
            c(i2{h}(j))  = 1;
            xCon.c       = c';
        end;
    end;
    clear i2 k c
end;

% Inactivate deletion if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(dec)
    dec = false;
end;

% Load mask image
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);

% Load beta images
%-------------------------------------------------------------------------%
B = MA_load_betas(SPM,m_ind);
v = length(m_ind);
d = floor(v/100);

% Load residual variance
%-------------------------------------------------------------------------%
s2_img = spm_read_vols(SPM.VResMS);
S2     = s2_img(m_ind);
clear s2_img


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','mma_cons: estimate');
spm_progress_bar('Init', 100, 'Calculate voxel-wise contrasts ...', '');

% Finalize contrast information ("spm_contrasts.m", ll. 82-89)
%-------------------------------------------------------------------------%
for k = 1:numel(xCon)
    if size(xCon(k).c,1) < size(SPM.xX.X,2)
        % augment contrast matrix (right-padding with zeros), if necessary
        xCon(k).c = [xCon(k).c; zeros(size(SPM.xX.X,2)-size(xCon(k).c,1),size(xCon(k).c,2))];
    end;
    [X1o, X0]     = spm_SpUtil('c->TSp',SPM.xX.X,xCon(k).c);
    [trMV,trMVMV] = spm_SpUtil('trMV',X1o,SPM.xX.V);
    xCon(k).X0    = X0;
    xCon(k).X1o   = X1o;
    xCon(k).eidf  = trMV^2/trMVMV;
end

% Prepare voxel-wise inference
%-------------------------------------------------------------------------%
covB = MMA.xX.covB;             % parameter covariance
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

% Delete existing contrasts
%-------------------------------------------------------------------------%
if dec
    for k = 1:numel(SPM.xCon)
        SPM.xCon(k).Vcon.fname
        SPM.xCon(k).Vspm.fname
    end;
    SPM = rmfield(SPM,{'xCon'});
    k0  = 0;
else
    k0  = numel(SPM.xCon);
end;

% Save contrast maps
%-------------------------------------------------------------------------%
for k = 1:numel(xCon)
    % contrast value image
    if strcmp(xCon(k).STAT,'T')
        H.fname   = strcat('con_',MF_int2str0(k0+k,4),'.nii');
        H.descrip = sprintf('mma_mma: contrast value image, contrast vector %d', k0+k);
    else
        H.fname   = strcat('ess_',MF_int2str0(k0+k,4),'.nii');
        H.descrip = sprintf('mma_mma: extra sums of squares, contrast matrix %d', k0+k);
    end;
    I(m_ind)  = CONs(k,:);
    H = spm_data_hdr_write(H);
    H = spm_data_write(H,I);
    xCon(k).Vcon = H;
    % statistical parametric map
    H.fname   = strcat('spm',xCon(k).STAT,'_',MF_int2str0(k0+k,4),'.nii');
    H.descrip = sprintf('mma_mma: statistical parametric map, contrast %d', k0+k);
    I(m_ind)  = SPMs(k,:);
    H = spm_data_hdr_write(H);
    H = spm_data_write(H,I);
    xCon(k).Vspm = H;
end;
SPM.xCon(k0+[1:numel(xCon)]) = xCon;
clear xCon;

% Save MMA and SPM structures
%-------------------------------------------------------------------------%
save(strcat(SPM.swd,'/','MMA.mat'),'MMA');
save(strcat(SPM.swd,'/','SPM.mat'),'SPM');
cd(orig_dir);