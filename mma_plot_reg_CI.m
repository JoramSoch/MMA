function [x, y, xr, yr, CI] = mma_plot_reg_CI(MMA, reg, xyz, alpha, reg_plot)
% _
% Plot Regression with Confidence Interval for Multi-Modal Analysis
% FORMAT [x, y, xr, yr, CI] = mma_plot_reg_CI(MMA, reg, xyz, alpha, reg_plot)
% 
%     MMA      - a structure specifying an estimated MMA model
%     reg      - a 1 x C vector indexing the regressors to be used
%     xyz      - a 1 x 3 vector of MNI coordinates [mm]
%     alpha    - the significance level, CIs are (1-alpha)
%     reg_plot - logical indicating whether to plot regression
% 
%     x, y, xr, yr, CI - 1 x C cell arrays containing:
%     x        - an nj x 1 vector of independent variable values
%     y        - an nj x 1 vector of dependent variable values
%     xr       - a 100 x 1 vector of regression line x-values
%     yr       - a 100 x 1 vector of regression line y-values
%     CI       - a 2 x 100 matrix of regression confidence bands
% 
% FORMAT [x, y, xr, yr, CI] = mma_plot_reg_CI(MMA, reg, xyz, alpha, reg_plot)
% computes and displays regression line as well as (1-alpha) confidence
% bands [1,2] for a linear relationship of measured responses with the
% regressors indexed by reg at selected coordinates xyz [1].
% 
% References:
% [1] Soch J (2020). Plot Regression with Confidence Interval. GitHub; URL:
%     https://github.com/JoramSoch/spm_helper/blob/master/spm_plot_reg_CI.m
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/02/2022, 14:42
%  Last edit: 10/02/2022, 14:42


% Set defaults
%-------------------------------------------------------------------------%
if isempty(reg)      || nargin < 2, reg      = 1;      end;
if isempty(xyz)      || nargin < 3, xyz      =[0,0,0]; end;
if isempty(alpha)    || nargin < 4, alpha    = 0.1;    end;
if isempty(reg_plot) || nargin < 5, reg_plot = true;   end;

% Change directory
%-------------------------------------------------------------------------%
SPM_mat = strcat(MMA.swd,'/','SPM.mat');
load(SPM_mat);
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

% Load mask image
%-------------------------------------------------------------------------%
 m_dim           = SPM.VM.dim;
[m_img, m_xyz]   = spm_read_vols(SPM.VM);
 m_ind           = find(m_img~=0);
 d_img           = sqrt(sum((m_xyz(:,m_ind) - repmat(xyz',[1 numel(m_ind)])).^2, 1));
[d_min, xyz_ind] = min(d_img);
clear m_dim m_img m_xyz d_img d_min

% Load all models
%-------------------------------------------------------------------------%
load(MMA.design.SPM_mat_DM);                % dependent modality SPM.mat
SPM1  = SPM;
nIMs  = numel(MMA.design.SPM_mat_IM);       % number of independent modalities
for h = 1:nIMs
    load(MMA.design.SPM_mat_IM{h});         % independent modality SPM.mat
    SPM2(h) = SPM;
end;
clear SPM
load(SPM_mat);

% Load all data
%-------------------------------------------------------------------------%
Finter = spm('FigName','mma_plot_reg_CI: load');
spm_progress_bar('Init',100,'Load in-mask time series...','');
n  = size(SPM1.xX.X,1);
d  = floor(n/100);
y1 = zeros(n,1);
y2 = zeros(n,nIMs);
for i = 1:n
    y_img = spm_read_vols(MMA.xY.DM.VY(i));
    y1(i) = y_img(m_ind(xyz_ind));
    for h = 1:nIMs
        y_img   = spm_read_vols(MMA.xY.IM(h).VY(i));
        y2(i,h) = y_img(m_ind(xyz_ind));
    end;
    if mod(i,d) == 0, spm_progress_bar('Set',(i/n)*100); end;
end;
clear y_img
spm_progress_bar('Clear');

% Get design matrices
%-------------------------------------------------------------------------%
con_mat_IM = MMA.design.con_mat_IM;
X2 = cell(1,nIMs);              % additional design matrices
p2 = zeros(1,nIMs);             % numbers of additional regressors
for h = 1:nIMs
    X2{h} = SPM2(h).xX.X;       % design matrices for
    X2{h} = X2{h}*con_mat_IM{h};% independent modalities,
    p2(h) = size(X2{h},2);      % reduced by contrast matrix
end;
clear SPM1 SPM2

% Create design matrix
%-------------------------------------------------------------------------%
i2 = MMA.xX.i2;
Xi = SPM.xX.X;
for h = 1:nIMs
    Xh = X2{h};
    for j = 1:p2(h)
        ij = find(Xh(:,j)==1);
        Xh(ij,j) = y2(ij,h);
    end;
    Xi(:,i2{h}) = Xh;
    clear Xh ij
end;

% Get signal and regressor
%-------------------------------------------------------------------------%
C  = numel(reg);
x  = cell(1,C);
y  = cell(1,C);
xL = cell(1,C);
for j = 1:C
    x{j}  = Xi( abs(Xi(:,reg(j)))>exp(-23), reg(j) );
    y{j}  = y1( abs(Xi(:,reg(j)))>exp(-23) );
    xL{j} = SPM.xX.name{reg(j)};
end;
clear X

% Perform linear regression
%-------------------------------------------------------------------------%
b = zeros(2,C);
for j = 1:C
    X{j}   = [ones(numel(x{j}),1), x{j}];
    b(:,j) = (X{j}'*X{j})^(-1) * X{j}'*y{j};
end;
clear Xj

% Find variable ranges
%-------------------------------------------------------------------------%
y_min = min(y1)-1/20*range(y1);
y_max = max(y1)+1/20*range(y1);
x_min = zeros(1,C);
x_max = zeros(1,C);
for j = 1:C
    x_min(j) = min(x{j})-1/20*range(x{j});
    x_max(j) = max(x{j})+1/20*range(x{j});
end;

% Create regression line
%-------------------------------------------------------------------------%
m  = 100;
xr = cell(1,C);
yr = cell(1,C);
for j = 1:C
    xr{j} = [x_min(j):((x_max(j)-x_min(j))/(m-1)):x_max(j)]';
    Xr    = [ones(m,1), xr{j}];
    yr{j} = Xr*b(:,j);
end;
clear Xr

% Create confidence intervals
%-------------------------------------------------------------------------%
CI = cell(1,C);
for j = 1:C
    nj    = numel(y{j});
    SSX   = sum((x{j} - mean(x{j})).^2);
    SSE   = sum((y{j} - X{j}*b(:,j)).^2);
    SE_yr = sqrt(SSE/(nj-2)) * sqrt( 1/nj + 1/SSX*((xr{j}-mean(x{j})).^2) );
    CI_yr = SE_yr * tinv(1-alpha, nj-2);
    CI{j} = [yr{j} - CI_yr, yr{j} + CI_yr]';
end;

% Plot estimates and intervals
%-------------------------------------------------------------------------%
if reg_plot
    
    % open figure and prepare color map
    figure;
    colormap('default');
    cmap = colormap;
    cols = zeros(C,3);
    for j = 1:C
        cols(j,:) = cmap(1+round((size(cmap,1)-1)*((j-1)/C)),:);
    end;
    
    % plot regressors as separate lines
    hold on;
    for j = 1:C
        plot(min(x_min)-1, y_min-1, '-', 'Color', cols(j,:), 'LineWidth', 2);
    end;
    for j = 1:C
        plot(x{j}, y{j}, '.', 'Color', cols(j,:), 'MarkerSize', 15);
        plot(xr{j}, yr{j}, '-', 'Color', cols(j,:), 'LineWidth', 2);
        plot(xr{j}', CI{j}(1,:), '--', 'Color', cols(j,:), 'LineWidth', 2);
        plot(xr{j}', CI{j}(2,:), '--', 'Color', cols(j,:), 'LineWidth', 2);
    end;
    axis([min(x_min), max(x_max), y_min y_max]);
    set(gca,'Box','On');
    legend(xL, 'Location', 'Best');
    xlabel('regressor value', 'FontSize', 12);
    ylabel(sprintf('measured response at [%d, %d, %d]', xyz(1), xyz(2), xyz(3)), 'FontSize', 12);
    title(sprintf('Predicted responses and %d%% confidence bands', round((1-alpha)*100)), 'FontSize', 16);

end;