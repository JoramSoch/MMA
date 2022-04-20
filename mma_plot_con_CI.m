function [cb, CI] = mma_plot_con_CI(MMA, con, xyz, alpha, CI_plot)
% _
% Plot Contrast with Confidence Intervals for Multi-Modal Analysis
% FORMAT [cb, CI] = spm_plot_con_CI(MMA, con, xyz, alpha, CI_plot)
% 
%     MMA     - a structure specifying an estimated MMA model
%     con     - an integer indexing the contrast to be used
%     xyz     - a 1 x 3 vector of MNI coordinates [mm]
%     alpha   - the significance level, CIs are (1-alpha)
%     CI_plot - logical indicating whether to plot CIs
% 
%     cb      - a 1 x q vector of contrast estimates
%     CI      - a 2 x q vector of confidence intervals
%               where q is the 2nd dim of the contrast matrix and
%               where 1st row is lower and 2nd row is upper end of CI
% 
% FORMAT [cb, CI] = mma_plot_con_CI(MMA, con, xyz, alpha, CI_plot)
% computes and displays contrast estimates and (1-alpha) confidence
% intervals of a contrast indexed by con at selected coordinates xyz [1].
% 
% References:
% [1] Soch J (2020). Plot Contrast with Confidence Intervals. GitHub; URL:
%     https://github.com/JoramSoch/spm_helper/blob/master/spm_plot_reg_CI.m
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/02/2022, 10:52
%  Last edit: 10/02/2022, 10:52


% Set defaults
%-------------------------------------------------------------------------%
if isempty(con)     || nargin < 2, con     = 1;      end;
if isempty(xyz)     || nargin < 3, xyz     =[0,0,0]; end;
if isempty(alpha)   || nargin < 4, alpha   = 0.1;    end;
if isempty(CI_plot) || nargin < 5, CI_plot = true;   end;

% Change directory
%-------------------------------------------------------------------------%
SPM_mat = strcat(MMA.swd,'/','SPM.mat');
load(SPM_mat);
try
    cd(SPM.swd)
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

% Load betas and sigma^2
%-------------------------------------------------------------------------%
p = numel(SPM.Vbeta);
b = zeros(p,1);
for j = 1:p
    b_img = spm_read_vols(SPM.Vbeta(j));
    b(j)  = b_img(m_ind(xyz_ind));
end;
s2_img = spm_read_vols(SPM.VResMS);
s2     = s2_img(m_ind(xyz_ind));
covB   = MMA.xX.covB(:,:,xyz_ind);
clear b_img s2_img

% Get contrast estimates
%-------------------------------------------------------------------------%
c  = SPM.xCon(con).c;
q  = size(c,2);
cb = (c'*b)';

% Get confidence intervals
%-------------------------------------------------------------------------%
SE = diag(sqrt(s2 * c'*covB*c))';
CI = SE * norminv(1-alpha/2, 0, 1);
CI = [cb - CI; cb + CI];

% Plot estimates and intervals
%-------------------------------------------------------------------------%
if CI_plot

    figure;
    hold on;
    bar([1:q], cb, 'b');
    errorbar([1:q], cb, cb-CI(1,:), CI(2,:)-cb, '.k', 'LineWidth', 2, 'CapSize', 20);
    axis([(1-0.5), (q+0.5), min([(11/10)*min(CI(1,:)), -(1/10)*max(CI(2,:))]), max([(11/10)*max(CI(2,:)), -(1/10)*min(CI(1,:))])]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:q]);
    xlabel('contrast row/column', 'FontSize', 12);
    ylabel(sprintf('contrast estimate at [%d, %d, %d]', xyz(1), xyz(2), xyz(3)), 'FontSize', 12);
    title(sprintf('Contrast estimates and %d%% confidence intervals: %s', round((1-alpha)*100), SPM.xCon(con).name), 'FontSize', 16);

end;