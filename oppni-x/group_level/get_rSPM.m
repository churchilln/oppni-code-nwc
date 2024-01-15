function [ rep, rSPM ] = get_rSPM( vect1, vect2, keepMean )
%
%    [ rep, rSPM ] = get_rSPM( vect1, vect2, keepMean )
%
%  * This script takes in 2 vectors, returns reproducibility (rep)
%    and reproducible SPM (rSPM)
%
%  * keepMean: reinsert the mean offset present in vectors
%
%
% ------------------------------------------------------------------------%
% Author: Nathan Churchill, University of Toronto
%  email: nathan.churchill@rotman.baycrest.on.ca
% ------------------------------------------------------------------------%
% version history: March 15 2012
% ------------------------------------------------------------------------%

rep = corr(vect1, vect2);

%(1) getting the mean offsets (normed by SD)
normedMean1 = mean(vect1)./std(vect1);
normedMean2 = mean(vect2)./std(vect2);
%    and rotating means into signal/noise axes
sigMean = (normedMean1 + normedMean2)/sqrt(2);
%noiMean = (normedMean1 - normedMean2)/sqrt(2);
%(2) getting  signal/noise axis projections of (zscored) betamaps
sigProj = ( zscore(vect1) + zscore(vect2) ) / sqrt(2);
noiProj = ( zscore(vect1) - zscore(vect2) ) / sqrt(2);
% noise-axis SD
noiStd = std(noiProj);
%(3) norming by noise SD:
%     ...getting the (re-normed) mean offsets
sigMean = sigMean./noiStd;
%noiMean = noiMean./noiStd; 
size(sigProj),
size(noiStd),
%  getting the normed signal/noise projection maps
sigProj = sigProj ./ noiStd;
%noiProj = noiProj ./ noiStd;

% Produce the rSPM:
if    ( keepMean == 1 )   rSPM = sigProj + sigMean;
elseif( keepMean == 0 )   rSPM = sigProj;
end
