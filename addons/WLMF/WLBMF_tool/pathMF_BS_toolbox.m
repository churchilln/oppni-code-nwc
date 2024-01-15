%%%%%% SETTING THE PATHS FOR WLBMF_tool
%
% Herwig Wendt, Lyon, 2006 - 2008

tmp=which('pathMF_BS_toolbox'); index=strfind(tmp,'\');
p=tmp(1:index(end));

addpath([p,'core'],[p,'misc'],[p,'core/constancy'],[p,'misc/constancy']);

% clear p tmp index
