function [paramBS, paramEst] = MF_BS_test_constancy_param(B1,B2,Block,Method,T_S,Alpha,Fun,j1,j2,wtype,Jflag,MomNul,methodWT,Cum,q)
% function [paramBS, paramEst] = MF_BS_test_constancy_param(B1,B2,Block,Method,T_S,Alpha,Fun,j1,j2,wtype ,Jflag,MomNul,methodWT,Cum,q)
%
% Write parameters for MF_test_constancy in structures
%
% Herwig Wendt, Lyon, 2006 - 2008


paramBS.Nresamp1=B1;
paramBS.Nresamp2=B2 ;
paramBS.blocklength=Block;
paramBS.Method=Method;
paramBS.T_S=T_S;
paramBS.Alpha=Alpha;

paramEst.EstFun=Fun;
paramEst.j1=j1;
paramEst.j2=j2;
paramEst.wtype=wtype;
paramEst.Jflag=Jflag;
paramEst.MomNul=MomNul;
paramEst.methodWT=methodWT;
paramEst.Cum=Cum;
paramEst.q=q;