function [paramEST, paramBS, paramTest]=MF_BS_tool_param(MomNul,gamint,j1,j2,wtype,Fun,Cum,q,B1,B2,Block,Method,Alpha,Jflag,Type,Tnull,T_S,CI,TEST);
% function [paramEST, paramBS, paramTest]=MFA_BS_light_param(MomNul,Norm,j1,j2,wtype,Fun,Cum,q,B1,B2,Block,Method,Alpha,Jflag,Type,Tnull,T_S);
% 
% Herwig Wendt, Lyon, 2006 - 2008

%%%% Estimation
paramEST.MomNul=MomNul;
paramEST.gamint=gamint;
paramEST.sym=0;
paramEST.j1=j1;
paramEST.j2=j2;
paramEST.wtype=wtype;
paramEST.Fun=Fun; % 0: only cp; 1: zeta(q), D(h), cp
paramEST.Cum=Cum;
paramEST.q=q;
%%%% Bootstrap
paramBS.Nresamp1=B1;
paramBS.Nresamp2=B2;
paramBS.blocklength=Block;
paramBS.Method=Method;
paramBS.Alpha=Alpha;
paramBS.Jflag=Jflag;
paramBS.T_S=T_S;
paramBS.CI=CI;
paramBS.TEST=TEST;
paramTest.Type=Type;
paramTest.Tnull=Tnull;
