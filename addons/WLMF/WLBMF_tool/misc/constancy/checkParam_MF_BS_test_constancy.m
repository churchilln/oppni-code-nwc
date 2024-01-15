% checkParam_MF_BS_test_constancy.m
%
% Check if all necessary parameters for 
% WriteParam_MF_test_constancy.m are defined
% In case of, give them default values
%
% Herwig Wendt, Lyon, 2006 - 2008
 
try 
    B1; 
catch; 
    B1=49; disp(['Set B1=',num2str(B1),' (default)']); end;
try 
    B2; 
catch; 
    B2=49; disp(['Set B2=',num2str(B2),' (default)']); end; 
try 
    methodWT; 
catch; 
    methodWT=[2]; disp(['Set methodWT=',num2str(methodWT),' (default)']);  end;
try 
    M; 
catch; 
    M=2; disp(['Set M=',num2str(M),' (default)']);  end;
try 
    T_S; 
catch; 
    T_S=0; disp(['Set T_S=',num2str(T_S),' (default)']);  end;
try 
    MomNul; 
catch; 
    MomNul=3; disp(['Set MomNul=',num2str(MomNul),' (default)']);  end;
try 
    Norm; 
catch; 
    Norm=1; disp(['Set Norm=',num2str(Norm),' (default)']);  end;
try 
    verbose; 
catch; 
    verbose=2; disp(['Set verbose=',num2str(verbose),' (default)']);  end;
try 
    FigNum; 
catch; 
    FigNum=0; disp(['Set FigNum=',num2str(FigNum),' (default)']);  end;
try 
    Alpha; 
catch; 
    Alpha=0.1; disp(['Set Alpha=',num2str(Alpha),' (default)']);  end;
try 
    Method; 
catch; 
    Method=3; disp(['Set Method=',num2str(Method),' (default)']);  end;
try 
    q; 
catch; 
    q=2; disp(['Set q=',num2str(q),' (default)']);  end;
try 
    Cum; 
catch; 
    Cum=2; disp(['Set Cum=',num2str(Cum),' (default)']);  end;
try 
    Jflag; 
catch; 
    Jflag=1; disp(['Set Jflag=',num2str(Jflag),' (default)']);  end;
try 
    Fun; 
catch; 
    Fun=101; disp(['Set Fun=',num2str(Fun),' (default)']);  end;
try 
    wtype; 
catch; 
    wtype=1; disp(['Set wtype=',num2str(wtype),' (default)']);  end;
try 
    j1; 
catch; 
    j1=3; disp(['Set j1=',num2str(j1),' (default)']);  end;

if ~T_S;
    try 
        Block; 
    catch; 
        Block=2*MomNul; end;
end;

try N>0;
    try 
        j2; 
    catch; 
        j2=log2(N)+1-log2(M)-5; disp(['Set j2=',num2str(j2),' (default)']);  end;
    try 
        Block; 
    catch; 
        Block=floor(N/32); end;
catch; 
    error('Cannot determine j2 or Blocklength: Need data length');
end;