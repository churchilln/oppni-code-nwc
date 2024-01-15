% checkParam_MF_BS_tool.m
%
% Check if all necessary parameters for 
% MFA_BS_light_param.m are defined
% In case of, give them default values
%
% Herwig Wendt, Lyon, 2006 - 2008

try 
    CI;
catch;
    CI=0; end
try 
    TEST;
catch;
    TEST=0; end
try 
    Type;
catch;
    Type=0; TEST=0; end
try 
    Tnull;
catch;
    Tnull=0; TEST=0; end
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
    T_S; 
catch; 
    T_S=0; disp(['Set T_S=',num2str(T_S),' (default)']);  end;
try 
    MomNul; 
catch; 
    MomNul=3; disp(['Set MomNul=',num2str(MomNul),' (default)']);  end;
try 
    gamint; 
catch; 
    gamint=0; disp(['Set gamint=',num2str(gamint),' (default)']);  end;
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
    sym; if sym~=0&sym~=1; sym=1; disp(['Set sym=',num2str(sym)]); end
catch; 
    sym=0; disp(['Set sym=',num2str(sym),' (default)']);  end;
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

try isempty(j2);
    if j2<=j1+2;  j2=j1+2; end
catch;
    try N>0;
        try
            j2;
        catch;
            j2=log2(N)-4; disp(['Set j2=',num2str(j2),' (default)']);  end;
        try
            Block;
        catch;
            Block=floor(N/32); end;
    catch;
        error('Cannot determine j2 or Blocklength: Need data length');
    end;
end


