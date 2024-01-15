function [coef, leaders, nj] = DxLx1d(data, Nwt, gamint,symm)

%
%  [coef, leaders, nj] = DxLx1d_newNan_wavelab(data, Nwt [,gamint,sym])
%
% compute Dyadic wavelet and leader coefficients  
%
% gamint -- fractional integration of order gamint is applied
% sym=0  -- use Daubechies Wavelet of order Nwt.
% sym=1  -- use Daubechies Symmetrized Wavelet of order Nwt.
%
%%%
% SR, ens-Lyon, 11/2013
%


Norm=1;

%-- paramater for  fractional integration
if nargin<3; gamint=0; end;
if isempty(gamint); gamint=0; end;

if nargin<4; symm=0; end;
%-- Initialize the wavelet filters
n = length(data) ;         % data length
if symm == 0 % Daubechies Wavelet
    h = rlistcoefdaub(Nwt) ;   % filter
    nl = length(h) ;           % length of filter, store to manage edge effect later
    gg1 = -1*(-1).^(1:nl).*h ; % wavelet filter
    hh1=fliplr(h);             % scaling filter
    % parameter for the centering of the wavelet
    x0=2; 
    x0Appro=2*Nwt; 

else % Daubechies Symmetrized Wavelet
    Nwt=abs(Nwt);
    h = rlistcoefdaub(Nwt) ;   % filter
    nl = length(h) ;           % length of filter, store to manage edge effect later
    gg1c = -1*(-1).^(1:nl).*h ; % wavelet filter
    hh1c=fliplr(h);             % scaling filter
    tmp=conv(h,hh1c)/sqrt(2);
    nu=-2*Nwt+1;
%     [v,nv] = tildeurflippeur(-tmp,nu);
    [v,nv] = tildeurflippeur(tmp,nu);
    
    [hh1,nh1] = flippeur(tmp,nu) ; nl=length(hh1);
    [gg1,ng1] = flippeur(v,nv) ;   

    
    % parameter for the centering of the wavelet
    x0=2*Nwt; 
    x0Appro=2*Nwt; 

end

%--- Predict the max # of octaves available given Nwt, and take the min with
nbvoies= fix( log2(length(data)) );
% nbvoies = min( fix(log2(n/(2*Nwt+1)))  , nbvoies); %   safer, casadestime having problems
nbvoies = min( fix(log2(n/(nl+1)))  , nbvoies);

%--- Compute the WT, calculate statistics
approW=data;
for j=1:nbvoies         % Loop Scales
    
    %-- Phase 1a: get the wavelet coefficients/appro at this scale
    njtemp = length(approW) ;
    conv_gg1 = conv(approW,gg1) ; conv_gg1(isnan(conv_gg1))=Inf;
    conv_hh1 = conv(approW,hh1) ; conv_hh1(isnan(conv_hh1))=Inf;
      
    
    % border effect
    fp=nl-1; % index of first good value
    lp=njtemp; % index of last good value
    
    % replace border with Inf     
    conv_gg1(1:fp-1)=Inf;
    conv_gg1(lp+1:end)=Inf;
    
    conv_hh1(1:fp-1)=Inf;
    conv_hh1(lp+1:end)=Inf;

    %-- centering and decimation 
    approW=conv_hh1((1:2:njtemp)+x0Appro-1);
    decime=conv_gg1((1:2:njtemp)+x0-1);     
             
    %-- passage Norme L1
    sig=sign(decime);
    AbsdqkW = abs(decime)*2^(j/2)/2^(j/Norm);   
     
    %-- max before integration
    lesi=find(isfinite(AbsdqkW));
    
    %% HW
    if length(lesi)<3; return; end
    
    coef(j).supcoefnoint=max(AbsdqkW(lesi));
    
    %-- fractional integration and max
    AbsdqkW = AbsdqkW*2^(gamint*j);
    [coef(j).supcoef, coef(j).supcoefid]=max(AbsdqkW(lesi));
    
    %-- store results
    coef(j).allvalue=AbsdqkW;
    coef(j).allxpos=1:2^j:n;
    coef(j).value=AbsdqkW(lesi);
    coef(j).value_noabs=AbsdqkW(lesi).*sig(lesi);
    coef(j).xpos=coef(j).allxpos(lesi);
    coef(j).gamma=gamint;
    nj.W(j) = length(lesi);            % number of coefficients
    
   
    % for leaders 
    if j ==1
        %-- compute and store leaders sans voisin  
        leaders(j).sans_voisin.allvalue = AbsdqkW;
        leaders(j).sans_voisin.allxpos=coef(j).allxpos;
        
        %% useless
        leaders(j).sans_voisin.value = AbsdqkW(lesi);
        leaders(j).sans_voisin.xpos=coef(j).allxpos(lesi);
        %%
        
        %-- compute leaders  
        nc=floor(length(AbsdqkW)/2);
        voisin = max([AbsdqkW(1:end-2); AbsdqkW(2:end-1);AbsdqkW(3:end)]);
       
        
        %-- store results
        ifini=find(isfinite(voisin));
        leaders(j).allvalue=voisin;        
        leaders(j).allxpos=coef(j).allxpos(2:end-1);    
        leaders(j).value=voisin(ifini);
        leaders(j).xpos=leaders(j).allxpos(ifini);
        
        
        nj.L(j) = length(leaders(j).value); % number of leaders
 
    else
        %-- compute and store leaders sans voisin
        nc=floor(length(leaders(j-1).sans_voisin.allvalue)/2);
        leaders(j).sans_voisin.allvalue = max([coef(j).allvalue(1:nc); leaders(j-1).sans_voisin.allvalue(1:2:2*nc) ; leaders(j-1).sans_voisin.allvalue(2:2:2*nc)]);      
        leaders(j).sans_voisin.allxpos = coef(j).allxpos(1:nc);
  
        %% useless
        ifini=find(isfinite(leaders(j).sans_voisin.allvalue));
        leaders(j).sans_voisin.value = leaders(j).sans_voisin.allvalue(ifini);
        leaders(j).sans_voisin.xpos = leaders(j).sans_voisin.allxpos(ifini);
        %% 
        
        %-- compute leaders
        voisin = max([leaders(j).sans_voisin.allvalue(1:end-2);leaders(j).sans_voisin.allvalue(2:end-1);leaders(j).sans_voisin.allvalue(3:end)]);
        
        %-- store results
        ifini=find(isfinite(voisin));
        leaders(j).allvalue=voisin;
        leaders(j).allxpos=coef(j).allxpos(2:nc-1);
        leaders(j).value=voisin(ifini);
        leaders(j).xpos=leaders(j).allxpos(ifini);
         
        nj.L(j) = length(leaders(j).value); % number of leaders
        
    end
 
     leaders(j).gamma=gamint;    
     [leaders(j).mincoef, leaders(j).mincoefid]=min(leaders(j).value);
     [leaders(j).supcoefL, leaders(j).supcoefidL]=max(leaders(j).value);
     leaders(j).supcoefnoint=coef(j).supcoefnoint; leaders(j).supcoef=coef(j).supcoef; leaders(j).supcoefid=coef(j).supcoefid;
     coef(j).mincoef=leaders(j).mincoef; coef(j).mincoefid=leaders(j).mincoefid;
     coef(j).supcoefL=leaders(j).supcoefL; coef(j).supcoefidL=leaders(j).supcoefidL;
         

end

