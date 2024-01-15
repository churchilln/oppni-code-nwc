% quickly get max scale...
%
NT = 247;

M=1;%number of windows
data=randn(NT,1);
MomNul=3;

% --- Calculate Coef Lead
[coef, lead, nj] = DxLx1d(data', MomNul);
j2max=find(nj.L>=4*M);
j2max=j2max(end); 
