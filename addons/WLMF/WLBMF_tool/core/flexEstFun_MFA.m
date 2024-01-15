function [varargout] = flexEstFun_MFA(Absdqk, param)
% function [varargout] = flexEstFun_MFA(Absdqk, param)
%
% Calculates Structure Functions (SF) and quantities for D(q), h(q) for a vector q (length(q)>1)
% Log Cumulants 1-5
%
% -- INPUT
%   Absdqk  -   |d_X(j,.)| or L_X(j,.)
%   param   -   structure with parameters
%       param.q      : vector with moments q to calculate
%       param.Cum    : highest order of Cumulant to calculate 
%       param.EstFun : what is to calculate
%                       number xyz 
%                       x : SF zeta(q) [0 or 1]
%                       y : SF D(h) [0 or 1]
%                       y : SF Cp [0 or 1]
%                       e.g. param.EstFun=001 calculates only Cp, 
%                            param.EstFun=110 calculates SF zeta(q) and D(h)
%
% Herwig Wendt, Lyon, 2006 - 2008

% SORTOUT=1; % sort out small values
SORTOUT=0;
thresh=eps;%1e-10;
if SORTOUT
    Absdqk=Absdqk((Absdqk>=thresh));
end
% try 
%     param.Fun;
%     Fun=param.Fun;
% catch;
%    Fun=bin2dec(num2str(param.EstFun));
    Fun=0;
    EstFun=param.EstFun;
    if EstFun>=100; Fun=Fun+4; EstFun=EstFun-100; end
    if EstFun>=10; Fun=Fun+2; EstFun=EstFun-10; end
    Fun=Fun+EstFun;
%end

% % structure functions and Moments. f_Dq and f_hq
% for  kq=1:length(q)    %  Loop on the values of q
%     Elogmuqj(kq) = log2(mean(Absdqk.^q(kq)));    % Calculate  log[ S_q(j) ]    ( in L1 language )
%     detail_kq = Absdqk.^q(kq); sum_detail_kq = sum(detail_kq);
%     f_Dq(kq) = sum(detail_kq .* log2(detail_kq / sum_detail_kq )) / sum_detail_kq + log2(length(Absdqk));
%     f_hq(kq) = sum(detail_kq .* log2(Absdqk)) / sum_detail_kq;
% end
% Matrix Calculation -> fast
argcount=1;
if Fun~=1; % if not only Cumulants
    q=param.q;
    %% structure functions and Moments. f_Dq and f_hq WITHOUT LOOP
    %% This version is faster than the loop for length(q)>3 and length(Absdqk)>10 (~0.7 times the loop)
    lenA=length(Absdqk);
    lenq=length(q);
    S=repmat(Absdqk,lenq,1);
    Q=repmat(q', 1, lenA);
    detail_kq=S.^Q;
    % zeta(q): Fun = 4, 5, 6, 7
    if Fun>=4; 
        Elogmuqj=(log2(mean(detail_kq,2)))';
        varargout{argcount}=Elogmuqj; argcount=argcount+1;
    end
    % D(h): Fun = 2, 3, 6, 7
    if (Fun~=4) && (Fun~=5) 
        sum_detail_kq = sum(detail_kq,2);
        f_Dq=(sum(detail_kq .* log2(detail_kq ./ repmat(sum_detail_kq,1,lenA) ),2)./ sum_detail_kq + log2(lenA))';
        f_hq=(sum(detail_kq .* log2(S),2) ./ sum_detail_kq)';
        varargout{argcount}=f_Dq; argcount=argcount+1;
        varargout{argcount}=f_hq; argcount=argcount+1;
    end
else
    Elogmuqj=NaN;
    f_Dq=NaN;
    f_hq=NaN;
end
% Cumulants: Fun = 1, 3, 5, 7
% if rem(Fun,2)==1
%     Cum=param.Cum;
%     %-- Cumulants
%     Absdqk=log(Absdqk);
%     Cp(1) = mean(Absdqk) ;
%     if Cum>1
%         Cp(2) = mean(Absdqk.^2) - Cp(1)^2 ;
%         if Cum>2
%             Cp(3) = mean(Absdqk.^3) - 3*Cp(2)*Cp(1) - Cp(1)^3 ;
%             if Cum>3
%                 Cp(4) = mean(Absdqk.^4) - 4*Cp(3)*Cp(1) - 3*Cp(2)^2 - 6*Cp(2)*Cp(1)^2 - Cp(1)^4;
%                 if Cum>4
%                     Cp(5) = mean(Absdqk.^5) - 5*Cp(4)*Cp(1) - 10*Cp(3)*Cp(2) - 10*Cp(3)*Cp(1)^2 - 15*Cp(2)^2*Cp(1) - 10*Cp(2)*Cp(1)^3 - Cp(1)^5 ;
%                 end;
%             end;
%         end;
%     end
%     varargout{argcount}=Cp;
% end
%% HW 26/05/2010
if rem(Fun,2)==1
    Cum=param.Cum;
    %-- Cumulants
    Absdqk0=log(Absdqk);
    Cp(1) = mean(Absdqk0) ;
    if Cum>1
        Absdqk=Absdqk0.^2;
        Cp(2) = mean(Absdqk) - Cp(1)^2 ;
        if Cum>2
            Absdqk=Absdqk.*Absdqk0;
            Cp(3) = mean(Absdqk) - 3*Cp(2)*Cp(1) - Cp(1)^3 ;
            if Cum>3
                Absdqk=Absdqk.*Absdqk0;
                Cp(4) = mean(Absdqk) - 4*Cp(3)*Cp(1) - 3*Cp(2)^2 - 6*Cp(2)*Cp(1)^2 - Cp(1)^4;
                if Cum>4
                    Absdqk=Absdqk.*Absdqk0;
                    Cp(5) = mean(Absdqk) - 5*Cp(4)*Cp(1) - 10*Cp(3)*Cp(2) - 10*Cp(3)*Cp(1)^2 - 15*Cp(2)^2*Cp(1) - 10*Cp(2)*Cp(1)^3 - Cp(1)^5 ;
                end;
            end;
        end;
    end
    varargout{argcount}=Cp;
end
