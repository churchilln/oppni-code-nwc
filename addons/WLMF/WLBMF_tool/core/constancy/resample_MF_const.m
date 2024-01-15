function [Xstar]=resample_MF_const(X,Block, B1, B2);
% function [Xstar]=resample_MF_const(X,Block, B1, B2);
%
% Herwig Wendt, Lyon, 2006 - 2008

if nargin==3; B2=1; end

N=length(X);
if Block>1
    % initialize blocks of indices
    bx=0:N-1;
    try % fast but memory intensive
        bx=0:N-1;
        BX=repmat(bx,Block,1);
        addit=repmat([0:Block-1]', 1,N);
        BX=BX+addit;
        BX=mod(BX,N)+1;
    catch % slow but memory save
        BX=[];bx=1:N;
        for bl=1:Block
            BX=[BX; bx];
            bx=[bx(end) bx(1:end-1)];   % create overlapping blocks
        end
    end
end

CIRCULAR=1;
if CIRCULAR
    N_range=N;
else
    N_range=max(N-Block+1,1);
end
N_BS=ceil(N/Block);
% N_BS=max(1,floor(N/Block)); % CHANGED: 23/08/2007
N_resample=N_BS*Block;

for b=1:B1  
    if Block==1
        index = fix(rand(1,N)*N)+1 ;
        Xstar{b}.t=X(index);
    else
        index = fix(rand(1,N_BS)*N_range)+1 ;
        tempid=reshape(BX(:,index),1,[]);
        Xstar{b}.t=X(tempid);
    end   
    
    if B2>1
    for bb=1:B2     
        if Block==1
            index = fix(rand(1,N)*N)+1 ;
            Xstar{b}.T(bb,:)=Xstar{b}(index);
        else  % block bootstrap
            if 0
                index2 = fix(rand(1,N_BS)*N_range)+1 ;
                tempid=reshape(BX(:,index2),1,[]);
                Xstar{b}.T(bb,:)=Xstar{b}.t(tempid);
            else
                % ---> new: resample from blocks
                index2 = fix(rand(1,N_BS)*N_BS)+1 ;
                tempid=reshape(BX(:,index(index2)),1,[]);
                Xstar{b}.T(bb,:)=X(tempid);
            end            
        end            
    end
    end
end



