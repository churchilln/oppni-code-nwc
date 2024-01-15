function out = orientation_checker( A,varargin )

% quick function - check for neurologic/radiologic orientation
% and whether 2 images have the same orientation

if    ( ischar(A)   && exist(A,'file')   ) VA=load_untouch_niiz(A);
elseif( isstruct(A) && isfield(A,'hdr')  ) VA=A;
elseif( isstruct(A) && isfield(A,'hist') ) VA.hdr=A;
else,  error('file cannot be found / not recognized nii struct');
end

%%
out.orient_A(1,1) = sign(VA.hdr.hist.srow_x(1));
out.orient_A(2,1) = sign(VA.hdr.hist.srow_y(2));
out.orient_A(3,1) = sign(VA.hdr.hist.srow_z(3));

str=[];
if(out.orient_A(1)>0) str = 'L->R';
else                  str = 'R->L';
end
if(out.orient_A(2)>0) str = [str,',P->A'];
else                  str = [str,',A->P'];
end
if(out.orient_A(2)>0) str = [str,',I->S'];
else                  str = [str,',S->I'];
end
out.orientStr_A = str;

out.neuro_A  = prod( out.orient_A == [1 1 1]' );
out.radio_A  = prod( out.orient_A == [-1 1 1]' );

if(nargin>1)
  
    B = varargin{1};

    if    ( ischar(B)   && exist(B,'file')   ) VB=load_untouch_niiz(B);
    elseif( isstruct(B) && isfield(B,'hdr')  ) VB=B;
    elseif( isstruct(B) && isfield(B,'hist') ) VB.hdr=B;
    else,  error('file cannot be found / not recognized nii struct');
    end

    %%
    out.orient_B(1,1) = sign(VB.hdr.hist.srow_x(1));
    out.orient_B(2,1) = sign(VB.hdr.hist.srow_y(2));
    out.orient_B(3,1) = sign(VB.hdr.hist.srow_z(3));
    
    str=[];
    if(out.orient_A(1)>0) str = 'L->R';
    else                  str = 'R->L';
    end
    if(out.orient_A(2)>0) str = [str,',P->A'];
    else                  str = [str,',A->P'];
    end
    if(out.orient_A(2)>0) str = [str,',I->S'];
    else                  str = [str,',S->I'];
    end
    out.orientStr_B = str;

    out.neuro_B  = prod( out.orient_B == [1 1 1]' );
    out.radio_B  = prod( out.orient_B == [-1 1 1]' );

    out.consist_AB = prod(out.orient_A == out.orient_B);
    
end

