function out = PCORboot_ph( datamat, x_in, x_out )

if nargin<3
    x_out=[];
end

econd = eps;
NBOOT = 2000;

    disp('partial correlation, bootstrapped...');

    % parameters
    n    = size(datamat,2);
    k    = size(x_in, 2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat,1), k, NBOOT );
    for(bsr=1:NBOOT)
        [bsr NBOOT],
        list = ceil(n*rand(n,1));
        % boot sample
        D = datamat(:,list);
        y = x_in(list,:);
        %
        if ~isempty(x_out)
            % other boot sample
            z = [ones(n,1) x_out(list,:)]; % aug w/ intercept
            % residualizeing on z
            bd = D * (z / (z'*z));
            D  = D - (bd * z');
            by = y' * (z / (z'*z));
            y  = y - (z * by');
        else
            D = bsxfun(@minus,D,mean(D,2));
            y = bsxfun(@minus,y,mean(y,1));
        end
        % unit norming
        D = bsxfun(@rdivide,D,sqrt(sum(D.^2,2))+eps);
        y = bsxfun(@rdivide,y,sqrt(sum(y.^2,1))+eps);
        % correlating - just take inner product now
        pcrmat(:,:,bsr) = (D * y); % inner product
    end
    out.bsr     = mean(pcrmat,3)./std(pcrmat,0,3);
    out.bsr_p   = 2*min(cat(3,sum(pcrmat>0,3), sum(pcrmat<0,3)),[],3)./bsr;
   
    out.testname = 'pcorr_bootstrap';
