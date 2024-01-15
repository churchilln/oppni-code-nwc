function out = Bootstrap_ph( datamat, design, estim, verb )
%
% . Bootstrapped analysis, compares 1-group rel. 0, or between 2 groups
%

if nargin<4
    verb=1;
end

NBOOT = 1000;

if(nargin<3) estim='mean'; end

if    ( isempty(design) || numel(unique(design))==1 )
    
    disp('Bootstrap, 1-sample...');
    
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat,1), NBOOT );
    for(bsr=1:NBOOT)
        if verb>0
        [bsr NBOOT],
        end

        list = ceil(n*rand(n,1));
        
        switch estim % select bootstrapping statistic
            case 'mean'
                bsrmat(:,bsr) = mean(datamat(:,list),2);
            case 'stdev'
                bsrmat(:,bsr) = std(datamat(:,list),0,2);
            case 'robmean'
                bsrmat(:,bsr) = location_mest( datamat(:,list), 1.5, 50 );
            case 'robstdev'
                bsrmat(:,bsr) = 1.49*median( abs(bsxfun(@minus,datamat(:,list),median(datamat(:,list),2))),2);
            case 'safemean'
                bsrmat(:,bsr) = mean(datamat(:,list),2,'omitnan');
            case 'safedev'
                bsrmat(:,bsr) = std(datamat(:,list),0,2,'omitnan');
            otherwise
                error(['undefined parameter input:', estim]);
        end
    end
    out.bsr     = mean(bsrmat,2)./std(bsrmat,0,2);
    out.bsr_p   = 2*min([mean(bsrmat>0,2), mean(bsrmat<0,2)],[],2);
    
    out.testname = 'bootstrap_1samp';

elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==1 )
    
    ix = unique(design);
    datamat1 = datamat(:,design==ix(1));
    datamat2 = datamat(:,design==ix(2));
    disp('Bootstrap, 2-sample (unpaired)...');

    % parameters
    n1 = size(datamat1,2);
    n2 = size(datamat2,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat1,1), NBOOT );
    for(bsr=1:NBOOT)
        
        if verb>0
        [bsr NBOOT],
        end

        list1 = ceil(n1*rand(n1,1));
        list2 = ceil(n2*rand(n2,1));  
        
        switch estim % select bootstrapping statistic
            case 'mean'
                bsrmat(:,bsr) = mean(datamat2(:,list2),2) - mean(datamat1(:,list1),2);
            case 'stdev'
                bsrmat(:,bsr) = std(datamat2(:,list2),0,2) - std(datamat1(:,list1),0,2);
            case 'robmean'
                bsrmat(:,bsr) = location_mest( datamat2(:,list2), 1.5, 50 ) - location_mest( datamat1(:,list1), 1.5, 50 );
            case 'robstdev'
                bsrmat(:,bsr) = 1.49* ( median( abs(bsxfun(@minus,datamat2(:,list2),median(datamat2(:,list2),2))),2)   -   median( abs(bsxfun(@minus,datamat1(:,list1),median(datamat1(:,list1),2))),2) );
            case 'safemean'
                bsrmat(:,bsr) = mean(datamat2(:,list2),2,'omitnan') - mean(datamat1(:,list1),2,'omitnan');
            case 'safedev'
                bsrmat(:,bsr) = std(datamat2(:,list2),0,2,'omitnan') - std(datamat1(:,list1),0,2,'omitnan');
            otherwise
                error(['undefined parameter input:', estim]);
        end
    end
    out.bsr     = mean(bsrmat,2)./std(bsrmat,0,2);
    out.bsr_p   = 2*min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;    

    out.testname = 'bootstrap_2samp_unpair';

elseif( numel(unique(design))>1 && size(design,1)>1 && size(design,2)==2 )
    
    ix = unique(design(:,1));
    datamat1 = datamat(:,design(:,1)==ix(1)); 
    datamat2 = datamat(:,design(:,1)==ix(2)); 
    
    clear datamat; %% wipe out old datamat
   
    des1     = design(design(:,1)==ix(1),2);
    des2     = design(design(:,1)==ix(2),2);
    
    if( size(datamat1,2) ~= size(datamat2,2) ) 
        error('unequal splits, cannot match');
    else
        datamat = zeros( size(datamat1) );
    end
    
    tmp=zeros(size(datamat2));
    for(i=1:length(des1))
        tmp(:,i) = datamat2(:,des2==des1(i)); 
    end
    datamat2 = tmp; clear tmp; % reorder to match datamat1
    
    disp('Bootstrap, 2-sample (paired)...');
    
    % parameters
    n    = size(datamat,2);
    disp('running resampling...');
    bsrmat = zeros( size(datamat1,1), NBOOT );
    for(bsr=1:NBOOT)
        
        if verb>0
        [bsr NBOOT],
        end
        
        list = ceil(n*rand(n,1));
        
        switch estim % select bootstrapping statistic
            case 'mean'
                bsrmat(:,bsr) = mean(datamat2(:,list),2) - mean(datamat1(:,list),2);
            case 'stdev'
                bsrmat(:,bsr) = std(datamat2(:,list),0,2) - std(datamat1(:,list),0,2);
            case 'robmean'
                bsrmat(:,bsr) = location_mest( datamat2(:,list), 1.5, 50 ) - location_mest( datamat1(:,list), 1.5, 50 );
            case 'robstdev'
                bsrmat(:,bsr) = 1.49* ( median( abs(bsxfun(@minus,datamat2(:,list),median(datamat2(:,list),2))),2)   -   median( abs(bsxfun(@minus,datamat1(:,list),median(datamat1(:,list),2))),2) );
            case 'safemean'
                bsrmat(:,bsr) = mean(datamat2(:,list),2,'omitnan') - mean(datamat1(:,list),2,'omitnan');
            case 'safedev'
                bsrmat(:,bsr) = std(datamat2(:,list),0,2,'omitnan') - std(datamat1(:,list),0,2,'omitnan');
            otherwise
                error(['undefined parameter input:', estim]);
        end
    end
    out.bsr     = mean(bsrmat,2)./std(bsrmat,0,2);
    out.bsr_p   = 2*min([sum(bsrmat>0,2) sum(bsrmat<0,2)],[],2)./bsr;    

    out.testname = 'bootstrap_2samp_pair';

end

