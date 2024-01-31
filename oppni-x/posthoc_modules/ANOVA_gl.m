function out = ANOVA_gl( datamat, design )

if( size(design,2)==1 ) %% 1-way anova

    %% convert into cell array for different classes
    ix=unique(design);
    for(i=1:length(ix)) tmpcell{i} = datamat(:,design==i); end
    datamat = tmpcell; clear tmpcell;

    %%
    des =[]; kk = length(datamat);
    Xmat=[];
    for(i=1:kk) 
        des = [des; i*ones(size(datamat{i},2),1)]; 
        Xmat=[Xmat datamat{i}];
        nn(i) = size(datamat{i},2); 
    end
    % grand mean
    GM = mean(Xmat,2);
    % total sum of squares
    SST= sum( bsxfun(@minus, Xmat,GM).^2, 2);
    SSW= 0;
    for(i=1:kk) 
        % treatment means
        TM(:,i) = mean(Xmat(:,des==i),2);
        % add to within-class variance
        SSW = SSW + sum(bsxfun(@minus,Xmat(:,des==i),TM(:,i)).^2,2);
    end
    % between-class variance
    SSB= (bsxfun(@minus, TM,GM).^2)*nn(:);
    % degrees of freedom
    df_b = kk-1;
    df_w = sum(nn) - kk;
    % mean squared errors
    MSB= SSB./df_b;
    MSW= SSW./df_w;

    % group effect
    out.fstat_G    = MSB./MSW;
    out.fstat_G_p  = 1-fcdf(out.fstat_G,(kk-1),sum(nn) - kk);

    out.testname = 'anova1';

elseif( size(design,2)==2 ) %% 2-way anova

    NF=2;
    for(f=1:2) NL(f) = numel( unique(design(:,f)) ); end
    for(i=1:NL(1))
    for(j=1:NL(2))
        rr(i,j) = sum( design(:,1)==i & design(:,2)==j );
    end
    end
    for(f=1:2) 
    for(i=1:NL(f))
        rf{f}(i)= sum( design(:,f)==i );
    end
    end

    ssw = 0;
    for(i=1:NL(1))
    for(j=1:NL(2))

        yijk = datamat(:, design(:,1)==i & design(:,2)==j );
        yij  = mean(yijk,2);
        ssw  = ssw+ sum( (yijk-mean(yij)).^2,2 );
    end
    end
    dfw = sum(rr(:)) - prod(NL);

    y_gm = mean(datamat,2);

    ssg1=0;
    for(i=1:NL(1))
        ssg1 = ssg1 + rf{1}(i)*( mean(datamat(:,design(:,1)==i),2) - mean(datamat,2) ).^2;
    end
    dfg1 = NL(1)-1;

    ssg2=0;
    for(j=1:NL(2))
        ssg2 = ssg2 + rf{2}(j)*( mean(datamat(:,design(:,2)==j),2) - mean(datamat,2) ).^2;
    end
    dfg2 = NL(2)-1;

    ssi=0;
    for(i=1:NL(1))
    for(j=1:NL(2))
        yij = mean( datamat( :, design(:,1)==i & design(:,2)==j ),2 );
        yi  = mean( datamat( :, design(:,1)==i ),2 );
        yj  = mean( datamat( :, design(:,2)==j ),2 );

        ssi = ssi + rr(i,j)*( yij - yi - yj +y_gm ).^2;
    end
    end
    dfi = (NL(1)-1)*(NL(2)-1);

    ms=[ ssg1./dfg1, ssg2./dfg2, ssi./dfi, ssw./dfw ];
    ff=bsxfun(@rdivide, ms(:,1:3), ms(:,4) );
    PP = [ 1-fcdf(ff(:,1),dfg1,dfw), 1-fcdf(ff(:,2),dfg2,dfw), 1-fcdf(ff(:,3),dfi,dfw)];

    %% statistics - treat1 / treat2 / interaction

    % group1 effect
    out.fstat_G1    = ff(:,1);
    out.fstat_G1_p  = PP(:,1);
    % group2 effect
    out.fstat_G2    = ff(:,2);
    out.fstat_G2_p  = PP(:,2);
    % interaction effect
    out.fstat_G12   = ff(:,3);
    out.fstat_G12_p = PP(:,3);

    out.testname = 'anova2';

else
    disp('can only run 1-way or 2-way anova. Only 1 or 2 design factors allowed.');
end