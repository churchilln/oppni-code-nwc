function out = glca_kde( I, v_range, d_set, model, b_width, makefigs, kd_nsamp )
%
% =========================================================================
% GLCA_KDE: script that takes in 3d image array and generates a set of
% Haralick texture features
% =========================================================================
%
% Syntax:
%            out = glca_kde( I, v_range, d_set, model, b_width, makefigs, kd_nsamp )
%
% Input:
%
%           I : (mandatory) 3d array of input values to analyze. elements with NaN values are excluded from analysis
%     v_range : (optional) specifies the min-max range to rescale values over. can be either 2d vector of numeric values or
%               A special string argument:
%                      'PctRange-<N>', where <N> denotes an integer between 1 and 50
%                      'StDev-<N>', where <N> denotes a number between 2 and 20
%                      'RobStDev-<N>', where <N> denotes a number between 2 and 20
%               default is 'PctRange-10'
%       d_set : (optional) set of 3d directions to measure texture. Specified as a
%               stack of row vectors, e.g., d_set = [ [1 0 0]; [0 1 0]; ...];
%               default is the full set of 13 indepdendent directions
%       model : (mandatory) specifies either histogram-based ('HDE') or kernel density based ('KDE') analysis
%               The HDE approach is *much* faster to compute but potentially less robust
%     b_width : (mandatory) bandwidth parameter specifying coarseness of histogram / kernel plot
%               can either be a single numeric value <1 (directly estimate features) or 
%               A special string argument:
%                       'estim', just estimates the optimal bandwidth
%                       'optim', selects optimal bandwidth and gets corresponding texture features
%    makefigs : (optional) special "display mode" that just generates plots and terminates early. 
%               0=off, 1=on. Default is off
%    kd_nsamp : (optional) specifies density of grid to sample on for KDE. Only change it if you are encountering errors.
%               Default should be 200
%
% Output:
%
%          out.grid_info : information about grid sampling density if running KDE 
%              out.TFset : (bw-values x directions) array of "truncation factor" values (optim) 
%     out.Copt_err(id,1) : (bw-values x 1) vector of optimal CV errors
%     out.Copt_idx(id,1) : (bw-values x 1) vector of optimal CV indices
%     out.Copt_bvl(id,1) : (bw-values x 1) vector of optimal bandwidths
%          out.Copt_bref : bandwidth chosen by plug-in estimator
%             out.CVfull : (bw-values x directions) array of cross-validation errors
%             out.SZfull : (bw-values x directions) array of bandwidth values
%       out.Bvl_for_fitt : best average bandwidth OR predefined value if provided 
%              out.TFopt : optimal bandwidth truncation factor
%         out.metrics_av : vector of average texture metrics
%         out.metrics_sd : vector of st-dev on texture metrics
%           out.bscalset : optimal/preselected bandwidth scale relative to plug-in estimator
%

scaled_noise = 0;
scaled_shift = 0;

%--- safeguard against discretization artifact
g = sort(I(isfinite(I)));
a = g(2:end)-g(1:end-1); 
a=a(a>eps);
amin=min(a); % smallest step-size detectable

% 1 = fit glca and get estimates / 2 = estimate optimal bandwidth
if ischar(b_width) 
    if strcmpi( b_width, 'estim')
        gl_mode = 2; % gets cv curve on bandwidth (no TEX analysis!)
    elseif strcmpi( b_width, 'optim')
        gl_mode = 3; % gets cfv curve on bandwidth, THEN tex-analysis on best
    else
        error('unrecognized bw arg')
    end
else
    gl_mode = 1; % tex-analysis using predefined bandwidth value
end
out.nelem = sum( isfinite(I(:)) );

vlist = I(isfinite(I));
if nargin<2 || isempty(v_range)
    vmin=min(vlist) - 0.10*range(vlist);
    vmax=max(vlist) + 0.10*range(vlist);

elseif ischar(v_range) && contains(v_range,'PctRange-')
    pct_cut = str2num(v_range(10:end));
    if isempty(pct_cut) || pct_cut<1 || pct_cut>50
        error('bad pct-cut');
    end
    vmin=min(vlist) - (pct_cut/100)*range(vlist);
    vmax=max(vlist) + (pct_cut/100)*range(vlist);

elseif ischar(v_range) && contains(v_range,'StDev-')
    std_cut = str2num(v_range(7:end));
    if isempty(std_cut) || std_cut<2 || std_cut>20
        error('bad std-cut');
    end
    vmin=mean(vlist) - std_cut*std(vlist);
    vmax=mean(vlist) + std_cut*std(vlist);

elseif ischar(v_range) && contains(v_range,'RobStDev-')
    std_cut = str2num(v_range(7:end));
    if isempty(std_cut) || std_cut<2 || std_cut>20
        error('bad std-cut');
    end
    vmin=median(vlist) - std_cut*iqr(vlist)/1.349;
    vmax=median(vlist) + std_cut*iqr(vlist)/1.349;
else
    vmin = v_range(1);
    vmax = v_range(2);
end
% quick check to make sure range bounds dont have major floor/ceiling effects
if mean(vlist<vmin) > 0.10
    error('big floor effect - more than 10% of data below your min-cut');
elseif mean(vlist<vmin) > 0.05
    warning('some floor effect - 5-10% of data below your min-cut');
elseif mean(vlist<vmin) > 0.01
    warning('possible floor effect - 1-5% of data below your min-cut');
end
if mean(vlist>vmax) > 0.10
    error('big ceiling effect - more than 10% of data above your max-cut');
elseif mean(vlist>vmax) > 0.05
    warning('some ceiling effect - 5-10% of data above your max-cut');
elseif mean(vlist>vmax) > 0.01
    warning('possible ceiling effect - 1-5% of data above your max-cut');
end

% now adjusting and bounding values on interval [0, 1]
I = (I-vmin)./(vmax-vmin);
I(I<0) = 0;
I(I>1) = 1;
% same for minstep, to check
amin2 = amin./(vmax-vmin);

% if no argin, will set...
if nargin<3 || isempty(d_set)
    d_set = [ 1  0  0;  0  1  0;  0  0  1; ... 
              1  1  0;  1  0  1;  0  1  1; ...
              1  1  1; ...
              1 -1  0;  1  0 -1;  0  1 -1; ...
             -1  1  1;  1 -1  1;  1  1 -1];
elseif numel(d_set)==1
    scal = d_set;
    d_set = [ 1  0  0;  0  1  0;  0  0  1; ... 
              1  1  0;  1  0  1;  0  1  1; ...
              1  1  1; ...
              1 -1  0;  1  0 -1;  0  1 -1; ...
             -1  1  1;  1 -1  1;  1  1 -1] .* scal;
end

if nargin<6
    makefigs = 0;
end
if nargin<7
    kd_nsamp = 'default';
end

% dimensions of array
nn=size(I);

if makefigs==1
    figure; 
    tiledlayout(2,4, 'Padding', 'none', 'TileSpacing', 'compact');
    for (i=1:7) 
        nexttile; imagesc( I(:,:,i),[0 1] ); colormap bone; set(gca,'XTick',[], 'YTick', []); 
    end
end

if strcmpi(model,'HDE')
    % no predefinitions here, since terms all depend on final choice of binwidth 
elseif strcmpi(model,'KDEU')
    % for gridding later...

    % note that for the kernel approach, del2 / mu0 / vr0 are predefined by
    % density of sampling grid (and must be pretty high)
    grange = [-0.5, 1.5];
    if isnumeric(kd_nsamp)
        NS=kd_nsamp;
    elseif strcmpi(kd_nsamp,'default')
        NS=300;
    end
    GRx  = repmat( linspace(grange(1),grange(2),NS), NS,1 );% (repmat( (1:NS),  NS,1 )-0.5)./NS;
    GRy  = repmat( linspace(grange(1),grange(2),NS)',1,NS );% (repmat( (1:NS)', 1,NS )-0.5)./NS;
    del2 = (GRx(1,2)-GRx(1,1))^2; 
    %
    mu_0 = mean(GRx(1,:)); % (2- -1)/2
    vr_0 = var( GRx(1,:));
    out.grid_info = [grange, NS, del2];

elseif strcmpi(model,'KDE')
    % for gridding later...

    % note that for the kernel approach, del2 / mu0 / vr0 are predefined by
    % density of sampling grid (and must be pretty high)
    grange = [0.0, 1.0];
    if isnumeric(kd_nsamp)
        NS=kd_nsamp;
    elseif strcmpi(kd_nsamp,'default')
        NS=200; % 175
    end
    % sampling at the center of NS "bins", each of width 1/NS, covering 0 to 1
    GRx  = repmat( linspace(grange(1)+1/(2*NS),grange(2)-1/(2*NS),NS), NS,1 );% (repmat( (1:NS),  NS,1 )-0.5)./NS;
    GRy  = repmat( linspace(grange(1)+1/(2*NS),grange(2)-1/(2*NS),NS)',1,NS );% (repmat( (1:NS)', 1,NS )-0.5)./NS;
    del2 = (GRx(1,2)-GRx(1,1))^2; 
    %
    mu_0 = mean(GRx(1,:)); % (2- -1)/2
    vr_0 = var( GRx(1,:));
    out.grid_info = [grange, NS, del2];
else
    error('unrecognized model!')
end

% ===> optimization loop (mode=1 or 3)
if gl_mode == 2 || gl_mode==3
    
    for id = 1:size(d_set,1)
    
        d = d_set(id,:);
        % shifted array(s) of matched size
        clear I12;
        % ---
        for k=1:3
            if d(k)>=0
                run1xyz(k,:) = [1,  (nn(k)-d(k))];
                run2xyz(k,:) = [(1+d(k)),  nn(k)];
            else
                run2xyz(k,:) = [1,  (nn(k)+d(k))];
                run1xyz(k,:) = [(1-d(k)),  nn(k)];
            end
        end
        I12(:,1) = reshape( I( run1xyz(1,1):run1xyz(1,2), run1xyz(2,1):run1xyz(2,2), run1xyz(3,1):run1xyz(3,2) )  ,[],1 );
        I12(:,2) = reshape( I( run2xyz(1,1):run2xyz(1,2), run2xyz(2,1):run2xyz(2,2), run2xyz(3,1):run2xyz(3,2) )  ,[],1 );
        % drop non-finite entries
        I12( ~isfinite(prod(I12,2)), : ) = [];
        xy = I12;
        N  = size(xy,1);
    
        % ====>> ====> display mode
        if makefigs==1 && id==1
            figure;
        end
        % ====>> ====> display mode

        % proceeding with kernelizations...
    
        if     strcmpi(model,'Hist' ) % "classic" histogram model
            
        elseif strcmpi(model,'HDE'  ) % histo density estimaqtion
            
            vbig = 3.5 * sqrt(mean(var(xy))) * size(xy,1)^(-1/4); % scot normal ref
            vlist = exp( linspace( log(vbig/5), log(2*vbig), 20 ) );
            disp('optimizing bandwidth...');
            for iv=1:numel(vlist)
                [id iv],
                vtmp=vlist(iv);      % temporary binwidth
                del2 = vtmp^2;
                NB  =ceil( 1/vtmp ); % num-bins, stretchf by +1 bin if spill-over
                %
                P0 = zeros(NB,NB);
                for bx=1:NB
                    Ix = (xy(:,1) >= vtmp*(bx-1)) & (xy(:,1) < vtmp*(bx));
                    for by=1:NB
                        Iy = (xy(:,2) >= vtmp*(by-1)) & (xy(:,2) < vtmp*(by));
                        P0(by,bx) = sum( Ix & Iy )/N;
                    end
                end
                KD = P0./del2; % --only used for display mode here!

                % ====>> ====> display mode
                if makefigs==1 && id==1
                    if sum(iv==[3 13 17])>0
                        if iv==5
                            subplot(2,3,1);
                            plot( xy(:,1), xy(:,2),'.k');
                            xlim([0 1]); ylim([0 1]);
                        end
                        subplot(2,3,round(iv/5)+3);
                        imagesc( KD );
                        set(gca,'YDir','normal');
                        runst= ((1:NB)*vtmp - vtmp/2);
                        xlim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
                        ylim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
                        linx = round(linspace( 2, numel(runst)-1, 4 ));
                        finx = round(runst(linx)*100)./100;
                        set(gca,'Xtick',linx);
                        set(gca,'XtickLabel',finx);
                        set(gca,'Ytick',linx);
                        set(gca,'YtickLabel',finx);
                    end
                end
                % ====>> ====> display mode

                CV(iv,id) = 2/(N-1)/del2 - ((N+1)/((N-1) * del2)) .* sum( P0(:).^2 );
            end
    
        elseif strcmpi(model,'KDEU') % kernel density estimation - gaussian
            
            vbig  = 1.06 * sqrt(mean(var(xy))) * size(xy,1)^(-1/6); % stdev - normal reference rule
            vlist = exp( linspace( log(vbig/5), log(2*vbig), 20 ) );
            pdist = (xy(:,1) - xy(:,1)').^2 + (xy(:,2) - xy(:,2)').^2;
    
            disp('optimizing bandwidth...');
            for iv=1:numel(vlist)
                [id iv],
                vtmp=vlist(iv);
                KD=0;
                for s=1:N
                    KD = KD + (1/(2*pi*(vtmp^2)*N)) * exp( -(1/(2*(vtmp^2))).*( (xy(s,1)-GRx).^2 + (xy(s,2)-GRy).^2 ) );
                end
                if abs((sum(KD(:))*del2) - 1) > 0.0005 % .05% deviation
                    abs((sum(KD(:))*del2) - 1),
                    error('KD bounding failed -- use a finer grid')
                end
                % ====>> ====> display mode
                if makefigs==1 && id==1 
                    if sum(iv==[5 10 15])>0
                        if iv==5
                            subplot(2,3,1);
                            plot( xy(:,1), xy(:,2),'.k');
                            xlim([0 1]); ylim([0 1]);
                        end
                        subplot(2,3,round(iv/5)+3);
                        imagesc( KD );
                        set(gca,'YDir','normal');
%                         runst= ((1:NB)*vtmp - vtmp/2);
%                         xlim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
%                         ylim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
%                         linx = round(linspace( 2, numel(runst)-1, 4 ));
%                         finx = round(runst(linx)*100)./100;
%                         set(gca,'Xtick',linx);
%                         set(gca,'XtickLabel',finx);
%                         set(gca,'Ytick',linx);
%                         set(gca,'YtickLabel',finx);
                    end
                end
                % ====>> ====> display mode

                S1 = sum( KD(:).^2 ) * del2;
    
                KD = (1/(2*pi*(vtmp^2)*(N-1))) * exp( -(1/(2*(vtmp^2))).*pdist );
                KD = KD .* (1-eye(size(KD)));
                S2 = sum(sum(KD,2)) * (2/N);
    
                CV(iv,id) = S1 - S2;
            end

        elseif strcmpi(model,'KDE') % kernel density estimation - gaussian

            vbig  = 1.06 * sqrt(mean(var(xy))) * size(xy,1)^(-1/6); % stdev - normal reference rule
            vlist = exp( linspace( log(vbig/5), log(2*vbig), 20 ) );
            pdist = (xy(:,1) - xy(:,1)').^2 + (xy(:,2) - xy(:,2)').^2;
            subaset = ones(N,1);
            TFACsub = ones(N,1);

            disp('optimizing bandwidth...');
            for iv=1:numel(vlist)
                [id iv],
                vtmp=vlist(iv);
                KD=0;
                for s=1:N
                    subaset(s,1) = (normcdf(1,xy(s,1),vtmp)-normcdf(0,xy(s,1),vtmp))*(normcdf(1,xy(s,2),vtmp)-normcdf(0,xy(s,2),vtmp));
                    KD = KD + (  (1/(2*pi*(vtmp^2)*N)) * exp( -(1/(2*(vtmp^2))).*( (xy(s,1)-GRx).^2 + (xy(s,2)-GRy).^2 ) )  );
                end
                TFAC = sum(subaset)/N; % truncation factor
                KD   = KD./TFAC;
                for s=1:N
                    TFACsub(s,1) = (sum(subaset) - subaset(s))/(N-1);
                end

                if abs((sum(KD(:))*del2) - 1) > 0.0005 % .05% deviation
                    abs((sum(KD(:))*del2) - 1),
                    error('KD bounding failed -- use a finer grid')
                end
                % ====>> ====> display mode
                if makefigs==1 && id==1 
                    if sum(iv==[5 10 15])>0
                        if iv==5
                            subplot(2,3,1);
                            plot( xy(:,1), xy(:,2),'.k');
                            xlim([0 1]); ylim([0 1]);
                        end
                        subplot(2,3,round(iv/5)+3);
                        imagesc( KD );
                        set(gca,'YDir','normal');
%                         runst= ((1:NB)*vtmp - vtmp/2);
%                         xlim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
%                         ylim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
%                         linx = round(linspace( 2, numel(runst)-1, 4 ));
%                         finx = round(runst(linx)*100)./100;
%                         set(gca,'Xtick',linx);
%                         set(gca,'XtickLabel',finx);
%                         set(gca,'Ytick',linx);
%                         set(gca,'YtickLabel',finx);
                    end
                end
                % ====>> ====> display mode

                S1 = sum( KD(:).^2 ) * del2;
    
                KD = (1/(2*pi*(vtmp^2)*(N-1))) * exp( -(1/(2*(vtmp^2))).*pdist );
                KD = KD .* (1-eye(size(KD)));
                S2 = sum(sum(KD,2)./TFACsub) * (2/N);
    
                CV(iv,id) = S1 - S2;

                %%%
                out.TFset(iv,id) =TFAC;
                %%%
            end
        end

        % store scale-list
        SZ(:,id) = vlist;
        % optima
        [vx,ixxv]=min(CV(:,id));
        out.Copt_err(id,1)  = vx;
        out.Copt_idx(id,1)  = ixxv;
        out.Copt_bvl(id,1)  = vlist(ixxv);
        out.Copt_bref(id,1) = vbig;
        out.CVfull = CV;
        out.SZfull = SZ;

        % ====>> ====> display mode
        if makefigs==1 && id==1 
            sfh1=subplot(2,3,2);
            cv_smo = (CV(:,id) + [0; CV(1:end-1,id)] + [CV(2:end,id); 0])./3;
            plot( vlist,cv_smo,'o-k','markerfacecolor',[0.5 0.5 0.5] );
            [vxt,ixt]=min(cv_smo);
            hold on; plot( vlist(ixt), cv_smo(ixt,1),'xr','markersize',10 );
            sfh1.Position = sfh1.Position + [0 0 0.28 0];
            xlabel('kernel scale');
            ylabel('cross-validation error');
        end
        % ====>> ====> display mode
    end
    
    % for now, collapse across directions
    out.Bvl_for_fitt = mean(out.Copt_bvl);
else
    % -- predefined bandwidth set
    out.Bvl_for_fitt = b_width;
end

%% ===> estimation loop (mode=2 or 3)

if gl_mode == 1 || gl_mode==3
    
    for id = 1:size(d_set,1)
    
        d = d_set(id,:);
        % shifted array(s) of matched size
        clear I12;
    
        % ---
        for k=1:3
            if d(k)>=0
                run1xyz(k,:) = [1,  (nn(k)-d(k))];
                run2xyz(k,:) = [(1+d(k)),  nn(k)];
            else
                run2xyz(k,:) = [1,  (nn(k)+d(k))];
                run1xyz(k,:) = [(1-d(k)),  nn(k)];
            end
        end
        I12(:,1) = reshape( I( run1xyz(1,1):run1xyz(1,2), run1xyz(2,1):run1xyz(2,2), run1xyz(3,1):run1xyz(3,2) )  ,[],1 );
        I12(:,2) = reshape( I( run2xyz(1,1):run2xyz(1,2), run2xyz(2,1):run2xyz(2,2), run2xyz(3,1):run2xyz(3,2) )  ,[],1 );
        % drop non-finite entries
        I12( ~isfinite(prod(I12,2)), : ) = [];
        xy = I12;
        N  = size(xy,1); % get sample size once non-finite dropped

        if scaled_noise>0 && scaled_shift>0
            error('cannot have both tested')
        end
        if scaled_noise>0
            sdsig = (std(xy,0,1));
            noise = randn(size(xy));
            sdnoi = (std(noise,0,1));
            mnnoi = (mean(noise,1));
            xy_n_noise = xy + ((noise-mnnoi).*(scaled_noise.*sdsig./sdnoi));
            xy_n_noise(xy_n_noise>(1-eps))=(1-eps);
            xy_n_noise(xy_n_noise<eps)=eps;
            xy = xy_n_noise; clear xy_n_noise; % replace with noisy version!
        end
        if scaled_shift>0
            xy_n_noise = xy + [abs(scaled_shift) abs(scaled_shift)];
            xy_n_noise(xy_n_noise>(1-eps))=(1-eps);
            xy_n_noise(xy_n_noise<eps)=eps;
            xy = xy_n_noise; clear xy_n_noise; % replace with shifted version!
        end

        % ===> proceeing with kernelizations
    
        if     strcmpi(model,'Hist' ) % "classic" histogram model
            
        elseif strcmpi(model,'HDE'  ) % histo density estimaqtion
    
            vbig = 3.5 * sqrt(mean(var(xy))) * size(xy,1)^(-1/4); % scot normal ref
            vtmp = out.Bvl_for_fitt;
            del2 = vtmp^2;
            NB  =ceil( 1/vtmp ); % stretch by +1 bin if spill-over
            %
            P0 = zeros(NB,NB);
            for bx=1:NB
                Ix = (xy(:,1) >= vtmp*(bx-1)) & (xy(:,1) < vtmp*(bx));
                for by=1:NB
                    Iy = (xy(:,2) >= vtmp*(by-1)) & (xy(:,2) < vtmp*(by));
                    P0(by,bx) = sum( Ix & Iy )/N;
                end
            end
            KD = P0./del2;

            %> increments in steps of ~0.001
            %> 
            
            % here, they grids are bin midpointses
            GRx  = repmat( ((1:NB)*vtmp - vtmp/2) ,NB,1 );% (repmat( (1:NS),  NS,1 )-0.5)./NS;
            GRy  = repmat( ((1:NB)*vtmp - vtmp/2)',1,NB );% (repmat( (1:NS)', 1,NS )-0.5)./NS;
            %
            mu_0 = mean(GRx(1,:)); % (2- -1)/2
            vr_0 = var( GRx(1,:));
    
        elseif strcmpi(model,'KDEU') % kernel density estimation - gaussian
            
            vbig = 1.06 * sqrt(mean(var(xy))) * size(xy,1)^(-1/6); % stdev - normal reference rule

            vtmp = out.Bvl_for_fitt;
            % fully gridding out samples
            KD=0;
            for s=1:N
                KD = KD + (1/(2*pi*(vtmp^2)*N)) * exp( -(1/(2*(vtmp^2))).*( (xy(s,1)-GRx).^2 + (xy(s,2)-GRy).^2 ) );
            end
            if abs((sum(KD(:))*del2) - 1) > 0.0005 % .05% deviation
                error('KD bounding failed -- use a finer grid')
            end

        elseif strcmpi(model,'KDE') % kernel density estimation - gaussian
            
            vbig = 1.06 * sqrt(mean(var(xy))) * size(xy,1)^(-1/6); % stdev - normal reference rule
            subaset = ones(N,1);

            vtmp = out.Bvl_for_fitt;
            % fully gridding out samples
            KD=0;
            for s=1:N
                subaset(s,1) = (normcdf(1,xy(s,1),vtmp)-normcdf(0,xy(s,1),vtmp))*(normcdf(1,xy(s,2),vtmp)-normcdf(0,xy(s,2),vtmp));
                KD = KD + (1/(2*pi*(vtmp^2)*N)) * exp( -(1/(2*(vtmp^2))).*( (xy(s,1)-GRx).^2 + (xy(s,2)-GRy).^2 ) );
            end
            TFAC = sum(subaset)/N; % truncation factor
            KD   = KD./TFAC;
            if abs((sum(KD(:))*del2) - 1) > 0.0005 % .05% deviation
                error('KD bounding failed -- use a finer grid')
            end
            out.TFopt = TFAC;
        end
        if makefigs==1 && id==1
            figure;
            subplot(1,2,1);
            plot( xy(:,1), xy(:,2),'.k');
            xlim([0 1]); ylim([0 1]);
            subplot(1,2,2);
            imagesc( KD );
            set(gca,'YDir','normal');
%             runst= ((1:NB)*vtmp - vtmp/2);
%             xlim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
%             ylim([1-(runst(1)./vtmp), NB+(1-runst(end))./vtmp])
%             linx = round(linspace( 2, numel(runst)-1, 4 ));
%             finx = round(runst(linx)*100)./100;
%             set(gca,'Xtick',linx)
%             set(gca,'XtickLabel',finx)
%             set(gca,'Ytick',linx)
%             set(gca,'YtickLabel',finx)
        end
    
        %--- start extracting features! ---%
        %
        % energy
        tmp_metrics(id,1) = sum( KD(:).^2 ) .* del2;
        % entropy
        lkd = log( KD(:) ); lkd(~isfinite(lkd))=0;
        tmp_metrics(id,2) = -sum( KD(:) .* lkd ) .* del2;
        % contrast
        tmp_metrics(id,3) = sum( KD(:) .* (GRx(:)-GRy(:)).^2 ) .* del2;
        % homogeneity
        tmp_metrics(id,4) = sum( KD(:) ./ (1 + (GRx(:)-GRy(:)).^2 ) ) .* del2;
        % correlation
        tmp_metrics(id,5) = sum( KD(:) .* ((GRx(:) - mu_0) .* (GRy(:) - mu_0))./vr_0 ) .* del2;
        % clust shade
        tmp_metrics(id,6) = abs(sum( KD(:) .* (GRx(:) + GRy(:) - 2*mu_0 ).^3 )) .* del2;
        % clust promin
        tmp_metrics(id,7) = sum( KD(:) .* (GRx(:) + GRy(:) - 2*mu_0 ).^4 ) .* del2;
        
        %-------
        bscalset(id,1) = vtmp/vbig;
    end
    
    % consolidate
    out.metrics_av = mean(tmp_metrics,  1);
    out.metrics_sd =  std(tmp_metrics,0,1);
    out.bscalset   = bscalset;
end

