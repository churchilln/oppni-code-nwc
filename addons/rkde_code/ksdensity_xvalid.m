function out = ksdensity_xvalid( Xin, U0, type )

if(nargin<3) type='LOO'; end

Nf = round(numel(Xin)/10);
Nh = round(numel(Xin)/2);

Urng = exp(linspace(log(U0/10), log(U0*10), 30 ));
%Urng = linspace(U0/10, U0*10, 100 );

% define range for numerical integration
 MIN=min(Xin(:))-range(Xin(:))/10;
 MAX=max(Xin(:))+range(Xin(:))/10;
erang=linspace(MIN,MAX,2000);

for(j=1:numel(Urng)) % for each kernel of interest
    
    [P2,~] = ksdensity(Xin(:), erang, 'Bandwidth', Urng(j) );% 'Function','cdf' );
    f2a(j,1) = sum( P2.^2 .* (erang(2)-erang(1)) );
    
    if(j==1)
        erang2=linspace(MIN,MAX,numel(erang)*2);
        [P2,~] = ksdensity(Xin(:), erang, 'Bandwidth', Urng(j) );% 'Function','cdf' );
        f2a_chk = sum( P2.^2 .* (erang(2)-erang(1)) );
        if( abs(f2a(j,1) - f2a_chk)/f2a(j,1) > 1E-6 ) 
            error('choose a finer sampling range for numeric integration'); 
        else
            %disp('numeric grid ok');
        end
    end
    
    if(strcmpi(type,'LOO'))

        for(i=1:numel(Xin))
            %[i,j],
            Xtrn =Xin; Xtrn(i)=[];
            % get
            [Px(i,j),~] = ksdensity(Xtrn(:), Xin(i), 'Bandwidth', Urng(j) );% 'Function','cdf' );
        end
        
    elseif(strcmpi(type,'10FOLD'))

        for(i=1:10)
            ixhold = (1:Nf) + (i-1)*Nf;
            Xtrn = Xin; Xtrn(ixhold)=[];
            [Pxtmp,~] = ksdensity(Xtrn(:), Xin(ixhold), 'Bandwidth', Urng(j) );% 'Function','cdf' );
            Px(i,j) = mean(Pxtmp);
        end
    elseif(strcmpi(type,'SPLITHALF'))
        
        for(i=1:10)
            ixhold = randperm(numel(Xin));
            ixhold1= ixhold(1:Nh);
            ixhold2= ixhold(Nh+1:end);
            %
            Xtrn = Xin; Xtrn(ixhold1)=[];
            [Pxtmp,~] = ksdensity(Xtrn(:), Xin(ixhold1), 'Bandwidth', Urng(j) );% 'Function','cdf' );
            Px(i,j) = mean(Pxtmp);
            %
            Xtrn = Xin; Xtrn(ixhold2)=[];
            [Pxtmp,~] = ksdensity(Xtrn(:), Xin(ixhold2), 'Bandwidth', Urng(j) );% 'Function','cdf' );
            Px(i+10,j) = mean(Pxtmp);
        end
    end
end

% figure, imagesc( Px );
% figure, plot( Px' );

[vx(1),ix(1)]=min( f2a - 2*mean( Px,1)' );
[vx(2),ix(2)]=min( f2a - 2*median( Px,1)' );
[ix(1) ix(2)],
% figure, plot( [f2a - 2*mean( Px,1)', f2a - 2*median( Px,1)'] );

if(ix(2)==1 || ix(2)==numel(Urng) )
    error('BW optimizes at the end of the explored range -- expand range limits!');
end
out.BW = Urng(ix);
