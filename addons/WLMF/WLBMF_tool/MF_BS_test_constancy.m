function [test, estm, logstatm] = MF_BS_test_constancy(data, M, paramEst, paramBS, verbose, FigNum);
% function [test, estm, logstatm] = MF_BS_test_constancy(data, M, paramEst, paramBS, verbose, FigNum);
%
%   Wavelet domain double block bootstrap based test of hypothesis of time constancy of multifractal attributes
%
% SEE demo_MF_BS_test_constancy for usage, outputs and parameters
%
% Herwig Wendt, Lyon, 2006 - 2008
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this tool, please quote:
%
%     @ARTICLE{WendtSPM2007,
%       author = {Herwig Wendt and Patrice Abry and Stéphane Jaffard},
%       title = {Bootstrap for Empirical Multifractal Analysis},
%       journal = {IEEE Signal Proc. Mag.},
%       volume = {24},
%       number = {4},
%       pages = {38--48},
%       year = {2007},
%     }
%
% and:
%
%     @INPROCEEDINGS{WendtICASSP2008,
%       author = {Herwig Wendt and Patrice Abry},
%       title = {Bootstrap tests for the time constancy of multifractal attributes},
%       booktitle = {Proc. IEEE Int. Conf. Acoust., Speech, and Signal Proc. (ICASSP)},
%       address = {Las Vegas, USA},
%       year = {2008},
%     }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


try 
    verbose;
catch
    verbose=0; end;
try 
    FigNum;
catch
    FigNum=0; end;

if FigNum==0;  PLOT=0; else;  PLOT=1; end;

MomNul=paramEst.MomNul;
j1=paramEst.j1;
j2=paramEst.j2;
methodWT=paramEst.methodWT;
EstFun=paramEst.EstFun;
B1=paramBS.Nresamp1;
T_S=paramBS.T_S;
Alpha=paramBS.Alpha;
bid=min(B1, ceil((B1+1)*(1-Alpha)));

Fun=bin2dec(num2str(EstFun));
poscount=1;
if Fun>=4; ZQ=1; zqpos=1; poscount=poscount+1; else; ZQ=0; end
if (Fun~=1) && (Fun~=4) && (Fun~=5); DH=1; Dqpos=poscount;  poscount=poscount+1;hqpos=poscount;  poscount=poscount+1; else;  DH=0; end
if rem(Fun,2)==1; CP=1; Cppos=poscount; else; CP=0;  end

if length(find(methodWT==1)); DWT=1;  else DWT=0;  end
if length(find(methodWT==2)); LWT=1;  else LWT=0; end


% --- Calculate Coef Lead
[coef, lead, nj] = DxLx1d(data, MomNul);
if T_S; j2max=find(nj.L>=8*M); else; j2max=find(nj.L>=4*M); end
j2max=j2max(end); j2=min(j2, j2max);

% --- Calculate Estimates and BS estimates
if DWT
    [est_m, logstat_m, est_all] = MF_const_BS(coef, nj.W, M, paramBS, paramEst);
    LEST=length(est_m{1}.t);
    % Convenience Variable
    for m=1:M
        tD(m,:)=est_m{m}.t;        % estimates(m)
        tBD(m,:,:)=est_m{m}.T';    % Bootstrap estimates(m)
        % STD estimates
        std_tBD(m,:)=est_m{m}.stdT';
        std_tBBD(m,:,:)=est_m{m}.stdTT';
        tVtheo(m)=est_m{m}.Vtheo;
        nj_m(m,:)=est_m{m}.nj;
        logstatm.DWT.Sqj(m,:,:)=logstat_m{m}.Sqj;
        logstatm.DWT.SqjBcut(m,:,:,:)=logstat_m{m}.SqjBcut;
        logstatm.DWT.SqjBall(m,:,:,:)=logstat_m{m}.SqjBall;
    end
    logstatm.DWT.Sqj=squeeze(logstatm.DWT.Sqj);
    logstatm.DWT.SqjBcut=squeeze(logstatm.DWT.SqjBcut);
    logstatm.DWT.SqjBall=squeeze(logstatm.DWT.SqjBall);
    % --- Abry Veitch Test, Graybill and Deal estimator for common mean, bootstrap
    %     using BS var est and GDE for common mean
    mtD_GD = sum(tD./std_tBD.^2)./sum(1./std_tBD.^2);
    Theta_GD = sum( 1./std_tBD.^2 .* ( tD  -  repmat(mtD_GD, M, 1)).^2);
    for bsid=1:B1
        mstartmp=squeeze(tBD(:,bsid,:));
        sstartmp=squeeze(std_tBBD(:,bsid,:));
        mmstar_GD=sum(mstartmp./sstartmp.^2)./sum(1./sstartmp.^2);
        ThetaB_GD(bsid,:)=sum(1./sstartmp.^2 .* (mstartmp-repmat(mmstar_GD, M, 1)).^2);
    end
    ThetaB_GD=sort(squeeze(ThetaB_GD));
    ThetaBcrit_GD=squeeze(ThetaB_GD(bid,:)); % test critical value
    for Estid=1:LEST
        if Theta_GD(Estid)>ThetaBcrit_GD(Estid); hh_GD( Estid)=1; else; hh_GD( Estid)=0; end % test decision
        ltmp=length(find(squeeze(ThetaB_GD(:,Estid))>=Theta_GD(Estid))); pp_GD( Estid)=ltmp/(B1+1); % p value
    end

    test.DWT.t=Theta_GD;
    test.DWT.T=ThetaB_GD;
    test.DWT.d=hh_GD;
    test.DWT.p=pp_GD;

    estm.DWT.LEst=est_m{m}.LEst;
    estm.DWT.t=tD;
    estm.DWT.T=tBD;
    estm.DWT.stdt=std_tBD;
    estm.DWT.stdT=std_tBBD;
    estm.DWT.Vtheo=tVtheo;
    estm.DWT.nj=nj_m;
    estm.DWT.j1=est_m{m}.j1;
    estm.DWT.j2=est_m{m}.j2;
    estm.DWT.t_all=est_all.t;
    estm.DWT.j1all=est_all.j1;
    estm.DWT.j2all=est_all.j2;
    estm.DWT.njall=est_all.nj;
   %logstatm;
end

if LWT
    [est_m, logstat_m, est_all] = MF_const_BS(lead, nj.L, M, paramBS, paramEst);
    LEST=length(est_m{1}.t);
    % Convenience Variable
    for m=1:M
        tD(m,:)=est_m{m}.t;        % estimates(m)
        tBD(m,:,:)=est_m{m}.T';    % Bootstrap estimates(m)
        % STD estimates
        std_tBD(m,:)=est_m{m}.stdT';
        std_tBBD(m,:,:)=est_m{m}.stdTT';
        nj_m(m,:)=est_m{m}.nj;
        logstatm.LWT.Sqj(m,:,:)=logstat_m{m}.Sqj;
        logstatm.LWT.SqjBcut(m,:,:,:)=logstat_m{m}.SqjBcut;
        logstatm.LWT.SqjBall(m,:,:,:)=logstat_m{m}.SqjBall;
    end
    logstatm.LWT.Sqj=squeeze(logstatm.LWT.Sqj);
    logstatm.LWT.SqjBcut=squeeze(logstatm.LWT.SqjBcut);
    logstatm.LWT.SqjBall=squeeze(logstatm.LWT.SqjBall);
    % --- Abry Veitch Test, Graybill and Deal estimator for common mean, bootstrap
    %     using BS var est and GDE for common mean
    mtD_GD = sum(tD./std_tBD.^2)./sum(1./std_tBD.^2);
    Theta_GD = sum( 1./std_tBD.^2 .* ( tD  -  repmat(mtD_GD, M, 1)).^2);
    for bsid=1:B1
        mstartmp=squeeze(tBD(:,bsid,:));
        sstartmp=squeeze(std_tBBD(:,bsid,:));
        mmstar_GD=sum(mstartmp./sstartmp.^2)./sum(1./sstartmp.^2);
        ThetaB_GD(bsid,:)=sum(1./sstartmp.^2 .* (mstartmp-repmat(mmstar_GD, M, 1)).^2);
    end
    ThetaB_GD=sort(squeeze(ThetaB_GD));
    ThetaBcrit_GD=squeeze(ThetaB_GD(bid,:)); % test critical value
    for Estid=1:LEST
        if Theta_GD(Estid)>ThetaBcrit_GD(Estid); hh_GD( Estid)=1; else; hh_GD( Estid)=0; end % test decision
        ltmp=length(find(squeeze(ThetaB_GD(:,Estid))>=Theta_GD(Estid))); pp_GD( Estid)=ltmp/(B1+1); % p value
    end

    test.LWT.t=Theta_GD;
    test.LWT.T=ThetaB_GD;
    test.LWT.d=hh_GD;
    test.LWT.p=pp_GD;

    estm.LWT.LEst=est_m{m}.LEst;
    estm.LWT.t=tD;
    estm.LWT.T=tBD;
    estm.LWT.stdt=std_tBD;
    estm.LWT.stdT=std_tBBD;
    estm.LWT.nj=nj_m;
    estm.LWT.j1=est_m{m}.j1;
    estm.LWT.j2=est_m{m}.j2;
    estm.LWT.t_all=est_all.t;
    estm.LWT.j1all=est_all.j1;
    estm.LWT.j2all=est_all.j2;
    estm.LWT.njall=est_all.nj;
end

if verbose>0
    q=paramEst.q;
    % ESTIMATION RESULTS
    if verbose>2
        fprintf(' ******************************************* \n');
        fprintf(' *************** ESTIMATION **************** \n');
        fprintf(' ******************************************* \n');
        if DWT
            LW=estm.DWT.LEst;
            fprintf(' ------------------------------------------- \n');
            fprintf(' WAVELET COEFFICIENTS ESTIMATION : \n');
            fprintf(' ------------------------------------------- \n');
            fprintf('       ');
            for id=1:LEST
                if ZQ&(id<=LW(2))
                    fprintf(' z(q=%1.1g) ',q(id));
                end
                if DH
                    if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                        fprintf(' D(q=%1.1g) ',q(id-LW(1+ZQ)));
                    end
                    if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                        fprintf(' h(q=%1.1g) ',q(id-LW(2+ZQ)));
                    end
                end
                if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                    fprintf(' c%1.1g     ',(id-LW(1+ZQ+2*DH)));
                end
            end
            fprintf('\n');
            fprintf(' est:  ');
            for id=1:LEST
                if ZQ&(id<=LW(2))
                    fprintf(' %1.3f  ',estm.DWT.t_all(id));
                end
                if DH
                    if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                        fprintf(' %1.3f  ',estm.DWT.t_all(id));
                    end
                    if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                        fprintf(' %1.3f  ',estm.DWT.t_all(id));
                    end
                end
                if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                    fprintf(' %1.3f  ',estm.DWT.t_all(id));
                end
            end            
            fprintf('\n');
            for m=1:M
                fprintf(' m=%2g: ',m);
                for id=1:LEST
                    if ZQ&(id<=LW(2))
                        fprintf(' %1.3f  ',estm.DWT.t(m,id));
                    end
                    if DH
                        if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                            fprintf(' %1.3f  ',estm.DWT.t(m,id));
                        end
                        if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                            fprintf(' %1.3f  ',estm.DWT.t(m,id));
                        end
                    end
                    if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                        fprintf(' %1.3f  ',estm.DWT.t(m,id));
                    end
                end
                fprintf('\n');
            end
        end
        if LWT
            LW=estm.LWT.LEst;
            fprintf(' ------------------------------------------- \n');
            fprintf(' WAVELET LEADERS ESTIMATION : \n');
            fprintf(' ------------------------------------------- \n');
            fprintf('       ');
            for id=1:LEST
                if ZQ&(id<=LW(2))
                    fprintf(' z(q=%1.1g) ',q(id));
                end
                if DH
                    if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                        fprintf(' D(q=%1.1g) ',q(id-LW(1+ZQ)));
                    end
                    if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                        fprintf(' h(q=%1.1g) ',q(id-LW(2+ZQ)));
                    end
                end
                if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                    fprintf(' c%1.1g     ',(id-LW(1+ZQ+2*DH)));
                end
            end
            fprintf('\n');
            fprintf(' est:  ');
            for id=1:LEST
                if ZQ&(id<=LW(2))
                    fprintf(' %1.3f  ',estm.LWT.t_all(id));
                end
                if DH
                    if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                        fprintf(' %1.3f  ',estm.LWT.t_all(id));
                    end
                    if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                        fprintf(' %1.3f  ',estm.LWT.t_all(id));
                    end
                end
                if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                    fprintf(' %1.3f  ',estm.LWT.t_all(id));
                end
            end            
            fprintf('\n');
            for m=1:M
                fprintf(' m=%2g: ',m);
                for id=1:LEST
                    if ZQ&(id<=LW(2))
                        fprintf(' %1.3f  ',estm.LWT.t(m,id));
                    end
                    if DH
                        if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                            fprintf(' %1.3f  ',estm.LWT.t(m,id));
                        end
                        if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                            fprintf(' %1.3f  ',estm.LWT.t(m,id));
                        end
                    end
                    if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                        fprintf(' %1.3f  ',estm.LWT.t(m,id));
                    end
                end
                fprintf('\n');
            end
            
        end        
        fprintf('\n ******************************************* \n');
        fprintf(' ********** TIME CONSTANCY TESTS *********** \n');
        fprintf(' ******************************************* \n\n');
    end
    % TEST RESULTS
    if DWT
        LW=estm.DWT.LEst;
        fprintf(' ------------------------------------------- \n');
        fprintf(' WAVELET COEFFICIENTS TIME CONSTANCY TESTS : \n');
        fprintf(' ------------------------------------------- \n');
        fprintf(' Estimate  decision  pvalue ');
        if verbose>1
            fprintf(' | statistic  crit.value \n');
        else
            fprintf('\n');
        end
        for id=1:LEST
            if ZQ&(id<=LW(2))
                fprintf('  z(q=%1.1g):     %1g      %1.3f ',q(id), test.DWT.d(id), test.DWT.p(id));
                if verbose>1
                    fprintf('  |  %1.3f      %1.3f\n',test.DWT.t(id), test.DWT.T(bid,id));
                else
                    fprintf('\n');
                end
                if id==LW(2); disp('  -------'); end
            end
            if DH
                if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                    fprintf('  D(q=%1.1g):     %1g      %1.3f ',q(id-LW(1+ZQ)), test.DWT.d(id), test.DWT.p(id));
                    if verbose>1
                        fprintf('  |  %1.3f      %1.3f\n',test.DWT.t(id), test.DWT.T(bid,id));
                    else
                        fprintf('\n');
                    end
                    if id==LW(2+ZQ); disp('  -------'); end
                end
                if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                    fprintf('  h(q=%1.1g):     %1g      %1.3f ',q(id-LW(2+ZQ)), test.DWT.d(id), test.DWT.p(id));
                    if verbose>1                        
                        fprintf('  |  %1.3f      %1.3f\n',test.DWT.t(id), test.DWT.T(bid,id));
                    else
                        fprintf('\n');
                    end
                    if id==LW(3+ZQ); disp('  -------'); end
                end
            end
            if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                fprintf('    c%1.1g:       %1g      %1.3f ',(id-LW(1+ZQ+2*DH)), test.DWT.d(id), test.DWT.p(id));
                if verbose>1
                    fprintf('  |  %1.3f      %1.3f\n',test.DWT.t(id), test.DWT.T(bid,id));
                else
                    fprintf('\n');
                end
                if id==LW(2+ZQ+2*DH); disp('  -------'); end
            end
        end
    end
    if LWT
        LW=estm.LWT.LEst;
        fprintf(' ------------------------------------------- \n');
        fprintf(' WAVELET LEADERS TIME CONSTANCY TESTS : \n');
        fprintf(' ------------------------------------------- \n');
        fprintf(' Estimate  decision  pvalue ');
        if verbose>1
            fprintf(' | statistic  crit.value \n');
        else
            fprintf('\n');
        end
        for id=1:LEST
            if ZQ&(id<=LW(2))
                fprintf('  z(q=%1.1g):     %1g      %1.3f ',q(id), test.LWT.d(id), test.LWT.p(id));
                if verbose>1
                    fprintf('  |  %1.3f      %1.3f\n',test.LWT.t(id), test.LWT.T(bid,id));
                else
                    fprintf('\n');
                end
                if id==LW(2); disp('  -------'); end
            end
            if DH
                if (id>LW(1+ZQ))&(id<=LW(2+ZQ))
                    fprintf('  D(q=%1.1g)      %1g      %1.3f ',q(id-LW(1+ZQ)), test.LWT.d(id), test.LWT.p(id));
                    if verbose>1
                        fprintf('  |  %1.3f      %1.3f\n',test.LWT.t(id), test.LWT.T(bid,id));
                    else
                        fprintf('\n');
                    end
                    if id==LW(2+ZQ); disp('  -------'); end
                end
                if (id>LW(2+ZQ))&(id<=LW(3+ZQ))
                    fprintf('  h(q=%1.1g):     %1g      %1.3f ',(id-LW(2+ZQ)), test.LWT.d(id), test.LWT.p(id));
                    if verbose>1
                        fprintf('  |  %1.3f      %1.3f\n',test.LWT.t(id), test.LWT.T(bid,id));
                    else
                        fprintf('\n');
                    end
                    if id==LW(3+ZQ); disp('  -------'); end
                end
            end
            if CP&(id>LW(1+ZQ+2*DH))&(id<=LW(2+ZQ+2*DH))
                fprintf('    c%1.1g:       %1g      %1.3f ',(id-LW(1+ZQ+2*DH)), test.LWT.d(id), test.LWT.p(id));
                if verbose>1
                    fprintf('  |  %1.3f      %1.3f\n',test.LWT.t(id), test.LWT.T(bid,id));
                else
                    fprintf('\n');
                end
                if id==LW(2+ZQ+2*DH); disp('  -------'); end
            end
        end
    end

end

if PLOT
    FontSize=12;
    FontSize1=10;
    if DWT
        LW=estm.DWT.LEst;
        if ZQ
            figure(FigNum); clf;
            Ltmp=LW(2)-LW(1);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for id=1:Ltmp;
                sb=subplot(nrow,ncol,id);
                dtmp=test.DWT.d(id); ptmp=test.DWT.p(id);
                ttmp=estm.DWT.t(:,id); sttmp=estm.DWT.stdt(:,id);
                mTtmp=mean(estm.DWT.T(:,:,id),2); sTtmp=std(estm.DWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['\zeta(q=',num2str(q(id)),') : coefficients'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
        end
        if DH
            figure(FigNum+ZQ); clf;
            Ltmp=LW(2+ZQ)-LW(1+ZQ);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for i=1:Ltmp;
                id=i+LW(1+ZQ);
                sb=subplot(nrow,ncol,i);
                dtmp=test.DWT.d(id); ptmp=test.DWT.p(id);
                ttmp=estm.DWT.t(:,id); sttmp=estm.DWT.stdt(:,id);
                mTtmp=mean(estm.DWT.T(:,:,id),2); sTtmp=std(estm.DWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['D(q=',num2str(q(i)),') : coefficients'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
            figure(FigNum+ZQ+DH); clf;
            Ltmp=LW(2+ZQ+DH)-LW(1+ZQ+DH);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for i=1:Ltmp;
                id=i+LW(1+ZQ+DH);
                sb=subplot(nrow,ncol,i);
                dtmp=test.DWT.d(id); ptmp=test.DWT.p(id);
                ttmp=estm.DWT.t(:,id); sttmp=estm.DWT.stdt(:,id);
                mTtmp=mean(estm.DWT.T(:,:,id),2); sTtmp=std(estm.DWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['h(q=',num2str(q(i)),') : coefficients'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
        end
        if CP
            figure(FigNum+ZQ+2*DH); clf;
            Ltmp=LW(2+ZQ+2*DH)-LW(1+ZQ+2*DH);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for i=(1:Ltmp);
                id=i+LW(1+ZQ+2*DH);
                sb=subplot(nrow,ncol,i);
                dtmp=test.DWT.d(id); ptmp=test.DWT.p(id);
                ttmp=estm.DWT.t(:,id); sttmp=estm.DWT.stdt(:,id);
                mTtmp=mean(estm.DWT.T(:,:,id),2); sTtmp=std(estm.DWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['c_',num2str((i)),' : coefficients'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
        end
    end
    if LWT
        LW=estm.LWT.LEst;
        if DWT; Figshift=ZQ+CP+2*DH; else; Figshift=0; end;
        if ZQ
            figure(FigNum+Figshift); clf;
            Ltmp=LW(2)-LW(1);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for id=1:Ltmp;
                sb=subplot(nrow,ncol,id);
                dtmp=test.LWT.d(id); ptmp=test.LWT.p(id);
                ttmp=estm.LWT.t(:,id); sttmp=estm.LWT.stdt(:,id);
                mTtmp=mean(estm.LWT.T(:,:,id),2); sTtmp=std(estm.LWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['\zeta(q=',num2str(q(id)),') : Leaders'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
        end
        if DH
            figure(FigNum+ZQ+Figshift); clf;
            Ltmp=LW(2+ZQ)-LW(1+ZQ);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for i=1:Ltmp;
                id=i+LW(1+ZQ);
                sb=subplot(nrow,ncol,i);
                dtmp=test.LWT.d(id); ptmp=test.LWT.p(id);
                ttmp=estm.LWT.t(:,id); sttmp=estm.LWT.stdt(:,id);
                mTtmp=mean(estm.LWT.T(:,:,id),2); sTtmp=std(estm.LWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['D(q=',num2str(q(i)),') : Leaders'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
            figure(FigNum+ZQ+DH+Figshift); clf;
            Ltmp=LW(2+ZQ+DH)-LW(1+ZQ+DH);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for i=1:Ltmp;
                id=i+LW(1+ZQ+DH);
                sb=subplot(nrow,ncol,i);
                dtmp=test.LWT.d(id); ptmp=test.LWT.p(id);
                ttmp=estm.LWT.t(:,id); sttmp=estm.LWT.stdt(:,id);
                mTtmp=mean(estm.LWT.T(:,:,id),2); sTtmp=std(estm.LWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['h(q=',num2str(q(i)),') : Leaders'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
        end
        if CP
            figure(FigNum+ZQ+2*DH+Figshift); clf;
            Ltmp=LW(2+ZQ+2*DH)-LW(1+ZQ+2*DH);
            if Ltmp>1; ncol=2; else; ncol=1; end
            nrow=ceil(Ltmp/ncol);
            m1=1:M; m=m1-0.05; mref=m+0.1;
            for i=(1:Ltmp);
                id=i+LW(1+ZQ+2*DH);
                sb=subplot(nrow,ncol,i);
                dtmp=test.LWT.d(id); ptmp=test.LWT.p(id);
                ttmp=estm.LWT.t(:,id); sttmp=estm.LWT.stdt(:,id);
                mTtmp=mean(estm.LWT.T(:,:,id),2); sTtmp=std(estm.LWT.T(:,:,id),1,2);
                errorbar(m, ttmp, 1.96*sttmp, 'ro-'); hold on;
                errorbar(mref, mTtmp, 1.96*sTtmp, 'bx-'); hold off; grid on;
                set(gca, 'XTick', m1); set(gca, 'FontSize', FontSize1);
                xlabel('m'); legend('\theta_m','\theta_m^*');
                text('string',['c_',num2str((i)),' : Leaders'], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.03 0.95]);
                text('string',['d=',num2str(dtmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.87]);
                text('string',['p=',num2str(ptmp)], 'fontsize',FontSize1,'FontName', 'Helvetica', 'units','normalized', 'pos',[0.06 0.79]);
                if Ltmp>1; sb_fils=get(sb,'Children'); set(sb_fils,'hittest','off'); set(sb,'ButtonDownFcn','my_click_and_plot(gcbo)'); end;
            end
        end
    end
end
