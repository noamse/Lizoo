addpath('/home/ocs/matlab/Lizoo');
Doc='/home/ocs/.claude/jobs/010f4e84/tmp/'; AxN='XY'; SidMonth=27.321661;
RUN=getenv('RUN'); FLD=getenv('FLD'); if isempty(RUN)||isempty(FLD), error('set RUN and FLD'); end
MinN=10; MaxSecz=1.3;   % bins with fewer epochs, and epochs above this airmass, are dropped
for K=1
  f=FLD;
  V=load(sprintf('/home/ocs/KMTdata/Results/%s/IFfinal_%s.mat',RUN,f));
  T=load(sprintf('/home/ocs/KMTdata/Results/%s/Tie_%s.mat',RUN,f));
  IF=V.IFsys; NF=V.Info; PMabs=T.PMabs;
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  P8=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',f)); Ie=find(P8.SrcIdx(NF.SrcInd(:))==P8.TargetMsc);
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  secz=IF.getTimeSeriesField(1,{'secz'});
  X0=IF.ParS(1,Ie); Y0=IF.ParS(2,Ie);
  List={Ie,'target','the target'};
  for L=1:size(List,1)
    Ii=List{L,1}; tag=List{L,2}; what=List{L,3};
    ok=W(:,Ii)>0 & secz<MaxSecz; if sum(ok)<50, continue; end
    JD=IF.JD(ok); t=(JD-IF.JD0)/365.25;
    dist=sqrt((IF.ParS(1,Ii)-X0)^2+(IF.ParS(2,Ii)-Y0)^2);
    pos=cell(1,2); linc=cell(1,2); res=cell(1,2); MuM=nan(1,2);
    for A=1:2
      if A==1, r=400*Rx(ok,Ii); p0=IF.ParS(1,Ii); mu=IF.ParS(3,Ii);
      else,    r=400*Ry(ok,Ii); p0=IF.ParS(2,Ii); mu=IF.ParS(4,Ii); end
      lin=400*(p0+mu.*t); pp=lin+r; off=mean(pp);
      pos{A}=pp-off; linc{A}=lin-off; res{A}=r; MuM(A)=400*mu;
    end
    lab=floor((JD-min(JD))/SidMonth); u=unique(lab);
    Tt=cell(1,2);Pp=cell(1,2);Ep=cell(1,2);Rr=cell(1,2);Er=cell(1,2);Ll=cell(1,2); Fe=nan(1,2);
    for A=1:2
      TT=nan(numel(u),1);MP=TT;EP=TT;MR=TT;ER=TT;ML=TT;
      for q=1:numel(u)
        s=lab==u(q); if sum(s)<MinN, continue; end
        TT(q)=mean(JD(s)); MP(q)=mean(pos{A}(s));
        MR(q)=mean(res{A}(s)); ER(q)=std(res{A}(s))/sqrt(sum(s)); ML(q)=mean(linc{A}(s));
      end
      g=isfinite(TT); Tt{A}=TT(g);Pp{A}=MP(g);Rr{A}=MR(g);Er{A}=ER(g);Ll{A}=ML(g);
      % s.e. scaled by the empirical factor sqrt(chi2/dof) of the bin means about zero;
      % the same bars serve the PM-retained row, whose bins hold the same epochs
      Fe(A)=max(1,sqrt(sum((Rr{A}./Er{A}).^2)/(numel(Rr{A})-1)));
      Er{A}=Er{A}*Fe(A); Ep{A}=Er{A};
    end
    hi=max([quantile(abs(Pp{1})+Ep{1},0.90),quantile(abs(Rr{1})+Er{1},0.90), ...
            quantile(abs(Pp{2})+Ep{2},0.90),quantile(abs(Rr{2})+Er{2},0.90)]);
    YL=max(20,10*ceil(1.1*hi/10));
    Fh=figure('Visible','off','Position',[40 40 1500 900]);
    TL=tiledlayout(Fh,2,2,'TileSpacing','compact','Padding','compact');
    for A=1:2
      D=datetime(Tt{A},'convertfrom','juliandate');
      nexttile(TL,A); hold on; box on; grid on; set(gca,'FontSize',13);
      errorbar(D,Pp{A},Ep{A},'o','Color',[0.10 0.35 0.65],'MarkerFaceColor',[0.35 0.60 0.85], ...
               'MarkerSize',4,'LineWidth',0.8,'CapSize',0);
      hL=plot(D,Ll{A},'-','Color',[0.75 0.10 0.10],'LineWidth',2);
      ylim([-YL YL]); ylabel(sprintf('%s position [mas]',AxN(A)),'FontSize',14);
      if A==1, cv='pixel X, which runs OPPOSITE to RA'; else, cv='pixel Y, aligned with Dec'; end
      legend(hL,sprintf('fitted PM %+.2f mas/yr along %s',MuM(A),cv), ...
             'Location','northwest','FontSize',11);
      title(sprintf('%s   position, proper motion RETAINED   (bars = s.e. \\times %.1f)',AxN(A),Fe(A)),'FontWeight','bold','FontSize',14);
      nexttile(TL,A+2); hold on; box on; grid on; set(gca,'FontSize',13);
      errorbar(D,Rr{A},Er{A},'o','Color',[0.65 0.10 0.10],'MarkerFaceColor',[0.95 0.45 0.35], ...
               'MarkerSize',4,'LineWidth',0.8,'CapSize',0);
      yline(0,'k-','LineWidth',1); ylim([-YL YL]);
      xlabel('date','FontSize',14); ylabel(sprintf('\\Delta%s residual [mas]',AxN(A)),'FontSize',14);
      title(sprintf('\\Delta%s   position and proper motion REMOVED   (bars = s.e. \\times %.1f)',AxN(A),Fe(A)),'FontWeight','bold','FontSize',14);
    end
    if Ii==Ie, wh=''; else, wh=sprintf(', %.0f pix from the target',dist); end
    title(TL,sprintf([RUN '  %s  %s  —  I = %.2f%s, at (%.0f,%.0f)   |   ', ...
        'ABSOLUTE sky PM %+.2f / %+.2f mas/yr   |   bins \\geq%d epochs, sec z < %.1f'], f, what, I(Ii), wh, ...
        IF.ParS(1,Ii), IF.ParS(2,Ii), PMabs(1,Ii), PMabs(2,Ii), MinN, MaxSecz), ...
        'FontWeight','bold','FontSize',14);
    exportgraphics(Fh,sprintf('%s%s_motion_%s_%s.png',Doc,RUN,f,tag),'Resolution',150);
    close(Fh);
    fprintf('%s %-9s I %.2f d %5.1f | PM %+6.2f/%+6.2f | factor %.2f / %.2f\n', f, tag, I(Ii), dist, PMabs(1,Ii), PMabs(2,Ii), Fe);
  end
  clear IF V
end
fprintf('RUNFIGS DONE\n');
