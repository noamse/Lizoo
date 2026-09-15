addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/'; AxN='XY';
RAd=celestial.coo.convertdms('17:52:38.09','gH','d');
Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120);
g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec');
Fs={'BLG41','BLG01'};

%% ---- 1. RMS vs magnitude -------------------------------------------------
Fh=figure('Visible','off','Position',[40 40 1400 1000]);
TL=tiledlayout(Fh,2,2,'TileSpacing','compact','Padding','compact');
for K=1:2
  f=Fs{K};
  V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v13/Tie_%s.mat',f));
  IF=V.IFsys; NF=V.Info; I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  [rx,ry]=IF.calculateRstd; Ie=IF.findClosestSource([150 150]);
  for A=1:2
    nexttile(TL,(A-1)*2+K); hold on; box on; grid on; set(gca,'FontSize',12);
    if A==1, R=rx; else, R=ry; end
    h1=plot(I,R,'o','MarkerSize',4,'MarkerFaceColor',[0.35 0.60 0.85],'MarkerEdgeColor',[0.1 0.3 0.6]);
    h2=plot(I(T.tieSet),R(T.tieSet),'o','MarkerSize',5,'MarkerFaceColor',[0.2 0.7 0.35],'MarkerEdgeColor','k');
    h3=plot(I(Ie),R(Ie),'p','MarkerSize',18,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
    % running median
    ed=14:0.5:19; cc=[]; mm=[];
    for q=1:numel(ed)-1
      s=I>=ed(q)&I<ed(q+1)&isfinite(R);
      if sum(s)>=5, cc(end+1)=mean(ed(q:q+1)); mm(end+1)=median(R(s)); end
    end
    h4=plot(cc,mm,'k-','LineWidth',2);
    set(gca,'YScale','log'); ylim([2 200]); xlim([13.5 19.5]);
    xlabel('OGLE I [mag]','FontSize',13);
    ylabel(sprintf('\\Delta%s residual RMS [mas]',AxN(A)),'FontSize',13);
    title(sprintf('%s   \\Delta%s',f,AxN(A)),'FontWeight','bold','FontSize',14);
    if K==1&&A==1
      legend([h1 h2 h3 h4],{'calibration stars (all in the fit)', ...
        'the Gaia tie set: RUWE<1.4 and 15<I<17','the target','running median'}, ...
        'Location','northwest','FontSize',10);
    end
  end
  clear IF V
end
title(TL,'Residual RMS after position and proper motion, per axis and field', ...
   'FontWeight','bold','FontSize',15);
exportgraphics(Fh,[Doc 'report_v13_RMS.png'],'Resolution',150); close(Fh);

%% ---- 2. our PM against Gaia ----------------------------------------------
Fh=figure('Visible','off','Position',[40 40 1400 1000]);
TL=tiledlayout(Fh,2,2,'TileSpacing','compact','Padding','compact');
for K=1:2
  f=Fs{K};
  V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v13/Tie_%s.mat',f));
  IF=V.IFsys; cl=T.clean; J=T.J; Ie=IF.findClosestSource([150 150]);
  Ours={T.PMabs(1,:)',T.PMabs(2,:)'}; Gaia={PMRA(J),PMDec(J)};
  nm={'\mu_\alpha cos\delta','\mu_\delta'};
  for A=1:2
    nexttile(TL,(A-1)*2+K); hold on; box on; grid on; set(gca,'FontSize',12);
    plot([-20 20],[-20 20],'k--','LineWidth',1);
    plot(Gaia{A}(cl),Ours{A}(cl),'o','MarkerSize',4, ...
        'MarkerFaceColor',[0.35 0.60 0.85],'MarkerEdgeColor',[0.1 0.3 0.6]);
    plot(Gaia{A}(T.tieSet),Ours{A}(T.tieSet),'o','MarkerSize',5, ...
        'MarkerFaceColor',[0.2 0.7 0.35],'MarkerEdgeColor','k');
    if cl(Ie)
      plot(Gaia{A}(Ie),Ours{A}(Ie),'p','MarkerSize',18,'MarkerFaceColor',[0.85 0.1 0.1],'MarkerEdgeColor','k');
    else
      plot(Gaia{A}(Ie),Ours{A}(Ie),'p','MarkerSize',18,'MarkerFaceColor',[0.95 0.6 0.1],'MarkerEdgeColor','k');
    end
    d=Ours{A}(cl)-Gaia{A}(cl);
    axis([-18 12 -18 12]);
    xlabel(sprintf('Gaia %s [mas/yr]',nm{A}),'FontSize',13);
    ylabel(sprintf('v13 %s [mas/yr]',nm{A}),'FontSize',13);
    dt_=Ours{A}(T.tieSet)-Gaia{A}(T.tieSet);
    title(sprintf('%s   %s   —  scatter %.2f (all %d), %.2f (tie set %d)', f, nm{A}, ...
        tools.math.stat.rstd(d), sum(cl), tools.math.stat.rstd(dt_), sum(T.tieSet)), ...
        'FontWeight','bold','FontSize',12);
  end
  clear IF V
end
title(TL,'Absolute proper motion against Gaia, after the affine tie (full gauge removed)', ...
   'FontWeight','bold','FontSize',15);
exportgraphics(Fh,[Doc 'report_v13_pm_vs_gaia.png'],'Resolution',150); close(Fh);
fprintf('V13F1 DONE\n');
