% Exactly Figure 1's points and colour code, on (mu-mu_Gaia)_alpha vs (mu-mu_Gaia)_delta
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120);
g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec');
Fs={'BLG41','BLG01'};
Fh=figure('Visible','off','Position',[40 40 1400 700]);
TL=tiledlayout(Fh,1,2,'TileSpacing','compact','Padding','compact');
for K=1:2
  f=Fs{K}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys');
  T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f)); IF=V.IFsys; cl=T.clean; J=T.J; Ie=IF.findClosestSource([150 150]);
  dA=T.PMabs(1,:)'-PMRA(J); dD=T.PMabs(2,:)'-PMDec(J);
  nexttile(TL,K); hold on; box on; grid on; axis equal; set(gca,'FontSize',12);
  plot([-8 8],[0 0],'k--','LineWidth',0.8); plot([0 0],[-8 8],'k--','LineWidth',0.8);
  h1=plot(dA(cl),dD(cl),'o','MarkerSize',4,'MarkerFaceColor',[0.35 0.60 0.85],'MarkerEdgeColor',[0.1 0.3 0.6]);
  h2=plot(dA(T.tieSet),dD(T.tieSet),'o','MarkerSize',5,'MarkerFaceColor',[0.2 0.7 0.35],'MarkerEdgeColor','k');
  if cl(Ie), fc=[0.85 0.1 0.1]; else, fc=[0.95 0.6 0.1]; end
  h3=plot(dA(Ie),dD(Ie),'p','MarkerSize',18,'MarkerFaceColor',fc,'MarkerEdgeColor','k');
  axis([-6 6 -6 6]);
  xlabel('(\mu - \mu_{Gaia})_\alpha cos\delta [mas/yr]','FontSize',13); ylabel('(\mu - \mu_{Gaia})_\delta [mas/yr]','FontSize',13);
  title(sprintf('%s   —  scatter %.2f / %.2f (all %d), %.2f / %.2f (tie set %d)', f, ...
     tools.math.stat.rstd(dA(cl)), tools.math.stat.rstd(dD(cl)), sum(cl), tools.math.stat.rstd(dA(T.tieSet)), tools.math.stat.rstd(dD(T.tieSet)), sum(T.tieSet)),'FontWeight','bold','FontSize',12);
  if K==1, legend([h1 h2 h3],{'RUWE < 1.4 (all stars of Figure 1)','the Gaia tie set: RUWE<1.4 and 15<I<17','the target'},'Location','southwest','FontSize',10); end
end
title(TL,'v8: absolute proper motion minus Gaia DR3, per star, after the affine tie (full gauge removed)','FontWeight','bold','FontSize',15);
exportgraphics(Fh,[Doc 'qa_v8_pmdiff_vs_gaia_plain.png'],'Resolution',150); close(Fh);
fprintf('PLAIN DONE\n');
