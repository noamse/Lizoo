% (mu - mu_Gaia)_alpha against (mu - mu_Gaia)_delta, same colour code as the report's Figure 1
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120);
g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec'); ePMRA=g('ErrPMRA'); ePMDec=g('ErrPMDec'); Gmag=g('phot_g_mean_mag');
Rw=load('/home/ocs/KMTdata/GaiaRef/ruwe_cone120.mat'); RUWE=Rw.RUWE;
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat'); Fs={'BLG41','BLG01'};
SigInt=struct('BLG41',[0.347 0.358],'BLG01',[0.227 0.376]);
Fh=figure('Visible','off','Position',[40 40 1400 700]);
TL=tiledlayout(Fh,1,2,'TileSpacing','compact','Padding','compact');
for K=1:2
  f=Fs{K}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info');
  T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f)); Ie=C.Sel.(['Ie_' f]);
  I=V.Info.SrcData.I_ogle(:); I(I>=99)=NaN;
  J=T.J; gm=T.Tie.Matched>0&isfinite(J); cl=T.clean; ts=T.tieSet;
  dA=nan(V.IFsys.Nsrc,1); dD=dA; dA(gm)=T.PMabs(1,gm)'-PMRA(J(gm)); dD(gm)=T.PMabs(2,gm)'-PMDec(J(gm));
  nexttile(TL,K); hold on; box on; grid on; axis equal; set(gca,'FontSize',12);
  sa=std(dA(cl),'omitnan'); sd=std(dD(cl),'omitnan'); th=linspace(0,2*pi,200);   % plain std, no outlier rejection
  for k=[1 2 3], plot(k*sa*cos(th),k*sd*sin(th),'-','Color',[0.5 0.5 0.5],'LineWidth',0.8); end
  h0=plot(dA(gm&~cl),dD(gm&~cl),'x','MarkerSize',5,'Color',[0.6 0.6 0.6]);
  h1=plot(dA(cl),dD(cl),'o','MarkerSize',4,'MarkerFaceColor',[0.35 0.60 0.85],'MarkerEdgeColor',[0.1 0.3 0.6]);
  h2=plot(dA(ts),dD(ts),'o','MarkerSize',5,'MarkerFaceColor',[0.2 0.7 0.35],'MarkerEdgeColor','k');
  xline(0,'k-'); yline(0,'k-');
  si=SigInt.(f); eg=[ePMRA(J(Ie)) ePMDec(J(Ie))]; et=sqrt(si.^2+eg.^2);
  plot(dA(Ie)+et(1)*cos(th),dD(Ie)+et(2)*sin(th),'-','Color',[0.85 0.1 0.1],'LineWidth',1.2);
  col=[0.85 0.1 0.1]; if ~cl(Ie), col=[0.95 0.6 0.1]; end
  h3=plot(dA(Ie),dD(Ie),'p','MarkerSize',18,'MarkerFaceColor',col,'MarkerEdgeColor','k');
  % where does the target sit in the clean population? normalised radius
  r=sqrt((dA/sa).^2+(dD/sd).^2); pc=100*mean(r(cl)<r(Ie));
  m=sum(cl); wm=cl&abs(I-I(Ie))<0.5; pcm=100*mean(r(wm)<r(Ie));
  xlabel('(\mu - \mu_{Gaia})_\alpha cos\delta [mas/yr]','FontSize',13); ylabel('(\mu - \mu_{Gaia})_\delta [mas/yr]','FontSize',13);
  L=ceil(max([3.2*max(sa,sd); max(abs(dA(cl))); max(abs(dD(cl)))])); axis([-L L -L L]);
  o=find(cl&(abs(dA)>3*sa|abs(dD)>3*sd)); fprintf('   beyond 3 std: %s\n', strjoin(arrayfun(@(i) sprintf('I=%.1f (%+.1f,%+.1f)',I(i),dA(i),dD(i)),o,'uni',0),'; '));
  title(sprintf(['%s: target \\Delta\\mu = %+.2f / %+.2f, Gaia G=%.1f RUWE=%.2f\n', ...
      'population std %.2f / %.2f; target at %.0f%% of %d clean stars, %.0f%% of %d within 0.5 mag'], ...
      f, dA(Ie), dD(Ie), Gmag(J(Ie)), RUWE(J(Ie)), sa, sd, pc, m, pcm, sum(wm)),'FontWeight','bold','FontSize',12);
  if K==1, legend([h1 h2 h3 h0],{'RUWE<1.4 (gauge fitted on these)','tie set: RUWE<1.4 and 15<I<17', ...
     'the target (circle: internal \oplus Gaia error)','Gaia-matched, RUWE\geq1.4'},'Location','southwest','FontSize',10); end
  fprintf('%s: target dmu %+.3f / %+.3f, Gaia err %.3f / %.3f, internal %.3f / %.3f -> %.1f / %.1f sigma; G %.2f RUWE %.2f; clean %d; pop std %.3f/%.3f; pct %.0f%% (all) %.0f%% (within 0.5 mag, n=%d)\n', ...
     f, dA(Ie), dD(Ie), eg, si, dA(Ie)/et(1), dD(Ie)/et(2), Gmag(J(Ie)), RUWE(J(Ie)), m, sa, sd, pc, pcm, sum(wm));
  fprintf('   Gaia PM of the target: %+.3f +- %.3f / %+.3f +- %.3f ; ours %+.3f / %+.3f\n', PMRA(J(Ie)), eg(1), PMDec(J(Ie)), eg(2), T.PMabs(1,Ie), T.PMabs(2,Ie));
end
title(TL,'v8: our absolute proper motion minus Gaia DR3, per star (grey ellipses: 1, 2, 3 plain standard deviations of the clean population, no outlier rejection)','FontWeight','bold','FontSize',14);
exportgraphics(Fh,[Doc 'qa_v8_pmdiff_vs_gaia_std.png'],'Resolution',150); close(Fh);
