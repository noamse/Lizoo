addpath('/home/ocs/matlab/Lizoo');
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120);
g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec'); ePMRA=g('ErrPMRA'); ePMDec=g('ErrPMDec'); Gm=g('phot_g_mean_mag');
Xi=(g('RA')*180/pi-RAd).*cosd(Decd).*3600; Eta=(g('Dec')*180/pi-Decd).*3600;
Rw=load('/home/ocs/KMTdata/GaiaRef/ruwe_cone120.mat'); RUWE=Rw.RUWE;
% Gaia sources within 4 arcsec of the target
d=sqrt(Xi.^2+Eta.^2); near=find(d<4); [~,o]=sort(d(near)); near=near(o);
fprintf('Gaia DR3 sources within 4 arcsec of the target position:\n');
for k=near', fprintf('   d=%.2f arcsec  dxi %+.2f deta %+.2f  G=%.2f  PM %+.2f+-%.2f / %+.2f+-%.2f  RUWE %.2f\n', d(k), Xi(k), Eta(k), Gm(k), PMRA(k), ePMRA(k), PMDec(k), ePMDec(k), RUWE(k)); end
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f)); V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info');
  Ie=C.Sel.(['Ie_' f]); Ic=C.Sel.(f)(strcmp(C.Sel.tags,'m18_d04')); IF=V.IFsys; I=V.Info.SrcData.I_ogle(:);
  M=IF.medianFieldSource({'MAG_PSF'}); M=M(:);
  fprintf('%s: target ours %+.2f/%+.2f (I_ogle %.2f, KMT mag %.2f) Gaia idx %d | companion ours %+.2f/%+.2f (I_ogle %.2f, KMT mag %.2f) Gaia idx %d, sep %.1f pix = %.1f arcsec | fwhm median %.1f pix\n', ...
    f, T.PMabs(1,Ie), T.PMabs(2,Ie), I(Ie), M(Ie), T.J(Ie), T.PMabs(1,Ic), T.PMabs(2,Ic), I(Ic), M(Ic), T.J(Ic), ...
    norm(IF.ParS(1:2,Ie)-IF.ParS(1:2,Ic)), 0.4*norm(IF.ParS(1:2,Ie)-IF.ParS(1:2,Ic)), median(IF.Data.fwhm(:),'omitnan'));
  if isfinite(T.J(Ic)), fprintf('   companion Gaia PM %+.2f+-%.2f / %+.2f+-%.2f, G %.2f\n', PMRA(T.J(Ic)), ePMRA(T.J(Ic)), PMDec(T.J(Ic)), ePMDec(T.J(Ic)), Gm(T.J(Ic))); end
  % if both KMT centroids are flux-weighted blends: ours_t = (1-w) mu_t + w mu_c ; solve w from Dec
  if isfinite(T.J(Ic))
    mt=PMDec(T.J(Ie)); mc=PMDec(T.J(Ic)); w=(T.PMabs(2,Ie)-mt)/(mc-mt); wc=(T.PMabs(2,Ic)-mc)/(mt-mc);
    fprintf('   Dec: blend weight needed to explain ours: target pulled %.0f%% toward companion; companion pulled %.0f%% toward target\n', 100*w, 100*wc);
    fprintf('   RA check with those weights: predicted %+.2f (ours %+.2f) and %+.2f (ours %+.2f)\n', (1-w)*PMRA(T.J(Ie))+w*PMRA(T.J(Ic)), T.PMabs(1,Ie), (1-wc)*PMRA(T.J(Ic))+wc*PMRA(T.J(Ie)), T.PMabs(1,Ic));
  end
end
