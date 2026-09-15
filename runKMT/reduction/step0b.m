addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
T=readtable('/home/ocs/KMTdata/GaiaRef/gaia_dr3_ruwe_cone3arcmin.csv');
% the SAME cone gaiaTie uses: 120 arcsec, so GaiaInd indexes these rows
[GC,Col]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120);
cra=GC(:,strcmp(Col,'RA'))*180/pi; cdec=GC(:,strcmp(Col,'Dec'))*180/pi;
D=sqrt(((cra-T.ra').*cosd(cdec)).^2+(cdec-T.dec').^2)*3600;
[d,j]=min(D,[],2); ok=d<0.05;
RUWE=nan(size(GC,1),1); RUWE(ok)=T.ruwe(j(ok));
fprintf('gaiaTie cone (120"): %d rows | matched to the archive %d (%.0f%%) | RUWE finite %d\n', ...
    size(GC,1), sum(ok), 100*mean(ok), sum(isfinite(RUWE)));
save('/home/ocs/KMTdata/GaiaRef/ruwe_cone120.mat','RUWE','-v7.3');
for F={'BLG41','BLG01'}
  f=F{1};
  Cat=ml.util.ogleCompanionCat([Doc(1:end-4) 'OGLEdata/OB260058/OB160058.mat'],'Field',f);
  S=load(sprintf('/home/ocs/KMTdata/Results/v3_Final/IFfinal_260058_CTIO_%s.mat',f));
  G2=load(sprintf('/home/ocs/KMTdata/Results/v3_Final/GaiaTie_%s.mat',f));
  IF=S.IFsys; NF=S.Info; Tie=G2.Tie;
  Iog=NF.SrcData.I_ogle(:); Iog(Iog>=99)=NaN; Ie=IF.findClosestSource([150 150]);
  Cal=logical(ml.util.selectRefSources(IF,'RefMag',Iog,'MagRange',[14 17], ...
      'CompanionRadius',5,'CompanionMaxMag',18,'CompanionCat',Cat)); Cal=Cal(:); Cal(Ie)=true;
  gm=Tie.Matched>0 & isfinite(Tie.GaiaInd);
  r=nan(IF.Nsrc,1); r(gm)=RUWE(Tie.GaiaInd(gm));
  clean=gm & isfinite(r) & r<1.4;
  fprintf('%s: calib %d | Gaia-matched %d | RUWE present %d | RUWE<1.4 %d -> FRAME STARS %d\n', ...
      f, sum(Cal), sum(Cal&gm), sum(Cal&gm&isfinite(r)), sum(Cal&clean), sum(Cal&clean));
  clear IF S
end
