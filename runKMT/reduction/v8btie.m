% Gaia tie for v8: fitted on RUWE<1.4 AND 15<I<17, full six-parameter gauge.
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
RUN=getenv('RUN'); if isempty(RUN), error('set RUN'); end
OGf=[Doc(1:end-4) 'OGLEdata/OB260058/OB160058.mat'];
RAd=celestial.coo.convertdms('17:52:38.09','gH','d');
Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
Rw=load('/home/ocs/KMTdata/GaiaRef/ruwe_cone120.mat'); RUWE=Rw.RUWE;
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120);
g=@(n) GC(:,strcmp(GCol,n));
Xi=(g('RA')*180/pi-RAd).*cosd(Decd).*3600; Eta=(g('Dec')*180/pi-Decd).*3600;
PMRA=g('PMRA'); PMDec=g('PMDec'); dt=2019.4127-2016.0;
for F={'BLG41','BLG01'}
  f=F{1};
  V=load(sprintf('/home/ocs/KMTdata/Results/%s/IFfinal_%s.mat',RUN,f));
  IF=V.IFsys; NF=V.Info;
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  Tie=ml.util.gaiaTie(IF,'Field',f,'OgleFile',OGf,'Verbosity',0);
  J=Tie.GaiaInd; gm=Tie.Matched>0 & isfinite(J);
  r=nan(IF.Nsrc,1); r(gm)=RUWE(J(gm));
  cleanAll = gm & isfinite(r) & r<1.4;
  tieSet   = cleanAll & I>=15 & I<=17;
  ii=find(tieSet);
  xi=Xi(J(ii))+PMRA(J(ii))/1000*dt; eta=Eta(J(ii))+PMDec(J(ii))/1000*dt;
  Xk=IF.ParS(1,ii)'; Yk=IF.ParS(2,ii)'; Gm=[ones(numel(ii),1) xi eta];
  ax=Gm\Xk; ay=Gm\Yk;
  for it=1:3
     rx=Xk-Gm*ax; ry=Yk-Gm*ay;
     good=abs(rx)<3*tools.math.stat.rstd(rx)&abs(ry)<3*tools.math.stat.rstd(ry);
     ax=Gm(good,:)\Xk(good); ay=Gm(good,:)\Yk(good);
  end
  Amap=[ax(2) ax(3); ay(2) ay(3)];
  mres=median(sqrt((Xk-Gm*ax).^2+(Yk-Gm*ay).^2)*400);
  sky=1000*(Amap\IF.ParS(3:4,:));
  Xc=IF.ParS(1,:)'-150; Yc=IF.ParS(2,:)'-150;
  Ag=[ones(numel(ii),1) Xc(ii) Yc(ii)];
  dxg=sky(1,ii)'-PMRA(J(ii)); dyg=sky(2,ii)'-PMDec(J(ii));
  cx=Ag\dxg; cy=Ag\dyg;
  for it=1:3
     rx=dxg-Ag*cx; ry=dyg-Ag*cy;
     kg=abs(rx)<3*tools.math.stat.rstd(rx)&abs(ry)<3*tools.math.stat.rstd(ry);
     cx=Ag(kg,:)\dxg(kg); cy=Ag(kg,:)\dyg(kg);
  end
  Aall=[ones(IF.Nsrc,1) Xc Yc];
  PMabs=[sky(1,:)'-Aall*cx, sky(2,:)'-Aall*cy]';
  Gauge=[cx cy]; clean=cleanAll;
  save(sprintf('/home/ocs/KMTdata/Results/%s/Tie_%s.mat',RUN,f), ...
       'Tie','Amap','ax','ay','PMabs','clean','tieSet','Gauge','J','-v7.3');
  te=find(cleanAll);
  fprintf('%s: Gaia-matched %d of %d | RUWE<1.4 %d | tie set (15<I<17) %d\n', ...
      f, sum(gm), IF.Nsrc, sum(cleanAll), sum(tieSet));
  fprintf('   map residual %.1f mas, scale %.4f pix/arcsec | gauge const %+.3f/%+.3f\n', ...
      mres, sqrt(abs(det(Amap))), cx(1), cy(1));
  fprintf('   scatter about Gaia over all RUWE-clean: %.3f / %.3f mas/yr\n', ...
      tools.math.stat.rstd(PMabs(1,te)'-PMRA(J(te))), tools.math.stat.rstd(PMabs(2,te)'-PMDec(J(te))));
  clear IF V
end
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/%s/IFfinal_%s.mat',RUN,f));
  T=load(sprintf('/home/ocs/KMTdata/Results/%s/Tie_%s.mat',RUN,f));
  Ie=V.IFsys.findClosestSource([150 150]);
  fprintf('%s %s TARGET PM (Gaia frame): %+.3f / %+.3f mas/yr\n', RUN, f, T.PMabs(1,Ie), T.PMabs(2,Ie));
end
