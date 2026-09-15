addpath('/home/ocs/matlab/Lizoo');
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120); g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec'); eD=g('ErrPMDec');
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f)); V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'Info'); Ie=C.Sel.(['Ie_' f]);
  I=V.Info.SrcData.I_ogle(:); I(I>=99)=NaN; cl=T.clean; J=T.J;
  dA=nan(size(I)); dD=dA; dA(cl)=T.PMabs(1,cl)'-PMRA(J(cl)); dD(cl)=T.PMabs(2,cl)'-PMDec(J(cl));
  for e=[14 17;17 18;18 19]'
    s=cl&I>=e(1)&I<e(2);
    fprintf('%s I %d-%d n=%3d: rstd dRA %.2f dDec %.2f | median Gaia ErrPMDec %.2f | frac |dDec|>2.0: %.0f%%\n', f, e, sum(s), tools.math.stat.rstd(dA(s)), tools.math.stat.rstd(dD(s)), median(eD(J(s))), 100*mean(abs(dD(s))>2.0));
  end
  s=cl&abs(I-I(Ie))<0.5; fprintf('%s within 0.5 mag n=%d: rstd dDec %.2f, target dDec %+.2f = %.1f sigma of that group; frac |dDec|>=target %.0f%%\n', f, sum(s), tools.math.stat.rstd(dD(s)), dD(Ie), dD(Ie)/tools.math.stat.rstd(dD(s)), 100*mean(abs(dD(s))>=abs(dD(Ie))));
end
