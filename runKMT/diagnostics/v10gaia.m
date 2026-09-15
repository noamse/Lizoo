addpath('/home/ocs/matlab/Lizoo');
RAd=celestial.coo.convertdms('17:52:38.09','gH','d'); Decd=celestial.coo.convertdms('-31:47:36.1','gD','d');
[GC,GCol]=catsHTM.cone_search('GAIADR3',RAd/180*pi,Decd/180*pi,120); g=@(n) GC(:,strcmp(GCol,n)); PMRA=g('PMRA'); PMDec=g('PMDec');
for F={'BLG41','BLG01'}
  f=F{1};
  P8=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',f)); P10=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep10_%s.mat',f));
  V8=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'Info'); V10=load(sprintf('/home/ocs/KMTdata/Results/v10/IFfinal_%s.mat',f),'Info');
  T8=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f)); T10=load(sprintf('/home/ocs/KMTdata/Results/v10/Tie_%s.mat',f));
  m8=P8.SrcIdx(V8.Info.SrcInd(:)); m10=P10.SrcIdx(V10.Info.SrcInd(:)); [~,i8,i10]=intersect(m8,m10);
  c=T8.clean(i8)&T10.clean(i10); a=i8(c); b=i10(c);
  fprintf('%s: %d common RUWE-clean stars | v8 %.3f/%.3f -> v10 %.3f/%.3f mas/yr\n', f, sum(c), ...
    tools.math.stat.rstd(T8.PMabs(1,a)'-PMRA(T8.J(a))), tools.math.stat.rstd(T8.PMabs(2,a)'-PMDec(T8.J(a))), ...
    tools.math.stat.rstd(T10.PMabs(1,b)'-PMRA(T10.J(b))), tools.math.stat.rstd(T10.PMabs(2,b)'-PMDec(T10.J(b))));
end
