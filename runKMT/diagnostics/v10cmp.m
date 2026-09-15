% v8 vs v10 on the SAME stars (v10's set), so the RMS comparison is not a selection effect
addpath('/home/ocs/matlab/Lizoo');
for F={'BLG41','BLG01'}
  f=F{1};
  P8=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',f)); P10=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep10_%s.mat',f));
  V8=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info');
  V10=load(sprintf('/home/ocs/KMTdata/Results/v10/IFfinal_%s.mat',f),'IFsys','Info');
  m8=P8.SrcIdx(V8.Info.SrcInd(:)); m10=P10.SrcIdx(V10.Info.SrcInd(:));
  [~,i8,i10]=intersect(m8,m10);
  I=V10.Info.SrcData.I_ogle(:); I(I>=99)=NaN; I=I(i10);
  [x8,y8]=V8.IFsys.calculateRstd; [x10,y10]=V10.IFsys.calculateRstd;
  x8=x8(i8); y8=y8(i8); x10=x10(i10); y10=y10(i10);
  fprintf('%s: %d common stars\n', f, numel(i8));
  for e=[14 17; 17 18; 18 19]'
    s=I>=e(1)&I<e(2);
    fprintf('   I %2d-%2d (n=%3d): v8 %.3f/%.3f -> v10 %.3f/%.3f mas  (%+.1f%% / %+.1f%%)\n', e(1),e(2),sum(s), ...
      median(x8(s)),median(y8(s)),median(x10(s)),median(y10(s)), 100*(median(x10(s))/median(x8(s))-1), 100*(median(y10(s))/median(y8(s))-1));
  end
  s=I<17; fprintf('   stars improved in X: %.0f%%, in Y: %.0f%% (I<17)\n', 100*mean(x10(s)<x8(s)), 100*mean(y10(s)<y8(s)));
end
