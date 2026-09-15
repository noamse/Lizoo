C=load('/home/ocs/KMTdata/Results/v13/cmp.mat'); S=C.Sel;
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f),'IFsys'); [a,b]=V.IFsys.calculateRstd;
  Ie=S.(['Ie_' f]); fprintf('%s target %.2f / %.2f\n', f, a(Ie), b(Ie));
  for q=1:numel(S.tags), i=double(S.(f)(q)); if isfinite(i), fprintf('   %-8s %.2f / %.2f\n', S.tags{q}, a(i), b(i)); end; end
end
