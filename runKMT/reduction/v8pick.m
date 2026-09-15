% The same six comparison stars, located in v8 by OGLE identity so they are the
% same physical objects in both fields and comparable with the v6 report.
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
OGf=[Doc(1:end-4) 'OGLEdata/OB260058/OB160058.mat'];
C6=load('/home/ocs/KMTdata/Results/v6/cmp.mat'); S6=C6.Sel;
Sel=struct('BLG41',[],'BLG01',[],'tags',{S6.tags});
for F={'BLG41','BLG01'}
  f=F{1};
  V6=load(sprintf('/home/ocs/KMTdata/Results/v6/IFfinal_%s.mat',f));
  V8=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f));
  P8=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',f));
  I8=V8.Info.SrcData.I_ogle(:); I8(I8>=99)=NaN;
  X6=V6.IFsys.ParS(1,:)'; Y6=V6.IFsys.ParS(2,:)';
  X8=V8.IFsys.ParS(1,:)'; Y8=V8.IFsys.ParS(2,:)';
  idx=nan(1,numel(S6.tags));
  fprintf('\n=== %s ===\n',f);
  for q=1:numel(S6.tags)
    i6=S6.(f)(q);
    d=sqrt((X8-X6(i6)).^2+(Y8-Y6(i6)).^2);
    [dm,i8]=min(d);
    if dm<1.0, idx(q)=i8;
      fprintf('  %-9s v6 idx %3d (%5.1f,%5.1f) -> v8 idx %3d (%5.1f,%5.1f)  d=%.2f pix  I=%.2f\n', ...
          S6.tags{q}, i6, X6(i6), Y6(i6), i8, X8(i8), Y8(i8), dm, I8(i8));
    else
      fprintf('  %-9s NOT FOUND in v8 (nearest %.2f pix) -- removed by the outlier cut\n', S6.tags{q}, dm);
    end
  end
  Sel.(f)=idx;
  % the target
  Ie=find(P8.SrcIdx(V8.Info.SrcInd(:))==P8.TargetMsc);
  Sel.(['Ie_' f])=Ie;
  fprintf('  target v8 idx %d at (%.1f,%.1f)\n', Ie, X8(Ie), Y8(Ie));
  clear V6 V8
end
save('/home/ocs/KMTdata/Results/v8/cmp.mat','Sel','-v7.3');
