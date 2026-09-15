% Iterate the calibrator selection: drop stars whose residual RMS sits more than
% 2 sigma above the running median for their magnitude, in log space so the cut
% is multiplicative and does not simply shave off the faint end.
%
% NOTE on indexing. v6.Info.SrcInd indexes v6's OWN INPUT LIST (the prep6
% selection), not the full MSc. The MSc index of v6 source i is therefore
% prep6.SrcIdx(v6.Info.SrcInd(i)). Using SrcInd directly as an MSc index picks
% the wrong stars entirely.
addpath('/home/ocs/matlab/Lizoo');
for F={'BLG41','BLG01'}
  f=F{1};
  P6=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep6_%s.mat',f));
  V=load(sprintf('/home/ocs/KMTdata/Results/v6/IFfinal_%s.mat',f));
  IF=V.IFsys; NF=V.Info; Ie=IF.findClosestSource([150 150]);
  mscAll = P6.SrcIdx(NF.SrcInd(:));        % true MSc index of every v6 source
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  [rx,ry]=IF.calculateRstd; R=sqrt(rx(:).^2+ry(:).^2);
  ok=isfinite(I)&isfinite(R)&R>0; L=log10(R);
  keep=true(IF.Nsrc,1); keep(~ok)=false;
  ed=14:0.5:19; nrej=0;
  for q=1:numel(ed)-1
    s=ok&I>=ed(q)&I<ed(q+1);
    if sum(s)<6, continue; end
    m=median(L(s)); sd=tools.math.stat.rstd(L(s));
    bad=s&(L>m+2*sd); keep(bad)=false; nrej=nrej+sum(bad);
  end
  keep(Ie)=true;                            % never drop the target
  SrcIdx=mscAll(keep); SrcIdx=SrcIdx(:);
  TargetMsc=mscAll(Ie);
  fprintf('\n=== %s ===\n',f);
  fprintf('  %d in, %d rejected above 2 sigma, %d kept\n', IF.Nsrc, nrej, sum(keep));
  fprintf('  MSc indices span %d..%d (full MSc, as required)\n', min(SrcIdx), max(SrcIdx));
  fprintf('  target MSc index %d, at entry %d of the new list, RMS %.1f mas\n', ...
      TargetMsc, find(SrcIdx==TargetMsc), R(Ie));
  save(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',f),'SrcIdx','keep','TargetMsc','-v7.3');
  clear IF V
end
