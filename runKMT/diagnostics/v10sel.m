% Second clipping iteration, on the season scale: drop stars whose season-offset
% chi2/DoF -- sum over seasons and axes of (season mean / s.e.)^2, divided by
% 2*Nseasons -- sits more than 2 sigma above the per-magnitude median in log space.
% Applied to the v8 solution (which is itself the decade-scale 2-sigma cut on v6).
% v8's Info.SrcInd indexes prep8's list, so the MSc index of v8 source i is
% prep8.SrcIdx(Info.SrcInd(i)).
addpath('/home/ocs/matlab/Lizoo');
for F={'BLG41','BLG01'}
  f=F{1};
  P8=load(sprintf('/home/ocs/KMTdata/GaiaRef/prep8_%s.mat',f));
  V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info');
  IF=V.IFsys; NF=V.Info; N=IF.Nsrc;
  mscAll=P8.SrcIdx(NF.SrcInd(:)); Ie=find(mscAll==P8.TargetMsc);
  if numel(Ie)~=1, error('v10sel:Target','target not found once in %s',f); end
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup;
  Ch=nan(N,nb);
  for i=1:N
    for b=1:nb
      s=W(:,i)>0&G==b; n=sum(s); if n<30, continue; end
      r=400*[Rx(s,i) Ry(s,i)]; mo=mean(r,1); se=std(r,0,1)/sqrt(n);
      Ch(i,b)=sum((mo./se).^2);
    end
  end
  chi=sum(Ch,2,'omitnan')./(2*sum(isfinite(Ch),2));
  L=log10(chi); ok=isfinite(I)&isfinite(L);
  keep=true(N,1); keep(~ok)=false; z=nan(N,1);
  ed=14:0.5:19; nrej=0;
  for q=1:numel(ed)-1
    s=ok&I>=ed(q)&I<ed(q+1);
    if sum(s)<6, continue; end
    m=median(L(s)); sd=tools.math.stat.rstd(L(s)); z(s)=(L(s)-m)/sd;
    bad=s&(L>m+2*sd); keep(bad)=false; nrej=nrej+sum(bad);
  end
  keep(Ie)=true;                           % never drop the target
  SrcIdx=mscAll(keep); SrcIdx=SrcIdx(:); TargetMsc=P8.TargetMsc;
  fprintf('\n=== %s ===\n',f);
  fprintf('  %d in, %d rejected above 2 sigma in season-offset chi2/DoF, %d kept\n', N, nrej, sum(keep));
  fprintf('  target: chi2/DoF %.1f, z %+.2f, entry %d of the new list\n', chi(Ie), z(Ie), find(SrcIdx==TargetMsc));
  rj=find(~keep); [~,o]=sort(z(rj),'descend'); rj=rj(o);
  fprintf('  rejected (I, chi2/DoF, z): %s\n', strjoin(arrayfun(@(i) sprintf('%.1f/%.0f/%.1f',I(i),chi(i),z(i)),rj,'uni',0),'  '));
  save(sprintf('/home/ocs/KMTdata/GaiaRef/prep10_%s.mat',f),'SrcIdx','keep','TargetMsc','chi','z','-v7.3');
  clear IF V
end
