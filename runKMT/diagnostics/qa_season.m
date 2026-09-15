addpath('/home/ocs/matlab/Lizoo'); C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys'); IF=V.IFsys; Ie=C.Sel.(['Ie_' f]);
  T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f),'Amap');
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD)|isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup; t=(IF.JD-IF.JD0)/365.25; ok=W(:,Ie)>0;
  S=1000*(T.Amap\[Rx(:,Ie)';Ry(:,Ie)']); ra=S(1,:)'; de=S(2,:)';   % sky, mas
  tb=nan(nb,1); mr=tb; md=tb; nn=tb; sr=tb; sd=tb;
  for b=1:nb, s=ok&G==b; nn(b)=sum(s); if nn(b)<20, continue; end
    tb(b)=mean(t(s)); mr(b)=mean(ra(s)); md(b)=mean(de(s)); sr(b)=std(ra(s))/sqrt(nn(b)); sd(b)=std(de(s))/sqrt(nn(b)); end
  k=isfinite(tb); A=[ones(sum(k),1) tb(k)]; er=mr(k)-A*(A\mr(k)); ed=md(k)-A*(A\md(k)); n=sum(k); Sxx=sum((tb(k)-mean(tb(k))).^2);
  fprintf('\n%s target, sky axes. seasons used %d, sqrt(Sxx) %.2f yr, span %.1f yr\n', f, n, sqrt(Sxx), max(tb)-min(tb));
  fprintf('   %-6s %5s %8s %8s %6s %6s\n','year','n','RA mean','Dec mean','seRA','seDec');
  for b=find(k)', fprintf('   %6.2f %5d %+8.2f %+8.2f %6.2f %6.2f\n', 2016+tb(b), nn(b), mr(b), md(b), sr(b), sd(b)); end
  fprintf('   scatter about the line: %.2f / %.2f mas | white-noise s.e. median %.2f / %.2f | ratio %.1f / %.1f\n', std(er,1)*sqrt(n/(n-2)), std(ed,1)*sqrt(n/(n-2)), median(sr(k)), median(sd(k)), std(er,1)*sqrt(n/(n-2))/median(sr(k)), std(ed,1)*sqrt(n/(n-2))/median(sd(k)));
  fprintf('   => PM error %.3f / %.3f mas/yr ; formula sigma/T*sqrt(12/N) with T=%.1f: %.3f / %.3f\n', std(er,1)*sqrt(n/(n-2))/sqrt(Sxx), std(ed,1)*sqrt(n/(n-2))/sqrt(Sxx), max(tb)-min(tb), std(er,1)*sqrt(n/(n-2))/(max(tb)-min(tb))*sqrt(12/n), std(ed,1)*sqrt(n/(n-2))/(max(tb)-min(tb))*sqrt(12/n));
end
