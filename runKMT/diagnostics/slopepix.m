% Pixel-phase test: does a star's within-season slope track the within-season
% drift of its own sub-pixel phase? The fit removes a GLOBAL quintic in phase,
% so any star-to-star variation in intra-pixel response survives.
addpath('/home/ocs/matlab/Lizoo');
for F={'BLG41','BLG01'}
  f=F{1};
  V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f));
  IF=V.IFsys;
  hasP = isfield(IF.Data,'Xphase') && ~isempty(IF.Data.Xphase);
  fprintf('\n=== %s ===  pixel phase stored: %d\n', f, hasP);
  if ~hasP, continue; end
  [Rx,~]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx))=0;
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup;
  Xp=IF.Data.Xphase;
  sl=[]; dph=[];
  for i=1:IF.Nsrc
    for b=1:nb
      s=W(:,i)>0 & G==b & isfinite(Xp(:,i)); n=sum(s);
      if n<50, continue; end
      t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t)); A=[ones(n,1) t];
      c1=A\(400*Rx(s,i)); c2=A\Xp(s,i);
      sl(end+1,1)=c1(2); dph(end+1,1)=c2(2);
    end
  end
  g=isfinite(sl)&isfinite(dph);
  fprintf('  %d star-season pairs\n', sum(g));
  fprintf('  corr(within-season slope, within-season pixel-phase drift) = %+.3f (Spearman %+.3f)\n', ...
      corr(sl(g),dph(g)), corr(sl(g),dph(g),'type','Spearman'));
  c=polyfit(dph(g),sl(g),1);
  r=sl(g)-polyval(c,dph(g));
  fprintf('  slope of the relation %.2f mas per unit phase drift; variance explained %.1f%%\n', ...
      c(1), 100*(1-var(r)/var(sl(g))));
  clear IF V
end
