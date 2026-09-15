addpath('/home/ocs/matlab/Lizoo'); SidMonth=27.321661;
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys'); IF=V.IFsys; Ie=C.Sel.(['Ie_' f]);
  secz=IF.getTimeSeriesField(1,{'secz'});
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD)|isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0; ok=W(:,Ie)>0;
  fprintf('%s: secz quantiles 50/90/95/99%%: %s ; fraction >1.2 %.1f%% >1.25 %.1f%% >1.3 %.1f%%\n', f, mat2str(quantile(secz(ok),[.5 .9 .95 .99]),3), 100*mean(secz(ok)>1.2),100*mean(secz(ok)>1.25),100*mean(secz(ok)>1.3));
  JD=IF.JD; T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f),'Amap');
  S=1000*(T.Amap\[Rx(:,Ie)';Ry(:,Ie)']);
  for cap=[Inf 1.3 1.25 1.2]
    o=ok&secz<cap; lab=floor((JD(o)-min(JD(ok)))/SidMonth); u=unique(lab); n=arrayfun(@(k) sum(lab==k),u);
    ra=S(1,o)'; de=S(2,o)'; MR=arrayfun(@(k) mean(ra(lab==k)),u); MD=arrayfun(@(k) mean(de(lab==k)),u);
    g=n>=10; dev=sqrt((MR(g)-mean(MR(g))).^2+(MD(g)-mean(MD(g))).^2);
    fprintf('   cap %4.2f n>=10: %d bins (%d dropped for n<10) | bins >15 mas: %d | rms RA/Dec %.2f/%.2f\n', cap, sum(g), sum(~g), sum(dev>15), std(MR(g)), std(MD(g)));
  end
end
