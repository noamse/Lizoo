% remaining outlier bins after the n>=10, secz<1.3 cuts: what are they made of?
addpath('/home/ocs/matlab/Lizoo'); SidMonth=27.321661; MinN=10; MaxSecz=1.3;
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info');
  T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f),'Amap');
  IF=V.IFsys; Ie=C.Sel.(['Ie_' f]); I=V.Info.SrcData.I_ogle(:); I(I>=99)=NaN;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD)|isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  secz=IF.getTimeSeriesField(1,{'secz'}); JD=IF.JD; night=floor(JD-0.5);
  ok=W(:,Ie)>0 & secz<MaxSecz;
  lab=floor((JD-min(JD(ok)))/SidMonth);
  % field-wide per-epoch median residual of stars 14<I<17 (frame-level systematic left after the fit)
  br=find(I>14&I<17); 
  Wb=W(:,br)>0; rxb=Rx(:,br); ryb=Ry(:,br); rxb(~Wb)=NaN; ryb(~Wb)=NaN;
  fx=400*median(rxb,2,'omitnan'); fy=400*median(ryb,2,'omitnan');
  px=400*Rx(:,Ie); py=400*Ry(:,Ie);
  u=unique(lab(ok)); R=[];
  for q=1:numel(u)
    s=ok&lab==u(q); n=sum(s); if n<MinN, continue; end
    mx=mean(px(s)); my=mean(py(s)); mdx=median(px(s)); mdy=median(py(s));
    R(end+1,:)=[mean(JD(s)) n numel(unique(night(s))) mx my mdx mdy std(px(s)) std(py(s)) mean(fx(s)) mean(fy(s))];
  end
  dev=sqrt(R(:,4).^2+R(:,5).^2); out=dev>12;
  fprintf('\n%s target: %d bins, %d with |mean|>12 mas\n', f, size(R,1), sum(out));
  fprintf('   %-10s %3s %3s | %6s %6s | %6s %6s | %5s %5s | %6s %6s\n','date','n','nts','meanX','meanY','medX','medY','sdX','sdY','fldX','fldY');
  for k=find(out)'
    fprintf('   %-10s %3d %3d | %+6.1f %+6.1f | %+6.1f %+6.1f | %5.1f %5.1f | %+6.1f %+6.1f\n', datestr(R(k,1)-1721058.5,'yyyy-mm-dd'), R(k,2:3), R(k,4:5), R(k,6:7), R(k,8:9), R(k,10:11));
  end
  fprintf('   all bins: median nights/bin %d, median within-bin sd %.1f/%.1f | outliers: nights %d, sd %.1f/%.1f\n', ...
     median(R(:,3)), median(R(:,8)), median(R(:,9)), median(R(out,3)), median(R(out,8)), median(R(out,9)));
  fprintf('   corr(target bin mean, field median in same bin): X %+.2f  Y %+.2f\n', corr(R(:,4),R(:,10)), corr(R(:,5),R(:,11)));
  % how many stars share the outlier bins? per bin, fraction of 16-18 mag stars with |bin mean| > 3x their own bin rms
  mid=find(I>16&I<19); nb=size(R,1); frac=nan(nb,1);
  for q=1:nb
    s=ok&abs(JD-R(q,1))<SidMonth/2; z=nan(numel(mid),1);
    for j=1:numel(mid), w=s&W(:,mid(j))>0; if sum(w)>=MinN, z(j)=abs(mean(400*Rx(w,mid(j))))/(std(400*Rx(w,mid(j)))/sqrt(sum(w))); end; end
    frac(q)=mean(z>4,'omitnan');
  end
  fprintf('   fraction of 16-19 mag stars with |bin mean|>4 s.e. in X: all bins median %.2f, in target-outlier bins %s\n', median(frac), mat2str(round(frac(out)',2)));
end
