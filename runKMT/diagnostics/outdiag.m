% why are some sidereal-month bins strong outliers? n, airmass, season edge, seeing
addpath('/home/ocs/matlab/Lizoo'); SidMonth=27.321661;
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v8/Tie_%s.mat',f));
  IF=V.IFsys; Ie=C.Sel.(['Ie_' f]); Amap=T.Amap;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  secz=IF.getTimeSeriesField(1,{'secz'}); fw=median(IF.Data.fwhm,2,'omitnan');
  ok=W(:,Ie)>0; JD=IF.JD(ok);
  S=1000*(Amap\[Rx(ok,Ie)'; Ry(ok,Ie)']); ra=S(1,:)'-mean(S(1,:)); de=S(2,:)'-mean(S(2,:));
  sz=secz(ok); fh=fw(ok);
  lab=floor((JD-min(JD))/SidMonth); u=unique(lab);
  R=[];
  for q=1:numel(u)
    s=lab==u(q); n=sum(s); if n<3, continue; end
    doy=day(datetime(mean(JD(s)),'convertfrom','juliandate'),'dayofyear');
    R(end+1,:)=[mean(JD(s)) n mean(ra(s)) std(ra(s))/sqrt(n) mean(de(s)) std(de(s))/sqrt(n) median(sz(s)) median(fh(s)) doy];
  end
  dev=sqrt(R(:,3).^2+R(:,5).^2); out=dev>15; 
  fprintf('\n%s: %d bins, %d with |resid|>15 mas\n', f, size(R,1), sum(out));
  fprintf('   median n per bin: all %d | outliers %d\n', median(R(:,2)), median(R(out,2)));
  fprintf('   median secz     : all %.2f | outliers %.2f\n', median(R(:,7)), median(R(out,7)));
  fprintf('   median fwhm     : all %.2f | outliers %.2f\n', median(R(:,8)), median(R(out,8)));
  fprintf('   day-of-year     : all %d..%d | outliers: %s\n', min(R(:,9)), max(R(:,9)), mat2str(R(out,9)'));
  fprintf('   %-10s %4s %8s %8s %6s %6s %4s\n','date','n','dRA','dDec','secz','fwhm','doy');
  for k=find(out)'
    fprintf('   %-10s %4d %+8.1f %+8.1f %6.2f %6.2f %4d\n', datestr(R(k,1)-1721058.5,'yyyy-mm-dd'), R(k,2), R(k,3), R(k,5), R(k,7), R(k,8), R(k,9));
  end
  % correlation of |deviation| with n
  fprintf('   Spearman |dev| vs n: %+.2f ; vs secz %+.2f ; vs fwhm %+.2f\n', ...
     corr(dev,R(:,2),'type','Spearman'), corr(dev,R(:,7),'type','Spearman'), corr(dev,R(:,8),'type','Spearman'));
  fprintf('   bins with n<10: %d of %d ; their median |dev| %.1f vs %.1f for n>=10\n', sum(R(:,2)<10), size(R,1), median(dev(R(:,2)<10)), median(dev(R(:,2)>=10)));
end
