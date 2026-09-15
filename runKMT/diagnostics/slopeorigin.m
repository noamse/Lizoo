% Where do the within-season slopes come from?
% Per star: the rms of its own within-season slopes across the seasons, then
% correlate that against every candidate cause.
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
OGf=[Doc(1:end-4) 'OGLEdata/OB260058/OB160058.mat'];
RUN=getenv('RUN'); if isempty(RUN), RUN='v8'; end
Q=struct();
for F={'BLG41','BLG01'}
  f=F{1};
  V=load(sprintf('/home/ocs/KMTdata/Results/%s/IFfinal_%s.mat',RUN,f));
  IF=V.IFsys; NF=V.Info; Ie=IF.findClosestSource([150 150]);
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup;
  secz=IF.getTimeSeriesField(1,{'secz'}); pa=IF.getTimeSeriesField(1,{'pa'});
  SL=nan(IF.Nsrc,nb,2); SE=nan(IF.Nsrc,nb);
  % within-season drift of the observing geometry, the suspected driver
  dsecz=nan(nb,1); dpa=nan(nb,1);
  for b=1:nb
    s=G==b; if sum(s)<50, continue; end
    t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t)); A=[ones(sum(s),1) t];
    c=A\secz(s); dsecz(b)=c(2);
    c=A\sin(pa(s)); dpa(b)=c(2);
  end
  for i=1:IF.Nsrc
    for b=1:nb
      s=W(:,i)>0 & G==b; n=sum(s); if n<50, continue; end
      t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t)); A=[ones(n,1) t];
      cx=A\(400*Rx(s,i)); cy=A\(400*Ry(s,i));
      SL(i,b,1)=cx(2); SL(i,b,2)=cy(2);
      rr=400*Rx(s,i)-A*cx; SE(i,b)=std(rr)/sqrt(n)/std(t);
    end
  end
  srms=[std(SL(:,:,1),0,2,'omitnan') std(SL(:,:,2),0,2,'omitnan')];
  eint=median(SE,2,'omitnan');           % expected slope error from white noise
  Cat=ml.util.ogleCompanionCat(OGf,'Field',f);
  X=IF.ParS(1,:)'; Y=IF.ParS(2,:)';
  D=sqrt((X-Cat(:,1)').^2+(Y-Cat(:,2)').^2); Ds=sort(D,2); nnd=Ds(:,2);
  C=IF.Data.C; col=median(C,1,'omitnan')';
  rad=sqrt((X-150).^2+(Y-150).^2);
  g=isfinite(srms(:,1))&isfinite(I);
  fprintf('\n===== %s (%s) =====  %d stars\n', f, RUN, sum(g));
  fprintf('  per-star season-slope rms: median %.2f / %.2f mas (X/Y)\n', ...
      median(srms(g,1)), median(srms(g,2)));
  fprintf('  expected from white noise : median %.2f mas  -> EXCESS factor %.1f\n', ...
      median(eint(g)), median(srms(g,1))./median(eint(g)));
  fprintf('  Spearman of slope rms with:\n');
  fprintf('     I_OGLE            %+.3f\n', corr(I(g),srms(g,1),'type','Spearman'));
  fprintf('     white-noise expect %+.3f\n', corr(eint(g),srms(g,1),'type','Spearman'));
  fprintf('     colour V-I        %+.3f\n', corr(col(g),srms(g,1),'type','Spearman','rows','complete'));
  fprintf('     nearest nbr dist  %+.3f\n', corr(nnd(g),srms(g,1),'type','Spearman'));
  fprintf('     radius from centre %+.3f\n', corr(rad(g),srms(g,1),'type','Spearman'));
  % excess over white noise, per star
  exc=srms(:,1)./eint;
  fprintf('  excess over white noise: median %.2f, 90th pct %.2f, target %.2f\n', ...
      median(exc(g)), quantile(exc(g),0.9), exc(Ie));
  fprintf('  TARGET slope rms %.2f / %.2f mas, percentile among |I-Itgt|<0.4: %.0f%% / %.0f%%\n', ...
      srms(Ie,1), srms(Ie,2), ...
      100*mean(srms(g&abs(I-I(Ie))<0.4,1)<srms(Ie,1)), ...
      100*mean(srms(g&abs(I-I(Ie))<0.4,2)<srms(Ie,2)));
  % does the SEASON pattern track the geometry drift?
  for ax=1:2
    mp=squeeze(median(SL(g,:,ax),1,'omitnan'))';
    k=isfinite(mp)&isfinite(dsecz);
    fprintf('  field-median season slope (axis %d) vs within-season d(secz): corr %+.3f ; vs d(sin pa) %+.3f\n', ...
        ax, corr(mp(k),dsecz(k)), corr(mp(k),dpa(k)));
  end
  Q.(f)=struct('SL',SL,'srms',srms,'I',I,'Ie',Ie,'col',col,'nnd',nnd,'eint',eint,'X',X,'Y',Y);
  clear IF V
end
% per-star reproduction between fields: the decisive test
save('/home/ocs/KMTdata/Results/slopeQ.mat','Q','-v7.3');
fprintf('\n(per-star cross-field matching done separately)\n');
