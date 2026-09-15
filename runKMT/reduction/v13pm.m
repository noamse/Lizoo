% Proper motions with uncertainties. The error is taken from the scatter of the
% SEASON means about the fitted line, not from the per-epoch count: the nightly
% noise is strongly correlated, so a per-epoch error would be far too small.
addpath('/home/ocs/matlab/Lizoo');
Res=struct();
for F={'BLG41','BLG01'}
  f=F{1};
  V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v13/Tie_%s.mat',f));
  IF=V.IFsys; Amap=T.Amap; Gauge=T.Gauge;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup;
  t=(IF.JD-IF.JD0)/365.25;
  SigPix=nan(IF.Nsrc,2);
  for i=1:IF.Nsrc
    ok=W(:,i)>0;
    px=400*(IF.ParS(1,i)+IF.ParS(3,i)*t)+400*Rx(:,i);
    py=400*(IF.ParS(2,i)+IF.ParS(4,i)*t)+400*Ry(:,i);
    tb=nan(nb,1); mx=nan(nb,1); my=nan(nb,1);
    for b=1:nb
      s=ok&(G==b); if sum(s)<20, continue; end
      tb(b)=mean(t(s)); mx(b)=mean(px(s)); my(b)=mean(py(s));
    end
    k=isfinite(tb)&isfinite(mx);
    if sum(k)<4, continue; end
    A=[ones(sum(k),1) tb(k)];
    cx=A\mx(k); cy=A\my(k);
    n=sum(k); Sxx=sum((tb(k)-mean(tb(k))).^2);
    sx=std(mx(k)-A*cx,1)*sqrt(n/(n-2))/sqrt(Sxx);
    sy=std(my(k)-A*cy,1)*sqrt(n/(n-2))/sqrt(Sxx);
    SigPix(i,:)=[sx sy];      % mas/yr, in pixel axes
  end
  % propagate the pixel-axis errors onto the sky through the same linear map
  % Amap maps arcsec -> pix, and sky motion = 1000*(Amap \ pixel motion).
  % Propagate the pixel-axis variances through that same linear map.
  Jac = 1000*(Amap\eye(2));
  Sig=nan(IF.Nsrc,2);
  for i=1:IF.Nsrc
     if ~isfinite(SigPix(i,1)), continue; end
     C=diag((SigPix(i,:)/400).^2);          % mas/yr on pixel axes -> pix/yr variance
     Sig(i,:)=sqrt(diag(Jac*C*Jac'))';
  end
  Ie=IF.findClosestSource([150 150]);
  Res.(f)=struct('PMabs',T.PMabs,'Sig',Sig,'Ie',Ie,'clean',T.clean);
  fprintf('\n=== %s ===\n',f);
  fprintf('  target  mu_a*cosd %+7.3f +- %5.3f   mu_d %+7.3f +- %5.3f  mas/yr\n', ...
      T.PMabs(1,Ie), Sig(Ie,1), T.PMabs(2,Ie), Sig(Ie,2));
  ok=isfinite(Sig(:,1));
  fprintf('  median PM uncertainty over all %d sources: %.3f / %.3f mas/yr\n', ...
      sum(ok), median(Sig(ok,1)), median(Sig(ok,2)));
  fprintf('  median over the %d RUWE-clean tie stars:   %.3f / %.3f mas/yr\n', ...
      sum(T.clean&ok), median(Sig(T.clean&ok,1)), median(Sig(T.clean&ok,2)));
  clear IF V
end
d=Res.BLG41.PMabs(:,Res.BLG41.Ie)'-Res.BLG01.PMabs(:,Res.BLG01.Ie)';
s=sqrt(Res.BLG41.Sig(Res.BLG41.Ie,:).^2+Res.BLG01.Sig(Res.BLG01.Ie,:).^2);
fprintf('\ntarget, difference between fields: %+.3f +- %.3f / %+.3f +- %.3f mas/yr\n', ...
    d(1), s(1), d(2), s(2));
