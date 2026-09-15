% Generalise: each star's residual is regressed on ALL per-epoch observing
% conditions; the within-season drift of those conditions then predicts the
% star's within-season slope.  slope_pred(i,b) = sum_k beta_k(i) * drift_k(b)
addpath('/home/ocs/matlab/Lizoo');
for F={'BLG41','BLG01'}
  f=F{1};
  V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f));
  IF=V.IFsys;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  fw  = median(IF.Data.fwhm,2,'omitnan');
  dpsf= median(IF.Data.DeltaPSFXY,2,'omitnan');
  secz= IF.Data.secz(:,1); pa=IF.Data.pa(:,1);
  chi = median(IF.Data.PSF_CHI2DOF,2,'omitnan');
  C = [fw, secz, secz.*sin(pa), secz.*cos(pa), dpsf, chi];
  nm= {'fwhm','secz','secz.sinpa','secz.cospa','dPSFXY','chi2dof'};
  good=all(isfinite(C),2);
  C=C-median(C(good,:),1,'omitnan');
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup;
  % within-season drift of each condition
  D=nan(nb,size(C,2));
  for b=1:nb
    s=G==b & good; if sum(s)<50, continue; end
    t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t)); A=[ones(sum(s),1) t];
    cc=A\C(s,:); D(b,:)=cc(2,:);
  end
  fprintf('\n===== %s =====\n', f);
  for ax=1:2
    R = (ax==1)*400.*Rx + (ax==2)*400.*Ry;
    pred=[]; obs=[];
    for i=1:IF.Nsrc
      ok=W(:,i)>0 & good; if sum(ok)<300, continue; end
      A=[ones(sum(ok),1) C(ok,:)];
      beta=A\R(ok,i); beta=beta(2:end);
      for b=1:nb
        s=W(:,i)>0 & G==b; if sum(s)<50 || any(~isfinite(D(b,:))), continue; end
        t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t));
        a=[ones(sum(s),1) t]\R(s,i);
        obs(end+1,1)=a(2); pred(end+1,1)=D(b,:)*beta;
      end
    end
    k=polyfit(pred,obs,1); res=obs-polyval(k,pred);
    fprintf('  axis %d: corr %+.3f over %d star-seasons; variance explained %.0f%%\n', ...
        ax, corr(pred,obs), numel(obs), 100*(1-var(res)/var(obs)));
    fprintf('           slope rms %.2f -> %.2f mas after removing the conditions model\n', std(obs), std(res));
  end
  clear IF V
end
