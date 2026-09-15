addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
OGf=[Doc(1:end-4) 'OGLEdata/OB260058/OB160058.mat'];
nm={'fwhm','secz','secz.sinpa','secz.cospa','dPSFXY','chi2dof'};

Keep=struct();
for F={'BLG41','BLG01'}
  f=F{1};
  V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f));
  IF=V.IFsys;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  fw=median(IF.Data.fwhm,2,'omitnan'); dpsf=median(IF.Data.DeltaPSFXY,2,'omitnan');
  secz=IF.Data.secz(:,1); pa=IF.Data.pa(:,1); chi=median(IF.Data.PSF_CHI2DOF,2,'omitnan');
  C=[fw, secz, secz.*sin(pa), secz.*cos(pa), dpsf, chi];
  good=all(isfinite(C),2); C=C-median(C(good,:),1,'omitnan');
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup;
  D=nan(nb,size(C,2));
  for b=1:nb
    s=G==b&good; if sum(s)<50, continue; end
    t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t)); cc=[ones(sum(s),1) t]\C(s,:); D(b,:)=cc(2,:);
  end
  fprintf('\n===== %s =====  variance explained, leaving ONE variable out (axis 1)\n', f);
  for drop=0:numel(nm)
    use=setdiff(1:numel(nm),drop);
    pred=[];obs=[];
    for i=1:IF.Nsrc
      ok=W(:,i)>0&good; if sum(ok)<300, continue; end
      A=[ones(sum(ok),1) C(ok,use)]; beta=A\(400*Rx(ok,i)); beta=beta(2:end);
      for b=1:nb
        s=W(:,i)>0&G==b; if sum(s)<50||any(~isfinite(D(b,use))), continue; end
        t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t));
        a=[ones(sum(s),1) t]\(400*Rx(s,i));
        obs(end+1,1)=a(2); pred(end+1,1)=D(b,use)*beta;
      end
    end
    k=polyfit(pred,obs,1); r=obs-polyval(k,pred); ve=100*(1-var(r)/var(obs));
    if drop==0, fprintf('   all six variables            : %5.1f%%\n', ve); base=ve;
    else,       fprintf('   without %-14s       : %5.1f%%   (costs %4.1f)\n', nm{drop}, ve, base-ve); end
  end
  % residual slopes after removing the model, for the cross-field test
  SLr=nan(IF.Nsrc,nb,2); SLo=nan(IF.Nsrc,nb,2);
  for i=1:IF.Nsrc
    ok=W(:,i)>0&good; if sum(ok)<300, continue; end
    A=[ones(sum(ok),1) C(ok,:)];
    bx=A\(400*Rx(ok,i)); by=A\(400*Ry(ok,i)); bx=bx(2:end); by=by(2:end);
    for b=1:nb
      s=W(:,i)>0&G==b; if sum(s)<50||any(~isfinite(D(b,:))), continue; end
      t=IF.JD(s); t=(t-mean(t))/(max(t)-min(t)); B=[ones(sum(s),1) t];
      a1=B\(400*Rx(s,i)); a2=B\(400*Ry(s,i));
      SLo(i,b,1)=a1(2); SLo(i,b,2)=a2(2);
      SLr(i,b,1)=a1(2)-D(b,:)*bx; SLr(i,b,2)=a2(2)-D(b,:)*by;
    end
  end
  Keep.(f)=struct('SLo',SLo,'SLr',SLr,'X',IF.ParS(1,:)','Y',IF.ParS(2,:)');
  clear IF V
end
save('/home/ocs/KMTdata/Results/slopemodel.mat','Keep','-v7.3');
% cross-field reproduction, before and after removing the conditions model
% match the two fields directly: the cut-outs differ by a fixed small offset
sh=[-1.113 0.241];
P1=[Keep.BLG41.X Keep.BLG41.Y]; P2=[Keep.BLG01.X Keep.BLG01.Y]+sh;
Dm=sqrt((P1(:,1)-P2(:,1)').^2+(P1(:,2)-P2(:,2)').^2);
[d12,j12]=min(Dm,[],2); [~,j21]=min(Dm,[],1);
a=[]; b=[];
for i=1:numel(d12)
   if d12(i)<1.0 && j21(j12(i))==i, a(end+1)=i; b(end+1)=j12(i); end
end
fprintf('\n===== cross-field reproduction over %d matched stars =====\n', numel(a));
for ax=1:2
  for which={'SLo','before removing the model'; 'SLr','after removing the model'}'
    r=nan(numel(a),1);
    for z=1:numel(a)
      p1=squeeze(Keep.BLG41.(which{1})(a(z),:,ax))'; p2=squeeze(Keep.BLG01.(which{1})(b(z),:,ax))';
      g=isfinite(p1)&isfinite(p2); if sum(g)>=6, r(z)=corr(p1(g),p2(g)); end
    end
    fprintf('  axis %d, %-28s : median corr %+.3f, %.0f%% positive\n', ...
        ax, which{2}, median(r,'omitnan'), 100*mean(r(isfinite(r))>0));
  end
end
