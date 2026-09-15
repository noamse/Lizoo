% season-scale clipping candidates on the v8 solution: how many stars go under each definition
addpath('/home/ocs/matlab/Lizoo');
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info');
  IF=V.IFsys; NF=V.Info; Ie=C.Sel.(['Ie_' f]); I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  i31=C.Sel.(f)(strcmp(C.Sel.tags,'m16_d31'));
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD)|isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  [G,GI]=ml.util.splitEpochGroups(IF.JD); nb=GI.Ngroup; N=IF.Nsrc;
  Rs=nan(N,nb); Off=nan(N,nb); Ch=nan(N,nb);
  for i=1:N, for b=1:nb
    s=W(:,i)>0&G==b; n=sum(s); if n<30, continue; end
    r=400*[Rx(s,i) Ry(s,i)];
    Rs(i,b)=sqrt(mean(r(:).^2));                      % rms about zero, both axes
    mo=mean(r,1); se=std(r,0,1)/sqrt(n); Off(i,b)=norm(mo); Ch(i,b)=sum((mo./se).^2);
  end, end
  Rdec=sqrt(mean(Rs.^2,2,'omitnan'));                 % decade rms (for reference)
  A_max=max(Rs,[],2);                                 % worst-season rms
  B_chi=sum(Ch,2,'omitnan')./(2*sum(isfinite(Ch),2)); % season-offset chi2/dof
  ed=14:0.5:19;
  stat={log10(A_max),'worst-season rms (log)'; log10(B_chi),'season-offset chi2/dof (log)'; log10(max(Rs./Rdec,[],2)),'worst season / decade rms (log)'};
  fprintf('\n=== %s (%d stars) ===\n',f,N); lbl={'kept','rejected'};
  for k=1:size(stat,1)
    L=stat{k,1}; ok=isfinite(I)&isfinite(L); rej=false(N,1); z=nan(N,1);
    for q=1:numel(ed)-1
      s=ok&I>=ed(q)&I<ed(q+1); if sum(s)<6, continue; end
      m=median(L(s)); sd=tools.math.stat.rstd(L(s)); z(s)=(L(s)-m)/sd; rej(s)=L(s)>m+2*sd;
    end
    rej(Ie)=false;
    fprintf('  %-36s: %2d rejected | target z %+.2f | m16_d31 z %+.2f %s\n', stat{k,2}, sum(rej), z(Ie), z(i31), lbl{rej(i31)+1});
    % also: per-season 2 sigma, flagged in >=2 seasons
  end
  % per-season clipping, count seasons flagged
  nfl=zeros(N,1);
  for b=1:nb
    L=log10(Rs(:,b)); ok=isfinite(I)&isfinite(L);
    for q=1:numel(ed)-1
      s=ok&I>=ed(q)&I<ed(q+1); if sum(s)<6, continue; end
      m=median(L(s)); sd=tools.math.stat.rstd(L(s)); nfl(s)=nfl(s)+(L(s)>m+2*sd);
    end
  end
  fprintf('  %-36s: >=1 season %d, >=2 seasons %d, >=3 seasons %d | target %d, m16_d31 %d seasons\n', 'per-season 2-sigma, count', sum(nfl>=1&(1:N)'~=Ie), sum(nfl>=2&(1:N)'~=Ie), sum(nfl>=3&(1:N)'~=Ie), nfl(Ie), nfl(i31));
  fprintf('  m16_d31 season rms [mas]: %s\n', mat2str(round(Rs(i31,:),1)));
  fprintf('  m16_d31 season offsets  : %s\n', mat2str(round(Off(i31,:),1)));
end
