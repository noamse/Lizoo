% Where does the 9 mas monthly scatter come from? rms of bin means against bin
% length, compared with the white-noise expectation sigma_pt/sqrt(n).
addpath('/home/ocs/matlab/Lizoo'); SidMonth=27.321661;
C=load('/home/ocs/KMTdata/Results/v8/cmp.mat');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys','Info');
  IF=V.IFsys; Ie=C.Sel.(['Ie_' f]); I=V.Info.SrcData.I_ogle(:); I(I>=99)=NaN;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD)|isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  JD=IF.JD; [G,~]=ml.util.splitEpochGroups(JD);
  stars=[Ie; find(I>16&I<17&(1:IF.Nsrc)'~=Ie)];
  fprintf('\n===== %s =====\n',f);
  for k=1:2
    if k==1, list=Ie; nm='TARGET'; else, list=stars(2:end); nm=sprintf('median of %d stars 16<I<17',numel(list)); end
    out=[];
    for i=list'
      ok=W(:,i)>0; r=400*Rx(ok,i); t=JD(ok); g=G(ok); n=sum(ok); sp=std(r);
      lab={floor(t-0.5), floor((t-min(t))/7), floor((t-min(t))/SidMonth), g};
      row=[n sp];
      for L=1:4
        u=unique(lab{L}); m=nan(numel(u),1); nn=m; se=m;
        for q=1:numel(u), s=lab{L}==u(q); nn(q)=sum(s); if nn(q)<3, continue; end; m(q)=mean(r(s)); se(q)=std(r(s))/sqrt(nn(q)); end
        v=isfinite(m); m=m(v); nn=nn(v); se=se(v);
        if L==4, A=[ones(sum(v),1) t(1)*0+arrayfun(@(q) mean(t(lab{L}==q)),u(v))]; m=m-A*(A\m); end   % seasons: about the fitted line
        row=[row numel(m) median(nn) rms(m) median(sp./sqrt(nn)) sqrt(sum((m./se).^2)/numel(m))];
      end
      out(end+1,:)=row;
    end
    o=median(out,1);
    fprintf('%s: %d valid epochs, per-epoch sigma_pt %.1f mas (X)\n', nm, o(1), o(2));
    fprintf('   %-8s %6s %8s | %10s %11s %7s | %s\n','bin','Nbins','n/bin','rms(means)','sig/sqrt(n)','ratio','sqrt(chi2/dof)');
    nms={'night','week','month','season*'};
    for L=1:4, c=2+(L-1)*5; fprintf('   %-8s %6d %8d | %10.2f %11.2f %7.1f | %.1f\n', nms{L}, o(c+1), o(c+2), o(c+3), o(c+4), o(c+3)/o(c+4), o(c+5)); end
  end
  % monthly bins of the target: unweighted rms vs weighted, and by population
  ok=W(:,Ie)>0; r=400*Rx(ok,Ie); t=JD(ok); lab=floor((t-min(t))/SidMonth); u=unique(lab);
  m=nan(numel(u),1); nn=m; se=m;
  for q=1:numel(u), s=lab==u(q); nn(q)=sum(s); if nn(q)<3, continue; end; m(q)=mean(r(s)); se(q)=std(r(s))/sqrt(nn(q)); end
  v=isfinite(m); m=m(v); nn=nn(v); se=se(v);
  fprintf('TARGET monthly bins, X: unweighted rms of means %.1f | weighted rms %.1f | rms for n>=50 %.1f (%d bins) | n<50 %.1f (%d bins) | median s.e. %.2f | sqrt(chi2/dof) %.2f\n', ...
     rms(m), sqrt(sum(m.^2./se.^2)/sum(1./se.^2)), rms(m(nn>=50)), sum(nn>=50), rms(m(nn<50)), sum(nn<50), median(se), sqrt(mean((m./se).^2)));
end
