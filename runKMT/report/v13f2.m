% chi^2 on sidereal-month bins, for every star; plus RA/Dec tracks.
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
SidMonth=27.321661; Fs={'BLG41','BLG01'};
C=load('/home/ocs/KMTdata/Results/v13/cmp.mat'); Sel=C.Sel;
Chi=struct();
for K=1:2
  f=Fs{K};
  V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v13/Tie_%s.mat',f));
  IF=V.IFsys; NF=V.Info; Amap=T.Amap;
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN; Ie=IF.findClosestSource([150 150]);
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  lab0=floor((IF.JD-min(IF.JD))/SidMonth);
  X2=nan(IF.Nsrc,1); DoF=nan(IF.Nsrc,1);
  for i=1:IF.Nsrc
    ok=W(:,i)>0; if sum(ok)<100, continue; end
    lab=lab0(ok); rx=400*Rx(ok,i); ry=400*Ry(ok,i);
    u=unique(lab); s2=0; nb=0;
    for q=1:numel(u)
      s=lab==u(q); n=sum(s); if n<5, continue; end
      mx=mean(rx(s)); my=mean(ry(s));
      ex=std(rx(s))/sqrt(n); ey=std(ry(s))/sqrt(n);
      if ex<=0||ey<=0, continue; end
      s2=s2+(mx/ex)^2+(my/ey)^2; nb=nb+1;
    end
    if nb>=6, X2(i)=s2; DoF(i)=2*nb-4; end
  end
  Chi.(f)=struct('X2',X2,'DoF',DoF,'Ie',Ie,'I',I);
  fprintf('%s: chi2 on %d stars | target chi2 %.0f, DoF %d, chi2/DoF %.1f | field median chi2/DoF %.1f\n', ...
      f, sum(isfinite(X2)), X2(Ie), DoF(Ie), X2(Ie)/DoF(Ie), median(X2./DoF,'omitnan'));
  clear IF V
end
save('/home/ocs/KMTdata/Results/v13/chi2.mat','Chi','-v7.3');

Fh=figure('Visible','off','Position',[40 40 1400 560]);
TL=tiledlayout(Fh,1,2,'TileSpacing','compact','Padding','compact');
for K=1:2
  f=Fs{K}; c=Chi.(f);
  nexttile(TL,K); hold on; box on; grid on; set(gca,'FontSize',12);
  v=c.X2(isfinite(c.X2)); ed=logspace(log10(max(10,min(v))),log10(max(v)*1.1),40);
  histogram(v,ed,'FaceColor',[0.35 0.60 0.85],'EdgeColor',[0.1 0.3 0.6]);
  set(gca,'XScale','log');
  yl=ylim; plot([c.X2(c.Ie) c.X2(c.Ie)],yl,'r-','LineWidth',3);
  text(c.X2(c.Ie),yl(2)*0.92,sprintf('  target  \\chi^2 = %.0f (DoF %d)',c.X2(c.Ie),c.DoF(c.Ie)), ...
      'Color','r','FontSize',12,'FontWeight','bold');
  xlabel('\chi^2 on sidereal-month bins (not divided by DoF)','FontSize',13);
  ylabel('number of stars','FontSize',13);
  title(sprintf('%s   —  %d stars, median DoF %d', f, sum(isfinite(c.X2)), ...
      round(median(c.DoF,'omitnan'))),'FontWeight','bold','FontSize',14);
end
title(TL,'v13: \chi^2 of the monthly binned residuals, all calibration stars','FontWeight','bold','FontSize',15);
exportgraphics(Fh,[Doc 'report_v13_chi2.png'],'Resolution',150); close(Fh);
fprintf('V13F2 DONE\n');
