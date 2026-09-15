% Target, before and after proper-motion detrending, both fields: 2x2.
% Columns are the fields, rows are PM retained / PM removed. Each panel carries
% both sky axes so the whole target behaviour is on one page.
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
SidMonth=27.321661; Fs={'BLG41','BLG01'};
C=load('/home/ocs/KMTdata/Results/v13/cmp.mat');
Fh=figure('Visible','off','Position',[40 40 1500 900]);
TL=tiledlayout(Fh,2,2,'TileSpacing','compact','Padding','compact');
for K=1:2
  f=Fs{K};
  V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v13/Tie_%s.mat',f));
  IF=V.IFsys; NF=V.Info; Amap=T.Amap; Ie=C.Sel.(['Ie_' f]);
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN;
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  ok=W(:,Ie)>0; JD=IF.JD(ok); t=(JD-IF.JD0)/365.25;
  pxA=400*(IF.ParS(1,Ie)+IF.ParS(3,Ie)*t)+400*Rx(ok,Ie);   % PM retained
  pyA=400*(IF.ParS(2,Ie)+IF.ParS(4,Ie)*t)+400*Ry(ok,Ie);
  pxB=400*Rx(ok,Ie); pyB=400*Ry(ok,Ie);                     % PM removed
  lab=floor((JD-min(JD))/SidMonth); u=unique(lab);
  for R=1:2
    if R==1, px=pxA; py=pyA; else, px=pxB; py=pyB; end
    S=1000*(Amap\[px'/400; py'/400]);
    ra=S(1,:)'-mean(S(1,:)); de=S(2,:)'-mean(S(2,:));
    TT=nan(numel(u),1); MR=TT;ER=TT;MD=TT;ED=TT;
    for q=1:numel(u)
      s=lab==u(q); n=sum(s); if n<3, continue; end
      TT(q)=mean(JD(s)); MR(q)=mean(ra(s)); ER(q)=std(ra(s))/sqrt(n);
      MD(q)=mean(de(s)); ED(q)=std(de(s))/sqrt(n);
    end
    g=isfinite(TT); D=datetime(TT(g),'convertfrom','juliandate');
    nexttile(TL,(R-1)*2+K); hold on; box on; grid on; set(gca,'FontSize',12);
    h1=errorbar(D,MR(g),ER(g),'o','Color',[0.10 0.35 0.65],'MarkerFaceColor',[0.35 0.60 0.85], ...
        'MarkerSize',4,'CapSize',0,'LineWidth',0.8);
    h2=errorbar(D,MD(g),ED(g),'s','Color',[0.65 0.10 0.10],'MarkerFaceColor',[0.95 0.45 0.35], ...
        'MarkerSize',4,'CapSize',0,'LineWidth',0.8);
    yline(0,'k-','LineWidth',1);
    ylim([-55 55]); xlabel('date','FontSize',13); ylabel('offset [mas]','FontSize',13);
    if R==1
      title(sprintf('%s   proper motion RETAINED',f),'FontWeight','bold','FontSize',14);
    else
      title(sprintf('%s   proper motion REMOVED   (rms %.1f / %.1f mas)', f, ...
          std(MR(g)), std(MD(g))),'FontWeight','bold','FontSize',14);
    end
    if K==1&&R==1
      legend([h1 h2],{'\Delta\alpha cos\delta','\Delta\delta'},'Location','southwest','FontSize',11);
    end
  end
  clear IF V
end
title(TL,'The target, before and after proper-motion detrending, sidereal-month bins', ...
   'FontWeight','bold','FontSize',15);
exportgraphics(Fh,[Doc 'report_v13_target_beforeafter.png'],'Resolution',150); close(Fh);
fprintf('BA DONE\n');
