% RA vs time, Dec vs time, and RA vs Dec colour-coded by time, for the target
% and two well-measured comparison stars.
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
SidMonth=27.321661; Fs={'BLG41','BLG01'};
MinN=10; MaxSecz=1.3;   % bins with fewer epochs, and epochs above this airmass, are dropped
C=load('/home/ocs/KMTdata/Results/v13/cmp.mat'); Sel=C.Sel;
Want={'m16_d24','m16_d30'};
for K=1:2
  f=Fs{K};
  V=load(sprintf('/home/ocs/KMTdata/Results/v13/IFfinal_%s.mat',f));
  T=load(sprintf('/home/ocs/KMTdata/Results/v13/Tie_%s.mat',f));
  IF=V.IFsys; NF=V.Info; Amap=T.Amap;
  I=NF.SrcData.I_ogle(:); I(I>=99)=NaN; Ie=C.Sel.(['Ie_' f]);
  [Rx,Ry]=IF.calculateResiduals; W=IF.calculateWes;
  Bad=isoutlier(Rx,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD) ...
    | isoutlier(Ry,'movmedian',30,'ThresholdFactor',1.5,'SamplePoints',IF.JD);
  W(Bad)=0; W(~isfinite(Rx)|~isfinite(Ry))=0;
  secz=IF.getTimeSeriesField(1,{'secz'});
  Lst={Ie,'target'};
  for q=1:numel(Sel.tags)
     if any(strcmp(Sel.tags{q},Want)) && isfinite(Sel.(f)(q))
        Lst(end+1,:)={Sel.(f)(q),Sel.tags{q}};
     end
  end
  for L=1:size(Lst,1)
    Ii=Lst{L,1}; tag=Lst{L,2};
    ok=W(:,Ii)>0 & secz<MaxSecz; if sum(ok)<100, continue; end
    JD=IF.JD(ok); t=(JD-IF.JD0)/365.25;
    % residual only: position AND proper motion removed
    px=400*Rx(ok,Ii);
    py=400*Ry(ok,Ii);
    S=1000*(Amap\[px'/400; py'/400]);          % mas on the sky
    ra=S(1,:)'-mean(S(1,:)); de=S(2,:)'-mean(S(2,:));
    lab=floor((JD-min(JD))/SidMonth); u=unique(lab);
    TT=nan(numel(u),1); MR=TT; ER=TT; MD=TT; ED=TT;
    for q=1:numel(u)
      s=lab==u(q); n=sum(s); if n<MinN, continue; end
      TT(q)=mean(JD(s)); MR(q)=mean(ra(s)); ER(q)=std(ra(s))/sqrt(n);
      MD(q)=mean(de(s)); ED(q)=std(de(s))/sqrt(n);
    end
    g=isfinite(TT); TT=TT(g);MR=MR(g);ER=ER(g);MD=MD(g);ED=ED(g);
    % empirical factor: sqrt(chi2/dof) of the bin means about zero, per axis,
    % so the bars carry the real bin-to-bin scatter rather than the white-noise s.e.
    Fr=max(1,sqrt(sum((MR./ER).^2)/(numel(MR)-1))); Fd=max(1,sqrt(sum((MD./ED).^2)/(numel(MD)-1)));
    ER=ER*Fr; ED=ED*Fd;
    D=datetime(TT,'convertfrom','juliandate'); yr=year(D)+(day(D,'dayofyear')-1)/365.25;
    Fh=figure('Visible','off','Position',[40 40 1500 520]);
    TL=tiledlayout(Fh,1,3,'TileSpacing','compact','Padding','compact');
    nexttile(TL,1); hold on; box on; grid on; set(gca,'FontSize',12);
    errorbar(D,MR,ER,'o','Color',[0.10 0.35 0.65],'MarkerFaceColor',[0.35 0.60 0.85], ...
        'MarkerSize',4,'CapSize',0,'LineWidth',0.8);
    xlabel('date','FontSize',13); ylabel('\Delta\alpha cos\delta [mas]','FontSize',13);
    yline(0,'k-','LineWidth',1);
    title(sprintf('\\Delta RA against time   (bars = s.e. \\times %.1f)',Fr),'FontWeight','bold','FontSize',13);
    nexttile(TL,2); hold on; box on; grid on; set(gca,'FontSize',12);
    errorbar(D,MD,ED,'o','Color',[0.65 0.10 0.10],'MarkerFaceColor',[0.95 0.45 0.35], ...
        'MarkerSize',4,'CapSize',0,'LineWidth',0.8);
    xlabel('date','FontSize',13); ylabel('\Delta\delta [mas]','FontSize',13);
    yline(0,'k-','LineWidth',1);
    title(sprintf('\\Delta Dec against time   (bars = s.e. \\times %.1f)',Fd),'FontWeight','bold','FontSize',13);
    nexttile(TL,3); hold on; box on; grid on; axis equal; set(gca,'FontSize',12);
    plot(MR,MD,'-','Color',[0.55 0.55 0.55],'LineWidth',0.5);   % consecutive months joined
    scatter(MR,MD,44,yr,'filled','MarkerEdgeColor','k');
    cb=colorbar; cb.Label.String='year'; colormap(gca,parula);
    xlabel('\Delta\alpha cos\delta [mas]','FontSize',13); ylabel('\Delta\delta [mas]','FontSize',13);
    plot(0,0,'k+','MarkerSize',14,'LineWidth',1.5);
    title('\Delta RA against \Delta Dec, colour = time','FontWeight','bold','FontSize',13);
    title(TL,sprintf(['v13  %s  %s  —  I = %.2f, position and proper motion ', ...
        'REMOVED, sidereal-month bins (\\geq%d epochs, sec z < %.1f)'], f, strrep(tag,'_','-'), I(Ii), MinN, MaxSecz),'FontWeight','bold','FontSize',15);
    exportgraphics(Fh,sprintf('%sreport_v13_radec_%s_%s.png',Doc,f,tag),'Resolution',150);
    close(Fh);
    fprintf('%s %-9s: %d bins | residual rms %.2f / %.2f mas | span %.0f / %.0f mas | factor %.2f / %.2f\n', ...
        f, tag, numel(TT), std(MR), std(MD), max(MR)-min(MR), max(MD)-min(MD), Fr, Fd);
  end
  clear IF V
end
fprintf('V13F3 DONE\n');
