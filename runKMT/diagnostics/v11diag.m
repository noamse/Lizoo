% What does trimming 14 nights at each end of every season remove, and what airmass cut would do the same?
addpath('/home/ocs/matlab/Lizoo');
for F={'BLG41','BLG01'}
  f=F{1}; V=load(sprintf('/home/ocs/KMTdata/Results/v8/IFfinal_%s.mat',f),'IFsys'); IF=V.IFsys;
  JD=IF.JD; secz=IF.getTimeSeriesField(1,{'secz'}); [G,GI]=ml.util.splitEpochGroups(JD);
  fprintf('\n%s: %d epochs, %d seasons; secz quantiles 50/75/90/95: %s\n', f, numel(JD), GI.Ngroup, mat2str(quantile(secz,[.5 .75 .9 .95]),3));
  fprintf('   %-6s %5s %5s %6s | %-24s | %-24s | %s\n','season','Nep','days','nights','first 14 nights: Nep, secz','last 14 nights: Nep, secz','middle secz med/90%');
  keep14=true(size(JD)); keepZ=true(size(JD));
  tot=0;
  for b=1:GI.Ngroup
    s=find(G==b); t=JD(s); t0=min(t); t1=max(t); nights=numel(unique(floor(t-0.5)));
    e=s(t<t0+14); l=s(t>t1-14); m=s(t>=t0+14&t<=t1-14);
    keep14(e)=false; keep14(l)=false; tot=tot+numel(e)+numel(l);
    D=datetime(t0,'convertfrom','juliandate');
    fprintf('   %-6s %5d %5.0f %6d | %5d  %.2f/%.2f (med/90)   | %5d  %.2f/%.2f (med/90)   | %.2f / %.2f\n', datestr(D,'yyyy'), numel(s), t1-t0, nights, ...
      numel(e), median(secz(e)), quantile(secz(e),.9), numel(l), median(secz(l)), quantile(secz(l),.9), median(secz(m)), quantile(secz(m),.9));
  end
  fprintf('   14-night trim removes %d epochs (%.1f%%)\n', tot, 100*tot/numel(JD));
  for zc=[1.5 1.4 1.3 1.25 1.2]
    r=secz>zc; fprintf('   secz>%.2f removes %d epochs (%.1f%%); overlap with the 14-night trim: %.0f%% of the trimmed epochs have secz>%.2f, %.0f%% of secz>%.2f epochs are in the trim\n', zc, sum(r), 100*mean(r), 100*mean(r(~keep14)), zc, 100*mean(~keep14(r)), zc);
  end
  % day-of-year of the trimmed epochs vs airmass: is airmass the actual variable?
  doy=day(datetime(JD,'convertfrom','juliandate'),'dayofyear');
  fprintf('   trimmed epochs: day-of-year ranges %d-%d (start) / %d-%d (end)\n', min(doy(~keep14&doy<180)), max(doy(~keep14&doy<180)), min(doy(~keep14&doy>=180)), max(doy(~keep14&doy>=180)));
end
