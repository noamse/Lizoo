function [keep, note] = v11flag(JD, Mode)
% Epoch selection for v11. Mode 'trim': drop the first 14 and last 14 nights of
% every season (the last 14 only where the season ends after day-of-year 250,
% so the truncated 2026 season keeps its final, mid-season, epochs).
% Mode 'seczNN': drop epochs with sec z > NN/100 (target's airmass at CTIO).
JD=JD(:);
switch Mode(1:4)
  case 'trim'
    [G,GI]=ml.util.splitEpochGroups(JD); keep=true(size(JD)); nS=0; nE=0;
    for b=1:GI.Ngroup
      s=G==b; t0=min(JD(s)); t1=max(JD(s));
      e=s&JD<t0+14; keep(e)=false; nS=nS+sum(e);
      if day(datetime(t1,'convertfrom','juliandate'),'dayofyear')>250
        l=s&JD>t1-14; keep(l)=false; nE=nE+sum(l);
      end
    end
    note=sprintf('trim: %d season starts, %d season ends dropped', nS, nE);
  case 'secz'
    zc=str2double(Mode(5:end))/100;
    Lat=-(30+10/60+1.84/3600); Lon=-(70+48/60+14.39/3600);
    RA=celestial.coo.convertdms('17:52:38.09','gH','d'); Dec=celestial.coo.convertdms('-31:47:36.1','gD','d');
    [~,Alt]=celestial.coo.radec2azalt(JD, RA, Dec, 'GeoCoo',[Lon Lat],'InUnits','deg','OutUnits','deg');
    secz=1./sind(Alt); keep=secz<=zc;
    note=sprintf('secz<=%.2f: %d dropped', zc, sum(~keep));
  otherwise, error('v11flag:Mode','unknown mode %s',Mode);
end
end
