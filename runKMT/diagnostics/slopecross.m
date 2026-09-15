% Do INDIVIDUAL stars show the same season-slope pattern in both fields?
addpath('/home/ocs/matlab/Lizoo'); Doc='/home/ocs/matlab/Lizoo/doc/';
OGf=[Doc(1:end-4) 'OGLEdata/OB260058/OB160058.mat'];
Q=load('/home/ocs/KMTdata/Results/slopeQ.mat'); Q=Q.Q;
% match the two fields' stars through the OGLE catalogue
id=cell(1,2); Fs={'BLG41','BLG01'};
for k=1:2
  f=Fs{k}; Cat=ml.util.ogleCompanionCat(OGf,'Field',f);
  X=Q.(f).X; Y=Q.(f).Y;
  D=sqrt((X-Cat(:,1)').^2+(Y-Cat(:,2)').^2);
  [d,j]=min(D,[],2); j(d>=1.5)=nan; id{k}=j;
end
[c,i1,i2]=intersect(id{1}(isfinite(id{1})),id{2}(isfinite(id{2})));
a=find(isfinite(id{1})); b=find(isfinite(id{2}));
a=a(i1); b=b(i2);
fprintf('%d stars matched between the fields through OGLE\n\n', numel(a));
for ax=1:2
  r=nan(numel(a),1);
  for z=1:numel(a)
    p1=squeeze(Q.BLG41.SL(a(z),:,ax))'; p2=squeeze(Q.BLG01.SL(b(z),:,ax))';
    g=isfinite(p1)&isfinite(p2);
    if sum(g)>=6, r(z)=corr(p1(g),p2(g)); end
  end
  gg=isfinite(r);
  fprintf('axis %d: per-star season-pattern correlation between fields\n', ax);
  fprintf('   median %+.3f over %d stars ; fraction positive %.0f%%\n', ...
      median(r(gg)), sum(gg), 100*mean(r(gg)>0));
  I=Q.BLG41.I(a);
  for e=[14 16; 16 17; 17 18; 18 19]'
    s=gg&I>=e(1)&I<e(2);
    if sum(s)<6, continue; end
    fprintf('     %4.1f-%4.1f n=%3d  median corr %+.3f\n', e(1),e(2),sum(s),median(r(s)));
  end
  % the target
  it=find(a==Q.BLG41.Ie);
  if ~isempty(it), fprintf('     TARGET corr %+.3f (percentile %.0f%%)\n', r(it), 100*mean(r(gg)<r(it))); end
end
