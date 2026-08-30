function [ParS] = initiateParS(IF,Args)

arguments
    IF;
    Args.Plx = true;
end

if IF.Plx
    ParS = zeros(5,IF.Nsrc);
else
    ParS = zeros(4,IF.Nsrc);
%     if IF.FakePlx 
%         ParS = zeros(9,IF.Nsrc);
%     else
%         ParS = zeros(5,IF.Nsrc);
%     end
% elseif IF.FakePlx 
%     ParS = zeros(8,IF.Nsrc);
% else
%     ParS = zeros(4,IF.Nsrc);
end

if ~isempty(IF.InitialXYGuess)
    Xguess = IF.InitialXYGuess(:,1);
    Yguess = IF.InitialXYGuess(:,2);

else
    Xguess = IF.medianFieldSource({'X'});
    Yguess = IF.medianFieldSource({'Y'});
end
ParS([1,2],:)= [Xguess';Yguess'];

% Parameters that are held rather than fitted start at the value they are held
% at, since the solver will never move them off it.
if ~isempty(IF.ParSFixed)
    if ~isequal(size(IF.ParSFixed),size(ParS))
        error('IterFit:initiateParS:BadParSFixed', ...
              'ParSFixed is %dx%d but ParS is %dx%d', ...
              size(IF.ParSFixed,1), size(IF.ParSFixed,2), size(ParS,1), size(ParS,2));
    end
    Held = isfinite(IF.ParSFixed);
    ParS(Held) = IF.ParSFixed(Held);
end

end