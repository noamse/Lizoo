function Bs = calculateBs(IF)

%[Ax,Ay] = IF.generateSourceDesignMat;
Ax = IF.AsX;
Ay = IF.AsY;

[Rx,Ry]     = IF.calculateResiduals;
Rx(isnan(Rx))= 0;
Ry(isnan(Ry))= 0;
Wes = calculateWes(IF);

Bmat = Ax'*(Rx.*Wes) + Ay'*(Ry.*Wes);

% Clear the right hand side of every held parameter, to go with the row and
% column calculateNss clears in the normal matrix. Together they make the
% increment of a held parameter exactly zero.
if ~isempty(IF.ParSFixed)
    Bmat(isfinite(IF.ParSFixed)) = 0;
end

Bs = reshape(Bmat ,[],1);





end