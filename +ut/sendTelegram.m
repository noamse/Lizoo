function varargout = sendTelegram(varargin)
% sendTelegram  (public stub) No-op for release. Prints message to console.
fprintf('[sendTelegram stub] ');
for i=1:nargin
    if ischar(varargin{i}) || isstring(varargin{i})
        fprintf('%s ', string(varargin{i}));
    else
        fprintf('<arg%d>', i);
    end
end
fprintf('\n');
if nargout>0
    varargout = cell(1,nargout);
    [varargout{:}] = deal([]);
end
end
