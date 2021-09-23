function varargout = nanall(varargin)


for i = 1:nargout
    varargout{i} = nan(varargin{:});
end