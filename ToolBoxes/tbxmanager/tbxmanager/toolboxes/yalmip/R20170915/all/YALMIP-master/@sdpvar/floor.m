function varargout=floor(varargin)
%FLOOR (overloaded)

switch class(varargin{1})
    
    case 'sdpvar' % Overloaded operator for SDPVAR objects. Pass on args and save them.
        
        x = varargin{1};
        dim = size(x);
        x = reshape(x,prod(dim),1);
        y = [];
        for i = 1:prod(dim)
            y = [y;yalmip('define',mfilename,extsubsref(x,i))];
        end
        y = reshape(y,dim);
        varargout{1} = y;
        
    case 'char'
        
        t = varargin{2};
        X = varargin{3};
        
        F = ([X-1 <= t <= X]) + (integer(t));
        
        varargout{1} = F;
        varargout{2} = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
        varargout{3} = X;
    otherwise
        error('Strange type on first argument in SDPVAR/FLOOR');
end
