function h=hypergeom2F1(a,b,c,z,tol)
%HYPERGEOM2F1  Gaussian or ordinary hypergeometric function for ABS(Z) < 1
%   H = HYPERGEOM2F1(A,B,C,Z) returns the hypergeometric function 2F1 for scalar
%   parameters A, B, C, and complex inputs Z. A finite number of terms from the
%   infinite summation representation of the function are used until the
%   absolute tolerance EPS is met.
%
%   H = HYPERGEOM2F1(...,TOL) specifies the absolute tolerance used to terminate
%   the summation of terms.
%   
%   Note:
%       Unless, C = A or C = B, if C is an integer <= 0, NaN is returned.
%       Additionally, the simple method of computation used can be very
%       inaccurate when ABS(Z) is close to 1 for some parameter combinations.
%
%   See also HYPERGEOM.

%   Andrew D. Horchler, horchler @ gmail . com, Created 5-12-13
%   Revision: 1.0, 5-19-14


% Check four required inputs
if isempty(a) || ~isscalar(a) || ~isfloat(a) || ~isreal(a) || ~isfinite(a)
    error('hypergeom2F1:AInvalid',...
          'A must be a non-empty finite real floating-point scalar.');
end
if isempty(b) || ~isscalar(b) || ~isfloat(b) || ~isreal(b) || ~isfinite(b)
    error('hypergeom2F1:BInvalid',...
          'B must be a non-empty finite real floating-point scalar.');
end
if isempty(c) || ~isscalar(c) || ~isfloat(c) || ~isreal(c) || ~isfinite(c)
    error('hypergeom2F1:CInvalid',...
          'C must be a non-empty finite real floating-point scalar.');
end
if ~isfloat(z) || ~all(isfinite(z))
    error('hypergeom2F1:ZInvalid','Z must be a finite floating-point array.');
end
if any(abs(z)>=1)
    error('hypergeom2F1:ZUndefined',...
         ['The standard Gaussian hypergeometric function is only defined ' ...
          'for ABS(Z) < 1']);
end

% Calculate Gaussian hypergeometric function via summation of terms
if isempty(z)
	h = z;
else
    % Set relative tolerance and maximum number of iterations
    dtype = superiorfloat(a,b,c,z);
    if nargin < 5
        tol = eps(dtype);
    end
    itermax = 2^15;
    
    h = zeros(size(z),dtype);
    for j = 1:numel(z)
        Z = z(j);
        
        if a == 1
            if b == c
                yi = Z;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*Z;
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            else
                yi = b*Z/c;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*(b+i)*Z/(c+i);
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            end
        elseif b == 1
            if a == c
                yi = Z;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*Z;
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            else
                yi = a*Z/c;
                y = 1+yi;
                for i = 1:itermax
                    yi = yi*(a+i)*Z/(c+i);
                    y = y+yi;
                    if abs(yi) < tol
                        break;
                    end
                end
            end
        elseif a == c
            yi = b*Z;
            y = 1+yi;
            for i = 1:itermax
                yi = yi*(b+i)*Z/(i+1);
                y = y+yi;
                if abs(yi) < tol
                    break;
                end
            end
        elseif b == c
            yi = a*Z;
            y = 1+yi;
            for i = 1:itermax
                yi = yi*(a+i)*Z/(i+1);
                y = y+yi;
                if abs(yi) < tol
                    break;
                end
            end
        else 
            yi = a*b*Z/c;
            y = 1+yi;
            for i = 1:itermax
                yi = yi*(a+i)*(b+i)*Z/((i+1)*(c+i));
                y = y+yi;
                if abs(yi) < tol
                    break;
                end
            end
        end
        
        h(j) = y;
    end
end