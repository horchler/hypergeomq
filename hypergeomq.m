function h=hypergeomq(n,d,z)
%HYPERGEOMQ  Fast evaluation of the generalized hypergeometric function.
%   H = HYPERGEOMQ(N,D,Z) evaluates the generalized hypergeometric function for
%   the vector parameters N and D at the values in the array Z. N, D, and Z can
%   be any numeric or logical datatype and may be complex. The output, H, will
%   have the same dimensions as Z and will be double-precision unless Z is
%   single-precision, in which case, H will be as well.
%
%   HYPERGEOMQ uses the same low-level function as HYPERGEOM, but implements
%   several optimizations to achieve a performance boost of approximately an
%   order of magnitude. Results from the two functions should agree to the
%   precision of the inputs. Additionaly, HYPERGEOMQ can avoid some errors due
%   to singularities that can occur when HYPERGEOM is evaluated with numeric Z
%   values.
%
%   Examples:
%       % An identity for the hypergeomentric function 1F0(a;;z)
%       z = -3:0.1:3; hz = hypergeomq(1,[],z); lz = 1./(1-z);
%       figure; plot(z,hz,'b',z,lz,'r.'); grid on; axis([-3 3 -10 10]);
%       xlabel('z'); ylabel('1F0(a;;z) == 1/(1-z)');
%
%       % Relationship between Confluent and Gaussian Hypergeometric Functions
%       z = logspace(1,2,20); h1 = hypergeomq(1,3,z); a2 = [2e2 5e2 2e3];
%       for i = 1:3; h2(i,:) = hypergeomq([1 a2(i)],3,z/a2(i)); end;
%       figure; loglog(z,abs(h1),'k',z,abs(h2)); grid on;
%       xlabel('z'); ylabel('|1F1(a;b;z)|, |2F1(a,b;c;z)|');
%       legend('_1F_1(a_1;b;z)','_2F_1(a_1,200;b;z/200)',...
%           '_2F_1(a_1,500;b;z/500)','_2F_1(a_1,2000;b;z/2000)',...
%           'Location','NorthWest');
%
%   Class support for ouput H:
%       float: double, single
%
%   See also: HYPERGEOM, SYM/HYPERGEOM.

%   Andrew D. Horchler, horchler @ gmail . com, Created 7-22-12
%   Revision: 1.1, 4-21-16


% Check zeroes
if ~isvector(n) && ~isempty(n)
    error('hypergeomq:NonVectorN','The first input must be a vector.');
end
if ~isnumeric(n) && ~islogical(n)
    error('hypergeomq:NonNumericN',...
          'The first input argument must be a numeric or logical vector.');
end
if ~all(isfinite(n))
    error('hypergeomq:NonFiniteN',...
          'The first input must be a vector of finite values.');
end

% Check poles
if ~isvector(d) && ~isempty(d)
    error('hypergeomq:NonVectorD','The second input must be a vector.');
end
if ~isnumeric(d) && ~islogical(d)
    error('hypergeomq:NonNumericD',...
          'The second input argument must be a numeric or logical vector.');
end
if ~all(isfinite(d))
    error('hypergeomq:NonFiniteD',...
          'The second input must be a vector of finite values.');
end

% Check Z values
if ~isnumeric(z) && ~islogical(z)
    error('hypergeomq:NonNumericD',...
          'The third input argument must be a numeric or logical array.');
end
if ~all(isfinite(z))
    error('hypergeomq:NonFiniteZ',...
          'The third input must be an array of finite values.');
end

if isa(z,'single')
    dtype = 'single';
else
    dtype = 'double';
end
if isempty(z)
    h = zeros(size(z),dtype);
else
    % Convert zeroes to a string
    N = float2str(n);
    
    % Convert poles to a string
    D = float2str(d);
    
    % Try getting an analytic expression as a function of z
    mupadmexExists = (exist('mupadmex','file') == 3);
    
    if mupadmexExists
        % Low-level function with pre-formatted arguments, no validation
        sy = mupadmex('symobj::hypergeom','z',N,D);
        sc = mupadmex('symobj::char',charcmd(sy),0);
        if ~isempty(strncmp(sc,'"_symans_',9))
            sc = char(sy);
        end
    else
        % Pass string arguments and convert output to string
        sc = char(hypergeom(N,D,'z'));
    end
    
    % Evaluate vectorized analytic expression if not function of hypergeometrics
    if isempty(strfind(sc,'hypergeom('))
        if ~isscalar(z)
            % Vectorize analytic expression
            sc = strrep(sc,'*','.*');
            sc = strrep(sc,'/','./');
            sc = strrep(sc,'^','.^');
        end
        
        % Analytical expression can return errors, e.g., hypergeomq(1,0.5,-1)
        try
            h = eval(sc);
            return;
        catch	%#ok<CTCH>  
        end
    end
    
    % Convert Z values to a string
    Z = float2str(z);

    % Evaluate numerically
    try
        if mupadmexExists
            % Low-level function with pre-formatted arguments and no validation
            h = mupadmex('map',Z,'symobj::hypergeom',N,D);
        else
            h = hypergeom(N,D,Z);
        end
    catch ME
        if strcmp(ME.identifier,'symbolic:mupadmex:CommandError')...
                && (~isempty(strfind(ME.message,'Singularity'))...
                || ~isempty(strfind(ME.message,'Division by zero')))
            error('hypergeomq:mupadmexSingularityError',...
                  'One or more singularities exist at or near the Z values.');
        else
            rethrow(ME);
        end
    end
    
    % Convert output to floating point - much faster than calling sym/double
	h = sym2float(h,dtype);
    
    % Reshape matrix and multi-dimensional array inputs
    if ~isvector(z)
        h = reshape(h,size(z));
    elseif size(z,1) > 1
        h = h.';
    end
end



function s=float2str(z)
%FLOAT2STR  Fast conversion of floating-point arrays to symbolic strings.
%   S = FLOAT2STR(Z) returns a compact string representation of the
%   floating-point array Z that can be used by low-level symbolic functions. The
%   input Z may include complex and non-finite values. The representation used
%   only describes the elements of Z, not the dimensions.
%
%   See also SPRINTF, SYM.

%   Andrew D. Horchler, adh9 @ case . edu, Created 11-11-13
%   Revision: 1.0, 5-19-14


if ~(isnumeric(z) || islogical(z))
    error('hypergeomq:float2str:NonNumeric',...
          'Input must be a numeric or logical array.')
end

z = z(:);
if isreal(z) || all(imag(z) == 0)
    if isinteger(z) || islogical(z)
        printStr = '%d.0,';
    else
        isInt = (floor(z) == z);
        if all(isInt)
            printStr = '%.1f,';
        elseif any(isInt)
            if isa(z,'double')
                printStrs = {'%.17g,','%0.1f,'};
            else
                printStrs = {'%.9g,','%.1f,'};
            end
            printStr = reshape(char(printStrs(isInt+1)).',1,[]);
        else
            if isa(z,'double')
                printStr = '%.17g,';
            else
                printStr = '%.9g,';
            end
        end
    end
    
    if isscalar(z)
        s = sprintf(printStr(1:end-1),z);
    else
        s = ['[' sprintf(printStr,z)];
        s(end) = ']';
    end
else
    isRealInt = (floor(real(z)) == real(z));
    isImagInt = (floor(imag(z)) == imag(z));
    if all(isRealInt) && all(isImagInt)
        printStr = '%.1f+%.1f*i,';
    elseif any(isRealInt) || any(isImagInt)
        if isa(z,'double')
            printStrs = {'%.17g+%.17g*i,','%0.1f+%.17g*i,',...
                         '%.17g+%0.1f*i,','%0.1f+%0.1f*i,'};
        else
            printStrs = {'%.9g+%.9g*i,','%.1f+%.9g*i,',...
                         '%.9g+%.1f*i,','%.1f+%.1f*i,'};
        end
        printStr = reshape(char(printStrs(isRealInt+2*isImagInt+1)).',1,[]);
    else
        if isa(z,'double')
            printStr = '%.17g+%.17g*i,';
        else
            printStr = '%.9g+%.9g*i,';
        end
    end
    
    if isscalar(z)
        s = sprintf(printStr(1:end-1),real(z),imag(z));
    else
        s = ['[' sprintf(printStr,real(z),imag(z))];
        s(end) = ']';
    end
end


function z=sym2float(s,dtype)
%SYM2FLOAT  Fast conversion of symbolic matrices and arrays to floating-point.
%   Z = SYM2FLOAT(S) converts the array of symbolic values, S, and returns them
%   as a column vector of doubles. S may include complex and non-finite values.
%   In general, this conversion is faster than evaluating Z = DOUBLE(S) or
%   similar for symbolic matrices and 1-D and 2-D symbolic arrays.
%
%   Z = SYM2FLOAT(S,DATATYPE) specifies an optional floating-point datatype,
%   'single' or 'double' (the default).
%
%   See also DOUBLE, SINGLE, CAST, STR2DOUBLE, SYM.

%   Andrew D. Horchler, adh9 @ case . edu, Created 11-11-13
%   Revision: 1.0, 5-19-14


if ~isa(s,'sym')
    error('hypergeomq:sym2float:InvalidSymbol',...
          'Input data must be of class ''sym''.');
end

if nargin > 1
    if strcmp(dtype,'double')
        scanStr = '%f64';
    elseif strcmp(dtype,'single')
        scanStr = '%f32';
    else
        error('hypergeomq:sym2float:InvalidDatatype',...
              'Optional datatype must be ''single'' or ''double'' (default).');
    end
else
    dtype = 'double';
    scanStr = '%f64';
end

% Fast conversion from sym to char
try
    if exist('mupadmex','file') == 3 && ndims(s) < 3	%#ok<ISMAT>
        % Convert from symbolic to string
        s = mupadmex('symobj::char',charcmd(s),0);
        if ~isempty(strncmp(s,'"_symans_',9))
            s = char(s);
        end
    else
        s = char(s);
    end
catch ME
    if strcmp(ME.identifier,'symbolic:mupadmex:CommandError')...
            && (~isempty(strfind(ME.message,'Singularity'))...
            || ~isempty(strfind(ME.message,'Division by zero')))
        error('hypergeomq:sym2float:mupadmexSingularityError',...
              'One or more singularities exist.');
    else
        rethrow(ME);
    end
end

% Replace special constants and format for conversion
s = strrep(s(2:end-1),'RD_INF','Inf');
s = strrep(s,'RD_NINF','-Inf');
s = strrep(s,'RD_NAN','NaN');
s = strrep(s,' ','');

if s(1) == 'm'
    % Count rows and trim matrix([[ ... ] ... [ ... ]])
    m = length(strfind(s,'],['))+1;
    s = s(find(s == '[',1):end-2);
    s = strrep(s,'[','');
    s = strrep(s,']','');

    % Format complex values as A + Bi
    if ~isempty(find(s == 'i',1))
        s = strrep(s,'*i','i');
        s = regexprep(s,'([^,]*i)([^,]*)','$2+$1');
        s = strrep(s,'+-','-');
    end

    s = textscan(s,scanStr,'Delimiter',',');
    z = [s{:}];
    d = [m numel(z)/m];
elseif s(1) == 'a'
    % Find dimensions and trim array( ... )
    idx = find(s == '=',1)+1;
    d = textscan(s(7:idx),'%*d64%*[..]%d64,');
    d = [d{:}].';
    s = s(idx:end-1);

    % Format complex values as A + Bi
    if ~isempty(find(s == 'i',1))
        s = strrep(s,'*i','i');
        s = regexprep(s,'([^,]*i)([^,]*)','$2+$1');
        s = strrep(s,'+-','-');
    end

    if any(d == 0)
        z = [];
    else
        s = textscan(s,...
            [scanStr ',(%*d64' repmat(',%*d64',1,length(d)-1) ')=']);
        z = [s{:}];
    end
else
    % Scalar case
    z = cast(str2double(s),dtype);
    d = 1;
end

% Reshape matrix inputs
if d(1) > 1
    if d(2) == 1
        z = reshape(z,d);
    else
        z = reshape(z,d).';
    end
else
    z = z.';
end