function hypergeomq_test
%
clear all
m = 5;
n = 202;
x = logspace(-154,0,n);
%r = digits(16);
h1 = hypergeom([0.5 0.5],[1.5 1.5],x);
tic
for i = 1:m
    h1 = hypergeom([0.5 0.5],[1.5 1.5],x);
end
t1 = toc;
%digits(r);
disp([num2str(t1)  ' seconds, ' num2str(t1/(m*n)) ' seconds per evaluation.'])
%{
clear all
m = 5;
n = 400;
x = logspace(-154,0,n);
r = digits(16);
h2 = double(mupadmex('symobj::map',['matrix('...
    sprintf('%d,1,[%.15g',n,x(1))...
    sprintf(',%.15g',x(2:end)) '])'],...
    'symobj::hypergeom','matrix(1,2,[0.5,0.5])',...
    'matrix(1,2,[1.5,1.5])'));
tic
for i = 1:m
    h2 = double(mupadmex('symobj::map',['matrix('...
    sprintf('%d,1,[%.15g',n,x(1))...
    sprintf(',%.15g',x(2:end)) '])'],...
    'symobj::hypergeom','matrix(1,2,[0.5,0.5])',...
    'matrix(1,2,[1.5,1.5])'));
end
t2 = toc;
digits(r);
disp([num2str(t2)  ' seconds, ' num2str(t2/(m*n)) ' seconds per evaluation.'])
%
clear all
m = 5;
n = 400;
x = logspace(-154,0,n);
r = digits(16);
h3 = double(mupadmex('symobj::map',[...
    sprintf('[%.15g',x(1))...
    sprintf(',%.15g',x(2:end)) ']'],...
    'symobj::hypergeom','[0.5,0.5]',...
    '[1.5,1.5]'))';
tic
for i = 1:m
    h3 = double(mupadmex('symobj::map',[...
    sprintf('[%.15g',x(1))...
    sprintf(',%.15g',x(2:end)) ']'],...
    'symobj::hypergeom','[0.5,0.5]',...
    '[1.5,1.5]'))';
end
t3 = toc;
digits(r);
disp([num2str(t3)  ' seconds, ' num2str(t3/(m*n)) ' seconds per evaluation.'])
%}
clear all
m = 5;
n = 202;
x = logspace(-154,0,n);
h4 = hypergeomq([0.5 0.5],[1.5 1.5],x);
tic
for i = 1:m
    h4 = hypergeomq([0.5 0.5],[1.5 1.5],x);
end
t4 = toc;
disp([num2str(t4)  ' seconds, ' num2str(t4/(m*n)) ' seconds per evaluation.'])