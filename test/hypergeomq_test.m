function hypergeomq_test
%HYPERGEOM_TEST  

%   Andrew D. Horchler, adh9 @ case . edu, Created 7-23-12
%   Revision: 1.1, 4-8-16


clear all;  %#ok<CLALL>
m = 5;
n = 202;
x = logspace(-154,0,n);
h1 = hypergeom([0.5 0.5],[1.5 1.5],x);      %#ok<NASGU>
tic
for i = 1:m
    h1 = hypergeom([0.5 0.5],[1.5 1.5],x);  %#ok<NASGU>
end
t1 = toc;
disp(['hypergeom: ' num2str(t1)  ' seconds, ' num2str(t1/(m*n)) ' seconds per evaluation.']);

clear all;  %#ok<CLALL>
m = 5;
n = 202;
x = logspace(-154,0,n);
h2 = hypergeomq([0.5 0.5],[1.5 1.5],x);     %#ok<NASGU>
tic
for i = 1:m
    h2 = hypergeomq([0.5 0.5],[1.5 1.5],x); %#ok<NASGU>
end
t2 = toc;
disp(['hypergeomq: ' num2str(t2)  ' seconds, ' num2str(t2/(m*n)) ' seconds per evaluation.']);