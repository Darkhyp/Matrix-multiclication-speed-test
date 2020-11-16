% matrix multiplivation Tests
% N = 1000; % 0.354141 seconds.
N = 2000; % 0.324714s (old0.591440) (gcc+openblas 20 times 0.377805) seconds.
% N = 5000; % 6.170928 seconds.
% N = 10000; % 45.787484 seconds.

DIR = 'D:\_install\_Programming\parallel_computing\_matrix_product_tests\';
isLoad = true;
% isLoad = false;
version -blas

tic
if (isLoad)
    fprintf('\nLoading matrix A...')
    A       = load([DIR,'A.dat']);
    fprintf('\nLoading matrix B...')
    B       = load([DIR,'B.dat']);
    fprintf('\nLoading matrix Ctmp...')
    Ctmp	= load([DIR,'C.dat']);
    N  = size(A,1);
else
    A = rand(N);
    B = rand(N);
%{
    A = rand(n)+1i*rand(n);
    B = rand(n)+1i*rand(n);
%}
end
toc

Ncount = 20;
fprintf('\nMatrix product in matlab: c = a*b for the size (%ix%i)',N,N)
tic;
for i=1:Ncount
    C = A*B;
end
fprintf('\nCalculation time is %gs', toc/Ncount)

if (isLoad)
    tic
    fprintf('\n norm(C-Ctmp)=%g',norm(C-Ctmp))
    toc
end


%{
gA = gpuArray(A);
gB = gpuArray(B);

gputimeit(@()gA*gB)
%}