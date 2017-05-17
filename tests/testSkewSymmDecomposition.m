%
absErr = 1e-12;

%% 1x1
A = 6;
[skew, symm] = SkewSymmDecomposition(A);

assert(skew == 0);
assert(symm == A);

%% Symm
A = [1 2; 2 1];
[skew, symm] = SkewSymmDecomposition(A);

assert(sum(abs(skew(:))) == 0);
assert(sum(abs(symm(:) - A(:))) == 0);

%% Skew
A = [0 2; -2 0];
[skew, symm] = SkewSymmDecomposition(A);

assert(sum(abs(symm(:))) == 0);
assert(sum(abs(skew(:) - A(:))) == 0);

%% Mixed
skewPart = [0 1 2; -1 0 3; -2 -3 0];
symmPart = [1 2 3; 2 4 5; 3 5 6];
A = skewPart + symmPart;
[skew, symm] = SkewSymmDecomposition(A);

assert(sum(abs(symm(:) - symmPart(:))) == 0);
assert(sum(abs(skew(:) - skewPart(:))) == 0);