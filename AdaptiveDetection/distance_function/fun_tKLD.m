function [ distance ] = fun_SKLD( A,B )
%%¡¶Matrix CFAR detectors based on symmetrized Kullback¨CLeibler and
%%%total Kullback¨CLeibler divergences¡·
%%%symmetrized  Kullback¨CLeibler divergences (14)Ê½
[N,~] = size(A);
distance = 0.5 * trac(inv(B) * A + inv(A) * B - 2 * eye(N));
end