function [ logA ] = fun_Logm( A )
[UA, LA] = svd(A);
logA = UA * diag(log(diag(LA))) * UA';
end