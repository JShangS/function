function [ tfr ] = JS_RXC( x ,q)
%JS_RXC 此处显示有关此函数的摘要
%   此处显示详细说明
%%计算x的对称瞬时自相关函数
%%根据tfrwv函数改变
[xrow,xcol] = size(x);
t=1:xrow; N=xrow ; 
[trow,tcol] = size(t);
if (nargin == 1)
    q = 0;
end
tfr= zeros (N,tcol);  
for icol=1:tcol
%     disp(['第',num2str(icol),'次循环']);
    ti= t(icol)+q/2;
    taumax=min([ti-1,xrow-ti,round(N/2)-1]) ;
    tau=-taumax:taumax;
    indices= rem(N+tau,N)+1;
    tfr(indices,icol) = x(ti+tau,1) .* conj(x(ti-tau,xcol));
    tau=round(N/2 - q/2)  ;
    if (ti<=xrow-tau)&(ti>=tau+1)
        tfr(tau+1,icol) = 0.5 * (x(ti+tau,1) * conj(x(ti-tau,xcol))  + ...
            x(ti-tau,1) * conj(x(ti+tau,xcol))) ;
    end;
    %  if trace, disprog(icol,tcol,10); end;
end; 
% tfr = circshift(tfr,[0,a]);
% tfr= fft(tfr); 
% if (xcol==1), tfr=real(tfr); end ;

% if (nargout==0),
%  tfrqview(tfr,x,t,'tfrwv');
% elseif (nargout==3),
%  f=(0.5*(0:N-1)/N)';
% end;

end

