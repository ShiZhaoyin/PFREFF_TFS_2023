function [varargout] = MyFCM(X,k,m)
[~,N] = size(X); J = inf; err = inf; iter = 1; Itermax = 100;
U = rand(N,k); thr =1e-4;
label = ceil(k*rand(1,N));
U = full(sparse(1:N,label,1,N,k,N));
while(err(iter)>thr && iter<Itermax)
    V = X*U.^m./sum(U.^m,1);
    G = pdist2(X',V');
    U = bsxfun(@rdivide,G.^(-2/(m-1)),sum(G.^(-2/(m-1)),2));
    J(iter+1) = sum(sum(U.^m.*G.^2,1),2);
    err(iter+1) = abs(J(iter)-J(iter+1))/J(iter+1);
    iter = iter+1;
end
J = J(2:end);
err = err(2:end);
[~,lable] = max(U,[],2);
varargout{1} = U;
varargout{2} = lable;
varargout{3} = V;
varargout{4} = J;
varargout{5} = err;
end

