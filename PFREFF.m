function [varargout] = PFREFF(H)
%% Global Spectral feature
[N,k,M] = size(H); X = [];
for i = 1:M
    X = [X,H(:,:,i)];
end
X = X./sum(X,2); [EL,Ev,~] = svd(X/sqrt(diag(sum(X))));
[~,ind] = sort(diag(Ev),'descend'); X = EL(:,ind(1:k)); X = X';
%%
F = mean(H,3);
Alpha = repmat(1/(M+1),[1,1,M+1]); Q = zeros(N,k,M); q = zeros(1,1,M);
Itermax = 300; iter = 1;  err(iter,1) = inf; J(iter,1) = inf;
tol = 1e-4; T  = repmat(eye(k),1,1,M); w = q; D = zeros(size(Alpha)); 
%%
while (err(iter,1)>tol && iter<Itermax)
    O = X*F/diag(sum(F));
    Dfcm = pdist2(X',O').^2;
    D(:,:,M+1) = trace(F'*Dfcm);
    for v = 1:M
        [A,~,B] = svd(F'*H(:,:,v));
        T(:,:,v) = B*A';
        Q(:,:,v) = H(:,:,v)*T(:,:,v);
        D(:,:,v) = log(norm(F - H(:,:,v)*T(:,:,v) ,'fro').^2+1);
        q(:,:,v) = 1/(norm(F - H(:,:,v)*T(:,:,v) ,'fro').^2+1);
    end
    d = reshape(D,1,M+1); m = 2;
    alpha = bsxfun(@rdivide,d.^(-1/(m-1)),sum(d.^(-1/(m-1)),2));
    Alpha = reshape(alpha,[1,1,M+1]);
    for v = 1:M
        w(:,:,v) = Alpha(:,:,v).^m*q(:,:,v);
    end
    for i = 1:N
        f = MY_NSS_Com([reshape(Q(i,:,:),k,M),reshape(Dfcm(i,:),k,1)],w,Alpha(:,:,M+1).^m);
        F(i,:) = f';
    end
    %%
    J(iter+1,1) = d(1:M)*reshape(Alpha(1:M).^m,M,1) + Alpha(:,:,M+1).^m*D(:,:,M+1);
    err(iter+1,1) = abs(J(iter+1,1) - J(iter,1))/J(iter+1,1);
    iter = iter + 1;    
end
J = J(2:end,1);
err = err(3:end,1);
[~,label] = max(F,[],2);
varargout{1} = F;
varargout{2} = label;
varargout{3} = T;
varargout{4} = Alpha;
varargout{5} = J; varargout{6} = err; 
function [x] = MY_NSS_Com(V,w,a)
%% solve 
%  V = [v(1),v(2),...,v(p)];
%  min  1/2 sum_(p)=1^M w(p) * || x - v(p)||^2
%  s.t. x>=0, 1'x=1
%  Transform it into 
%  f(lambda_m) = 1/n*sum_{i=1}^n max(lambada_m - u_j,0) - lambda_m = 0;
%  v = sum(w(p)*v(p));
%  u = v - ones(n,n)/n*v + sum(w(p))/n*ones(n,1);
%  x_j = max(u - lambda_m, 0)/sum(w(p));
%  if umin > 0, lambda_m = 0; x = u/sum(w(p));
%  else : Newton Method.
%%
[n,M] = size(V); M = M-1;
w = reshape(w,1,M); 
v = sum(V(:,1:M).*repmat(w,n,1),2)-a*V(:,M+1);
u = v - mean(v) + sum(w)/n;
umin = min(u);
if umin >= 0
    x = u/sum(w);
else
    f = 1;
    iter = 1;
    lambda_m = 0;
    while abs(f) > 10^-10
        p = lambda_m - u;
        k = p>0;
        g = sum(k)/n - 1;
        f = sum(p(k))/n - lambda_m;
        lambda_m = lambda_m - f/g;
        iter = iter + 1;
        if iter>100
            break;
        end
    end
    x = max(-p,0)/sum(w);
end