function dxdt = bivirus_hoi(t,x)
%BIVIRUS Differential equation of bivirus dynamics
%   The global variables are the parameters of the bivirus model. Also, dfx
%   gives the Jacobian of the bivirus system. The vector x = [x^1', x^2']'
%   is a column vector where x^i is the n-vector of infection probability
%   for all nodes, for virus i


global D1 D2 A1 A2 B1 B2 b1 b2 n

x1 = x(1:n); x2 = x(n+1:end);   %The two vectors of virus infection prob


%This part is specifically for calculating the HOI terms, in a way that
%makes the final equation more compact.
x1_hoi = zeros(n,1);
x2_hoi = zeros(n,1);

for i = 1:n
    x1_hoi(i) = x1'*B1(:,:,i)*x1;
    x2_hoi(i) = x2'*B2(:,:,i)*x2;
end

%The networked bivirus equations
dx1dt = (-D1 + (eye(n) - diag(x1) - diag(x2))*A1)*x1 + b1*(eye(n) - diag(x1) - diag(x2))*x1_hoi;
dx2dt = (-D2 + (eye(n) - diag(x1) - diag(x2))*A2)*x2+ b2*(eye(n) - diag(x1) - diag(x2))*x2_hoi;

dxdt = [dx1dt', dx2dt']';
end

