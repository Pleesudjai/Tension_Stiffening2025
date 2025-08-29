function s = Linear_Equation(bc,P,h,psi,Af,Ef,k,g,nodeCRK)

n = length(k);          % total number of nodes
a = zeros(n,1);         % coefficient matrix
b = zeros(n,1);         % coefficient matrix
c = zeros(n,1);         % coefficient matrix
d = zeros(n,1);         % right hand side vector
x = zeros(n,1);         % unknown slip

% FD equation for finite length WITHOUT spring
for i=2:n-1
    a(i) =  1/2*Af*Ef(i-1)+1/2*Af*Ef(i);
    b(i) = -1/2*Af*Ef(i+1)-1/2*Af*Ef(i-1)-Af*Ef(i)-psi*k(i)*h^2;
    c(i) =  1/2*Af*Ef(i)+1/2*Af*Ef(i+1);
end

% % impose spring equation if they exists at node i
% % FD equation for finite length WITH spring
% for i=1:length(g)
%     j = nodeCRK(i);
%     % consider spring contribution if they are not at the end nodes
%     if j>2 || j<n-2
%         b(j) = b(j)-g(i)*h;
%     end
% end

% impose boundary conditions on both end
if isempty(g) == 1
    switch bc
        case 1  % bc=1 F(1)=0, F(L)=P
            b(1) = -1/2*psi*k(1)*h^2-1/2*Af*Ef(1)-1/2*Af*Ef(2);
            c(1) =  1/2*Af*Ef(2)+1/2*Af*Ef(1);
            a(n) =  1/2*Af*Ef(n-1)+1/2*Af*Ef(n);
            b(n) = -1/2*Af*Ef(n-1)-1/2*Af*Ef(n)-1/2*psi*k(n)*h^2;
            d(n) = -P*h;
        case 2  % bc=2 F(1)=P, F(L)=0
            b(1) =  1/2*Af*Ef(1)+1/2*Af*Ef(2)+1/2*psi*k(1)*h^2;
            c(1) = -1/2*Af*Ef(1)-1/2*Af*Ef(2);
            a(n) = -1/2*Af*Ef(n-1)-1/2*Af*Ef(n);
            b(n) =  1/2*Af*Ef(n-1)+1/2*Af*Ef(n)+1/2*psi*k(n)*h^2;
            d(1) = -P*h;       
        case 3  % bc=3 F(1)=P, F(L)=P
            b(1) =  1/2*Af*Ef(1)+1/2*Af*Ef(2)+1/2*psi*k(1)*h^2;
            c(1) = -1/2*Af*Ef(1)-1/2*Af*Ef(2);
            a(n) =  1/2*Af*Ef(n-1)+1/2*Af*Ef(n);
            b(n) = -1/2*Af*Ef(n-1)-1/2*Af*Ef(n)-1/2*psi*k(n)*h^2;
            d(1) = -P*h;
            d(n) = -P*h;
    end
else
    switch bc
        case 1  % bc=1 F(1)=0, F(L)=P
            b(1) = -1/2*psi*k(1)*h^2-1/2*Af*Ef(1)-1/2*Af*Ef(2);
            c(1) =  1/2*Af*Ef(2)+1/2*Af*Ef(1);
            a(n) =  1/2*Af*Ef(n-1)+1/2*Af*Ef(n);
            b(n) = -1/2*Af*Ef(n-1)-1/2*Af*Ef(n)-1/2*psi*k(n)*h^2-g(n)*h;
            d(n) = -P*h;
        case 2  % bc=2 F(1)=P, F(L)=0
            b(1) =  1/2*Af*Ef(1)+1/2*Af*Ef(2)+1/2*psi*k(1)*h^2+g(1)*h;
            c(1) = -1/2*Af*Ef(1)-1/2*Af*Ef(2);
            a(n) = -1/2*Af*Ef(n-1)-1/2*Af*Ef(n);
            b(n) =  1/2*Af*Ef(n-1)+1/2*Af*Ef(n)+1/2*psi*k(n)*h^2+g(n)*h;
            d(1) = -P*h;       
        case 3  % bc=3 F(1)=P, F(L)=P
            b(1) =  1/2*Af*Ef(1)+1/2*Af*Ef(2)+1/2*psi*k(1)*h^2+g(1)*h;
            c(1) = -1/2*Af*Ef(1)-1/2*Af*Ef(2);
            a(n) =  1/2*Af*Ef(n-1)+1/2*Af*Ef(n);
            b(n) = -1/2*Af*Ef(n-1)-1/2*Af*Ef(n)-1/2*psi*k(n)*h^2;
            d(1) = -P*h;
            d(n) = -P*h;
    end
end


% solve linear system equations using Thomas Algorithm
% http://www.cfd-online.com/Wiki/Tridiagonal_matrix_algorithm_-_TDMA_(Thomas_algorithm)
% check conditions
% cond = abs(b)-(abs(a)+abs(c))

for i = 2:n 
    m    = a(i)/b(i-1);
    b(i) = b(i) - m*c(i-1);
    d(i) = d(i) - m*d(i-1);
end

x(n) = d(n)/b(n);
for i = n-1:-1:1
    x(i) =  (d(i)-c(i)*x(i+1))/b(i);
end

s = x';