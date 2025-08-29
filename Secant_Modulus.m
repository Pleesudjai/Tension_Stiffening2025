function ks = Secant_Modulus(X,Y,K, x)

for n=1:length(x)
    if x(n) <= X(2)
        ks(n) = K(1);
    elseif x(n) > X(end)
        ks(n) = Y(end)/x(n);
    else
        for i=2:length(X)-1
            if x(n)>X(i) && x(n)<=X(i+1)
                ks(n) = (Y(i) + K(i)*(x(n)-X(i)))/x(n);
            end
        end
    end
end