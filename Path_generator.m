function P = Path_generator(xs, xf, eta, B)
    k = 1;
    eps = 0.1;
    P(1,k) = xs(1);
    P(2,k) = xs(2);
    [m,~] = size(B);
    while sqrt((P(1,k)-xf(1))^2+(P(2,k)-xf(2))^2) >= 0.1
        Fa = -eta*[P(1,k)-xf(1);P(2,k)-xf(2)];
        Fre = [0;0];
        for i = 1:m
            M = [B(i,1);B(i,2)];
            L = [B(i,3);B(i,4)];
            alpha = inv([L-M,-[0 -1;1 0]*(L-M)])*[P(1,k)-M(1);P(2,k)-M(2)];
            H = alpha(1)*(L-M) + M;
            if alpha(1)>=0 && alpha(1) <=1
                b = H;
            elseif alpha(1)<0
                b = M;
            else
                b = L;
            end
            R = sqrt((b(1)-P(1,k))^2+(b(2)-P(2,k))^2);
            if R <= B(i,5)
                FFr = B(i,6)*(1/(R^3))*((1/R)-(1/B(i,5)))*[P(1,k)-b(1);P(2,k)-b(2)];
                Fr = [FFr(1);FFr(2)];
            else
                Fr = [0;0];
            end
            Fre = Fre + Fr;
        end
        F = Fre + Fa;
        P(1,k+1) = P(1,k)+(eps/sqrt((F(1))^2+(F(2))^2))*F(1);
        P(2,k+1) = P(2,k)+(eps/sqrt((F(1))^2+(F(2))^2))*F(2);
        k = k+1;
    end
    P = P';
end

