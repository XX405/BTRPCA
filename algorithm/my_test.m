function E= my_test(W,e)

    % E = zeros(size(t));
    E=W;
    % n = size(W,2);
    % t = pagenorm(W);
    for i = 42
         t = pagenorm(E);
        if abs(t(i)) < e
            E(i) = -((abs(t(i)) - e)^2 / 2) + e^2 / 2;
        else
            E(i) = e^2 / 2;
        end
    end
end
