function z = fun(x,y)
    
    z = 0.01*ones(size(x));
    z(abs(x-1)  < 0.5) = 0.005;
    z(abs(x-1)  < 0.25) = 0.003;