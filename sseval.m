function sse = sseval(x,tdata,ydata)
    a = x(1);
    b = x(2);
    
    sse = sum((ydata - a + a*exp(-tdata/b)).^2);
end
