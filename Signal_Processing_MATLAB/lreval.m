function sse = lreval(x,tdata,ydata)
    k = x(1);
    sse = sum((ydata - k*tdata).^2);
end