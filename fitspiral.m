function s()
    xdata = [399 347 292 235 167 109  76  96 197 262 323 410 463 559 566 557 527];
    ydata = [413 434 440 435 410 367 325 184 106  70  70  80 128 278 338 381 458];
    xydata = [xdata; ydata]';
    th=linspace(0, 2.2*pi, length(xdata));
    const=[ 196 156 301 266 5.40 10 .8];
    param = lsqcurvefit(@fitspiral, const, th, xydata)
    [xyfit] = fitspiral(param,th);
    figure(1)
    plot(xdata, ydata, 'bp')
    hold on
    plot(xyfit(:,1), xyfit(:,2), '-m')
    hold off
    grid
end 
%%
   function [xy] = fitspiral(param, th)
    a = param(1);
    b = param(2);
    c = param(3);
    d = param(4);
    e = param(5);
    f = param(5);
    g = param(5);
    R=((a+(f+g).*th).^2.*(b+(f+g).*th).^2./((a+(f+g).*th).^2.*sin(th).^2+(b+(f+g).*th).^2.*cos(th).^2)).^(.5);
    x=R.*cos(th)+c;
    y=R.*sin(th)+d+(th/pi).^e;
    xy = [x; y]';
   end