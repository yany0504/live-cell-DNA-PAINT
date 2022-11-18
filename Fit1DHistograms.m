function [x,resnorm,residual,exitflag]=Fit1DHistograms(Z,xdata,x0);

lb = [0,min(xdata),0,0,min(xdata),0];
ub = [realmax('double'),0,5,realmax('double'),0,5];
opts=optimoptions('lsqcurvefit');
opts.MaxFunctionEvaluations = 100000;
opts.MaxIterations = 100000;
[x,resnorm,residual,exitflag] = lsqcurvefit(@Two1DGaussFunction,x0,xdata,Z,lb,ub,opts);