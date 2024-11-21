function [fitresult,gof,gauss2,WB] = gauss2fit(x,y,z)
%gauss2fit(X, Y, Z) 二维高斯拟合
%  二维高斯拟合
%
%  https://blog.csdn.net/u012366767/article/details/90743083
%
%  f(x,y) = A*exp(-(x-X0)^2/(2*sigmaX2)-(y-Y0)^2/(2*sigmaY2))
%
%  ln(f) = lnA - (x-X0)^2/(2*sigmaX2) - (y-Y0)^2/(2*sigmaY2)
%  ln(f) = p00 + p10*x + p01*y + p20*x^2 + p02*y^2;
%  p00 = lnA - X0^2/(2*sigmaX2) - Y0^2/(2*sigmaY2)
%  p10 = X0/sigmaX2
%  p01 = Y0/sigmaY2
%  p20 = -0.5/sigmaX2
%  p02 = -0.5/sigmaY2
%
%

% 利用多项式拟合
[xData,yData,zData] = prepareSurfaceData(x,y,log(z));
ft = fittype('poly22');
opts = fitoptions('Method','LinearLeastSquares');
opts.Lower = [-Inf 0 0 -Inf 0 -Inf];
opts.Upper = [Inf Inf Inf 0 0 0];
% opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf];
% opts.Upper = [Inf Inf Inf Inf Inf Inf];
[fitresult,gof] = fit([xData,yData],zData,ft,opts);
sigmaX2 = -0.5/fitresult.p20;
sigmaY2 = -0.5/fitresult.p02;
% 将多项式拟合结果转换为二维高斯结果
ft = fittype('A*exp((-(x-X0).^2/sigmaX2-(y-Y0).^2/sigmaY2)/2)', ...
    independent=["x" "y"],dependent='img', ...
    coefficients=["A" "X0" "Y0" "sigmaX2" "sigmaY2"]);
X0 = fitresult.p10*sigmaX2;
Y0 = fitresult.p01*sigmaY2;
A = exp(fitresult.p00+X0^2/(2*sigmaX2)+Y0^2/(2*sigmaY2));
gauss2 = sfit(ft,A,X0,Y0,sigmaX2,sigmaY2);
% 将多项式拟合结果转换为全光束拟合(Whole Beam Fit Equations)结果
WBft = fittype('A*exp(-2*(((x-X0)./(wX/2)).^2+((y-Y0)./(wY/2)).^2))', ...
    independent=["x" "y"],dependent='img', ...
    coefficients=["A" "X0" "Y0" "wX" "wY"]);
wX=sqrt(sigmaX2)*4;
wY=sqrt(sigmaY2)*4;
A = exp(fitresult.p00+X0^2/(2*sigmaX2)+Y0^2/(2*sigmaY2));
WB = sfit(WBft,A,X0,Y0,wX,wY);
end