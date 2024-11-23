%% 
% input:
%   - fiberPath  fiber对应的.bgData路径
%   - WGPath     波导.bgData文件夹路径
%
% Return:
%   请参看函数<modefield_overlap_bgData>help内容
%   
%   给"debug"赋值1后，默认在命令行中显示计算结果并绘图
%

fiberPath = '.\test\光纤.bgData';
WGPath = '.\test\V500KHz-50mW-160μm-13mms.bgData';

warning('off','curvefit:prepareFittingData:removingNaNAndInf');
[Int2,Int2_Sim,gof] = modefield_overlap_bgData(WGPath,fiberPath,"Llimit",0.1,"debug",1,"Colormap",'OSI_rainbow');
