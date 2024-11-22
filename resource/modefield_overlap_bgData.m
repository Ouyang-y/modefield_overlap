function [Int2,Int2_Sim,gof] = modefield_overlap_bgData(file_WG,file_fiber,options)
% modefield_overlap_bgData  calc overlap int from 2 .bgData files.
%   [Int2,Int2_Sim,gof] = modefield_overlap_bgData(file_WG,file_fiber)
%   [Int2,Int2_Sim,gof] = modefield_overlap_bgData(file_WG,file_fiber,"Llimit",Llimit)
%   [Int2,Int2_Sim,gof] = modefield_overlap_bgData(file_WG,file_fiber,"Llimit",Llimit,"debug",debug)
%   [Int2,Int2_Sim,gof] = modefield_overlap_bgData(".\test\1.bgData",".\test\fiber.bgData","Llimit",0.1,"debug",1)
%
% Syntax: (这里添加函数的调用格式, `[]`的内容表示可选参数)
%	[Int2,Int2_Sim,gof] = bgDataRead(file_WG,file_fiber ...
%							[, 'Llimit', 0.1 ...
%							 , 'Colormap', 'OSI_rainbow'...
%							 , 'debug', 1...
%							 , 'debug_figure', 0...
%							 , 'debug_txt', 0]);
%
% Params:
%   - file_WG       [required]  [file=*.bgData] 波导.bgData文件路径
%   - file_fiber    [required]  [file=*.bgData] 光纤.bgData文件路径
%   - Llimit        [namevalue] [numeric; >=0; <1] BeamGage 颜色→z轴刻度下限
%   - Colormap      [namevalue] [char] 展示图像的colormap，默认OSI_rainbow，'jet'、'gray'、'OSI_rainbow'
%   - debug         [namevalue] [logical] 是否给出计算过程图像和文字
%   - debug_figure  [namevalue] [logical] 可关闭图像
%   - debug_txt     [namevalue] [logical] 可关闭文字
%
% Return:
%   - Int2      计算实验重叠积分
%   - Int2_Sim  计算拟合重叠积分
%   - gof       WG拟合符合度
%
% Matlab Version: R2024b
%
% Author: oyy
arguments
    file_WG {mustBeFile}
    file_fiber {mustBeFile}
    options.Llimit (1,1) {mustBeLessThan(options.Llimit,1),mustBeGreaterThanOrEqual(options.Llimit,0)} = 0.1
    options.Colormap {mustBeNonzeroLengthText} = 'OSI_rainbow'
    options.debug {mustBeNumericOrLogical} = 0
    options.debug_figure {mustBeNumericOrLogical} = 1
    options.debug_txt {mustBeNumericOrLogical} = 1
end
%% 图片处理
%   mf_WG_original    .bgData读取数据
%   mf_WG             裁剪掉无光强部分并归一化光强
%   im_WG             计算实验重叠积分图片
% 读图
[mf_WG_original,~] = bgDataRead(file_WG,Llimit=options.Llimit);
[mf_fiber_original,~] = bgDataRead(file_fiber,Llimit=options.Llimit);
% 去除无光强部分
mf_WG = beamGageGray64ImgPrepare(mf_WG_original);
mf_fiber = beamGageGray64ImgPrepare(mf_fiber_original);

% 统一中心点位置
% 对比图片大小设置计算范围
% 在大小之间取最大值作为新图进行计算
% 找出二维高斯中心并移动到相同大小的图片中
% warning('off','curvefit:prepareFittingData:removingNaNAndInf');
[size_WG(1),size_WG(2)] = size(mf_WG);
[~,gof,gauss2_WG,~] = gauss2fit(1:size_WG(2),1:size_WG(1),mf_WG);
mid_WG = round([gauss2_WG.Y0,gauss2_WG.X0]);

[size_fiber(1),size_fiber(2)] = size(mf_fiber);
[~,~,gauss2_fiber,~] = gauss2fit(1:size_fiber(2),1:size_fiber(1),mf_fiber);
mid_fiber = round([gauss2_fiber.Y0,gauss2_fiber.X0]);

size0 = max([size(mf_WG);size(mf_fiber)]);

[im_WG,im_fiber] = deal(zeros(size0));
shift = round(size0/2-mid_WG);
im_WG(1:size_WG(1),1:size_WG(2))=mf_WG;im_WG = circshift(im_WG,shift);
shift = round(size0/2-mid_fiber);
im_fiber(1:size_fiber(1),1:size_fiber(2))=mf_fiber;im_fiber = circshift(im_fiber,shift);
%% 计算重叠积分
% 实验图片归一化后重叠积分
Int2 = sum(sum((im_WG.^0.5).*(im_fiber.^0.5)))^2/sum(sum(im_WG))/sum(sum(im_fiber));

% 验证统一性
[fitr_WG,~,gauss_WG,~] = gauss2fit(1:size0(2),1:size0(1),im_WG);
[fitr_fiber,~,gauss_fiber,~] = gauss2fit(1:size0(2),1:size0(1),im_fiber);

% 拟合图片重叠积分
% warning('on','curvefit:prepareFittingData:removingNaNAndInf');
[x,y] = meshgrid(1:size0(2),1:size0(1));
fitr_WG_im = exp(feval(fitr_WG,x,y));
fitr_fiber_im = exp(feval(fitr_fiber,x,y));
Int2_Sim = sum(sum((fitr_WG_im.^0.5).*(fitr_fiber_im.^0.5)))^2/sum(sum(fitr_WG_im))/sum(sum(fitr_fiber_im));

%% 绘图
if options.debug
    if options.debug_txt
        fprintf('WG.X0=%f,fiber.X0=%f\nWG.Y0=%f,fiber.Y0=%f\n',gauss_WG.X0,gauss_fiber.X0,gauss_WG.Y0,gauss_fiber.Y0)
        fprintf('实验重叠积分：%.12f\t损耗：%.4f dB\n',Int2,-10*log10(Int2))
        fprintf('拟合重叠积分：%.12f\t损耗：%.4f dB\n',Int2_Sim,-10*log10(Int2_Sim))
    end
    if options.debug_figure
        if strcmp(options.Colormap,'OSI_rainbow'),load("OSI_rainbow.mat");options.Colormap=OSI_rainbow;end

        figure(1),clf
        tiledlayout(2,2,"TileSpacing","tight");
        % fig1-(1,1)实验波导图
        nexttile;imagesc(mf_WG);title('WG');set(gca,'XTickLabel',[],'YColor','b','XColor','b')
        % fig1-(1,2)实验光纤图，与(1,1)尺寸不同
        nexttile;imagesc(mf_fiber);title('Fiber');set(gca,'YAxisLocation','right','XTickLabel',[],'YColor','b','XColor','b')
        % fig1-(2,1)归一化后实验波导图
        nexttile;imagesc(im_WG);set(gca,'YColor','r','XColor','r');
        % fig1-(2,2)归一化后实验光纤图，与(2,1)尺寸相同
        nexttile;imagesc(im_fiber);set(gca,'YAxisLocation','right','YColor','r','XColor','r');
        colormap(options.Colormap);cb = colorbar;cb.Layout.Tile = 'east';cb.Ticks=linspace(cb.Limits(1),cb.Limits(2),6);

        fit_erro_WG = abs(im_WG-fitr_WG_im)/255;
        fit_erro_fiber = abs(im_fiber-fitr_fiber_im)/255;
        figure(2),clf
        tiledlayout(3,2,"TileSpacing","tight");
        % fig2-(1,1)归一化后实验波导图，同fig1-(2,1)
        nexttile;mesh(im_WG);title('WG');set(gca,'YColor','b','XColor','b','ZColor','b')
        % fig2-(1,2)归一化后实验光纤图，同fig1-(2,2)
        nexttile;mesh(im_fiber);title('Fiber');set(gca,'YColor','r','XColor','r','ZColor','r')
        colormap(options.Colormap);cb = colorbar;cb.Location = 'eastoutside';cb.Ticks=linspace(cb.Limits(1),cb.Limits(2),6);
        % fig2-(2,1)拟合波导图
        nexttile;mesh(fitr_WG_im);set(gca,'YColor','b','XColor','b','ZColor','b');
        % fig2-(2,2)拟合光纤图
        nexttile;mesh(fitr_fiber_im);set(gca,'YColor','r','XColor','r','ZColor','r');
        % fig2-(3,1)拟合波导误差图
        nexttile;imagesc(fit_erro_WG);title('fit error (%)');set(gca,'YColor','b','XColor','b','ZColor','b')
        colormap(options.Colormap);cb = colorbar;cb.Location = 'southoutside';
        cb.Limits(2)=max(max(fit_erro_WG));
        cb.Limits(2)=max(max(fit_erro_WG(fit_erro_WG<0.99*cb.Limits(2))));
        cb.Ticks=linspace(cb.Limits(1),cb.Limits(2),6);cb.Color='b';
        cb.TickLabels=split(num2str(ceil(cb.Ticks*1000)/1000));
        % fig2-(3,1)拟合光纤误差图
        nexttile;imagesc(fit_erro_fiber);title('fit error (%)');set(gca,'YAxisLocation','right','YColor','r','XColor','r','ZColor','r')
        colormap(options.Colormap);cb = colorbar;cb.Location = 'southoutside';
        cb.Limits(2)=max(max(fit_erro_fiber));
        cb.Limits(2)=max(max(fit_erro_fiber(fit_erro_fiber<0.99*cb.Limits(2))));
        cb.Ticks=linspace(cb.Limits(1),cb.Limits(2),6);cb.Color='r';
        cb.TickLabels=split(num2str(ceil(cb.Ticks*1000)/1000));
    end
end
end

%% function
function img = beamGageGray64ImgPrepare(img0)
% 找到光斑列范围
C = sum(img0);
[~,XlocM,w,~]=findpeaks(C,MinPeakHeight=max(C)/4,MinPeakDistance=length(C)/10);
zl = C==0;
XlocL=find(zl(1:XlocM(1)),1,"last");
XlocR=find(zl(XlocM(end):end),1)+XlocM(end)-1;
if isempty(XlocR)||isempty(XlocL)
    XlocL=floor(XlocM-w*1.5);
    XlocR=ceil(XlocM+w*1.5);
end
% 找到光斑行范围
R = sum(img0,2);
[~,YlocM,w]=findpeaks(R,MinPeakHeight=max(R)/4,MinPeakDistance=length(R)/10);
zl = R==0;
YlocL=find(zl(1:YlocM(1)),1,"last");
YlocR=find(zl(YlocM(end):end),1)+YlocM(end)-1;
if isempty(XlocR)||isempty(XlocL)
    YlocL=floor(YlocM-w*1.5);
    YlocR=ceil(YlocM+w*1.5);
end
img = img0(YlocL:YlocR,XlocL:XlocR);
img = padarray(img,[10 10],0,'both');
% % 依据最大值归一化
maxA = double(max(max(img)));
img = double(img)/maxA*255;
end
