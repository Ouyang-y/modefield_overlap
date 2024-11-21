%% 
% 输入：
%   - bgDataPath    .bgData文件夹路径
%   - fiberPath     fiber对应的.bgData路径
%
clear
bgDataPath = '';
fiberPath = '';
% 处理光纤

[mf_WG_original,~] = bgDataRead("1.bgData");
[mf_fiber_original,~] = bgDataRead("2.bgData");
mf_WG = beamGageGray64ImgPrepare(mf_WG_original);
mf_fiber = beamGageGray64ImgPrepare(mf_fiber_original);

%% 统一中心点位置
warning('off','curvefit:prepareFittingData:removingNaNAndInf');
% 找出二维高斯中心并移动到相同大小的图片中
[size_WG(1),size_WG(2)] = size(mf_WG);
[f,g,gauss2_WG] = gauss2fit(1:size_WG(2),1:size_WG(1),log(double(mf_WG)));
midrc = round([size_WG(1),size_WG(2)]/2);
mid_WG = round([gauss2_WG.X0,gauss2_WG.Y0]);
size1 = max([mid_WG;midrc])*2;

[size_fiber(1),size_fiber(2)] = size(mf_fiber);
[~,~,gauss2_fiber] = gauss2fit(1:size_fiber(2),1:size_fiber(1),log(double(mf_fiber)));
midrc = round([size_fiber(1),size_fiber(2)]/2);
mid_fiber = round([gauss2_fiber.X0,gauss2_fiber.Y0]);
size2 = max([mid_fiber;midrc])*2;

size0 = max([size1;size2]);

im_WG = zeros(size0);im_fiber = im_WG;
shift = size0/2-mid_WG;
im_WG(1:size_WG(1),1:size_WG(2))=mf_WG;im_WG = circshift(im_WG,shift(end:-1:1));
shift = size0/2-mid_fiber;
im_fiber(1:size_fiber(1),1:size_fiber(2))=mf_fiber;im_fiber = circshift(im_fiber,shift(end:-1:1));

%% 计算重叠积分
Int2 = sum(sum((im_WG.^0.5).*(im_fiber.^0.5)))^2/sum(sum(im_WG))/sum(sum(im_fiber));
fprintf('重叠积分：%.12f\n损耗：%.4f dB\n',Int2,-10*log10(Int2))

