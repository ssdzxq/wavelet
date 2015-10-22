%% tracking
clear
avi = VideoReader('TestVideo_1.avi');
numFrames = avi.NumberOfFrames;

% vmin = 0.5*255;
% smin = 0.5*255;
% hmin = 220/360*255;
% hmax = 260/360*255;
% 找出来的蓝色不准，只识别了墙边蓝色的框子

for i = 1:100
    frame = read(avi,i);
    % 存储第一帧图片中感兴趣区域的颜色
    if i == 1
        figure(1);imshow(frame);
        rect = imrect;
        pos = getPosition(rect);
        x_leftdn = pos(1);
        y_leftdn = pos(2);
        width = pos(3);
        height = pos(4);
        sample = frame(x_leftdn:x_leftdn+width, y_leftdn:y_leftdn+height,:);
        hsvSample = rgb2hsv(sample);
        hSam = hsvSample(:,:,1)*255;
        sSam = hsvSample(:,:,2)*255;
        vSam = hsvSample(:,:,3)*255;
        % 如何设置标准值？
    end
        
    % 获取直方图之后如何在后续帧中进行识别？ 
    % 1.逐点判断，置零？找符合的点块，取一个矩形包含这些点，矩形（或者矩形的中心）即为目标
    % 2.利用meanshift(particle filter),在当前目标附近产生一些样本？？（还是一些点）
    % 代入状态方程进行迭代，再检测下一帧中这些点的hsv值，舍去相差超过阈值的点。
    % 选择差距最小的点作为新的的目标，即下一次迭代的起始位置
    hsvImg = rgb2hsv(frame);
    h = hsvImg(:,:,1)*255;
    s = hsvImg(:,:,2)*255;
    v = hsvImg(:,:,3)*255;
    
    index = find(s<smin | v<vmin);
    s(double(index)) = 0;
    v(double(index)) = 0;
    h(double(index)) = 0;
    
    index = find(h<hmin | h>hmax);
    s(double(index)) = 0;
    v(double(index)) = 0;
    h(double(index)) = 0;
    
    result(:,:,1) = h;
    result(:,:,2) = s;
    result(:,:,3) = v;
    result = hsv2rgb(result/255);
    figure(2);imshow(result);
    
    axis image off
    drawnow;
end;

%% 蓝色试值
% 红0 黄60 绿120 青180 蓝240 品红300
% hz = (240/360*255)*ones(24,24);
% sz = (0.5*255)*ones(24,24);
% vz = (0.5*255)*ones(24,24);
% rz(:,:,1) = hz; rz(:,:,2) = sz; rz(:,:,3) = vz;
% I = hsv2rgb(rz/255);
% imshow(I)

%% 找出第一帧中感兴趣的部分，存取其HSV直方图。可以用光流法找吗

%% 利用鼠标，手动选择感兴趣的区域



