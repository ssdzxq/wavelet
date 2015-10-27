%% tracking
clear

% 追踪视频中的目标
% avi = VideoReader('TestVideo_1.avi');
% numFrames = avi.NumberOfFrames;
% vmin = 0.5*255;
% smin = 0.5*255;
% hmin = 220/360*255;
% hmax = 260/360*255;
% 找出来的蓝色不准，只识别了墙边蓝色的框子

% particle filter 相关参数的初始化
N = 20;
sigma_noise = 3;

for i = 1:20
    
%     frame = read(avi,i);
    
    % 模拟目标运动
    frame = ones(300,300,3);
    frame(10+(i-1)*10:40+(i-1)*10,10+(i-1)*10:40+(i-1)*10,3) = 0.5;
    
    % 存储第一帧图片中感兴趣区域的颜色
    if i == 1
        figure(1);imshow(frame);
        % 利用鼠标，手动选择感兴趣的区域
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
        
        % 用于逐点判断的标准值
        vmin = mean(mean(vSam)); smin = mean(mean(hSam)); 
        hmin = mean(mean(hSam))-20; hmax = mean(mean(hSam))+20;
        
        % 用于particle filter的参数设置
        pos_center_x = x_leftdn + width/2;
        pos_center_y = y_leftdn + height/2;
    end
        
    % 更新 particle filter 粒子的坐标
    particles.position = [pos_center_x*ones(1,N); pos_center_y*ones(1,N)] + normrnd(0,sigma_noise,2,N);
        
    % 获取直方图之后如何在后续帧中进行识别？ 
    % 1.逐点判断，置零？找符合的点块，取一个矩形包含这些点，矩形（或者矩形的中心）即为目标
    % 2.利用meanshift(particle filter),在当前目标附近产生一些样本？？（还是一些点？）
    % 代入状态方程进行迭代，再检测下一帧中这些点的hsv值，舍去相差超过阈值的点。
    % 选择差距最小的点作为新的的目标，即下一次迭代的起始位置
    
    % 逐点判断-------------------------------------------------------------
%     hsvImg = rgb2hsv(frame);
%     h = hsvImg(:,:,1)*255; s = hsvImg(:,:,2)*255; v = hsvImg(:,:,3)*255;
%     
%     index = find(s<smin | v<vmin);
%     s(double(index)) = 0; v(double(index)) = 0; h(double(index)) = 0;
%     index = find(h<hmin | h>hmax);
%     s(double(index)) = 0; v(double(index)) = 0; h(double(index)) = 0;
%     result(:,:,1) = h; result(:,:,2) = s; result(:,:,3) = v;
%     result = hsv2rgb(result/255);
%     figure(2);imshow(result); 

    % particle filter------------------------------------------------------
    figure(2);imshow(frame);
    hold on; scatter(particles.position(1,:),particles.position(2,:),'filled');
    
    axis image off
    drawnow;
    
    % particle filter------------------------------------------------------
    % 经过第一轮的初始化设置，得到N个与目标HSV相同的点，进行运动预测
    particlesNew.position = particles.position + [10*ones(1,N);10*ones(1,N)] + normrnd(0,sigma_noise,2,N);
    cnt = 0;
    for k = 1:N
        frameHsv = rgb2hsv(frame);
        checkH = frameHsv(round(particlesNew.position(1,k)),round(particlesNew.position(2,k)),1) * 255;
        checkS = frameHsv(round(particlesNew.position(1,k)),round(particlesNew.position(2,k)),2) * 255;
        checkV = frameHsv(round(particlesNew.position(1,k)),round(particlesNew.position(2,k)),3) * 255;
        sign = checkH >= hmin && checkH <= hmax && checkS >= smin && checkV >= vmin;
        if sign == 0
            particlesNew.position(1,k) = 0;
            particlesNew.position(2,k) = 0;
        else
            cnt = cnt + 1;
        end
    end
    % 未考虑粒子的权重
    pos_center_x = sum(particlesNew.position(1,:))/cnt; 
    pos_center_y = sum(particlesNew.position(2,:))/cnt;
end;

%% 蓝色试值
% 红0 黄60 绿120 青180 蓝240 品红300
% hz = (240/360*255)*ones(24,24);
% sz = (0.5*255)*ones(24,24);
% vz = (0.5*255)*ones(24,24);
% rz(:,:,1) = hz; rz(:,:,2) = sz; rz(:,:,3) = vz;
% I = hsv2rgb(rz/255);
% imshow(I)





