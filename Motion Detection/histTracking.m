%% histogram-based tracking & edge-based tracking
clear

% particle filter 相关参数的初始化
N = 10;
sigma_noise = 3;
originalPos = [11,40]; % 默认object为方形
speed_x = 10; speed_y = 10;

for i = 1:20
    % 模拟目标运动
    frame = ones(300,300,3); % 注意：图片像素值已转化为double！
    frame(originalPos(1)+(i-1)*speed_x:originalPos(2)+(i-1)*speed_x,originalPos(1)+(i-1)*speed_y:originalPos(2)+(i-1)*speed_y,3) = 0.5; % 取一色块作为运动目标
    obj = frame(originalPos(1)+(i-1)*speed_x:originalPos(2)+(i-1)*speed_x,originalPos(1)+(i-1)*speed_y:originalPos(2)+(i-1)*speed_y,:);    
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
        sample = frame(int16(x_leftdn):int16(x_leftdn+width), int16(y_leftdn):int16(y_leftdn+height),:);
        hsvSample = rgb2hsv(sample);
        hSam = hsvSample(:,:,1)*255;
        sSam = hsvSample(:,:,2)*255;
        vSam = hsvSample(:,:,3)*255;
        
        % 检测方式：hsv (1)设置标准值
        vmin = min(min(vSam)); smin = min(min(hSam)); 
        hmin = min(min(hSam))-20; hmax = max(max(hSam))+20;

        pos_center_x = x_leftdn + width/2;
        pos_center_y = y_leftdn + height/2;
    end
        
    % 更新 particle filter 粒子的中心点坐标
    particles.center = [pos_center_x*ones(1,N); pos_center_y*ones(1,N)] + normrnd(0,sigma_noise,2,N);
    particles.width = width; % width和height可以改变
    particles.height = height;

    figure(2);imshow(frame);
    hold on; scatter(particles.center(1,:),particles.center(2,:),'filled');
    
    axis image off
    drawnow;
%-----------------------------------------------------------------------------------------------------------------
    % 检测方式：hsv---------------------------------------------------------
    % 经过第一轮的初始化设置，得到N个与目标HSV相同的点，进行运动预测.<只关注中心点>
    speed_x = 10; speed_y = 10;
    particlesNew.center = particles.center + [speed_x*ones(1,N);speed_y*ones(1,N)] + normrnd(0,sigma_noise,2,N);
    particlesNew.width = particles.width;
    particlesNew.height = particles.height;
    
%     cnt = 0;
%     for k = 1:N
%         % 检测方式：hsv (2)找到目标范围内的粒子
%         frameHsv = rgb2hsv(frame);
%         checkH = frameHsv(round(particlesNew.center(1,k)),round(particlesNew.center(2,k)),1) * 255;
%         checkS = frameHsv(round(particlesNew.center(1,k)),round(particlesNew.center(2,k)),2) * 255;
%         checkV = frameHsv(round(particlesNew.center(1,k)),round(particlesNew.center(2,k)),3) * 255;
%         sign = checkH >= hmin && checkH <= hmax && checkS >= smin && checkV >= vmin;
%         
%         if sign == 0
%             particlesNew.center(1,k) = 0;
%             particlesNew.center(2,k) = 0;
%         else
%             cnt = cnt + 1;
%         end
%     end
%     % 未考虑粒子的权重
%     pos_center_x = sum(particlesNew.center(1,:))/cnt;
%     pos_center_y = sum(particlesNew.center(2,:))/cnt;
%-------------------------------------------------------------------------------------------------------
    
    % 检测方式：RGB histogram---------------------------------------------
    % (1) 存储object的RGB histogram
%    彩色图像直方图统计：彩色图像一般由RGB三个通道构成，每一个通道由8位构成，最大为255，
% 如果直接根据三个通道每一个不同的值构造直方图很先得很庞大，为256*256*256=256三次方个bins，
% 为简单起见，每一个通道设置8个bins，这样一来，每一个通道最大值256/8=32,
% 即每一个通道划分8bins，每一个bins里面可以存放32个数，0-31,32-63,64-127，.....,224-255等。
% 彩色RGB转化为一维总共8*8*8=512个bins，
%    r = image[(y*W+x)*3] >> R_SHIFT;  
%    g = image[(y*W+x)*3+1] >> G_SHIFT;
%    b = image[(y*W+x)*3+2] >> B_SHIFT;
%   这里R_SHIFT， G_SHIFT，B_SHIFT宏定义5，右移5位，每一个R，G,B值除以32映射到相对应的8个bins中。
% 0-31映射到bins1,32-63映射到bins2中......224-255映射到bins8中.
%   总结：对于每一个RGB像素值，通过计算都可以映射到唯一的index，根据index累加，将相应的核密度权值累加，统计出直方图

    % 设置8*8*8bins存储目标的RGB直方图信息
    R_Shift = -5; G_Shift = -5; B_Shift = -5;
    R_Bins = 8; G_Bins = 8; B_Bins = 8;
%     [obj_height,obj_width,channels] = size(sample);
    obj_r = bitshift(round(sample(:,:,1)*255),R_Shift);
    obj_g = bitshift(round(sample(:,:,2)*255),G_Shift);
    obj_b = bitshift(round(sample(:,:,3)*255),B_Shift);
    obj_index = obj_r*G_Bins*B_Bins + obj_g*B_Bins + obj_b;
    dim = R_Bins*G_Bins*B_Bins;
    
    % 计算颜色分布color distribution pu
    k_r = getWeight(height,width); % 已归一化
    obj_hist = zeros(dim,1);
    for j = 1:width
        for k = 1:height
            % 如果区域内所有的点权重相同，则 +1
            % obj_hist(obj_index(j,k)) = obj_hist(obj_index(j,k)) + 1;
            obj_hist(obj_index(j,k)) = obj_hist(obj_index(j,k)) + k_r(k,j);
        end
    end
    obj_hist = obj_hist/sum(sum(k_r));
    
    % (2) 设置搜索窗大小和目标大小相同，搜索窗为以粒子为中心的矩形，得到N个搜索窗
    searchWin.center = particlesNew.center;
    searchWin.width = particlesNew.width;
    searchWin.height = particlesNew.height;

    % (3) 计算搜索窗的color distribution qu and distance
    dist =ones(N,1);
    for k = 1:N
        % 计算直方图
        rows = int16(searchWin.center(2,k)-searchWin.height/2+1):int16(searchWin.center(2,k)+searchWin.height/2);
        cols = int16(searchWin.center(1,k)-searchWin.width/2+1):int16(searchWin.center(1,k)+searchWin.width/2);
        searchArea = frame(rows,cols,:);
%         [search_height,search_width,channels] = size(searchArea);
        candi_r = bitshift(round(searchArea(:,:,1)*255),R_Shift);
        candi_g = bitshift(round(searchArea(:,:,2)*255),G_Shift);
        candi_b = bitshift(round(searchArea(:,:,3)*255),B_Shift);
        candi_index = candi_r*G_Bins*B_Bins + candi_g*B_Bins + candi_b;
        k_r = getWeight(searchWin.height,searchWin.width);
        candi_hist = zeros(dim,1);
        for j = 1:searchWin.width
            for p = 1:searchWin.height
                candi_hist(candi_index(p,j)) = candi_hist(candi_index(p,j)) + k_r(p,j);
            end
        end
        candi_hist = candi_hist/sum(sum(k_r));
        rou = sum(sqrt(obj_hist.*candi_hist));
        dist(k) = sqrt(1-rou);
    end
    
    % 更新目标中心位置
    index = find(dist==min(dist));
    pos_center_x = searchWin.center(1,index(1));
    pos_center_y = searchWin.center(2,index(1));
    
%-------------------------------------------------------------------------------------------
    % 检测方式：object和candidate的similarity，借助sparse representation
    
    
    
    
 
end;




