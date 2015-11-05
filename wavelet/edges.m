clear
% img = imread('/home/zxq/github/wavelet/01-1m.jpg');
img = imread('/home/zxq/Pictures/screenshot.png');
img = rgb2gray(img);
w = 1/25*ones(5,5);
img = imfilter(img,w);
img = imfilter(img,w);
% img  = img(1:4:end,1:4:end);
[c,s] = wavedec2(img,2,'db2');
dx = detcoef2('h',c,s,1);
dy = detcoef2('v',c,s,1);
grad = abs(dx) + abs(dy);
ang = atan((dy+1e-6)./(dx+1e-6));
dir = ang;
index = ang>=0 & ang<pi/4;
dir(index) = 1;
index = ang>=-pi & ang<-3*pi/4;
dir(index) = 1;
index = ang>=pi/4 & ang<pi/2;
dir(index) = 2;
index = ang>=-3*pi/4 & ang<-pi/2;
dir(index) = 2;
index = ang>=pi/2 & ang<3*pi/4;
dir(index) = 3;
index = ang>=-pi/2 & ang<-pi/4;
dir(index) = 3;
index = ang>=3*pi/4 & ang<=pi;
dir(index) = 4;
index = ang>=-pi/4 & ang<0;
dir(index) = 4;

[m,n] = size(dir);
result = zeros(m,n);

for i = 2:m-1
    for j = 2:n-1
        if dir(i,j) == 1
            if grad(i,j)>=grad(i,j-1) && grad(i,j)>=grad(i,j+1)
                result(i,j) = 1;
            end
        elseif dir(i,j) == 2
            if grad(i,j)>=grad(i-1,j+1) && grad(i,j)>=grad(i+1,j-1)
                result(i,j) = 1;
            end
        elseif dir(i,j) == 3
            if grad(i,j)>=grad(i-1,j) && grad(i,j)>=grad(i+1,j)
                result(i,j) = 1;
            end
        else 
            if grad(i,j)>=grad(i-1,j-1) && grad(i,j)>=grad(i+1,j+1)
                result(i,j) = 1;
            end
        end
    end
end

result = imfilter(result,w);

figure(1)
subplot(2,2,1);imshow(img);title('原图');
subplot(2,2,2);imshow(dx);title('水平高频分量');
subplot(2,2,3);imshow(dy);title('垂直高频分量');
subplot(2,2,4);imshow(result);title('模极大提取');
