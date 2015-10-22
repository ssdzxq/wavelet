clear

file_path = '.\imm_face_db\';
img_path_list = dir(strcat(file_path,'*.jpg'));
img_num = length(img_path_list);

row_delta = 40;
line_delta = 40;
n = img_num; % total of samples
k = 40; % training samples include 40 classes, each class includes 6 samples
numPerClass = n/k;
m = 16*12; % dimension of a single sample
A = zeros(m,n);

if img_num > 0
    for j = 1:img_num
        image_name = img_path_list(j).name;
        image = imread(strcat(file_path,image_name));
        
        image = rgb2gray(image);
        [row,line] = size(image);
        data_new = image(1:row_delta:row,:); % downsampling in row
        data_new = data_new(:,1:line_delta:line); % downsampling in column
%         [filename,pathname] = uiputfile('*.jpg','enter name:');
%         pathfile = [pathname,filename];
%         imwrite(data_new,pathfile,'jpeg');
%         imwrite(data_new,['.\im_downsample\',j],'JPEG'); % save the downsampled image
        A(:,j) = data_new( : );
    end
end

testImage = imread('01-6m.jpg');
testImage = rgb2gray(testImage);
testImage = im2double(testImage);
testImage = imadjust(testImage); % as a test sample
testImage_ds = testImage(1:row_delta:row,:);
testImage_ds = testImage_ds(:,1:line_delta:line);
image_before = testImage_ds;
y = testImage_ds( : );

A = A./repmat(sqrt(sum(A.^2,1)),size(A,1),1); % normalize the columns of A to have unit l2-norm

% (1) without taking occlusion into consideration, so the problem is y = A*x
x0 = A'*y; % x0 is positive, but xp contains negative elements

% (2) for occlusion, y = [A I]*[x e0]


% (3) for occlusion and nonnegative constraint, y = [A I -I]*[x e0+ e0-] = B*c,
% c >= 0




tic
xp = l1eq_pd(x0, A, [], y, 1e-3);
toc

image_after = reshape(A*xp, 12, 16);

figure(1);plot(x0);
figure(2);plot(xp);

% the result of recovery is good, but it can't be classfied to any class.
figure(3);imshow(image_before);
figure(4);imshow(image_after);

% compute the residuals
res = zeros(k,1);
deltol = 1e-3;
for i = 1:k
    delta = zeros(n,1);
    delta(((i-1)*6+1):i*6) = xp(((i-1)*6+1):i*6);
    index = find(abs(delta)<deltol);
    delta(index) = 0;
    temp = y - A*delta;
    res(i) = norm(temp);
end

figure(5);stem(res);
