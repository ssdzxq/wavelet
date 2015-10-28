clear
load woman
[c,s] = wavedec2(X,2,'db1');
[m,n] = size(X);
J = c(1:(m/4*n/4));
P = reshape(J,m/4,n/4);
imshow(uint8(P),[]);

I = X;
[m,n] = size(X);
[LO_D,HI_D,LO_R,HI_R] = wfilters('db1');
c1 = c(1:m/2*n/2);
I_cols = reshape(I,1,m*n);
J1 = downsample(I,2);
J1 = downsample(J1',2);
J1 = J1';
J1 = reshape(J1,1,m/2*n/2);
[c2,s2] = wavedec(J1,1,'db1');

% I_cols = reshape(downsample(I,2),1,m*n/2); % 行采样，按列一维化
% I_rows = reshape(downsample(I',2),1,m*n/2); % 列采样，按行一维化
% [c11,s11] = wavedec(I_cols,2,'db1'); % y方向的低频和高频
% [c12,s12] = wavedec(I_rows,2,'db2'); % x方向的低频和高频
% 
% I11 = reshape(c11(1:m/2*n/2),m/2,n/2); % y方向的低频
% I12 = reshape(c11(m/2*n/2+1:m*n/2),m/2,n/2); % y方向的高频
% I13 = reshape(c12(1:m/2*n/2),m/2,n/2); % x方向的低频
% I14 = reshape(c12(m/2*n/2+1:m*n/2),m/2,n/2); % y方向的高频
% figure(1);
% subplot(2,2,1);imshow(I11);
% subplot(2,2,2);imshow(I12);
% subplot(2,2,3);imshow(I13);
% subplot(2,2,4);imshow(I14);
% 
% f11 = filter(LO_D,double(I_cols)); % y方向的低频
% f12 = filter(HI_D,double(I_cols)); % y方向的高频
% f13 = filter(LO_D,double(I_rows)); % x方向的低频
% f14 = filter(HI_D,double(I_rows)); % x方向的高频
% J11 = reshape(f11,m/2,n/2);
% J12 = reshape(f12,m/2,n/2);
% J13 = reshape(f13,m/2,n/2);
% J14 = reshape(f14,m/2,n/2);
% figure(2);
% subplot(2,2,1);imshow(J11);
% subplot(2,2,2);imshow(J12);
% subplot(2,2,3);imshow(J13);
% subplot(2,2,4);imshow(J14);