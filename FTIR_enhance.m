%% Image Read
I = im2double(imread('C:\Users\hjjiang\Desktop\大三上\数字图像处理\作业\综合作业1\FTIR.bmp'));

%% Image Segmentation
pad_I = padarray(I(:,1:387),[32,32],0.9961); % 388列的数据有一定干扰，扩充后可以看到一条黑线，因此去除 
imgsize = size(pad_I);
height = imgsize(1);
width = imgsize(2);
xcenter = 36:8:(width-32);
ycenter = 36:8:(height-32);
xnum = size(xcenter);
ynum = size(ycenter);
block_32 = cell(ynum(2), xnum(2));
block_32_DFT = cell(ynum(2), xnum(2));
block_32_DFT_amplitude = cell(ynum(2), xnum(2));
block_32_DFT_mean_amplitude =zeros(ynum(2), xnum(2));
for i = 1:xnum(2)
    for j = 1:ynum(2)
        block_32{j,i} = pad_I(ycenter(j)-15:ycenter(j)+16,xcenter(i)-15:xcenter(i)+16);  % 第一个维度是高
        block_32{j,i} = histeq(block_32{j,i});
        % Quick DFT
        vector = 0:1:31;
        matrix = 2*pi.*(vector'*vector)./32;
        G1 = cos(matrix) - 1i*sin(matrix);
        G2 = G1;
        block_DFT = G1*block_32{j,i}*G2;
        block_DFT_shift = fftshift(block_DFT);
        block_32_DFT{j,i} = block_DFT_shift;  
        block_32_DFT_amplitude{j,i} = abs(block_32_DFT{j,i});
        block_32_DFT_mean_amplitude(j,i) = mean(mean(block_32_DFT_amplitude{j,i}));
    end
end

% 有指纹的区域的中心坐标
coordinate = block_32_DFT_mean_amplitude > 1;

%% Orientation and Frequency Compute
y0 = 0; x0 = 0; y1 = 0; x1 = 0; y2 = 0; x2 = 0;
v = zeros(47,48);
u = zeros(47,48);
block_degree = zeros(47,48);
block_frequency = zeros(47,48);
for i = 1:xnum(2)
    for j = 1:ynum(2)
        if coordinate(j,i) == 1
            temp_DFT_amplitude = block_32_DFT_amplitude{j,i};
            temp_DFT_amplitude(16:17,16:17) = 0;
            [y1,x1] = find(temp_DFT_amplitude == max(max(temp_DFT_amplitude)),1);
            temp_DFT_amplitude(y1,x1) = 0;
            [y2,x2] = find(temp_DFT_amplitude == max(max(temp_DFT_amplitude)),1);
            
            % 频率估计 根据小作业4中的正弦波得到 f与两个峰值成正比
            block_frequency(j,i) = sqrt((y1-y2)^2+(x1-x2)^2);
            
            % 方向的方程为 y = ax + b
            if y1 ~= y2
                block_degree(j,i) = atand((x2-x1)/(y2-y1));
                v(j,i) = sin((block_degree(j,i)/180)*pi);
                u(j,i) = cos((block_degree(j,i)/180)*pi);
            else % y1 = y2
                if block_degree(j-1,i) > 0 || block_degree(j,i+1) > 0 
                    block_degree(j,i) = 90;
                else
                    block_degree(j,i) = 270;
                end
                v(j,i) = sin((block_degree(j,i)/180)*pi);
                u(j,i) = cos((block_degree(j,i)/180)*pi);
            end
        end
    end
end

%% Smooth Orientation Block

W = fspecial('gaussian', 9, 5);
size_w = (9-1)/2;
phix_org = cos(2*(block_degree/180)*pi);
phiy_org = sin(2*(block_degree/180)*pi);

phi_matrix = phix_org + 1i*phiy_org;
phi_matrix = padarray(phi_matrix,[size_w,size_w]);

smooth_block_degree = block_degree;
smooth_v = v;
smooth_u = u;

for i = 1+size_w:xnum(2)+size_w
    for j = 1+size_w:ynum(2)+size_w
        if coordinate(j-size_w,i-size_w) == 1
            phi_matrix(j-size_w,i-size_w) = sum(sum(W.*phi_matrix(j-size_w:j+size_w, i-size_w:i+size_w)));
            smooth_block_degree(j-size_w,i-size_w) = 0.5*180*angle(phi_matrix(j-size_w,i-size_w))/pi;
            smooth_v(j-size_w,i-size_w) = sin((smooth_block_degree(j-size_w,i-size_w)/180)*pi);
            smooth_u(j-size_w,i-size_w) = cos((smooth_block_degree(j-size_w,i-size_w)/180)*pi);
        end
    end
end

%% Smooth Frequency Block
W = fspecial('gaussian',5,5);
smooth_block_frequency = imfilter(block_frequency, W);

%% Gabor Filtering
enhanced_block= cell(ynum(2),xnum(2));
enhanced_image = pad_I;

for i = 1:xnum(2)
    for j = 1:ynum(2)
        if coordinate(j,i) == 1
            orientation = smooth_block_degree(j,i) + 90;
            wavelength = 20*pi/smooth_block_frequency(j,i);
            [mag, phase] = imgaborfilt(block_32{j,i},wavelength,orientation);
            enhanced_block{j,i} = mag.*cos(phase+90);
            enhanced_image(ycenter(j)-3:ycenter(j)+4,xcenter(i)-3:xcenter(i)+4) = enhanced_block{j,i}(13:20,13:20);
        end
    end
end


%% Show Result
smooth_block_frequency_image = zeros(438,451);
block_frequency_image = zeros(438,451);
fingerprint_area = zeros(438,451);
for i = 1:xnum(2)
    for j = 1:ynum(2)
        if coordinate(j,i) == 1
            temp = zeros(8);
            temp(:) = smooth_block_frequency(j,i);
            smooth_block_frequency_image(ycenter(j)-3:ycenter(j)+4,xcenter(i)-3:xcenter(i)+4) = temp;
            temp(:) = block_frequency(j,i);
            block_frequency_image(ycenter(j)-3:ycenter(j)+4,xcenter(i)-3:xcenter(i)+4) = temp;
            temp(:) = coordinate(j,i);
            fingerprint_area(ycenter(j)-3:ycenter(j)+4,xcenter(i)-3:xcenter(i)+4) = temp;
        end
    end
end

figure(1);
subplot(2,4,1),imshow(pad_I),title('Original Image')
subplot(2,4,2),imshow(fingerprint_area),title('Fingerprint Area')
[x, y] = meshgrid(36:8:451-32, 438-32:-8:36);
subplot(2,4,3),quiver(x,y,u,v),title('Orientation Image Before Smooth')
subplot(2,4,4),quiver(x,y,smooth_u,smooth_v),title('Orientation Image After Smooth')
subplot(2,4,5),imshow(block_frequency_image,[]),title('Frequency Image Before Smooth')
subplot(2,4,6),imshow(smooth_block_frequency_image,[]),title('Frequency Image After Smooth')
subplot(2,4,7),imshow(enhanced_image),title('Enhanced Image')
ax(1) = subplot(2,4,1);
ax(2) = subplot(2,4,2);
ax(3) = subplot(2,4,3);
ax(4) = subplot(2,4,4);
ax(5) = subplot(2,4,5);
ax(6) = subplot(2,4,6);
ax(7) = subplot(2,4,7);
linkaxes(ax, 'xy');
