%% Image Read
I = im2double(imread('C:\Users\hjjiang\Desktop\大三上\数字图像处理\作业\综合作业1\phone.bmp'));

%% Block-wise Spatial Enhance
pad_size = 12;
pad_I = padarray(I,[12,12],1);
imgsize = size(pad_I);
height = imgsize(1);
width = imgsize(2);
enhance_block_size = 12;
enhance_block_half_size = 6;
xcenter = enhance_block_half_size:8:(width-enhance_block_half_size);
ycenter = enhance_block_half_size:8:(height-enhance_block_half_size);
xnum = size(xcenter);
ynum = size(ycenter);
block_enhance = cell(ynum(2), xnum(2));
enhance_image = zeros(height,width);
for i = 1:xnum(2)
    for j = 1:ynum(2)
        block_enhance{j,i} = pad_I(ycenter(j)-enhance_block_half_size+1:ycenter(j)+enhance_block_half_size,xcenter(i)-enhance_block_half_size+1:xcenter(i)+enhance_block_half_size);  % 第一个维度是高
        block_enhance{j,i} = 1 - block_enhance{j,i}; % 反色
        block_enhance{j,i} = histeq(block_enhance{j,i}); % 直方图增强
        block_enhance{j,i} = block_enhance{j,i}.^1.5;
        enhance_image(ycenter(j)-enhance_block_half_size+1:ycenter(j)+enhance_block_half_size,xcenter(i)-enhance_block_half_size+1:xcenter(i)+enhance_block_half_size) = block_enhance{j,i};
    end
end

%% DFT Compute
block_size = 32;
block_half_size = 16;
xcenter = block_half_size:8:(width-block_half_size);
ycenter = block_half_size:8:(height-block_half_size);
xnum = size(xcenter);
ynum = size(ycenter);
block_32 = cell(ynum(2), xnum(2));
block_32_DFT = cell(ynum(2), xnum(2));
block_32_DFT_amplitude = cell(ynum(2), xnum(2));
block_32_DFT_mean_amplitude = zeros(ynum(2), xnum(2));
for i = 1:xnum(2)
    for j = 1:ynum(2)
        block_32{j,i} = enhance_image(ycenter(j)-block_half_size+1:ycenter(j)+block_half_size,xcenter(i)-block_half_size+1:xcenter(i)+block_half_size);  % 第一个维度是高
        
        % Quick DFT
        vector = 0:1:block_size-1;
        matrix = 2*pi.*(vector'*vector)./block_size;
        G1 = cos(matrix) - 1i*sin(matrix);
        G2 = G1;
        block_DFT = G1*block_32{j,i}*G2;
        block_DFT_shift = fftshift(block_DFT);
        block_32_DFT{j,i} = block_DFT_shift;  
        block_32_DFT_amplitude{j,i} = abs(block_32_DFT{j,i});
        block_32_DFT_mean_amplitude(j,i) = mean(mean(block_32_DFT_amplitude{j,i}));
    end
end
coordinate = block_32_DFT_mean_amplitude > 4; % 有指纹的区域的中心坐标

%% Orientation and Frequency Compute
y0 = 0; x0 = 0; y1 = 0; x1 = 0; y2 = 0; x2 = 0;
v = zeros(ynum(2), xnum(2));
u = zeros(ynum(2), xnum(2));
block_degree = ones(ynum(2), xnum(2));
block_frequency = zeros(ynum(2), xnum(2));
for i = 1:xnum(2)
    for j = 1:ynum(2)
        if coordinate(j,i) == 1
            temp_DFT_amplitude = block_32_DFT_amplitude{j,i};
            temp_DFT_amplitude(16:17,16:17) = 0;

            [y1,x1] = find(temp_DFT_amplitude == max(max(temp_DFT_amplitude)),1);
            temp_DFT_amplitude(y1,x1) = 0;
            [y2,x2] = find(temp_DFT_amplitude == max(max(temp_DFT_amplitude)),1);

            % 频率估计 根据小作业4中的正弦波得到 f=(两个峰值距离)/块边长2倍
            block_frequency(j,i) = sqrt((y1-y2)^2+(x1-x2)^2);
            %block_frequency(j,i) = sqrt((y1)^2+(x1)^2)/32;

            % 方向的方程为 y = ax + b
            if y1 ~= y2
                %a = (x2 - x1)/(y2 - y1);
                %orientaion_32_block = PaintBlockOrientation(a);
                block_degree(j,i) = atand((x2-x1)/(y2-y1));
                v(j,i) = sin((block_degree(j,i)/180)*pi);
                u(j,i) = cos((block_degree(j,i)/180)*pi);
            else % y1 = y2
                block_degree(j,i) = 270;
                v(j,i) = sin((block_degree(j,i)/180)*pi);
                u(j,i) = cos((block_degree(j,i)/180)*pi);
            end
        end
    end
end

%% Smooth Orientation Block
W = fspecial('gaussian', 5, 1.5);
size_w = (5-1)/2;
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
W = fspecial('gaussian', 5, 1);
smooth_block_frequency = imfilter(block_frequency,W);

%% Gabor Filtering
enhanced_block= cell(ynum(2),xnum(2));
enhanced_image = zeros(628,384);
for i = 1:xnum(2)
    for j = 1:ynum(2)
        if coordinate(j,i) == 1
            orientation = smooth_block_degree(j,i)+90;
            wavelength = 17.5*pi/smooth_block_frequency(j,i);
            [mag, phase] = imgaborfilt(block_32{j,i},wavelength,orientation);
            enhanced_block{j,i} = mag.*cos(phase+90);
            enhanced_image(ycenter(j)-3:ycenter(j)+4,xcenter(i)-3:xcenter(i)+4) = enhanced_block{j,i}(13:20,13:20);
        end
    end
end

%% Show Result
smooth_block_frequency_image = zeros(628,384);
block_frequency_image = zeros(628,384);
for i = 1:xnum(2)
    for j = 1:ynum(2)
        if coordinate(j,i) == 1
            temp = zeros(8);
            temp(:) = smooth_block_frequency(j,i);
            smooth_block_frequency_image(ycenter(j)-3:ycenter(j)+4,xcenter(i)-3:xcenter(i)+4) = temp;
            temp(:) = block_frequency(j,i);
            block_frequency_image(ycenter(j)-3:ycenter(j)+4,xcenter(i)-3:xcenter(i)+4) = temp;
        end
    end
end

figure(2);
subplot(2,4,1),imshow(pad_I,[]),title('Original Image')
subplot(2,4,2),imshow(enhance_image,[]),title('Spatial Enhanced Image')
[x, y] = meshgrid(block_half_size:8:width-block_half_size, height-block_half_size:-8:block_half_size);
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
