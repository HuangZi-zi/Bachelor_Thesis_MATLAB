% 使用tum数据集
PATCH_SIZE = 20;
% 校准数据
K_rgb=[5.4886723733696215e+02 0. 3.1649655835885483e+02;
    0. 5.4958402532237187e+02 2.2923873484682150e+02;
    0. 0. 1.];
dist_coeffs_rgb=[7.5390330809303213e-02 -2.8252164126047641e-02 ...
    -9.5892915347463477e-03 -4.6479123130854875e-04 ...
    -1.4030650062583919e-01];
K_ir=[5.7592685448804468e+02 0. 3.1515026356388171e+02;
    0. 5.7640791601093247e+02 2.3058580662101753e+02;
    0. 0. 1.];
dist_coeffs_ir=[-4.5096285298902479e-02 2.3953450975051985e-01 ...
    -9.5873062860576726e-03 -2.9427017385906013e-03 ...
    -1.0147872415919272e+00];
R_stereo=[9.9993556622059065e-01 -1.1131100247700660e-02 2.2275579414240366e-03;
    1.1108891478543052e-02 9.9989080121799856e-01 9.7456744121138719e-03;
    -2.3357947736726680e-03 -9.7203007620451660e-03 9.9995002865936788e-01];
t_stereo=[1.1497548441022023e+01; 3.5139088879273231e+01; 2.1887459420807019e+01];
fx_ir = K_ir(1,1);
fy_ir = K_ir(2,2);
cx_ir = K_ir(1,3);
cy_ir = K_ir(2,3);
fx_rgb = K_rgb(1,1);
fy_rgb = K_rgb(2,2);
cx_rgb = K_rgb(1,3);
cy_rgb = K_rgb(2,3);

% 维度数据
rgb_img=imread("C:\Users\YawnFun\Downloads\rgbd_dataset_freiburg1_floor\rgb\1305033539.743168.png");
d_img=imread("C:\Users\YawnFun\Downloads\rgbd_dataset_freiburg1_floor\depth\1305033539.743179.png");
[height,width,~]=size(rgb_img);
nr_horizontal_cells = width/PATCH_SIZE;
nr_vertical_cells = height/PATCH_SIZE;

% 预计算
X_pre=cast(zeros(height,width),'uint16');
Y_pre=cast(zeros(height,width),'uint16');
for r=0:1:height-1
    for c=0:1:width-1
        X_pre(r+1,c+1) = (c-cx_ir)/fx_ir;
        Y_pre(r+1,c+1) = (r-cy_ir)/fy_ir;
    end
end

cell_map=zeros(height,width);
for r=1:1:height
    cell_r = floor((r-1)/PATCH_SIZE);
    local_r = rem((r-1),PATCH_SIZE);
    for c=1:1:width
        cell_c = floor((c-1)/PATCH_SIZE);
        local_c = rem((c-1),PATCH_SIZE);
        cell_map(r,c) = (cell_r*nr_horizontal_cells+cell_c)*PATCH_SIZE*PATCH_SIZE + local_r*PATCH_SIZE + local_c;
    end
end

% X=zeros(height,width);
% Y=zeros(height,width);
% X_t=zeros(height,width);
% Y_t=zeros(height,width);
% cloud_array=zeros(width*height,3);
% cloud_array_organized=zeros(width*height,3);

color_code = rand(500, 3);

% 将深度图反投影成点云
X = X_pre.*d_img;
Y = Y_pre.*d_img;
X_t = (R_stereo(1,1))*X+(R_stereo(1,2))*Y+(R_stereo(1,3))*d_img + t_stereo(1);
Y_t = (R_stereo(2,1))*X+(R_stereo(2,2))*Y+(R_stereo(2,3))*d_img + t_stereo(2);
d_img = (R_stereo(3,1))*X+(R_stereo(3,2))*Y+(R_stereo(3,3))*d_img + t_stereo(3);
[U,V,cloud_array]=projectPointCloud(X_t, Y_t, d_img, fx_rgb, fy_rgb, cx_rgb, cy_rgb, t_stereo(3));

seg_rz = zeros(height,width,3);
seg_output = zeros(height,width);

% Run CAPE
% int nr_planes, nr_cylinders;
% vector<PlaneSeg> plane_params;
% vector<CylinderSeg> cylinder_params;

cloud_array_organized=organizePointCloudByCell(cloud_array,  cell_map);
% plane_detector->process(cloud_array_organized, nr_planes, nr_cylinders, seg_output, plane_params, cylinder_params);


%% Function projectPointCloud
function [U, V, cloud_array] = projectPointCloud(X, Y, Z, fx_rgb, fy_rgb, cx_rgb, cy_rgb, z_min)
    [height, width] = size(X);

    % Project to image coordinates
    U = X ./ Z;
    V = Y ./ Z;
    U = U * fx_rgb + cx_rgb;
    V = V * fy_rgb + cy_rgb;

    % Reusing U as cloud index
    %U = V * width + U + 0.5;

    cloud_array = zeros(height * width, 3);

    it = 1;
    for r = 1:height
        for c = 1:width
            z = Z(r, c);
            u = U(r, c);
            v = V(r, c);
            if (z > z_min && u > 0 && v > 0 && u < width && v < height)
                id = floor(v) * width + u;
                cloud_array(id, 1) = X(r, c);
                cloud_array(id, 2) = Y(r, c);
                cloud_array(id, 3) = z;
            end
            it = it + 1;
        end
    end
end

% function [U, V, cloud_array] = projectPointCloud(X, Y, Z,  fx_rgb, fy_rgb, cx_rgb, cy_rgb, z_min)
% [height,width,~]=size(X);
% % Project to image coordinates
% U=X./Z;
% V=Y./Z;
% U = U*fx_rgb + cx_rgb;
% V = V*fy_rgb + cy_rgb;
% % Reusing U as cloud index
% % U = V*width + U + 0.5;
% 
% for r=0:1:height-1
%     sx = (r+1);
%     sy = (r+1);
%     sz = (r+1);
%     u_ptr = U(r+1);
%     v_ptr = V(r+1);
%     for c=0:1:width-1
%         z = Z(sz+c+1);
%         u = u_ptr(c+1);
%         v = v_ptr(c+1);
%         if z>z_min && u>0 && v>0 && u<width && v<height
%             id = floor(v)*width + u;
%             cloud_array(id,1) = sx(c+1);
%             cloud_array(id,2) = sy(c+1);
%             cloud_array(id,3) = z;
%         end
%     end
% end
% end
%% Function organizePointCloudByCell
function cloud_out=organizePointCloudByCell(cloud_in, cell_map)
    width = size(cell_map, 2);
    height = size(cell_map, 1);
    mxn = width * height;
    mxn2 = 2 * mxn;

    cloud_out = zeros(mxn, size(cloud_in, 2)); % Initialize output cloud

    it = 1; % iterate
    for r = 1:height
        for c = 1:width
            id = cell_map(r, c)+1; % index
            cloud_out(id, :) = cloud_in(it, :);
            cloud_out(mxn + id, :) = cloud_in(mxn + it, :);
            cloud_out(mxn2 + id, :) = cloud_in(mxn2 + it, :);
            it = it + 1;
        end
    end
end
% function organizePointCloudByCell(cloud_in, cloud_out,cell_map)
% 
% width = cell_map.cols;
% height = cell_map.rows;
% mxn = width*height;
% mxn2 = 2*mxn;
% 
%     int id, it(0);
%     int * cell_map_ptr;
%     for r=0:1:height-1
%         cell_map_ptr = cell_map(r+1);
%         for c=0:1:width-1
%             id = cell_map_ptr(c+1);
%             cloud_out.data() + id = *(cloud_in.data() + it);
%             cloud_out.data() + mxn + id) = *(cloud_in.data() + mxn + it);
%             *(cloud_out.data() + mxn2 + id) = *(cloud_in.data() + mxn2 + it);
%             it++;
%         }
%     }
% }
