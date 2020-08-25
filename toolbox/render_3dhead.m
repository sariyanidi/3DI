function render_3dhead(data, opts, fit_params)

% I is the image -- grascaled and resized
I = data.Is{1};

% 3D rotation parameters
u = fit_params.us(1,:);

% 3D translation parameters
taux = fit_params.taus(1,1);
tauy = fit_params.taus(1,2);
tauz = fit_params.taus(1,3);

% These are the non-rigid parameters of the 3D model fitting
% alpha is shape parameters for identity and epsilon is expression
% parameter
alpha = fit_params.alpha;
epsilon = fit_params.epsilons(1,:);

% rotation matrix
[R] = get_rot_mat(u);

% p is the 3D shape of face -- without view transformation
Deltap = data.model.A*alpha' + data.model.E*epsilon';
p = reshape(data.p0+Deltap, [3, data.N]);

X = p(1,:)';
Y = p(2,:)';
Z = p(3,:)';
Z = -Z;

% The final 3D face shape -- it contains all variations (rigid & non-rigid)
vpose = R*p+[taux; tauy; tauz];

% The 2D location of the dense face shape. 
% Warning: those are for the resized image
xp = data.phix*vpose(1,:)./vpose(3,:) + data.cx;
yp = data.phiy*vpose(2,:)./vpose(3,:) + data.cy;

% These are the 2D location of dense face shape for the original image
xp_orig = xp/data.RESIZE_COEF;
yp_orig = yp/data.RESIZE_COEF;

% These are the landmarks as detected by 3DFAN -- not fitting of our method
px = data.pxs{1};
py = data.pys{1};

% the commands below construct the face rectangle in the image
px0 = min(px);
py0 = min(py);
facewidth = max(px)-min(px);
faceheight = max(py)-min(py);
facesize = (facewidth+faceheight)/2;
px0 = px0-facesize;
py0 = py0-facesize;
rect = [px0 py0 facesize*3 facesize*3];

% here we do a correction that is needed only for visualization purposes
R1 = eul2rotm([ 0 0 deg2rad(180) ]');
vpose_vis = R1*vpose;
vpose_vis(1,:) = vpose_vis(1,:)-mean(vpose_vis(1,:));
vpose_vis(2,:) = vpose_vis(2,:)-mean(vpose_vis(2,:));



% The rest of the code does only visualization
% width and height of the figure that will be plotted
width = 410*4;  
height = 410;

clf
subplot(1,4,1);
fig_handle = sfigure(1);
set(fig_handle, 'position', [0, 0, width, height])

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
set(gca,'LooseInset',get(gca,'TightInset'));
axes('OuterPosition', [0,0,1,1]);

tl = data.model.tl;

objpose.vertices = vpose_vis';
objpose.faces = tl;

rect(1) = max(rect(1), 1);
rect(2) = max(rect(2), 1);


% Panel 1
subplot(1,4,1);
imshow(imcrop(I,rect));


% Panel 2
subplot(1,4,2);
patch(objpose, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'phong'); light ;   axis off;% axis equal;
xlim([-2.5 2.5]);
ylim([-3.5 3.5]);

% Panel 3
subplot(1,4,3);
rect_orig = rect/data.RESIZE_COEF;
imshow(imcrop(data.Iorigs{1},rect_orig));
hold on
plot(xp_orig(1:5:end)-rect_orig(1), yp_orig(1:5:end)-rect_orig(2), '.', 'MarkerSize', 5);
plot(xp_orig(opts.li)-rect_orig(1), yp_orig(opts.li)-rect_orig(2), 'g.', 'MarkerSize', 12, 'LineWidth', 2.4);

% Panel 4
subplot(1,4,4);
rect_orig = rect/data.RESIZE_COEF;
imshow(imcrop(data.Iorigs{1},rect_orig));
hold on
plot(xp_orig(opts.li)-rect_orig(1), yp_orig(opts.li)-rect_orig(2), 'g.', 'MarkerSize', 15, 'LineWidth', 2.4);

results_dir = sprintf('output');
if ~isdir(results_dir)
mkdir(results_dir);
end



set(fig_handle, 'position', [0, 0, width, height])
grid on

set(gca,'xtick',[]);
set(gca,'xticklabel',[]);
set(gca,'ytick',[]);
set(gca,'yticklabel',[]);
set(gca,'LooseInset',get(gca,'TightInset'));

tmp = getframe(gcf);
oI = tmp.cdata;
oI = imresize(oI, 0.5);
imwrite(oI, opts.oimpath);

end


function [R] = get_rot_mat(u)

if size(u,1) == 1
    u = u';
end

theta = norm(u);
unorm = u/theta;
unorm_skew = skew(unorm);
u_skew = skew(u);

I = eye(3);
R = I+sin(theta)*unorm_skew +(1-cos(theta))*(unorm_skew *unorm_skew );


end


function Ux = skew(u)
Ux = [0     -u(3)   u(2);
      u(3)  0       -u(1);
      -u(2) u(1)    0];
end


function rotm = eul2rotm(eul, sequence)
    if (size(eul,1) ~= 3)
        eul = eul';
%         error('eul2rotm: %s', WBM.wbmErrorMsg.WRONG_VEC_DIM);
    end

    if ~exist('sequence', 'var')
        % use the default axis sequence ...
        sequence = 'ZYX';
    end
    rotm = zeros(3,3);

    s_1 = sin(eul(1,1)); % theta_z or theta_z1
    c_1 = cos(eul(1,1));
    s_2 = sin(eul(2,1)); % theta_y
    c_2 = cos(eul(2,1));
    s_3 = sin(eul(3,1)); % theta_x or theta_z2
    c_3 = cos(eul(3,1));

    %% Convert the given Euler angles theta for the x, y and z-axis into the corresponding
    %  direction cosine rotation matrix R, in dependency of the axis rotation sequence for
    %  the multiplication order of the rotation factors:
    % For further details see:
    %   [1] Geometric Tools Engine, Documentation: <http://www.geometrictools.com/Documentation/EulerAngles.pdf>, p. 9 & 16.
    %   [2] MATLAB TOOLBOX FOR RIGID BODY KINEMATICS, Hanspeter Schaub & John L. Junkins,
    %       9th AAS/AIAA Astrodynamics Specialist Conference, AAS 99-139, 1999, <http://hanspeterschaub.info/Papers/mtb1.pdf>, p. 4.
    %   [3] GitHub: ShoolCode/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox,
    %       <https://github.com/ZachDischner/SchoolCode/tree/master/ASEN 5010-Spacecraft Attitude Dynamics and Control/AIAA Software (2nd)/Matlab Toolbox/>
    %   [4] Modelling and Control of Robot Manipulators, L. Sciavicco & B. Siciliano, 2nd Edition, Springer, 2008,
    %       pp. 31-32, formulas (2.18) and (2.20).
    switch sequence
        case 'ZYX'
            %            |c_1*c_2    c_1*s_2*s_3 - s_1*c_3    c_1*s_2*c_3 + s_1*s_3|
            % R(Theta) = |s_1*c_2    s_1*s_2*s_3 + c_1*c_3    s_1*s_2*c_3 - c_1*s_3|
            %            |   -s_2                  c_2*s_3                  c_2*c_3|
            rotm(1,1) =  c_1*c_2;
            rotm(1,2) =  c_1*s_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2*c_3 + s_1*s_3;

            rotm(2,1) =  s_1*c_2;
            rotm(2,2) =  s_1*s_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2*c_3 - c_1*s_3;

            rotm(3,1) = -s_2;
            rotm(3,2) =  c_2*s_3;
            rotm(3,3) =  c_2*c_3;
        case 'ZYZ'
            %            |c_1*c_2*c_3 - s_1*s_3   -c_1*c_2*s_3 - s_1*c_3    c_1*s_2|
            % R(Theta) = |s_1*c_2*c_3 + c_1*s_3   -s_1*c_2*s_3 + c_1*c_3    s_1*s_2|
            %            |             -s_2*c_3                  s_2*s_3        c_2|
            rotm(1,1) =  c_1*c_2*c_3 - s_1*s_3;
            rotm(1,2) = -c_1*c_2*s_3 - s_1*c_3;
            rotm(1,3) =  c_1*s_2;

            rotm(2,1) =  s_1*c_2*c_3 + c_1*s_3;
            rotm(2,2) = -s_1*c_2*s_3 + c_1*c_3;
            rotm(2,3) =  s_1*s_2;

            rotm(3,1) = -s_2*c_3;
            rotm(3,2) =  s_2*s_3;
            rotm(3,3) =  c_2;
        otherwise
            error('eul2rotm: %s', WBM.wbmErrorMsg.UNKNOWN_AXIS_SEQ);
    end
end


