% 1) Obtain the Basel 2009 face model: https://faces.dmi.unibas.ch/bfm/index.php?nav=1-0&id=basel_face_model
% 2) Copy the 01_MorphableModel.mat file of the Basel model in this folder
% 3) Obtain the expression model from
% http://www.cbsr.ia.ac.cn/users/xiangyuzhu/projects/3DDFA/Database/300W-3D/main.htm
% (it's in the first of the two links, i.e.
% https://drive.google.com/file/d/0B7OEHD3T4eCkRFRPSXdFWEhRdlE/view?usp=sharing )
% 4) Copy the Model_Exp.mat in this folder (Model_Exp.mat is in the .zip file you downloaded
% in step 3)
% 5) Run this file: create_combined_model.m


w = load('w.dat');
tl = load('tl.dat');

model2 = load('01_MorphableModel.mat');

vertex_idx = load('./BaselIDS', '-ASCII')+1;

lix = (vertex_idx-1)*3+1;
liy = (vertex_idx-1)*3+2;
liz = (vertex_idx-1)*3+3;
li = [lix'; liy'; liz'];
li = li(:);

model2.texPC = model2.texPC(li,:);
model2.texMU = model2.texMU(li);

model2.shapePC = model2.shapePC(li,:);
model2.shapeMU = model2.shapeMU(li,:);
       
model.texPC = double(model2.texPC(w,:));
model.texMU = double(model2.texMU(w));
model.texEV = double(model2.texEV);

model.shapePC = double(model2.shapePC(w,:));
model.shapeMU = double(model2.shapeMU(w));
model.shapeEV = double(model2.shapeEV);
model.tl = tl;

ST = model.texPC;
STR = ST(1:3:end,:);
STG = ST(2:3:end,:);
STB = ST(3:3:end,:);

ST = 0.2989 * double(STR) + 0.5870 * double(STG)+ 0.1140 *double(STB);
model.texPC = ST;

texMU = model.texMU;
model.texMU = 0.2989 * texMU(1:3:end) + 0.5870 * texMU(2:3:end) + 0.1140 * texMU(3:3:end);

load Model_Exp.mat

model.expMU = mu_exp;%(w);
model.expPC = w_exp;%(w,:);
model.expEV = sigma_exp;


K = size(w_exp,2);

norms = zeros(K,1);

for k=1:K
    norms(k) = norm(model.expPC(:,k));
    model.expEV(k) = model.expEV(k)*norms(k);
    model.expPC(:,k) = model.expPC(:,k)/norms(k);
end


model.expPC = model.expPC(w,:);

save('Basel2009_zreduced_grayscale_3DFFA.mat', 'model');







