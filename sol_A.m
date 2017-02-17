clear;
% CIGs parameters 
off0   =-500;
offmax = 500;
doff   = 10;

noffs = (offmax-off0)/doff+1;

nCIGs = 21;

% Model parameters
nx = 501;
nz = 301;

len_lmd = (nx-1)*(nz-1); 

% Size of matrix A and vector b

A_row = noffs*nCIGs;

A_col = len_lmd;

b_row = A_row;

Av_fileID = fopen('Av_file');
Av = fread(Av_fileID,'float');
fclose(Av_fileID);

Ai_fileID = fopen('Ai_file');
Ai = fread(Ai_fileID,'int');
fclose(Ai_fileID);

Aj_fileID = fopen('Aj_file');
Aj = fread(Aj_fileID,'int');
fclose(Aj_fileID);

b_fileID = fopen('b_file');
b = fread(b_fileID,'float');
fclose(b_fileID);

A_full = sparse(Ai,Aj,Av,A_row,A_col);
%options = optimset('LargeScale','on','Display','iter','MaxIter', 1);

%[xx, RESNORM, RESIDUAL, EXITFLAG] = lsqlin(A_full,b,[],[],[],[],[],[],[],options);
% xx = A_full \ b;
% 
% nx = nx-1;
% nz = nz-1;
% dvp = reshape(xx(1:nx*nz),[nz nx]);
% 
% figure;imagesc(1:10:10*nx,1:10:10*nz,dvp);
% colorbar; xlabel('x (m)'); ylabel('z (m)'); title('Vp0 update (m/s)');
% 
% dvh = reshape(xx(nx*nz+1:2*nx*nz),[nz nx]);
% 
% figure;imagesc(1:10:10*nx,1:10:10*nz,dvh);
% colorbar; xlabel('x (m)'); ylabel('z (m)'); title('Vh update (m/s)');
% 
% dvn = reshape(xx(2*nx*nz+1:3*nx*nz),[nz nx]);
% 
% figure;imagesc(1:10:10*nx,1:10:10*nz,dvn);
% colorbar; xlabel('x (m)'); ylabel('z (m)'); title('Vn update (m/s)');


% filename='dvp_file';
% fid=fopen(filename,'w');
% fwrite(fid,dvp,'float');
% fclose(fid);

%b_full = sparse([1:b_row]',1,b,1);