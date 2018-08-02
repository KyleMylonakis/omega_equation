clear;
clc;
%load('LTER04_mid_channel_3D_solver_test_data_w_zeta.mat');
%load('LTER04_LongSouth_grid_3D_solver_input_data.mat');
%load('LTER04_mid_grid_3D_solver_input_data.mat');
%load('LTER04_MaxNS_grid_3D_solver_input_data.mat');
load('LTER04_MaxNS_grid_3D_solver_input_data_for_neumann.mat');
interp_type  = 'spline';
global divQgeo_pruned_interp N_sqrd_pruned f zeta;
global g_Qx g_Qy;
divQgeo_pruned_interp = fillmissing(fillmissing(divQgeo_pruned,interp_type),interp_type,2);
g_Qx = fillmissing(fillmissing(Qx_pruned,interp_type),interp_type,2);
g_Qy = fillmissing(fillmissing(Qy_pruned,interp_type),interp_type,2);

global a b c
global g_dx g_dy g_dz
global g_H

g_H = H;
%g_dx = dx/1000;
%g_dy = dy/1000;
g_dx = 2000;
g_dy = 2000;
dz = 2; %dz is 2 meters
g_dz = dz;
a = length(divQgeo_pruned(:,1,1));
b = length(divQgeo_pruned(1,:,1));
c = length(divQgeo_pruned(1,1,:));

RHS = zeros((a-2)*(b-2)*(c-2),1);

for k = 2:c-1
    for j = 2:b-1
        for i = 2:a-1
            RHS((k-2)*(b-2)*(a-2) + (j-2)*(a-2) + i-1) = RHS_f(i,j,k);
        end
    end
end


%fft_z(RHS);
%N_sqrd_pruned = N_sqrd_pruned/(1000*1000);
%w = gmres(@A,RHS,10,10^(-6),10,@fft_z);
%w = gmres(@A,RHS,[],[],100,@ifft_z);
%w = gmres(@A,RHS,[],10^(-3),(a-2)*(b-2)*(c-2),@ifft_z);
w = gmres(@A,RHS,[],[],100);
%w = gmres(@A,RHS,[],[],(a-2)*(b-2)*(c-2));

%norm(RHS - ifft_z(fft_z(RHS)))


omega = zeros(a-2,b-2,c-2);
omega = reshape(w,[a-2,b-2,c-2]);
omega_ft = fftn(omega);
%omega_ft(1,1,1) = 0;
%omega = ifftn(omega_ft);
%for k = 2:c-1
%    for j = 2:b-1
%        for i = 2:a-1
%            omega(i,j,k) = w( (k-2)*(b-2)*(a-2) + (j-2)*(a-2) + i-1 );
%        end
%    end
%end

plt = contourf(mean(omega,3),20);
colorbar()
title('Vertically Averaged Vertical Velocity')

save('omega_eqn_zero_neumann_soln');




function output = fft_z(w)
    global a b c;
    omega = reshape(w,[a-2,b-2,c-2]);
    odd_omega = zeros(a-2,b-2,2*c -2);
    odd_omega(:,:,2:c-1) = omega;
    for kk = 1:c-2
        odd_omega(:,:,c+kk) = -omega(:,:, c-kk-1 );
    end
    %for i = 1:a-1
    %    for j = 1:b-1
    %        for k = 1:c-1
    %            omega(i,j,k) = w(I(i-1,j-1,k-1));
    %        end
    %    end
    %end
    odd_omega = -0.5*imag(fft(omega,2*c - 2,3));
    omega = odd_omega(:,:,2:c-1);
    output = omega(:);
    
end

function output = ifft_z(w)
    global a b c;
    global g_dz;
    omega = reshape(w,[a-2,b-2,c-2]);
    odd_omega = zeros(a-2,b-2,2*c -2);
    odd_omega(:,:,2:c-1) = omega;
    for kk = 1:c-2
        odd_omega(:,:,c+kk) = -omega(:,:, c-kk-1 );
    end
    %for i = 1:a-1
    %    for j = 1:b-1
    %        for k = 1:c-1
    %            omega(i,j,k) = w(I(i-1,j-1,k-1));
    %        end
    %    end
    %end
    odd_omega = 2*imag(ifft(omega,2*c - 2,3));
    omega = odd_omega(:,:,2:c-1);
    for kk = 1:c-2
        omega(:,:,kk) = omega(:,:,kk)./(-4*(sin(pi*(kk)/(2*(c-2) + 1))).^2/(g_dz^2));
    end
    output = omega(:);
    
end




function output = f1_f(i,j,k)
    global N_sqrd_pruned
    output = N_sqrd_pruned(i,j,k);
    %output = 1;
end

function output = f2_f(i,j,k)
    global f zeta g_H
    output = f.*(f+(zeta(i,j)./g_H));
    %output = 1;
end



function output = RHS_f(i,j,k)
    global divQgeo_pruned_interp
    output = 2 * divQgeo_pruned_interp(i,j,k);
end

function output = I(i,j,k,a,b)
    % a = length(x), b = length(y)
    output = (k-1)*(a-2)*(b-2) + (j-1)*(a-2) + i;
end

function output = J(x)
    global a b c;
    output = zeros(3,1);
    output(1) = mod(x-1,a-2) + 1;
    output(2) = mod((x - output(1))/(a-2),b-2) + 1;
    output(3) = (x - output(1) - (output(2)-1)*(a-2))/((a-2)*(b-2)) + 1;
end

function output = A(w)
    %global x y z;
    global g_dx g_dy g_dz;
    global g_Qx g_Qy;
    omega = zeros(length(w),1);
    output = zeros(length(w),1);
    %length(output)
    global a b c;
    for k = 2:c-1
        for j = 2:b-1
            for i = 2:a-1
                temp = 0;
                if i ~= 2
                    temp = temp...
                        + f1_f((i-1),(j),(k))*w(I(i-2,j-1,k-1,a,b));
                else
                    temp = temp...
                          + f1_f((i+1),(j),(k))*w(I(i,j-1,k-1,a,b))...
                          - 4*g_dy*g_Qy(i,j,k);
                end
                if i ~=a-1
                    temp = temp...
                        + f1_f((i+1),(j),(k))*w(I(i,j-1,k-1,a,b));
                else
                    temp = temp...
                    + f1_f((i-1),(j),(k))*w(I(i-2,j-1,k-1,a,b))...
                    + 4*g_dy*g_Qy(i,j,k);
                    
                end
                temp = temp ...
                    -2 * f1_f((i),(j),(k))*w(I(i-1,j-1,k-1,a,b));
                
                output(I(i-1,j-1,k-1,a,b)) = output(I(i-1,j-1,k-1,a,b)) ...
                    + temp/g_dy^2;
                
                temp = 0;
                if j ~= 2
                    temp = temp...
                        + f1_f((i),(j-1),(k))*w(I(i-1,j-2,k-1,a,b));
                else
                    temp = temp...
                        + f1_f((i),(j+1),(k))*w(I(i-1,j,k-1,a,b))...
                        - 4*g_dx*g_Qx(i,j,k);
                end
                
                if j ~=b-1
                    temp = temp...
                        + f1_f((i),(j+1),(k))*w(I(i-1,j,k-1,a,b));
                else
                    temp = temp...
                         + f1_f((i),(j-1),(k))*w(I(i-1,j-2,k-1,a,b))...
                         + 4*g_dx*g_Qx(i,j,k);
                end
                
                temp = temp -2 * f1_f((i),(j),(k))*w(I(i-1,j-1,k-1,a,b));
                output(I(i-1,j-1,k-1,a,b)) = output(I(i-1,j-1,k-1,a,b))...
                    + temp/g_dx^2;
                
                temp = 0;
                
                if k ~= 2
                    temp = temp...
                        + f2_f((i),(j),(k))*w(I(i-1,j-1,k-2,a,b));
                else
                    % Comment out for Dirichlet, uncomment for Neumann
                   %  temp = temp...
                    %    + f2_f((i),(j),(k))*w(I(i-1,j-1,k,a,b));
                end
                if k ~=c-1
                    temp = temp...
                        + f2_f((i),(j),(k))*w(I(i-1,j-1,k,a,b));
                else
                    temp = temp...
                         + f2_f((i),(j),(k))*w(I(i-1,j-1,k-2,a,b)); 
                end
                    
                temp = temp...
                    -2 * f2_f((i),(j),(k))*w(I(i-1,j-1,k-1,a,b));

                output(I(i-1,j-1,k-1,a,b)) = output(I(i-1,j-1,k-1,a,b))...
                    + temp/g_dz^2;
            end
        end
    end
    
end
