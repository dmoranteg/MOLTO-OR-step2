
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Obtain the guess for the Ampl shadow Betts problem %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
testcase = num2str(tfinal);
%
%% LOAD THE INPUT DATA STRUCTURE FROM qlaw.
%
% %
trj        = load(initialGuess);
%
Re = 6.378140000e+03;
%
% STATE
%
t        = (trj(:,1)-trj(1,1))/(3600*24*365);
p        = trj(:,2)/Re;
f        = trj(:,3);
g        = trj(:,4);
h        = trj(:,5);
k        = trj(:,6);
L        = trj(:,7);
m        = trj(:,8)/trj(1,8);%
% CONTROL
%
alpha    = trj(:,9);
beta     = trj(:,10);
T = ones(size(alpha));
%
n = ceil(L(end)/(2*pi)*n);

N = n+1;
%
% INTERPOLATE DATA FOR AT THE REQUIRED NODES.
%
Linterp = linspace(L(1),L(end),N)';
%
p_interp = interp1(L,p,Linterp,'pchip','extrap');
f_interp = interp1(L,f,Linterp,'pchip','extrap');
g_interp = interp1(L,g,Linterp,'pchip','extrap');
h_interp = interp1(L,h,Linterp,'pchip','extrap');
k_interp = interp1(L,k,Linterp,'pchip','extrap');
t_interp = interp1(L,t,Linterp,'pchip','extrap');
m_interp = interp1(L,m,Linterp,'pchip','extrap');
T_interp = interp1(L,double(T),Linterp,'linear','extrap');

%signoid =  ecpse<0;
%alpha = alpha(signoid>0 & T>0);
%beta  = beta(signoid>0 & T>0);
%L  = L(signoid>0 & T>0);%

c_alpha_interp = interp1(L,cos(alpha),Linterp,'linear','extrap');
s_alpha_interp = interp1(L,sin(alpha),Linterp,'linear','extrap');
alpha_interp   = atan2(s_alpha_interp,c_alpha_interp);

c_beta_interp  = interp1(L,cos(beta),Linterp,'linear','extrap');
s_beta_interp  = interp1(L,sin(beta),Linterp,'linear','extrap');
beta_interp    = atan2(s_beta_interp,c_beta_interp);

%T_interp = T_interp > 0.5;

% figure
% plot(Linterp,t_interp)
% stop

fprintf(fid2,'let nodes:= %i;\n',n);
fprintf(fid ,'let tfinal:= %f;\n', tfinal/365);

%
for i = 1:N
    %p
    fprintf(fid,'let p[%i] := %15.13f;\n',i-1,p_interp(i));
    %f
    fprintf(fid,'let f[%i] := %15.13f;\n',i-1,f_interp(i));
    %g
    fprintf(fid,'let g[%i] := %15.13f;\n',i-1,g_interp(i));
    %h
    fprintf(fid,'let h[%i] := %15.13f;\n',i-1,h_interp(i));
    %k
    fprintf(fid,'let k[%i] := %15.13f;\n',i-1,k_interp(i));
    %t
    fprintf(fid,'let t[%i] := %15.13f;\n',i-1,t_interp(i));
    %m
    fprintf(fid,'let m[%i] := %15.13f;\n',i-1,m_interp(i));
    %alpha
    fprintf(fid,'let alpha[%i] := %15.13f;\n',i-1,alpha_interp(i));
    %beta
    fprintf(fid,'let beta[%i] := %15.13f;\n',i-1 ,beta_interp(i));
    %T_magnitude
    fprintf(fid,'let Thr[%i] := %15.13f;\n',i-1 ,T_interp(i));
    %
end

fprintf(fid,'let Lf := %15.13f;\n',(L(end))/(2*pi)-2);

%% PROVIDE PARAMETERS FOR THE THE INTERPOLATION OF THE EARTH-SUN VECTOR

[A]= load('ephemerid2.mat');
RAN_sun  = A.ephemerid(:,1)*pi/180;
i_sun    = A.ephemerid(:,2)*pi/180;

r1 = cos(i_sun).*cos(RAN_sun);
r2 = cos(i_sun).*sin(RAN_sun);
r3 = sin(i_sun);

t = linspace (0,1,numel( A.ephemerid(:,2)));
nu = 6;
[s1,p1] = polyfit(t',r1,nu);
[s2,p2] = polyfit(t',r2,nu);
[s3,p3] = polyfit(t',r3,nu);

for i = 1:nu+1
    fprintf(fid,'let s1[%i] := %15.13f;\n',i ,s1(i));
    fprintf(fid,'let s2[%i] := %15.13f;\n',i ,s2(i));
    fprintf(fid,'let s3[%i] := %15.13f;\n',i ,s3(i));
end

fprintf(fid3,['cd ','./',dir,';\n ']);

fprintf(fid3,'model "../AMPL/betts_shadow_collocation_cone.mod";');


fclose(fid);
fclose(fid2);
fclose(fid3);


