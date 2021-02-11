
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
Re = 6.378140000e+06;
%
% STATE
%
t        = (trj(:,1)-trj(1,1))/365;
a        = trj(:,8)/Re;
e        = trj(:,9);
in       = trj(:,10);
Om       = trj(:,11);
w        = trj(:,12);
v        = trj(:,13);
m        = trj(:,21)/trj(1,21);%
% CONTROL
%
alpha    = trj(:,14);
beta     = trj(:,15);
T        = trj(:,16)/max(trj(:,16));
T(end)   = 1;
ecpse    = trj(:,19);
%
% CONVERT STATE TO MODIFIED EQUINOCCIAL ELEMENTS
%
mee  = oe2mee([ a , e , in , Om , w , v ]);
%
p = mee(:,1); f = mee(:,2); g = mee(:,3) ; h = mee(:,4) ; k = mee(:,5); L = mee(:,6);
%
L = wrapTo2Pi(L);
%
for i = 2:numel(L)
    while (L(i)-L(i-1)) < 0
        L(i) = L(i) + 2*pi;
    end
end
%
n = ceil(L(end)/(2*pi)*n);

N = n+1;
%
% INTERPOLATE DATA FOR AT THE REQUIRED NODES.
%
Linterp = linspace(L(1),L(end),N)';
%
p_interp = interp1(L,p,Linterp,'linear','extrap');
f_interp = interp1(L,f,Linterp,'linear','extrap');
g_interp = interp1(L,g,Linterp,'linear','extrap');
h_interp = interp1(L,h,Linterp,'linear','extrap');
k_interp = interp1(L,k,Linterp,'linear','extrap');
t_interp = interp1(L,t,Linterp,'linear','extrap');
m_interp = interp1(L,m,Linterp,'linear','extrap');
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

fprintf(fid,'let Lf := %15.13f;\n',(L(end))/(2*pi));

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


