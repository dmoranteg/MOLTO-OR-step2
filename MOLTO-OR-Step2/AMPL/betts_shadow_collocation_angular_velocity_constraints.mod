###############################################################
# Low-Thrust Earth Orbit Optimal Control Problem
# Transcription Method using Hermite-Collocation 
# Using Nonlinear Programming IPOPT Solver.
# Implemented by David Morante, 2015
###############################################################
#
reset;
#
# UNITY CONVERTORS 
#
param pi  := 4*atan(1.);
param d2r := pi/180;
param r2d := 180/pi;
#
# CONSTANT PARAMETERS 
#
param g0  := 9.80665;
param Re  := 6378.145e3;
param mu  := 398601e9;
param tc  := sqrt(Re^3/mu);
param ac  := mu/(Re*Re);
param J2  := 1082.6e-06;
#
##############################################################################
#                 EPHEMERIDES INTERPOLATION PARAMETERS                       #
##############################################################################
param s1{1..7};
param s2{1..7};
param s3{1..7};
##############################################################################
#                          MISSION PARAMETERS                                #
##############################################################################
#
# INITIAL ORBIT 
#
param a0   := 3.82;
param e0   := 0.731;
param in0  := 27*d2r;
param w0   := 0;
param Om0  := 99*d2r;
#
# FINAL ORBIT 
#
param af   := 6.6107;
param ef   := 1e-4;
param incf := 1e-4*d2r;
#
# SC PARAMETERS 
#
param T_W  := 4.55e-05; # THRUST TO MASS RATIO
param T    := T_W*g0 /ac;
param m0   := 450; # INITIAL MASS (kg)
param Isp  := 3300;# SPECIFIC IMPULSE (s)
#
##############################################################################
#                              VARIABLES                                     #
##############################################################################
#
# DEFINE SETS
#
param nodes;
set N   := {0 .. nodes};
set N2  := {0 .. nodes-1};
#
# LOAD PARAMETERS FOR PHASES 
#
include model_param.dat
#
# STATES - modified equinoccial elements 
#
var p   {N};
var f   {N};
var g   {N};
var h   {N};
var k   {N};
var t   {N};
var m   {N};
var Thr {N};
#
# CONTROL - inplane(alpha, wrt transversal direction), out-of-plane (beta) 
#
var alpha {N};
var beta  {N};
#
# FINAL TIME 
#
param tfinal;
#
# NUMBER OF REVOLUTIONS 
#
var Lf;
#
# LOAD INITIAL GUESS 
#
include ampl_guess_coll.dat;
#
# DEPENDENT VARIABLE
#
param step := 1 /(nodes); # vector of steps
#
param L0:= (Om0 + w0)/(2*pi);
var L{i in N} = L0*2*pi + i*step*(Lf-L0)*2*pi;
#
var scale = (Lf-L0)*2*pi;
#
# DEFINE ADDITONAL VARIABLES
#
var w{i in N} = 1 + f[i]*cos(L[i]) + g[i]*sin(L[i]);
var r{i in N} = p[i]/w[i];
var s{i in N} = sqrt(1 + h[i]^2 + k[i]^2);
#
##############################################################################
#                            SHADOW CONSTRAINTS                              #
##############################################################################
#
var alpha_aux_2{i in N} = h[i]^2 - k[i]^2;
#
var sun_vector1{i in N} = s1[7] + s1[6]*t[i] +s1[5]*t[i]^2 +s1[4]*t[i]^3+ s1[3]*t[i]^4 +s1[2]*t[i]^5 +s1[1]*t[i]^6;
var sun_vector2{i in N} = s2[7] + s2[6]*t[i] +s2[5]*t[i]^2 +s2[4]*t[i]^3+ s2[3]*t[i]^4 +s2[2]*t[i]^5 +s2[1]*t[i]^6;
var sun_vector3{i in N} = s3[7] + s3[6]*t[i] +s3[5]*t[i]^2 +s3[4]*t[i]^3+ s3[3]*t[i]^4 +s3[2]*t[i]^5 +s3[1]*t[i]^6;
#
var r_vector1{i in N}   = r[i]/s[i]^2 *( cos(L[i])+alpha_aux_2[i]*cos(L[i])+2*h[i]*k[i]*sin(L[i]));
var r_vector2{i in N}   = r[i]/s[i]^2 *( sin(L[i])-alpha_aux_2[i]*sin(L[i])+2*h[i]*k[i]*cos(L[i]));
var r_vector3{i in N}   = r[i]/s[i]^2 * 2*(h[i]*sin(L[i])-k[i]*cos(L[i]));
#
var elem1{i in N}       = sun_vector1[i]*r_vector1[i] +sun_vector2[i]*r_vector2[i]+sun_vector3[i]*r_vector3[i];
var elem2{i in N}       = r_vector1[i]^2 + r_vector2[i]^2 + r_vector3[i]^2;
#
C_elem2{i in N} : r_vector1[i]^2 + r_vector2[i]^2 + r_vector3[i]^2 >=1;
#
var delta{i in N}       = -(elem1[i] + sqrt(elem2[i]-1));
#
var signoid{i in N}    = 1 / (1+exp(delta[i]/0.01));
#
##############################################################################
#                       PERTURBING ACCELERATION                              #
##############################################################################
#
# THRUST PERTURBATION
#
var aT  {i in N} = signoid[i]*Thr[i]*T/m[i];
#
var ar_T{i in N} = aT[i]*cos(beta[i])*sin(alpha[i]);
var ao_T{i in N} = aT[i]*cos(beta[i])*cos(alpha[i]);
var ah_T{i in N} = aT[i]*sin(beta[i]);
#
# PERTURBATIONS DUE TO THE EARTH OBLATENESS
#
var ar_j2{i in N} =  -  3* J2/(2*r[i]^4)*(1 - 12*( h[i]*sin(L[i])-k[i]*cos(L[i]) )^2/s[i]^4);
var ao_j2{i in N} =  - 12* J2/r[i]^4 *(h[i]*sin(L[i]) - k[i]*cos(L[i]))*(h[i]*cos(L[i])+ k[i]*sin(L[i]))/s[i]^4;
var ah_j2{i in N} =  -  6* J2/r[i]^4* (h[i]*sin(L[i]) - k[i]*cos(L[i]))*(1 - h[i]^2 - k[i]^2)/s[i]^4;
#
# TOTAL PERTURBATION 
#
var ar{i in N} = ar_T[i] + ar_j2[i];
var ao{i in N} = ao_T[i] + ao_j2[i];
var ah{i in N} = ah_T[i] + ah_j2[i];
#
##############################################################################
#                           DYNAMICAL EQUATIONS                              #
##############################################################################
#
var A1{i in N}  = sqrt(p[i])/w[i]*(               0*ar[i]                        +2*p[i]*ao[i]                                        +  0*ah[i] );
var A2{i in N}  = sqrt(p[i])/w[i]*(  w[i]*sin(L[i])*ar[i]  + ((w[i]+1)*cos(L[i]) + f[i])*ao[i]     - g[i]*(h[i]*sin(L[i]) -k[i]*cos(L[i]))*ah[i] );
var A3{i in N}  = sqrt(p[i])/w[i]*( -w[i]*cos(L[i])*ar[i]  + ((w[i]+1)*sin(L[i]) + g[i])*ao[i]     + f[i]*(h[i]*sin(L[i]) -k[i]*cos(L[i]))*ah[i] );
var A4{i in N}  = sqrt(p[i])/w[i]*(               0*ar[i]                             +0*ao[i]                       +  s[i]^2*cos(L[i])/2*ah[i] );
var A5{i in N}  = sqrt(p[i])/w[i]*(               0*ar[i]                             +0*ao[i]                       +  s[i]^2*sin(L[i])/2*ah[i] );
var A6{i in N}  = sqrt(p[i])/w[i]*(               0*ar[i]                             +0*ao[i]         +  (h[i]*sin(L[i]) -k[i]*cos(L[i]))*ah[i] );
var A7{i in N}  = -aT[i]*m[i]/(g0*Isp) ;
#
var sigma{i in N} = A6[i] + sqrt(p[i])*(w[i]/p[i])^2 ;
#
var pdot {i in N}  = A1[i]/sigma[i];
var fdot {i in N}  = A2[i]/sigma[i];
var gdot {i in N}  = A3[i]/sigma[i];
var hdot {i in N}  = A4[i]/sigma[i];
var kdot {i in N}  = A5[i]/sigma[i];
var tdot {i in N}  = 1/sigma[i]*tc/(365*24*3600);
var mdot {i in N}  = A7[i]/sigma[i]*tc*ac;
#
##############################################################################
#                       DEPENDENT VARIABLES AT MIDPOINT                      #
##############################################################################
#
var p_mid{i in N2} = 0.5*(p[i]+p[i+1]) + step*scale/8*(pdot[i]-pdot[i+1]);
var f_mid{i in N2} = 0.5*(f[i]+f[i+1]) + step*scale/8*(fdot[i]-fdot[i+1]);
var g_mid{i in N2} = 0.5*(g[i]+g[i+1]) + step*scale/8*(gdot[i]-gdot[i+1]);
var h_mid{i in N2} = 0.5*(h[i]+h[i+1]) + step*scale/8*(hdot[i]-hdot[i+1]);
var k_mid{i in N2} = 0.5*(k[i]+k[i+1]) + step*scale/8*(kdot[i]-kdot[i+1]);
var t_mid{i in N2} = 0.5*(t[i]+t[i+1]) + step*scale/8*(tdot[i]-tdot[i+1]);
var m_mid{i in N2} = 0.5*(m[i]+m[i+1]) + step*scale/8*(mdot[i]-mdot[i+1]);
#
# INDEPENDENT VARIABLE AND CONTROL AT THE MIDPOINT
#
var L_mid    {i in N2} = 0.5*(L[i+1]+L[i]);
var alpha_mid{i in N2} = 0.5*(alpha[i]+alpha[i+1]);
var beta_mid {i in N2} = 0.5*(beta[i+1]+beta[i+1]);
var Thr_mid  {i in N2} = Thr[i];
#
# ADDITIONAL VARIABLES AT MIDPOINT
#
var w_mid{i in N2}  = 1 + f_mid[i]*cos(L_mid[i]) + g_mid[i]*sin(L_mid[i]);
var r_mid{i in N2}  = p_mid[i]/w_mid[i];
var s_mid{i in N2}  = sqrt(1 + h_mid[i]^2 + k_mid[i]^2);
#
##############################################################################
#                    SHADOW CONSTRAINTS AT MIDPOINT                          #
##############################################################################
#
var alpha_aux_2_mid{i in N2} = h_mid[i]^2 - k_mid[i]^2;
#
var sun_vector1_mid{i in N2} = s1[7] + s1[6]*t_mid[i] +s1[5]*t_mid[i]^2 +s1[4]*t_mid[i]^3+ s1[3]*t_mid[i]^4 +s1[2]*t_mid[i]^5 +s1[1]*t_mid[i]^6;
var sun_vector2_mid{i in N2} = s2[7] + s2[6]*t_mid[i] +s2[5]*t_mid[i]^2 +s2[4]*t_mid[i]^3+ s2[3]*t_mid[i]^4 +s2[2]*t_mid[i]^5 +s2[1]*t_mid[i]^6;
var sun_vector3_mid{i in N2} = s3[7] + s3[6]*t_mid[i] +s3[5]*t_mid[i]^2 +s3[4]*t_mid[i]^3+ s3[3]*t_mid[i]^4 +s3[2]*t_mid[i]^5 +s3[1]*t_mid[i]^6;
#
var r_vector1_mid{i in N2}   = r_mid[i]/s_mid[i]^2 *( cos(L_mid[i])+alpha_aux_2_mid[i]*cos(L_mid[i])+2*h_mid[i]*k_mid[i]*sin(L_mid[i]));
var r_vector2_mid{i in N2}   = r_mid[i]/s_mid[i]^2 *( sin(L_mid[i])-alpha_aux_2_mid[i]*sin(L_mid[i])+2*h_mid[i]*k_mid[i]*cos(L_mid[i]));
var r_vector3_mid{i in N2}   = r_mid[i]/s_mid[i]^2 * 2*(h_mid[i]*sin(L_mid[i])-k_mid[i]*cos(L_mid[i]));
#
var elem1_mid{i in N2}       = sun_vector1_mid[i]*r_vector1_mid[i] +sun_vector2_mid[i]*r_vector2_mid[i]+sun_vector3_mid[i]*r_vector3_mid[i];
var elem2_mid{i in N2}       = r_vector1_mid[i]^2 + r_vector2_mid[i]^2 + r_vector3_mid[i]^2;
#
C_elem2_mid{i in N2} : r_vector1_mid[i]^2 + r_vector2_mid[i]^2 + r_vector3_mid[i]^2 >=1;
#
var delta_mid{i in N2}       = - (elem1_mid[i] + sqrt(elem2_mid[i]-1));
#
var signoid_mid{i in N2}     = 1 / (1+exp(delta_mid[i]/0.01));
#
##############################################################################
#              PERTURBING ACCELERATION AT MIDPOINT                           #
##############################################################################
#
# THRUST PERTURBATION
#
var aT_mid{i in N2}    = signoid_mid[i]*Thr_mid[i]*T/m_mid[i];
#
var ar_T_mid{i in N2}  = aT_mid[i]*cos(beta_mid[i])*sin(alpha_mid[i]);
var ao_T_mid{i in N2}  = aT_mid[i]*cos(beta_mid[i])*cos(alpha_mid[i]);
var ah_T_mid{i in N2}  = aT_mid[i]*sin(beta_mid[i]);
#
# PERTURBATIONS DUE TO THE EARTH OBLATENESS
#
var ar_j2_mid{i in N2} =  -  3* J2/(2*r_mid[i]^4)*( 1 - 12*( h_mid[i]*sin(L_mid[i])-k_mid[i]*cos(L_mid[i]) )^2/s_mid[i]^4);
var ao_j2_mid{i in N2} =  - 12* J2/r_mid[i]^4 *( h_mid[i]*sin(L_mid[i])-k_mid[i]*cos(L_mid[i]) )*( h_mid[i]*cos(L_mid[i])+ k_mid[i]*sin(L_mid[i]) )/s_mid[i]^4;
var ah_j2_mid{i in N2} =  -  6* J2/r_mid[i]^4* ( h_mid[i]*sin(L_mid[i])-k_mid[i]*cos(L_mid[i]) )*(1 - h_mid[i]^2 - k_mid[i]^2)/s_mid[i]^4;
#
# TOTAL PERTURBATION 
#
var ar_mid{i in N2} = ar_T_mid[i] + ar_j2_mid[i];
var ao_mid{i in N2} = ao_T_mid[i] + ao_j2_mid[i];
var ah_mid{i in N2} = ah_T_mid[i] + ah_j2_mid[i];
#
##############################################################################
#                 DYNAMICAL EQUATIONS AT MIDPOINT                            #
##############################################################################
#
var A1_mid{i in N2} = sqrt(p_mid[i])/w_mid[i]*(               0*ar_mid[i]                        +2*p_mid[i]*ao_mid[i]                                     +  0*ah_mid[i] );
var A2_mid{i in N2} = sqrt(p_mid[i])/w_mid[i]*(  w_mid[i]*sin(L_mid[i])*ar_mid[i]  + ((w_mid[i]+1)*cos(L_mid[i]) + f_mid[i])*ao_mid[i]  - g_mid[i]*(h_mid[i]*sin(L_mid[i]) -k_mid[i]*cos(L_mid[i]))*ah_mid[i] );
var A3_mid{i in N2} = sqrt(p_mid[i])/w_mid[i]*( -w_mid[i]*cos(L_mid[i])*ar_mid[i]  + ((w_mid[i]+1)*sin(L_mid[i]) + g_mid[i])*ao_mid[i]    + f_mid[i]*(h_mid[i]*sin(L_mid[i]) -k_mid[i]*cos(L_mid[i]))*ah_mid[i] );
var A4_mid{i in N2} = sqrt(p_mid[i])/w_mid[i]*(               0*ar_mid[i]                             +0*ao_mid[i]                    +  s_mid[i]^2*cos(L_mid[i])/2*ah_mid[i] );
var A5_mid{i in N2} = sqrt(p_mid[i])/w_mid[i]*(               0*ar_mid[i]                             +0*ao_mid[i]                   +  s_mid[i]^2*sin(L_mid[i])/2*ah_mid[i] );
var A6_mid{i in N2} = sqrt(p_mid[i])/w_mid[i]*(               0*ar_mid[i]                             +0*ao_mid[i]     +  (h_mid[i]*sin(L_mid[i]) -k_mid[i]*cos(L_mid[i]))*ah_mid[i] );
var A7_mid{i in N2} = - aT_mid[i]*m_mid[i]/(g0*Isp) ;
#
var sigma_mid{i in N2} = A6_mid[i] + sqrt(p_mid[i])*(w_mid[i]/p_mid[i])^2 ;
#
var pdot_mid {i in N2} = A1_mid[i]/sigma_mid[i];
var fdot_mid {i in N2} = A2_mid[i]/sigma_mid[i];
var gdot_mid {i in N2} = A3_mid[i]/sigma_mid[i];
var hdot_mid {i in N2} = A4_mid[i]/sigma_mid[i];
var kdot_mid {i in N2} = A5_mid[i]/sigma_mid[i];
var tdot_mid {i in N2} = 1/sigma_mid[i]*tc/(365*24*3600);
var mdot_mid {i in N2} = A7_mid[i]/sigma_mid[i]*tc*ac;
#
##################################################################################################
##################################################################################################
#                            HERMITE-SIMPSON CONSTRAINTS                                         #
##################################################################################################
##################################################################################################
#
state_1 {i in N2}: ( p[i]-p[i+1] +step*scale/6*( pdot[i]+pdot[i+1] + 4*pdot_mid[i]) )  = 0;   
state_2 {i in N2}: ( f[i]-f[i+1] +step*scale/6*( fdot[i]+fdot[i+1] + 4*fdot_mid[i]) )  = 0;        
state_3 {i in N2}: ( g[i]-g[i+1] +step*scale/6*( gdot[i]+gdot[i+1] + 4*gdot_mid[i]) )  = 0;          
state_4 {i in N2}: ( h[i]-h[i+1] +step*scale/6*( hdot[i]+hdot[i+1] + 4*hdot_mid[i]) )  = 0;        
state_5 {i in N2}: ( k[i]-k[i+1] +step*scale/6*( kdot[i]+kdot[i+1] + 4*kdot_mid[i]) )  = 0;      
state_6 {i in N2}: ( t[i]-t[i+1] +step*scale/6*( tdot[i]+tdot[i+1] + 4*tdot_mid[i]) )  = 0;
state_7 {i in N2}: ( m[i]-m[i+1] +step*scale/6*( mdot[i]+mdot[i+1] + 4*mdot_mid[i]) )  = 0;
##################################################################################################
#
##############################################################################
#                 INEQUALITY CONSTRAINTS                                     #
##############################################################################
#
# CONTROL
#
ALPHA_C1{i in N}: -pi   <= alpha[i] <= pi;
BETA_C1 {i in N}: -pi/2 <= beta[i]  <= pi/2;
THR_C1{i in N}:       0 <= Thr[i] <= 1;
#
# STATES
#
p_C1{i in N}:  0.5  <= p[i] <= 10 ;
f_C1{i in N}: -2 <= f[i] <= 2;
g_C1{i in N}: -2  <= g[i] <= 2;
h_C1{i in N}: -2 <= h[i] <= 2;
k_C1{i in N}: -2  <= k[i] <= 2;
t_C1{i in N}:  0    <= t[i] <= tfinal;
m_C1{i in N}:  0.1  <= m[i] <= 1; 
#
p_C1_mid{i in N2}:  0.5  <= p_mid[i] <= 10 ;
f_C1_mid{i in N2}: -2    <= f_mid[i] <= 2;
g_C1_mid{i in N2}: -2    <= g_mid[i] <= 2;
h_C1_mid{i in N2}: -2    <= h_mid[i] <= 2;
k_C1_mid{i in N2}: -2    <= k_mid[i] <= 2;
t_C1_mid{i in N2}:  0    <= t_mid[i] <= tfinal;
m_C1_mid{i in N2}:  0.1  <= m_mid[i] <= 1; 
#
# ROTATIONAL CONSTRAINTS
#
#var dot_product{i in N2} = sin(alpha[i+1])*sin(alpha[i])*cos(beta[i])*cos(beta[i+1]) + cos(alpha[i+1])*cos(alpha[i])*cos(beta[i])*cos(beta[i+1]) + sin(beta[i])*sin(beta[i+1]);
#
#var alpha_rate{i in N2}  = acos(dot_product[i])/(t[i+1]-t[i]);
#
#var alpha_rate2{i in N2} = alpha_rate[i]*180/pi/(365*24);
#
#rot_cons{i in N2}: 0 <= alpha_rate2[i]/1000000 <= 10/1000000;
#
##############################################################################
#                 EQUALITY CONSTRAINTS                                      #
##############################################################################
#
# INITIAL ORBIT
#
p_0:            p[0] = a0*(1-e0^2);
f_0:	        f[0] = e0*cos(w0+Om0);
g_0:         	g[0] = e0*sin(w0+Om0);
h_0:         	h[0] = tan(in0/2)*cos(Om0);
k_0:         	k[0] = tan(in0/2)*sin(Om0);
t_0:         	t[0] = 0;
m_0:         	m[0] = 1;
#
# FINAL ORBIT
#
a_f:         	af*(1-ef^2)   <= p[nodes]     <= af*(1-ef^2);
e_f:            ef^2          <=f[nodes]^2+g[nodes]^2<=ef^2;
i_f:            tan(incf/2)^2 <= h[nodes]^2+k[nodes]^2<=tan(incf/2)^2;
#
############################################################################
#                         COST FUNCTION PROBLEM                            #
############################################################################

minimize final_time: t[nodes]*365;

#############################################################################
#                            CALLING SOLVER                                 #
#############################################################################
option solver ipopt;
option presolve 0;
option ipopt_options "max_iter 10000 halt_on_ampl_error yes mu_init 1e-5 output_file ipoptlog.txt";
solve;
###############################################################
#                        PRINT OUTPUT                         #
###############################################################

printf {i in N}: "%17.15f %17.15f %17.15f %17.15f %17.15f %17.15f %17.15f\n",
 t[i], p[i], f[i], g[i] , h[i], k[i], m[i]> output.out;
 
printf {i in N}: "%17.15f  \n",
alpha[i] > alpha.out;

printf {i in N}: "%17.15f  \n",
Thr[i] > T.out;
 
printf {i in N}: "%17.15f  \n",
beta[i] > beta.out;

printf{i in N}:"%17.15f  \n",
aT[i] > at.out;

printf{i in N}:"%19.15f  \n",
L[i] > L.out;

printf{i in N}:"%19.15f  \n",
signoid[i] > signoid.out;


printf{i in N}:"%17.15f  \n",
T_magnitude[i] > T_magnitude.out;

