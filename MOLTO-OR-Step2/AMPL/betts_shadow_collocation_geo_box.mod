# Betts optimization problem
# Optimal Control using collocation 
# Using Nonlinear Programming; IPOPT Solver.
# Implemented by David Morante, 2015

###############################################################
reset;
################# UNITY CONVERTORS ######################
param pi  := 4*atan(1.);
param d2r := pi/180;
param r2d := 180/pi;
################# CONSTANT PARAMETERS ######################
param g0  := 9.80665;
param Re  := 6378.145e3;
param mu  := 398601e9;
param tc  := sqrt(Re^3/mu);
param ac  := mu/(Re*Re);
param J2  := 1082.6e-06;
########################## EPHEMERIDES INTERPOLATION PARAMETERS ######################
param s1{1..7};
param s2{1..7};
param s3{1..7};
############################### BOUNDARY VALUES ######################################
param af    := 6.6107;
#param ef   := 1e-4;
param ef   := 0.00169;
#param ef    := 1e-4;
param incf := 0.2*d2r;
#param incf := 0.001778768484826; # phi
#param incf  := 1e-4*d2r;
param T_W   := 1.515006944995779e-05; 
param T     := T_W*g0 /ac;
param m0    := 3500;
param Isp   := 1600;
#################################### SETS #############################################

param nodes;
param nphases;
param phase_up;
param phase_down;

set iset  := {0 ..nphases-1};
set jset  := {0 .. nodes-1};
set jset2 := {0 .. nodes-2};

######################  PARAMETERS FOR PHASES ##########################

include model_param.dat

######################   INDEPENDENT VARIABLES ##########################
# states

param L0:= (Om0 + w0 +v0);

var p {iset,jset};
var f {iset,jset};
var g {iset,jset};
var h {iset,jset};
var k {iset,jset};
var t {iset,jset};
var m {iset,jset};
var T_magnitude{iset,jset};
var L_Om{iset,jset};

### CONTROL #######

var alpha{iset,jset};
var beta {iset,jset};

###################   DEPENDENT VARIABLE    ###########################
#
include ampl_guess_coll.dat;
#
param step = 1/(nodes-1); # vector of steps
#
var L{i in iset,j in jset} = L_Om[i,0] + j*step*(L_Om[i,nodes-1]-L_Om[i,0]);
#
var scale {i in iset} = (L_Om[i,nodes-1]-L_Om[i,0]);
#
c_scale {i in iset}: scale[i] >=0;
#
### DEFINE ADDITONAL VARIABLES ###
#
var w {i in iset,j in jset} = 1 + f[i,j]*cos(L[i,j]) + g[i,j]*sin(L[i,j]);
var r {i in iset,j in jset} = p[i,j]/w[i,j];
var s {i in iset,j in jset} = sqrt(1 + h[i,j]^2 + k[i,j]^2);
#
##############################################################################
#                             SHADOW CONSTRAINTS                             #
##############################################################################
var alpha_aux_2{i in iset,j in jset} = h[i,j]^2 - k[i,j]^2;
#
var sun_vector1{i in iset,j in jset} = s1[7] + s1[6]*t[i,j] +s1[5]*t[i,j]^2 +s1[4]*t[i,j]^3+ s1[3]*t[i,j]^4 +s1[2]*t[i,j]^5 +s1[1]*t[i,j]^6;
var sun_vector2{i in iset,j in jset} = s2[7] + s2[6]*t[i,j] +s2[5]*t[i,j]^2 +s2[4]*t[i,j]^3+ s2[3]*t[i,j]^4 +s2[2]*t[i,j]^5 +s2[1]*t[i,j]^6;
var sun_vector3{i in iset,j in jset} = s3[7] + s3[6]*t[i,j] +s3[5]*t[i,j]^2 +s3[4]*t[i,j]^3+ s3[3]*t[i,j]^4 +s3[2]*t[i,j]^5 +s3[1]*t[i,j]^6;
#
var r_vector1{i in iset,j in jset}   = r[i,j]/s[i,j]^2 *( cos(L[i,j])+alpha_aux_2[i,j]*cos(L[i,j])+2*h[i,j]*k[i,j]*sin(L[i,j]));
var r_vector2{i in iset,j in jset}   = r[i,j]/s[i,j]^2 *( sin(L[i,j])-alpha_aux_2[i,j]*sin(L[i,j])+2*h[i,j]*k[i,j]*cos(L[i,j]));
var r_vector3{i in iset,j in jset}   = r[i,j]/s[i,j]^2 * 2*(h[i,j]*sin(L[i,j])-k[i,j]*cos(L[i,j]));
#
var elem1{i in iset,j in jset}       = sun_vector1[i,j]*r_vector1[i,j] +sun_vector2[i,j]*r_vector2[i,j]+sun_vector3[i,j]*r_vector3[i,j];
var elem2{i in iset,j in jset}       = r_vector1[i,j]^2 + r_vector2[i,j]^2 + r_vector3[i,j]^2;

#C_elem2{i in iset,j in jset} : r_vector1[i,j]^2 + r_vector2[i,j]^2 + r_vector3[i,j]^2 >=1;

var delta{i in iset,j in jset}       = -(elem1[i,j] + sqrt(elem2[i,j]-1));
#
#var signoid{i in iset,j in jset}   = 1 / (1+exp(delta[i,j]/0.01));
param signoid{i in iset,j in jset}:= 1;
#
########## COMPUTE THRUST ACCELERATION #######################
#
var aT  {i in iset,j in jset} = signoid[i,j]*T_magnitude[i,j]*T/m[i,j];
#
var ar_T{i in iset,j in jset} = aT[i,j]*cos(beta[i,j])*sin(alpha[i,j]);
var ao_T{i in iset,j in jset} = aT[i,j]*cos(beta[i,j])*cos(alpha[i,j]);
var ah_T{i in iset,j in jset} = aT[i,j]*sin(beta[i,j]);

### COMPUTE THE PERTURBATIONS DUE TO THE EARTH OBLATENESS ####
var ar_j2{i in iset,j in jset}=  -  3* J2/(2*r[i,j]^4)*(1 - 12*( h[i,j]*sin(L[i,j])-k[i,j]*cos(L[i,j]) )^2/s[i,j]^4);
var ao_j2{i in iset,j in jset} =  - 12* J2/r[i,j]^4 *(h[i,j]*sin(L[i,j]) - k[i,j]*cos(L[i,j]))*(h[i,j]*cos(L[i,j])+ k[i,j]*sin(L[i,j]))/s[i,j]^4;
var ah_j2{i in iset,j in jset} =  -  6* J2/r[i,j]^4* (h[i,j]*sin(L[i,j]) - k[i,j]*cos(L[i,j]))*(1 - h[i,j]^2 - k[i,j]^2)/s[i,j]^4;

### COMPUTE THE WHOLE PERTURBATION ######
var ar{i in iset,j in jset} = ar_T[i,j] + ar_j2[i,j];
var ao{i in iset,j in jset} = ao_T[i,j] + ao_j2[i,j];
var ah{i in iset,j in jset} = ah_T[i,j] + ah_j2[i,j];

### COMPUTE THE DYNAMICAL EQUATIONS #####

var A1{i in iset,j in jset}  = sqrt(p[i,j])/w[i,j]*(               0*ar[i,j]                        +2*p[i,j]*ao[i,j]                                        +  0*ah[i,j] );
var A2{i in iset,j in jset}  = sqrt(p[i,j])/w[i,j]*(  w[i,j]*sin(L[i,j])*ar[i,j]  + ((w[i,j]+1)*cos(L[i,j]) + f[i,j])*ao[i,j]     - g[i,j]*(h[i,j]*sin(L[i,j]) -k[i,j]*cos(L[i,j]))*ah[i,j] );
var A3{i in iset,j in jset}  = sqrt(p[i,j])/w[i,j]*( -w[i,j]*cos(L[i,j])*ar[i,j]  + ((w[i,j]+1)*sin(L[i,j]) + g[i,j])*ao[i,j]     + f[i,j]*(h[i,j]*sin(L[i,j]) -k[i,j]*cos(L[i,j]))*ah[i,j] );
var A4{i in iset,j in jset}  = sqrt(p[i,j])/w[i,j]*(               0*ar[i,j]                             +0*ao[i,j]                       +  s[i,j]^2*cos(L[i,j])/2*ah[i,j] );
var A5{i in iset,j in jset}  = sqrt(p[i,j])/w[i,j]*(               0*ar[i,j]                             +0*ao[i,j]                       +  s[i,j]^2*sin(L[i,j])/2*ah[i,j] );
var A6{i in iset,j in jset}  = sqrt(p[i,j])/w[i,j]*(               0*ar[i,j]                             +0*ao[i,j]         +  (h[i,j]*sin(L[i,j]) -k[i,j]*cos(L[i,j]))*ah[i,j] );
var A7{i in iset,j in jset}  = -aT[i,j]*m[i,j]/(g0*Isp) ;

var sigma{i in iset,j in jset} = A6[i,j] + sqrt(p[i,j])*(w[i,j]/p[i,j])^2 ;

############### OPTMAL CONTROL PROBLEM #######################
# Cost function

minimize tiempo_final: t[nphases-1,nodes-1]*365;
#minimize tiempo_final: 1;
#minimize masa_final: -m[nphases-1,nodes-1]*0.1;

##############################################################
 
var pdot {i in iset,j in jset}  = A1[i,j]/sigma[i,j];
var fdot {i in iset,j in jset}  = A2[i,j]/sigma[i,j];
var gdot {i in iset,j in jset}  = A3[i,j]/sigma[i,j];
var hdot {i in iset,j in jset}  = A4[i,j]/sigma[i,j];
var kdot {i in iset,j in jset}  = A5[i,j]/sigma[i,j];
var tdot {i in iset,j in jset}  = 1/sigma[i,j]*tc/(365*24*3600);
var mdot {i in iset,j in jset}  = A7[i,j]/sigma[i,j]*tc*ac;

### COMPUTE DEPENDENT VARIABLES AT THE MIDPOINT ####

var p_mid{i in iset,j in jset2} = 0.5*(p[i,j]+p[i,j+1]) + step*scale[i]/8*(pdot[i,j]-pdot[i,j+1]);
var f_mid{i in iset,j in jset2} = 0.5*(f[i,j]+f[i,j+1]) + step*scale[i]/8*(fdot[i,j]-fdot[i,j+1]);
var g_mid{i in iset,j in jset2} = 0.5*(g[i,j]+g[i,j+1]) + step*scale[i]/8*(gdot[i,j]-gdot[i,j+1]);
var h_mid{i in iset,j in jset2} = 0.5*(h[i,j]+h[i,j+1]) + step*scale[i]/8*(hdot[i,j]-hdot[i,j+1]);
var k_mid{i in iset,j in jset2} = 0.5*(k[i,j]+k[i,j+1]) + step*scale[i]/8*(kdot[i,j]-kdot[i,j+1]);
var t_mid{i in iset,j in jset2} = 0.5*(t[i,j]+t[i,j+1]) + step*scale[i]/8*(tdot[i,j]-tdot[i,j+1]);
var m_mid{i in iset,j in jset2} = 0.5*(m[i,j]+m[i,j+1]) + step*scale[i]/8*(mdot[i,j]-mdot[i,j+1]);

### COMPUTE INDEPENDENT VARIABLE AND CONTROL AT THE MIDPOINT ####
var L_mid    {i in iset,j in jset2} = 0.5*(L[i,j+1]+L[i,j]);
var alpha_mid{i in iset,j in jset2} = alpha[i,j];
var beta_mid {i in iset,j in jset2} = 0.5*(beta[i,j]+beta[i,j+1]) ;
var T_magnitude_mid {i in iset,j in jset2}  = T_magnitude[i,j];

### COMPUTE DIFFERENTIAL EQUATIONS AT MIDPOINT ###
var w_mid{i in iset,j in jset2} = 1 + f_mid[i,j]*cos(L_mid[i,j]) + g_mid[i,j]*sin(L_mid[i,j]);
var r_mid{i in iset,j in jset2}  = p_mid[i,j]/w_mid[i,j];
var s_mid{i in iset,j in jset2}  = sqrt(1 + h_mid[i,j]^2 + k_mid[i,j]^2);

##############################################################################
#                    SHADOW CONSTRAINTS AT MIDPOINT                          #
##############################################################################
var alpha_aux_2_mid{i in iset,j in jset2} = h_mid[i,j]^2 - k_mid[i,j]^2;
#
var sun_vector1_mid{i in iset,j in jset2} = s1[7] + s1[6]*t_mid[i,j] +s1[5]*t_mid[i,j]^2 +s1[4]*t_mid[i,j]^3+ s1[3]*t_mid[i,j]^4 +s1[2]*t_mid[i,j]^5 +s1[1]*t_mid[i,j]^6;
var sun_vector2_mid{i in iset,j in jset2} = s2[7] + s2[6]*t_mid[i,j] +s2[5]*t_mid[i,j]^2 +s2[4]*t_mid[i,j]^3+ s2[3]*t_mid[i,j]^4 +s2[2]*t_mid[i,j]^5 +s2[1]*t_mid[i,j]^6;
var sun_vector3_mid{i in iset,j in jset2} = s3[7] + s3[6]*t_mid[i,j] +s3[5]*t_mid[i,j]^2 +s3[4]*t_mid[i,j]^3+ s3[3]*t_mid[i,j]^4 +s3[2]*t_mid[i,j]^5 +s3[1]*t_mid[i,j]^6;
#
var r_vector1_mid{i in iset,j in jset2}   = r_mid[i,j]/s_mid[i,j]^2 *( cos(L_mid[i,j])+alpha_aux_2_mid[i,j]*cos(L_mid[i,j])+2*h_mid[i,j]*k_mid[i,j]*sin(L_mid[i,j]));
var r_vector2_mid{i in iset,j in jset2}   = r_mid[i,j]/s_mid[i,j]^2 *( sin(L_mid[i,j])-alpha_aux_2_mid[i,j]*sin(L_mid[i,j])+2*h_mid[i,j]*k_mid[i,j]*cos(L_mid[i,j]));
var r_vector3_mid{i in iset,j in jset2}   = r_mid[i,j]/s_mid[i,j]^2 * 2*(h_mid[i,j]*sin(L_mid[i,j])-k_mid[i,j]*cos(L_mid[i,j]));
#
var elem1_mid{i in iset,j in jset2}       = sun_vector1_mid[i,j]*r_vector1_mid[i,j] +sun_vector2_mid[i,j]*r_vector2_mid[i,j]+sun_vector3_mid[i,j]*r_vector3_mid[i,j];
var elem2_mid{i in iset,j in jset2}       = r_vector1_mid[i,j]^2 + r_vector2_mid[i,j]^2 + r_vector3_mid[i,j]^2;
var delta_mid{i in iset,j in jset2}       = - (elem1_mid[i,j] + sqrt(elem2_mid[i,j]-1));
#
#var signoid_mid{i in iset,j in jset2}     = 1 / (1+exp(delta_mid[i,j]/0.01));
#
param signoid_mid{i in iset,j in jset2}    := 1;

#C_elem2_mid{i in iset,j in jset2} : r_vector1_mid[i,j]^2 + r_vector2_mid[i,j]^2 + r_vector3_mid[i,j]^2 >=1;

### COMPUTE THE PERTURBATION DUE TO THE THRUST ##########
var aT_mid{i in iset,j in jset2}    = signoid_mid[i,j]*T_magnitude_mid[i,j]*T/m_mid[i,j];

var ar_T_mid{i in iset,j in jset2}  = aT_mid[i,j]*cos(beta_mid[i,j])*sin(alpha_mid[i,j]);
var ao_T_mid{i in iset,j in jset2}  = aT_mid[i,j]*cos(beta_mid[i,j])*cos(alpha_mid[i,j]);
var ah_T_mid{i in iset,j in jset2}  = aT_mid[i,j]*sin(beta_mid[i,j]);

### COMPUTE THE PEETURBATIONS DUE TO THE EARTH OBLATENESS ####
var ar_j2_mid{i in iset,j in jset2} =  -  3* J2/(2*r_mid[i,j]^4)*( 1 - 12*( h_mid[i,j]*sin(L_mid[i,j])-k_mid[i,j]*cos(L_mid[i,j]) )^2/s_mid[i,j]^4);
var ao_j2_mid{i in iset,j in jset2} =  - 12* J2/r_mid[i,j]^4 *( h_mid[i,j]*sin(L_mid[i,j])-k_mid[i,j]*cos(L_mid[i,j]) )*( h_mid[i,j]*cos(L_mid[i,j])+ k_mid[i,j]*sin(L_mid[i,j]) )/s_mid[i,j]^4;
var ah_j2_mid{i in iset,j in jset2} =  -  6* J2/r_mid[i,j]^4* ( h_mid[i,j]*sin(L_mid[i,j])-k_mid[i,j]*cos(L_mid[i,j]) )*(1 - h_mid[i,j]^2 - k_mid[i,j]^2)/s_mid[i,j]^4;

### COMPUTE THE WHOLE PERTURBATION #############
var ar_mid{i in iset,j in jset2} = ar_T_mid[i,j] + ar_j2_mid[i,j];
var ao_mid{i in iset,j in jset2} = ao_T_mid[i,j] + ao_j2_mid[i,j];
var ah_mid{i in iset,j in jset2} = ah_T_mid[i,j] + ah_j2_mid[i,j];

### COMPUTE THE DYNAMICAL EQUATIONS #####
var A1_mid{i in iset,j in jset2} = sqrt(p_mid[i,j])/w_mid[i,j]*(               0*ar_mid[i,j]                        +2*p_mid[i,j]*ao_mid[i,j]                                     +  0*ah_mid[i,j] );
var A2_mid{i in iset,j in jset2} = sqrt(p_mid[i,j])/w_mid[i,j]*(  w_mid[i,j]*sin(L_mid[i,j])*ar_mid[i,j]  + ((w_mid[i,j]+1)*cos(L_mid[i,j]) + f_mid[i,j])*ao_mid[i,j]  - g_mid[i,j]*(h_mid[i,j]*sin(L_mid[i,j]) -k_mid[i,j]*cos(L_mid[i,j]))*ah_mid[i,j] );
var A3_mid{i in iset,j in jset2}= sqrt(p_mid[i,j])/w_mid[i,j]*( -w_mid[i,j]*cos(L_mid[i,j])*ar_mid[i,j]  + ((w_mid[i,j]+1)*sin(L_mid[i,j]) + g_mid[i,j])*ao_mid[i,j]    + f_mid[i,j]*(h_mid[i,j]*sin(L_mid[i,j]) -k_mid[i,j]*cos(L_mid[i,j]))*ah_mid[i,j] );
var A4_mid{i in iset,j in jset2}= sqrt(p_mid[i,j])/w_mid[i,j]*(               0*ar_mid[i,j]                             +0*ao_mid[i,j]                    +  s_mid[i,j]^2*cos(L_mid[i,j])/2*ah_mid[i,j] );
var A5_mid{i in iset,j in jset2} = sqrt(p_mid[i,j])/w_mid[i,j]*(               0*ar_mid[i,j]                             +0*ao_mid[i,j]                   +  s_mid[i,j]^2*sin(L_mid[i,j])/2*ah_mid[i,j] );
var A6_mid{i in iset,j in jset2} = sqrt(p_mid[i,j])/w_mid[i,j]*(               0*ar_mid[i,j]                             +0*ao_mid[i,j]     +  (h_mid[i,j]*sin(L_mid[i,j]) -k_mid[i,j]*cos(L_mid[i,j]))*ah_mid[i,j] );
var A7_mid{i in iset,j in jset2} = - aT_mid[i,j]*m_mid[i,j]/(g0*Isp) ;

var sigma_mid{i in iset,j in jset2} = A6_mid[i,j] + sqrt(p_mid[i,j])*(w_mid[i,j]/p_mid[i,j])^2 ;

var pdot_mid {i in iset,j in jset2} = A1_mid[i,j]/sigma_mid[i,j];
var fdot_mid {i in iset,j in jset2} = A2_mid[i,j]/sigma_mid[i,j];
var gdot_mid {i in iset,j in jset2} = A3_mid[i,j]/sigma_mid[i,j];
var hdot_mid {i in iset,j in jset2} = A4_mid[i,j]/sigma_mid[i,j];
var kdot_mid {i in iset,j in jset2} = A5_mid[i,j]/sigma_mid[i,j];
var tdot_mid {i in iset,j in jset2} = 1/sigma_mid[i,j]*tc/(365*24*3600);
var mdot_mid {i in iset,j in jset2} = A7_mid[i,j]/sigma_mid[i,j]*tc*ac;

##################################################################################################
##################################################################################################
#                                  DIFFERENTIAL EQUATIONS                                        #
##################################################################################################
##################################################################################################
state_1 {i in iset,j in jset2}: ( p[i,j]-p[i,j+1] +step*scale[i]/6*( pdot[i,j]+pdot[i,j+1] + 4*pdot_mid[i,j]) )*1  = 0;   
state_2 {i in iset,j in jset2}: ( f[i,j]-f[i,j+1] +step*scale[i]/6*( fdot[i,j]+fdot[i,j+1] + 4*fdot_mid[i,j]) )*1  = 0;        
state_3 {i in iset,j in jset2}: ( g[i,j]-g[i,j+1] +step*scale[i]/6*( gdot[i,j]+gdot[i,j+1] + 4*gdot_mid[i,j]) )*1  = 0;          
state_4 {i in iset,j in jset2}: ( h[i,j]-h[i,j+1] +step*scale[i]/6*( hdot[i,j]+hdot[i,j+1] + 4*hdot_mid[i,j]) )*1  = 0;        
state_5 {i in iset,j in jset2}: ( k[i,j]-k[i,j+1] +step*scale[i]/6*( kdot[i,j]+kdot[i,j+1] + 4*kdot_mid[i,j]) )*1  = 0;      
state_6 {i in iset,j in jset2}: ( t[i,j]-t[i,j+1] +step*scale[i]/6*( tdot[i,j]+tdot[i,j+1] + 4*tdot_mid[i,j]) )*1  = 0;
state_7 {i in iset,j in jset2}: ( m[i,j]-m[i,j+1] +step*scale[i]/6*( mdot[i,j]+mdot[i,j+1] + 4*mdot_mid[i,j]) )*1  = 0;
#
##################################################################################################
## INCLUDE THE GEOBOX CONSTRAINT
##################################################################################################
########
# CONSTRAINT AT THE ASCENDING DESCENDING NODE
########
param Hgeo    := 35786e3;
param rgeo    := (Re + Hgeo)/Re;
param boxlimx1 := 75e+3/Re/2;
param boxlimx := 75e+3/Re/2;
param boxlimy := 75e+3/Re;
set N  := {0 .. nphases-4 };
set N3 := {phase_up .. nphases-4 by 2};
set N4 := {phase_down .. nphases-4 by 2};
#
# ASCENDING NODE CONSTRAINT
#
#z_c{i in 0 ..nphases-2} :r_vector3[i,nodes-1]/1000=0;
z_c{i in 0 ..nphases-2} : -1e-7 <= r_vector3[i,nodes-1] <= 1e-7;
#z_c{i in 0 ..nphases-2} : r_vector3[i,nodes-1] = 0;
z_c2{i in N3, j in 0..nodes-1}:r_vector3[i,j]>= -1e-7;
z_c3{i in N4, j in 0..nodes-1}:r_vector3[i,j]<=  1e-7;
#
var r_inplane{i in 0 ..nphases-4} = sqrt(r_vector1[i,nodes-1]^2+r_vector2[i,nodes-1]^2);
#
gbox_cons1  {i in 0 ..nphases-4}: (r_inplane[i]-rgeo)^2 >= boxlimx1^2  ;
#
############################################
# OPTION 1 
############################################
#
param phi_up   :=  boxlimy/rgeo;
param phi_down := -boxlimy/rgeo;
#
var inc        {i in N} = 2*atan(sqrt(h[i,nodes-1]^2+k[i,nodes-1]^2));
c_inc          {i in N} : (h[i,nodes-1]^2+k[i,nodes-1]^2)*10000 >= (tan(0.2*d2r/2))^2*100000 ;
#c_inc2        : h[nphases-4,nodes-1]^2+k[nphases-4,nodes-1]^2 = (tan(0.2*d2r/2))^2 ;

var a          {i in N} = p[i,nodes-1]/(1-f[i,nodes-1]^2-g[i,nodes-1]^2);
var e          {i in N} = sqrt(f[i,nodes-1]^2 + g[i,nodes-1]^2);
var ww         {i in N} = atan2(g[i,nodes-1]*h[i,nodes-1]-f[i,nodes-1]*k[i,nodes-1],f[i,nodes-1]*h[i,nodes-1] + g[i,nodes-1]*k[i,nodes-1]);
#
var theta_up   {i in N} = asin( sin( phi_up )/sin( inc[i] ) );
var theta_down {i in N} = asin( sin( phi_down )/sin( inc[i] ) );

var r_up_des   {i in  N3} = a[i] * ( 1 - e[i]^2 ) / (1 + e[i] * cos(pi -(ww[i]-theta_up[i])));
var r_down_des {i in  N3} = a[i] * ( 1 - e[i]^2 ) / (1 + e[i] * cos(pi -(ww[i]-theta_down[i])));

var r_up_as    {i in  N4}  = a[i] * ( 1 - e[i]^2 ) / (1 + e[i] * cos((ww[i]-theta_up[i])));
var r_down_as  {i in  N4}  = a[i] * ( 1 - e[i]^2 ) / (1 + e[i] * cos((ww[i]-theta_down[i])));

#c_r_up_as      {i in  N4}: r_down_as[i]*1  >= (rgeo+boxlimx)*1;
#c_r_up_as      {i in  N4}: r_up_as[i]*1    <= (rgeo-boxlimx)*1 ;
#c_r_down_as    {i in  N4}: r_down_as[i]*1  <= (rgeo-boxlimx)*1 ;
#c_r_up_des     {i in  N3}: r_up_des[i]*1   <= (rgeo-boxlimx)*1 ;
#c_r_down_des   {i in  N3}: r_down_des[i]*1 <= (rgeo-boxlimx)*1 ;

#c_r_down_as   {i in  N4}: r_down_as[i] <= r_up_as[i];
#
#geo_box_cons1{i in N3}: ( (rgeo-boxlimx) - r_down_des[i])*((rgeo-boxlimx)-r_up_des[i])/10000000 >= 0;
#geo_box_cons2{i in N3}: ( (rgeo+boxlimx) - r_down_des[i])*((rgeo+boxlimx)-r_up_des[i])/10000000 >= 0;
#
#geo_box_cons3{i in N4}: ( r_up_as[i]-(rgeo-boxlimx))*(r_down_as[i]-(rgeo-boxlimx)) >= 0;
#geo_box_cons4{i in N4}: ( r_down_as[i]-(rgeo+boxlimx))*(r_up_as[i]-(rgeo+boxlimx)) >= 0;
#
#geo_box_cons5{i in N4}: (r_up_as[i]-rgeo)^2/100   >= boxlimx1^2/100;
#geo_box_cons6{i in N4}: (r_down_as[i]-rgeo)^2/100 >= boxlimx1^2/100;

geo_box_cons7{i in N4}: (r_up_as[i]-rgeo  )^2 >= boxlimx1^2  ;
geo_box_cons8{i in N4}: (r_down_as[i]-rgeo)^2 >= boxlimx1^2  ;
geo_box_cons3{i in N4}: ( r_up_as[i]-rgeo)*(r_down_as[i]-rgeo) >= 0;
#
#geo_box_cons5{i in N4}: ( (rgeo-boxlimx/2) - r_down_as[i])*((rgeo-boxlimx/2)-r_up_as[i])/10000 >= 0;
#geo_box_cons6{i in N4}: ( (rgeo+boxlimx/2) - r_down_as[i])*((rgeo+boxlimx/2)-r_up_as[i])/10000 >= 0;
#
#geo_box_cons7{i in N4}: ( (rgeo-1*boxlimx/4) - r_down_as[i])*((rgeo-1*boxlimx/4)-r_up_as[i])/10000 >= 0;
#geo_box_cons8{i in N4}: ( (rgeo+1*boxlimx/4) - r_down_as[i])*((rgeo+1*boxlimx/4)-r_up_as[i])/10000 >= 0;
#
#geo_box_cons9{i in N4}: ( (rgeo-3*boxlimx/4) - r_down_as[i])*((rgeo-3*boxlimx/4)-r_up_as[i])/10000 >= 0;
#geo_box_cons10{i in N4}: ( (rgeo+3*boxlimx/4)- r_down_as[i])*((rgeo+3*boxlimx/4)-r_up_as[i])/10000 >= 0;
#
#geo_box_con11{i in N4}: (r_up_as[i]-rgeo)^2/10000 >= boxlimx^2/10000;;
#geo_box_cons12{i in N4}: (r_down_as[i]-rgeo)^2/10000 >= boxlimx^2/10000;;

##################################################################################################
## INCLUDE THE PHASING CONSTRAINT
##################################################################################################

#param desired_anomaly :=  3*pi/2;

#var true_anomaly   = atan2(r_vector2[nphases-1,nodes-1], r_vector1[nphases-1,nodes-1]);

#var true_anomaly    = L[nphases-1,nodes-1];

#c_phasing1: cos(desired_anomaly)/10000 = cos(true_anomaly)/10000;
#c_phasing2: sin(desired_anomaly)/10000 = sin(true_anomaly)/10000;

## OPTION2 ####
#c_phasing3:(desired_anomaly - 1*pi/180)/10000 <= true_anomaly/10000 <= (desired_anomaly + 1*pi/180)/10000;

##################################################################################################
#
# inequality constraints for the controls
#
ALPHA_C1{i in iset,j in jset}: -1   <= alpha[i,j]/pi     <= 1;
BETA_C1 {i in iset,j in jset}: -1   <= beta[i,j]/(pi/2)  <= 1;
#
# inequality constraints for the states
#
p_C1{i in iset,j in jset}     :  0    <= p[i,j]/9 <= 1 ;
f_C1{i in iset,j in jset}     : -1    <= f[i,j]/2 <= 1;
g_C1{i in iset,j in jset}     : -1    <= g[i,j]/2 <= 1;
h_C1{i in iset,j in jset}     : -1    <= h[i,j]/2 <= 1;
k_C1{i in iset,j in jset}     : -1    <= k[i,j]/2 <= 1;
t_C1{i in iset,j in jset}     :  0    <= t[i,j]/1 <= 1;
m_C1{i in iset,j in jset}     :  0    <= m[i,j]/1 <= 1;  
#
p_C1_mid{i in iset,j in jset2}:  0    <= p_mid[i,j]/9 <= 1 ;
f_C1_mid{i in iset,j in jset2}: -1    <= f_mid[i,j]/2 <= 1;
g_C1_mid{i in iset,j in jset2}: -1    <= g_mid[i,j]/2 <= 1;
h_C1_mid{i in iset,j in jset2}: -1    <= h_mid[i,j]/2 <= 1;
k_C1_mid{i in iset,j in jset2}: -1    <= k_mid[i,j]/2 <= 1;
t_C1_mid{i in iset,j in jset2}:  0    <= t_mid[i,j]/1 <= 1;
m_C1_mid{i in iset,j in jset2}:  0    <= m_mid[i,j]/1 <= 1; 
#
T_magnitude_C1{i in iset,j in jset} :  1  <= T_magnitude[i,j] <= 1;
#
# Continuity conditions
#
p_C2{i in 0 ..nphases-2}:   p[i,nodes-1]     = p[i+1,0];
f_C2{i in 0 ..nphases-2}:   f[i,nodes-1]     = f[i+1,0];
g_C2{i in 0 ..nphases-2}:   g[i,nodes-1]     = g[i+1,0];
h_C2{i in 0 ..nphases-2}:   h[i,nodes-1]     = h[i+1,0];
k_C2{i in 0 ..nphases-2}:   k[i,nodes-1]     = k[i+1,0];
t_C2{i in 0 ..nphases-2}:   t[i,nodes-1]     = t[i+1,0];
m_C2{i in 0 ..nphases-2}:   m[i,nodes-1]     = m[i+1,0]; 
a_C2{i in 0 ..nphases-2}:   alpha[i,nodes-1] = alpha[i+1,0]; 
b_C2{i in 0 ..nphases-2}:   beta[i,nodes-1]  = beta [i+1,0]; 
T_C2{i in 0 ..nphases-2}:   T_magnitude[i,nodes-1] = T_magnitude[i+1,0] ; 
L_C2{i in 0 ..nphases-2}:   L[i,nodes-1] = L[i+1,0] ; 
#
# Boundary values at the initial point
#
p_0:            p[0,0] = a0*(1-e0^2);
f_0:	        f[0,0] = e0*cos(w0+Om0);
g_0:         	g[0,0] = e0*sin(w0+Om0);
h_0:         	h[0,0] = tan(in0/2)*cos(Om0);
k_0:         	k[0,0] = tan(in0/2)*sin(Om0);
t_0:         	t[0,0] = 0;
m_0:         	m[0,0] = 1;
L_0:         	L_Om[0,0] = L0;
#
# boundary values at the final point
#
a_f:            p[nphases-1,nodes-1] = af*(1-ef^2);
e_f:            0<=f[nphases-1,nodes-1]^2+g[nphases-1,nodes-1]^2 <= ef^2;
i_f:            0<=h[nphases-1,nodes-1]^2+k[nphases-1,nodes-1]^2 <= tan(incf/2)^2;
#
########################### CALLING SOLVER #######################
option solver ipopt;
option presolve 0;
#mu_strategy adaptive
option ipopt_options "max_iter 10000 halt_on_ampl_error yes mu_init 1e-5 output_file ipoptlog.txt";
solve;
############################   DISPLAYS    #######################


printf {i in iset,j in 1 ..nodes-1}: "%17.15f %17.15f %17.15f %17.15f %17.15f %17.15f %17.15f\n",
 t[i,j], p[i,j], f[i,j], g[i,j] , h[i,j], k[i,j], m[i,j]> output.out;

printf {i in iset,j in 1 ..nodes-1}: "%17.15f  \n",
alpha[i,j] > alpha.out;

printf {i in iset,j in 1 ..nodes-1}: "%17.15f  \n",
T_magnitude[i,j] > T.out;
 
printf {i in iset,j in 1 ..nodes-1}: "%17.15f  \n",
beta[i,j] > beta.out;

printf {i in iset,j in 1 ..nodes-1}:"%17.15f  \n",
aT[i,j] > at.out;

printf {i in iset,j in 1 ..nodes-1}:"%19.15f  \n",
L[i,j] > L.out;

printf {i in iset,j in 1 ..nodes-1}:"%19.15f  \n",
signoid[i,j] > signoid.out;

printf {i in iset,j in 1 ..nodes-1}:"%17.15f  \n",
T_magnitude[i,j] > T_magnitude.out;





