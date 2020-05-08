clear;
clc;
% Main Dimension
rpm=input('Motor speed in rpm=')
f=input('Frequency of machine=')
pf=input('Desired power factor=')
W=input('Desired power(in KW)=')
e=input('Desired efficiency=')
Es=input('Stator voltage per phase=')
phase=input('phase=')
disp('Synchronus speed(in rps)')
Ns=rpm/60
disp('Poles,')
p=(2*f)/Ns
p=input('Assume poles=')
bav=input('Assume Bav=')
ac=input('Assume ac=')
kw=input('Assume Kw=')
disp('KVA input,')
Q=W/(pf*e)
disp('Output co-efficient,')
C0=W*kw*bav*ac*10^(-3)
disp('Product (D^2)*l,')
x=Q/(C0*Ns)
disp('For minimum cost L/t should be between 1.25 to 2. For good power factor L/t should be between 1 to 1.25')
l_t=input('Assume the value for L/t=')
D=((x*p)/(3.1416*l_t))^(1/3)
L=l_t*((3.1416*D)/p)
t=L/l_t
L=input('Assume L=')
D=input('Assume D=')
dn=input('Number of duct=')
dw=input('Duct width=')
Ki=input('Value of Ki=')
disp('Net iron length,')
Li=Ki*(L-dw*dn)
disp('Leminations of Lohys steel .5mm thick are used.')
% Stator Design
% Winding
disp('The stator winding is Delta connected in order that the machine has to be started by a star delta starter.')
disp('Flux per pole,')
f_p=bav*t*L
disp('Stator turns per pole,')
Ts=Es/(4.44*f*f_p*kw)
Ts=input('Assume stator turns per pole=')
q=input('Assume slots per pole per phase=')
disp('Stator slots,')
Ss=phase*q*p
disp('Stator slot pitch,')
yss=(3.1416*D)/Ss
disp('Total stator conductor,')
Tsc=p*Ts
disp('Conductor per slot,')
Zss=Tsc/Ss
Zss=input('Assume conductor per slot=')
disp('Total conductors used,')
Tc_1=Zss*Ss
disp('Stator turns per phase,')
Ts2=Tc_1/p
disp('The number of conductors used is very near to the calculated value and therefore the value of flux density need not be modified. Mush winding is parallel sided semi-enclosed slots is used.')
disp('Coil span,')
Cs=Ss/p
if rem(Cs,2)~=0
    disp('As coil span is odd,')
    kp=1
else
    disp('But even value of coil span is not expected, so coil span,')
    k_1=Cs;
    Cs=Cs-1
    disp('Angle of chording,')
    a=(1/k_1)*180
    disp('Pitch factor,')
    kp=cos(deg2rad(a/2))
end
disp('Distribution factor,')
kd=sin(deg2rad(60/2))/(q*sin(deg2rad(60/(2*q))))
disp('Stator winding factor,')
kws=kd*kp
% Conductor Size
disp('Stator current per phase,')
Is=(W*1000)/(3*Es*pf*e)
CD=input('Assume current density(A/m^2)=')
disp('Area of stator conductor,')
as=Is/CD
n_c=input('Number of conductor in parallel=')
d_cb=input('Diameter of conductor(bare) in metre=')
d_ci=input('Diameter of conductor(insulated) in metre=')
disp('Area of conductor provided,')
ad=n_c*(3.1416/4)*d_cb^2
% Slot Dimension
s_wc=input('Assume wire number along the width of the slot=')
disp('Slot width,')
Wss=(s_wc*d_ci)+(1/1000)+(1.48/1000)
s_dc=input('Assume wire number along the depth of the slot=')
disp('Slot depth,')
Wdss=(s_dc*d_ci)+(1.5/1000)+(3/1000)+(1/1000)+(2.43/1000)
disp('Slot pitch at,')
AA=(3.1416*(D*1000+8))/Ss
disp('Tooth width at AA,')
Wt=AA-Wss
disp('Flux density in teeth at AA,')
f_dt=f_p/((Ss/p)*Wt*Li)
disp('This is within limits. The dimension of stator slot are,')
Wss=input('Assume Wss in metre=')
Wdss=input('Assume Wdss in metre=')
disp('Length of mean turn,')
Lmis=2*L+2.3*t+.24
% Stator core
disp('Flux in core,')
f_c=f_p/2
f_ds=input('Assume flux density for the core in Wb/m^2=')
disp('Area of stator core,')
Acs=f_c/f_ds
disp('Depth of stator core,')
dcs=Acs/Li
dcs=input('The depth of core is used in metre=')
disp('Outer diameter of stator laminations,')
Do=D+2*(Wdss+dcs)
% Rotor Design
% Air Gap
disp('Length of air gap,')
Lg=(.2+2*(D*L)^(1/2))*10^-3
disp('Diameter of rotor,')
Dr=D-2*Lg
% Rotor Slots
disp('Number of rotor slots must be greater than stator slot.')
disp('Number of rotor slots,')
Sr=Ss+p/2
disp('Rotor slot pitch,')
ysr=(3.1416*Dr)/Sr
% Rotor Bars
disp('Current in each rotor bar,')
Ib=((2*3*kws*Ts)/Sr)*Is*pf
CD1=input('Assume current density in the bar=')
disp('Area of each bar,')
ab=Ib/CD1
disp('Assuming a rotor bar is used of')
Wrb=input('Width in metre=')
Drb=input('Depth in metre=')
disp('Area of bar used,')
ab=Wrb*Drb
disp('Dimension of rotor slot,')
Wsr=input('Width in metre=')
Dsr=input('Depth in metre=')
disp('Slot pitch at the root of teeth,')
Sprt=(3.1416*(Dr-2*Dsr))/Sr
disp('Slot width at the root of teeth,')
Wt1=Sprt-Wsr
disp('Flux density at the root of teeth,')
f_drt=f_p/((Sr/p)*Li*Wt1)
disp('This is within limits and therefore mmf required for rotor teeth would not be excessive.')
Ebar=input('Assuming extending the bars beyond the core on each side in metre is=')
Lsc=input('Assuming increase the length by scewing=')
disp('Length of each bar,')
Lb=L+2*Ebar+Lsc
SR=input('Specific resistance of rotor bar=')
disp('Resistance of each bar,')
rb=(SR*Lb)/(ab*1000000)
disp('Total Cu loss in bar,')
CL=Sr*(Ib^2)*rb
% End Rings
disp('End ring current,')
Ie=(Sr*Ib)/(3.1416*p)
disp('Area of each end ring,')
ae=Ie/CD1
disp('Using end ring dimension according to area,')
de=input('Depth of ring=')
te=input('Thickness of ring=')
disp('Area of ring provided,')
ae=de*te
d_or=Dr-2*Dsr;
d_er=d_or-2*de;
disp('Mean diameter of ring,')
De=(d_or+d_er)/2
disp('Resistance of each end ring,')
r_er=((SR*3.1416*De)/ae)*(10^(-6))
disp('Cu loss in 2 end rings,')
Cul_er=(2*(Ie^2)*r_er)
disp('Total rotor Cu loss,')
Trcl=Cul_er+CL
disp('Slip for induction motor,')
sim=Trcl/(Trcl+(W*1000))
% Rotor Core
disp('Depth of rotor core,')
dcr=dcs
disp('Inner diameter of rotor lamination,')
D_idrl=Dr-2*Dsr-2*dcr
% No Load Current
% Magnetising Current
% 1. Air Gap
s_o=input('Stator slot openning=')
disp('Ratio slot openning/gap length,')
rso=s_o/(Lg*1000)
C_C1=input('Taking Carter"s co-efficient for stator=')
disp('Gap contraction factor for stator slot,')
kgss=(yss*1000)/((yss*1000)-C_C1*s_o)
r_o=input('Rotor slot openning=')
disp('Ratio slot openning/gap length,')
rso1=r_o/(Lg*1000)
C_C2=input('Taking Carter"s co-efficient for rotor=')
disp('Gap contraction factor for rotor slot,')
kgsr=(ysr*1000)/((ysr*1000)-C_C2*r_o)
disp('Gap contraction factor for rotor slot,')
kgs=kgss*kgsr
disp('Ratio duct width/(1/2)*gap length,')
R_dwgl=dw/(1/2)*Lg
C_C3=input('Taking Carter"s co-efficient for duct=')
disp('Gap contraction factor for duct,')
kgd=L/(L-dn*C_C3*dw)
disp('Total gap contraction factor,')
kg=kgs*kgd
disp('Effective length of air gap,')
Lgs=kg*Lg
disp('Area of air gap,')
Ag=t*L
disp('Flux density in air gap at 60 degree interpolar axis,')
b60=1.36*bav
disp('mmf for air gap,')
ATg=800000*b60*Lgs
% Stator Teeth
disp('Width of stator teeth at 1/3 height from narrow end,')
W_st=((3.1416*(Ts2+(2*Wdss)/3))/Ss)-Wss
disp('Area of stator teeth per pole at 1/3 height,')
Ats=(Ss/p)*Li*W_st
disp('Flux density of stator teeth,')
Bts=f_p/Ats
Bts60=1.36*Bts
atts=input('Corresponding to this flux density, atts=')
disp('MMF required for stator teeth,')
ATts=atts*Wdss
% Stator Core
disp('Area of stator core,')
Acs=dcs*Li
disp('Flux density in stator core,')
Bcs=(f_p/2)/Acs
atcs=input('Corresponding to this flux density, atcs=')
disp('Length of path through stator core,')
lcs=(3.1416*(Ts2+dcs))/(phase*p)
disp('MMF required for stator core,')
ATcs=lcs*atcs
% Rotor Teeth
disp('Width of rotor teeth at 1/3 height from narrow end,')
W_tr=((3.1416*(Dr-(4/3)*Dsr))/Sr)-Wsr
disp('Area of rotor teeth per pole,')
a_rtp=(Sr/p)*Li*W_tr
disp('Flux density in rotor teeth at 1/3 height,')
Btr_3=f_p/a_rtp
Btr_60=1.36*Btr_3
attr=input('Corresponding to this flux density, attr=')
disp('MMF required for rotor teeth,')
ATtr=attr*Dsr
% Rotor Core
disp('Area of rotor core,')
Acr=Acs
disp('Flux density in rotor core,')
Bcr=Bcs
atcr=input('Corresponding to this flux density, atcr=')
disp('Length of path through rotor core,')
lcr=(3.1416*(Ts2-2*Dsr-dcs))/(phase*p)
disp('MMF required for rotor core,')
ATcr=lcr*atcr
disp('Magnetising current per phase,')
Im=(0.427*p*(ATg+ATts+ATcs+ATtr+ATcr))/(kws*Ts2)
% % Loss Component 
% % Iron Loss in Stator Teeth
% disp('Volume of stator teeth,')
% disp('Weight of stator teeth,')
% disp('Maximum flux density in teeth,')
% disp('Using 0.5 mm thick (Lohys steel) lamination.')
% =input('Corresponding to this flux density, loss per kg=')
% disp('Iron loss in stator teeth,')
% % Iron Loss in Stator Core
% disp('Volume of stator core,')
% disp('Weight of stator core,')
% disp('Flux density in core,')
% disp('Using 0.5 mm thick (Lohys steel) lamination.')
% =input('Corresponding to this flux density, loss per kg=')
% disp('Iron loss in stator core,')
% disp('Total iron loss,')
% disp('As the actual iron loss will be 2 times of this loss, total iron loss,')
% Short Circuit Current
% Circle Diagram
% Stator Tempurature Rise