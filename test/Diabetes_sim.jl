#Initalizations

##system states
Qsto1 = 0
Qsto2 = 0
Qgut = 0
#Gp = Gpb
#Gt = Gtb
#Gs = Gpb
Isc1 = 0
#Il = Ilb
#Ip = Ipb
X = 0
#I_prime = Ib
#XL = Ib
Isc2 = 0
y = [Qsto1;Qsto2;Qgut;Gp;Gt;Gs;Isc1;Il;Ip;X;I_prime;XL;Isc2]

##hyperparameters
###non-diabetic values
ke1 = 5e-4
ke2 = 339
kp1 = 2.7
kp2 = 2.1e-3
kp3 = 9e-3
kp4 = 6.18e-2
ki = 7.9e-3
Fsnc = 1
Vm0 = 2.5
Vmx = 4.7e-2
km0 = 225.59
p2U = 3.31e-2
Vg = 1.88
k1 = 6.5e-2
k2 = 7.9e-2
Vi = 5e-2
m1 = 0.19
m2 = 0.484
m4 = 0.194
m5 = 3.04e-2
m6 = 0.6471
HEb = 0.6
m3 = (HEb * m1)/(1-HEb)
kmax = 5.58e-2
kmin = 8e-3
kabs = 5.7e-2
f = 0.9
a = 1.3e-4
b = 0.82
d = 2.36e-3
ka1 = 1.8e-3
ka2 = 1.82e-2
kd = 1.64e-2
kgut = 3.78e-2

#Renal Extraction
if y[4] > ke2

    escr = ke1 * (y[4] - k2)

else

    escr = 0

end

#Plasma Glucose
dery[4] = max((kp1 - kp2*y[4]-kp3*y[12]),0) + max((kabs*y[3]*f/BW),0) - Fsnc - K1*y[4] + K2*y[5] - escr
dery[5] = K1*y[4] - K2*y[5] - (Vm0 + Vmx*y[10])/(km0 + y[5])*y[5]

#Absorption
dery[1] = -kmax*y[1]
Qsto = y[1] + y[2]

if (dose == zeros(5,1)) || conta == 1

    dery[2] = 0
    dery[3] = 0

else

    aa = 5/2/(1-b)/dosekempt
    cc = 5/2/d/dosekempt
    kgut = kmin + (kmax - kmin)/2*(tanh(aa*(Qsto-b*dosekempt)) - tanh(cc*(Qsto - d*dosekempt))+2)
    dery[2] = kmax*y[1] - y[2]*kgut
    dery[3] = kgut*y[2] - kabs*y[3]
 
end

#Subcutaneous glucose kinetics
dery[6] = -0.1*y[6] + 0.1*y[4]

#Plasma Insulin and PID
P = Kp*(y[6]/Vg - G_target)
D = Kp*TD*dery[6]/Vg
dery[13] = P/TI
PID = (P + D + y[13])/BW

if PID > 0

    dery[7] = PID - (kd + ka1)*y[7]

else

    dery[7] = -(kd + ka1)*y[7]

end

dery[14] = kd*y[7] - ka2*y[14]

#Plasma Insulin
dery[8] = -(m1 + m3)*y[8] + m2*y[9]
dery[9] = m1*y[8] - (m2 + m4)*y[9] + ka1*y[7] + ka2*y[14]

#Insulin Action
dery[10] = -p2U*y[10] + p2U*(y[9]/Vi-Ib)

#Production
dery[11] = ki*(y[9]/Vi-y[11])
dery[12] = ki*(y[11]-y[12])