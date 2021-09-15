import matplotlib.pyplot as plt
import sympy as sym
import numpy as np
import time



def load_v3():
    import sympy as sym
    t = sym.Symbol('t')
    delta,Psi,alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

    V = (1 / 2) * (beta + (-1) * delta) ** (-1) * (beta ** 2 + ((-1) + beta) * delta) ** (-1) * (4 * beta + (- 1) * (1 + delta) ** 2) ** (-1) * Psi * (2 * sym.exp(1) ** (((-1) * t + t0) * beta) * IdA0 * ((-1) + beta) * beta * (beta + (-1) * delta) * Psi + (-2) * (alpha + (-1) * beta + delta) * (beta ** 2 + ((- 1) + beta) * delta) * Psi + sym.exp(1) ** ((1 / 2) * (t + (-1) * t0) * ((-1) + delta + (-1) * Psi)) * (IdA0 * beta * (beta + (-1) * delta) * ((-1) + (-1) * delta + beta * (3 + delta + (-1) * Psi) + Psi)+ (-1) * (beta ** 2 + (-1) * delta + beta * delta) * (alpha * (1 + (-2) * beta + delta + (-1) * Psi) + (beta + (-1) * delta) * ((-1) + 2 * IaA0 * beta + (-1) * delta + Psi + V0 * ((-1) + (-1) * delta + Psi))))+ sym.exp(1) ** ((1 / 2) * (t + (-1) * t0) * ((-1) + delta + Psi)) * ((-1) * IdA0 * beta * (beta+(-1) * delta) * ((-1) + (-1) * delta + (-1) * Psi + beta * (3 + delta + Psi)) + ( beta ** 2 + (-1) * delta+beta * delta) * (alpha * (1 + (-2) * beta + delta + Psi) + (beta + (-1) * delta) * ((-1) + 2 * IaA0 *  beta+(-1) * delta + (-1) * Psi + (-1) * V0 * (1 + delta + Psi)))))
    Iadap = (1 / 2) * sym.exp(1) ** (t0 + (-1) * t0 * delta + (-1 / 2) * t * ((-1) + 2 * beta + delta + Psi)) * (beta+(-1) * delta) ** (-1) * (beta ** 2 + ((-1) + beta) * delta) ** (-1) * ( 4 * beta + (-1) * ( 1 + delta) ** 2) ** (-1) * ( 2 * sym.exp(1) ** (t0 * ((-1) + beta + delta) + ( 1 / 2) * t * (( -1) + delta + Psi)) * IdA0 * beta * (beta + ( -1) * delta) * (4 * beta + (-1) * (1 + delta) ** 2) + (-2) * sym.exp(1) ** (t0 * ((-1) + delta) + (1 / 2) * t * ( (-1) + 2 * beta + delta + Psi)) * alpha * (beta ** 2 + ((-1) + beta) * delta) * ((-4) * beta + ( 1 + delta) ** 2) + sym.exp(1) ** ((1 / 2) * t0 * ((-1) + delta + (-1)* Psi) + t * ((-1) + beta + delta + Psi)) * ((-1) * IdA0 * beta * (beta + (-1) * delta) * ((-1) * (1 + delta) ** 2 + ((-1) + delta) * Psi + 2 * beta * (2 + Psi)) + (beta ** 2 + ((-1) + beta) * delta) * (alpha * (1+(-4) * beta + delta * ( 2 + delta + (-1) * Psi) + Psi) + (beta + (-1) * delta) * (4 * IaA0 * beta+(-2) * (1 + V0) * Psi + IaA0 * (1 + delta) * ((-1) + (-1) * delta + Psi))))+sym.exp(1) ** (t * ((-1) + beta + delta)+(1 / 2) * t0 * ((-1) + delta + Psi))*(IdA0 * beta * (beta + (-1) * delta) * ((1 + delta) ** 2 + 2 * beta * ((-2) + Psi) + ((-1) + delta) * Psi) + (beta ** 2 + ((-1) + beta) * delta)* (alpha * ((-4) * beta + (1 + delta) ** 2 + ((-1) + delta) * Psi) + (beta + (-1) * delta) * (4 * IaA0 * beta+2 * (1+V0) * Psi + (-1) * IaA0 * (1 + delta) * (1 + delta + Psi)))))
    Idep = sym.exp(1) ** (((-1) * t + t0) * beta) * IdA0

    return [V, Iadap, Idep]

def exp_cum(x, a, b):
    return a * (1 - np.exp(-b * x))

def monod(x, a, b, c, alp):
    return c + (a * np.exp(b) * x) / (alp + x)

tic = time.perf_counter()
EL=-65.88080929487171

#vres=-57.0
#vtm=-45.8

Vconvfact=-EL
vth=-0.7745615818727873#-vth*EL
vrm=-0.9521167423453459#vres/Vconvfact

t0_val=0

delta,Psi,alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

ts=np.inf

#neurone 958137
cost_idep_ini=0.3074538609199179
Cm=1996.500183108281
tao_m=98786.96214851207
sc=149.40116888490152
bet=0.22081378942479013
delta1=0.008911203404208803
cost_idep_ini=0.3074538609199179
Idep_ini_vr= 95.026900922085
psi1=0.36694230971007974
a=144.39555377980298
b=-0.015746877152471795
c=-29.665468077998163
alp=8.171843459232637
istim_min_spikinig_exp=200
istim_max_spikinig_exp=1000
#a=14.7314
#b=0.9767
#c=-10.
#alp=1.1144
time_scale=1 / (-sc / (Cm * EL))

d_dt=0.1
dt=d_dt/time_scale
sim_lenght=450
ref_t=2
t,delta,Psi,alpha, beta, gamma, IaA0, IdA0, t0, V0 = sym.symbols('t,delta,Psi,alpha,beta,gamma,IaA0,IdA0,t0,V0')

t0_val=0
[V,Iadap,Idep]=load_v3()
psi1=((-4)*bet+((1+delta1)**2))**(0.5)

Idep_ini=0
Iadap_ini=0
out=[]
t_out=[]
Istim=1800
Istim2=180
Istim3=1600
Istim4=180

t_final=t0_val+dt
v_ini=-1

mul=15

cor=np.ones(int(sim_lenght/d_dt))*Istim
change_cur=130
cor[int(change_cur/d_dt):int(sim_lenght/d_dt)+1]=np.ones(len(cor[int(change_cur/d_dt):int(sim_lenght/d_dt)+1]))*Istim2
beg_sign=30.0
cor[0:int(beg_sign/d_dt)]=0

change_cur=250
cor[int(change_cur/d_dt):int(sim_lenght/d_dt)+1]=np.ones(len(cor[int(change_cur/d_dt):int(sim_lenght/d_dt)+1]))*Istim3

change_cur=300
cor[int(change_cur/d_dt):int(sim_lenght/d_dt)+1]=np.ones(len(cor[int(change_cur/d_dt):int(sim_lenght/d_dt)+1]))*Istim4

#corr_list=[]
#fh = open('G1_1000.txt','r')
#for i in enumerate(fh):
#    corr_list.append(np.float64(i[1]))
#fh.close()
#cor=np.array(corr_list)
i=0


soglia_sign=10
Ide=[]
Iada=[]
Ide2=[]
Iada2=[]
ith=1.0969699757829885
tetalist=[]
tspk = []

curr_conv_fact = (bet*Cm*((1/(delta1*tao_m))**2))*(-EL) / (1/(delta1*tao_m));


while(t_final*time_scale<sim_lenght):
    #print(i)
    if cor[i] < ith:
        # print(i)
        Idep_ini = 0
        Iadap_ini = (cor[i] / sc) / (bet - delta1)
        if ((cor[i] / sc) / (bet - delta1) - 1)<v_ini:
            v_ini = ((cor[i] / sc) / (bet - delta1) - 1)
        out.append(v_ini)

    else:
        if cor[i] < cor[i-1]:
            # print(i)
            #Idep_ini =Iadap_ini-100*((cor[i-1]-cor[i])/(cor[i]))*(cor[i] / sc) / bet
            #versione originale
            #teta>-1
            #Idep_ini = Iadap_ini + teta *(cor[i] / sc) / bet
            #teta=(((out[i-1] - out[i - 2])/dt)-delta1*(1+out[i-1]))/(cor[i-1] / sc)-1
            teta=(out[i-1]/(cor[i-1] / sc))*(1/dt-delta1)-(out[i-2]/((cor[i-1] / sc)*dt))-delta1/(cor[i-1] / sc)-1
            Idep_ini = Iadap_ini + teta * (cor[i] / sc) / bet
            tetalist.append(teta)
            print("teta ", teta, (out[i-1]/(cor[i-1] / sc)), (1/dt-delta1), out[i-2], ((cor[i-1] / sc)*dt), delta1/(cor[i-1] / sc) )

    #print("adap nell'else")
    #print(Iadap_ini)
        out.append(V.subs(alpha, cor[i]/sc).subs(beta, bet).subs(delta, delta1).subs(t0, t0_val).subs(V0,v_ini ).subs(IaA0,Iadap_ini).subs(IdA0, Idep_ini).subs(Psi, psi1).subs(t,t_final))
        Iadap_ini = Iadap.subs(alpha, cor[i] / sc).subs(beta, bet).subs(delta, delta1).subs(t0, t0_val).subs(V0,v_ini).subs(IaA0, Iadap_ini).subs(IdA0, Idep_ini).subs(Psi, psi1).subs(t, t_final)
        Idep_ini = Idep.subs(beta, bet).subs(IdA0, Idep_ini).subs(t0, t0_val).subs(t, t_final)
        v_ini = out[i]
        #Iada.append(Iadap_ini)
        #Ide.append(Idep_ini)


    t_out.append(t_final)
    Iada.append(Iadap_ini)
    Ide.append(Idep_ini)

    if cor[i]>0:
        if cor[i-1]<soglia_sign or i==0:
            init_sign=t_final
            Idep_ini = cost_idep_ini*cor[i]
    if v_ini>vth:
        v_ini=vrm


        print('************')
        print('val_ist V')
        print(v_ini * Vconvfact)
        print('adap')
        print(Iadap_ini)
        print('t_fin')
        print(t_final)
        print('t_ini')
        print(init_sign)
        print('************')
        tspk.append(t_final*time_scale)

        if cor[i]<istim_min_spikinig_exp or cor[i]>istim_max_spikinig_exp:

            c_aux=0.8*Idep_ini_vr + (cor[i]/(sc)) / bet+(delta1/bet)*(1+vrm)-a*np.exp(b*cor[i]/1000)
            Iadap_ini = monod((t_final - init_sign) * time_scale, a, b * cor[i] / 1000, c_aux, alp)
        else:
            Iadap_ini=monod((t_final-init_sign)*time_scale,a,b*cor[i]/1000,c,alp)
        if cor[i]<1:
            Idep_ini=0
            Iadap_ini =0
        else:
            Idep_ini=Idep_ini_vr
        print('adap')
        print(Iadap_ini)
        #print(v_ini)
        #print(Iadap_ini)
        #print(Idep_ini)
        for k in range(int(ref_t / d_dt)):
            out.append(v_ini)
            t_final = t_final + dt
            t_out.append(t_final)
            Iada.append(Iadap_ini)
            Ide.append(Idep_ini)
            i = i + 1
        #print('*************init*************')
        #print(t_final)
        #print(v_ini)
        #print(Iadap_ini)
        #print(Idep_ini)
        #print('***************fin******************')



    i = i + 1
    t0_val = t_final
    t_final = t0_val + dt
plt.figure();
plt.scatter(np.array(t_out)*time_scale,np.array(out)*Vconvfact,s=2)
plt.title('voltage')
plt.figure();
plt.scatter(np.array(t_out)*time_scale,np.array(cor[0:len(t_out)]),s=2)
plt.title('corrente')
print("lengths ",len(t_out), len(out), len(Iada), len(Ide))
plt.figure();
plt.scatter(np.array(t_out)*time_scale,np.array(Iada)*curr_conv_fact,s=2)
plt.title('Iadap')
plt.figure();
plt.scatter(np.array(t_out)*time_scale,np.array(Ide)*curr_conv_fact,s=2)
plt.title('Idep')


toc = time.perf_counter()
print(f"time: {toc - tic:0.4f} seconds")
plt.show()

print("tspk", tspk)

with open('variables_neuron_simulator.dat', 'w') as f:
    for i in range(len(t_out)):
        f.write(str(t_out[i]*time_scale)+" "+str(out[i]*Vconvfact)+" "+str(Iada[i]*curr_conv_fact)+" "+str(Ide[i]*curr_conv_fact)+" "+str(cor[i])+"  \n")
