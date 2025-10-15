# This script helps build an analytic PMSM or SYRM for a given requirements
# Copyright (C) 2025 Sai Gopal Rachakonda gopalsai20021909@gmail.com
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import math
import femm
import cmath 
import numpy as np

# === INITIAL MOTOR PARAMETERS ===
Pm =  2.2  # Mechanical power kW
Vl = math.sqrt(3)*230  # Line-to-Line Voltage in Volts
PF = 0.8  # Power Factor
Ns = 1500 #Speed
n_eff = 0.8


# === MOTOR TOPOLOGY ===
q = 48 #Number of slots
p = 8  #Number of poles
Np = 4 #Parallel paths
m = 3 #Number of phases
g = 0.45 #Air-gap

# == WINDING PARAMETERS ==
kw = 0.975 #Winding Factor
connection  = "STAR" #Winding connection
winding_layers = "DOUBLE" #Winding layers
gamma_emf = 0.91 #Eph/Vph

# == SLOT PARAMETERS ==
bs0 = 0.1 #Slot opening width
hs0 = 0.1 #Slot opening depth
hs1 = 0.1 #Slot opening dd

B_avg = 0.5  #Specific Magnetic Loading  (T)  ~~~~~~~2/pi * (0.8 to 1.05)
ac = 30  #Specific Electric Loading   (kA/m)  ~~~~~~~Z*Iz/(pi D)
Jsw = 5 #A/mm2
ar_m = 0.888 #L/D

Bst = 1.7 #Maximum tooth flux density
Bsy = 1.4 #Maximum stator yoke flux density
Bry = 1.4 #Maximum rotor yoke flux density
ki = 0.95 #Stacking factor
kf = 0.35 #Fill factor

# == Mechanical constrains
tr_t = 1 #Tangential rib thickness


# == INPUTS FEMM ==
filename = "Machine.fem"
BM = np.array([q, p])/math.gcd(q, p)
winding = np.array([
            [1, 1, -3, -3, 2, 2],
            [1, -3, -3, 2, 2, -1]
          ])

ar = ar_m*p/math.pi #aspect ratio L/pole_pitch

n_layers = 1 if winding_layers == "SINGLE" else 2 if winding_layers == "DOUBLE" else 0 #Number of coils per phase
nrps = (Ns/60)
G = 1.11*(math.pi**2)*kw*B_avg*ac #Machine constant ~~~~~Torque per unit volume i.e., P = G*(D^2)*L*rps

Pe = Pm/n_eff #Electrical Power 
P_g = Pe/PF #Airgap power "Q"
f = p*Ns/120 #Frequency
Il = 1000*Pe/(math.sqrt(3)*Vl*PF)  # Line Current in Amps
Vp = (Vl/math.sqrt(3)) if connection == "STAR" else Vl if connection == "DELTA" else 0 #Phase voltage
Ip = Il if connection == "STAR" else (Il/math.sqrt(3)) if connection == "DELTA" else 0 #Phase Current
Nc = (q/m/2) if winding_layers == "SINGLE" else (q/m) if winding_layers == "DOUBLE" else 0 #Number of coils per phase
Nc_p = Nc/Np #Number of coils in series per parallel-path
Ic = Ip/Np #Coil current
Vc = Vp/Nc_p #Coil voltage

#tau_p = math.pi*
D2L = P_g/G/nrps*1000000000 #mm^3

#Derived Dimension
D = math.ceil((D2L/ar_m)**(1/3))
L = math.ceil(D*ar_m)
Dr = D - g*2 #Rotot Outer Diameter


phi = B_avg*math.pi*D*L/1000000*1000 #Total flux mWb
phi_p = phi/p #Flux per pole mWb

phi_st_max = phi_p/2*math.sin(math.pi*p/q) #Maximum stator tooth flux

wt = phi_st_max/(Bst*L*ki)*1000 #Tooth width
bs1 = 2*(((math.tan(math.pi/q))*(D/2+hs0+hs1))-((wt/2)/(math.cos(math.pi/q))))

wsy = (phi_p/2)/(L/1000*ki*Bsy)


Ep = gamma_emf*Vp #Emf of a phase

Nt_ph = Ep/(4.44*f*kw*phi_p)*1000 #Number of series turns per phase
Nt_c = Nt_ph/Nc_p #Number of turns per coil

cAsc = Ic/Jsw #Single conductor area
cAca = cAsc*Nt_c #Copper area per coil
gAca = cAca/kf #Gross copper area i.e., slot area

sA = n_layers*gAca #Slot area for all the layers

d = (-bs1 + math.sqrt((bs1**2) + (4*math.tan(math.pi/q)*sA)))/(2*math.tan(math.pi/q))
bs2 = bs1 + 2*d*(math.tan(math.pi/q))


Ds = D + 2*(hs0+hs1+d+wsy)


addd = {
    "Bore Diameter (D)": D,
    "Stack Length (L)": L,
    "Flux per Pole (phi_p)": phi_p,
    "Max Stator Tooth Flux (phi_st_max)": phi_st_max,
    "Tooth Width (wt)": wt,
    "Slot Opening Width (bs1)": bs1,
    "Slot Depth (d)": d,
    "Slot Width (bs2)": bs2,
    "Stator Yoke Thickness (wsy)": wsy,
    "Stator Outer Diameter (Ds)": Ds
}

print(addd)


def drawStator(D, Ds, q, hs0, hs1, d, bs0, bs1, bs2, wsy, BM, winding, Nt_c, filename):
    femm.openfemm()
    femm.opendocument(filename)
    
    bore_radius = D / 2
    stator_outer_radius = Ds/2
    slot_angle = 360 / q    

    angle_slot_opening = math.asin(bs0 / (2 * bore_radius))
    
    # Adjusted slot points with shift applied
    temp0 = complex(bore_radius, 0)* (cmath.exp(1j * angle_slot_opening))
    temp1 = complex(bore_radius, 0)* (cmath.exp(1j * slot_angle/2 * math.pi/180))
    temp2 = complex(stator_outer_radius, 0)* (cmath.exp(1j * slot_angle/2 * math.pi/180))
    
    outer = (temp0.real,temp0.imag)
    top = (temp0.real + hs0 , 0)
    mid = (temp0.real + hs0 + hs1, 0)
    bottom1 = (temp0.real + hs0 + hs1 + d/2, 0)
    bottom2 = (temp0.real + hs0 + hs1 + d/2 + d/2, 0)
    
    
    Rlayer1 = temp0.real + hs0 + hs1 + d*0.25
    Rlayer2 = temp0.real + hs0 + hs1 + d*0.75
    
    # Add nodes for one side
    femm.mi_addnode(outer[0], bs0 / 2)
    femm.mi_addnode(top[0], bs0 / 2)
    femm.mi_addnode(mid[0], bs1 / 2)
    
    femm.mi_addnode(bottom1[0], bs1/2 + (bs2-bs1)/2/2)
    femm.mi_addnode(bottom2[0], bs1/2 + 2*(bs2-bs1)/2/2)
    
    femm.mi_addnode(temp1.real, temp1.imag)
    femm.mi_addnode(temp2.real, temp2.imag)
    
    # Mirror to create the other side
    femm.mi_addnode(outer[0], -bs0 / 2)
    femm.mi_addnode(top[0], -bs0 / 2)
    femm.mi_addnode(mid[0], -bs1 / 2)
    
    
    femm.mi_addnode(bottom1[0], -(bs1/2 + (bs2-bs1)/2/2))
    femm.mi_addnode(bottom2[0], -(bs1/2 + 2*(bs2-bs1)/2/2))
    
    
    femm.mi_addnode(temp1.real, -temp1.imag)
    femm.mi_addnode(temp2.real, -temp2.imag)
    
    # Add segments
    femm.mi_addsegment(outer[0], bs0 / 2, top[0], bs0 / 2)
    femm.mi_addsegment(outer[0], -bs0 / 2, top[0], -bs0 / 2)
    
    femm.mi_addsegment(top[0], bs0 / 2, mid[0], bs1 / 2)
    femm.mi_addsegment(top[0], -bs0 / 2, mid[0], -bs1 / 2)
    
    femm.mi_addsegment(mid[0], -bs1 / 2, mid[0], bs1 / 2)
    
    femm.mi_addsegment(mid[0], bs1 / 2, bottom1[0], (bs1/2 + (bs2-bs1)/2/2))
    femm.mi_addsegment(mid[0], -bs1 / 2, bottom1[0], -(bs1/2 + (bs2-bs1)/2/2))
    
    femm.mi_addsegment(bottom1[0], (bs1/2 + (bs2-bs1)/2/2), bottom1[0], -(bs1/2 + (bs2-bs1)/2/2))
    
    femm.mi_addsegment(bottom1[0], (bs1/2 + (bs2-bs1)/2/2), bottom2[0], bs2 / 2)
    femm.mi_addsegment(bottom1[0], -(bs1/2 + (bs2-bs1)/2/2), bottom2[0], -bs2 / 2)
        
    
    femm.mi_addsegment(bottom2[0], bs2 / 2, bottom2[0], -bs2 / 2)
    
    
    femm.mi_addsegment(temp1.real, temp1.imag, temp2.real, temp2.imag)
    femm.mi_addsegment(temp1.real, -temp1.imag, temp2.real, -temp2.imag)
    
    
    # Use an arc instead of a segment for the outer boundary
    #femm.mi_addarc(outer[0], -bs0 / 2, outer[0], bs0 / 2, angle_slot_opening*2*180/math.pi, 1)
    femm.mi_addarc(outer[0], bs0 / 2, temp1.real, temp1.imag, slot_angle/2 - angle_slot_opening*2*180/2/math.pi, 1)
    femm.mi_addarc(temp1.real, -temp1.imag, outer[0], -bs0 / 2, slot_angle/2 - angle_slot_opening*2*180/2/math.pi, 1)
    femm.mi_addarc(temp2.real, -temp1.imag, temp2.real, temp1.imag, slot_angle, 1)
    
    
       
    
    
    femm.mi_selectcircle(0,0,stator_outer_radius*1.5,4)
    femm.mi_setgroup(1001)
    
    femm.mi_selectgroup(1001)
    femm.mi_moverotate(0,0,slot_angle/2)
    femm.mi_clearselected()
    femm.mi_selectgroup(1001)
    femm.mi_copyrotate(0,0,slot_angle, BM[0]-1)
    
    femm.mi_addnode((bore_radius - g/3), 0)
    
    femm.mi_addsegment(bore_radius, 0, (bore_radius - g/3),0)
    
    
    femm.mi_addnode((bore_radius - g/3) * np.cos(np.radians(BM[0]*slot_angle)),
                    (bore_radius - g/3) * np.sin(np.radians(BM[0]*slot_angle)))
    
    femm.mi_addsegment(
                    (bore_radius) * np.cos(np.radians(BM[0]*slot_angle)),
                    (bore_radius) * np.sin(np.radians(BM[0]*slot_angle)),
                    (bore_radius - g/3) * np.cos(np.radians(BM[0]*slot_angle)),
                    (bore_radius - g/3) * np.sin(np.radians(BM[0]*slot_angle))
                    )
    
    femm.mi_addarc(bore_radius - g/3, 0, 
                   (bore_radius - g/3) * np.cos(np.radians(BM[0]*slot_angle)),
                   (bore_radius - g/3) * np.sin(np.radians(BM[0]*slot_angle)), BM[0]*slot_angle, 1)
        
    for i in range(int(BM[0])):
        yokeLabel = (stator_outer_radius - wsy/2)*np.exp(1j*np.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(yokeLabel.real, yokeLabel.imag)
        
        femm.mi_selectlabel(yokeLabel.real, yokeLabel.imag)
        femm.mi_setblockprop("US Steel Type 2-S 0.024 inch thickness", 0, wsy/5)
        femm.mi_clearselected()
        
        Lay1 = Rlayer1*np.exp(1j*np.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(Lay1.real, Lay1.imag)
        
        femm.mi_selectlabel(Lay1.real, Lay1.imag)
        femm.mi_setblockprop("24 SWG", 0, d/2 , f"fase{np.abs(winding[0,i])}",0,1001,Nt_c*winding[0,i]/np.abs(winding[0,i]))
        femm.mi_clearselected()
        
        Lay2 = Rlayer2*np.exp(1j*np.radians(slot_angle/2 + slot_angle*i))
        femm.mi_addblocklabel(Lay2.real, Lay2.imag)
        
        femm.mi_selectlabel(Lay2.real, Lay2.imag)
        femm.mi_setblockprop("24 SWG", 0, d/2 , f"fase{np.abs(winding[1,i])}",0,1001,Nt_c*winding[1,i]/np.abs(winding[1,i]))
        femm.mi_clearselected()
        
        
    femm.mi_selectcircle(0,0,stator_outer_radius*1.5,4)
    femm.mi_setgroup(1001)
    femm.mi_clearselected()
    
    femm.mi_saveas(filename)
    femm.mi_close()
    femm.closefemm()

def drawRotor(Dr, tr_t, g, p, BM, n_r_samples, n_th_samples, filename):
    femm.openfemm()
    femm.opendocument(filename)
    
    rotor_radius = Dr / 2
    agap_radius = rotor_radius + g/2
    shaft_radius = rotor_radius * 0.6
    pole_angle = 360 / p
    
    #n_r_samples = 10
    #n_th_samples = 5
    
    rotor_r_samples = np.append(0, np.linspace(shaft_radius, rotor_radius - tr_t, n_r_samples))
    rotor_th_samples = np.linspace(-pole_angle / 2, pole_angle / 2, n_th_samples)
    
    femm.mi_addnode(0, 0)
     
    for r_temp in range(n_r_samples):
        a0 = rotor_r_samples[r_temp]
        a1 = rotor_r_samples[r_temp + 1]
        
        for t_temp in range(n_th_samples - 1):
            b0 = rotor_th_samples[t_temp]
            b1 = rotor_th_samples[t_temp + 1]
            
            femm.mi_addnode(a1 * np.cos(np.radians(b0)), a1 * np.sin(np.radians(b0)))
            femm.mi_addsegment(a0 * np.cos(np.radians(b0)), a0 * np.sin(np.radians(b0)), 
                               a1 * np.cos(np.radians(b0)), a1 * np.sin(np.radians(b0)))
            femm.mi_addnode(a1 * np.cos(np.radians(b1)), a1 * np.sin(np.radians(b1)))
            femm.mi_addsegment(a0 * np.cos(np.radians(b1)), a0 * np.sin(np.radians(b1)), 
                               a1 * np.cos(np.radians(b1)), a1 * np.sin(np.radians(b1)))
            femm.mi_addarc(a1 * np.cos(np.radians(b0)), a1 * np.sin(np.radians(b0)), 
                           a1 * np.cos(np.radians(b1)), a1 * np.sin(np.radians(b1)), 
                           b1 - b0, 1)
    

    femm.mi_addnode(rotor_radius * np.cos(np.radians(rotor_th_samples[0])), 
                    rotor_radius * np.sin(np.radians(rotor_th_samples[0])))
    
    femm.mi_addsegment(rotor_radius * np.cos(np.radians(rotor_th_samples[0])), 
                       rotor_radius * np.sin(np.radians(rotor_th_samples[0])), 
                       rotor_r_samples[-1] * np.cos(np.radians(rotor_th_samples[0])), 
                       rotor_r_samples[-1] * np.sin(np.radians(rotor_th_samples[0])))
    
    
    
    femm.mi_addnode((rotor_radius + g/3) * np.cos(np.radians(rotor_th_samples[0])), 
                    (rotor_radius + g/3) * np.sin(np.radians(rotor_th_samples[0])))
    
    femm.mi_addsegment(rotor_radius * np.cos(np.radians(rotor_th_samples[0])), 
                       rotor_radius * np.sin(np.radians(rotor_th_samples[0])), 
                       (rotor_radius + g/3) * np.cos(np.radians(rotor_th_samples[0])), 
                       (rotor_radius + g/3) * np.sin(np.radians(rotor_th_samples[0])))
    
    
    #print(rotor_radius, ((rotor_th_samples[0])))
         
    femm.mi_addnode(rotor_radius * np.cos(np.radians(rotor_th_samples[-1])),
                    rotor_radius * np.sin(np.radians(rotor_th_samples[-1])))
    
    femm.mi_addsegment(rotor_radius * np.cos(np.radians(rotor_th_samples[-1])), 
                       rotor_radius * np.sin(np.radians(rotor_th_samples[-1])), 
                       rotor_r_samples[-1] * np.cos(np.radians(rotor_th_samples[-1])), 
                       rotor_r_samples[-1] * np.sin(np.radians(rotor_th_samples[-1])))

    femm.mi_addnode((rotor_radius + g/3) * np.cos(np.radians(rotor_th_samples[-1])), 
                    (rotor_radius + g/3) * np.sin(np.radians(rotor_th_samples[-1])))
    
    femm.mi_addsegment(rotor_radius * np.cos(np.radians(rotor_th_samples[-1])), 
                       rotor_radius * np.sin(np.radians(rotor_th_samples[-1])), 
                       (rotor_radius + g/3) * np.cos(np.radians(rotor_th_samples[-1])), 
                       (rotor_radius + g/3) * np.sin(np.radians(rotor_th_samples[-1])))

    
    femm.mi_addarc(rotor_radius * np.cos(np.radians(rotor_th_samples[0])), 
                   rotor_radius * np.sin(np.radians(rotor_th_samples[0])), 
                   rotor_radius * np.cos(np.radians(rotor_th_samples[-1])), 
                   rotor_radius * np.sin(np.radians(rotor_th_samples[-1])), 
                   rotor_th_samples[-1] - rotor_th_samples[0], 1)

    femm.mi_addarc((rotor_radius + g/3) * np.cos(np.radians(rotor_th_samples[0])), 
                   (rotor_radius + g/3) * np.sin(np.radians(rotor_th_samples[0])), 
                   (rotor_radius + g/3) * np.cos(np.radians(rotor_th_samples[-1])), 
                   (rotor_radius + g/3) * np.sin(np.radians(rotor_th_samples[-1])), 
                   rotor_th_samples[-1] - rotor_th_samples[0], 1)
                
    
    femm.mi_selectcircle(0,0,agap_radius,4)
    femm.mi_setgroup(2001)
    
    femm.mi_selectgroup(2001)
    femm.mi_moverotate(0,0,pole_angle/2)
    femm.mi_clearselected()
    femm.mi_selectgroup(2001)
    femm.mi_copyrotate(0,0,pole_angle, BM[1]-1)
    
    femm.mi_saveas(filename)
    femm.mi_close()
    femm.closefemm()



def rotorLables(Dr, tr_t, g, p, BM, n_r_samples, n_th_samples, filename):
    femm.openfemm()
    femm.opendocument(filename)
    
    rotor_radius = Dr / 2
    shaft_radius = rotor_radius * 0.6
    pole_angle = 360 / p
    
    ribring_radius = rotor_radius - tr_t/2
    
    
    n_r_samples = n_r_samples
    n_th_samples2 = int((n_th_samples-1)/2 + 1)
    
    
    shaft_mesh_size = shaft_radius
    rotor_mesh_size = (rotor_radius - shaft_radius)/n_r_samples*2
    rib_mesh_size = tr_t

    
    rotor_r_samples = np.append(0, np.linspace(shaft_radius, rotor_radius - tr_t, n_r_samples))
    rotor_th_samples = np.linspace(-pole_angle/2, 0, n_th_samples2)
    
    rotor_block_r = (rotor_r_samples[:-1] + np.diff(rotor_r_samples)/2).T
    rotor_block_th = (rotor_th_samples[:-1] + np.diff(rotor_th_samples)/2) #+ pole_angle/2
  
    core_block_nodes =  ribring_radius
    nomesh_block_nodes =  (rotor_block_r[0]*np.exp(1j*np.radians(rotor_block_th))).reshape(-1,1)
    
    opt_block_nodes = (((rotor_block_r[1:]).reshape(-1,1))*(np.exp(1j*np.radians(rotor_block_th)))).reshape(-1,1)


    temp_orig = core_block_nodes*np.exp(1j*np.radians(pole_angle/2))
    
    femm.mi_addblocklabel(temp_orig.real, temp_orig.imag)
    femm.mi_selectlabel(temp_orig.real, temp_orig.imag)
    femm.mi_setblockprop("US Steel Type 2-S 0.024 inch thickness", 0, rib_mesh_size)
    femm.mi_clearselected()
    
    for i in (nomesh_block_nodes):
        
        temp_orig = i*np.exp(1j*np.radians(pole_angle/2))
        temp_conj = np.conj(i)*np.exp(1j*np.radians(pole_angle/2))
        
        femm.mi_addblocklabel(*temp_orig.real, *temp_orig.imag)
        femm.mi_addblocklabel(*temp_conj.real, *temp_conj.imag)
        
        femm.mi_selectlabel(*temp_orig.real, *temp_orig.imag)
        femm.mi_setblockprop("Air", 0, shaft_mesh_size)
        femm.mi_clearselected()
        
        femm.mi_selectlabel(*temp_conj.real, *temp_conj.imag)
        femm.mi_setblockprop("Air", 0, shaft_mesh_size)
        femm.mi_clearselected()
        
        

    for i in (opt_block_nodes):
        
        temp_orig = i*np.exp(1j*np.radians(pole_angle/2))
        temp_conj = np.conj(i)*np.exp(1j*np.radians(pole_angle/2))
        
        femm.mi_addblocklabel(*temp_orig.real, *temp_orig.imag)
        femm.mi_addblocklabel(*temp_conj.real, *temp_conj.imag)
    
        femm.mi_selectlabel(*temp_orig.real, *temp_orig.imag)
        femm.mi_setblockprop("Air", 0, rotor_mesh_size)
        femm.mi_clearselected()
        
        femm.mi_selectlabel(*temp_conj.real, *temp_conj.imag)
        femm.mi_setblockprop("Air", 0, rotor_mesh_size)
        femm.mi_clearselected()
        

femm.openfemm()
femm.newdocument(0)
femm.mi_probdef(0, "millimeters", "planar", 1e-8, 0, 30)
femm.mi_getmaterial("US Steel Type 2-S 0.024 inch thickness")
femm.mi_getmaterial("Air")
femm.mi_getmaterial("24 SWG")
femm.mi_smartmesh(0)

for i in range(m):
    femm.mi_addcircprop(f"fase{i+1}", 0, 1)
    
femm.mi_saveas(filename)


n_r_samples = 10
n_th_samples = 11 #always be an odd number 
    
drawStator(D, Ds, q, hs0, hs1, d, bs0, bs1, bs2, wsy, BM, winding, Nt_c, filename)
drawRotor(Dr, tr_t, g, p, BM, n_r_samples, n_th_samples, filename)


rotorLables(Dr, tr_t, g, p, BM, n_r_samples, n_th_samples, filename)

