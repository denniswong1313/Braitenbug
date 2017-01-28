from brian2 import *
import matplotlib.pyplot as plt
map_size = 100
global foodx, foody, food_count, bug_plot, food_plot, sr_plot, sl_plot
foodx = -50
foody = 50
food_count = 0


# Sensor neurons
a = 0.02/ms
b = 0.2/ms
c = -65*mV
d = 0.5*mV/ms
Erev = -80*mV
I0 = 750
I1 = 1000
tau = 5*ms

# --- *** Small bug sensor equations *** ---
sensor_eqsL = '''
dvL/dt = (0.04/ms/mV)*vL**2 + (5/ms)*vL + (140*mV/ms) - uL + iL :  volt
duL/dt = a*(b*vL - uL): volt / second

iL : volt/second

x : 1
y : 1
I :1
x_disp : 1
y_disp : 1
foodx : 1
foody : 1
smell: 1

'''

sensor_eqsR = '''
dvR/dt = (0.04/ms/mV)*vR**2 + (5/ms)*vR + (140*mV/ms) - uR + iR :  volt
duR/dt = a*(b*vR - uR): volt / second

iR : volt/second

x : 1
y : 1
I :1
x_disp : 1
y_disp : 1
foodx : 1
foody : 1
smell: 1
'''

# --- *** Alpha bug additional sensor equations *** ---
sensorA_eqsL = '''
dvL/dt = (0.04/ms/mV)*vL**2 + (5/ms)*vL + (140*mV/ms) - uL + iLA + iLA2 :  volt
duL/dt = a*(b*vL - uL): volt / second

iLA : volt/second
iLA2 : volt/second

x : 1
y : 1
I :1
x_disp : 1
y_disp : 1
foodx : 1
foody : 1
smell : 1
'''

sensorA_eqsR = '''
dvR/dt = (0.04/ms/mV)*vR**2 + (5/ms)*vR + (140*mV/ms) - uR + iRA + iRA2 :  volt
duR/dt = a*(b*vR - uR): volt / second

iRA : volt/second
iRA2 : volt/second

x : 1
y : 1
I :1
x_disp : 1
y_disp : 1
foodx : 1
foody : 1
smell:1
'''

# --- *** Create groups for small bug sensors *** ---
#             ****** Bug 1 ******
sr = NeuronGroup(1, sensor_eqsR,
                    threshold =
                    'vR > 30*mV',
                    reset = '''
                    vR = c
                    uR+=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
sr.x_disp = 5
sr.y_disp = 5
sr.x = sr.x_disp
sr.y = sr.y_disp
sr.foodx = foodx
sr.foody = foody
sr.vR = c
sr.uR = b*c

sl = NeuronGroup(1, sensor_eqsL,
                    threshold =
                    'vL > 30*mV',
                    reset = '''
                    vL = c
                    uL +=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
sl.x_disp = -5
sl.y_disp = 5
sl.x = sl.x_disp
sl.y = sl.y_disp
sl.foodx = foodx
sl.foody = foody
sl.vL = c
sl.uL = b*c

#             ****** Bug 2 ******
sr2 = NeuronGroup(1, sensor_eqsR,
                    threshold =
                    'vR > 30*mV',
                    reset = '''
                    vR = c
                    uR+=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
sr2.x_disp = 5
sr2.y_disp = 5
sr2.x = sr2.x_disp
sr2.y = sr2.y_disp
sr2.foodx = foodx
sr2.foody = foody
sr2.vR = c
sr2.uR = b*c

sl2 = NeuronGroup(1, sensor_eqsL,
                    threshold =
                    'vL > 30*mV',
                    reset = '''
                    vL = c
                    uL +=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
sl2.x_disp = -5
sl2.y_disp = 5
sl2.x = sl2.x_disp
sl2.y = sl2.y_disp
sl2.foodx = foodx
sl2.foody = foody
sl2.vL = c
sl2.uL = b*c

#             ****** Alpha bug ******
srA = NeuronGroup(1, sensor_eqsR,
                    threshold =
                    'vR > 30*mV',
                    reset = '''
                    vR = c
                    uR+=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
srA.x_disp = 5
srA.y_disp = 5
srA.x = sr.x_disp
srA.y = sr.y_disp
srA.foodx = foodx
srA.foody = foody
srA.vR = c
srA.uR = b*c

slA = NeuronGroup(1, sensor_eqsL,
                    threshold =
                    'vL > 30*mV',
                    reset = '''
                    vL = c
                    uL +=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
slA.x_disp = -5
slA.y_disp = 5
slA.x = sl.x_disp
slA.y = sl.y_disp
slA.foodx = foodx
slA.foody = foody
slA.vL = c
slA.uL = b*c

srA_pred = NeuronGroup(1, sensorA_eqsR,
                    threshold =
                    'vR > 30*mV',
                    reset = '''
                    vR = c
                    uR+=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
srA_pred.x_disp = 5
srA_pred.y_disp = 5
srA_pred.x = sr.x_disp
srA_pred.y = sr.y_disp
srA_pred.foodx = foodx
srA_pred.foody = foody
srA_pred.vR = c
srA_pred.uR = b*c

slA_pred = NeuronGroup(1, sensorA_eqsL,
                    threshold =
                    'vL > 30*mV',
                    reset = '''
                    vL = c
                    uL +=d
                    ''',
                    clock=Clock(0.2*ms),
                    method = 'Euler')
slA_pred.x_disp = -5
slA_pred.y_disp = 5
slA_pred.x = sl.x_disp
slA_pred.y = sl.y_disp
slA_pred.foodx = foodx
slA_pred.foody = foody
slA_pred.vL = c
slA_pred.uL = b*c

# --- *** Small bug motoneuron equations *** ---
moto_eqsL = '''
dvmL/dt = (0.04/ms/mV)*vmL**2 + (5/ms)*vmL + (140*mV/ms) - umL + imL :  volt
dumL/dt = a*(b*vmL - umL): volt / second
dgL/dt = -gL/tau + zL/ms : siemens
dzL/dt = -zL/tau : siemens
imL= -gL*(vmL)/(nsiemens*ms) : volt/second
'''

moto_eqsR = '''
dvmR/dt = (0.04/ms/mV)*vmR**2 + (5/ms)*vmR + (140*mV/ms) - umR + imR :  volt
dumR/dt = a*(b*vmR - umR): volt / second
dgR/dt = -gR/tau + zR/ms : siemens
dzR/dt = -zR/tau :siemens
imR= -gR*(vmR)/(nsiemens*ms) : volt/second
'''

# --- *** Create groups for small bug motoneurons *** ---
#             ****** Bug 1 ******
mL = NeuronGroup(1,moto_eqsL,
                    threshold=
                    'vmL>30*mV',
                    reset=
                    '''
                    vmL = c
                    umL +=d
                    ''',
                    clock = Clock(0.2*ms),
                    method = 'Euler')
                    
mR = NeuronGroup(1,moto_eqsR,
                    threshold=
                    'vmR>30*mV',
                    reset=
                    '''
                    vmR = c
                    umR +=d
                    ''',
                    clock = Clock(0.2*ms),
                    method = 'Euler')
                    
mL.vmL = c
mL.umL = b*c
mR.vmR = c
mR.umR = b*c

#             ****** Bug 2 ******
mL2 = NeuronGroup(1,moto_eqsL,
                    threshold=
                    'vmL>30*mV',
                    reset=
                    '''
                    vmL = c
                    umL +=d
                    ''',
                    clock = Clock(0.2*ms),
                    method = 'Euler')
                    
mR2 = NeuronGroup(1,moto_eqsR,
                    threshold=
                    'vmR>30*mV',
                    reset=
                    '''
                    vmR = c
                    umR +=d
                    ''',
                    clock = Clock(0.2*ms),
                    method = 'Euler')
                    
mL2.vmL = c
mL2.umL = b*c
mR2.vmR = c
mR2.umR = b*c

#             ****** Alpha Bug ******
mLA = NeuronGroup(1,moto_eqsL,
                    threshold=
                    'vmL>30*mV',
                    reset=
                    '''
                    vmL = c
                    umL +=d
                    ''',
                    clock = Clock(0.2*ms),
                    method = 'Euler')
                    
mRA = NeuronGroup(1,moto_eqsR,
                    threshold=
                    'vmR>30*mV',
                    reset=
                    '''
                    vmR = c
                    umR +=d
                    ''',
                    clock = Clock(0.2*ms),
                    method = 'Euler')
                    
mLA.vmL = c
mLA.umL = b*c
mRA.vmR = c
mRA.umR = b*c

# --- *** Create synapses between sensors and motoneurons for small bugs *** ---
#             ****** Bug 1 ******
syn_intr1 = Synapses(sl, mR, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zR += w
                    ''')
syn_intr1.connect(i=[0],j=[0])
syn_intr1.w= (0.05*ms/(tau*exp(-1))) * nsiemens

syn_intl1 = Synapses(sr, mL, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zL += w
                    ''')
syn_intl1.connect(i=[0],j=[0])
syn_intl1.w= (0.05*ms/(tau*exp(-1))) * nsiemens

#             ****** Bug 2 ******
syn_intr2 = Synapses(sl2, mR2, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zR += w
                    ''')
syn_intr2.connect(i=[0],j=[0])
syn_intr2.w= (0.05*ms/(tau*exp(-1))) * nsiemens

syn_intl2 = Synapses(sr2, mL2, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zL += w
                    ''')
syn_intl2.connect(i=[0],j=[0])
syn_intl2.w= (0.05*ms/(tau*exp(-1))) * nsiemens

#             ****** Alpha Bug ******
syn_intrA = Synapses(srA, mRA, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zR += w
                    ''')
syn_intrA.connect(i=[0],j=[0])
syn_intrA.w= (0.05*ms/(tau*exp(-1))) * nsiemens

syn_intlA = Synapses(slA, mLA, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zL += w
                    ''')
syn_intlA.connect(i=[0],j=[0])
syn_intlA.w= (0.05*ms/(tau*exp(-1))) * nsiemens

# --- *** The small bugs *** ---

taum = 4*ms
base_speed = 40
turn_rate = 3000*Hz

bug_eqs = '''
dvl/dt=  (-vl/taum) : 1
dvr/dt = (-vr/taum) : 1
speed = (vl + vr)/2 + base_speed : 1
dangle/dt = (vr - vl)*turn_rate : 1
dx/dt = speed*cos(angle)*15*Hz : 1
dy/dt = speed*sin(angle)*15*Hz : 1
tl:1
tr:1
'''

bug1 = NeuronGroup(1, bug_eqs, clock=Clock(0.2*ms), method = 'Euler')

bug1.angle = pi/2
bug1.x = -50
bug1.y = -50

bug2 = NeuronGroup(1, bug_eqs, clock=Clock(0.2*ms), method = 'Euler')

bug2.angle = pi/2
bug2.x = -50
bug2.y = -55

# --- *** The Alpha Bug *** ---
taumA = 10*ms
base_speedA = 1
turn_rateA = 10000*Hz
baseline = 150

# Remove base speed dependency
# Bug speed determined by baseline
bugA_eqs = '''
dvl/dt=  (-vl+baseline)/taumA : 1
dvr/dt = (-vr+baseline)/taumA : 1
speed = (vl + vr)/2 : 1 
dangle/dt = (vr - vl)*turn_rateA : 1
dx/dt = speed*cos(angle)*15*Hz : 1
dy/dt = speed*sin(angle)*15*Hz : 1
tl:1
tr:1
'''

bugA = NeuronGroup(1, bugA_eqs, clock=Clock(0.2*ms), method = 'Euler')

bugA.angle = pi/2
bugA.x = 0
bugA.y = 0

# --- *** Synapses (motoneurons communicate with bug motor) *** ---
#             ****** Bug 1 ******
syn_r1 = Synapses(mR, bug1, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vr += w
                    ''')
syn_r1.connect(i=[0],j=[0])
syn_r1.w= (0.1*ms/(taum*exp(-1)))

syn_l1 = Synapses(mL, bug1, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vl += w
                    ''')
syn_l1.connect(i=[0],j=[0])
syn_l1.w= (0.1*ms/(taum*exp(-1)))

#             ****** Bug 2 ******
syn_r2 = Synapses(mR2, bug2, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vr += w
                    ''')
syn_r2.connect(i=[0],j=[0])
syn_r2.w= (0.1*ms/(taum*exp(-1)))

syn_l2 = Synapses(mL2, bug2, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vl += w
                    ''')
syn_l2.connect(i=[0],j=[0])
syn_l2.w= (0.1*ms/(taum*exp(-1)))

#             ****** Alpha Bug ******
syn_rA = Synapses(mRA, bugA, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vr -= w
                    ''')
syn_rA.connect(i=[0],j=[0])
syn_rA.w= (0.05*ms/(taumA*exp(-1)))

syn_lA = Synapses(mLA, bugA, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vl -= w
                    ''')
syn_lA.connect(i=[0],j=[0])
syn_lA.w= (0.05*ms/(taumA*exp(-1)))

syn_intrA_pred = Synapses(slA_pred, bugA, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vr += w
                    ''')
syn_intrA_pred.connect(i=[0],j=[0])
syn_intrA_pred.w= (0.05*ms/(taumA*exp(-1)))

syn_intlA_pred = Synapses(srA_pred, bugA, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vl += w
                    ''')
syn_intlA_pred.connect(i=[0],j=[0])
syn_intlA_pred.w= (0.05*ms/(taumA*exp(-1)))
                    


# --- *** Monitor Spikes *** ---
#MSR = StateMonitor(srA_pred,('vR'),record=True)
#MSL = StateMonitor(slA_pred,('vL'),record=True)
#
##MMR = StateMonitor(mR,('vmR'),record=True)
##MML = StateMonitor(mL,('vmL'),record=True)
#
#MBR = StateMonitor(bugA,('vr','x','y'),record=True)
#MBL = StateMonitor(bugA,('vl'),record=True)
# --- -*** End monitor spikes *** ---

f = figure(1)
bug1_plot = plot(bug1.x, bug1.y, 'ko')
bug2_plot = plot(bug2.x, bug2.y, 'ko')
bugA_plot = plot(bugA.x, bugA.y, 'ko')
food_plot = plot(foodx, foody, 'b*')
sr1_plot = plot([0], [0], 'w')
sl1_plot = plot([0], [0], 'w')
sr2_plot = plot([0], [0], 'w')   
sl2_plot = plot([0], [0], 'w')
srA_plot = plot([0], [0], 'w')   
slA_plot = plot([0], [0], 'w')
# Additional update rules (not covered/possible in above eqns)
# Initialize the target for the guardian bug
target = bug1

#This block of code updates the position of the bug and food and makes the appropriate rotations
@network_operation()
def update_positions(dt=1*ms):
    global foodx, foody, food_count, target

#            --- *** Sensors detect distance to food *** ---
#                           ****** Bug 1 ******
    sr.iR = I0 / sqrt(((sr.x-foodx)**2+(sr.y-foody)**2)) * mV/ms
    sl.iL = I0 / sqrt(((sl.x-foodx)**2+(sl.y-foody)**2)) * mV/ms
    
    sr.x = bug1.x + sr.x_disp*cos(bug1.angle-pi/2) - sr.y_disp*sin(bug1.angle-pi/2)
    sr.y = bug1.y + sr.x_disp*sin(bug1.angle-pi/2) + sr.y_disp*cos(bug1.angle-pi/2)

    sl.x = bug1.x + sl.x_disp*cos(bug1.angle-pi/2) - sl.y_disp*sin(bug1.angle-pi/2)
    sl.y = bug1.y + sl.x_disp*sin(bug1.angle-pi/2) + sl.y_disp*cos(bug1.angle-pi/2)

#                           ****** Bug 2 ******
    sr2.iR = I0 / sqrt(((sr2.x-foodx)**2+(sr2.y-foody)**2)) * mV/ms
    sl2.iL = I0 / sqrt(((sl2.x-foodx)**2+(sl2.y-foody)**2)) * mV/ms
    
    sr2.x = bug2.x + sr2.x_disp*cos(bug2.angle-pi/2) - sr2.y_disp*sin(bug2.angle-pi/2)
    sr2.y = bug2.y + sr2.x_disp*sin(bug2.angle-pi/2) + sr2.y_disp*cos(bug2.angle-pi/2)

    sl2.x = bug2.x + sl2.x_disp*cos(bug2.angle-pi/2) - sl2.y_disp*sin(bug2.angle-pi/2)
    sl2.y = bug2.y + sl2.x_disp*sin(bug2.angle-pi/2) + sl2.y_disp*cos(bug2.angle-pi/2)
    
#                           ****** Alpha Bug ******
    srA.iR = I0 / sqrt(((srA.x-foodx)**2+(srA.y-foody)**2)) * mV/ms
    slA.iL = I0 / sqrt(((slA.x-foodx)**2+(slA.y-foody)**2)) * mV/ms
    
    srA_pred.iRA = I1 / (sqrt((srA.x-target.x)**2+(srA.y-target.y)**2))**6 * mV/ms
    slA_pred.iLA = I1 / (sqrt((slA.x-target.x)**2+(slA.y-target.y)**2))**6 * mV/ms
    
    srA.x = bugA.x + srA.x_disp*cos(bugA.angle-pi/2) - srA.y_disp*sin(bugA.angle-pi/2)
    srA.y = bugA.y + srA.x_disp*sin(bugA.angle-pi/2) + srA.y_disp*cos(bugA.angle-pi/2)

    slA.x = bugA.x + slA.x_disp*cos(bugA.angle-pi/2) - slA.y_disp*sin(bugA.angle-pi/2)
    slA.y = bugA.y + slA.x_disp*sin(bugA.angle-pi/2) + slA.y_disp*cos(bugA.angle-pi/2)

#    --- *** What to do when food gets eaten *** ---
    if ((bug1.x-foodx)**2+(bug1.y-foody)**2) < 16:
	food_count += 1
	foodx = randint(-map_size+10, map_size-10)
	foody = randint(-map_size+10, map_size-10)
	
    if ((bug2.x-foodx)**2+(bug2.y-foody)**2) < 16:
	food_count += 1
	foodx = randint(-map_size+10, map_size-10)
	foody = randint(-map_size+10, map_size-10)
	
    if ((bugA.x-bug1.x)**2+(bugA.y-bug1.y)**2) < 16:
	bug1.x = randint(-map_size+10, map_size-10)
	bug1.y = randint(-map_size+10, map_size-10)
	if ((foodx-bug1.x)**2 + (foody-bug1.y)**2) > ((foodx-bug2.x)**2 + (foody-bug2.y)**2):
	    target = bug2
	
    if ((bugA.x-bug2.x)**2+(bugA.y-bug2.y)**2) < 16:
	bug2.x = randint(-map_size+10, map_size-10)
	bug2.y = randint(-map_size+10, map_size-10)
	if ((foodx-bug1.x)**2 + (foody-bug1.y)**2) < ((foodx-bug2.x)**2 + (foody-bug2.y)**2):
	    target = bug1
	    
# Update target if one bug becomes closer to the food than the other    
    if ((foodx-bug1.x)**2 + (foody-bug1.y)**2) < ((foodx-bug2.x)**2 + (foody-bug2.y)**2):
        target = bug1

    if ((foodx-bug1.x)**2 + (foody-bug1.y)**2) > ((foodx-bug2.x)**2 + (foody-bug2.y)**2):
	target = bug2
# --- *** What to do if bug reaches edge of map *** ---
    if (bug1.x < -map_size):
        bug1.x = -map_size
        bug1.angle = pi - bug1.angle
    if (bug1.x > map_size):
	bug1.x = map_size
	bug1.angle = pi - bug1.angle
    if (bug1.y < -map_size):
	bug1.y = -map_size
	bug1.angle = -bug1.angle
    if (bug1.y > map_size):
	bug1.y = map_size
	bug1.angle = -bug1.angle
	
    if (bug2.x < -map_size):
        bug2.x = -map_size
        bug2.angle = pi - bug2.angle
    if (bug2.x > map_size):
	bug2.x = map_size
	bug2.angle = pi - bug2.angle
    if (bug2.y < -map_size):
	bug2.y = -map_size
	bug2.angle = -bug2.angle
    if (bug2.y > map_size):
	bug2.y = map_size
	bug2.angle = -bug2.angle
	
    if (bugA.x < -map_size):
        bugA.x = -map_size
        bugA.angle = pi - bugA.angle
    if (bugA.x > map_size):
	bugA.x = map_size
	bugA.angle = pi - bugA.angle
    if (bugA.y < -map_size):
	bugA.y = -map_size
	bugA.angle = -bugA.angle
    if (bugA.y > map_size):
	bugA.y = map_size
	bugA.angle = -bugA.angle

# --- *** Update where the food is relative to the bug *** ---
    sr.foodx = foodx
    sr.foody = foody
    sl.foodx = foodx
    sl.foody = foody
    
    sr2.foodx = foodx
    sr2.foody = foody
    sl2.foodx = foodx
    sl2.foody = foody
    
    srA.foodx = foodx
    srA.foody = foody
    slA.foodx = foodx
    slA.foody = foody

# --- *** Update the plots so the bug and food move *** ---
@network_operation(dt=1*ms)
def update_plot():
    global foodx, foody, bug1_plot, bug2_plot, bugA_plot, food_plot, sr1_plot, sl1_plot, sr2_plot, sl2_plot, srA_plot, slA_plot
    bug1_plot[0].remove()
    bug2_plot[0].remove()
    bugA_plot[0].remove()
    food_plot[0].remove()
    sr1_plot[0].remove()
    sl1_plot[0].remove()
    sr2_plot[0].remove()
    sl2_plot[0].remove()
    srA_plot[0].remove()
    slA_plot[0].remove()
    
#   --- *** Bug 1 *** ---    
    bug1_x_coords = [bug1.x, bug1.x-2*cos(bug1.angle), bug1.x-4*cos(bug1.angle)]    # ant-like body
    bug1_y_coords = [bug1.y, bug1.y-2*sin(bug1.angle), bug1.y-4*sin(bug1.angle)]
    bug1_plot = plot(bug1_x_coords, bug1_y_coords, 'ko')     # Plot the bug's current position
    sr1_plot = plot([bug1.x, sr.x], [bug1.y, sr.y], 'b')
    sl1_plot = plot([bug1.x, sl.x], [bug1.y, sl.y], 'r')
    
#   --- *** Bug 2 *** ---
    bug2_x_coords = [bug2.x, bug2.x-2*cos(bug2.angle), bug2.x-4*cos(bug2.angle)]    # ant-like body
    bug2_y_coords = [bug2.y, bug2.y-2*sin(bug2.angle), bug2.y-4*sin(bug2.angle)]
    bug2_plot = plot(bug2_x_coords, bug2_y_coords, 'ko')     # Plot the bug's current position
    sr2_plot = plot([bug2.x, sr2.x], [bug2.y, sr2.y], 'b')
    sl2_plot = plot([bug2.x, sl2.x], [bug2.y, sl2.y], 'r')
    
#   --- *** Alpha Bug *** ---
    bugA_x_coords = [bugA.x, bugA.x-2*cos(bugA.angle), bugA.x-4*cos(bugA.angle)]    # ant-like body
    bugA_y_coords = [bugA.y, bugA.y-2*sin(bugA.angle), bugA.y-4*sin(bugA.angle)]
    bugA_plot = plot(bugA_x_coords, bugA_y_coords, 'ko')     # Plot the bug's current position
    srA_plot = plot([bugA.x, srA.x], [bugA.y, srA.y], 'b')
    slA_plot = plot([bugA.x, slA.x], [bugA.y, slA.y], 'r')

    food_plot = plot(foodx, foody, 'b*')
    axis([-100,100,-100,100])
    draw()
    #print "."
    pause(0.01)


run(5000*ms,report='text')

#figure(2)
#subplot(2,1,1)
#plot(MSL.t/ms, MSL.vL[0],'r')
#subplot(2,1,2)
#plot(MSR.t/ms, MSR.vR[0],'b')
#title('Sensor')
#figure(3)
#subplot(2,1,1)
#plot(MML.t/ms, MML.vmL[0],'r')
#subplot(2,1,2)
#plot(MMR.t/ms, MMR.vmR[0],'b')
#title('Motoneuron')
#figure(3)
#plot(MBL.t/ms, MBL.vl[0],'r')
#plot(MBR.t/ms, MBR.vr[0],'b')
#title('Motor')
#show()


#plt.clf()
#plt.plot(MBR.x[0], MBR.y[0])
#plt.plot(foodx, foody, 'b*')
#axis([-100,100,-100,100])
#title('Path')
#show()





    