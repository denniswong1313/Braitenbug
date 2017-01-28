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
tau = 10*ms

# You will enter your Izhikevich neuron eqns here- you can also introduce the input current that is a function of distance
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

'''

# for each neuron group you will also need to initialize parameters associated with membrane model, and assign thresholds and resets as needed- sveral neurons can be in each group
sr = NeuronGroup(1, sensor_eqsR,
                    threshold =
                    'vR > 30*mV',
                    reset = '''
                    vR = c
                    uR+=d
                    ''',
                    clock=Clock(0.2*ms))
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
                    clock=Clock(0.2*ms))
sl.x_disp = -5
sl.y_disp = 5
sl.x = sl.x_disp
sl.y = sl.y_disp
sl.foodx = foodx
sl.foody = foody
sl.vL = c
sl.uL = b*c

moto_eqsL = '''
dvmL/dt = (0.04/ms/mV)*vmL**2 + (5/ms)*vmL + (140*mV/ms) - umL + imL  :  volt
dumL/dt = a*(b*vmL - umL): volt / second
dgL/dt = -gL/tau + zL/ms : siemens
dzL/dt = -zL/tau : siemens
imL= -gL*(vmL)/(nsiemens*ms) : volt/second
'''

moto_eqsR = '''
dvmR/dt = (0.04/ms/mV)*vmR**2 + (5/ms)*vmR + (140*mV/ms) - umR + imR  :  volt
dumR/dt = a*(b*vmR - umR): volt / second
dgR/dt = -gR/tau + zR/ms : siemens
dzR/dt = -zR/tau :siemens
imR= -gR*(vmR)/(nsiemens*ms) : volt/second
'''

mL = NeuronGroup(1,moto_eqsL,
                    threshold=
                    'vmL>30*mV',
                    reset=
                    '''
                    vmL = c
                    umL +=d
                    ''',
                    clock = Clock(0.2*ms))
                    
mR = NeuronGroup(1,moto_eqsR,
                    threshold=
                    'vmR>30*mV',
                    reset=
                    '''
                    vmR = c
                    umR +=d
                    ''',
                    clock = Clock(0.2*ms))
                    
mL.vmL = c
mL.umL = b*c
mR.vmR = c
mR.umR = b*c

syn_intr = Synapses(sr, mR, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zR += w
                    ''')
syn_intr.connect(i=[0],j=[0])
syn_intr.w= (0.05*ms/(tau*exp(-1))) * nsiemens

syn_intl = Synapses(sl, mL, 
                    model = '''
                    w:siemens
                    ''',
                    on_pre = '''
                    zL += w
                    ''')
syn_intl.connect(i=[0],j=[0])
syn_intl.w= (0.05*ms/(tau*exp(-1))) * nsiemens
                    
# The virtual bug

taum = 10*ms
base_speed = 1
turn_rate = 8000*Hz
eta = 10
baseline = 150

# the motor equations will need to be differential equations that are incremented each time they receive an input spike from the pre-synaptic neurons

bug_eqs = '''
dvl/dt=  (-vl+baseline)/taum : 1
dvr/dt = (-vr+baseline)/taum : 1
speed = (vl + vr)/2 : 1
dangle/dt = (vr - vl)*turn_rate : 1
dx/dt = speed*cos(angle)*15*Hz : 1
dy/dt = speed*sin(angle)*15*Hz : 1
tl:1
tr:1
'''

bug = NeuronGroup(1, bug_eqs, clock=Clock(0.2*ms))

bug.angle = pi/2
bug.x = 0
bug.y = 0
#bug.vl = 10
#bug.vr = 10

# Synapses (sensors communicate with bug motor)
# You will need to add the Synapse pieces here to connect the sensors with the bug motors 

syn_r = Synapses(mR, bug, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vr -= w
                    ''')
syn_r.connect(i=[0],j=[0])
syn_r.w= (0.05*ms/(taum*exp(-1)))

syn_l = Synapses(mL, bug, 
                    model = '''
                    w:1
                    ''',
                    on_pre = '''
                    vl -= w
                    ''')
syn_l.connect(i=[0],j=[0])
syn_l.w= (0.05*ms/(taum*exp(-1)))

MSR = StateMonitor(sr,('vR'),record=True)
MSL = StateMonitor(sl,('vL'),record=True)

MMR = StateMonitor(mR,('vmR'),record=True)
MML = StateMonitor(mL,('vmL'),record=True)

MBR = StateMonitor(bug,('vr','x','y'),record=True)
MBL = StateMonitor(bug,('vl'),record=True)

f = figure(1)
bug_plot = plot(bug.x, bug.y, 'ko')
food_plot = plot(foodx, foody, 'b*')
sr_plot = plot([0], [0], 'w')   # Just leaving it blank for now
sl_plot = plot([0], [0], 'w')
# Additional update rules (not covered/possible in above eqns)


#This block of code updates the position of the bug and food and makes the appropriate rotations
@network_operation()
def update_positions(dt=1*ms):
    global foodx, foody, food_count

#note that in the version I1 and I2 are computed here and drive the motors directly- you will need to remove and add the currents to the sensor equations and use the neuron outputs to drive the motors.

    sr.iR = I0 / sqrt(((sr.x-foodx)**2+(sr.y-foody)**2)) * mV/ms
    sl.iL = I0 / sqrt(((sl.x-foodx)**2+(sl.y-foody)**2)) * mV/ms
    
    sr.x = bug.x + sr.x_disp*cos(bug.angle-pi/2) - sr.y_disp*sin(bug.angle-pi/2)
    sr.y = bug.y + sr.x_disp*sin(bug.angle-pi/2) + sr.y_disp*cos(bug.angle-pi/2)

    sl.x = bug.x + sl.x_disp*cos(bug.angle-pi/2) - sl.y_disp*sin(bug.angle-pi/2)
    sl.y = bug.y + sl.x_disp*sin(bug.angle-pi/2) + sl.y_disp*cos(bug.angle-pi/2)


    if ((bug.x-foodx)**2+(bug.y-foody)**2) < 16:
	food_count += 1
	foodx = randint(-map_size+20, map_size-20)
	foody = randint(-map_size+20, map_size-20)

    if (bug.x < -map_size):
        bug.x = -map_size
        bug.angle = pi - bug.angle
    if (bug.x > map_size):
	bug.x = map_size
	bug.angle = pi - bug.angle
    if (bug.y < -map_size):
	bug.y = -map_size
	bug.angle = -bug.angle
    if (bug.y > map_size):
	bug.y = map_size
	bug.angle = -bug.angle

    sr.foodx = foodx
    sr.foody = foody
    sl.foodx = foodx
    sl.foody = foody


# this block of code updates the plots so you can see the bug and food move
@network_operation(dt=1*ms)
def update_plot():
    global foodx, foody, bug_plot, food_plot, sr_plot, sl_plot
    bug_plot[0].remove()
    food_plot[0].remove()
    sr_plot[0].remove()
    sl_plot[0].remove()
    bug_x_coords = [bug.x, bug.x-2*cos(bug.angle), bug.x-4*cos(bug.angle)]    # ant-like body
    bug_y_coords = [bug.y, bug.y-2*sin(bug.angle), bug.y-4*sin(bug.angle)]
    bug_plot = plot(bug_x_coords, bug_y_coords, 'ko')     # Plot the bug's current position
    sr_plot = plot([bug.x, sr.x], [bug.y, sr.y], 'b')
    sl_plot = plot([bug.x, sl.x], [bug.y, sl.y], 'r')
    food_plot = plot(foodx, foody, 'b*')
    axis([-100,100,-100,100])
    draw()
    #print "."
    pause(0.01)


run(1000*ms,report='text')

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
#figure(4)
#plot(MBL.t/ms, MBL.vl[0],'r')
#plot(MBR.t/ms, MBR.vr[0],'b')
#title('Motor')
#show()


plt.clf()
plt.plot(MBR.x[0], MBR.y[0])
plt.plot(foodx, foody, 'b*')
axis([-100,100,-100,100])
title('Guardian Bug Path')
show()





    