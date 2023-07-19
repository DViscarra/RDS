import matplotlib.pyplot as plt
import numpy as np
from random import choice

class chemical:
    default_id = 64

    def __init__(self,diff_const=1,id="",color=""):
        if id == "":
            chemical.default_id = chemical.default_id + 1
            self.id = chr(chemical.default_id)
        else:
            self.id = id
        self.diff_const = diff_const

        if color == "":
            self.color = '#' + ''.join([choice('0123456789ABCDEF') for j in range(6)])
        else:
            self.color = color

    def __repr__(self):
        return self.id

class reaction:
    def __init__(self,chems_in:dict,chems_out:dict,reaction_const=1):
        self.chems_in = chems_in.copy()
        self.chems_out = chems_out.copy()
        self.chems_list = set(list(chems_in) + list(chems_out))
        self.reaction_const = reaction_const

        self.id = ""
        first = True
        if self.chems_in:
            for chem in self.chems_in:
                if first:
                    self.id = self.id + str(self.chems_in[chem]) + chem.id + " "
                    first = False
                else:
                    self.id = self.id + "+ " + str(self.chems_in[chem]) + chem.id + " "
        else:
            self.id = self.id + "Ø "

        self.id = self.id + "-(" + str(reaction_const) + ")> "

        first = True
        if self.chems_out:
            for chem in self.chems_out:
                if first:
                    self.id = self.id + str(self.chems_out[chem]) + chem.id + " "
                    first = False
                else:
                    self.id = self.id + "+ " + str(self.chems_out[chem]) + chem.id + " "
        else:
            self.id = self.id + "Ø"

    def __repr__(self):
        return self.id

    def has_chem(self,chem:chemical):
        return chem in self.chems_list

    def pred_prey(self,chems_val:dict):
        c = self.reaction_const
        for chem in self.chems_in:
            c = c * chems_val[chem]**self.chems_in[chem]
    
        # consumption
        consumption = {}
        # if not self.chems_out:
        #     for chem in self.chems_in:
        #         consumption[chem] = self.reaction_const*self.chems_in[chem]
        for chem in self.chems_list:
            if chem in self.chems_in:
                if chem in self.chems_out:
                    if self.chems_in[chem] > self.chems_out[chem]:
                        val = self.chems_in[chem] - self.chems_out[chem]
                    else:
                        val = 0
                else:
                    val = self.chems_in[chem]
            else:
                val = 0
            consumption[chem] = c * val
            
        # creation
        creation = {}
        if not self.chems_in:
            for chem in self.chems_out:
                creation[chem] = self.reaction_const*self.chems_out[chem]
        for chem in self.chems_list:
            if chem in self.chems_out:
                if chem in self.chems_in:
                    if self.chems_in[chem] < self.chems_out[chem]:
                        val = self.chems_out[chem] - self.chems_in[chem]
                    else:
                        val = 0
                else:
                    val = self.chems_out[chem]
            else:
                val = 0
            creation[chem] = c * val

        return consumption, creation

class rod:
    rod_num = 0

    def __init__(self,dt,dx,length=1,temp=1):
        self.chems_list = set()
        self.chems = {}
        self.chems_next = {}
        self.reacts = []
        self.dx = dx
        self.dt = dt
        self.nx = int(length/dx)
        self.length = length
        self.temp = temp
        self.coef = dt/(dx**2)
        self.x = []
        for i in range(self.nx):
            self.x.append(i*self.dx)

        self.fig, self.ax = plt.subplots()
        rod.rod_num += 1
        self.fig.canvas.manager.set_window_title("Rod " + str(rod.rod_num))
        plt.ion()

    def __repr__(self):
        out = "=================================\nChemicals\n=================================\n"
        # chemicals
        for chem in self.chems_list:
            out = out + chem.id + " "
            for i in self.chems[chem]:
                out = out + str(i) +" "
            out = out + "\n"

        out = out + "=================================\nReactions\n=================================\n"
        #reactions
        for react in self.reacts:
            out = out + react.id + "\n"

        return out

    def add_reaction(self,reacts):
        if type(reacts) == reaction:
            self.reacts.append(reacts)
        elif type(reacts) == list:
            self.reacts = self.reacts + reacts

    def add_chem(self,chem:chemical,region:list,vals:list):
        if len(region) != len(vals):
            raise Exception("Provide region and vals of the same size")
        
        if not chem in self.chems_list:
            self.chems_list.add(chem)
            self.chems[chem] = [0] * self.nx
            self.chems_next[chem] = [0] * self.nx

        for i in range(len(region)):
            self.chems[chem][region[i]] = self.chems[chem][region[i]] + vals[i]
        
    def transfer_to(self,space,chem,range_in=[],range_out=[]):
        # to do: add surface transfer support
        
        if range_in == []:
            range_in = list(range(self.nx))

        if range_out == []:
            range_out = range_in.copy()

        transfer_val = 0
        for i in range(len(range_in)):
            transfer_val = transfer_val + self.chems[chem][i]
            self.chems[chem][range_in[i]] = 0

        transfer_val = [transfer_val/len(range_out)] * len(range_out)
        space.add_chem(chem,range_out,transfer_val)

    def integrate_chems(self,chems=None) -> dict:
        integrated = {}
        if chems == None or type(chems) == list:
            if chems == None:
                chems = self.chems_list

            for chem in chems:
                integrated[chem] = 0
                for i in self.chems[chem]:
                    integrated[chem] = integrated[chem] + i*self.dx

        elif type(chems) == chemical:
            chem = chems
            integrated[chem] = 0
            for i in self.chems[chem]:
                integrated[chem] = integrated[chem] + i*self.dx

        else:
            raise Exception("Give no argumant for all chemicals, list of chemicals wanted, or single chemical")
        
        return integrated

    def update(self):
        for i in range(1,self.nx-1):
            cons = {}
            cre = {}
            for react in self.reacts:
                if set(react.chems_list).issubset(set(self.chems_list)):
                    chems_val = {}
                    for chem in react.chems_in:
                        chems_val[chem] = self.chems[chem][i]
                    cons_temp, cre_temp = react.pred_prey(chems_val)
                    for j in cons_temp:
                        if j in cons:
                            cons[j] = cons[j] + cons_temp[j]
                            cre[j] = cre[j] + cre_temp[j]
                        else:
                            cons[j] = cons_temp[j]
                            cre[j] = cre_temp[j]
                else:
                    raise Exception("Missing chemicals in rod for reaction")
            for chem in self.chems_list:
                if not chem in cons:
                    cons[chem] = 0
                if not chem in cre:
                    cre[chem] = 0
                
                self.chems_next[chem][i] = self.coef*chem.diff_const*(self.chems[chem][i-1] - 2*self.chems[chem][i] + self.chems[chem][i+1])\
                    + self.chems[chem][i]\
                    + self.dt*(cre[chem]- cons[chem])

        for chem in self.chems_list:
            self.chems_next[chem][0] = self.chems_next[chem][1]
            self.chems_next[chem][self.nx-1] = self.chems_next[chem][self.nx-2]

        for chem in self.chems_list:
            self.chems[chem] = self.chems_next[chem].copy()

    def show(self,time=None,ylim=True,pause=0.00000001):
        self.ax.clear()
        if time != None:
            self.fig.suptitle("Time: " + str(time))
        else:
            self.fig.suptitle("")
        self.ax.set_xlabel("x")
        self.ax.set_ylabel("Chemical Concentration")
        for chem in self.chems_list:
            self.ax.plot(self.x,self.chems[chem],color=chem.color,label=repr(chem))
        self.fig.legend()
        if type(ylim) != bool:
            self.ax.set_ylim(ylim[0],ylim[1])
        self.fig.canvas.draw()
        self.fig.show()
        plt.pause(pause)

class surface:
    surface_num = 0

    def __init__(self,dx,dy,dt,length_x=1,length_y=1,temp=1):
        self.chems_list = set()
        self.chems = {}
        self.chems_next = {}
        self.reacts = []
        self.dx = dx
        self.dy = dy
        self.dt = dt
        self.nx = int(length_x/dx)
        self.ny = int(length_y/dy)
        self.length_x = length_x
        self.length_y = length_y
        self.temp = temp 
        self.coef = None #####################################################
        self.x = np.arange(0,1,self.dx)
        self.y = np.arange(0,1,self.dy)
        self.x, self.y = np.meshgrid(self.x,self.y)

        self.fig, self.ax = plt.subplots()
        surface.surface_num += 1
        self.fig.canvas.manager.set_window_title("Surface " + str(surface.rod_num))
        plt.ion()

    def __repr__(self):
        pass

    def add_reaction(self,react:reaction):
        self.reacts.append(react)

    def add_chem(self,chem:chemical,region:list,vals:list):
        if len(region) != len(vals):
            raise Exception("Provide region and vals of the same size")
        
        if not chem in self.chems_list:
            self.chems_list.add(chem)
            self.chems[chem] = [0] * self.nx
            self.chems_next[chem] = [0] * self.nx

        for i in range(len(region)):
            self.chems[chem][region[i]] = self.chems[chem][region[i]] + vals[i]

    def transfer_to(self,space,chem,range_in=[],range_out=[]):
        pass

    def integrate_chems(self,chems=None) -> dict:
        pass

    def update(self):
        pass

    def show(self,time=None):
        pass