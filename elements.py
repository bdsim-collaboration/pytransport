import numpy as _np
from scipy import constants as _con
import string as _string
from _General import functions
    
class elements(functions):
    
    def drift(self,line):
        label = self._get_label(line)
        for ele in line[1:]:
            if len(ele) > 0: #I.E. Not a blank space
                endofline = self._endofline(ele)
                if endofline == -1:
                    driftlen = ele
                    break
                elif endofline != -1:
                    driftlen = ele[:endofline]
                    break
        
        if self.units['element_length'] != 'm':
            length_in_metres = _np.float(driftlen) * self.scale[self.units['element_length'][0]]  #Convert to metres.
        else:
            length_in_metres = _np.float(driftlen)

        elementid = 'DR'+_np.str(self.elementprops.drifts)
        self.elementprops.drifts += 1

        self.machine.AddDrift(name=elementid,length=length_in_metres)


        
    def dipole(self,line):        
        label = self._get_label(line)
        dipoledata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    dipoledata.append(_np.float(ele))
                except ValueError:
                    dipoledata.append(ele)
        length = dipoledata[0]          # First two non-blanks must be the entries in a specific order.
        if self.benddef:
            bfield = dipoledata[1]
            field_in_Gauss = bfield * self.scale[self.units['magnetic_fields'][0]]     # Scale to Gauss
            field_in_Tesla = field_in_Gauss * 1e-4                                      # Convert to Tesla
            rho = self.beamprops.brho / (_np.float(field_in_Tesla))             # Calculate bending radius.
            angle = (_np.float(length) / rho) * self.elementprops.bending       # for direction of bend
        elif not self.benddef:
            angle_in_deg = dipoledata[1]
            angle = angle_in_deg * (_np.pi/180.)
        
        if self.units['element_length'] != 'm':
            length_in_metres = length * self.scale[self.units['element_length'][0]]
        else:
            length_in_metres = length
        
        elementid = 'BM'+_np.str(self.elementprops.dipoles)
        self.elementprops.dipoles += 1
        
        self.machine.AddDipole(name=elementid,length=length_in_metres,angle=_np.round(angle,4))
        



    def change_bend(self,line):
        '''Function to change the direction of the dipole bend. Can be a direction other than horizontal (i.e != n*pi).
            '''
        ## NOT FULLY TESTED.
        angle = 0
        for index,ele in enumerate(line[1:]): #For loops are iterating over blank space (delimiter)
            if ele != ' ':
                endofline = self._endofline(ele)
                if len(ele) > 0 and endofline == -1: #I.E endofline is not in this element
                    angle = ele
                    break
                elif len(ele) > 0 and endofline != -1: #endofline is in this element
                    angle = _np.str(ele[:endofline])
                    break
        self.elementprops.angle = _np.float(angle)
        if self.elementprops.angle >= 360:
            self.elementprops.angle = _np.mod(self.elementprops.angle,360)
        if self.elementprops.angle <= -360:
            self.elementprops.angle = _np.mod(self.elementprops.angle,-360)

        if self.elementprops.angle == 180 or self.elementprops.angle == -180: #If 180 degrees, switch bending angle
            self.elementprops.bending *= -1
        elif self.elementprops.angle != 0:                        #If not 180 degrees, use transform3d.      
            self.elementprops.angle *= -1                         #For conversion to correct direction. Eg in TRANSPORT -90 is upwards, in BDSIM, 90 is upwards.  
            anginrad = self.elementprops.angle * (_np.pi / 180)
            elementid = 't'+_np.str(self.elementprops.transforms)
            self.elementprops.transforms += 1
            self.machine.AddTransform3D(name=elementid,psi=anginrad)
        self.angle=0    




    def quadrupole(self,line):
        label = self._get_label(line)
        quaddata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    quaddata.append(_np.float(ele))
                except ValueError:
                    quaddata.append(ele)
        length = quaddata[0]        # First three non-blanks must be the entries in a specific order.
        field_at_tip = quaddata[1]  # Field in TRANSPORT units 
        pipe_rad = quaddata[2]      # Pipe Radius In TRANSPORT units
        
        field_in_Gauss = field_at_tip * self.scale[self.units['magnetic_fields'][0]] #Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  #Convert to Tesla
        
        if self.units['bend_vert_gap'] != 'm':
            pipe_in_metres = pipe_rad * self.scale[self.units['bend_vert_gap'][0]]  #Scale to meters
        else:
            pipe_in_metres = pipe_rad
    
        if self.units['element_length'] != 'm':
            length_in_metres = length * self.scale[self.units['element_length'][0]] #Scale to meters
        else:
            length_in_metres = length
        
        field_gradient = (field_in_Tesla / pipe_in_metres) / self.beamprops.brho    #K1 in correct units
        
        if label is not None: #Write to file
            if field_gradient > 0:
                elementid = 'QF'+_np.str(self.elementprops.quads)
            elif field_gradient < 0:
                elementid = 'QD'+_np.str(self.elementprops.quads)
            else:
                elementid = 'NULLQUAD'+_np.str(self.elementprops.quads)  #For K1 = 0. 
        
        self.elementprops.quads += 1

        self.machine.AddQuadrupole(name=elementid,length=length_in_metres,k1=_np.round(field_gradient,4))



    def collimator(self,line):
        ### Was used to write the location of a collimator as a string, redundant if file writing done with pybdsim.
        label = self._get_label(line)
        collstarted = False
        for index in self._collindex: #Look for existing collimator elements of the same name
            if index == label:
                collstarted = True      #If one already exists, that must be the start of the collimator
                break
        
        if collstarted == True:          
            #If the start already exists, input line must be for collimator end
            coll = 'ends'
        else:
            coll = 'starts'
        self._collindex.append(label) #Add to collimator list

        collidata = []
        for index,ele in enumerate(line[1:]): # Iterate over line to extract data
            if ele != '':
                try:
                    collidata.append(_np.float(ele))
                except ValueError:
                    dummy=1

        horwidth = 'Unknown'
        verwidth = 'Unknown'
        # Determine which entry is for horiz. and vert.
        if collidata[0] == 1.0:
            horwidth = _np.str(collidata[1])
        elif collidata[0] == 3.0:
            verwidth = _np.str(collidata[1])

        if len(collidata) > 2:
            if collidata[2] == 1.0:
                horwidth = _np.str(collidata[3])
            elif collidata[2] == 3.0:
                verwidth = _np.str(collidata[3])

        collline  = '! A collimator labelled ' + label +' ' + coll + ' here'
        collline2 = '! with slit size half widths of x = '+horwidth+' '+self.units['x']+' and y = '+verwidth+' '+self.units['x']+'.'




    def accelerator(self,line):
        ''' A Function that writes the properties of an acceleration element
            as a gmad comment. Assumed only 1 acceleration element in the machine.
            '''
        # Redundant function until comments and /or acceleration components can be handled
        
        label = self._get_label(line)
        accdata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    accdata.append(_np.float(ele))
                except ValueError:
                    dummy = 1
        if len(accdata) == 2:       # Newer case with multiple elements
            self._acc_sequence(line)
        elif len(accdata) == 4:     # Older case for single element
            acclen = accdata[0]
            e_gain = accdata[1]
            phase_lag = accdata[2]
            wavel = accdata[3]
            
            #Write to file
            accline =  '! An accelerator element goes here of length '+_np.str(acclen)+' '+self.units['element_length']+', \n'
            accline2 = '! with an energy gain of '+_np.str(e_gain)+' '+self.units['p_egain']+', phase lag of '+_np.str(phase_lag)+' degrees, \n'
            accline3 = '! and a wavelength of '+_np.str(wavel)+' '+_self.units['bunch_length']+'. \n'


    def _acc_sequence(self,inputline):
        '''Function to calulate the total length of a sequence of accelerator components.
            '''
        ## UNTESTED ##
        
        # Redundant function until comments and /or acceleration components can be handled.
        concat=''
        for ele in inputline:                       #Concat string back together
            concat += ele
            if ele != inputline[-1]:
                concat += ' '
        ele=0
        for linenum,line in enumerate(self.data):   #Find start of accelerator element sequence
            if line == concat:
                seq_start = linenum
                break
        accparts=[]
        accelements = 0
        isacc = True
        while isacc:                                #Find remainder of acc sequence
            a = self.data[seq_start+accelements].split(' ')
            if a[0] == '11.':
                accparts.append(a)
                accelements += 1
            else:
                isacc = False
        acccopy = accparts #temp copy of data
        accparts=[]
        for part in acccopy:               #Calculate total length of accelerator part.
            accdata=[]
            for index,ele in enumerate(part[1:]):
                if ele != '':
                    try:
                        accdata.append(_np.float(ele))
                    except ValueError:
                        pass
            accparts.append(accdata)
        accarray = _np.array(accparts)
        tot_len = _np.sum(accarray[2:,0])
        #Write to file
        accline = '! An electrostatic accelerator section goes here, of length '+_np.str(tot_len)+' '+self.units['element_length']+','
        accline2= '!split into '+_np.str(accelements-1)+' elements, with a total voltage of '+_np.str(_np.round(accarray[0,1],3))+' GV.'



    def sextupole(self,line):
        label = self._get_label(line)
        sextudata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    sextudata.append(_np.float(ele))
                except ValueError:
                    sextudata.append(ele)
        length = sextudata[0]        # First three non-blanks must be the entries in a specific order.
        field_at_tip = sextudata[1]  # Field in TRANSPORT units
        pipe_rad = sextudata[2]      # Pipe Radius In TRANSPORT units
        
        field_in_Gauss = field_at_tip * self.scale[self.units['magnetic_fields'][0]] #Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  #Convert to Tesla
        
        if self.units['bend_vert_gap'] != 'm':
            pipe_in_metres = pipe_rad * self.scale[self.units['bend_vert_gap'][0]]  #Scale to meters
        else:
            pipe_in_metres = pipe_rad
        
        if self.units['element_length'] != 'm':
            length_in_metres = length * self.scale[self.units['element_length'][0]] #Scale to meters
        else:
            length_in_metres = length
        
        field_gradient = (2*field_in_Tesla / pipe_in_metres**2) / self.beamprops.brho    #K2 in correct units
        
        elementid = 'SEXT'+_np.str(self.elementprops.sextus)
        
        self.elementprops.sextus += 1
        
        self.machine.AddSextupole(name=elementid,length=length_in_metres,k2=_np.round(field_gradient,4))



    def solenoid(self,line):
        label = self._get_label(line)
        soledata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    soledata.append(_np.float(ele))
                except ValueError:
                    soledata.append(ele)
        length = soledata[0]        # First three non-blanks must be the entries in a specific order.
        field = soledata[1]         # Field in TRANSPORT units
        
        field_in_Gauss = field * self.scale[self.units['magnetic_fields'][0]] #Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  #Convert to Tesla
        
        if self.units['element_length'] != 'm':
            length_in_metres = length * self.scale[self.units['element_length'][0]] #Scale to meters
        else:
            length_in_metres = length
                
        elementid = 'SOLE'+_np.str(self.elementprops.solenoids)

        self.elementprops.solenoids += 1
        
        self.machine.AddSolenoid(name=elementid,length=length_in_metres,ks=_np.round(field_in_Tesla,4))



    def printline(self,line):
        label = self._get_label(line)
        for ele in line[1:]:
            if ele == '48.':
                self.benddef = False
            if ele == '47.':
                self.benddef = True

