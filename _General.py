import numpy as _np
from scipy import constants as _con
import string as _string


class functions():
    def _is_addition(self,line):
        '''Function to check if a BEAM line of TRANSPORT code is a beam definition or r.m.s addition.
            '''
#        concat=''
#        for ele in line:                       #Concat string back together
#            concat += ele
#            if ele != line[-1]:
#                concat += ' '
#        ele=0
#
#        endindex = self._endofline(line)
#        if concat[endindex-1] == '0' or concat[endindex-2] == '0':
#            return True
#        else:
#            return False
        if len(line) > 8:
            if line[8] == '0.':
                return True
            else:
                return False
        else:
            print('Incorrect number of entries for a beam definition')
            return None
    
    def _is_sentinel(self,line):
        for element in line:
            if element[:8] == 'SENTINEL':
                return True
        return False


    def _facerotation(self,line,linenum):
        faceline = self.data[linenum]
        splitline = faceline.split(' ')
        splitline = self._remove_label(splitline)
        splitline = self._remove_spaces(splitline)
        try:
            if _np.float(splitline[0]) == 2.0:
                endof = self._endofline(splitline[1])
                if endof != -1:
                    angle = _np.round(_np.float(splitline[1][:endof]),4)
                else:
                    angle = _np.round(_np.float(splitline[1]),4)
            else:
                angle = 0
        except ValueError:
            angle = 0

        return angle
    
    def _remove_spaces(self,line):
        elementlist=[]
        for i in line:
            if (i != '') and (i != ' '):
                elementlist.append(i)
        line = _np.array(elementlist)
        return line


    def _endofline(self,line):   #Find the end of the line of code
        endpos = -1
        breakloop=False
        if isinstance(line,_np.str):
            for charnum,char in enumerate(line):
                if char == ';':
                    endpos = charnum
                    break
        elif isinstance(line,_np.ndarray):
            for index,ele in enumerate(line):
                for char in ele:
                    if char == ';':
                        endpos = index
                        breakloop=True
                        break
                if breakloop:
                    break
        return endpos
                    


    def _get_label(self,line):
        '''Function to get element label from code line.'''
        for elenum,ele in enumerate(line):
            startslash = _string.find(ele,"/")
            startquote = _string.find(ele,"'")
            startequal = _string.find(ele,"=")
            if startslash != -1:
                end = 1 + startslash + _string.find(ele[(startslash+1):],"/")
                label = ele[startslash+1:end]
                break
            elif startquote != -1:
                end = 1 + startquote + _string.find(ele[(startslash+1):],"'")
                label = ele[startslash+1:end]
                break
            elif startequal != -1:
                end = 1 + startequal + _string.find(ele[(startslash+1):],"=")
                label = ele[startslash+1:end]
                break
            else:
                label = None
        return label,elenum
            
            
    def _remove_label(self,line):
        '''Function to remove the label from a line.
            '''
        label,elenum = self._get_label(line)
        if label != None:
            element = line[elenum]
            lablen = len(label)
            newval = element
            for index in range(len(element)):
                if element[index:index+lablen] == label:
                    prelabel = element[:index-1]
                    postlabel = element[index+lablen+1:]
                    newval = prelabel+' '+postlabel
                    break
            line[elenum] = newval 
        return line 
            


    def _get_indic(self):
        '''Function to read the indicator number. Must be 0, 1, or 2, where:
                0 is a new lattice
                1 is for fitting with the first lattice (from a 0 indicator file)
                2 if for a second fitting which suppresses the first fitting.
            '''
        for linenum,line in enumerate(self.data):
            if line[0] == '0':
                if line[1] == '\r' or '\n':
                    indc = 0
                    break
            if line[0] == '1':
                if line[1] == '\r' or '\n':
                    indc = 1    
                    break
            if line[0] == '2':
                if line[1] == '\r' or '\n':
                    indc = 2    
                    break
        return indc,linenum


    def _get_comment(self,line):
        '''Function to extract a comment from a line. 
            Will only find one comment per line.'''
        concat=''
        for ele in line:
            concat += ele
            concat += ' '
        commstart = _string.find(concat,'(')
        commend = _string.find(concat,')')
        if commstart != -1 and commend != -1:
            comment = concat[commstart+1 : commend]
            gmadcomment = '! '+comment
        else:
            gmadcomment = None
        return gmadcomment
    
    def _bunch_length_convert(self,bunch_length):
        ### Function to convert bunch length unit in TRANSPORT into seconds.
        scale = self.scale[self.units['bunch_length'][0]]   
        blmeters = bunch_length * scale     #Bunch length scaled to metres
        blseconds = blmeters / (self.beamprops.beta*_con.c)       #Length converted to seconds
        return blseconds

    def unit_change(self,line):
        '''Function to change the units (scaling) of various parameters'''
        label,elenum = self._get_label(line)
        
        if label=='CM' or label=='MM' or label=='UM' or label=='NM':
            label = label.lower()
        ### Convert Energy Unit Cases:
        if label=='EV':
            label='eV'
        if label=='KEV':
            label='keV'
        if label=='MEV':
            label='MeV'
        if label=='GEV':
            label='GeV'
        if label=='TEV':
            label='TeV'
        
        if line[1] == '1.0':    #Horizontal and vertical beam size
            self.units['x'] = label
            self.units['y'] = label
        if line[1] == '2.0':    #Horizontal and vertical divergence
            self.units['xp'] = label
            self.units['yp'] = label
        if line[1] == '3.0':    #Bending Magnet Gap
            self.units['bend_vert_gap'] = label
        if line[1] == '4.0':    #Vertical Divergence ONLY
            self.units['yp'] = label
        if line[1] == '5.0':    #Pulsed Beam Length 
            self.units['bunch_length'] = label
        if line[1] == '6.0':    #Momentum Spread
            self.units['momentum_spread'] = label   ## Percent
        if line[1] == '7.0':    #Bend/pole face rotation
            pass
        if line[1] == '8.0':    #Element Length
            self.units['element_length'] = label
        if line[1] == '9.0':    #Magnetic Field
            self.units['magnetic_fields'] = label
        if line[1] == '10.0':   #Mass
            print('Cannot change mass scale.')
        if line[1] == '11.0':   #Momentum / energy gain during acc.
            self.units['p_egain'] = label



    def _calculate_energy(self,momentum):
        '''Function to calculate:
                Total Energy
                Kinetic Energy
                Momentum
                Lorentz factor (gamma)
                Velocity (beta)
                Magnetic rigidity (brho)
            '''
        
        momentum = _np.float(momentum)
        self.beamprops.momentum = momentum
        p_mass = self.beamprops.mass ## Particle rest mass (in GeV)
        
        mom_unit = self.units['p_egain']
        if mom_unit != 'eV':
            scaling = 1e9 / self.scale[mom_unit[0]]     ## Scaling relative to mom. unit
            mom_in_ev = momentum * self.scale[mom_unit[0]]
        elif mom_unit == 'eV':
            scaling = 1e9                               ## Scaling relative to 1 eV
            mom_in_ev = momentum
        p_mass *= scaling                               ## Scale particle rest mass
        energy = _np.sqrt((p_mass**2 * _con.c**2) + (momentum**2 * _con.c**2)) / _con.c
        self.beamprops.tot_energy = energy
        self.beamprops.k_energy = energy - p_mass
        self.beamprops.gamma = energy / p_mass
        self.beamprops.beta = _np.sqrt((1.0 - (1.0 / self.beamprops.gamma**2)))
        self.beamprops.brho = mom_in_ev / _con.c
    
  
