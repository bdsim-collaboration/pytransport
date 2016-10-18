import numpy as _np
from scipy import constants as _con
import string as _string
import glob as _glob
import os as _os
import reader as _reader
import sys

class functions():
    def _is_addition(self,line,type='input'):
        '''Function to check if a BEAM line of TRANSPORT code is a beam definition or r.m.s addition.
            '''
        ## Output file is a standard format, any RMS addition line should always be 10 long.
        if type == 'output':
            if len(line) == 10:
                return True
    
        elif type == 'input':
            if len(line) > 8:
                if (line[8] == '0.') or (line[8] == '0'):
                    return True
                else:
                    return False
        else:
            raise ValueError("File type can only be input or output")

    
    def _is_sentinel(self,line):
        for element in line:
            if element[:8] == 'SENTINEL':
                return True
        return False


    def _facerotation(self,line,linenum):
        anglein=0
        angleout=0
        
        linelist=[]
        for i in self.data[(linenum-5):linenum]:
            linelist.append(i)
        linelist.reverse()   #Search for poleface in reverse line order

        for line in linelist:
            try:
                if _np.float(line[0]) == 4.0:
                    break
                elif _np.float(line[0]) == 2.0:
                    endof = self._endofline(line[1])
                    if endof != -1:
                        anglein = _np.round(_np.float(line[1][:endof]),4)
                    else:
                        anglein = _np.round(_np.float(line[1]),4)
                    break
                else:
                    pass
            except ValueError:
                pass
        for line in self.data[linenum+1:(linenum+6)]:
            try:
                if _np.float(line[0]) == 4.0:
                    break
                elif _np.float(line[0]) == 2.0:
                    endof = self._endofline(line[1])
                    if endof != -1:
                        angleout = _np.round(_np.float(line[1][:endof]),4)
                    else:
                        angleout = _np.round(_np.float(line[1]),4)
                    break
                else:
                    pass
            except ValueError:
                pass
        return anglein,angleout
    
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
                    
    def _load_file(self,input):
        '''Load file to be converted into gmad format. 
           Some processing here too (removal of blank spaces etc)
            '''
        data=[]
        filedata=[]
        temp = _reader.reader()

        if not isinstance(input, _np.str):
            raise TypeError("Input must be a string")
        
        infile = input.split('/')[-1] #Remove filepath, leave just filename
        self._file = infile[:-4]         #Remove extension
        f = open(input)
        try:
            numlines=0
            for inputline in f:
                # If a line in the file is equal to one of these, then it's
                # likely to be an output file, so instead, use the reader
                if (inputline == '0    0\n') or (inputline == '0    0\r\n'):
                    #reset lists
                    data = []
                    filedata = []
                    flist = temp._file_to_list(input)
                    lattice,output=temp._get_latticeandoutput(flist)
                    for linenum, latticeline in enumerate(lattice):
                        latticeline = latticeline.replace(';','')
                        line = _np.array(latticeline.split(' '),dtype=_np.str)
                        line = self._remove_blanks(line)
                        
                        # Method of dealing with split lines in the output
                        # Should only be applicable to type 12 entry (up to 15 variables)
                        # It is assumed that the line is always split, so be careful.
                        prevline = lattice[linenum-1].replace(';','')
                        prevline = _np.array(prevline.split(' '),dtype=_np.str)
                        prevline = self._remove_blanks(prevline)
                        
                        #Catch for any lines that are comments
                        try:
                            if (linenum > 0) and _np.abs(_np.float(line[0])) == 12.0:
                                latticeline, line = self._joinsplitlines(linenum,lattice)
                            # Ignore line after type 12 entry (second part of split line)
                            if (linenum > 1) and _np.abs(_np.float(prevline[0])) == 12.0:
                                pass
                            else:
                                data.append(line)
                                filedata.append(latticeline)
                        except ValueError:
                            data.append(line)
                            filedata.append(latticeline)
                    break
                else:
                    endoflinepos = self._endofline(inputline)
                    templine = inputline
                    if endoflinepos > 0:
                        templine = inputline[:endoflinepos]
                    line = _np.array(templine.split(' '),dtype=_np.str)
                    line = self._remove_blanks(line)
                    data.append(line)
                    filedata.append(inputline)
            self._fileloaded = True
        except IOError:
            print 'Cannot open file.'
        f.close()
        self.data=data
        self.filedata=filedata


    def _joinsplitlines(self,linenum,lattice):
        firstline = lattice[linenum].replace(';','')
        latticeline = firstline #Copy for later
        firstline = _np.array(firstline.split(' '),dtype=_np.str)
        firstline = self._remove_blanks(firstline)
        numericals = []
        
        #Keep entries that are strings of numbers
        for i in firstline:
            try:
                number = _np.float(i)
                numericals.append(_np.str(number))
            except ValueError:
                pass
    
        #Number of numerical elements minus the first which should be the entry type number.
        #This is bascially a way of extracting any label or comments.
        numelements = len(numericals) - 1

        secline = lattice[linenum+1].replace(';','')
        secline = _np.array(secline.split(' '),dtype=_np.str)
        secline = self._remove_blanks(secline)
        secnumericals = []

        for i in secline:
            try:
                number = _np.float(i)
                secnumericals.append("%.4f" %number)
            except ValueError:
                pass

        #Second line should be 15 minus number of numerical elements from prev line.
        #This is done to skip erroneous numbers in the line such as '000' which have
        #appeared when lines have been split.
        secline = secnumericals[-15+numelements:]
        numericals.extend(secline)

        #Add to latticeline so as to appear like one single line in the file
        seclinetxt=""
        for i in secline:
            newline = "     " + i
            seclinetxt += newline
        latticeline += seclinetxt

#        if line[0] == '000':
#            prevline = data[-1]
#            prevfileline = filedata[-1]
#            data.pop()
#            filedata.pop()
#            templine = latticeline
#            latticeline = prevfileline[:-1] + templine
#            prevline = list(prevline)
#            prevline.extend(list(line[1:]))
#            line = _np.array(prevline)
        line = _np.array(numericals)
        return latticeline,line
                        
    def _remove_blanks(self,line):
        ''' Function to remove '' from lines.
            '''
        linelist=[]
        for element in line:
            if element != '' and element != '"':
                linelist.append(element)
        line = _np.array(linelist)
        return line
    

    def _get_preamble(self):        #Redundant until pybdsim can handle comments.
        '''Function to read any preamble at the start of the TRANSPORT file.
            '''
        indc,linenum = self._get_indic()
        gmadpreamble=[]
        for line in self.data[:linenum-1]:
            if line == '\r\n':
                pass
            else:
                gmadline = '!' + line
                gmadpreamble.append(gmadline)
        return gmadpreamble


    def _get_label(self,line):
        '''Function to get element label from code line.'''
        #if isinstance(line,_np.str):
        for elenum,ele in enumerate(line):
            startslash   = _string.find(ele,"/")
            startquote   = _string.find(ele,"'")
            startequal   = _string.find(ele,"=")
            startdbquote = _string.find(ele,'"')
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
            elif startdbquote != -1:
                end = 1 + startdbquote + _string.find(ele[(startdbquote+1):],'"')
                label = ele[startdbquote+1:end]
                break
            else:
                label = None
        return label
    
#        elif isinstance(line,_np.ndarray):
#            label=''
#            for element in line:
#                if element[0] == '"':
#                    label = element[1:-1]
#                    break
#                if element[0] == "/":
#                    label = element[1:-1]
#                    break
#                if element[0] == "=":
#                    label = element[1:-1]
#                    break
#                if element[0] == "'":
#                    label = element[1:-1]
#                    break
#            return label


    def _dir_exists(self,dir):
        dirs = _glob.glob('*/')
        if dir[-1] != '/':
            dir += '/'
        if dirs.__contains__(dir):
            return True
        return False

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
            
    def _get_elementdata(self,line):
        data = []
        name=''
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    data.append(_np.float(ele))
                except ValueError:
                    name = ele
        data.append(name)
        return data



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

    def _printout(self,line):
        sys.stdout.write(line+'\n')
        logfile = self._file + '_conversion.log'
        if self._outlog:
            self._logfile = open(logfile, 'a')
            self._logfile.write(line)
            self._logfile.write('\n')
            self._logfile.close()


    def _bunch_length_convert(self,bunch_length):
        ### Function to convert bunch length unit in TRANSPORT into seconds.
        scale = self.scale[self.units['bunch_length'][0]]   
        blmeters = bunch_length * scale     #Bunch length scaled to metres
        blseconds = blmeters / (self.beamprops.beta*_con.c)       #Length converted to seconds
        return blseconds

    def unit_change(self,line):
        '''Function to change the units (scaling) of various parameters'''
        #slabel,elenum = self._get_label(line)
        label = line['label']
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

        if _np.float(line['number']) == 1:    #Horizontal and vertical beam size
            self.units['x'] = label
            self.units['y'] = label
            self.units['bend_vert_gap'] = label
            #self.units['pipe_rad'] = label
        if _np.float(line['number']) == 2:    #Horizontal and vertical divergence
            self.units['xp'] = label
            self.units['yp'] = label
        if _np.float(line['number']) == 3:    #Bending Magnet Gap
            self.units['bend_vert_gap'] = label
        if _np.float(line['number']) == 4:    #Vertical Divergence ONLY
            self.units['yp'] = label
        if _np.float(line['number']) == 5:    #Pulsed Beam Length
            self.units['bunch_length'] = label
        if _np.float(line['number']) == 6:    #Momentum Spread
            self.units['momentum_spread'] = label   ## Percent
        if _np.float(line['number']) == 7:    #Bend/pole face rotation
            pass
        if _np.float(line['number']) == 8:    #Element Length
            self.units['element_length'] = label
        if _np.float(line['number']) == 9:    #Magnetic Field
            self.units['magnetic_fields'] = label
        if _np.float(line['number']) == 10:   #Mass
            self._printout('Cannot change mass scale.')
        if _np.float(line['number']) == 11:   #Momentum / energy gain during acc.
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
        energy = _np.sqrt((p_mass**2) + (momentum**2))
        self.beamprops.tot_energy = energy
        self.beamprops.tot_energy_current = energy
        self.beamprops.k_energy = energy - p_mass
        self.beamprops.gamma = energy / p_mass
        self.beamprops.beta = _np.sqrt((1.0 - (1.0 / self.beamprops.gamma**2)))
        self.beamprops.brho = mom_in_ev / _con.c

    def _calculate_momentum(self,k_energy):
        '''Function to calculate:
                Total Energy
                Kinetic Energy
                Momentum
                Lorentz factor (gamma)
                Velocity (beta)
                Magnetic rigidity (brho)
            '''
        
        k_energy = _np.float(k_energy)

        self.beamprops.k_energy = k_energy
        p_mass = self.beamprops.mass ## Particle rest mass (in GeV)
        
        e_unit = self.units['p_egain']
        if e_unit != 'eV':
            scaling = 1e9 / self.scale[e_unit[0]]     ## Scaling relative to mom. unit
        elif e_unit == 'eV':
            scaling = 1e9                               ## Scaling relative to 1 eV
        p_mass *= scaling                               ## Scale particle rest mass
        
        #energy = _np.sqrt((p_mass**2 * _con.c**2) + (momentum**2 * _con.c**2)) / _con.c

        self.beamprops.tot_energy_current = k_energy + p_mass
        self.beamprops.momentum = _np.sqrt((self.beamprops.tot_energy_current**2) - (p_mass**2))

        self.beamprops.gamma = self.beamprops.tot_energy_current / p_mass
        self.beamprops.beta = _np.sqrt((1.0 - (1.0 / self.beamprops.gamma**2)))

        if e_unit != 'eV':
            mom_in_ev = self.beamprops.momentum * self.scale[e_unit[0]]
        elif e_unit == 'eV':
            mom_in_ev = self.beamprops.momentum

        self.beamprops.brho = mom_in_ev / _con.c
