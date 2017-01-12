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
            elecode = 0
            try:
                elecode = _np.float(line[0])
            except ValueError:
                anglein = 0
                break
            
            if _np.float(line[0]) == 4.0:
                break
            elif _np.float(line[0]) == 2.0:
                endof = self._endofline(line[1])
                if endof != -1:
                    try:
                        anglein = _np.round(_np.float(line[1][:endof]),4)
                    except ValueError:
                        try:
                            anglein = _np.round(_np.float(line[2][:endof]),4)
                        except ValueError:
                            pass
                    else:
                        pass
                else:
                    try:
                        anglein = _np.round(_np.float(line[1]),4)
                    except ValueError:
                        try:
                            anglein = _np.round(_np.float(line[2]),4)
                        except ValueError:
                            pass
                    else:
                        pass
                break
            else:
                pass

        for line in self.data[linenum+1:(linenum+6)]:
            elecode = 0
            try:
                elecode = _np.float(line[0])
            except ValueError:
                angleout=0
                break
            
            if _np.float(line[0]) == 4.0:
                break
            elif _np.float(line[0]) == 2.0:
                endof = self._endofline(line[1])
                if endof != -1:
                    try:
                        angleout = _np.round(_np.float(line[1][:endof]),4)
                    except ValueError:
                        try:
                            angleout = _np.round(_np.float(line[2][:endof]),4)
                        except ValueError:
                            pass
                    else:
                        pass
                else:
                    try:
                        angleout = _np.round(_np.float(line[1]),4)
                    except ValueError:
                        try:
                            angleout = _np.round(_np.float(line[2]),4)
                        except ValueError:
                            pass
                    else:
                        pass
                break
            else:
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
        temp = _reader.reader()

        if not isinstance(input, _np.str):
            raise TypeError("Input must be a string")
        
        infile = input.split('/')[-1]       #Remove filepath, leave just filename
        self._file = infile[:-4]            #Remove extension
        self._filename = input
        isOutput = self._is_Output(input)   #Is a TRANSPORT standard output file.

        if isOutput:
            lattice,output=temp._getLatticeAndOptics(input)
            fits,fitres = temp._getFits(input)
            self._outputfits_to_registry(fitres)
            if self._debug:
                self._printout('\tAdding any fitting output to the fitting registry (self._fitReg)')
            for linenum, latticeline in enumerate(lattice):
                latticeline = latticeline.replace(';','')
                line = _np.array(latticeline.split(' '),dtype=_np.str)
                line = self._remove_illegals(line)
                
                # Method of dealing with split lines in the output
                # Should only be applicable to type 12 entry (up to 15 variables)
                # It is assumed that the line is always split, so be careful.
                prevline = lattice[linenum-1].replace(';','')
                prevline = _np.array(prevline.split(' '),dtype=_np.str)
                prevline = self._remove_illegals(prevline)
                
                try:
                    if (linenum > 0) and _np.abs(_np.float(line[0])) == 12.0:
                        latticeline, line = self._joinsplitlines(linenum,lattice)
                    # Ignore line after type 12 entry (second part of split line)
                    if (linenum > 1) and _np.abs(_np.float(prevline[0])) == 12.0:
                        pass
                    else:
                        self.data.append(line)
                        self.filedata.append(latticeline)
                except ValueError:
                    self.data.append(line)
                    self.filedata.append(latticeline)
                except IndexError:
                    pass
                    
        else:
            f = open(input)
            for inputline in f:
                endoflinepos = self._endofline(inputline)
                templine = inputline
                if endoflinepos > 0:
                    templine = inputline[:endoflinepos]
                line = _np.array(templine.split(' '),dtype=_np.str)
                # do not change comment lines
                if not line[0][0] == '(':
                    line = self._remove_illegals(line)
                self.data.append(line)
                self.filedata.append(inputline)
            f.close()

        self._fileloaded = True


    def _joinsplitlines(self,linenum,lattice):
        firstline = lattice[linenum].replace(';','')
        latticeline = firstline #Copy for later
        firstline = _np.array(firstline.split(' '),dtype=_np.str)
        firstline = self._remove_illegals(firstline)
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
        secline = self._remove_illegals(secline)
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
    

    def _remove_illegals(self,line):
        ''' Function to remove '' and stray characters from lines.
            '''
        illegal = ['"','','(',')']
        
        linelist=[]
        for element in line:
            if not illegal.__contains__(element):
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
                if end <= startslash:
                    label = ele[startslash+1:]
                else:
                    label = ele[startslash+1:end]
                break
            elif startquote != -1:
                end = 1 + startquote + _string.find(ele[(startslash+1):],"'")
                if end <= startquote:
                    label = ele[startquote+1:]
                else:
                    label = ele[startquote+1:end]
                break
            elif startequal != -1:
                end = 1 + startequal + _string.find(ele[(startslash+1):],"=")
                if end <= startequal:
                    label = ele[startequal+1:]
                else:
                    label = ele[startequal+1:end]
                break
            elif startdbquote != -1:
                end = 1 + startdbquote + _string.find(ele[(startdbquote+1):],'"')
                if end <= startdbquote:
                    label = ele[startdbquote+1:]
                else:
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
        #data.append(name)
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


    def _scale_to_meters(self,quantity):
        ''' Function to scale quantity (string) to meters, returns conversion factor'''
        if self.units[quantity] != 'm':
            conversionFactor = self.scale[self.units[quantity][0]]
        else:
            conversionFactor = 1
        return conversionFactor

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


    def _process_fits(self,fits):
        # First split the fitting output into its respective sections (input problem step).
        fitsections = []
        fitsstarts=[]
        # Start line of each section
        for linenum,line in enumerate(fits):
            if line.startswith('1'):
                fitsstarts.append(linenum)
            
        for secnum in range(len(fitsstarts)):
            if secnum+1 < len(fitsstarts):
                section = fits[fitsstarts[secnum]:fitsstarts[secnum+1]]
            else:
                section = fits[fitsstarts[secnum]:]
            lines=[]
            for line in section:
                lines.append(self._remove_illegals(line.split(' ')))
            fitsections.append(lines)
            
        magnetlines = []
        for section in fitsections:
            for line in section:
                if (len(line) > 0) and (line[0][0] == '*' and line[0][-1] == '*') and line[0] != '*FIT*':
                    magnetlines.append(line)


    def _outputfits_to_registry(self,outputdata):
        isLegal ={'*DRIFT*' : 3.0,
                  '*QUAD*'  : 5.0,
                  '*BEND*'  : 4.0}
                  
        for line in outputdata:
            append = False
            linedict = {'elementnum' : 0.0,
                        'name'       : '',
                        'length'     : 0.0}
            data    = self._remove_illegals(line.split(' '))
            eledata = self._get_elementdata(data)
            label   = self._get_label(data)
            if isLegal.__contains__(data[0]):
                linedict['elementnum']  = isLegal[data[0]]
                linedict['name']        = label
                linedict['data']        = eledata[1:] #first value is elementnum.
                linedict['length']      = eledata[1]
                append = True

            #Only add an element with a name to the fitting registry.
            #(Element has to be named to be varied in the fitting routine).
            #Otherwise update the total length of the machine.
            if append and (label is not None) and (label != ''):
                self._fitReg.AddToRegistry(linedict,line)
            else:
                self._fitReg.UpdateLength(linedict)


    def _is_Output(self,inputfile):
        ''' Function to check if a file is a standard TRANSPORT output file.
            Based upon existence of the lines:
                "0  XXX"
                
            being present, which represents the TRANSPORT indicator card line.
            X can be 0, 1, 2. Default is 0. 
            '''
        temp = _reader.reader()
        isOutput = False
        try:
            f = open(inputfile)
            for inputline in f:
                inputline = inputline.replace("\r",'')
                inputline = inputline.replace("\n",'')
                if (temp._allowedIndicatorLines.__contains__(inputline)):
                    isOutput = True
                    break
            f.close()
        except IOError:
            self._printout('Cannot open file.')
        return isOutput


    def _checkSingleLineOutputApplied(self,file):
        ''' Function to check if the control element that print element output in
            a single line was successfully applied. Check needed as not all versions
            of TRANSPORT can run this type code.
            '''
        reader = _reader.reader()
        flist  = reader._general._LoadFile(file)
        optics = reader.optics._getOptics(flist)
        for element in optics:
            if element == 'IO: UNDEFINED TYPE CODE 13. 19. ;':
                return False
        return True


    def _getTypeNum(self,line):
        ''' Function to extract the element type number (type code).
            Written because element types can contain alphabetical 
            characters when fits are used, e.g: 5.0A. Only the number 
            is required, the use of fitting does not need to be known.
            '''
        eleNum = line[0]
        if (len(eleNum) > 2):
            for characNum in range(len(eleNum[2:])):
                try:
                    converted = _np.float(eleNum[:characNum+2])
                except ValueError:
                    break
            typeNum = _np.float(eleNum[:characNum+2])
        else:
            typeNum = _np.float(eleNum)
        return typeNum

    def _transformUpdate(self,linedict):
        if self._debug and (linedict['elementnum'] == 6.0):
            errorline  = '\tElement is either a transform update or a collimator. The type code 6 definition'
            errorline2 = '\thas not been switched to collimators, therefore nothing will be done for this element.'
            self._printout(errorline)
            self._printout(errorline2)
