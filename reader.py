import numpy as _np
import pybdsim

class reader():
    def get_data(self,file,type=None):
        if isinstance(type,_np.str):
            if type == 'beam':
                transdata = self.get_beam_output(file)
                return transdata
            elif type == 'standard':
                transdata = self.get_standard_output(file)
                return transdata
    
        f = open(file)
        flist=[]
        transdata = None
        for line in f:
            #remove any carriage returns (both Mac and Unix)
            line = line.rstrip('\r\n')
            flist.append(line)
            splitline = self._remove_blanks(line.split(' '))
            if splitline and splitline[0] == '*BEAM*':     ## Is beam output
                #print "'*BEAM*' found in line " + _np.str(i+1)
                f.close()
                transdata = self.get_beam_output(file)
                break
            elif line == '0    0':
                #print "'0    0' found in line " + _np.str(i+1)
                f.close()
                transdata = self.get_standard_output(file)
                break
            else:
                pass
        if transdata == None:
            errorstring =  "Could not find an indicator in the input file for either\n"
            errorstring += "a beam output file (indicator = '*BEAM*), or a standard output file (indicator = '0    0').\n"
            errorstring += "Please check the input file or specify the input type with the type argument in the get_output function.\n"
            errorstring += "Note that the only accepted values for type are 'standard' or 'beam'."
            raise IOError(errorstring)
        return transdata

    def get_beam_output(self,file):
        ''' Returns a BDSAsciiData instance of parameters from the input file.
            The input file is assumed to contain the beam data as output 
            manually by TRANSPORT. As such, the structure of said output 
            is assumed to remain unchanged from lattice to lattice. This code 
            is designed to work with that output only. An example of that 
            output would be:
            
            *DRIFT*      z = 34.372 m
            *SIGMA*
             Center:         0.000 mm    0.000 mrad    0.000 mm    0.000 mrad
             horz. Par. :   18.361 mm   23.134 mrad    0.997
             vert. Par. :    8.447 mm   10.242 mrad   -0.995
            *TWISS PARAMETERS* (for dp/p = 1.000 % )
               alfax:     betax:         alfay:     betay:
             -13.51800   10.75807 m      9.55361    7.92270 m
            *TRANSFORM 1*
                horz:                              vert:
              -14.22000   -0.64752   -9.32539      0.13646    0.97529    0.00000
              -18.40864   -0.90858  -10.33185     -1.18959   -1.17399    0.00000

            Anything other than this format will be read incorrectly or will 
            not be read. There has to be a blank line between every element's 
            output in the file.
            
            Note:
                It is assumed that the element (3,2) in the displayed *TRANSFORM 1*
                matrices corresponds to the dispersion (-9.32539 and 0.0000 in the above
                example). Some output however appears to have the correct magnitude, but 
                incorrect sign. This doesn't affect the resulting beam size, but beware 
                that a direct dispersion comparison to another lattice may appear incorrect.
            '''
        flist = self._file_to_list(file)
        transdata = self.process_beam_output(flist)
        return transdata
  
  
    def get_standard_output(self,file):
        flist = self._file_to_list(file)
        lattice,output = self._get_latticeandoutput(flist)
        transdata = self._process_elements(output)
        return transdata
    

    def _file_to_list(self,file):
        ''' Converts the input file into a list. The data has to be in a format other than
            handling line by line (using next()). 
            
            The function for processing the file reads section by section rather than line
            by line. This may be an inefficent method but the input file should not be very 
            large so it should not require a large amount of memory.
            '''
        self.file = file
        flist=[]
        infile = open(file)
        ##Loop over lines and remove any carriage returns (both Mac and Unix)
        for line in infile:
            cleanline = line.rstrip('\r\n')
            flist.append(cleanline)
        infile.close()
        return flist


    def _get_lattice(self,flist):
        ''' Function to extract the lattice from a standard output file. No processing at the moment, but
            this function should identify the chunk that is the lattice.
            '''
        foundlatticestart = False
        foundlatticeend = False

        for linenum,line in enumerate(flist):
            if line == '0    0':
                if not foundlatticestart:
                    latticestart = linenum+1
                    foundlatticestart = True
            if line == '0SENTINEL':
                if not foundlatticeend:
                    latticeend = linenum
                    foundlatticeend = True
        if not foundlatticestart:
            if not foundlatticeend:
                raise IOError('No lattice found in '+self.file+'.')
            else:
                errorstring  = 'The end of a lattice (line = "0SENTINEL") was found at line '+ _np.str(latticeend+1)+',\n'
                errorstring += 'but the start of a lattice (line = "0    0") was not found. Please check the input file.'
                raise IOError(errorstring)
        elif not foundlatticeend:
                errorstring  = 'The start of a lattice (line = "0    0") was found at line '+ _np.str(latticestart-1)+',\n'
                errorstring += 'but the end of a lattice (line = "0SENTINEL") was not found. Please check the input file.'
                raise IOError(errorstring)
        else:
            lattice = flist[latticestart:latticeend]
        return lattice


    def _get_output(self,flist):
        ''' Function to extract the output from a standard output file. The output will be an list of the lines
            for each element which conatins the beam data. Each element should contain the R and TRANSPORT matrices which
            are necessary so the beam info can be calculated.
            '''
        foundoutputstart = False
        foundoutputend = False

        for linenum,line in enumerate(flist):
            splitline = self._remove_blanks(line.split(' '))
            if splitline and splitline[0] == '*BEAM*' and not foundoutputstart:
                outputstart=linenum
                foundoutputstart = True
            if splitline and splitline[0] == '0*LENGTH*':
                outputend = linenum
                foundoutputend = True
                break
        if not foundoutputstart:
            if not foundoutputend:
                raise IOError('No output found in '+self.file+'.')
            else:
                errorstring  = 'The end of a lattice (line containing "0*LENGTH*") was found at line '+ _np.str(outputend+1)+',\n'
                errorstring += 'but the start of a lattice (first line containing "*BEAM*") was not found. Please check the input file.'
                raise IOError(errorstring)
        elif not foundoutputend:
                errorstring  = 'The start of a lattice (first line containing "*BEAM*") was found at line '+ _np.str(outputstart-1)+',\n'
                errorstring += 'but the end of a lattice (line containing "0*LENGTH*") was not found. Please check the input file.'
                raise IOError(errorstring)
        else:
            output = flist[outputstart:outputend]

        #Split the list of all element data into their individual elements.
        elementlist=[]
        elementstart = False
        elementend = False
        finalelement = False
        for linenum,line in enumerate(output):
            if linenum == (len(output)-1):
                finalelement = True
            try:
                if line[1] == '*':
                    if line[2:11] == 'TRANSFORM':   #Is midway through element output
                        pass
                    elif elementstart == False: #Current line must be start of the element
                        elementstart = True
                        elementstartline = linenum
                    elif elementstart == True:
                        elementend = True       #Otherwise the line must be the start of the next element
                        elementendline = linenum
                if elementstart and finalelement:
                    element = output[elementstartline:]
                    elementlist.append(element)
                if elementstart and elementend: #If the start and end of the element are found, append and reset
                    element = output[elementstartline:elementendline]
                    elementlist.append(element)
                    elementend = False
                    elementstart = False
                if elementstart == False and elementend == False:   #Though if it's been reset, it must be because the
                    elementstart = True                             #current line is the start of next element, so set the
                    elementstartline = linenum                      #start line for the next element
            except IndexError:
                pass
        return elementlist


    def _get_latticeandoutput(self,flist):
        lattice = self._get_lattice(flist)
        output = self._get_output(flist)
        return lattice,output
    

    def _process_elements(self,elementlist):
        transdata = {
            'Sigma_x'   :[],
            'Sigma_xp'  :[],
            'Sigma_y'   :[],
            'Sigma_yp'  :[],
            'S'         :[],
            'Alph_x'    :[],
            'Alph_y'    :[],
            'Beta_x'    :[],
            'Beta_y'    :[],
            'Emitt_x'   :[],
            'Emitt_y'   :[],
            'Disp_x'    :[],
            'Disp_y'    :[],
            'Sigma_p'   :[],
            'Momentum'  :[],
            'E'         :[], # kinetic energy
            'Name'      :[],
            'Type'      :[],
            }
        num_elements = 0
        # initialise momentum/energy since not given for every element
        momentum = 0.0
        energy   = 0.0
        proton_mass = 938.272
        for element in elementlist:
            if len(element) > 1:  # I.e not a fit or matrix-modifying element
                # type is in between * can have a space (for space charge *SP CH*)
                type    = element[0].split('*')[1]
                # rest of first line split with spaces
                splitline = self._remove_blanks(element[0].split('*')[2].split(' '))
                name    = splitline[1].strip('"')
                if type=="BEAM" or type=="ACC":
                    momentum = _np.float(splitline[-2])
                    energy = _np.sqrt(proton_mass*proton_mass + momentum*momentum) - proton_mass
                s       = _np.float(self._remove_blanks(element[1].split(' '))[0])
                sigx    = _np.float(self._remove_blanks(element[1].split(' '))[3])
                sigxp   = _np.float(self._remove_blanks(element[2].split(' '))[1])
                sigy    = _np.float(self._remove_blanks(element[3].split(' '))[1])
                sigyp   = _np.float(self._remove_blanks(element[4].split(' '))[1])
                sigt    = _np.float(self._remove_blanks(element[5].split(' '))[1])
                sigp    = _np.float(self._remove_blanks(element[6].split(' '))[1])

                r21     = _np.float(self._remove_blanks(element[2].split(' '))[3])
                r43     = _np.float(self._remove_blanks(element[4].split(' '))[5])
                dx      = _np.float(self._split_negatives(self._remove_blanks(element[8].split(' ')))[5])
                dy      = _np.float(self._split_negatives(self._remove_blanks(element[10].split(' ')))[5])

                xpint   = _np.sqrt(sigxp**2 * (1 - r21**2))
                ypint   = _np.sqrt(sigyp**2 * (1 - r43**2))
                
                ex      = sigx*xpint
                ey      = sigy*ypint
                betx    = sigx**2 / ex
                bety    = sigy**2 / ey
                gammax  = sigxp**2 / ex
                gammay  = sigyp**2 / ey
                alfx    = _np.sqrt((gammax * betx)-1)
                alfy    = _np.sqrt((gammax * betx)-1)

                transdata['Sigma_x'].append(sigx/1000)
                transdata['Sigma_xp'].append(sigxp)
                transdata['Sigma_y'].append(sigy/1000)
                transdata['Sigma_yp'].append(sigyp)
                transdata['S'].append(s)
                transdata['Alph_x'].append(alfx)
                transdata['Alph_y'].append(alfy)
                transdata['Beta_x'].append(betx)
                transdata['Beta_y'].append(bety)
                transdata['Emitt_x'].append(ex)
                transdata['Emitt_y'].append(ey)
                transdata['Disp_x'].append(dx)
                transdata['Disp_y'].append(dy)
                transdata['Sigma_p'].append(sigp)
                transdata['Momentum'].append(momentum)
                transdata['E'].append(energy)
                transdata['Name'].append(name)
                transdata['Type'].append(type)
                num_elements += 1 

        def get_elementdata(index):             # Function to get the data for each element, rather than each key.
            elementlist2=[]                      # There's probably a better container type for this, but I'm familiar with
            for name in transdata.keys():   # dictionaries, and it works (for now).
                elementlist2.append(transdata[name][index])
        
            return elementlist2

        data = pybdsim.Data.BDSAsciiData()      # Now convert the dict into BDSAsciiData instance for final output.
        for name in transdata.keys():
            data._AddProperty(name)
        for i in range(num_elements):
            data.append(get_elementdata(i))

        return data



    def process_beam_output(self,flist):
        transdata = {
            'Sigma_x'   :[],
            'Sigma_xp'  :[],
            'Sigma_y'   :[],
            'Sigma_yp'  :[],
            'S'         :[],
            'Alph_x'    :[],
            'Alph_y'    :[],
            'Beta_x'    :[],
            'Beta_y'    :[],
            'Emitt_x'   :[],
            'Emitt_y'   :[],
            'Disp_x'    :[],
            'Disp_y'    :[],
            'Sigma_p'   :[],
            'Name'      :[],
            }
        num_elements = 0
        for elenum,element in enumerate(flist):
            if element == '':   ### The first line of the section should be a blank line.
                try:
                    section = flist[elenum+1:elenum+13]
                    
                    ### Get the element name, type and start position (in s)
                    line0 = section[0]
                    line0split = line0.split('*')
                    typestr = line0split[1]
                    reststr = line0split[2]
                    restsplit = reststr.split(' ')
                    restsplit = self._remove_blanks(restsplit)
                    transdata['S'].append(_np.float(restsplit[2]))
                    namestr = restsplit[4]
                    transdata['Name'].append(namestr)
                    #transdata['type'].append(typestr)

                    ### Get sigma_x and sigma_xp
                    line3 = section[3].split(' ')
                    line3 = self._remove_blanks(line3)
                    transdata['Sigma_x'].append(_np.float(line3[3])/1000)
                    transdata['Sigma_xp'].append(_np.float(line3[5])/1000)

                    ### Get sigma_y and sigma_yp
                    line4 = section[4].split(' ')
                    line4 = self._remove_blanks(line4)
                    transdata['Sigma_y'].append(_np.float(line4[3])/1000)
                    transdata['Sigma_yp'].append(_np.float(line4[5])/1000)
                    
                    ### Get momentum spread
                    line5 = section[5].split(' ')
                    transdata['Sigma_p'].append(_np.float(line5[5])/100)
                    
                    ### Get alfa and beta twiss transdata for x and y.
                    line7 = section[7].split(' ')
                    line7 = self._remove_blanks(line7)
                    transdata['Alph_x'].append(_np.float(line7[0]))
                    transdata['Beta_x'].append(_np.float(line7[1]))
                    transdata['Alph_y'].append(_np.float(line7[3]))
                    transdata['Beta_y'].append(_np.float(line7[4]))

                    ### Get horizontal and vertical dispersion
                    line10 = section[10].split(' ')
                    line10 = self._remove_blanks(line10)
                    transdata['Disp_x'].append(_np.float(line10[2])/10)
                    transdata['Disp_y'].append(_np.float(line10[5])/10)
                    
                    ### Terms for calculating the emittance.
                    term1x = _np.float(line3[3])**2
                    term2x = (_np.float(line10[2])*(_np.float(line5[5])/100))**2
                    term1y = _np.float(line4[3])**2
                    term2y = (_np.float(line10[5])*(_np.float(line5[5])/100))**2

                    ### Get horizontal and vertical emittance
                    emittx = (term1x - term2x)/ _np.float(line7[1])
                    emitty = (term1y - term2y)/ _np.float(line7[4])
                    transdata['Emitt_x'].append(emittx)
                    transdata['Emitt_y'].append(emitty)
                
                    num_elements += 1
                except ValueError:
                    errstr = "Could not process section beginning at line " + _np.str(elenum) + " : "
                    print(errstr)
                    print " "
                    print(element)

        def get_elementdata(index):             # Function to get the data for each element, rather than each key.
            elementlist=[]                      # There's probably a better container type for this, but I'm familiar with
            for name in transdata.keys():   # dictionaries, and it works (for now).
                elementlist.append(transdata[name][index])
            return elementlist

        data = pybdsim.Data.BDSAsciiData()      # Now convert the dict into BDSAsciiData instance for final output.
        for name in transdata.keys():
            data._AddProperty(name)
        for i in range(num_elements):
            data.append(get_elementdata(i))

        return data


    def _remove_blanks(self,line):
        ''' Removes any blanks from a string (blanks being '' and not white spaces).
            '''
        newline=''
        for element in line:
            if element != '':
                newline += element
                newline += ' '
        stripline = newline.split(' ')
        # remove last element as it will be blank (due to added space)
        return stripline[:-1]

    def _split_negatives(self,line):
        newline=[]
        for element in line:
            negpos = element.find('-')
            if negpos > 1:
                parts = element.split('-')
                newline.append(parts[0])
                newline.append('-'+parts[1])
            else:
                newline.append(element)
        return newline




