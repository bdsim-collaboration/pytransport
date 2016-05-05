import numpy as _np

class reader():
    def get_data(self,file):
        ''' Returns a dictionary of parameters from the input file.
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
        params = self._process_data(flist)
        return params

    def _file_to_list(self,file):
        ''' Converts the input file into a list. The data has to be in a format other than
            handling line by line (using next()). 
            
            The function for processing the file reads section by section rather than line
            by line. This may be an inefficent method but the input file should not be very 
            large so it should not require a large amount of memory.
            '''
        num_lines = _np.sum(1 for line in open(file))
        flist=[]
        infile = open(file)
        for i in range(num_lines):
            line = infile.next()
            if line[-1] == '\n':
                line = line[:-1]
            if line[-1] == '\r':
                line = line[:-1]
            flist.append(line)
        infile.close()
        return flist


    def _process_data(self,flist):
        params = {
            'sigx'      :[],
            'sigxp'     :[],
            'sigy'      :[],
            'sigyp'     :[],
            's'         :[],
            'alfx'      :[],
            'alfy'      :[],
            'betx'      :[],
            'bety'      :[],
            'ex'        :[],
            'ey'        :[],
            'dx'        :[],
            'dy'        :[],
            'sigp'      :[],
            'type'      :[],
            'name'      :[],
            'code'      :'transport'
            }

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
                    params['s'].append(_np.float(restsplit[2]))
                    namestr = restsplit[4]
                    params['name'].append(namestr)
                    params['type'].append(typestr)

                    ### Get sigma_x and sigma_xp
                    line3 = section[3].split(' ')
                    line3 = self._remove_blanks(line3)
                    params['sigx'].append(_np.float(line3[3])/1000)
                    params['sigxp'].append(_np.float(line3[5])/1000)

                    ### Get sigma_y and sigma_yp
                    line4 = section[4].split(' ')
                    line4 = self._remove_blanks(line4)
                    params['sigy'].append(_np.float(line4[3])/1000)
                    params['sigyp'].append(_np.float(line4[5])/1000)
                    
                    ### Get momentum spread
                    line5 = section[5].split(' ')
                    params['sigp'].append(_np.float(line5[5])/100)
                    
                    ### Get alfa and beta twiss params for x and y.
                    line7 = section[7].split(' ')
                    line7 = self._remove_blanks(line7)
                    params['alfx'].append(_np.float(line7[0]))
                    params['betx'].append(_np.float(line7[1]))
                    params['alfy'].append(_np.float(line7[3]))
                    params['bety'].append(_np.float(line7[4]))

                    ### Get horizontal and vertical dispersion
                    line10 = section[10].split(' ')
                    line10 = self._remove_blanks(line10)
                    params['dx'].append(_np.float(line10[2])/10)
                    params['dy'].append(_np.float(line10[5])/10)
                    
                    ### Terms for calculating the emittance.
                    term1x = _np.float(line3[3])**2
                    term2x = (_np.float(line10[2])*(_np.float(line5[5])/100))**2
                    term1y = _np.float(line4[3])**2
                    term2y = (_np.float(line10[5])*(_np.float(line5[5])/100))**2

                    ### Get horizontal and vertical emittance
                    emittx = (term1x - term2x)/ _np.float(line7[1])
                    emitty = (term1y - term2y)/ _np.float(line7[4])
                    params['ex'].append(emittx)
                    params['ey'].append(emitty)
                
                except ValueError:
                    errstr = "Could not process section beginning at line " + _np.str(elenum) + " : "
                    print(errstr)
                    print " "
                    print(element)
            
        for key in params.keys():
            params[key] = _np.array(params[key])

        return params


    def _remove_blanks(self,line):
        ''' Removes any blanks from a string (blanks being '' and not white spaces).
            '''
        newline=''
        for element in line:
            if element != '':
                newline += element
                newline += ' '
        stripline = newline.split(' ')
        return stripline







