import numpy as _np
from scipy import constants as _con
from pymadx import Builder as _mdBuilder
from pybdsim import Builder as _pyBuilder
import string as _string
from _General import functions
    
class elements(functions):
    def define_beam(self,line):
        if self._is_addition(line):
            if self._debug:
                print('\tIgnoring beam rms addition.')
            return
        if self._beamdefined:
            self._numberparts += 1
            self.write()
            print('Writing.')
            del self.gmadmachine
            del self.madxmachine
            self.gmadmachine = _pyBuilder.Machine()
            self.madxmachine = _mdBuilder.Machine()
            self._correctedbeamdef = False
            
            print('\tBeam redefinition found. Writing previous section to file.')
            print('\tSplitting into multiple machines.')
        
        line = self._remove_label(line)
        if len(line) < 8:
            raise IndexError("Incorrect number of beam parameters.")

        endofline = self._endofline(line[7])
        
        #Find momentum
        if endofline != -1:
            momentum = line[7][:endofline]
        else:
            momentum = line[7]

        self._beamdefined = True
       
        #Convert momentum to energy and set distribution params.
        self._calculate_energy(momentum)
        self.beamprops.SigmaX  = _np.float(line[1])
        self.beamprops.SigmaY  = _np.float(line[3])
        self.beamprops.SigmaXP = _np.float(line[2])
        self.beamprops.SigmaYP = _np.float(line[4])
        self.beamprops.SigmaE  = _np.float(line[6]) * 0.01 * (self.beamprops.beta**2) ## Convert from percentage mom spread to absolute espread
        self.beamprops.SigmaT  = self._bunch_length_convert(_np.float(line[5])) ## Get bunch length in seconds.

        
        #Calculate Initial Twiss params
        self.beamprops.betx = self.beamprops.SigmaX / self.beamprops.SigmaXP
        self.beamprops.bety = self.beamprops.SigmaY / self.beamprops.SigmaYP
        self.beamprops.emitx = self.beamprops.SigmaX * self.beamprops.SigmaXP / 1000.0
        self.beamprops.emity = self.beamprops.SigmaY * self.beamprops.SigmaYP / 1000.0

        if self._debug:
            print('\t Beam definition :')
            print('\t distrType = ' + self.beamprops.distrType)
            print('\t energy = ' + _np.str(_np.round(self.beamprops.tot_energy,3))+ ' ' +self.units['p_egain'])
            print('\t SigmaX = ' + _np.str(self.beamprops.SigmaX)  + ' ' +self.units['x'])
            print('\t SigmaXP = '+ _np.str(self.beamprops.SigmaXP) + ' ' +self.units['xp'])
            print('\t SigmaY = ' + _np.str(self.beamprops.SigmaY)  + ' ' +self.units['y'])
            print('\t SigmaYP = '+ _np.str(self.beamprops.SigmaYP) + ' ' +self.units['yp'])
            print('\t SigmaE = ' + _np.str(self.beamprops.SigmaE))
            print('\t SigmaT = ' + _np.str(self.beamprops.SigmaT))
            print('\t (brho = '  + _np.str(_np.round(self.beamprops.brho,2))+' Tm)')
            print('\t Twiss Params:')
            print('\t BetaX = ' +_np.str(self.beamprops.betx) + ' ' + self.units['beta_func'])
            print('\t BetaY = ' +_np.str(self.beamprops.bety) + ' ' + self.units['beta_func'])
            print('\t AlphaX = '+_np.str(self.beamprops.alfx))
            print('\t AlphaY = '+_np.str(self.beamprops.alfy))
            print('\t Emittx = '+_np.str(self.beamprops.emitx) + ' ' + self.units['emittance'])
            print('\t EmittY = '+_np.str(self.beamprops.emity) + ' ' + self.units['emittance'])
        
    
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

        elementid = 'DR'+_np.str(self.machineprops.drifts)
        self.machineprops.drifts += 1

        self.gmadmachine.AddDrift(name=elementid,length=length_in_metres)
        self.madxmachine.AddDrift(name=elementid,length=length_in_metres)
        
        if self._debug:
            print('\tConverted to:')
            debugstring = 'Drift '+elementid+', length '+_np.str(length_in_metres)+' m'
            print('\t'+debugstring)

        
    def dipole(self,line,linenum):
        label = self._get_label(line)
        dipoledata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    dipoledata.append(_np.float(ele))
                except ValueError:
                    dipoledata.append(ele)
        length = dipoledata[0]          # First two non-blanks must be the entries in a specific order.
        
        ## Get poleface rotationa
        e1 = self._facerotation(line,linenum-1) * (_np.pi / 180.0)   ## Entrance pole face rotation.
        e2 = self._facerotation(line,linenum+1) * (_np.pi / 180.0)   ## Exit pole face rotation.
        if self._debug:
            if e1 != 0:
                print('\tPreceding element ('+_np.str(linenum-1)+') provides an entrance poleface rotation of '+_np.str(_np.round(e1,4))+' rad.')
            if e2 != 0:
                print('\tFollowing element ('+_np.str(linenum+1)+') provides an exit poleface rotation of '+_np.str(_np.round(e2,4))+' rad.')
        
        ##Calculate bending angle
        if self.machineprops.benddef:
            bfield = dipoledata[1]
            field_in_Gauss = bfield * self.scale[self.units['magnetic_fields'][0]]  # Scale to Gauss
            field_in_Tesla = field_in_Gauss * 1e-4                                  # Convert to Tesla
            if field_in_Tesla == 0:
                angle = 0                                                           # zero field = zero angle
            else:
                rho = self.beamprops.brho / (_np.float(field_in_Tesla))             # Calculate bending radius.
                angle = (_np.float(length) / rho) * self.machineprops.bending       # for direction of bend
            if self._debug:
                print('\tbfield = '+_np.str(field_in_Gauss)+' kG')
                print('\tbfield = '+_np.str(field_in_Tesla)+' T')
                print('\tCorresponds to angle of '+_np.str(_np.round(angle,4)) + ' rad.')
        elif not self.machineprops.benddef:
            angle_in_deg = dipoledata[1]
            angle = angle_in_deg * (_np.pi/180.) * self.machineprops.bending
        
        ##Convert element length
        if self.units['element_length'] != 'm':
            length_in_metres = length * self.scale[self.units['element_length'][0]]
        else:
            length_in_metres = length
        
        elementid = 'BM'+_np.str(self.machineprops.dipoles)
        self.machineprops.dipoles += 1
        
        ##Check for non zero pole face rotation
        if (e1 != 0) and (e2 != 0):
            self.gmadmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4),e1=_np.round(e1,4),e2=_np.round(e2,4))
            self.madxmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4),e1=_np.round(e1,4),e2=_np.round(e2,4))
        elif (e1 != 0) and (e2 == 0):
            self.gmadmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4),e1=_np.round(e1,4))
            self.madxmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4),e1=_np.round(e1,4))
        elif (e1 == 0) and (e2 != 0):
            self.gmadmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4),e2=_np.round(e2,4))
            self.madxmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4),e2=_np.round(e2,4))
        else:
            self.gmadmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4))
            self.madxmachine.AddDipole(name=elementid,category='sbend',length=length_in_metres,angle=_np.round(angle,4))

        ## Debug output
        if self._debug:
            if (e1 != 0) and (e2 != 0):
                polefacestr = ', e1= '+_np.str(_np.round(e1,4))+' rad, e2= '+_np.str(_np.round(e2,4))+' rad'
            elif (e1 != 0) and (e2 == 0):
                polefacestr = ', e1= '+_np.str(_np.round(e1,4))+' rad'
            elif (e1 == 0) and (e2 != 0):
                polefacestr = ', e2= '+_np.str(_np.round(e2,4))+' rad'
            else:
                polefacestr = ''

            print('\tConverted to:')
            debugstring = 'Dipole '+elementid+', length= '+_np.str(length_in_metres)+' m, angle= '+_np.str(_np.round(angle,4))+' rad'+polefacestr
            print('\t'+debugstring)


    def change_bend(self,line):
        '''Function to change the direction of the dipole bend. Can be a direction other than horizontal (i.e != n*pi).
            '''
        ## NOT FULLY TESTED.
        angle = 0
        rotation = False
        for index,ele in enumerate(line[1:]): #For loops are iterating over blank space (delimiter)
            if ele != ' ':
                endofline = self._endofline(ele)
                if len(ele) > 0 and endofline == -1: #I.E endofline is not in this element
                    angle = ele
                    break
                elif len(ele) > 0 and endofline != -1: #endofline is in this element
                    angle = _np.str(ele[:endofline])
                    break
        self.machineprops.angle = _np.float(angle)
        if self.machineprops.angle >= 360:
            self.machineprops.angle = _np.mod(self.machineprops.angle,360)
        if self.machineprops.angle <= -360:
            self.machineprops.angle = _np.mod(self.machineprops.angle,-360)

        if self.machineprops.angle == 180 or self.machineprops.angle == -180: #If 180 degrees, switch bending angle
            self.machineprops.bending *= -1

        
        elif self.machineprops.angle != 0:                        #If not 180 degrees, use transform3d.      
            self.machineprops.angle *= -1                         #For conversion to correct direction. Eg in TRANSPORT -90 is upwards, in BDSIM, 90 is upwards.  
            anginrad = self.machineprops.angle * (_np.pi / 180)
            elementid = 't'+_np.str(self.machineprops.transforms)
            self.machineprops.transforms += 1
            self.gmadmachine.AddTransform3D(name=elementid,psi=anginrad)
            self.madxmachine.AddTransform3D(name=elementid,psi=anginrad)
            rotation = True
        
        if self._debug:
            if rotation:
                print('\tConverted to:')
                debugstring = 'Transform3D '+elementid+', angle '+_np.str(_np.round(angle,4))+' rad'
                print('\t'+debugstring)
            elif self.machineprops.angle == 1:
                print('Bending direction set to Right')
            elif self.machineprops.angle == -1:
                print('Bending direction set to Left')
        


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
                elementid = 'QF'+_np.str(self.machineprops.quads)
            elif field_gradient < 0:
                elementid = 'QD'+_np.str(self.machineprops.quads)
            else:
                elementid = 'NULLQUAD'+_np.str(self.machineprops.quads)  #For K1 = 0. 
        
        self.machineprops.quads += 1

        self.gmadmachine.AddQuadrupole(name=elementid,length=length_in_metres,k1=_np.round(field_gradient,4))
        self.madxmachine.AddQuadrupole(name=elementid,length=length_in_metres,k1=_np.round(field_gradient,4))
        
        if self._debug:
            string1 = '\tQuadrupole, field in gauss = ' + _np.str(field_in_Gauss) + ' KG, field in Tesla = ' + _np.str(field_in_Tesla) + ' T.'
            string2 = '\tBeampipe radius = ' + _np.str(pipe_in_metres) + ' m. Field gradient = '+ _np.str(field_in_Tesla/pipe_in_metres) + ' T/m.'
            string3 = '\tBrho = ' + _np.str(_np.round(self.beamprops.brho,4)) + ' Tm. K1 = ' +_np.str(_np.round(field_gradient,4)) + ' m^-2'
            print(string1)
            print(string2)
            print(string3)
            print('\tConverted to:')
            debugstring = 'Quadrupole '+elementid+', length= '+_np.str(length_in_metres)+' m, k1= '+_np.str(_np.round(field_gradient,4))+' T/m'
            print('\t'+debugstring)




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
        
        elementid = 'SEXT'+_np.str(self.machineprops.sextus)
        
        self.machineprops.sextus += 1
        
        self.gmadmachine.AddSextupole(name=elementid,length=length_in_metres,k2=_np.round(field_gradient,4))
        self.madxmachine.AddSextupole(name=elementid,length=length_in_metres,k2=_np.round(field_gradient,4))


        if self._debug:
            print('\tConverted to:')
            debugstring = 'Sextupole '+elementid+', length '+_np.str(length_in_metres)+' m, k2 '+_np.str(_np.round(field_gradient,4))+' T/m^2'
            print('\t'+debugstring)



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
                
        elementid = 'SOLE'+_np.str(self.machineprops.solenoids)

        self.machineprops.solenoids += 1
        
        self.gmadmachine.AddSolenoid(name=elementid,length=length_in_metres,ks=_np.round(field_in_Tesla,4))
        self.madxmachine.AddSolenoid(name=elementid,length=length_in_metres,ks=_np.round(field_in_Tesla,4))

        if self._debug:
            print('\tConverted to:')
            debugstring = 'Solenoid '+elementid+', length '+_np.str(length_in_metres)+' m, ks '+_np.str(_np.round(field_in_Tesla,4))+' T'
            print('\t'+debugstring)



    def printline(self,line):
        label = self._get_label(line)
        for ele in line[1:]:
            try:
                number = _np.float(ele)
                if number == 48:
                    self.machineprops.benddef = False
                    print('Switched Dipoles to Angle definition.')
                if number == 47:
                    self.machineprops.benddef = True
                    print('Switched Dipoles to field definition.')
            except ValueError:
                dummy=0



    def correction(self,line,linenum):
        if self._correctedbeamdef == True:
            print('\t Not Correction to original beam definition')
            return
        #Check if the previous line was the original beam definition and not an rms update
        prevline = self.data[linenum-1].split(' ')
        if _np.float(prevline[0]) == 1.0 and not self._is_addition(line) and self._beamdefined:
            self._correctedbeamdef = True
        
        label = self._get_label(line)
        correctiondata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    correctiondata.append(_np.float(ele))
                except ValueError:
                    correctiondata.append(ele)

        if len(correctiondata) > 15: #15 sigma elements
            sigma21 = correctiondata[0]
            sigma43 = correctiondata[5]
        else:
            print('\tLength of correction line is incorrect')
            return

        emittoverbeta = self.beamprops.SigmaXP**2 * (1 - sigma21**2)
        emittbeta = self.beamprops.SigmaX**2
        betx = _np.sqrt(emittbeta / emittoverbeta)
        emitx = emittbeta / betx
        slope = sigma21 * self.beamprops.SigmaXP / self.beamprops.SigmaX
        alfx = -1.0 * slope * betx
        
        self.beamprops.betx = betx
        self.beamprops.emitx = emitx / 1000.0
        self.beamprops.alfx = alfx
        
        emittoverbeta = self.beamprops.SigmaYP**2 * (1 - sigma43**2)
        emittbeta = self.beamprops.SigmaY**2
        bety = _np.sqrt(emittbeta / emittoverbeta)
        emity = emittbeta / bety
        slope = sigma43 * self.beamprops.SigmaYP / self.beamprops.SigmaY
        alfy = -1.0 * slope * bety

        self.beamprops.bety = bety
        self.beamprops.emity = emity / 1000.0
        self.beamprops.alfy = alfy

        self.beamprops.distrType = 'gausstwiss'

        if self._debug:
            print('\tConverted to:')
            print('\t Beam Correction. Sigma21 = ' + _np.str(sigma21) + ', Sigma43 = '  + _np.str(sigma43) + '.')
            print('\t Beam distribution type now switched to "gausstwiss":')
            print('\t Twiss Params:')
            print('\t BetaX = ' +_np.str(self.beamprops.betx) + ' ' + self.units['beta_func'])
            print('\t BetaY = ' +_np.str(self.beamprops.bety) + ' ' + self.units['beta_func'])
            print('\t AlphaX = '+_np.str(self.beamprops.alfx))
            print('\t AlphaY = '+_np.str(self.beamprops.alfy))
            print('\t Emittx = '+_np.str(self.beamprops.emitx) + ' ' + self.units['emittance'])
            print('\t EmittY = '+_np.str(self.beamprops.emity) + ' ' + self.units['emittance'])



    def special_input(self,line):
        label = self._get_label(line)
        specialdata = []
        for index,ele in enumerate(line[1:]):
            if ele != '':
                try:
                    specialdata.append(_np.float(ele))
                except ValueError:
                    specialdata.append(ele)
        if specialdata[0] == 16.0:  #X0 offset
            self.beamprops.X0 = specialdata[1]
        if specialdata[0] == 17.0:  #Y0 offset
            self.beamprops.Y0 = specialdata[1]
        if specialdata[0] == 18.0:  #Z0 offset
            self.beamprops.Z0 = specialdata[1]
        #if specialdata[0] == 5.0:   #beampiperadius (technically only vertical, but will apply a circle for now)
        #    self.machineprops.beampiperadius = specialdata[1]

        #if self._debug:
        #    print('\tConverted to:')
        #    print('\t'+_np.str(specialdata[2]))

