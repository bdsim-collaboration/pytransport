import numpy as _np
from scipy import constants as _con
import string as _string
from pybdsim import Options as _Options
from pybdsim import Builder as _pyBuilder
from pymadx import Builder as _mdBuilder
from elements import elements
import os as _os

class _beamprops():
    '''A class containing the properties of the inital beam distribution.
        '''
    def __init__(self,p_mass=938.272):
        # beam properties that are updated along the lattice
        self.momentum = 0
        self.k_energy = 0
        self.tot_energy_current = p_mass
        self.gamma = 1
        self.beta = 0
        self.brho = 0
        # beam properties that are from the initial beam and fixed
        self.mass = p_mass
        self.tot_energy = p_mass # initial energy
        self.SigmaX = 0
        self.SigmaY = 0
        self.SigmaXP = 0
        self.SigmaYP = 0
        self.SigmaE = 0
        self.SigmaT = 0
        self.X0 = 0
        self.Y0 = 0
        self.Z0 = 0
        self.T0 = 0
        self.Xp0 = 0
        self.Yp0 = 0
        self.betx = 0
        self.alfx = 0
        self.bety = 0
        self.alfy = 0
        self.emitx = 0
        self.emity = 0
        self.distrType = 'gauss'


class _machineprops():
    '''A class containing the number of elements and angular properties (i.e bending direction)
        '''
    def __init__(self):
        self.benddef = True # True = dipole defined by 4. L B n. False = dipole defined by 4. L angle n.
        self.bending = 1    # +VE = bends to the right for positive particles
        self.angle = 0      # dipole rotation angle
        self.drifts     = 0 # nr of drifts
        self.dipoles    = 0
        self.rf         = 0
        self.quads      = 0
        self.sextus     = 0
        self.transforms = 0
        self.solenoids  = 0
        self.beampiperadius = 20


class _Registry():
    def __init__(self):
        self.elements = []
        self.names    = []
        self.lines    = []

    def AddToRegistry(self,linedict,line):
        if not isinstance(linedict,dict):
            raise TypeError("Added element is not a Dictionary")
        self.elements.append(linedict)
        self.names.append(linedict['name'])
        self.lines.append(line)


class pytransport(elements):
    '''A module for converting a TRANSPORT file into gmad for use in BDSIM.
        
        To use:
            self = pytransport.convert.pytransport(inputfile)
            
        Will output the lattice in the appropriate format.

        Parameters
        -------------------------------
        
        particle: string
            The particle type, default = 'proton'.
        
        debug: boolean
            Output debug strings, default = False.
            
        distrType: string
            The distribution type of the beam, default = 'gauss'.
            Can only handle 'gauss' and 'gausstwiss'. If madx output is specified,
            the madx beam distribution is 'madx'.
            
        gmad: boolean
            Write the converted output into gmad format, default = True.
            
        madx: boolean
            write the converted outout into madx format, dafault = False.
            
        auto: boolean
            Automatically convert and output the file, default = True.

        '''
    def __init__(self,inputfile,
                 particle   = 'proton',
                 debug      = False,
                 distrType  = 'gauss',
                 gmad       = True,
                 gmadDir    = 'gmad',
                 madx       = False,
                 madxDir    = 'madx',
                 auto       = True,
                 outlog     = True):

        if particle == 'proton':
            p_mass = (_con.proton_mass) * (_con.c**2 / _con.e) / 1e9        ## Particle masses in same unit as TRANSPORT (GeV)
        elif particle == 'e-' or particle == 'e+':                          
            p_mass = (_con.electron_mass) * (_con.c**2 / _con.e) / 1e9
        self._elementReg  = _Registry()
        self._fitReg      = _Registry()
        self._particle    = particle
        self._beamdefined = False
        self._correctedbeamdef = False
        self._fileloaded  = False
        self._gmadoutput  = gmad
        self._gmadDir     = gmadDir
        self._madxoutput  = madx
        self._madxDir     = madxDir
        self._numberparts = -1
        self._collindex   = []  # An index of collimator labels
        self._accstart    = []  # An index of the start of acceleration elements.
        self.data         = []  # A list that will contain arrays of the element data
        self.filedata     = []  # A list that will contain the raw strings from the input file

        self.units        = {   # Default TRANSPORT units
        'x'                     :'cm',
        'xp'                    :'mrad',
        'y'                     :'cm',
        'yp'                    :'mrad',
        'bunch_length'          :'cm',
        'momentum_spread'       :'pc',
        'element_length'        :'m',
        'magnetic_fields'       :'kG',
        'p_egain'               :'GeV', # Momentum / energy gain during acceleration.
        'bend_vert_gap'         :'cm',  # Vertical Gap in dipoles
        'pipe_rad'              :'cm',
        'beta_func'             :'m',
        'emittance'             :'mm mrad'
            }
        self.scale={
        'p':1e-12,
        'n':1e-9,
        'u':1e-6,
        'm':1e-3,
        'c':1e-2,
        'k':1e+3,
        'K':1e+3, # Included both cases of k just in case.  
        'M':1e+6,
        'G':1e+9,
        'T':1e+12
            }
        self._debug  = debug
        #boolean for outputting stream to a log file.
        self._outlog = outlog

        #pytransport conversion classes
        self.beamprops = _beamprops(p_mass)
        self.beamprops.distrType = distrType
        self.machineprops = _machineprops()

        #a machine for both gmad and madx. Both created by default, input booleans only decide writing.
        self.gmadmachine = _pyBuilder.Machine()
        self.madxmachine = _mdBuilder.Machine()

        #load file automatically
        self._load_file(inputfile)
        if auto:
            self.transport2gmad()


    def write(self):
        if self._numberparts < 0:
            self._filename = self._file
        else:
            self._numberparts += 1
            self._filename = self._file+'_part'+_np.str((self._numberparts))
        self.create_beam()
        self.create_options()
        self.gmadmachine.AddSampler('all')
        self.madxmachine.AddSampler('all')
        if self._gmadoutput:
            if not self._dir_exists(self._gmadDir):
                _os.mkdir(self._gmadDir)
            _os.chdir(self._gmadDir)
            filename = self._filename + '.gmad'
            self.gmadmachine.Write(filename)
            _os.chdir('../')
        if self._madxoutput:
            if not self._dir_exists(self._madxDir):
                _os.mkdir(self._madxDir)
            _os.chdir(self._madxDir)
            filename = self._filename + '.madx'
            self.madxmachine.Write(filename)
            _os.chdir('../')


    def transport2gmad(self):
        '''Function to convert TRANSPORT file on a line by line basis.
            '''
        if not self._fileloaded:
            self._printout('No file loaded.')
            return
        for linenum,line in enumerate(self.data):
            if self._debug:
                self._printout('Processing line '+_np.str(linenum)+' :')
                self._printout('\t' + self.filedata[linenum])

            self._line = line
            self._linenum = linenum
            if self._is_sentinel(self._line):   # Checks if the SENTINEL line is found. SENTINEL relates to TRANSPORT
                if self._debug:                 # fitting routine and is only written after the lattice definition,
                    self._printout('Sentinel Found.')    # so there's no point reading lines beyond it.
                break
            ### Test for positive element, negative ones ignored in TRANSPORT so ignored here too.
            try: 
                if _np.float(self._line[0]) > 0:
                    if self.data[0][0] == 'OUTPUT':
                        self._element_prepper(self._line,linenum,'output')
                    else:
                        self._line = self._remove_illegals(self._line)
                        self._element_prepper(self._line,linenum,'input')
                else:
                    if self._debug:
                        self._printout('\tType code is 0 or negative, ignoring line.')
            except ValueError:
                if self._debug:
                    if line[0] == '(' or line[0] == '/':
                        errorline = '\tCannot process line '+_np.str(linenum)+', line is a comment.'
                    elif line[0] == 'S':
                        errorline = '\tCannot process line '+_np.str(linenum)+', line is for TRANSPORT fitting routine.'
                    elif line[0] == '\n':
                        errorline = '\tCannot process line '+_np.str(linenum)+', line is blank.'
                    else:
                        errorline = '\tCannot process line '+_np.str(linenum)+', reason unknown.'

                    self._printout(errorline)

        if self._debug:
            self._printout('Converting registry elements to pybdsim compatable format and adding to machine builder.')

        for element in self._elementReg.elements:
            self._get_type(element)

        self.write()
        

    def _element_prepper(self,line,linenum,filetype='input'):
        ''' Function to extract the data and prepare it for processing by each element function.
            This has been written as the lattice lines from an input file and output file are different,
            so it just a way of correctly ordering the information.
            '''
        linedict = {'elementnum' : 0.0,
                    'name'       : ''}
        
        if _np.float(line[0]) == 15.0:
            linedict['elementnum'] = 15.0
            label  = self._get_label(line)
            if filetype == 'output':
                linedict['label'] = line[2].strip('"')
            if filetype == 'input':
                linedict['label'] = label
            linedict['number'] = line[1]
            if self._debug:
                self._printout("\tEntry is a Unit Control, adding to the element registry.")
        
        if _np.float(line[0]) == 20.0:
            linedict['elementnum'] = 20.0
            angle = 0 ##Default
#                for index,ele in enumerate(line[1:]): #For loops are iterating over blank space (delimiter)
#                    if ele != ' ':
#                        endofline = self._endofline(ele)
#                        if len(ele) > 0 and endofline == -1: #I.E endofline is not in this element
#                            angle = ele
#                            break
#                        elif len(ele) > 0 and endofline != -1: #endofline is in this element
#                            angle = _np.str(ele[:endofline])
#                            break
            endofline = self._endofline(line[1])
            angle = line[1][:endofline]
            linedict['angle'] = angle
            if self._debug:
                self._printout("\tEntry is a coordinate rotation, adding to the element registry.")

    
        if _np.float(line[0]) == 1.0:
            linedict['elementnum'] = 1.0
            linedict['name'] = self._get_label(line)
            linedict['isAddition'] = False
            if self._is_addition(line,filetype):
                linedict['isAddition'] = True
            #line = self._remove_label(line)
            if len(line) < 8:
                raise IndexError("Incorrect number of beam parameters.")
    
            if filetype == 'input':
                n = 0
            elif filetype == 'output':
                n = 1
            
            #Find momentum
            #endofline = self._endofline(line[7+n])
            #if endofline != -1:
            #    linedict['momentum'] = line[7+n][:endofline]
            #else:
            linedict['momentum'] = line[7+n]
            linedict['Sigmax']  = line[1+n]
            linedict['Sigmay']  = line[3+n]
            linedict['Sigmaxp'] = line[2+n]
            linedict['Sigmayp'] = line[4+n]
            linedict['SigmaT']  = line[5+n]
            linedict['SigmaE']  = line[6+n]
            if self._debug:
                self._printout("\tEntry is a Beam definition or r.m.s addition, adding to the element registry.")
                
        
        if _np.float(line[0]) == 3.0:
            linedict['elementnum'] = 3.0
            linedict['name'] = self._get_label(line)
            data = self._get_elementdata(line)
            linedict['driftlen'] = data[0]
            if self._debug:
                self._printout("\tEntry is a drift tube, adding to the element registry.")
            
        if _np.float(line[0]) == 4.0:
            linedict['elementnum'] = 4.0
            linedict['name'] = self._get_label(line)
            linedict['linenum'] = linenum
            linedict['data'] = self._get_elementdata(line)
            e1,e2 = self._facerotation(line,linenum)
            linedict['e1'] = e1
            linedict['e2'] = e2
            if self._debug:
                self._printout("\tEntry is a dipole, adding to the element registry.")

        if _np.float(line[0]) == 5.0:
            linedict['elementnum'] = 5.0
            linedict['name'] = self._get_label(line)
            linedict['data'] = self._get_elementdata(line)
            if self._debug:
                self._printout("\tEntry is a quadrupole, adding to the element registry.")

        if _np.float(line[0]) == 6.0:
            linedict['elementnum'] = 6.0
            if self._debug:
                self._printout("\tEntry is a Transform update, adding to the element registry.")

        if _np.float(line[0]) == 12.0:
            linedict['elementnum'] = 12.0
            if filetype == 'input':
                linedict['name'] = self._get_label(line)
                linedict['data'] = self._get_elementdata(line)
            elif filetype == 'output':
                linedict['data'] = self._get_elementdata(line)[1:]
                linedict['name'] = linedict['data'][1]

            prevline = self.data[linenum-1]#.split(' ')
            linedict['prevlinenum'] = _np.float(prevline[0])
            linedict['isAddition'] = False
            if self._is_addition(line):
                linedict['isAddition'] = True
            if self._debug:
                self._printout("\tEntry is a beam rotation, adding to the element registry.")
                
        if _np.float(line[0]) == 11.0:
            linedict['elementnum'] = 11.0
            linedict['name'] = self._get_label(line)
            linedict['data'] = self._get_elementdata(line)
            if self._debug:
                self._printout("\tEntry is an acceleration element, adding to the element registry.")
        
        if _np.float(line[0]) == 13.0:
            linedict['elementnum'] = 13.0
            linedict['data'] = self._get_elementdata(line)
            if self._debug:
                self._printout("\tEntry is a Input/Output control, adding to the element registry.")
        
        if _np.float(line[0]) == 16.0:
            linedict['elementnum'] = 16.0
            linedict['data'] = self._get_elementdata(line)
            if self._debug:
                self._printout("\tEntry is a special input, adding to the element registry.")

        if _np.float(line[0]) == 18.0:
            linedict['elementnum'] = 18.0
            linedict['name'] = self._get_label(line)
            linedict['data'] = self._get_elementdata(line)
            if self._debug:
                self._printout("\tEntry is a sextupole, adding to the element registry.")
        
        if _np.float(line[0]) == 19.0:
            linedict['elementnum'] = 19.0
            linedict['name'] = self._get_label(line)
            linedict['data'] = self._get_elementdata(line)
            if self._debug:
                self._printout("\tEntry is a solenoid, adding to the element registry.")
        
        if _np.float(line[0]) == 9.0:
            linedict['elementnum'] = 9.0
            if self._debug:
                self._printout("\tEntry is a repetition control, adding to the element registry.")

        rawline = self.filedata[linenum]
        self._elementReg.AddToRegistry(linedict,rawline)


    def _get_type(self,linedict):
        '''Function to read element type.
            '''
        if linedict['elementnum'] == 15.0:
            self.unit_change(linedict)
        if linedict['elementnum'] == 20.0:
            self.change_bend(linedict)
        if linedict['elementnum'] == 1.0:
            self.define_beam(linedict)
        if linedict['elementnum'] == 3.0:
            self.drift(linedict)
        if linedict['elementnum'] == 4.0:
            self.dipole(linedict)
        if linedict['elementnum'] == 5.0:
            self.quadrupole(linedict)
        if linedict['elementnum'] == 6.0:
            pass
            #self.collimator(line)
        if linedict['elementnum'] == 12.0:
            self.correction(linedict)
        if linedict['elementnum'] == 11.0:
            self.acceleration(linedict)
        if linedict['elementnum'] == 13.0:
            self.printline(linedict)
        if linedict['elementnum'] == 16.0:
            self.special_input(linedict)
        if linedict['elementnum'] == 18.0:
            self.sextupole(linedict)
        if linedict['elementnum'] == 19.0:
            self.solenoid(line)

        # 9.  : 'Repetition' - for nesting elements
        if linedict['elementnum'] == 9.0:
            errorline = '\tWARNING Repetition Element not implemented in converter!' + _np.str(linenum) + '\n'
            self._printout(errorline)

        ### OTHER TYPES WHICH CAN BE IGNORED:
        # 2.  : Dipole poleface rotation (handled in dipole line).
        # 6.0.X : Update RX matrix used in TRANSPORT
        # 7.  : 'Shift beam centroid'
        # 8.  : Magnet alignment tolerances
        # 10. : Fitting constraint
        # 14. : Arbitrary transformation of TRANSPORT matrix
       
    def create_beam(self):
        '''Function to prepare the beam for writing.
            '''
        #different beam objects depending on output type
        self.madxbeam = self.madxmachine.beam
        self.gmadbeam = self.gmadmachine.beam
        
        #convert energy to GeV (madx only handles GeV)
        energy_in_gev = self.beamprops.tot_energy * self.scale[self.units['p_egain'][0]] / 1e9
        self.beamprops.tot_energy = energy_in_gev
        
        self.madxbeam.SetParticleType(self._particle)
        self.madxbeam.SetEnergy(energy=self.beamprops.tot_energy,unitsstring = 'GeV')

        self.gmadbeam.SetParticleType(self._particle)
        self.gmadbeam.SetEnergy(energy=self.beamprops.tot_energy,unitsstring = 'GeV')

        #set gmad parameters depending on distribution
        if self.beamprops.distrType == 'gausstwiss':
            self.gmadbeam.SetDistributionType(self.beamprops.distrType)
            self.gmadbeam.SetBetaX(self.beamprops.betx)
            self.gmadbeam.SetBetaY(self.beamprops.bety)
            self.gmadbeam.SetAlphaX(self.beamprops.alfx)
            self.gmadbeam.SetAlphaY(self.beamprops.alfy)
            self.gmadbeam.SetEmittanceX(self.beamprops.emitx,unitsstring='mm')
            self.gmadbeam.SetEmittanceY(self.beamprops.emity,unitsstring='mm')
            self.gmadbeam.SetSigmaE(self.beamprops.SigmaE)
            self.gmadbeam.SetSigmaT(self.beamprops.SigmaT)

        else:
            self.gmadbeam.SetDistributionType(self.beamprops.distrType)
            self.gmadbeam.SetSigmaX(self.beamprops.SigmaX,unitsstring=self.units['x'])
            self.gmadbeam.SetSigmaY(self.beamprops.SigmaY,unitsstring=self.units['y'])
            self.gmadbeam.SetSigmaXP(self.beamprops.SigmaXP,unitsstring=self.units['xp'])
            self.gmadbeam.SetSigmaYP(self.beamprops.SigmaYP,unitsstring=self.units['yp'])
            self.gmadbeam.SetSigmaE(self.beamprops.SigmaE)
            self.gmadbeam.SetSigmaT(self.beamprops.SigmaT)
            
            # calculate betas and emittances regardless for madx beam
            try:
                self.beamprops.betx = self.beamprops.SigmaX / self.beamprops.SigmaXP
            except ZeroDivisionError:
                self.beamprops.betx= 0
            try:
                self.beamprops.bety = self.beamprops.SigmaY / self.beamprops.SigmaYP
            except ZeroDivisionError:
                self.beamprops.bety= 0
            self.beamprops.emitx = self.beamprops.SigmaX * self.beamprops.SigmaXP / 1000.0
            self.beamprops.emity = self.beamprops.SigmaY * self.beamprops.SigmaYP / 1000.0
        
        #set madx beam
        self.madxbeam.SetDistributionType('madx')
        self.madxbeam.SetBetaX(self.beamprops.betx)
        self.madxbeam.SetBetaY(self.beamprops.bety)
        self.madxbeam.SetAlphaX(self.beamprops.alfx)
        self.madxbeam.SetAlphaY(self.beamprops.alfy)
        self.madxbeam.SetEmittanceX(self.beamprops.emitx/1000)
        self.madxbeam.SetEmittanceY(self.beamprops.emity/1000)
        self.madxbeam.SetSigmaE(self.beamprops.SigmaE)
        self.madxbeam.SetSigmaT(self.beamprops.SigmaT)
        
        #set beam offsets in gmad if non zero
        if self.beamprops.X0 != 0:
            self.gmadbeam.SetX0(self.beamprops.X0,unitsstring=self.units['x'])
        if self.beamprops.Y0 != 0:
            self.gmadbeam.SetY0(self.beamprops.Y0,unitsstring=self.units['y'])
        if self.beamprops.Z0 != 0:
            self.gmadbeam.SetZ0(self.beamprops.Z0,unitsstring=self.units['x'])

        
        if self._debug:
            self._printout('\t Beam definition :')
            self._printout('\t distrType = ' + self.beamprops.distrType)
            self._printout('\t energy = ' + _np.str(self.beamprops.tot_energy)+ ' GeV')
            self._printout('\t SigmaX = ' + _np.str(self.beamprops.SigmaX)  + ' ' +self.units['x'])
            self._printout('\t SigmaXP = '+ _np.str(self.beamprops.SigmaXP) + ' ' +self.units['xp'])
            self._printout('\t SigmaY = ' + _np.str(self.beamprops.SigmaY)  + ' ' +self.units['y'])
            self._printout('\t SigmaYP = '+ _np.str(self.beamprops.SigmaYP) + ' ' +self.units['yp'])
            self._printout('\t SigmaE = ' + _np.str(self.beamprops.SigmaE))
            self._printout('\t SigmaT = ' + _np.str(self.beamprops.SigmaT))
            self._printout('\t (Final brho = '+_np.str(_np.round(self.beamprops.brho,2))+' Tm)')
            self._printout('\t Twiss Params:')
            self._printout('\t BetaX = ' +_np.str(self.beamprops.betx) + ' ' + self.units['beta_func'])
            self._printout('\t BetaY = ' +_np.str(self.beamprops.bety) + ' ' + self.units['beta_func'])
            self._printout('\t AlphaX = '+_np.str(self.beamprops.alfx))
            self._printout('\t AlphaY = '+_np.str(self.beamprops.alfy))
            self._printout('\t Emittx = '+_np.str(self.beamprops.emitx) + ' ' + self.units['emittance'])
            self._printout('\t EmittY = '+_np.str(self.beamprops.emity) + ' ' + self.units['emittance'])

        self.gmadmachine.AddBeam(self.gmadbeam)
        self.madxmachine.AddBeam(self.madxbeam)


    def create_options(self):
        '''Function to set the Options for the BDSIM machine.'''
        self.options = _Options.Options()
        self.options.SetPhysicsList(physicslist='em')
        self.options.SetBeamPipeRadius(beampiperadius=self.machineprops.beampiperadius,unitsstring=self.units['pipe_rad'])
        self.options.SetOuterDiameter(outerdiameter=0.5,unitsstring='m')
        self.options.SetTunnelRadius(tunnelradius=1,unitsstring='m')
        self.options.SetBeamPipeThickness(bpt=5,unitsstring='mm')
        self.options.SetSamplerDiameter(radius=1,unitsstring='m')
        self.options.SetStopTracks(stop=True)
        
        self.gmadmachine.AddOptions(self.options)
        self.madxmachine.AddOptions(self.options)   #redundant

    



