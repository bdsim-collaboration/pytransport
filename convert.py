import numpy as _np
from scipy import constants as _con
import string as _string
from pybdsim import Options as _Options
from pybdsim import Builder as _pyBuilder
from pymadx import Builder as _mdBuilder
from elements import elements

class _beamprops():
    '''A class containing the properties of the inital beam distribution.
        '''
    def __init__(self,p_mass=938.272):
        self.momentum = 0
        self.mass = p_mass
        self.k_energy = 0
        self.tot_energy = p_mass
        self.gamma = 1
        self.beta = 0
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
        self.brho = 0
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
        self.benddef = True   #True = dipole defined by 4. L B n. False = dipole defined by 4. L angle n.
        self.bending = 1   #+VE = bends to the right for positive particles
        self.angle = 0      #dipole rotation angle
        self.drifts = 1
        self.dipoles = 1
        self.quads = 1
        self.sextus = 1
        self.transforms = 1
        self.solenoids = 1
        self.beampiperadius = 20


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
                 madx       = False,
                 auto       = True):

        if particle == 'proton':
            p_mass = (_con.proton_mass) * (_con.c**2 / _con.e) / 1e9        ## Particle masses in same unit as TRANSPORT (GeV)
        elif particle == 'e-' or particle == 'e+':                          
            p_mass = (_con.electron_mass) * (_con.c**2 / _con.e) / 1e9
        self._particle = particle
        self._beamdefined = False
        self._correctedbeamdef = False
        self._fileloaded = False
        self._gmadoutput = gmad
        self._madxoutput = madx
        self._numberparts = -1
        self._collindex=[]  # An index of collimator labels
        self._accstart=[]   # An index of the start of acceleration elements.
        self.units={    ### Default TRANSPORT units
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
        self._debug = False
        if debug:
            self._debug = True

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
            self._filename = self._file + '.gmad'
        else:
            self._numberparts += 1
            self._filename = self._file+'_part'+_np.str((self._numberparts))
        self.create_beam()
        self.create_options()
        self.gmadmachine.AddSampler('all')
        self.madxmachine.AddSampler('all')
        if self._gmadoutput:
            self.gmadmachine.Write(self._file)
        if self._madxoutput:
            self.madxmachine.Write(self._file)

    def transport2gmad(self):
        '''Function to convert TRANSPORT file on a line by line basis.
            '''
        if not self._fileloaded:
            print('No file loaded.')
            return
        for linenum,line in enumerate(self.data):
            if self._debug:
                print('Processing line '+_np.str(linenum)+' :')
                print('\t' + line)
            if len(line) > 1:   #i.e line isn't equal to escape sequence line.
                                #This is a bit slapdash at the moment, needs better implementation.
                self._line = _np.array(line.split(' '))
                if self._is_sentinel(self._line):
                    if self._debug:
                        print('Sentinel Found.')
                    break
                linelist=[]
                for element in self._line:
                    if element != '':
                        linelist.append(element)
                self._line = _np.array(linelist)
                ### Test for positive element, negative ones ignored in TRANSPORT so ignored here too.
                try: 
                    if _np.float(self._line[0]) > 0:
                        self._get_type(self._line,linenum)
                    #if self._debug:
                    #    print('\n')
                except ValueError:
                    if self._debug:
                        if line[0] == '(' or line[0] == '/':
                            errorline = '\tCannot process line '+_np.str(linenum)+', line is a comment.\n'
                        elif line[0] == 'S':
                            errorline = '\tCannot process line '+_np.str(linenum)+', line is for TRANSPORT fitting routine\n'
                        else:
                            errorline = '\tCannot process line '+_np.str(linenum)+' \n'

                        print(errorline)
        self.write()


    def _get_type(self,line,linenum):
        '''Function to read element type.
            '''
        if _np.float(line[0]) == 15.0:      
            self.unit_change(line)
        if _np.float(line[0]) == 20.0:    
            self.change_bend(line)
        if _np.float(line[0]) == 1.0:
            self.define_beam(line)
        if _np.float(line[0]) == 3.0:     
            self.drift(line)
        if _np.float(line[0]) == 4.0:     
            self.dipole(line,linenum)
        if _np.float(line[0]) == 5.0:     
            self.quadrupole(line)
        if _np.float(line[0]) == 6.0:     
            pass
            #self.collimator(line)
        if _np.float(line[0]) == 12.0:
            self.correction(line,linenum)
        if _np.float(line[0]) == 11.0:    
            self.acceleration(line)
        if _np.float(line[0]) == 13.0:
            self.printline(line)
        if _np.float(line[0]) == 16.0:
            self.special_input(line)
        if _np.float(line[0]) == 18.0:
            self.sextupole(line)
        if _np.float(line[0]) == 19.0:
            self.solenoid(line)

        ### OTHER TYPES WHICH CAN BE IGNORED:
        # 2.  : Dipole fringe fields.
        # 6.0.X : Update RX matrix used in TRANSPORT
        # 7.  : 'Shift beam centroid'
        # 8.  : Magnet alignment tolerances
        # 9.  : 'Repetition' - for nesting elements
        # 10. : Fitting constraint
        # 14. : Arbitrary transformation of TRANSPORT matrix
        
 

    def create_beam(self):
        '''Defines the beam. Will ALWAYS be a Gaussian. Must input the TRANSPORT line that defines the beam.
            '''
        self.beam = _Beam.Beam()
        self.beam.SetParticleType(self._particle)
        self.beam.SetEnergy(energy=_np.round(self.beamprops.tot_energy,3),unitsstring=self.units['p_egain'])

        if self.beamprops.distrType == 'gausstwiss':
            self.beam.SetDistributionType(self.beamprops.distrType)
            self.beam.SetBetaX(self.beamprops.betx)
            self.beam.SetBetaY(self.beamprops.bety)
            self.beam.SetAlphaX(self.beamprops.alfx)
            self.beam.SetAlphaY(self.beamprops.alfy)
            self.beam.SetEmittanceX(self.beamprops.emitx,unitsstring='mm')
            self.beam.SetEmittanceY(self.beamprops.emity,unitsstring='mm')
            self.beam.SetSigmaE(self.beamprops.SigmaE)
            self.beam.SetSigmaT(self.beamprops.SigmaT)
        
        else:
            self.beam.SetDistributionType(self.beamprops.distrType)
            self.beam.SetSigmaX(self.beamprops.SigmaX,unitsstring=self.units['x'])
            self.beam.SetSigmaY(self.beamprops.SigmaY,unitsstring=self.units['y'])
            self.beam.SetSigmaXP(self.beamprops.SigmaXP,unitsstring=self.units['xp'])
            self.beam.SetSigmaYP(self.beamprops.SigmaYP,unitsstring=self.units['yp'])
            self.beam.SetSigmaE(self.beamprops.SigmaE)
            self.beam.SetSigmaT(self.beamprops.SigmaT)
        
        self.beam.SetX0(self.beamprops.X0,unitsstring=self.units['x'])
        self.beam.SetY0(self.beamprops.Y0,unitsstring=self.units['y'])
        self.beam.SetZ0(self.beamprops.Z0,unitsstring=self.units['x'])
        
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
            print('\t (brho = '+_np.str(_np.round(self.beamprops.brho,2))+' Tm)')
            print('\t Twiss Params:')
            print('\t BetaX = ' +_np.str(self.beamprops.betx) + ' ' + self.units['beta_func'])
            print('\t BetaY = ' +_np.str(self.beamprops.bety) + ' ' + self.units['beta_func'])
            print('\t AlphaX = '+_np.str(self.beamprops.alfx))
            print('\t AlphaY = '+_np.str(self.beamprops.alfy))
            print('\t Emittx = '+_np.str(self.beamprops.emitx) + ' ' + self.units['emittance'])
            print('\t EmittY = '+_np.str(self.beamprops.emity) + ' ' + self.units['emittance'])

        self.machine.AddBeam(self.beam)



    def create_options(self):
        '''Function to create the Options gmad file.'''
        self.options = _Options.Options()
        self.options.SetPhysicsList(physicslist='em')
        self.options.SetBeamPipeRadius(beampiperadius=self.machineprops.beampiperadius,unitsstring=self.units['pipe_rad'])
        self.options.SetOuterDiameter(outerdiameter=0.5,unitsstring='m')
        self.options.SetTunnelRadius(tunnelradius=1,unitsstring='m')
        self.options.SetBeamPipeThickness(bpt=5,unitsstring='mm')
        self.options.SetSamplerDiameter(radius=1,unitsstring='m')
        self.machine.AddOptions(self.options)

    



