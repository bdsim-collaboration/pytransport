import numpy as _np
from scipy import constants as _con
import string as _string
from pybdsim import Options as _Options
from pybdsim import Beam as _Beam
from pybdsim import Builder as _Builder
from elements import elements

class _beamprops():      #Beam properties
    def __init__(self,p_mass=938.272):
        self.momentum = 0
        self.mass = p_mass
        self.k_energy = 0
        self.tot_energy = p_mass
        self.gamma = 1
        self.beta = 0
        self.brho = 0


class _elementprops():       #Number of elements and angular properties
    def __init__(self):
        self.bending = 1   #+VE = bends to the right for positive particles
        self.angle = 0      #dipole rotation angle
        self.drifts = 1
        self.dipoles = 1
        self.quads = 1
        self.sext = 1
        self.transforms = 1


class pytransport(elements):
    '''A module for converting a TRANSPORT file into gmad for use in BDSIM.
        
        To use:
            self = pytransport.convert.pytransport()
            self.load_file(TRANSPORTfile)
            self.convert()
            
        Will output:
        
        '''
    def __init__(self,particle='proton'):
        if particle == 'proton':
            p_mass = (_con.proton_mass) * (_con.c**2 / _con.e) / 1e9        ## Particle masses in same unit as TRANSPORT (GeV)
        elif particle == 'e-' or particle == 'e+':                          
            p_mass = (_con.electron_mass) * (_con.c**2 / _con.e) / 1e9
        self._particle = particle
        self._beamfilemade = False
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
        'pipe_rad'              :'cm'
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
        'T':1e+12}
        
        self.beamprops = _beamprops(p_mass)
        self.elementprops = _elementprops()
        self.machine = _Builder.Machine()

        
    def load_file(self,infile):
        '''Load file to be converted into gmad format.
            '''
        self.data=[]
        self._file = infile[:-4]
        self._filename = self._file+'.gmad'
        try:
            for line in open(infile):
                self.data.append(line)
        except IOError:
            print 'Cannot open file.'
                

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




    def convert(self):
        '''Function to convert TRANSPORT file on a line by line basis.
            '''
        for linenum,line in enumerate(self.data):
            if len(line) > 1:   #i.e line isn't equal to escape sequence line.
                                #This is a bit slapdash at the moment, needs better implementation.
                a = _np.array(line.split(' '))
                if self._is_sentinel(a):
                    print('Sentinel Found.')
                    break
                linelist=[]
                for element in a:
                    if element != '':
                        linelist.append(element)
                a = _np.array(linelist)
                ### Test for positive element, negative ones ignored in TRANSPORT so ignored here too.
                try: 
                    if _np.float(a[0]) > 0:
                        self._get_type(a)
                except ValueError:
                    if line[0] == '(' or line[0] == '/':
                        errorline = 'Cannot process line '+_np.str(linenum)+', line is a comment.'
                    elif line[0] == 'S':
                        errorline = 'Cannot process line '+_np.str(linenum)+', line is for TRANSPORT fitting routine'
                    else:
                        errorline = 'Cannot process line '+_np.str(linenum)+': \n'
                        errorline += '\t'+line[:-2]

                    print(errorline)

        self.create_options()
        self.machine.AddSampler('all')
        self.machine.Write(self._filename)
                


    def _get_type(self,line):
        '''Function to read element type.
            '''
        if _np.float(line[0]) == 15.0:      
            self.unit_change(line)
        if _np.float(line[0]) == 20.0:    
            self.change_bend(line)
        if _np.float(line[0]) == 1.0 and not self._is_addition(line):  
            if not self._beamfilemade:      ### To ignore beam r.m.s additions
                self.create_beam(line)
        if _np.float(line[0]) == 3.0:     
            self.drift(line)
        if _np.float(line[0]) == 4.0:     
            self.dipole(line)
        if _np.float(line[0]) == 5.0:     
            self.quadrupole(line)
        if _np.float(line[0]) == 6.0:     
            pass
            #self.collimator(line)
        if _np.float(line[0]) == 11.0:    
            self.acceleration(line)
                
        ### OTHER TYPES WHICH CAN BE IGNORED:
        # 2.  : Dipole fringe fields.
        # 6.0.X : Update RX matrix used in TRANSPORT
        # 7.  : 'Shift beam centroid'
        # 8.  : Magnet alignment tolerances
        # 9.  : 'Repetition' - for nesting elements
        # 10. : Fitting constraint
        # 12. : Something to do with the outputting the beam for use in another TRANSPORT system
        # 13. : Print output to terminal
        # 14. : Arbitrary transformation of TRANSPORT matrix
        
 

    def create_beam(self,line):
        '''Defines the beam. Will ALWAYS be a Gaussian. Must input the TRANSPORT line that defines the beam.
            '''
        line = self._remove_label(line)
        if len(line) < 8:
            raise IndexError("Incorrect number of beam parameters.")
        self.beam = _Beam.Beam()
        self.beam.SetParticleType(self._particle)
        endofline = self._endofline(line[7])
        
        #Find momentum
        if endofline != -1:
            momentum = line[7][:endofline]
        else:
            momentum = line[7]
        
        #Convert momentum to energy and set beam type/distribution.
        self._calculate_energy(momentum)  
        self.beam.SetEnergy(energy=_np.round(self.beamprops.tot_energy,3),unitsstring=self.units['p_egain'])
        self.beam.SetDistributionType('gauss')
            
        self.beam.SetSigmaX(sigmax=_np.float(line[1]),unitsstring=self.units['x'])
        self.beam.SetSigmaY(sigmay=_np.float(line[3]),unitsstring=self.units['y'])
        self.beam.SetSigmaXP(sigmaxp=_np.float(line[2]),unitsstring=self.units['xp'])
        self.beam.SetSigmaYP(sigmayp=_np.float(line[4]),unitsstring=self.units['yp'])
        self.beam.SetSigmaE(sigmae=_np.float(line[6]))
        
        bunchl = self._bunch_length_convert(_np.float(line[5])) ## Get bunch length in seconds.
        self.beam.SetSigmaT(sigmat=bunchl)     
        self._beamfilemade = True
        self.machine.AddBeam(self.beam)



    def create_options(self):
        '''Function to create the Options gmad file.'''
        self.options = _Options.Options()
        self.options.SetPhysicsList(physicslist='em_standard')
        self.options.SetBeamPipeRadius(beampiperadius=10,unitsstring='cm')
        self.options.SetOuterDiameter(outerdiameter=0.5,unitsstring='m')
        self.options.SetTunnelRadius(tunnelradius=1,unitsstring='m')
        self.options.SetBeamPipeThickness(bpt=5,unitsstring='mm')
        self.options.SetSamplerDiameter(radius=2,unitsstring='m')
        self.machine.AddOptions(self.options)

    



