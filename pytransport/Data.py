import numpy as _np
import os as _os
from scipy import constants as _con
import _General

useRootNumpy = True

try:
    import root_numpy as _rnp
except ImportError:
    useRootNumpy = False
    pass


def Load(filepath):
    extension = filepath.split('.')[-1]
    if not _os.path.isfile(filepath):
        raise IOError("File does not exist")
    if ("elosshist" in filepath) or (".hist" in filepath):
        return _LoadAsciiHistogram(filepath)
    elif extension == 'root':
        try:
            return _LoadRoot(filepath)
        except NameError:
            # raise error rather than return None, saves later scripting errors.
            raise IOError('Root loader not available.')
    else:
        raise IOError("Unknown file type - not BDSIM data")


def _LoadAsciiHistogram(filepath):
    data = BDSData()
    f = open(filepath, 'r')
    for i, line in enumerate(f):
        # first line is header (0 counting)
        if i == 1:
            names, units = _ParseHeaderLine(line)
            for name, unit in zip(names, units):
                data._AddProperty(name, unit)
        elif "nderflow" in line:
            data.underflow = float(line.strip().split()[1])
        elif "verflow" in line:
            data.overflow = float(line.strip().split()[1])
        elif i >= 4:
            data.append(tuple(map(float, line.split())))
    f.close()
    return data


def _LoadRoot(filepath):
    if not useRootNumpy:
        raise IOError("root_numpy not available - can't load ROOT file")
    data = BDSData()
    trees = _rnp.list_trees(filepath)

    if 'optics' in trees:
        branches = _rnp.list_branches(filepath, 'optics')
        treedata = _rnp.root2array(filepath, 'optics')
    elif 'orbit' in trees:
        branches = _rnp.list_branches(filepath, 'orbit')
        treedata = _rnp.root2array(filepath, 'orbit')
    else:
        raise IOError("This file doesn't have the required tree 'optics'.")
    for element in range(len(treedata[branches[0]])):
        elementlist = []
        for branch in branches:
            if element == 0:
                data._AddProperty(branch)
            elementlist.append(treedata[branch][element])
        data.append(elementlist)
    return data


def _ParseHeaderLine(line):
    names = []
    units = []
    for word in line.split():
        if word.count('[') > 0:
            names.append(word.split('[')[0])
            units.append(word.split('[')[1].strip(']'))
        else:
            names.append(word)
            units.append('NA')
    return names, units


class BDSData(list):
    """
    General class representing simple 2 column data.

    Inherits python list.  It's a list of tuples with extra columns of 'name' and 'units'.
    """
    def __init__(self, *args, **kwargs):
        list.__init__(self, *args, **kwargs)
        self.units   = []
        self.names   = []
        self.columns = self.names

    def __getitem__(self, index):
        return dict(zip(self.names, list.__getitem__(self, index)))

    def GetItemTuple(self, index):
        """
        Get a specific entry in the data as a tuple of values rather than a dictionary.
        """
        return list.__getitem__(self, index)
        
    def _AddMethod(self, variablename):
        """
        This is used to dynamically add a getter function for a variable name.
        """
        def GetAttribute():
            if self.names.count(variablename) == 0:
                raise KeyError(variablename+" is not a variable in this data")
            ind = self.names.index(variablename)
            return _np.array([event[ind] for event in self])
        setattr(self, variablename, GetAttribute)

    def ConcatenateMachine(self, *args):
        """
        This is used to concatenate machines.
        """
        # Get final position of the machine (different param for survey)
        lastSpos = self.GetColumn('S')[-1]
        
        for machine in args:
            if isinstance(machine, _np.str):
                machine = Load(machine)
        
            # check names sets are equal
            if len(set(self.names).difference(set(machine.names))) != 0:
                raise AttributeError("Cannot concatenate machine, variable names do not match")

            if self.names.count('S') != 0:
                sind = self.names.index('S')
            else:
                raise KeyError("S is not a variable in this data")
        
            # Have to convert each element to a list as tuples can't be modified
            for index in range(len(machine)):
                element = machine.GetItemTuple(index)
                elementlist = list(element)
                elementlist[sind] += lastSpos
                
                self.append(tuple(elementlist))

            lastSpos += machine.GetColumn('S')[-1]

    def _AddProperty(self, variablename, variableunit='NA'):
        """
        This is used to add a new variable and hence new getter function
        """
        self.names.append(variablename)
        self.units.append(variableunit)
        self._AddMethod(variablename)

    def _DuplicateNamesUnits(self, bdsdata2instance):
        d = bdsdata2instance
        for name, unit in zip(d.names, d.units):
            self._AddProperty(name, unit)

    def MatchValue(self, parametername, matchvalue, tolerance):
        """
        This is used to filter the instance of the class based on matching
        a parameter withing a certain tolerance.

        >>> a = pytransport.Data.Load("myfile.txt")
        >>> a.MatchValue("S",0.3,0.0004)
        
        this will match the "S" variable in instance "a" to the value of 0.3
        within +- 0.0004.

        You can therefore used to match any parameter.

        Return type is BDSAsciiData
        """
        if hasattr(self, parametername):
            a = BDSData()                 #build bdsdata2
            a._DuplicateNamesUnits(self)  #copy names and units
            pindex = a.names.index(parametername)
            filtereddata = [event for event in self if abs(event[pindex] - matchvalue) <= tolerance]
            a.extend(filtereddata)
            return a
        else:
            print "The parameter: ", parametername, " does not exist in this instance"

    def Filter(self, booleanarray):
        """
        Filter the data with a booleanarray.  Where true, will return
        that event in the data.

        Return type is BDSData
        """
        a = BDSData()
        a._DuplicateNamesUnits(self)
        a.extend([event for i, event in enumerate(self) if booleanarray[i]])
        return a

    def NameFromNearestS(self, S):
        i = self.IndexFromNearestS(S)
        if not hasattr(self, "Name"):
            raise ValueError("This file doesn't have the required column Name")
        return self.Name()[i]
    
    def IndexFromNearestS(self, S):
        """
        IndexFromNearestS(S) 

        return the index of the beamline element clostest to S 

        Only works if "SStart" column exists in data
        """
        # check this particular instance has the required columns for this function
        if not hasattr(self, "SStart"):
            raise ValueError("This file doesn't have the required column SStart")
        if not hasattr(self, "Arc_len"):
            raise ValueError("This file doesn't have the required column Arc_len")
        s = self.SStart()
        l = self.Arc_len()

        # iterate over beamline and record element if S is between the
        # sposition of that element and then next one
        # note madx S position is the end of the element by default
        ci = [i for i in range(len(self) - 1) if (S > s[i] and S < s[i]+l[i])]
        try:
            ci = ci[0] # return just the first match - should only be one
        except IndexError:
            # protect against S positions outside range of machine
            if S > s[-1]:
                ci -= 1
            else:
                ci = 0
        # check the absolute distance to each and return the closest one
        # make robust against s positions outside machine range
        return ci

    def GetColumn(self, columnstring):
        """
        Return a numpy array of the values in columnstring in order
        as they appear in the beamline
        """
        if columnstring not in self.columns:
            raise ValueError("Invalid column name")
        ind = self.names.index(columnstring)
        return _np.array([element[ind] for element in self])

    def __repr__(self):
        s = ''
        s += 'pytransport.Data.BDSData instance\n'
        s += str(len(self)) + ' entries'
        return s


class ConversionData:
    def __init__(self,
                 inputfile,
                 machine,
                 options       = None,  # None as madx has no options class. gmad needs pybdsim.Options passing in
                 particle      = 'proton',
                 debug         = False,
                 distrType     = 'gauss',
                 gmad          = True,
                 gmadDir       = 'gmad',
                 madx          = False,
                 madxDir       = 'madx',
                 auto          = True,
                 dontSplit     = False,
                 keepName      = False,
                 combineDrifts = False,
                 outlog        = True):

        if particle == 'proton':
            p_mass = _con.proton_mass * (_con.c ** 2 / _con.e) / 1e9  # Particle masses in same unit as TRANSPORT (GeV)
        elif particle == 'e-' or particle == 'e+':
            p_mass = _con.electron_mass * (_con.c ** 2 / _con.e) / 1e9
        else:
            p_mass = 1

        # pytransport data container classes
        self.convprops = _General._conversionProps(inputfile, particle, debug, gmad, gmadDir, madx, madxDir,
                                                   auto, dontSplit, keepName, combineDrifts, outlog)
        self.beamprops = _General._beamprops(p_mass)
        self.beamprops.distrType = distrType
        self.machineprops = _General._machineprops()
        self.options = options

        # the gmad/madx machine and beam that will be written.
        self.machine = machine
        self.beam = self.machine.beam

        # make a copy of the empty machine. Copy needed in case machine is split and a new machine is needed.
        self._machineCopy = machine

        # initialise registries
        self.ElementRegistry = _General._Registry()
        self.FitRegistry = _General._Registry()

        self.units = {  # Default TRANSPORT units
            'x': 'cm',
            'xp': 'mrad',
            'y': 'cm',
            'yp': 'mrad',
            'bunch_length': 'cm',
            'momentum_spread': 'pc',
            'element_length': 'm',
            'magnetic_fields': 'kG',
            'p_egain': 'GeV',  # Momentum / energy gain during acceleration.
            'bend_vert_gap': 'cm',  # Vertical half-gap in dipoles
            'pipe_rad': 'cm',
            'beta_func': 'm',
            'emittance': 'mm mrad'
        }
        self.scale = {
            'p': 1e-12,
            'n': 1e-9,
            'u': 1e-6,
            'm': 1e-3,
            'c': 1e-2,
            'k': 1e+3,
            'K': 1e+3,  # Included both cases of k just in case.
            'M': 1e+6,
            'G': 1e+9,
            'T': 1e+12
        }

        self.accstart = []  # An index of the start of acceleration elements.
        self.data = []  # A list that will contain arrays of the element data
        self.filedata = []  # A list that will contain the raw strings from the input file

    def AddOptions(self):
        """
        Function to set the Options for the BDSIM machine.
        """
        self.options.SetPhysicsList(physicslist='em')
        self.options.SetBeamPipeRadius(beampiperadius=self.machineprops.beampiperadius,
                                       unitsstring=self.units['pipe_rad'])
        self.options.SetOuterDiameter(outerdiameter=0.5, unitsstring='m')
        self.options.SetTunnelRadius(tunnelradius=1, unitsstring='m')
        self.options.SetBeamPipeThickness(bpt=5, unitsstring='mm')
        self.options.SetSamplerDiameter(radius=1, unitsstring='m')
        self.options.SetStopTracks(stop=True)
        self.options.SetIncludeFringeFields(on=True)

        self.machine.AddOptions(self.options)

    def AddBeam(self):
        """
        Function to prepare the beam for writing.
        """
        # convert energy to GeV (madx only handles GeV)
        energy_in_gev = self.beamprops.tot_energy * self.scale[self.units['p_egain'][0]] / 1e9
        self.beamprops.tot_energy = energy_in_gev

        self.beam.SetParticleType(self.convprops.particle)
        self.beam.SetEnergy(energy=self.beamprops.tot_energy, unitsstring='GeV')

        if self.convprops.gmadoutput:
            # set gmad parameters depending on distribution
            if self.beamprops.distrType == 'gausstwiss':
                self.beam.SetDistributionType(self.beamprops.distrType)
                self.beam.SetBetaX(self.beamprops.betx)
                self.beam.SetBetaY(self.beamprops.bety)
                self.beam.SetAlphaX(self.beamprops.alfx)
                self.beam.SetAlphaY(self.beamprops.alfy)
                self.beam.SetEmittanceX(self.beamprops.emitx, unitsstring='mm')
                self.beam.SetEmittanceY(self.beamprops.emity, unitsstring='mm')
                self.beam.SetSigmaE(self.beamprops.SigmaE)
                self.beam.SetSigmaT(self.beamprops.SigmaT)

            else:
                self.beam.SetDistributionType(self.beamprops.distrType)
                self.beam.SetSigmaX(self.beamprops.SigmaX, unitsstring=self.units['x'])
                self.beam.SetSigmaY(self.beamprops.SigmaY, unitsstring=self.units['y'])
                self.beam.SetSigmaXP(self.beamprops.SigmaXP, unitsstring=self.units['xp'])
                self.beam.SetSigmaYP(self.beamprops.SigmaYP, unitsstring=self.units['yp'])
                self.beam.SetSigmaE(self.beamprops.SigmaE)
                self.beam.SetSigmaT(self.beamprops.SigmaT)

            # set beam offsets in gmad if non zero
            if self.beamprops.X0 != 0:
                self.beam.SetX0(self.beamprops.X0, unitsstring=self.units['x'])
            if self.beamprops.Y0 != 0:
                self.beam.SetY0(self.beamprops.Y0, unitsstring=self.units['y'])
            if self.beamprops.Z0 != 0:
                self.beam.SetZ0(self.beamprops.Z0, unitsstring=self.units['z'])

        elif self.convprops.madxoutput:
            # calculate betas and emittances regardless for madx beam
            try:
                self.beamprops.betx = self.beamprops.SigmaX / self.beamprops.SigmaXP
            except ZeroDivisionError:
                self.beamprops.betx = 0
            try:
                self.beamprops.bety = self.beamprops.SigmaY / self.beamprops.SigmaYP
            except ZeroDivisionError:
                self.beamprops.bety = 0
                self.beamprops.emitx = self.beamprops.SigmaX * self.beamprops.SigmaXP / 1000.0
                self.beamprops.emity = self.beamprops.SigmaY * self.beamprops.SigmaYP / 1000.0

            # set madx beam
            self.beam.SetDistributionType('madx')
            self.beam.SetBetaX(self.beamprops.betx)
            self.beam.SetBetaY(self.beamprops.bety)
            self.beam.SetAlphaX(self.beamprops.alfx)
            self.beam.SetAlphaY(self.beamprops.alfy)
            self.beam.SetEmittanceX(self.beamprops.emitx / 1000)
            self.beam.SetEmittanceY(self.beamprops.emity / 1000)
            self.beam.SetSigmaE(self.beamprops.SigmaE)
            self.beam.SetSigmaT(self.beamprops.SigmaT)

        self.machine.AddBeam(self.beam)

    def _NewMachines(self):
        """
        Delete the machine and set to be the empty copied at class instantiation.
        """
        del self.machine
        self.machine = self._machineCopy
