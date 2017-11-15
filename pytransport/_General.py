import numpy as _np
from scipy import constants as _con
import sys as _sys
import os as _os
import string as _string
import glob as _glob
import Reader as _reader


class _beamprops:
    """
    A class containing the properties of the inital beam distribution.
    """
    def __init__(self, p_mass=938.272):
        # beam properties that are updated along the lattice
        self.momentum = 0
        self.k_energy = 0
        self.tot_energy_current = p_mass
        self.gamma = 1
        self.beta = 0
        self.brho = 0
        # beam properties that are from the initial beam and fixed
        self.mass = p_mass
        self.tot_energy = p_mass  # initial energy
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


class _machineprops:
    """
    A class containing the number of elements and angular properties (i.e bending direction)
    """
    def __init__(self):
        self.benddef        = True  # True = dipole defined by 4. L B n. False = dipole defined by 4. L angle n.
        self.bending        = 1     # +VE = bends to the right for positive particles
        self.angle          = 0     # dipole rotation angle
        self.drifts         = 0     # nr of drifts
        self.dipoles        = 0
        self.rf             = 0
        self.quads          = 0
        self.sextus         = 0
        self.transforms     = 0
        self.solenoids      = 0
        self.collimators    = 0
        self.beampiperadius = 20
        self.fringeIntegral = 0
        self.dipoleVertAper = 0
        self.apertureType   = 'circular'
        self._totalAccVoltage = 0
        self._e_gain_prev   = 0


class _conversionProps:
    """
    A class containing the settings for the conversion.
    """
    def __init__(self, inputfile,
                 particle      = 'proton',
                 debug         = False,
                 gmad          = True,
                 gmadDir       = 'gmad',
                 madx          = False,
                 madxDir       = 'madx',
                 auto          = True,
                 dontSplit     = False,
                 keepName      = False,
                 combineDrifts = False,
                 outlog        = True):

        self.debug = debug
        self.outlog = outlog
        self.typeCode6IsTransUpdate = True  # Definition of type code 6, true is transform update, false is collimator
        self.isAccSequence = False  # Definition of type code 11 is not an accelerator sequence

        # Automatic writing and machine splitting
        self.auto = auto
        self.dontSplit = dontSplit

        # beam definition
        self.particle = particle
        self.beamdefined = False
        self.correctedbeamdef = False

        # File input and output
        self.file = inputfile
        self.fileloaded = False
        self.gmadoutput = gmad
        self.gmadDir = gmadDir
        self.madxoutput = madx
        self.madxDir = madxDir
        self.numberparts = -1

        # transport optics output is modified to be a single line
        self.singleLineOptics = False
        self.keepName = keepName
        self.combineDrifts = combineDrifts


class _Registry:
    def __init__(self):
        self.elements = []
        self.names = []
        self.lines = []
        self.length = []
        self._uniquenames = []
        self._totalLength = 0

    def AddToRegistry(self, linedict, line):
        if not isinstance(linedict, dict):
            raise TypeError("Added element is not a Dictionary")
        self.elements.append(linedict)
        self.names.append(linedict['name'])
        if not linedict['name'] in self._uniquenames:
            self._uniquenames.append(linedict['name'])

        self.lines.append(line)
        # Cumulative length
        length = round(linedict['length'], 5)
        if len(self.length) > 0:
            self.length.append(length + self._totalLength)
        else:
            self.length.append(length)
        self._totalLength += length

    def GetElementIndex(self, name):
        elenums = []
        if name not in self.names:
            return elenums
        else:
            # Add all elements of the same name as a single element may be
            # used multiple times.
            for index, elename in enumerate(self.names):
                if elename == name:
                    elenums.append(index)
            return elenums

    def GetElement(self, name):
        elenum = self.GetElementIndex(name)
        if isinstance(elenum, list):
            elementList = []
            for num in elenum:
                elementList.append(self.elements[num])
            return elementList
        else:
            return self.elements[elenum]

    def GetElementEndSPosition(self, name):
        elenum = self.GetElementIndex(name)
        if isinstance(elenum, list):
            elementList = []
            for num in elenum:
                elementList.append(self.length[num])
            return elementList
        else:
            return self.length[elenum]

    def GetElementStartSPosition(self, name):
        elenum = self.GetElementIndex(name)
        endS = self.GetElementEndSPosition(name)

        if isinstance(elenum, list):
            elementList = []
            for index, num in enumerate(elenum):
                element = self.elements[num]
                length = element['length']
                startS = endS[index] - length
                elementList.append(round(startS, 5))
            return elementList
        else:
            element = self.elements[elenum]
            length = element['length']
            startS = endS - length
            return round(startS, 5)

    def UpdateLength(self, linedict):
        """
        Function to increases the machines length, but does not add element data.
        This is so the S positions of named elements in the fitting registry can
        be calculated correctly.
        """
        if not isinstance(linedict, dict):
            raise TypeError("Added element is not a Dictionary")
        self._totalLength += linedict['length']


class _Writer:
    def __init__(self, debugOutput=False, writeToLog=False, logfile=''):
        self.debug = debugOutput
        self.logfile = logfile
        self.outlog = writeToLog

    def Printout(self, line, outToTerminal=True):
        if outToTerminal:
            _sys.stdout.write(line+'\n')
        if self.outlog:
            logfile = open(self.logfile, 'a')
            logfile.write(line)
            logfile.write('\n')
            logfile.close()

    def DebugPrintout(self, line):
        if self.debug:
            self.Printout(line, outToTerminal=False)

    def BeamDebugPrintout(self, beam, units):
        if not isinstance(beam, _beamprops):
            raise TypeError("Beam must be pytransport _beamprops instance.")
        self.DebugPrintout('\tBeam definition :')
        self.DebugPrintout('\tdistrType = ' + beam.distrType)
        self.DebugPrintout('\tenergy = '  + _np.str(beam.tot_energy) + ' GeV')
        self.DebugPrintout('\tSigmaX = '  + _np.str(beam.SigmaX) + ' ' + units['x'])
        self.DebugPrintout('\tSigmaXP = ' + _np.str(beam.SigmaXP) + ' ' + units['xp'])
        self.DebugPrintout('\tSigmaY = '  + _np.str(beam.SigmaY) + ' ' + units['y'])
        self.DebugPrintout('\tSigmaYP = ' + _np.str(beam.SigmaYP) + ' ' + units['yp'])
        self.DebugPrintout('\tSigmaE = '  + _np.str(beam.SigmaE))
        self.DebugPrintout('\tSigmaT = '  + _np.str(beam.SigmaT))
        self.DebugPrintout('\t(Final brho = ' + _np.str(_np.round(beam.brho, 2)) + ' Tm)')
        self.DebugPrintout('\tTwiss Params:')
        self.DebugPrintout('\tBetaX = '  + _np.str(beam.betx) + ' ' + units['beta_func'])
        self.DebugPrintout('\tBetaY = '  + _np.str(beam.bety) + ' ' + units['beta_func'])
        self.DebugPrintout('\tAlphaX = ' + _np.str(beam.alfx))
        self.DebugPrintout('\tAlphaY = ' + _np.str(beam.alfy))
        self.DebugPrintout('\tEmittx = ' + _np.str(beam.emitx) + ' ' + units['emittance'])
        self.DebugPrintout('\tEmittY = ' + _np.str(beam.emity) + ' ' + units['emittance'])

    def ElementPrepDebugPrintout(self, elementType, numElements):
        debugString = "\tEntry is a " + elementType + ", adding to the element registry as element "
        debugString += numElements + "."
        self.DebugPrintout(debugString)

    def Write(self, transport, filename):
        """
        Write the converted TRANSPORT file to disk.
        """
        if not isinstance(filename, _np.str):
            raise TypeError("Filename must be a string")
        if transport.convprops.gmadoutput:
            fname = filename + '.gmad'
            dir = transport.convprops.gmadDir
            self.Printout('Writing to file: ' + dir + "/" + fname)
            if not CheckDirExists(dir):
                _os.mkdir(dir)
            _os.chdir(dir)
            transport.gmadmachine.Write(fname)
            _os.chdir('../')
        if transport.convprops.madxoutput:
            fname = filename + '.madx'
            dir = transport.convprops.madxDir
            self.Printout('Writing to file: ' + dir + "/" + fname)
            if not CheckDirExists(dir):
                _os.mkdir(dir)
            _os.chdir(dir)
            transport.madxmachine.Write(fname)
            _os.chdir('../')


def GetTypeNum(line):
    """
    Function to extract the element type number (type code).
    Written because element types can contain alphabetical
    characters when fits are used, e.g: 5.0A. Only the number
    is required, the use of fitting does not need to be known.
    """
    eleNum = line[0]
    characNum = 0
    if len(eleNum) > 2:
        for characNum in range(len(eleNum[2:])):
            try:
                converted = _np.float(eleNum[:characNum+2])
                del converted
            except ValueError:
                break
        typeNum = _np.float(eleNum[:characNum+2])
    else:
        typeNum = _np.float(eleNum)
    return typeNum


def CheckSingleLineOutputApplied(inputfile):
    """
    Function to check if the control element that print element output in
    a single line was successfully applied. Check needed as not all versions
    of TRANSPORT can run this type code.
    """
    reader = _reader.Reader()
    flist = _reader._LoadFile(inputfile)
    optics = reader.optics._getOptics(flist)
    for element in optics:
        if element == 'IO: UNDEFINED TYPE CODE 13. 19. ;':
            return False
    return True


def RemoveSpaces(line):
    elementlist = []
    for i in line:
        if (i != '') and (i != ' '):
            elementlist.append(i)
    line = _np.array(elementlist)
    return line


def FindEndOfLine(line):  # Find the end of the line of code
    endpos = -1
    breakloop = False
    if isinstance(line, _np.str):
        for charnum, char in enumerate(line):
            if char == ';':
                endpos = charnum
                break
    elif isinstance(line, _np.ndarray):
        for index, ele in enumerate(line):
            for char in ele:
                if char == ';':
                    endpos = index
                    breakloop = True
                    break
            if breakloop:
                break
    return endpos


def RemoveIllegals(line):
    """
    Function to remove '' and stray characters from lines.
    """
    illegal = ['"', '', '(', ')']

    linelist = [element for element in line if element not in illegal]
    line = _np.array(linelist)
    return line


def CheckIsSentinel(line):
    for element in line:
        if element[:8] == 'SENTINEL':
            return True
    return False


def CheckIsAddition(line, filetype='input'):
    """
    Function to check if a BEAM line of TRANSPORT code is a beam definition or r.m.s addition.
    """
    # Output file is a standard format, any RMS addition line should always be 10 long.
    if filetype == 'output':
        if len(line) == 10:
            return True

    elif filetype == 'input':
        if len(line) > 8:
            if (line[8] == '0.') or (line[8] == '0'):
                return True
            else:
                return False
    else:
        raise ValueError("File type can only be input or output")


def GetComment(line):
    """
    Function to extract a comment from a line.
    Will only find one comment per line.
    """
    concat = ''
    for ele in line:
        concat += ele
        concat += ' '
    commstart = _string.find(concat, '(')
    commend = _string.find(concat, ')')
    if commstart != -1 and commend != -1:
        comment = concat[commstart+1: commend]
        gmadcomment = '! '+comment
    else:
        gmadcomment = None
    return gmadcomment


def GetLabel(line):
    """
    Function to get element label from code line.
    """
    label = None
    for elenum, ele in enumerate(line):
        startslash = _string.find(ele, "/")
        startquote = _string.find(ele, "'")
        startequal = _string.find(ele, "=")
        startdbquote = _string.find(ele, '"')
        if startslash != -1:
            end = 1 + startslash + _string.find(ele[(startslash + 1):], "/")
            if end <= startslash:
                label = ele[startslash + 1:]
            else:
                label = ele[startslash + 1:end]
            break
        elif startquote != -1:
            end = 1 + startquote + _string.find(ele[(startslash + 1):], "'")
            if end <= startquote:
                label = ele[startquote + 1:]
            else:
                label = ele[startquote + 1:end]
            break
        elif startequal != -1:
            end = 1 + startequal + _string.find(ele[(startslash + 1):], "=")
            if end <= startequal:
                label = ele[startequal + 1:]
            else:
                label = ele[startequal + 1:end]
            break
        elif startdbquote != -1:
            end = 1 + startdbquote + _string.find(ele[(startdbquote + 1):], '"')
            if end <= startdbquote:
                label = ele[startdbquote + 1:]
            else:
                label = ele[startdbquote + 1:end]
            break
        else:
            label = None
    return label


def CheckDirExists(directory):
    dirs = _glob.glob('*/')
    if directory[-1] != '/':
        directory += '/'
    if directory in dirs:
        return True
    return False


def RemoveLabel(line):
    """
    Function to remove the label from a line.
    """
    label, elenum = GetLabel(line)
    if label is not None:
        element = line[elenum]
        lablen = len(label)
        newval = element
        for index in range(len(element)):
            if element[index:index + lablen] == label:
                prelabel = element[:index - 1]
                postlabel = element[index + lablen + 1:]
                newval = prelabel + ' ' + postlabel
                break
        line[elenum] = newval
    return line


def GetElementData(line):
    data = []
    for index, ele in enumerate(line[1:]):
        if ele != '':
            try:
                data.append(_np.float(ele))
            except ValueError:
                pass
    return data


def GetIndicator(data):
    """
    Function to read the indicator number. Must be 0, 1, or 2, where:
        0 is a new lattice
        1 is for fitting with the first lattice (from a 0 indicator file)
        2 if for a second fitting which suppresses the first fitting.
    """
    indc = 0
    linenum = 0
    for linenum, line in enumerate(data):
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
    return indc, linenum


def JoinSplitLines(linenum, lattice):
    firstline = lattice[linenum].replace(';', '')
    latticeline = firstline  # Copy for later
    firstline = _np.array(firstline.split(' '), dtype=_np.str)
    firstline = RemoveIllegals(firstline)
    numericals = []
    nonnumericals = []
    # Keep entries that are strings of numbers
    for i in firstline:
        try:
            number = _np.float(i)
            numericals.append(_np.str(number))
        except ValueError:
            nonnumericals.append(i)

    # Number of numerical elements minus the first which should be the entry type number.
    # This is bascially a way of extracting any label or comments.
    numelements = len(numericals) - 1

    secline = lattice[linenum + 1].replace(';', '')
    secline = _np.array(secline.split(' '), dtype=_np.str)
    secline = RemoveIllegals(secline)
    secnumericals = []

    for i in secline:
        try:
            number = _np.float(i)
            secnumericals.append("%.4f" % number)
        except ValueError:
            pass

    # Second line should be 15 minus number of numerical elements from prev line.
    # This is done to skip erroneous numbers in the line such as '000' which have
    # appeared when lines have been split.
    secline = secnumericals[-15 + numelements:]
    numericals.extend(secline)

    # Add to latticeline so as to appear like one single line in the file
    seclinetxt = ""
    for i in secline:
        newline = "     " + i
        seclinetxt += newline
    latticeline += seclinetxt

    # Add name to output
    if len(nonnumericals) == 1:
        numericals.append(nonnumericals[0])
    line = _np.array(numericals)
    return latticeline, line


def GetFaceRotationAngles(data, linenum):

    def searchForAngle(linelist):
        angle = 0
        for line in linelist:
            try:
                elecode = _np.float(line[0])
                del elecode
            except ValueError:
                angle = 0
                break

            if _np.float(line[0]) == 4.0:
                break
            elif _np.float(line[0]) == 2.0:
                endof = FindEndOfLine(line[1])
                if endof != -1:
                    try:
                        angle = _np.round(_np.float(line[1][:endof]), 4)
                    except ValueError:
                        try:
                            angle = _np.round(_np.float(line[2][:endof]), 4)
                        except ValueError:
                            pass
                    else:
                        pass
                else:
                    try:
                        angle = _np.round(_np.float(line[1]), 4)
                    except ValueError:
                        try:
                            angle = _np.round(_np.float(line[2]), 4)
                        except ValueError:
                            pass
                    else:
                        pass
                break
            else:
                pass
        return angle

    lineList = [i for i in data[(linenum - 5):linenum]]
    lineList.reverse()  # Search for input poleface in reverse line order

    anglein = searchForAngle(lineList)
    angleout = searchForAngle(data[linenum + 1:(linenum + 6)])

    return anglein, angleout


def _get_preamble(data):  # Redundant until pybdsim can handle comments.
    """
    Function to read any preamble at the start of the TRANSPORT file.
    """
    # indc, linenum = GetIndicator(data)
    gmadpreamble = []
    # for line in self.Transport.data[:linenum-1]:
    for line in data:
        if line == '\r\n':
            pass
        else:
            gmadline = '!' + line
            gmadpreamble.append(gmadline)
    return gmadpreamble


def _process_fits(fits):  # redundant
    # First split the fitting output into its respective sections (input problem step).
    fitsections = []
    fitsstarts = []
    # Start line of each section
    for linenum, line in enumerate(fits):
        if line.startswith('1'):
            fitsstarts.append(linenum)

    for secnum in range(len(fitsstarts)):
        if secnum + 1 < len(fitsstarts):
            section = fits[fitsstarts[secnum]:fitsstarts[secnum + 1]]
        else:
            section = fits[fitsstarts[secnum]:]
        lines = []
        for line in section:
            lines.append(RemoveIllegals(line.split(' ')))
        fitsections.append(lines)

    magnetlines = []
    for section in fitsections:
        for line in section:
            if (len(line) > 0) and (line[0][0] == '*' and line[0][-1] == '*') and line[0] != '*FIT*':
                magnetlines.append(line)


def OutputFitsToRegistry(transport, outputdata):
    isLegal = {'*DRIFT*': 3.0,
               '*QUAD*': 5.0,
               '*BEND*': 4.0}

    for line in outputdata:
        append = False
        linedict = {'elementnum': 0.0,
                    'name': '',
                    'length': 0.0}
        data = RemoveIllegals(line.split(' '))
        eledata = GetElementData(data)
        label = GetLabel(data)
        if data[0] in isLegal:
            linedict['elementnum'] = isLegal[data[0]]
            linedict['name'] = label
            linedict['data'] = eledata[1:]  # first value is elementnum.
            linedict['length'] = eledata[1]
            append = True

        # Only add an element with a name to the fitting registry.
        # (Element has to be named to be varied in the fitting routine).
        # Otherwise update the total length of the machine.
        if append and (label is not None) and (label != ''):
            transport.FitRegistry.AddToRegistry(linedict, line)
        else:
            transport.FitRegistry.UpdateLength(linedict)
    return transport


def ConvertBunchLength(transport, bunch_length):
    """
    Function to convert bunch length unit in TRANSPORT into seconds.
    """
    scale = transport.scale[transport.units['bunch_length'][0]]
    blmeters = bunch_length * scale  # Bunch length scaled to metres
    blseconds = blmeters / (transport.beamprops.beta*_con.c)  # Length converted to seconds
    return blseconds


def UpdateMomentumFromEnergy(transport, k_energy):
    """
    Function to calculate (from kinetic energy):
        Total Energy
        Kinetic Energy
        Momentum
        Lorentz factor (gamma)
        Velocity (beta)
        Magnetic rigidity (brho)
    """
    scaling = 1  # defaults
    mom_in_ev = 0
    k_energy = _np.float(k_energy)

    transport.beamprops.k_energy = k_energy
    p_mass = transport.beamprops.mass  # Particle rest mass (in GeV)

    e_unit = transport.units['p_egain']
    if e_unit != 'eV':
        scaling = 1e9 / transport.scale[e_unit[0]]  # Scaling relative to mom. unit
    elif e_unit == 'eV':
        scaling = 1e9  # Scaling relative to 1 eV
    p_mass *= scaling  # Scale particle rest mass

    # energy = _np.sqrt((p_mass**2 * _con.c**2) + (momentum**2 * _con.c**2)) / _con.c

    transport.beamprops.tot_energy_current = k_energy + p_mass
    transport.beamprops.momentum = _np.sqrt((transport.beamprops.tot_energy_current ** 2) - (p_mass ** 2))

    transport.beamprops.gamma = transport.beamprops.tot_energy_current / p_mass
    transport.beamprops.beta = _np.sqrt((1.0 - (1.0 / transport.beamprops.gamma ** 2)))

    if e_unit != 'eV':
        mom_in_ev = transport.beamprops.momentum * transport.scale[e_unit[0]]
    elif e_unit == 'eV':
        mom_in_ev = transport.beamprops.momentum

    transport.beamprops.brho = mom_in_ev / _con.c
    return transport


def UpdateEnergyFromMomentum(transport, momentum):
    """
    Function to calculate (from momentum):
        Total Energy
        Kinetic Energy
        Momentum
        Lorentz factor (gamma)
        Velocity (beta)
        Magnetic rigidity (brho)
    """
    momentum = _np.float(momentum)
    transport.beamprops.momentum = momentum
    p_mass = transport.beamprops.mass  # Particle rest mass (in GeV)
    scaling = 1
    mom_in_ev = momentum

    mom_unit = transport.units['p_egain']
    if mom_unit != 'eV':
        scaling = 1e9 / transport.scale[mom_unit[0]]     # Scaling relative to mom. unit
        mom_in_ev = momentum * transport.scale[mom_unit[0]]
    elif mom_unit == 'eV':
        scaling = 1e9                               # Scaling relative to 1 eV
        mom_in_ev = momentum
    p_mass *= scaling                               # Scale particle rest mass
    energy = _np.sqrt((p_mass**2) + (momentum**2))
    transport.beamprops.tot_energy = energy
    transport.beamprops.tot_energy_current = energy
    transport.beamprops.k_energy = energy - p_mass
    transport.beamprops.gamma = energy / p_mass
    transport.beamprops.beta = _np.sqrt((1.0 - (1.0 / transport.beamprops.gamma**2)))
    transport.beamprops.brho = mom_in_ev / _con.c
    return transport


def ScaleToMeters(transport, quantity):
    """
    Function to scale quantity (string) to meters, returns conversion factor.
    """
    if transport.units[quantity] != 'm':
        conversionFactor = transport.scale[transport.units[quantity][0]]
    else:
        conversionFactor = 1
    return conversionFactor


def CheckIsOutput(inputfile):
    """
    Function to check if a file is a standard TRANSPORT output file.
    Based upon existence of the lines:
        "0  XXX"
    being present, which represents the TRANSPORT indicator card line.
    X can be 0, 1, 2. Default is 0.
    """
    temp = _reader.reader()
    isOutput = False
    try:
        f = open(inputfile)
        for inputline in f:
            inputline = inputline.replace("\r", '')
            inputline = inputline.replace("\n", '')
            if inputline in temp._allowedIndicatorLines:
                isOutput = True
                break
        f.close()
    except IOError:
        raise IOError('Cannot open file.')
    return isOutput


def RemoveFileExt(inputfile):
    """
    Remove the file extension from the input file name. Only works on known extensions.
    """
    exts = [".txt", ".dat"]
    if inputfile[-4:] in exts:
        return inputfile[:-4]
    return inputfile
