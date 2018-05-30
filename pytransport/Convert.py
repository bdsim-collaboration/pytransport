# pytransport.Elements - tools for Transport element conversion
# Version 1.0
# W. Shields and J. Snuverink
# william.shields.2010@live.rhul.ac.uk

"""
Wrapper function for automatically converting to bdsim and/or madx

Class containing functions to help convert elements from Transport to gmad/madx.
Class instantiated in pybdsim.Convert.Transport2Gmad and pymadx.Convert.Transport2Madx.

Classes:
_Convert - a class used to convert different element types.

"""

import numpy as _np
import string as _str

import _General
from _General import _Writer
from Data import ConversionData as _convData
from pytransport import Reader as _Reader


def Convert(inputfile,
            particle='proton',
            distrType='gauss',
            output='bdsim',
            outputDir='',
            debug=False,
            dontSplit=False,
            keepName=False,
            combineDrifts=False,
            options=None,
            machine=None):
    """
    **Convert** convert a Transport input or output file into an input file for bdsim, madx, or both

    +-------------------------------+-------------------------------------------------------------------+
    | **inputfile**                 | dtype = string                                                    |
    |                               | path to the input file                                            |
    +-------------------------------+-------------------------------------------------------------------+
    | **particle**                  | dtype = string. Optional, default = "proton"                      |
    |                               | the particle species                                              |
    +-------------------------------+-------------------------------------------------------------------+
    | **distrType**                 | dtype = string. Optional, Default = "gauss".                      |
    |                               | the beam distribution type. Can be either gauss or gausstwiss.    |
    +-------------------------------+-------------------------------------------------------------------+
    | **output**                    | dtype=string. Optional, default = "bdsim"                         |
    |                               | the output type, can be "bdsim", "madx", or "both"                |
    +-------------------------------+-------------------------------------------------------------------+
    | **outputDir**                 | dtype=string. Optional, default = ""                              |
    |                               | the output directory where the files will be written              |
    |                               | if no input supplied, defaults to output type. if outputDir is    |
    |                               | supplied and output is "both", outputDir is appended with the     |
    |                               | output type, e.g. outputDir_bdsim                                 |
    +-------------------------------+-------------------------------------------------------------------+
    | **debug**                     | dtype = bool. Optional, default = False                           |
    |                               | output a log file (inputfile_conversion.log) detailing the        |
    |                               | conversion process, element by element                            |
    +-------------------------------+-------------------------------------------------------------------+
    | **dontSplit**                 | dtype = bool. Optional, default = False                           |
    |                               | the converter splits the machine into multiple parts when a beam  |
    |                               | is redefined in a Transport lattice. dontSplit overrides this and |
    |                               | forces the machine to be written to a single file                 |
    +-------------------------------+-------------------------------------------------------------------+
    | **keepName**                  | dtype = bool. Optional, default = False                           |
    |                               | keep the names of elements as defined in the Transport inputfile. |
    |                               | Appends element name with _N where N is an integer if the element |
    |                               | name has already been used                                        |
    +-------------------------------+-------------------------------------------------------------------+
    | **combineDrifts**             | dtype = bool. Optional, default = False                           |
    |                               | combine multiple consecutive drifts into a single drift           |
    +-------------------------------+-------------------------------------------------------------------+
    | **machine**                   | dtype = pybdsim.Builder.Machine or pymadx.Builder.Machine         |
    |                               | required, default = None                                          |
    |                               | machine instance required for conversion.                         |
    +-------------------------------+-------------------------------------------------------------------+
    | **options**                   | dtype = pybdsim.Options.Options. Optional, default = None         |
    |                               | options instance required to write options in bdsim conversion.   |
    |                               | Ignored if converting to "madx" format.                           |
    +-------------------------------+-------------------------------------------------------------------+

    Example:

    >>> Convert(inputfile, machine=pybdsim.Builder.Machine(), options=pybdsim.Builder.Options())
    >>> Convert(inputfile, output="madx", machine=pymadx.Builder.Machine())

    Writes converted machine to disk. Reader automatically detects if the supplied input file is a Transport input
    file or Transport output file.

    """
    outputType = _str.lower(output)

    if machine is None:
        raise TypeError("machine instance must be supplied")
    if ((outputType == 'bdsim') or (output == 'both')) and (options is None):
        raise TypeError("pybdsim.Options.Options must be supplied for bdsim conversion")
    if (outputType == 'madx') and (options is not None):
        print("Ignoring supplied options for madx conversion")

    # default for output='bdsim'
    gmadDir = outputDir
    madxDir = outputDir

    if outputType == 'bdsim':
        gmad = True
        madx = False
        if outputDir == '':
            gmadDir = 'bdsim'
            madxDir = ''
    elif outputType == 'madx':
        gmad = False
        madx = True
        if outputDir == '':
            gmadDir = ''
            madxDir = 'madx'
    elif outputType == 'both':
        gmad = True
        madx = True
        if outputDir == '':
            gmadDir = 'bdsim'
            madxDir = 'madx'
        else:
            gmadDir += '_bdsim'
            madxDir += '_madx'
    else:
        raise IOError("Unknown output type '"+output+"'")

    converter = _Convert(_convData(inputfile=inputfile,
                                   particle=particle,
                                   distrType=distrType,
                                   gmad=gmad,
                                   madx=madx,
                                   gmadDir=gmadDir,
                                   madxDir=madxDir,
                                   debug=debug,
                                   dontSplit=dontSplit,
                                   keepName=keepName,
                                   combineDrifts=combineDrifts,
                                   options=options,
                                   machine=machine))
    converter.Convert()


class _Convert:
    """
    Class for converting different Transport element types into gmad/madx format.

    Required: transportData, dtype = pytransport.Data.ConversionData

    """
    def __init__(self, transportData):
        if not isinstance(transportData, _convData):
            raise TypeError("transportData must be a pytransport.Data.ConversionData instance")
        self.Transport = transportData
        logfileName = _General.RemoveFileExt(self.Transport.convprops.file) + '_conversion.log'
        self.Writer = _Writer(debugOutput=self.Transport.convprops.debug,
                              writeToLog=self.Transport.convprops.outlog,
                              logfile=logfileName)

    def LoadFile(self, inputfile):
        """
        Load a Transport file.
        """
        if not isinstance(inputfile, _np.str):
            raise TypeError("Input must be a string")

        infile = inputfile.split('/')[-1]  # Remove filepath, leave just filename
        self.Transport._file = infile[:-4]  # Remove extension
        self.Transport._filename = inputfile
        isOutput = _General.CheckIsOutput(inputfile)  # Is a TRANSPORT standard output file.
        self.Writer.DebugPrintout("File Read.")
        if isOutput:
            lattice = _Reader.GetLattice(inputfile)
            fitres = _Reader.GetResultsFromFitting(inputfile)
            self.Transport = _General.OutputFitsToRegistry(self.Transport, fitres)
            self.Writer.DebugPrintout('Adding any fitting output to the fitting registry (self.FitRegistry)')

            self.Writer.DebugPrintout('Processing file and adding to Transport class.')
            for linenum, latticeline in enumerate(lattice):
                latticeline = latticeline.replace(';', '')
                line = _np.array(latticeline.split(' '), dtype=_np.str)
                line = _General.RemoveIllegals(line)

                # Method of dealing with split lines in the output
                # Should only be applicable to type 12 entry (up to 15 variables)
                # It is assumed that the line is always split, so be careful.
                prevline = lattice[linenum - 1].replace(';', '')
                prevline = _np.array(prevline.split(' '), dtype=_np.str)
                prevline = _General.RemoveIllegals(prevline)

                try:
                    if (linenum > 0) and _np.abs(_np.float(line[0])) == 12.0:
                        latticeline, line = _General.JoinSplitLines(linenum, lattice)
                    # Ignore line after type 12 entry (second part of split line)
                    if (linenum > 1) and _np.abs(_np.float(prevline[0])) == 12.0:
                        pass
                    else:
                        self.Transport.data.append(line)
                        self.Transport.filedata.append(latticeline)
                except ValueError:
                    self.Transport.data.append(line)
                    self.Transport.filedata.append(latticeline)
                except IndexError:
                    pass

        else:
            f = open(inputfile)
            for inputline in f:
                endoflinepos = _General.FindEndOfLine(inputline)
                templine = inputline
                if endoflinepos > 0:
                    templine = inputline[:endoflinepos]
                line = _np.array(templine.split(' '), dtype=_np.str)
                # do not change comment lines
                if not line[0][0] == '(':
                    line = _General.RemoveIllegals(line)
                self.Transport.data.append(line)
                self.Transport.filedata.append(inputline)
            f.close()
        self.Transport.convprops.fileloaded = True

    def Write(self):
        """
        Write the converted TRANSPORT file to disk.
        """
        self.Writer.DebugPrintout("Adding beam to gmad machine:")
        self.Transport.AddBeam()
        if self.Transport.convprops.gmadoutput:
            self.Transport.AddOptions()
        self.Transport.machine.AddSampler('all')
        self.Writer.BeamDebugPrintout(self.Transport.beamprops, self.Transport.units)
        fname = _General.RemoveFileExt(self.Transport.convprops.file)
        if self.Transport.convprops.numberparts < 0:
            filename = fname
        else:
            self.Transport.convprops.numberparts += 1
            filename = fname + '_part' + _np.str(self.Transport.convprops.numberparts)
        self.Writer.Write(self.Transport, filename)

    def Convert(self):
        """
        Convert, process, and write.
        """
        self.LoadFile(self.Transport.convprops.file)

        if not self.Transport.convprops.fileloaded:
            self.Writer.Printout('No file loaded.')
            return
        self.ProcessAndBuild()
        self.Write()

    def ProcessAndBuild(self):
        """
        Function that loops over the lattice, adds the elements to the element registry,
        and updates any elements that have fitted parameters.
        It then converts the registry elements and adds to the gmad machine.
        """
        self.Writer.DebugPrintout('Processing tokenised lines from input file and adding to element registry.\n')

        for linenum, line in enumerate(self.Transport.data):
            self.Writer.DebugPrintout('Processing tokenised line ' + _np.str(linenum) + ' :')
            self.Writer.DebugPrintout('\t' + str(line))
            self.Writer.DebugPrintout('\tOriginal :')
            self.Writer.DebugPrintout('\t' + self.Transport.filedata[linenum])

            # Checks if the SENTINEL line is found. SENTINEL relates to TRANSPORT fitting routine and is only written
            # after the lattice definition, so there's no point reading lines beyond it.
            if _General.CheckIsSentinel(line):
                self.Writer.DebugPrintout('Sentinel Found.')
                break
            # Test for positive element, negative ones ignored in TRANSPORT so ignored here too.
            try:
                typeNum = _General.GetTypeNum(line)
                if typeNum > 0:
                    if self.Transport.data[0][0] == 'OUTPUT':
                        self._ElementPrepper(line, linenum, 'output')
                    else:
                        line = _General.RemoveIllegals(line)
                        self._ElementPrepper(line, linenum, 'input')
                else:
                    self.Writer.DebugPrintout('\tType code is 0 or negative, ignoring line.')
            except ValueError:
                errorline = '\tCannot process line ' + _np.str(linenum) + ', '
                if line[0][0] == '(' or line[0][0] == '/':
                    errorline += 'line is a comment.'
                elif line[0][0] == 'S':  # S used as first character in SENTINEL command.
                    errorline = 'line is for TRANSPORT fitting routine.'
                elif line[0] == '\n':
                    errorline = 'line is blank.'
                else:
                    errorline = 'reason unknown.'
                self.Writer.DebugPrintout(errorline)
        self._UpdateElementsFromFits()
        self.Writer.DebugPrintout(
            'Converting registry elements to pybdsim compatible format and adding to machine builder.\n')

        skipNextDrift = False  # used for collimators
        lastElementWasADrift = True  # default value
        if self.Transport.convprops.combineDrifts:
            lastElementWasADrift = False
        for linenum, linedict in enumerate(self.Transport.ElementRegistry.elements):
            self.Writer.DebugPrintout('Converting element number ' + _np.str(linenum) + ':')
            convertline = '\t'
            for keynum, key in enumerate(linedict.keys()):
                if keynum != 0:
                    convertline += ', '
                if key == 'data':
                    convertline += 'element data:'
                    for ele in linedict[key]:
                        convertline += ('\t' + _np.str(ele))
                else:
                    convertline += (key + ': ' + _np.str(linedict[key]))
                if keynum == len(linedict.keys()):
                    convertline += '.'
            self.Writer.DebugPrintout(convertline)

            if self.Transport.convprops.combineDrifts:
                if lastElementWasADrift and linedict['elementnum'] != 3.0 and linedict['elementnum'] < 20.0:
                    # write possibly combined drift
                    self.Writer.DebugPrintout('\n\tConvert delayed drift(s)')
                    self.Drift(linedictDrift)
                    lastElementWasADrift = False
                    self.Writer.DebugPrintout('\n\tNow convert element number' + _np.str(linenum))

            if linedict['elementnum'] == 15.0:
                self.UnitChange(linedict)
            if linedict['elementnum'] == 20.0:
                self.ChangeBend(linedict)
            if linedict['elementnum'] == 1.0:  # Add beam on first definition
                if not self.Transport.convprops.beamdefined:
                    self.DefineBeam(linedict)
                elif not self.Transport.convprops.dontSplit:  # Only update beyond first definition if splitting is permitted
                    self.DefineBeam(linedict)
            if linedict['elementnum'] == 3.0:
                if skipNextDrift:
                    skipNextDrift = False
                    continue
                if self.Transport.convprops.combineDrifts:
                    self.Writer.DebugPrintout('\tDelay drift')
                    if lastElementWasADrift:
                        linedictDrift['length'] += linedict['length']  # update linedictDrift
                        if not linedictDrift['name']:
                            linedictDrift['name'] = linedict['name']  # keep first non-empty name
                    else:
                        linedictDrift = linedict  # new linedictDrift
                        lastElementWasADrift = True
                else:
                    self.Drift(linedict)
            if linedict['elementnum'] == 4.0:
                self.Dipole(linedict)
            if linedict['elementnum'] == 5.0:
                self.Quadrupole(linedict)
            if linedict['elementnum'] == 6.0:
                if not self.Transport.convprops.typeCode6IsTransUpdate:
                    self.Collimator(linedict)
                    # Length gotten from next drift
                    if linedict['length'] > 0.0:
                        skipNextDrift = True
                else:
                    self.TransformUpdate(linedict)
            if linedict['elementnum'] == 12.0:
                self.Correction(linedict)
            if linedict['elementnum'] == 11.0:
                self.Acceleration(linedict)
            if linedict['elementnum'] == 13.0:
                self.Printline(linedict)
            if linedict['elementnum'] == 16.0:
                self.SpecialInput(linedict)
            if linedict['elementnum'] == 18.0:
                self.Sextupole(linedict)
            if linedict['elementnum'] == 19.0:
                self.Solenoid(linedict)

            # 9.  : 'Repetition' - for nesting elements
            if linedict['elementnum'] == 9.0:
                self.Writer.DebugPrintout('\tWARNING Repetition Element not implemented in converter!')
            if linedict['elementnum'] == 2.0:
                errorline = '\tLine is a poleface rotation which is handled by the previous or next dipole element.'
                self.Writer.DebugPrintout(errorline)

            self.Writer.DebugPrintout('\n')

        # OTHER TYPES WHICH CAN BE IGNORED:
        # 6.0.X : Update RX matrix used in TRANSPORT
        # 7.  : 'Shift beam centroid'
        # 8.  : Magnet alignment tolerances
        # 10. : Fitting constraint
        # 14. : Arbitrary transformation of TRANSPORT matrix
        # 22. : Space charge element
        # 23. : RF Cavity (Buncher), changes bunch energy spread

        # Write also last drift
        if self.Transport.convprops.combineDrifts:
            if lastElementWasADrift:
                self.Writer.DebugPrintout('\tConvert delayed drift(s)')
                self.Drift(linedictDrift)

    def _ElementPrepper(self, line, linenum, filetype='input'):
        """
        Function to extract the data and prepare it for processing by each element function.
        This has been written as the lattice lines from an input file and output file are different,
        so it just a way of correctly ordering the information.
        """
        linedict = {'elementnum': 0.0,
                    'name': '',
                    'length': 0.0,
                    'isZeroLength': True}
        numElements = _np.str(len(self.Transport.ElementRegistry.elements))
        typeNum = _General.GetTypeNum(line)
        linedict['elementnum'] = typeNum

        if typeNum == 15.0:
            label = _General.GetLabel(line)
            if filetype == 'output':
                linedict['label'] = line[2].strip('"')
            if filetype == 'input':
                linedict['label'] = label
            linedict['number'] = line[1]
            self.Writer.ElementPrepDebugPrintout("Unit Control", numElements)

        if typeNum == 20.0:
            angle = 0  # Default
            if len(line) == 2:  # i.e has a label
                endofline = _General.FindEndOfLine(line[1])
                angle = line[1][:endofline]
            else:
                for index in line[1:]:
                    try:
                        angle = _np.float(index)
                        break
                    except ValueError:
                        pass
            linedict['angle'] = angle
            self.Writer.ElementPrepDebugPrintout("coordinate rotation", numElements)

        if typeNum == 1.0:
            linedict['name'] = _General.GetLabel(line)
            linedict['isAddition'] = False
            if _General.CheckIsAddition(line, filetype):
                linedict['isAddition'] = True
            # line = self._remove_label(line)
            if len(line) < 8:
                raise IndexError("Incorrect number of beam parameters.")
            n = 1
            if filetype == 'input':
                n = 0
            elif filetype == 'output':
                n = 1

            linedict['momentum'] = line[7 + n]
            linedict['Sigmax'] = line[1 + n]
            linedict['Sigmay'] = line[3 + n]
            linedict['Sigmaxp'] = line[2 + n]
            linedict['Sigmayp'] = line[4 + n]
            linedict['SigmaT'] = line[5 + n]
            linedict['SigmaE'] = line[6 + n]
            self.Writer.ElementPrepDebugPrintout("Beam definition or r.m.s addition", numElements)

        if typeNum == 2.0:
            linedict['name'] = _General.GetLabel(line)
            linedict['data'] = _General.GetElementData(line)
            self.Writer.ElementPrepDebugPrintout("poleface rotation", numElements)

        if typeNum == 3.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("drift", numElements)

        if typeNum == 4.0:
            linedict['name'] = _General.GetLabel(line)
            linedict['linenum'] = linenum
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            e1, e2 = _General.GetFaceRotationAngles(self.Transport.data, linenum)
            linedict['e1'] = e1
            linedict['e2'] = e2
            self.Writer.ElementPrepDebugPrintout("dipole", numElements)

        if typeNum == 5.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("quadrupole", numElements)

        if typeNum == 6.0:
            # element is a collimator or transform update
            # transform update is later ignored so only update linedict as if collimator

            physicalElements = [1.0, 3.0, 4.0, 5.0, 11.0, 18.0, 19.0]

            # Only iterate if not the last element
            if linenum == len(self.Transport.data):
                pass
            else:
                # Since collimators have zero length in TRANSPORT, chosen to use length of next drift instead if
                # present. Check all remaining elements for the next drift, following element(s) may be non-physical
                # element which can be ignored as it shouldnt affect the beamline. Ignore beam definition too in
                # the case where machine splitting is not permitted.
                for nextline in self.Transport.data[linenum + 1:]:
                    nextTypeNum = _General.GetTypeNum(nextline)
                    if nextTypeNum == 3.0:
                        nextData = _General.GetElementData(nextline)
                        linedict['length'] = nextData[0]
                        linedict['isZeroLength'] = False
                        linedict['name'] = _General.GetLabel(line)
                        data = _General.GetElementData(line)
                        linedict['data'] = data
                        break
                    # stop if physical element or beam redef if splitting permitted
                    elif nextTypeNum in physicalElements:
                        if (nextTypeNum == 1.0) and self.Transport.convprops.dontSplit:
                            pass
                        elif (nextTypeNum == 6.0) and self.Transport.convprops.typeCode6IsTransUpdate:
                            pass
                        else:
                            break
                    # ignore non-physical element
                    else:
                        pass
            # Can be either transform update or collimator, a 16. 14 entry changes the definition but is only
            # processed after ALL elements are added to the registry.
            self.Writer.ElementPrepDebugPrintout("transform update or collimator", numElements)

        if typeNum == 9.0:
            self.Writer.ElementPrepDebugPrintout("repetition control", numElements)

        if typeNum == 11.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['voltage'] = data[1]
            linedict['isZeroLength'] = False
            if len(data) == 4:  # Older case for single element
                linedict['phase_lag'] = data[2]
                linedict['wavel'] = data[3]
            self.Writer.ElementPrepDebugPrintout("acceleration element", numElements)

        if typeNum == 12.0:
            linedict['data'] = _General.GetElementData(line)
            linedict['name'] = _General.GetLabel(line)

            prevline = self.Transport.data[linenum - 1]  # .split(' ')
            linedict['prevlinenum'] = _np.float(prevline[0])
            linedict['isAddition'] = False
            if _General.CheckIsAddition(line):
                linedict['isAddition'] = True
            self.Writer.ElementPrepDebugPrintout("beam rotation", numElements)

        if typeNum == 13.0:
            linedict['data'] = _General.GetElementData(line)
            self.Writer.ElementPrepDebugPrintout("Input/Output control", numElements)

        if typeNum == 16.0:
            linedict['data'] = _General.GetElementData(line)
            self.Writer.ElementPrepDebugPrintout("special input", numElements)

        if typeNum == 18.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("sextupole", numElements)

        if typeNum == 19.0:
            linedict['name'] = _General.GetLabel(line)
            data = _General.GetElementData(line)
            linedict['data'] = data
            linedict['length'] = data[0]
            linedict['isZeroLength'] = False
            self.Writer.ElementPrepDebugPrintout("solenoid", numElements)

        if typeNum == 22.0:
            self.Writer.ElementPrepDebugPrintout("space charge element", numElements)

        if typeNum == 23.0:
            self.Writer.ElementPrepDebugPrintout("buncher", numElements)

        rawline = self.Transport.filedata[linenum]
        self.Transport.ElementRegistry.AddToRegistry(linedict, rawline)

    def DefineBeam(self, linedict):
        if linedict['isAddition']:
            if self.Transport.convprops.debug:
                self.Writer.DebugPrintout('\tIgnoring beam rms addition.')
            return
        if self.Transport.convprops.beamdefined and not self.Transport.convprops.dontSplit:
            self.Writer.Printout('Beam redefinition found. Writing previous section to file.')
            self.Writer.Printout('Splitting into multiple machines.')
            self.Transport.convprops.numberparts += 1
            self.Transport.AddBeam()
            if self.Transport.convprops.gmadoutput:
                self.Transport.AddOptions()
            self.Transport.machine.AddSampler('all')
            self.Writer.BeamDebugPrintout(self.Transport.beamprops, self.Transport.units)
            fname = _General.RemoveFileExt(self.Transport.convprops.file)
            if self.Transport.convprops.numberparts < 0:
                filename = fname
            else:
                self.Transport.convprops.numberparts += 1
                filename = fname + '_part' + _np.str(self.Transport.convprops.numberparts)
            self.Writer.Write(self.Transport, filename)
            self.Transport.ResetMachine()
            self.Transport.convprops.correctedbeamdef = False

        momentum = linedict['momentum']

        self.Transport.convprops.beamdefined = True

        # Convert momentum to energy and set distribution params.
        self.Transport = _General.UpdateEnergyFromMomentum(self.Transport, momentum)
        self.Transport.beamprops.SigmaX  = _np.float(linedict['Sigmax'])
        self.Transport.beamprops.SigmaY  = _np.float(linedict['Sigmay'])
        self.Transport.beamprops.SigmaXP = _np.float(linedict['Sigmaxp'])
        self.Transport.beamprops.SigmaYP = _np.float(linedict['Sigmayp'])
        self.Transport.beamprops.SigmaE  = _np.float(linedict['SigmaE']) * 0.01 * (self.Transport.beamprops.beta**2)  # Convert from percentage mom spread to absolute espread
        self.Transport.beamprops.SigmaT  = _General.ConvertBunchLength(self.Transport, _np.float(linedict['SigmaT']))  # Get bunch length in seconds.

        # Calculate Initial Twiss params
        try:
            self.Transport.beamprops.betx = self.Transport.beamprops.SigmaX / self.Transport.beamprops.SigmaXP
        except ZeroDivisionError:
            self.Transport.beamprops.betx = 0
        try:
            self.Transport.beamprops.bety = self.Transport.beamprops.SigmaY / self.Transport.beamprops.SigmaYP
        except ZeroDivisionError:
            self.Transport.beamprops.bety = 0
        self.Transport.beamprops.emitx = self.Transport.beamprops.SigmaX * self.Transport.beamprops.SigmaXP / 1000.0
        self.Transport.beamprops.emity = self.Transport.beamprops.SigmaY * self.Transport.beamprops.SigmaYP / 1000.0

        self.Writer.BeamDebugPrintout(self.Transport.beamprops, self.Transport.units)

    def Drift(self, linedict):
        driftlen = linedict['length']
        if driftlen <= 0:
            self.Writer.DebugPrintout('\tZero or negative length element, ignoring.')
            return

        lenInM = driftlen * _General.ScaleToMeters(self.Transport, 'element_length')  # length in metres

        self.Transport.machineprops.drifts += 1
        elementid = self._GetTransportElementName(linedict['name'])
        if not elementid:  # check on empty string
            elementid = 'DR'+_np.str(self.Transport.machineprops.drifts)

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddDrift(name=elementid, length=lenInM)

        self.Writer.DebugPrintout('\tConverted to:')
        self.Writer.DebugPrintout('\t' + 'Drift ' + elementid + ', length ' + _np.str(lenInM) + ' m')

    def Dipole(self, linedict):
        linenum = linedict['linenum']
        dipoledata = linedict['data']
        length = dipoledata[0]  # First two non-blanks must be the entries in a specific order.

        # Get poleface rotation
        e1 = linedict['e1'] * ((_np.pi / 180.0)*self.Transport.machineprops.bending)  # Entrance pole face rotation.
        e2 = linedict['e2'] * ((_np.pi / 180.0)*self.Transport.machineprops.bending)  # Exit pole face rotation.

        if e1 != 0:
            self.Writer.DebugPrintout('\tPreceding element on line (' + _np.str(linenum-1) + ') of the inout file provides an entrance poleface rotation of ' + _np.str(_np.round(e1, 4)) + ' rad.')
        if e2 != 0:
            self.Writer.DebugPrintout('\tFollowing element on line (' + _np.str(linenum+1) + ') of the input file provides an exit poleface rotation of ' + _np.str(_np.round(e2, 4)) + ' rad.')

        # Fringe Field Integrals
        fintVal = 0
        fintxVal = 0
        fintSecVal = 0
        fintxSecVal = 0
        if e1 != 0:
            fintVal = self.Transport.machineprops.fringeIntegral
            fintSecVal = self.Transport.machineprops.secondfringeInt
        if e2 != 0:
            fintxVal = self.Transport.machineprops.fringeIntegral
            fintxSecVal = self.Transport.machineprops.secondfringeInt

        if (fintVal != 0) or (fintxVal != 0):
            self.Writer.DebugPrintout('\tA previous entry set the fringe field integral K1=' + _np.str(self.Transport.machineprops.fringeIntegral) + '.')
            self.Writer.DebugPrintout('\tA previous entry set the second fringe field integral K2=' + _np.str(self.Transport.machineprops.secondfringeInt) + '.')

        hgap = self.Transport.machineprops.dipoleVertAper * self.Transport.scale[self.Transport.units['bend_vert_gap'][0]]

        # Poleface curvatures
        h1 = self.Transport.machineprops.bendInCurvature
        h2 = self.Transport.machineprops.bendOutCurvature
        if self.Transport.machineprops.bendInCurvature != 0:
            self.Writer.DebugPrintout('\tA previous entry set the dipole poleface entrance curvature H1=' + _np.str(self.Transport.machineprops.bendInCurvature) + '.')
        if self.Transport.machineprops.bendOutCurvature != 0:
            self.Writer.DebugPrintout('\tA previous entry set the dipole poleface entrance curvature H2=' + _np.str(self.Transport.machineprops.bendOutCurvature) + '.')

        # Calculate bending angle
        if self.Transport.machineprops.benddef:
            bfield = dipoledata[1]
            field_in_Gauss = bfield * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Scale to Gauss
            field_in_Tesla = field_in_Gauss * 1e-4                                  # Convert to Tesla
            if field_in_Tesla == 0:
                angle = 0                                                           # zero field = zero angle
            else:
                rho = self.Transport.beamprops.brho / (_np.float(field_in_Tesla))         # Calculate bending radius.
                angle = (_np.float(length) / rho) * self.Transport.machineprops.bending   # for direction of bend
            self.Writer.DebugPrintout('\tbfield = ' + _np.str(field_in_Gauss) + ' kG')
            self.Writer.DebugPrintout('\tbfield = ' + _np.str(field_in_Tesla) + ' T')
            self.Writer.DebugPrintout('\tCorresponds to angle of ' + _np.str(_np.round(angle, 4)) + ' rad.')
        else:
            angle_in_deg = dipoledata[1]
            angle = angle_in_deg * (_np.pi/180.) * self.Transport.machineprops.bending

        # Convert element length
        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        self.Transport.machineprops.dipoles += 1
        elementid = self._GetTransportElementName(linedict['name'])
        if not elementid:  # check on empty string
            elementid = 'BM' + _np.str(self.Transport.machineprops.dipoles)

        # pybdsim and pymadx set differently depending on second fringe field integral. Check for non zero pole face rotation.
        if (e1 != 0) and (e2 != 0):
            if self.Transport.convprops.madxoutput:
                self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM,
                                                 angle=_np.round(angle, 4), e1=_np.round(e1, 4), e2=_np.round(e2, 4),
                                                 fint=fintVal, fintx=fintxVal, hgap=hgap)
            else:
                self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM,
                                                 angle=_np.round(angle, 4), e1=_np.round(e1, 4), e2=_np.round(e2, 4),
                                                 fint=fintVal, fintx=fintxVal, fintK2=fintSecVal, fintxK2=fintxSecVal,
                                                 hgap=hgap)
        elif (e1 != 0) and (e2 == 0):
            if self.Transport.convprops.madxoutput:
                self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM,
                                                 angle=_np.round(angle, 4), e1=_np.round(e1, 4), fint=fintVal, fintx=0,
                                                 hgap=hgap)
            else:
                self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM,
                                                 angle=_np.round(angle, 4), e1=_np.round(e1, 4), fint=fintVal, fintx=0,
                                                 fintK2=fintSecVal, fintxK2=0, hgap=hgap)
        elif (e1 == 0) and (e2 != 0):
            if self.Transport.convprops.madxoutput:
                self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM,
                                                 angle=_np.round(angle, 4), e2=_np.round(e2, 4), fint=0, fintx=fintxVal,
                                                 hgap=hgap)
            else:
                self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM,
                                                 angle=_np.round(angle, 4), e2=_np.round(e2, 4), fint=0, fintx=fintxVal,
                                                 fintK2=0, fintxK2=fintxSecVal, hgap=hgap)
        else:
            self.Transport.machine.AddDipole(name=elementid, category='sbend', length=lenInM, angle=_np.round(angle, 4))

        # Debug output
        if (e1 != 0) and (e2 != 0):
            polefacestr = ', e1= ' + _np.str(_np.round(e1, 4)) + ' rad, e2= ' + _np.str(_np.round(e2, 4)) + ' rad'
        elif (e1 != 0) and (e2 == 0):
            polefacestr = ', e1= ' + _np.str(_np.round(e1, 4)) + ' rad'
        elif (e1 == 0) and (e2 != 0):
            polefacestr = ', e2= ' + _np.str(_np.round(e2, 4)) + ' rad'
        else:
            polefacestr = ''

        if (fintVal != 0) and (fintxVal != 0):
            fringestr = ' , fint= ' + _np.str(fintVal) + ', fintx= ' + _np.str(fintxVal)
        elif (fintVal != 0) and (fintxVal == 0):
            fringestr = ' , fint= ' + _np.str(fintVal)
        elif (fintVal == 0) and (fintxVal != 0):
            fringestr = ' , fintx= ' + _np.str(fintxVal)
        else:
            fringestr = ''

        if (fintSecVal != 0) and (fintxSecVal != 0):
            fringestr2 = ' , fintK2= ' + _np.str(fintSecVal) + ', fintxK2= ' + _np.str(fintxSecVal)
        elif (fintVal != 0) and (fintxSecVal == 0):
            fringestr2 = ' , fintK2= ' + _np.str(fintSecVal)
        elif (fintVal == 0) and (fintxSecVal != 0):
            fringestr2 = ' , fintxK2= ' + _np.str(fintxSecVal)
        else:
            fringestr2 = ''

        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Dipole ' + elementid + ', length= ' + _np.str(lenInM) + ' m, angle= ' + \
                      _np.str(_np.round(angle, 4)) + ' rad' + polefacestr + fringestr + fringestr2
        self.Writer.DebugPrintout('\t' + debugstring)

    def ChangeBend(self, linedict):
        """
        Function to change the direction of the dipole bend. Can be a direction other than horizontal (i.e != n*pi).
        """
        # NOT FULLY TESTED.
        angle = linedict['angle']
        rotation = False
        elementid = ''

        self.Transport.machineprops.angle = _np.float(angle)
        if self.Transport.machineprops.angle >= 360:
            self.Transport.machineprops.angle = _np.mod(self.Transport.machineprops.angle, 360)
        if self.Transport.machineprops.angle <= -360:
            self.Transport.machineprops.angle = _np.mod(self.Transport.machineprops.angle, -360)

        if self.Transport.machineprops.angle == 180 or self.Transport.machineprops.angle == -180:  # If 180 degrees, switch bending angle
            self.Transport.machineprops.bending *= -1

        elif self.Transport.machineprops.angle != 0:  # If not 180 degrees, use transform3d.
            # self.machineprops.angle *= -1
            # For conversion to correct direction. Eg in TRANSPORT -90 is upwards, in BDSIM, 90 is upwards.
            anginrad = self.Transport.machineprops.angle * (_np.pi / 180)
            self.Transport.machineprops.transforms += 1
            elementid = self._GetTransportElementName(linedict['name'])
            if not elementid:  # check on empty string
                elementid = 't' + _np.str(self.Transport.machineprops.transforms)

            # only call for gmad, warning for madx
            if self.Transport.convprops.gmadoutput:
                self.Transport.machine.AddTransform3D(name=elementid, psi=anginrad)
            elif self.Transport.convprops.madxoutput:
                self.Writer.DebugPrintout('\tWarning, MadX Builder does not have Transform 3D!')

            rotation = True

        if rotation:
            self.Writer.DebugPrintout('\tConverted to:')
            debugstring = '\tTransform3D ' + elementid + ', angle ' + _np.str(_np.round(self.Transport.machineprops.angle, 4)) + ' rad'
            self.Writer.DebugPrintout('\t'+debugstring)
        elif self.Transport.machineprops.angle == 180:
            self.Writer.DebugPrintout('\tBending direction set to Right')
        elif self.Transport.machineprops.angle == -180:
            self.Writer.DebugPrintout('\tBending direction set to Left')

    def Quadrupole(self, linedict):
        quaddata = linedict['data']
        length = quaddata[0]        # First three non-blanks must be the entries in a specific order.
        field_at_tip = quaddata[1]  # Field in TRANSPORT units
        pipe_rad = quaddata[2]      # Pipe Radius In TRANSPORT units

        field_in_Gauss = field_at_tip * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  # Convert to Tesla

        pipe_in_metres = pipe_rad * _General.ScaleToMeters(self.Transport, 'bend_vert_gap')
        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        field_gradient = (field_in_Tesla / pipe_in_metres) / self.Transport.beamprops.brho  # K1 in correct units

        self.Transport.machineprops.quads += 1

        elementid = self._GetTransportElementName(linedict['name'])
        if not elementid:  # check on empty string
            if field_gradient > 0:
                elementid = 'QF' + _np.str(self.Transport.machineprops.quads)
            elif field_gradient < 0:
                elementid = 'QD' + _np.str(self.Transport.machineprops.quads)
            else:
                elementid = 'NULLQUAD' + _np.str(self.Transport.machineprops.quads)  # For K1 = 0.

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddQuadrupole(name=elementid, length=lenInM, k1=_np.round(field_gradient, 4))

        string1 = '\tQuadrupole, field in gauss = ' + _np.str(field_in_Gauss) + ' G, field in Tesla = ' + _np.str(field_in_Tesla) + ' T.'
        string2 = '\tBeampipe radius = ' + _np.str(pipe_in_metres) + ' m. Field gradient = '+ _np.str(field_in_Tesla/pipe_in_metres) + ' T/m.'
        string3 = '\tBrho = ' + _np.str(_np.round(self.Transport.beamprops.brho, 4)) + ' Tm. K1 = ' +_np.str(_np.round(field_gradient, 4)) + ' m^-2'
        self.Writer.DebugPrintout(string1)
        self.Writer.DebugPrintout(string2)
        self.Writer.DebugPrintout(string3)
        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Quadrupole ' + elementid + ', length= ' + _np.str(lenInM) + ' m, k1= ' + _np.str(_np.round(field_gradient, 4)) + ' T/m'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Collimator(self, linedict):
        """
        A Function that writes the properties of a collimator element
        Only added for gmad, not for madx!
        """
        if linedict['length'] <= 0:
            self.Writer.DebugPrintout('\tZero or negative length element, ignoring.')
            return
        colldata = linedict['data']

        # Determine which entry is for horiz. and vert.
        aperx = self.Transport.machineprops.beampiperadius
        apery = self.Transport.machineprops.beampiperadius
        if _np.float(colldata[0]) == 1.0:
            aperx = colldata[1]
        elif _np.float(colldata[0]) == 3.0:
            apery = colldata[1]

        if len(colldata) > 2:
            if _np.float(colldata[2]) == 1.0:
                aperx = colldata[3]
            elif _np.float(colldata[2]) == 3.0:
                apery = colldata[3]
        aperx = _np.float(aperx)
        apery = _np.float(apery)

        lenInM = linedict['length'] * _General.ScaleToMeters(self.Transport, 'element_length')
        aperx_in_metres = aperx * _General.ScaleToMeters(self.Transport, 'x')
        apery_in_metres = apery * _General.ScaleToMeters(self.Transport, 'y')

        self.Transport.machineprops.collimators += 1
        elementid = self._GetTransportElementName(linedict['name'])
        if not elementid:  # check on empty string
            elementid = 'COL'+_np.str(self.Transport.machineprops.collimators)

        collimatorMaterial = 'copper'  # Default in BDSIM, added to prevent warnings
        # only call for gmad, warning for madx
        if self.Transport.convprops.gmadoutput:
            self.Transport.machine.AddRCol(name=elementid, length=lenInM, xsize=aperx_in_metres, ysize=apery_in_metres,
                                           material=collimatorMaterial)
        elif self.Transport.convprops.madxoutput:
            self.Writer.DebugPrintout('\tWarning, MadX Builder does not have RCOL')

        debugstring = '\tCollimator, x aperture = ' + _np.str(aperx_in_metres) \
                      + ' m, y aperture = ' + _np.str(apery_in_metres) + ' m.'
        self.Writer.DebugPrintout(debugstring)
        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Collimator ' + elementid + ', length= ' + _np.str(lenInM)\
                      + ' m, xsize= ' + _np.str(_np.round(aperx_in_metres, 4))
        debugstring += ' m, ysize= ' + _np.str(_np.round(apery_in_metres, 4)) + ' m.'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Acceleration(self, linedict):
        """
        A Function that writes the properties of an acceleration element
        Only RF added for gmad, not for madx!
        """
        # Redundant function until comments and /or acceleration components can be handled

        accdata = linedict['data']
        acclen = linedict['length']
        e_gain = linedict['voltage']

        # TODO add phase_lag and wavelength to BDSIM

        # If zero length then start of a sequence, save total accelerating voltage
        if acclen == 0.0:
            self.Transport.convprops.isAccSequence = True
            self.Transport.machineprops._totalAccVoltage = e_gain
            self.Transport.machineprops._e_gain_prev = 0.0  # start at 0
            return

        if self.Transport.convprops.isAccSequence:
            # voltage means voltage relative to the end of this segment
            e_rel_gain = e_gain - self.Transport.machineprops._e_gain_prev
            self.Transport.machineprops._e_gain_prev = e_gain  # store for next segment
            if e_gain == 1.0:  # end of sequence
                self.Transport.convprops.isAccSequence = False

            e_gain = e_rel_gain * self.Transport.machineprops._totalAccVoltage

        gradient = e_gain * (self.Transport.scale[self.Transport.units['p_egain'][0]] / 1e6)
        gradient /= (acclen * self.Transport.scale[self.Transport.units['element_length'][0]])  # gradient in MV/m

        self.Transport.machineprops.rf += 1
        elname = "ACC" + _np.str(self.Transport.machineprops.rf)

        # only call for gmad, warning for madx
        if self.Transport.convprops.gmadoutput:
            self.Transport.machine.AddRFCavity(name=elname, length=acclen, gradient=gradient)
        elif self.Transport.convprops.madxoutput:
            self.Writer.DebugPrintout('\tWarning, MadX Builder does not have RF Cavity')

        # Update beam parameters
        self.Transport = _General.UpdateMomentumFromEnergy(self.Transport, self.Transport.beamprops.k_energy + e_gain)

        # Commented out untested code
        # if len(accdata) == 2:  # Newer case with multiple elements
        # self._acc_sequence(line)
        if len(accdata) == 4:     # Older case for single element
            phase_lag = linedict['phase_lag']
            wavel = linedict['wavel']

            # Write to file
            accline = '! An accelerator element goes here of length ' + \
                      _np.str(acclen) + ' ' + self.Transport.units['element_length'] + ', \n'
            accline2 = '! with an energy gain of ' + _np.str(e_gain) + ' ' + \
                       self.Transport.units['p_egain'] + ', phase lag of ' + _np.str(phase_lag) + ' degrees, \n'
            accline3 = '! and a wavelength of ' + _np.str(wavel) + ' ' + self.Transport.units['bunch_length'] + '. \n'

    def Sextupole(self, linedict):
        sextudata = linedict['data']
        length = sextudata[0]        # First three non-blanks must be the entries in a specific order.
        field_at_tip = sextudata[1]  # Field in TRANSPORT units
        pipe_rad = sextudata[2]      # Pipe Radius In TRANSPORT units

        field_in_Gauss = field_at_tip * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  # Convert to Tesla

        pipe_in_metres = pipe_rad * _General.ScaleToMeters(self.Transport, 'bend_vert_gap')
        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        field_gradient = (2*field_in_Tesla / pipe_in_metres**2) / self.Transport.beamprops.brho  # K2 in correct units

        self.Transport.machineprops.sextus += 1
        elementid = self._GetTransportElementName(linedict['name'])
        if not elementid:  # check on empty string
            elementid = 'SEXT'+_np.str(self.Transport.machineprops.sextus)

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddSextupole(name=elementid, length=lenInM, k2=_np.round(field_gradient, 4))

        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Sextupole ' + elementid + ', length ' + _np.str(lenInM) + \
                      ' m, k2 ' + _np.str(_np.round(field_gradient, 4)) + ' T/m^2'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Solenoid(self, linedict):
        soledata = linedict['data']
        length = soledata[0]  # First three non-blanks must be the entries in a specific order.
        field = soledata[1]   # Field in TRANSPORT units

        field_in_Gauss = field * self.Transport.scale[self.Transport.units['magnetic_fields'][0]]  # Convert to Gauss
        field_in_Tesla = field_in_Gauss * 1e-4  # Convert to Tesla

        lenInM = length * _General.ScaleToMeters(self.Transport, 'element_length')

        self.Transport.machineprops.solenoids += 1
        elementid = self._GetTransportElementName(linedict['name'])
        if not elementid:  # check on empty string
            elementid = 'SOLE'+_np.str(self.Transport.machineprops.solenoids)

        # pybdsim and pymadx are the same.
        self.Transport.machine.AddSolenoid(name=elementid, length=lenInM, ks=_np.round(field_in_Tesla, 4))

        self.Writer.DebugPrintout('\tConverted to:')
        debugstring = 'Solenoid ' + elementid + ', length ' + _np.str(lenInM) + \
                      ' m, ks ' + _np.str(_np.round(field_in_Tesla, 4)) + ' T'
        self.Writer.DebugPrintout('\t' + debugstring)

    def Printline(self, linedict):
        number = linedict['data'][0]
        self.Writer.DebugPrintout('\tTRANSPORT control line,')
        try:
            number = _np.float(number)
            if number == 48:
                self.Transport.machineprops.benddef = False
                self.Writer.DebugPrintout('\t48. Switched Dipoles to Angle definition.')
            elif number == 47:
                self.Transport.machineprops.benddef = True
                self.Writer.DebugPrintout('\t47. Switched Dipoles to field definition.')
            elif number == 19:
                if _General.CheckSingleLineOutputApplied(self.Transport.convprops.filename):
                    self.Transport.convprops.singleLineOptics = True
                self.Writer.DebugPrintout('\t19. Optics output switched to single line per element.')
            else:
                self.Writer.DebugPrintout('\tCode 13. ' + _np.str(number) + ' handling not implemented.')
        except ValueError:
            pass

    def Correction(self, linedict):
        if self.Transport.convprops.correctedbeamdef:
            self.Writer.DebugPrintout('\tNot Correction to original beam definition')
            return
        # Check if the previous line was the original beam definition and not an rms update
        if linedict['prevlinenum'] == 1.0 and not linedict['isAddition'] and self.Transport.convprops.beamdefined:
            self.Transport.convprops.correctedbeamdef = True

        correctiondata = linedict['data']
        if len(correctiondata) >= 15:  # 15 sigma elements
            sigma21 = correctiondata[0]
            sigma43 = correctiondata[5]
        else:
            self.Writer.DebugPrintout('\tLength of correction line is incorrect')
            return

        emittoverbeta = self.Transport.beamprops.SigmaXP**2 * (1 - sigma21**2)
        emittbeta = self.Transport.beamprops.SigmaX**2
        betx = _np.sqrt(emittbeta / emittoverbeta)
        emitx = emittbeta / betx
        slope = sigma21 * self.Transport.beamprops.SigmaXP / self.Transport.beamprops.SigmaX
        alfx = -1.0 * slope * betx

        self.Transport.beamprops.betx = betx
        self.Transport.beamprops.emitx = emitx / 1000.0
        self.Transport.beamprops.alfx = alfx

        emittoverbeta = self.Transport.beamprops.SigmaYP**2 * (1 - sigma43**2)
        emittbeta = self.Transport.beamprops.SigmaY**2
        bety = _np.sqrt(emittbeta / emittoverbeta)
        emity = emittbeta / bety
        slope = sigma43 * self.Transport.beamprops.SigmaYP / self.Transport.beamprops.SigmaY
        alfy = -1.0 * slope * bety

        self.Transport.beamprops.bety = bety
        self.Transport.beamprops.emity = emity / 1000.0
        self.Transport.beamprops.alfy = alfy

        self.Transport.beamprops.distrType = 'gausstwiss'

        self.Writer.DebugPrintout('\tConverted to:')
        self.Writer.DebugPrintout('\t Beam Correction. Sigma21 = ' + _np.str(sigma21) + ', Sigma43 = ' + _np.str(sigma43) + '.')
        self.Writer.DebugPrintout('\t Beam distribution type now switched to "gausstwiss":')
        self.Writer.DebugPrintout('\t Twiss Params:')
        self.Writer.DebugPrintout('\t BetaX = ' + _np.str(self.Transport.beamprops.betx) + ' ' + self.Transport.units['beta_func'])
        self.Writer.DebugPrintout('\t BetaY = ' + _np.str(self.Transport.beamprops.bety) + ' ' + self.Transport.units['beta_func'])
        self.Writer.DebugPrintout('\t AlphaX = ' + _np.str(self.Transport.beamprops.alfx))
        self.Writer.DebugPrintout('\t AlphaY = ' + _np.str(self.Transport.beamprops.alfy))
        self.Writer.DebugPrintout('\t Emittx = ' + _np.str(self.Transport.beamprops.emitx) + ' ' + self.Transport.units['emittance'])
        self.Writer.DebugPrintout('\t EmittY = ' + _np.str(self.Transport.beamprops.emity) + ' ' + self.Transport.units['emittance'])

    def SpecialInput(self, linedict):
        specialdata = linedict['data']
        self.Writer.DebugPrintout('\tSpecial Input line:')

        if specialdata[0] == 5.0:  # beampiperadius (technically only vertical, but will apply a circle for now)
            self.Writer.DebugPrintout('\tType 5: Vertical half aperture,')
            self.Transport.machineprops.dipoleVertAper = specialdata[1]
            self.Writer.DebugPrintout('\tHalf aperture set to ' + _np.str(specialdata[1]) + '.')
            if self.Transport.machineprops.fringeIntegral == 0:
                self.Transport.machineprops.fringeIntegral = 0.5  # default if a vertical aperture is specified.
                self.Writer.DebugPrintout('Fringe field integral not set, setting to default of 0.5.')
        elif specialdata[0] == 7.0:  # Fringe Field integral
            self.Transport.machineprops.fringeIntegral = specialdata[1]
            self.Writer.DebugPrintout('\tType 7: K1 Fringe field integral,')
            self.Writer.DebugPrintout('\tIntegral set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 8.0:  # Second Fringe Field integral
            self.Transport.machineprops.secondfringeInt = specialdata[1]
            self.Writer.DebugPrintout('\tType 7: K2 Fringe field integral,')
            self.Writer.DebugPrintout('\tIntegral set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 14.0:  # Definition of element type code 6.
            if self.Transport.convprops.typeCode6IsTransUpdate:
                self.Transport.convprops.typeCode6IsTransUpdate = False
                typeCode6def = 'Collimator'
            else:
                self.Transport.convprops.typeCode6IsTransUpdate = True
                typeCode6def = 'Transform Update'
            self.Writer.DebugPrintout('\tType 14: Type code 6 definition,')
            self.Writer.DebugPrintout('\tDefinition set to ' + typeCode6def + '.')
        elif specialdata[0] == 12.0:  # entrance poleface curvature
            self.Transport.machineprops.bendInCurvature = specialdata[1]
            self.Writer.DebugPrintout('\tType 12: Dipole entrance poleface curvature,')
            self.Writer.DebugPrintout('\tCurvature set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 13.0:  # exit poleface curvature
            self.Transport.machineprops.bendOutCurvature = specialdata[1]
            self.Writer.DebugPrintout('\tType 13: Dipole exit poleface curvature,')
            self.Writer.DebugPrintout('\tCurvature set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 16.0:  # X0 offset
            self.Transport.beamprops.X0 = specialdata[1]
            self.Writer.DebugPrintout('\tType 16: X0 beam offset,')
            self.Writer.DebugPrintout('\tOffset set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 17.0:  # Y0 offset
            self.Transport.beamprops.Y0 = specialdata[1]
            self.Writer.DebugPrintout('\tType 17: Y0 beam offset,')
            self.Writer.DebugPrintout('\tOffset set to ' + _np.str(specialdata[1]) + '.')
        elif specialdata[0] == 18.0:  # Z0 offset
            self.Transport.beamprops.Z0 = specialdata[1]
            self.Writer.DebugPrintout('\tType 18: Z0 beam offset,')
            self.Writer.DebugPrintout('\tOffset set to ' + _np.str(specialdata[1]) + '.')
        else:
            self.Writer.DebugPrintout('\tCode type not yet supported, or unknown code type.')

    def UnitChange(self, linedict):
        """
        Function to change the units (scaling) of various parameters.
        """
        label = linedict['label']
        number = linedict['number']
        if label == 'CM' or label == 'MM' or label == 'UM' or label == 'NM':
            label = label.lower()
        # Convert Energy Unit Cases:
        if label == 'EV':
            label = 'eV'
        if label == 'KEV':
            label = 'keV'
        if label == 'MEV':
            label = 'MeV'
        if label == 'GEV':
            label = 'GeV'
        if label == 'TEV':
            label = 'TeV'

        debugstring2 = '\tConverted to ' + label

        if _np.float(number) == 1:  # Horizontal and vertical beam size
            self.Transport.units['x'] = label
            self.Transport.units['y'] = label
            self.Transport.units['bend_vert_gap'] = label
            # self.units['pipe_rad'] = label
            debugstring1 = '\tType 1: Horizontal and vertical beam extents, and magnet apertures,'

        elif _np.float(number) == 2:  # Horizontal and vertical divergence
            self.Transport.units['xp'] = label
            self.Transport.units['yp'] = label
            debugstring1 = '\tType 2: Horizontal and vertical angles,'

        elif _np.float(number) == 3:  # Bending Magnet Gap
            self.Transport.units['y'] = label
            self.Transport.units['bend_vert_gap'] = label
            debugstring1 = '\tType 3: Vertical (only) beam extent and magnet aperture,'

        elif _np.float(number) == 4:  # Vertical Divergence ONLY
            self.Transport.units['yp'] = label
            debugstring1 = '\tType 4: Vertical (only) beam angle,'

        elif _np.float(number) == 5:  # Pulsed Beam Length
            self.Transport.units['bunch_length'] = label
            debugstring1 = '\tType 5: Bunch length,'

        elif _np.float(number) == 6:  # Momentum Spread
            self.Transport.units['momentum_spread'] = label  # Percent
            debugstring1 = '\tType 6: Momentum spread,'

        elif _np.float(number) == 7:  # Bend/pole face rotation
            debugstring1 = '\tType 7: Bend and poleface rotation angles,'
            debugstring2 = '\tCONVERTION NOT IMPLEMENTED YET.'
            pass

        elif _np.float(number) == 8:  # Element Length
            self.Transport.units['element_length'] = label
            debugstring1 = '\tType 8: Element length,'

        elif _np.float(number) == 9:  # Magnetic Field
            self.Transport.units['magnetic_fields'] = label
            debugstring1 = '\tType 9: Magnetic Fields,'

        elif _np.float(number) == 10:  # Mass
            debugstring1 = '\tType 10: Mass,'
            debugstring2 = '\tCONVERTION NOT IMPLEMENTED YET.'
            pass

        elif _np.float(number) == 11:  # Momentum / energy gain during acc.
            self.Transport.units['p_egain'] = label
            debugstring1 = '\tType 11: Momentum and accelerator energy gain,'
        else:
            # default output
            debugstring1 = '\tCode type not yet supported, or unknown code type.'
            debugstring2 = ''

        self.Writer.DebugPrintout('\tUnit change line:')
        self.Writer.DebugPrintout(debugstring1)
        self.Writer.DebugPrintout(debugstring2)

    def TransformUpdate(self, linedict):
        if linedict['elementnum'] == 6.0:
            errorline = '\tElement is either a transform update or a collimator. The type code 6 definition'
            errorline2 = '\thas not been switched to collimators, therefore nothing will be done for this element.'
            self.Writer.DebugPrintout(errorline)
            self.Writer.DebugPrintout(errorline2)

    def _GetTransportElementName(self, elementName):
        """
        Checks if name already used by an element in the machine.
        Appends a _N to the name where N is the lowest integer not already
        used in a elementName_N name in the machine.
        """
        if not self.Transport.convprops.keepName:
            return ''

        currElements = self.Transport.machine.elementsd.keys()
        updatedName = elementName
        if elementName in currElements:
            nameNumber = 1
            while True:
                updatedName = elementName + "_" + _np.str(nameNumber)
                if updatedName not in currElements:
                    break
                else:
                    nameNumber += 1

        return updatedName

    def _UpdateElementsFromFits(self):
        # Functions that update the elements in the element registry.
        # For debugging purposes, they return dictionaries of the element type,
        # length change details, and which parameters were updated and the values in a list which
        # follows the pattern of [parameter name (e.g. 'field'),oldvalue,newvalue]

        # Length update common to nearly all elements, seperate function to prevent duplication
        # fitIndex was used in the past. Pass in anyway in case of future need, just delete for now.
        def _updateLength(eleIndex, fitIndex, element):
            del fitIndex
            oldlength = self.Transport.ElementRegistry.elements[eleIndex]['length']
            lengthDiff = self.Transport.ElementRegistry.elements[eleIndex]['length'] - element['length']
            self.Transport.ElementRegistry.elements[eleIndex]['length'] = element['length']  # Update length
            self.Transport.ElementRegistry.length[eleIndex:] += lengthDiff  # Update running length of subsequent elements.
            self.Transport.ElementRegistry._totalLength += lengthDiff  # Update total length
            lendict = {'old': _np.round(oldlength, 5),
                       'new': _np.round(element['length'], 5)}
            return lendict

        def _updateDrift(eleIndex, fitIndex, element):
            eleDict = {'updated': False,
                       'element': 'Drift',
                       'params': []}

            # Only length can be varied
            if self.Transport.ElementRegistry.elements[eleIndex]['length'] != element['length']:
                lendict = _updateLength(eleIndex, fitIndex, element)
                eleDict['updated'] = True
                eleDict['length'] = lendict
            return eleDict

        def _updateQuad(eleIndex, fitIndex, element):
            eleDict = {'updated': False,
                       'element': 'Quadrupole',
                       'params': []}

            if self.Transport.ElementRegistry.elements[eleIndex]['data'][1] != element['data'][1]:  # Field
                oldvalue = self.Transport.ElementRegistry.elements[eleIndex]['data'][1]
                self.Transport.ElementRegistry.elements[eleIndex]['data'][1] = element['data'][1]
                eleDict['updated'] = True
                data = ['field', oldvalue, element['data'][1]]
                eleDict['params'].append(data)


            if self.Transport.ElementRegistry.elements[eleIndex]['length'] != element['length']:
                self.Transport.ElementRegistry.elements[eleIndex]['data'][0] = element['data'][0]  # Length in data
                lendict = _updateLength(eleIndex, fitIndex, element)
                eleDict['length'] = lendict
                eleDict['updated'] = True

            return eleDict

        def _updateDipole(eleIndex, fitIndex, element):
            eleDict = {'updated': False,
                       'element': 'Dipole',
                       'params': []}

            # TODO: Need code in here to handle variation in poleface rotation. Not urgent for now.
            if self.Transport.ElementRegistry.elements[eleIndex]['data'][1] != element['data'][1]:  # Field
                oldvalue = self.Transport.ElementRegistry.elements[eleIndex]['data'][1]
                if not self.Transport.machineprops.benddef:  # Transport can switch dipole input definition
                    par = 'field'
                    self.Transport.ElementRegistry.elements[eleIndex]['data'][1] = element['data'][1]
                    eleDict['updated'] = True
                else:
                    par = 'angle'
                    self.Transport.ElementRegistry.elements[eleIndex]['data'][1] = element['data'][3]
                    eleDict['updated'] = True
                data = [par, oldvalue, element['data'][3]]
                eleDict['params'].append(data)
            if self.Transport.ElementRegistry.elements[eleIndex]['length'] != element['length']:
                self.Transport.ElementRegistry.elements[eleIndex]['data'][0] = element['data'][0]  # Length in data
                lendict = _updateLength(eleIndex, fitIndex, element)
                eleDict['updated'] = True
                eleDict['length'] = lendict
            return eleDict

        self.Writer.DebugPrintout("")
        self.Writer.DebugPrintout("Updating elements with results from fitting")
        for index, name in enumerate(self.Transport.FitRegistry._uniquenames):
            fitindex = self.Transport.FitRegistry.GetElementIndex(name)
            eleindex = self.Transport.ElementRegistry.GetElementIndex(name)

            for fitnum, fit in enumerate(fitindex):
                for elenum, ele in enumerate(eleindex):
                    fitelement = self.Transport.FitRegistry.elements[fitindex[fitnum]]
                    if fitelement['elementnum'] == 3:
                        eledict = _updateDrift(eleindex[elenum], fitindex[fitnum], fitelement)
                    elif fitelement['elementnum'] == 4:
                        eledict = _updateDipole(eleindex[elenum], fitindex[fitnum], fitelement)
                    elif fitelement['elementnum'] == 5:
                        eledict = _updateQuad(eleindex[elenum], fitindex[fitnum], fitelement)

                    if eledict['updated']:
                        self.Writer.DebugPrintout(
                            "Element " + _np.str(eleindex[elenum]) + " was updated from fitting:")
                        self.Writer.DebugPrintout("\tOptics Output line:")
                        self.Writer.DebugPrintout(
                            "\t\t'" + self.Transport.FitRegistry.lines[fitindex[fitnum]] + "'")
                        if eledict.has_key('length'):
                            lenline = "\t" + eledict['element'] + " length updated to "
                            lenline += _np.str(eledict['length']['new'])
                            lenline += " (from " + _np.str(eledict['length']['old']) + ")."
                            self.Writer.DebugPrintout(lenline)
                        for param in eledict['params']:
                            parline = "\t" + eledict['element'] + " " + param[0]
                            parline += " updated to " + _np.str(param[2]) + " (from " + _np.str(param[1]) + ")."
                            self.Writer.DebugPrintout(parline)
                            # self.Writer.DebugPrintout("\n")

                        break
