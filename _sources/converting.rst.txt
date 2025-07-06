==========
Conversion
==========

pytransport can convert a TRANSPORT "FOR001" file to BDSIM's `gmad` format::

  >>> from pytransport import Convert
  >>> import pybdsim.Builder
  >>> file = 'FOR001.DAT'
  >>> Convert.Convert(file, machine=pybdsim.Builder.Machine(), options=pybdsim.Options.Options())

  Writing to file: bdsim/FOR001.gmad
  Lattice written to:
  FOR001_components.gmad
  FOR001_sequence.gmad
  FOR001_beam.gmad
  FOR001_options.gmad
  All included in main file:
  FOR001.gmad
