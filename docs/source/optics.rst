======
Optics
======

pytransport can read a TRANSPORT `FOR002` output file that has sigma matrices.

  >>> import pytransport
  >>> optics = pytransport.Reader.GetOptics('FOR002.DAT')
  >>> import matplotlib.pyplot as plt
  >>> plt.plot(optics.S(), optics.Sigma_X())

Also with `pybdsim` the TRANSPORT optics can be directly compared with BDSIM::

  >>> pybdsim.Compare.TransportVsBDSIM('FOR002.DAT', 'bdsim_optics.root')