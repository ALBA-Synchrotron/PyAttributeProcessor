PyAttributeProcessor by srubio@cells.es,cpascual@cells.es

This class inherits from fandango.DynamicDS, it is based on PySignalSimulator but providing extra functionality 
like FFT, GaussPeak, PeakFit, all numpy and scipy methods and properties to load extra libraries on demand.

https://github.com/tango-controls/fandango/blob/documentation/doc/recipes/DynamicDS_and_Simulators.rst

Modules on ExtraModules property are just added like::

  PyTangoArchiving
  PyTangoArchiving.Reader as HDB

So you can use them in DynamicAttributes declaration::

  MeanPeriod = DevDouble(VAR('MeanPeriod',VALUE,default=3600)) # An R/W attribute
  LastValues = DevVarDoubleArray(t[1] for t in HDB().get_attribute_values("my/tango/attribute/name",- MeanPeriod))
  MeanValue = DevDouble(numpy.average(LastValues))



