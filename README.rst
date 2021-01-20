PyAttributeProcessor by srubio@cells.es,cpascual@cells.es

Development repository for this project is:

https://github.com/alba-synchrotron/pyattributeprocessor

This class inherits from fandango.DynamicDS, but providing extra libraries and functions:

https://github.com/tango-controls/fandango/blob/documentation/doc/recipes/DynamicDS_and_Simulators.rst

Extra functions/modules available are:

 - FFT
 - GaussPeak 
 - PeakFit 
 - math
 - time
 - random
 - numpy.*
 - scipy.*
 - scipy.signal.*
 
Other Modules can be loaded using ExtraModules property::

  PyTangoArchiving
  PyTangoArchiving.Reader as HDB

So you can use them in DynamicAttributes declaration::

  MeanPeriod = DevDouble(VAR('MeanPeriod',VALUE,default=3600)) # An R/W attribute
  LastValues = DevVarDoubleArray(t[1] for t in HDB().get_attribute_values("my/tango/attribute/name",- MeanPeriod))
  MeanValue = DevDouble(numpy.average(LastValues))



