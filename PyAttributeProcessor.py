#!/usr/bin/env python

#############################################################################
##
## file :    header.py
##
## description : Device Server that processes attributes.
##               Python source for the PyAttributeProcessor and its commands. 
##               The class is derived from Device. It represents the
##               CORBA servant object which will be accessed from the
##               network. All commands which can be executed on the
##               PyAttributeProcessor are implemented in this file.
##               Based on the code of PySignalSimulator by srubio
##
## project :    Tango-ds
##
## developers history: cpascual@cells.es
##                     (based on  PySignalSimulator by srubio)
##
## $Revision:  $
##
## $Log:  $
##
## copyleft :    Cells / Alba Synchrotron
##               Bellaterra
##               Spain
##
#############################################################################
##
## This file is part of Tango-ds.
##
## This is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This software is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, see <http://www.gnu.org/licenses/>.
###########################################################################


import sys
from scipy import *
from scipy.signal import *
from scipy.optimize import leastsq

import PyTango
from fandango.dynamic import DynamicDS,DynamicDSClass,DynamicAttribute
from fandango.interface import FullTangoInheritance


def FFT(a, power=True, phase=True, **kwargs):
    """returns full output of the fft of a by concatenating the power
    and the phase vectors for the calculated fft.
    if power==True, it returns the power spectrum
    if phase==True it returns the phase spectrum
    If both are true, it concatenates them.  @TODO: return a 2D attribute instead
    See scipy.fft for the **kwargs
    """    
    #When using scipy.fft(), each element of the output corresponds
    #to a frequency that is built as follows:
    #freq=[0,1,2,3,...,n/2,-n/2+1, -n/2+2,...,-2,-1] 
    #where n=a.size and the unit is [1/bin]
    #
    #For FFTs of real numbers (as is our case), the following is true:
    #   fft(x_i)=conj(fft(x_(N-i))
    #So the power spectrum is symetric (and the phase antisymetric) with respect 
    #to the 0 frequency, and therefore only the first half of the solution
    #from scipy.fft() is needed (we discard the higest frequency, the n/2 term).
    
    ##@TODO: include a vector with frequencies in the output (we'll do it when Image support is available)
    a=array(a)
    n=kwargs.get('n',a.size)
    result=fft(a,**kwargs)[:n/2]
    if power and phase: return concatenate((abs(result), arctan2(result.imag,result.real)))
    elif power: return abs(result)
    elif phase: return arctan2(result.imag,result.real)
    
def GaussPeak(x, *args):
    """Generates a gaussian peak as a function of its mean, FWHM and the
    maximum heigh. Optionally also an ordinate offset"""
    mean, ymax, fwhm =args
    if len(args)>3: yoffset=args[3]
    else: yoffset=0 
    return yoffset + ymax*power(2,-4*((x-mean)/fwhm)**2)

def PEAKFIT(yexp, x=None, roimin=None, roimax=None, fitfunc=None, Chi2Warning=2, **fitparameters):
    """It fits experimental data within a ROI to a peak (for the moment, to a gaussian)"""
    if roimin is None: roimin=0 ##@TODO:read it from an attribute before defaulting
    if roimax is None: roimax=yexp.size ##@TODO:read it from an attribute before defaulting
    yexp=yexp[roimin:roimax]
    if x is None: x=arange(roimin,roimax) ##@TODO:read it from an attribute before defaulting
    if fitfunc is None: fitfunc=GaussPeak

    if fitfunc is GaussPeak:
        ymax=fitparameters.get('ymax',yexp.max()) ##@TODO:read it from an attribute before defaulting
        mean=fitparameters.get('mean',x[argmax(yexp)]) ##@TODO:read it from an attribute before defaulting
        fwhm=fitparameters.get('fwhm',4) ##@TODO:read it from an attribute before defaulting
        yoffset=fitparameters.get('yoffset',None) ##@TODO:read it from an attribute before defaulting
        if yoffset is None:  params=(mean, ymax, fwhm)
        else: params=(mean, ymax, fwhm, yoffset)
    else:
        raise NotImplementedError
    
    residuals = lambda p, x, yexp: fitfunc(x, *p)-yexp
    params,ok=leastsq(residuals, params, args=(x,yexp))
    
    #check quality of result
    if ok:
        imin=max(roimin, int(mean-fwhm))
        imax=min(roimax, int(mean+fwhm))
        peakres=residuals(params,x[imin:imax],yexp[imin:imax])
        normChi2=(peakres*peakres/yexp[imin:imax]).sum()/(peakres.size-params.size)
        if normChi2 < Chi2Warning:
            quality=PyTango.AttrQuality.ATTR_VALID
        else:
            quality=PyTango.AttrQuality.ATTR_WARNING
    else:
        print "Fit failed"
        quality=PyTango.AttrQuality.ATTR_INVALID
        normChi2=99999
        
    return DynamicAttribute(value=concatenate(([normChi2],params)),quality=quality)
    
class PyAttributeProcessor(PyTango.Device_4Impl):

#--------- Add you global variables here --------------------------

#------------------------------------------------------------------
#    Device constructor
#------------------------------------------------------------------
    def __init__(self,cl, name):
        #PyTango.Device_4Impl.__init__(self,cl,name)
        print 'IN PYATTRIBUTEPROCESSOR.__INIT__'
        _locals = {}
        _locals.update(locals())
        _locals.update(globals())
        #print '_locals are:\n%s' % _locals
        DynamicDS.__init__(self,cl,name,_locals=_locals,useDynStates=True)
        PyAttributeProcessor.init_device(self)

#------------------------------------------------------------------
#    Device destructor
#------------------------------------------------------------------
    def delete_device(self):
        print "[Device delete_device method] for device",self.get_name()


#------------------------------------------------------------------
#    Device initialization
#------------------------------------------------------------------
    def init_device(self):
        print "In ", self.get_name(), "::init_device()"
        self.set_state(PyTango.DevState.ON)
        self.get_device_properties(self.get_device_class())
        if self.DynamicStates: self.set_state(PyTango.DevState.UNKNOWN)

#------------------------------------------------------------------
#    Always excuted hook method
#------------------------------------------------------------------
    def always_executed_hook(self):
        print "In ", self.get_name(), "::always_excuted_hook()"
        DynamicDS.always_executed_hook(self)

#==================================================================
#
#    PyAttributeProcessor read/write attribute methods
#
#==================================================================
#------------------------------------------------------------------
#    Read Attribute Hardware
#------------------------------------------------------------------
    def read_attr_hardware(self,data):
        print "In ", self.get_name(), "::read_attr_hardware()"
    
            


#==================================================================
#
#    PyAttributeProcessor command methods
#
#==================================================================
    
#==================================================================
#
#    PyAttributeProcessorClass class definition
#
#==================================================================
class PyAttributeProcessorClass(PyTango.DeviceClass):

    #    Class Properties
    class_property_list = {
        }


    #    Device Properties
    device_property_list = {
        'DynamicAttributes':
            [PyTango.DevVarStringArray,
            "Attributes and formulas to create for this device.\n<br/>\nThis Tango Attributes will be generated dynamically using this syntax:\n<br/>\nT3=int(SomeCommand(7007)/10.)\n\n<br/>\nSee the class description to know how to make any method available in attributes declaration.",
            [ ] ],
        'DynamicStates':
            [PyTango.DevVarStringArray,
            "This property will allow to declare new States dinamically based on\n<br/>\ndynamic attributes changes. The function Attr will allow to use the\n<br/>\nvalue of attributes in formulas.<br/>\n\n\n<br/>\nALARM=Attr(T1)>70<br/>\nOK=1",
            [ ] ],
        }


    #    Command definitions
    cmd_list = {
        }


    #    Attribute definitions
    attr_list = {
        }


#------------------------------------------------------------------
#    PyAttributeProcessorClass Constructor
#------------------------------------------------------------------
    def __init__(self, name):
        PyTango.DeviceClass.__init__(self, name)
        self.set_type(name);
        print "In PyAttributeProcessorClass  constructor"

#==================================================================
#
#    PyAttributeProcessor class main method
#
#==================================================================
if __name__ == '__main__':
    try:
        py = PyTango.Util(sys.argv)
        # Adding TRUE DeviceServer Inheritance
        from fandango.interface import FullTangoInheritance
        PyAttributeProcessor,PyAttributeProcessorClass = FullTangoInheritance('PyAttributeProcessor',PyAttributeProcessor,PyAttributeProcessorClass,DynamicDS,DynamicDSClass,ForceDevImpl=True)
        py.add_TgClass(PyAttributeProcessorClass,PyAttributeProcessor,'PyAttributeProcessor')

        U = PyTango.Util.instance()
        U.server_init()
        U.server_run()

    except PyTango.DevFailed,e:
        print '-------> Received a DevFailed exception:',e
    except Exception,e:
        print '-------> An unforeseen exception occured....',e
