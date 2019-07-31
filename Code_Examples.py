Correção NMO

""" OpenCPS_tool_v1
===begin_json_info===
{
    "name"       : "ApplyNmo",
    "purpose"    : "compute and apply NMO using a picked table",
    "categories" : ["Python"],
   
    "guiwidgets" : {
        "vel_tbl" : {
            "widget"     : "path",
            "label"      :"Velocity table:",
            "assettype"  : ["table"]
        }
    }
}
===end_json_info===
===begin_latex_help===
\thistool\ is an example showing how to use tables in python scripts. To apply normal moveout to data,
use \toolref{NormalMoveout}.
===end_latex_help===
"""
import math
import numpy

from opencps.tool import *

class Tool(basetool):
    def update(self):
        self.label = 'ApplyNMO'
        self.vels = spatialTableParam(self, "vel_tbl")
        self.readsHeaders("ILINE", "XLINE", "OFFSET")

    def execute(self, ntr, hds, trs):
        if ntr <= 0: return ntr

        ## figure out the inline and crossline coordinates
        il  = hds[0]["ILINE"]
        xl  = hds[0]["XLINE"]
        ## read in the velocity function
        vels = self.vels.getInterpolatedFunction(y=il, x=xl, m0=self.o1In, dm=self.d1In, nm=self.n1In)

        ## apply NMO to each trace in the gather
        for itr in xrange(0, ntr):
            self.applyNmo(vels, trs[itr], hds[itr]["OFFSET"])
           
        return ntr

    def applyNmo(self, vels, tr, off):
        t0 = self.o1In*.001+numpy.arange(self.n1In)*self.d1In*.001

        ## NMO equation
        tmap = (t0**2 + off**2 / vels**2)**.5

        ## convert to fractional sample indices
        map  = (tmap-self.o1In*.001) / self.d1In * 1000

        ## convert to integer sample indices and compute fractions
        imap = numpy.array(map,numpy.int32)
        frac = map - imap

        ## set out-of-range samples indices to 0
        imap[map<0] = 0; imap[map>=self.n1In-1] = 0

        ## linear interpolation
        tmp = tr[imap]*(1-frac)+tr[imap+1]*frac
 
        ## zero out out-of-range samples indices to 0
        tmp[map<0] = 0; tmp[map>=self.n1In-1] = 0

        ## set the result on the output trace
        tr[:self.n1In] = tmp[:]

____________________________________________________
FGain
""" OpenCPS_tool_v1
===begin_json_info===
{   "name"       : "FGain",
    "categories" : ["Python"],
    "purpose"    : "example for showing how to use FFTs in Python",
   
    "guiwidgets" : {
        "specopt" : {
            "widget"  : "radio",
            "label"   : "Specify:",
            "labels"  : ["frequency-amplitude pairs", "frequency-db pairs"],
            "options" : ["freqamp"                  , "freqdb"            ]            
        },
        "freqamp_pairs" : {
            "widget" : "string",
            "label"  : "freq-Amp pairs:",
            "indent" : 20,
            "width"  : 200,
            "tooltip" : "Specify frequency- Absolute amplitude pairs: Frequency1-Amplitude1 / Frequency2-Amplitude2 to spectrally shape the trace"
        },
        "freqdb_pairs" : {
            "widget" : "string",
            "label"  : "freq-DB pairs:",
            "indent" : 20,
            "width"  : 200,
            "tooltip": "Specify frequency- Db amplitude pairs: Frequency1 - Amplitude1 / Frequency2 - Amplitude2 to spectrally shape the trace"
        },
        "do_smooth" : {
            "widget"  : "checkbox",
            "label"   : "Smooth the filter",
            "tooltip" : "Enable to specify a smoothing filter that reduces/prevents ringing"            
        },
        "smooth_len" : {
            "widget" : "float",
            "label"  : "smoothing length (Hz):",
            "indent" : 20,
            "tooltip": "Specify a length of smoothing filter, in Hz, to apply to the filter before amplitudes are adjusted. This is used to reduce/prevent ringing"
        }
    }
}
===end_json_info===

===begin_latex_help===
\subsubsection*{License required}

This tool is available under the base license package.

\subsubsection*{MPI capability}

This tool supports limited MPI, the data will be split on the slowest rank and re-combined on Output. This tool does not support OpenMP multithreading.  

\subsubsection*{Description}

FGain allows the user to manipulate the power spectrum of the data. This is commonly referred to as spectral shaping.  
The shaping can be specified as frequency-amplitude pairs in either absolute amplitudes or dB log amplitudes.  A smoothing operator, specified in Hz, can be applied to the filter before application to the data.
Parameters follow a standard type-in format: Frequency1 - Amplitude1 / Frequency2 - Amplitude2 / Frequency3 - Amplitude3 / Frequency4 - Amplitude4
Note: Absolute amplitude pairs will multiply the amplitude of the frequency band linearly; i.e. an amplitude of 2 will double the resultant amplitudes. Db amplitudes follow the Decibel scale, so 6 dB will double the resultant amplitudes. The smoothing operator is specified in Hz and is the length of the averaging filter being applied.
\subsubsection*{Input}

This tool operates trace by trace so input order does not matter. Input data should normally be .seis formatted seismic traces.
\subsubsection*{Output}

Output data will be spectrally shaped by the FGain operation. Amplitudes will vary and frequency content will change. Seismic traces are output.

===end_latex_help===
"""
import math
import numpy as np
import re

import pythontool
from pythontool import *
import opencps_native_wrapper as wrapper

class Tool(basetool):
    def update(self):
        self.label = 'FGain'
        freqamp_pairs  = getStringParam("freqamp_pairs", self.locals, "0-0 / 10-1 / 40-1 / 50-0"   )
        freqdb_pairs   = getStringParam("freqdb_pairs" , self.locals, "30-0 / 40-10 / 50-10 / 60-0")
        self.specopt   = getStringParam("specopt"      , self.locals, "freqamp")
        self.dosmooth  = getBoolParam  ("do_smooth"    , self.locals, True)
        self.smoothlen = getFloatParam ("smooth_len"   , self.locals,   5.)
        setVisible(self.locals, "freqdb_pairs" , self.specopt=="freqdb" )
        setVisible(self.locals, "freqamp_pairs", self.specopt=="freqamp")
        setVisible(self.locals, "smooth_len"   , self.dosmooth)
       
        self.nfftr = fastFFTSize(self.n1In)
        self.df = 1/(self.d1In*.001*self.nfftr)
        self.nf = self.nfftr/2+1
       
        (self.fmult,err) = freqDbPairs2fmult(freqdb_pairs if self.specopt=="freqdb" else freqamp_pairs)
        if err: self.errors += [err]
           
    def startExecution(self):
        self.filt = interp(self.fmult, self.df, self.nf)
        if self.specopt=="freqdb":
            self.filt = np.power(10,.05*self.filt)
        if self.dosmooth:
            npts = clamp(int(.5+self.smoothlen/self.df),1,self.nfftr)
            self.filt = wrapper.mathutil_runningAverage(npts,self.filt)
            print npts
           
    def execute(self, ntr, hds, trs):
        if ntr <= 0: return ntr
       
        tr = np.zeros(self.nfftr)
        for itr in xrange(0,ntr):
            tr[:self.n1In] = trs[itr,:]
            ftr = numpy.fft.rfft(tr)
            ftr[:] *= self.filt
            trs[itr,:] = numpy.fft.irfft(ftr,self.nfftr)[:self.n1In]
        return ntr

def freqDbPairs2fmult(freqdb_pairs):
    try:
        map = {}
        for pair_ in freqdb_pairs.split("/"):
            pair = pair_.split("-",1)
            f = float(pair[0])
            v = float(pair[1])
            map[f] = v
        return (map,None)
    except:
        return (None,"sorry, freq-db pairs cannot be parsed")
   
def interp(map,df,nf):
    keys = np.array(map.keys())
    keys.sort()
    filt = np.zeros(nf,'float32')
    for ifreq in xrange(0,nf):
        f = df*ifreq
        r = np.searchsorted(keys,f)
        if   r==0        : vinterp = map[keys[r]]
        elif r>=len(keys): vinterp = map[keys[len(keys)-1]]
        else:
            lf = keys[r-1]; lv = map[lf]
            rf = keys[r  ]; rv = map[rf]
            frac = (f-lf)/(rf-lf)
            vinterp = lv*(1-frac)+rv*frac
        filt[ifreq] = vinterp
    return filt
 
def fastFFTSize(nt):
    while True:
        m = nt
        while m%2==0: m /= 2
        while m%3==0: m /= 3
        while m%5==0: m /= 5
        if m<=1: break
        nt+=1
    return nt


______________________________________________________-

fft

""" OpenCPS_tool_v1
===begin_json_info===
{
    "name"       : "FFT",
    "categories" : ["Python", "Transforms"],
    "purpose"    : "1D Fourier transform",
   
    "guiwidgets" : {
        "mode" : {
            "widget"  : "combo",
            "label"   : "Output mode",
            "options" : ["real+imag","real", "imag"]
        },
        "fast_nfftr" : {
            "widget"  : "checkbox",
            "label"   : "Use fast FFT size"
        },
        "__inv_lbl" : {
            "widget"  : "staticlabel",
            "label"   : "Input real+imag traces detected",
            "tooltip" : ["Input traces are assumed to be complex if ",
                         "the following four parameters are present ",
                         "in the input parameter list:\n",
                         "  fft_orig_d1\n",
                         "  fft_orig_n1\n",
                         "  fft_orig_l1\n",
                         "  fft_orig_nfftr"]
        }
    }                        
}
===end_json_info===
===begin_latex_help===
===end_latex_help===
"""
import math
import numpy as np
import numpy.fft

from opencps.tool import *
import opencps_native_wrapper as wrapper

class Tool(basetool):
    def update(self):
        self.inv = True
        for k in ['d1','n1','nfftr','l1']:
            if 'fft_orig_'+k not in self.inputs: self.inv = False
        setVisible(self.locals,'fast_nfftr', not self.inv)
        setVisible(self.locals,'mode'      , not self.inv)
        setVisible(self.locals,'__inv_lbl' ,     self.inv)
        if self.inv:
            self.label = 'iFFT'
            self.d1Out = self.inputs['fft_orig_d1'   ]
            self.n1Out = self.inputs['fft_orig_n1'   ]
            self.nfftr = self.inputs['fft_orig_nfftr']
            self.outputs["dims"][0] = self.inputs['fft_orig_l1']
            del self.inputs['fft_orig_d1'   ]
            del self.inputs['fft_orig_n1'   ]
            del self.inputs['fft_orig_nfftr']
            del self.inputs['fft_orig_l1'   ]
            self.nf = self.n1In/2
        else:
            self.label = 'FFT'
            self.isfast = getBoolParam  ("fast_nfftr", self.locals, False)
            self.mode   = getStringParam("mode"      , self.locals, "real+imag")
            self.nfftr = self.n1In
            if self.isfast:
                self.nfftr = np_fast_fft_size(self.nfftr)
            self.df = 1/(self.d1In*.001*self.nfftr)
            self.nf = self.nfftr/2+1
            self.d1Out = self.df
            self.outputs['dims'][0] = 'FREQ'
            if self.mode=="real+imag":
                self.n1Out = self.nf*2
                self.outputs['fft_orig_d1'   ] = self.d1In
                self.outputs['fft_orig_n1'   ] = self.n1In
                self.outputs['fft_orig_nfftr'] = self.nfftr
                self.outputs['fft_orig_l1'   ] = self.inputs["dims"][0]
            else:
                self.n1Out = self.nf
               
    def startExecution(self):
        self.tbuf = np.zeros(self.nfftr, 'float')
        self.fbuf = np.zeros(self.nf   ,'cfloat')
   
    def execute(self, ntr, hds, trs):
        for itr in xrange(0,ntr):
            tr = trs[itr,:]
            if self.inv:
                real = tr[       :self.nf]
                imag = tr[self.nf:       ]
                self.fbuf[:] = real+1j*imag
                tr_nfftr = np.fft.irfft(self.fbuf,self.nfftr)
                tr[:self.n1Out] = tr_nfftr[:self.n1Out]
            else:
                if self.nfftr!=self.n1In:
                    self.tbuf[:self.n1In] = tr[:self.n1In]
                    self.fbuf[:         ] = np.fft.rfft(self.tbuf)
                else:
                    self.fbuf[:] = np.fft.rfft(tr[:self.n1In])
                if self.n1Out==self.nf*2:
                    tr[       :self.nf] = np.real(self.fbuf)
                    tr[self.nf:       ] = np.imag(self.fbuf)
                elif self.mode=='real':
                    tr[       :self.nf] = np.real(self.fbuf)
                elif self.mode=='imag':
                    tr[       :self.nf] = np.imag(self.fbuf)
        return ntr


_____________________________________________________--
example

# OpenCPS_tool_v1
"""===begin_json_info===
{
    "name"       : "KitchenSinkPy",
    "purpose"    : "To demonstrate the Python API",
    "categories" : ["Examples"],
   
   
    "guiwidgets" : {
        "positive_number" : {
            "widget"  : "integer",
            "label"   : "Positive integer:",
            "tooltip" : "Enter a positive integer. If you enter a negative integer, the tool will show an error."
        },
        "header_to_add" : {
            "widget"  : "string",
            "label"   : "Add header:",
            "tooltip" : "Add a header by this name to the output header list."
        },
        "header_to_remove" : {
            "widget"  : "string",
            "label"   : "Remove header:",
            "tooltip" : ["Remove header by this name from the output header list. By default, all ",
                         "input headers are copied over to the output header list.\n",
                         "If no value is specified, then no header will be removed from the output ",
                         "header list."]
        },
        "header_to_require" : {
            "widget"  : "string",
            "label"   : "Require input header:",
            "tooltip" : ["Require this header in the input header list.\n",
                         "If the header does not exist in the input header list, an error will be shown."]
        },
        "do_extract_picks" : {
            "widget"  : "checkbox",
            "label"   : "Extract function from velocity table at center of grid"
        },
        "extraction_tbl" : {
            "widget"     : "path",
            "label"      :"picked (i.e. sparse, .tbl) table:",
            "assettype"  : ["table"],
            "indent"     : 20,
            "tooltip"    : "The picks for this table will be printed to the flow log."
        }
    }
}
===end_json_info==="""

from operator import itemgetter
from opencps.tool import *
   
class Tool(basetool):
    def update(self):
        self.label = "KitchenSinkPy"
       
        print("This is printed from the KitchenSinkPy tool's update() method")
       
        if getBoolParam("__initdefaults", self.locals, False):
            # default other params
            self.locals["header_to_add"    ] = "FOO"
            self.locals["header_to_remove" ] = "CMP"
            self.locals["header_to_require"] = "FFID"            
            self.locals["positive_number"  ] = 1
            self.locals["remove_dead"      ] = True        
       
        # Validate user-specified parameters in update().
        # A flow is not allowed to execute if any tool errors exist.
        if getFloatParam("positive_number", self.locals, 1)<=0:
            self.errors += ["specify a positive integer"]
                   
        # By calling readsHeaders() and writesHeaders(), we ensure that headers exist in the
        # input or output header list.
        #
        # If the header added has a standard header name (e.g. "OFFSET"), then it will have
        # the type of the standard header.
        # Similarly, if the header added has a non-standard header name but is provided as
        # an input header, then the header will retain the type of the input header.
        # If it is a non-standard name (e.g. "FOO") that doesn't already exist and we are writing it,
        # then we need to specify its type.
        self.readsHeaders(getStringParam("header_to_require", self.locals, ""))
        self.writesHeaders(getStringParam("header_to_add", self.locals, "")+":double")
       
        # Note that if you want to read header TRC_TYPE only if it is already on the input header list,
        # you do not need to call readsHeaders(TRC_TYPE, ...). Rather, just check whether header TRC_TYPE
        # already exists, and only read the header if it is present. If the header does not exist, the user
        # will not see any error message.
        self.readTrcType = self.canReadHeaders("TRC_TYPE")

        # Remove this header from the output header list. If it is not provided as an input
        # header, then this call will have no effect.
        self.removesHeaders(getStringParam("header_to_remove", self.locals, ""))
       
        # Set a "transient" parameter that is used in the set_il_xl parameter label
        is2d = getStringParam("is2d", self.globals, False)
        if is2d: self.locals["__ilxl_or_cmp"] = "CMP"
        else   : self.locals["__ilxl_or_cmp"] = "ILINE and XLINE"
       
        self.doExtractPicks = getBoolParam("do_extract_picks", self.locals, False)
        setVisible(self.locals, "extraction_tbl", self.doExtractPicks)
        if self.doExtractPicks:
            # Make sure that the grid geom is defined since we'll use it to determine
            # the center of the survey.
            if not self.gridgeom.defined():
                self.locals["__errors"] += ["Project grid geometry must be defined"]
           
            # Open the table to read        
            # If there's an error such as:
            # - extraction_tbl parameter is empty,
            # - read permission for extraction_tbl,
            # - or extraction_tbl does not exist in the file system,
            # this tool will flag an error and the flow will not execute.
            self.table = spatialTableParam(self, "extraction_tbl")
           
            ### This demonstrate how set warning messages.
            # If the table is initialized properly without any error,
            # set the warning flag if the dimensions does not match.
            # Note that setting a warning message will still allow the flow to run
            if self.table.ok():
                if self.table.getXName()!="XLINE":
                    self.locals["__warnings"] += ["expected XLINE but got "+self.table.getXName()+" from table"]
                if self.table.getYName()!="ILINE":
                    self.locals["__warnings"] += ["expected ILINE but got "+self.table.getYName()+" from table"]

    def startExecution(self):
        # You can check if a tool is running interactively by using the "interactive"
        # parameter.
        interactive = getBoolParam("interactive", self.globals, False)
        if not interactive and self.doExtractPicks:
            minxl, maxxl = self.gridgeom.minmaxXl()
            minil, maxil = self.gridgeom.minmaxIl()
            xl = (minxl+maxxl)/2
            il = (minil+maxil)/2
            print("Reading velocity function at ({0}, {1})".format(xl, il))
            # get the function at 20x coarser sampling than the traces
            vfn = self.table.getInterpolatedFunction(y=il, x=xl, m0=self.o1In, dm=20*self.d1In, nm=self.n1In/20)
            print("\t{0}".format(vfn))
   
    def finishExecution(self):
        # Don't write out our database table if the flow is running interactively.
        interactive = getBoolParam("interactive", self.globals, False)
        if interactive: return
       
    def execute(self, ntr, hds, trs):
        # This tool cannot provide traces if not provided with input traces,
        # so in that case just return immediately.
        #
        # A tool may be passed two values for ntr that are less than zero:
        #  - NO_MORE_TRACES (0) means that the flow will not provide the tool
        #    with any more traces. If the tool is unable to provide additional
        #    traces, it should return NO_MORE_TRACES.
        #  - NEED_TRACES (-2) means that the flow is asking the tool to provide
        #    traces (without receiving any traces as input). If the tool is
        #    unable to do so, it should return NEED_TRACES.
        if ntr<=0: return ntr
       
        # Loop over traces. For each trace, determine whether that trace is alive,
        # then copy it to the output buffer.
        #
        # Note that we do not always need to copy to the output buffer. This tool is
        # written as a "one set" tool, which means that the input and output of
        # KitchenSink both reside in the same array.
        #
        # The result of this is that a gather containing no dead traces will not require
        # that any headers or traces be copied. If the gather contains at least one dead
        # trace and the tool is removing dead traces, then all subsequent traces in the
        # gather will need to be copied into a lower index in the header and traces array.
        # See the usages of the itr_out variable in the loop below.
        itr_out = 0
        for itr in xrange(0, ntr):
            alive = hds[itr]["TRC_TYPE"]==1
            if alive and itr_out==itr:
                # Not much to do; the trace is already in the buffer in the correct position
                itr_out +=1
            elif alive:
                # We need to copy the trace and header from position itr to position itr_out
                trs[itr_out][:] = trs[itr][:]
                self.headersOut().copyBuffer(self.headersIn(), hds[itr], hds[itr_out])
                itr_out +=1
            else:
                # The trace is dead. Do nothing and don't bump itr_out
                pass

        # If the entire gather is filled with dead traces, then we have nothing to output. In that case,
        # we need to return -2, which asks the engine for more traces instead of 0, which would tell the
        # engine that we are finished.    
        if itr_out==0: return -2
        else         : return itr_out
