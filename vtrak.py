#! /usr/local/bin/python3

import os
import sys
import warnings
import yaml
import datetime
import math 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import skimage.io
import skimage.filters
from skimage import measure, data, color, util
from skimage.feature import canny
from skimage.draw import circle_perimeter
from skimage.transform import hough_circle, hough_circle_peaks
matplotlib.style.use('ggplot')

if sys.platform == "linux" or sys.platform == "darwin" \
or sys.platform == "linux2":
	try:
		from dialog import Dialog
		dlg = Dialog(dialog="dialog")
		rows, cols = os.popen("stty size","r").read().split()
		rows = int(rows); cols = int(cols)

		def fs():
			path = './'
			entries = os.listdir(path)
			tagtuples = []
			for i in entries:
				if '.tif' in i and '_overlay' not in i:
					tagtuples.append((i,i,"off"))
			code, paths = dlg.buildlist("Select Files",rows-10,cols-10,rows-14,tagtuples)
			if code == Dialog.OK:
				os.system('clear')
				return paths
			if code == Dialog.CANCEL:
				os.system('clear')
				return None
	except:
		print("Dialog module not found. Defaulting to CLI.")

if sys.platform == "win32":
    try:
        import wx
        def fs():
            app = wx.App(None)
            style = wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE
            dialog = wx.FileDialog(None, 'Open',wildcard='*.tif',style=style)
            mult = None
            if dialog.ShowModal() == wx.ID_OK:
                try:
                    paths = dialog.GetPath()
                    mult = False
                except:
                    paths = dialog.GetPaths()
                    mult = True
            else:
                paths = None
            dialog.Destroy()
            return paths        
    except:
        print("WX module not found. Defaulting to CLI.")

class exp_param:
    def __init__(self,root=None,params={},model=None,obj=None,exp=None,um2pix=None):
       self.root   = root
       self.model  = model
       self.obj    = obj
       self.exp    = exp
       self.um2pix = um2pix
       if self.root != None:
         print("> Importing YAML data ("+self.root+")...",end="")
         with open(self.root+'.yml','r') as yfile:
          self.params = yaml.load(yfile)
         if self.params['mscope']['model'].upper().lower() == 'nikon':
          self.model = 'nikon'
         if self.params['mscope']['obj'].upper().lower() == '40x':
          self.obj = '40x'
          self.um2pix = 3.9052  # Pixels per micron 
         if self.params['mscope']['obj'].upper().lower() == '20x':
          self.obj = '20x'
          self.um2pix = 3.9052  # Pixels per micron
         if self.params['mscope']['obj'].upper().lower() == '10x':
          self.obj = '10x'
          self.um2pix = 3.9052
         self.channels = {}
         for key, value in self.params['mscope']['channels'].items():
          newkey = key.upper().lower()
          self.channels[newkey] = value
         self.sol = {}
         for key, value in self.params['sol_content'].items():
          newkey = key.upper().lower()
          self.sol[newkey] = value
         self.pip = {}
         for key, value in self.params['pipette'].items():
          newkey = key.upper().lower()    
          self.pip[newkey] = value
         self.ves = {}
         for key, value in self.params['vesicle'].items():
          newkey = key.upper().lower()
          self.ves[newkey] = value
         self.calc = {}
         for key, value in self.params['calc'].items():
          newkey = key.upper().lower()
          self.calc[newkey] = value
         self.dt = self.params['mscope']['dt']
         self.date = str(self.params['date']).upper().lower()
         self.exp = str(self.params["type"]).upper().lower()
         try:
          self.tension = []
          tension_arr = self.params['tension']
          for i in tension_arr:
              for key, value in i.items():
                 self.tension.append((key,value))
         except:
          pass
         print("Done.")
       
    def output_params(self,output='stdout'):
       if output.upper().lower() == 'stdout':
         print("\nEXPERIMENTAL DETAILS")
         print("Date  >  "+self.date,end="\n")
         print("Name  >  "+self.root,end="\n")
         print("Type  >  "+self.exp,end="\n\n")
         print("MICROSCOPE")
         print("Model >  "+self.model,end="\n")
         print("Obj   >  "+self.obj,end="\n")
         print("Channels:",end="\n")
         for key, value in self.channels.items():
             print("  {} > {}".format(key,value),end="\n")
             print("\nSAMPLE")
             print("Solution:",end="\n")
             print("  In:")
         for key, value in self.sol['inside'].items():
             print("    {} {}".format(key,value),end="\n")      
             print("  Out:",end="\n")
         for key, value in self.sol['outside'].items():
             print("    {} {}".format(key,value),end="\n")      
    
class vesicleData:
    def __init__(self,imgstack=None,params=None):
       self.imgstack=imgstack
       self.params=params
       self.overlay=None
       self.circle_runs=[]
       self.prots=None
       self.time=None
       self.relative_area=None
       self.unit='um'

    def load_stack(self):
       self.imgstack = skimage.io.MultiImage(self.params.root+'.tif')[0]
    
    def gen_circles(self,bar=None):
       print("> Detecting vesicle contours...")
       self.overlay = self.imgstack
       if self.params.calc['pix_err'] == 0:
         pix_range = [0]
       else: 
         pix_range=[-self.params.calc['pix_err'],0,self.params.calc['pix_err']]
       for i, err in enumerate(pix_range):
         print("> Applying cicular Hough transform {}...".format(i+1))
         # Convert um to pix
         pip_r = math.ceil(self.params.pip['diameter']*self.params.um2pix/2)+err
         ves_r = math.ceil(self.params.ves['diameter']*self.params.um2pix/2)+err
         circles = {'vesicle':[],'protrusion':[]}
         for i in range(round(len(self.overlay))):
          if progress:
              percent = round(100*float(i/len(self.imgstack)))
              progress.update(pcent=percent)
          filt = canny(self.overlay[i],sigma=2,\
                 low_threshold=self.params.calc['canny_lo'],\
                 high_threshold=self.params.calc['canny_hi'])
          lrg_rd_range = np.array([ves_r,ves_r])
          sm_rd_range = np.array([pip_r,pip_r])
          lrg_res = hough_circle(filt,lrg_rd_range)
          acc, lx, ly, lrad = hough_circle_peaks(lrg_res,\
                        lrg_rd_range, total_num_peaks=1)
          sm_res = hough_circle(filt,sm_rd_range)
          acc, sx, sy, smrad = hough_circle_peaks(sm_res,\
                                     sm_rd_range, total_num_peaks=1)
          lcircy, lcircx = circle_perimeter(int(np.average(ly)),\
                     int(np.average(lx)),int(np.average(lrad)))
          self.overlay[i][lcircy, lcircx] = self.params.calc['c_int']
          smcircy, smcircx = circle_perimeter(int(np.average(sy)),\
                       int(np.average(sx)),int(np.average(smrad)))
          self.overlay[i][smcircy, smcircx] = self.params.calc['c_int']
          circles['vesicle'].append((np.average(ly),\
                       np.average(lx),np.average(lrad)))
          circles['protrusion'].append((np.average(sy),\
                          np.average(sx),np.average(smrad)))
          self.circle_runs.append(circles)
    
    def gen_prot_data(self):
       prot_array = []
       for i, circles in enumerate(self.circle_runs):    
         prots = []
         time  = []
         lrg_dat = circles['vesicle']
         sm_dat = circles['protrusion']
         for j, (lrg, sml) in enumerate(zip(lrg_dat,sm_dat)):
          d = math.sqrt((lrg[1]-sml[1])**2 + (lrg[0]-sml[0])**2)
          prot_length = d + sml[2]-lrg[2]
          prot_length = math.fabs(prot_length)    
          prots.append(prot_length/self.params.um2pix)
          if i == len(self.circle_runs):
              time.append((j*self.params.dt))
         prot_array.append(np.array(prots))
       self.prots = np.average(prot_array,axis=0)
       self.time = time
    
    def gen_prot_data(self):
       prot_arr = []
       for circles in self.circle_runs:
         prots = []
         time  = []
         lrg_dat = circles['vesicle']
         sm_dat = circles['protrusion']
         for i, (lrg, sml) in enumerate(zip(lrg_dat,sm_dat)):
          d = math.sqrt((lrg[1]-sml[1])**2 + (lrg[0]-sml[0])**2)
          prot_length = d + sml[2]-lrg[2]
          prot_length = math.fabs(prot_length)    
          prots.append(prot_length/self.params.um2pix)
          time.append((i*self.params.dt))
         prot_arr.append(prots)
       self.prots = np.average(prot_arr,axis=0)
       self.time = time

    def gen_area_data(self):
       dA_arr = []
       A0 = 0
       for circles in self.circle_runs:
         dA = []
         lrg_dat = circles['vesicle']
         sm_dat  = circles['protrusion']
         for i,(lrg,sml) in enumerate(zip(lrg_dat,sm_dat)):
          dp = sml[2]*2/self.params.um2pix
          dv = lrg[2]*2/self.params.um2pix
          dA.append(math.pi*dp*(1-(dp/dv))*(self.prots[i]-self.prots[0]))
         A0 = A0 + 4*math.pi*((lrg_dat[0][2]/self.params.um2pix)**2)
         dA_arr.append(dA)
       A0 = A0/len(self.circle_runs)
       dA_avg = 100*(np.average(dA_arr,axis=0))/A0
       self.relative_area = np.array(dA_avg)

    def conv_pix(self):
       if self.unit == 'um':
         self.prots *= self.params.um2pix
         self.realtive *= (self.params.um2pix)**2
         self.unit = 'pix'
       else:
         pass

    def conv_um(self): 
       if self.unit == 'pix':
         self.prots /= self.params.um2pix
         self.realtive /= (self.params.um2pix)**2
         self.unit = 'um'
       else:
         pass

    def log_prot(self):
       f = open(self.params.root+'_prot.dat','w')
       f.write("PROTRUSION LENGTH DATA\n")
       f.write("SOURCE: {}\n".format(self.params.root))
       f.write(str(datetime.datetime.now())+"\n\n")
       f.write("TIME    PROTRUSION LENGTH\n\n")
       for i in range(len(self.prots)):
         f.write("{:4}    {:8}\n".format(self.time[i],self.prots[i]))
       f.close() 

    def log_area(self):
       f = open(self.params.root+'_area.dat','w')
       f.write("# RELATIVE AREA CHANGE DATA\n")
       f.write("# SOURCE: {}\n".format(self.params.root))
       f.write("# "+str(datetime.datetime.now())+"\n\n")
       f.write("# TIME    RELATIVE AREA CHANGE (PERCENTAGE)\n\n")
       for i in range(len(self.prots)):
         f.write("{:4}    {:8}\n".format(self.time[i],self.relative_area[i]))
       f.close()

    def load_prot_dat(self):
       self.time, self.prots = np.loadtxt(self.root+'_prot.dat',\
              dtype='float',comments='#',usecols=(0,1),unpack=True)

    def load_area_dat(self):
        self.time, self.relative_area = np.loadtxt(self.root+'_area.dat',dtype='float',comments='#',usecols=(0,1),unpack=True)

class pbar:
    def __init__(self,length=30,marker=">",lm="|",rm="|"):
       self.length=length
       self.marker=marker
       self.lm=lm
       self.rm=rm

    if sys.platform == "linux" or "linux2" or "darwin":
       def col_reset(self):
         print('\x1b[0m')
       
       def update(self,pcent):
         red    = '\x1b[0;31;49m'
         yellow = '\x1b[0;33;49m'
         green  = '\x1b[0;32;49m' 
         sub = math.floor(pcent*self.length/100)
         if sub <= self.length/3:
          print(red+(self.marker*sub)+(" "*(self.length-sub))+\
              self.rm+" {}%".format(pcent),end="\r  "+self.lm)
         if sub > self.length/3 and sub <= 2*self.length/3:
          print(yellow+(self.marker*sub)+(" "*(self.length-sub))+\
              self.rm+" {}%".format(pcent),end="\r  "+self.lm)
         if sub > 2*self.length/3:
          print(green+(self.marker*sub)+(" "*(self.length-sub))+\
              self.rm+" {}%".format(pcent),end="\r  "+self.lm)
         if sub == self.length:
          self.col_reset()
    
    if sys.platform == "win32":
       def update(self,pcent):
         sub = round(pcent*self.length/100)
         print((self.marker*sub)+(" "*(self.length-sub))+self.rm,end="\r"+self.lm)

def area_routine(config):
	data = vesicleData(params=config)
	data.load_stack()
	data.gen_circles(bar=progress)
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		skimage.io.imsave(data.params.root+'_overlay.tif',data.overlay)
	data.gen_prot_data(); data.log_prot()
	data.gen_area_data(); data.log_area()
	print("> Data logged.")
	return data

#### MAIN ####

if '--fs' not in sys.argv:
    config = exp_param(root=str(sys.argv[-1]).split('.')[0])
    progress = pbar()
    data = area_routine(config)

    plt.plot(data.time,data.relative_area,'o',markersize=4)
    plt.title(data.params.root)
    plt.xlabel('Time [s]')
    plt.ylabel('Relative Area Change '+r'$\left(\frac{\Delta A}{A}\right) \%$')
    plt.savefig(data.params.root+"_area.pdf")
    plt.show()

if '--fs' in sys.argv:
    progress = pbar()
    inputfiles = fs()
    if len(inputfiles) != 1:
       datalist = []
       for i in range(len(inputfiles)):
         config = exp_param(root=str(inputfiles[i]).split('.')[0])
         data = area_routine(config)
         datalist.append(data)
       
       plt.figure('Multiple Vesicle Analysis')
       print("\nEnter Plot title: ",end="")
       TITLE = input()   
       plt.title(str(TITLE))
       for i,data in enumerate(datalist):
         plt.plot(data.time,data.relative_area,'o',markersize=4)
       plt.xlabel('Time [s]')    
       plt.ylabel('Relative Area Change '+r'$\left(\frac{\Delta A}{A}\right) \%$')
       if TITLE != '':
         plt.savefig(TITLE+"_area.pdf")
       plt.show()
    else:
       config = exp_param(root=str(inputfiles[0]).split('.')[0])
       data = area_routine(config)

       plt.plot(data.time,data.relative_area,'o',markersize=4)
       plt.title(data.params.root)
       plt.xlabel('Time [s]')
       plt.ylabel('Relative Area Change '+r'$\left(\frac{\Delta A}{A}\right) \%$')
       plt.savefig(data.params.root+"_area.pdf")
       plt.show()
