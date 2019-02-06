# Automated Ecologit GUI written using Tkinter, PIL (mostly) and matplotlib
# Author: James Beattie
#

########################################################################################################################
# Tkinter
from Tkinter import *
import ttk;

# Base Libraries
import os;
import pickle;
import csv;
import numpy as np
import sys;

# Image Processing Libraries
from PIL import Image, ImageTk;
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
import scipy.ndimage

########################################################################################################################

# Classes
########################################################################################################################

class DrawLabel:
    def __init__(self, line):
        self.line   = line

        # Coordinates of line segments and lengths of line segment vectors.
        self.xs     = list(line.get_xdata());
        self.ys     = list(line.get_ydata());
        self.coords = [];

        # Mouse click coordinates
        self.xclick = list(line.get_xdata());
        self.yclick = list(line.get_ydata());
        self.press  = None;

        # Store the number of lines
        self.NumberOfLines = 0;

        self.cidpress   = line.figure.canvas.mpl_connect('button_press_event', self.on_press);
        self.cidrelease = line.figure.canvas.mpl_connect('button_release_event', self.on_release);
        self.cidmotion  = line.figure.canvas.mpl_connect('motion_notify_event', self.on_motion);

    def on_press(self, event):
        if event.inaxes!=self.line.axes: return
        self.press = 1;
        self.xclick.append(event.xdata);
        self.yclick.append(event.ydata);
        self.NumberOfLines += 1;
        print('Button clicked at:{}'.format((event.xdata,event.ydata)))


    def on_motion(self, event):
        if event.inaxes!=self.line.axes: return
        if self.press is not None and self.NumberOfLines > 1:
            self.xs.append(event.xdata)
            self.ys.append(event.ydata)
            self.line.set_data(self.xs, self.ys)
            self.line.figure.canvas.draw()
            print((event.xdata,event.ydata), event)

    def on_release(self, event):
        print(event)
        self.press = None;
        self.line.figure.canvas.draw();

        if self.NumberOfLines > 1:
            self.coords.append(('Line Number:',self.NumberOfLines-1))
            self.coords.append(('xcoord','ycoord'))
            self.coords.append((self.xs,self.ys))

        self.xs = [];
        self.ys = [];



# Functions for buttons
########################################################################################################################
def line_self(line):

	# Coordinates of line segments and lengths of line segment vectors.
	xs     = list(line.get_xdata());
	ys     = list(line.get_ydata());
	coords = [];

	# Mouse click coordinates
	xclick = list(line.get_xdata());
	yclick = list(line.get_ydata());
	press  = None;

	# Store the number of lines
	NumberOfLines = 0;

	cidpress   = line.figure.canvas.mpl_connect('button_press_event', on_press);
	cidrelease = line.figure.canvas.mpl_connect('button_release_event', on_release);
	cidmotion  = line.figure.canvas.mpl_connect('motion_notify_event', on_motion);


def on_press(event):
	if event.inaxes!=line.axes: return
	press = 1;
	xclick.append(event.xdata);
	yclick.append(event.ydata);
	NumberOfLines += 1;
	print('Button clicked at:{}'.format((event.xdata,event.ydata)))


def on_motion(event):
	if event.inaxes!=line.axes: return
	if press is not None and NumberOfLines > 1:
		xs.append(event.xdata)
		ys.append(event.ydata)
		line.set_data(xs,ys)
		line.figure.canvas.draw()
		print((event.xdata,event.ydata), event)


def on_release(self, event):
	print(event)
	press = None;
	line.figure.canvas.draw();

	if NumberOfLines > 1:
		coords.append(('Line Number:',NumberOfLines-1))
		coords.append(('xcoord','ycoord'))
		coords.append((xs,ys))

	xs = [];
	ys = [];


def return_entry():
	global all_images

	"""Updates the directory based on the entry widget and clears the list"""
	try:
		content = directory_entry.get();
		directory.set(directory_entry.get());
		lbox.delete(0, END);
		files 	= os.listdir(directory.get());
		os.system('ls {}/*.png {}/*.JPG > imagefiles.txt'.format(directory.get(),directory.get()))

		for name in files:
		    lbox.insert('end', name);

		all_images = len(files)

	except OSError:
		print 'You have not chosen a valid directory.'
		pass


def onselect(evt):
	global index

	"""Save the items that are selected in the list wdiget and saves them to the image_var"""
	w 		= evt.widget;
	index 	= int(w.curselection()[0]);
	value 	= w.get(index);

	print ('You selected item %d: "%s"' % (index, value))

	listitem.set(value);

	image_var.set(directory.get() + '/' + value);


def image_press():
	"""If a new image is selected in the list, and then select image button is pressed it will
	change the image."""
	try:
		f 	= Figure(figsize=(6,6), dpi=150,tight_layout=True)
		a 	= f.add_subplot(111)

		global img

		img = Image.open(image_var.get())
		a.imshow(img)
		a.axis('off')

		# Matplotlib plot embedded in the canvas
		canvas = FigureCanvasTkAgg(f,ImageViewer)
		canvas.show()
		canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,E,W))
		#canvas.update()

		# Toolbar for plot
		toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer)
		toolbar.grid(column=1,row=2,sticky=(N,S,E,W))

	except ValueError:
		pass


def image_forward():
	"""The button that goes forward through the images"""
	try:
		f 			= open('imagefiles.txt','rb');
		lines 		= f.readlines();
		counter  	= 0;

		for files in lines:

			if files.rstrip() == image_var.get():

				try:
					image_var.set(lines[counter+1].rstrip());
				except IndexError:
					image_var.set(lines[0].rstrip())

			else:
				counter += 1;

		f 	= Figure(figsize=(6,6), dpi=150,tight_layout=True)
		a 	= f.add_subplot(111)

		global img

		img = Image.open(image_var.get())
		a.imshow(img)
		a.axis('off')

		# Matplotlib plot embedded in the canvas
		canvas = FigureCanvasTkAgg(f,ImageViewer)
		canvas.show()
		canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,E,W))

		# Toolbar for plot
		toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer)
		toolbar.grid(column=1,row=2,sticky=(N,S,E,W))

	except ValueError:
		pass


def image_backward():
	"""The button that goes backwards through the images"""
	try:
		f 			= open('imagefiles.txt','rb');
		lines 		= f.readlines();
		counter  	= 0;

		for files in lines:

			if files.rstrip() == image_var.get():

				try:
					image_var.set(lines[counter-1].rstrip());
				except IndexError:
					print('this works')
					image_var.set(lines[-1].rstrip());

			else:
				counter += 1;

		f 		= Figure(figsize=(6,6), dpi=150,tight_layout=True);
		a 		= f.add_subplot(111);

		global img

		img 	= Image.open(image_var.get());
		a.imshow(img);
		a.axis('off');

		# Matplotlib plot embedded in the canvas
		canvas = FigureCanvasTkAgg(f,ImageViewer);
		canvas.show();
		canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,E,W));

		# Toolbar for plot
		toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer);
		toolbar.grid(column=1,row=2,sticky=(N,S,E,W));

	except ValueError:
		pass


def ClassOnOff():
	""" Turn the classes on and off"""

	global img_downsample, img_d, img, PlotNumber

	PlotNumber 		= image_var.get().split('/')[-1]; 	#extract the file

	# Now remove the processed.JPG from the end
	PlotNumber 	= PlotNumber.replace('processed.JPG','',1)

	# Create the pickle file
	PickleFile 		= PlotNumber + '.JPG.pickle';

	# Create the directory file
	PickleDirectory = directory.get() + '/' + PickleFile;

	# ClassState check
	if ClassState.get() == 0:
		# If the ClassState is off, turn it on
		ClassState.set(1);

		# Update all button states
		SelectImageButton.config(state=DISABLED);
		NextImageButton.config(state=DISABLED);
		PreviousImageButton.config(state=DISABLED);
		NextClassButton.config(state=NORMAL);
		PreviousClassButton.config(state=NORMAL);

		print 'Trying to find pickle file...'

		try: # try to open the pickle file
			PFile 	= open(PickleDirectory,'rb');
			print('Found pickle file {}, now unpacking.'.format(PickleFile));
			PData 	= pickle.load(PFile);

			# Change the total class number
			ClassNumber.set(len(np.unique(PData[1])));
			ttk.Label(Class, text="Class {} / {}:".format(CurrentClass.get(),ClassNumber.get())).grid(column=1, row=1, sticky=W)
			print 'Unpacking finshed, there are {} classes.'.format(ClassNumber.get())

			# Find the size of the image and downsample by 5
			maxpix_d 	= [(map(max, zip(*PData[0]))[0]/5),(map(max, zip(*PData[0]))[1])/5];

			# Set this is a global variable so that I can act on it when I change classes
			img_d 		= np.zeros([maxpix_d[0],maxpix_d[1]]);


			# Read in the image from coordinates, populate img_d
			n = 0;
			for i in PData[0]:
				j 			= tuple(map(lambda x: x/5-1, i));
				img_d[j] 	= PData[1][n]+1;
				n 			+= 1;

		except IOError:
			print "There doesn't exist a pickle file for the image you have selected."
			pass
			ClassState.set(0);

	elif ClassState.get() == 1:
		ClassState.set(0);

		# Update all button states
		SelectImageButton.config(state=NORMAL);
		NextImageButton.config(state=NORMAL);
		PreviousImageButton.config(state=NORMAL);
		NextClassButton.config(state=DISABLED);
		PreviousClassButton.config(state=DISABLED);
		ClassDrawButton.config(state=DISABLED);

		ClassNumber.set('-');	# reset the total number of classes
		CurrentClass.set(0);	# reset the current class to 0
		CoverageVar.set(0.0);	# reset the coverage to -

		# Update the coverage environment
		CoverageLab = ttk.Label(Class, text="{}".format(CoverageVar.get()),width=20)
		CoverageLab.grid(column=1,row=12,sticky=W)

		# Update the total number of classses
		ttk.Label(Class, text="Class {} / {} :".format(CurrentClass.get(),ClassNumber.get())).grid(column=1, row=1, sticky=W);

		# Update the image so that it returns to the downsampled image.
		img = Image.open(image_var.get())
		a.imshow(img);
		a.axis('off');

		# Matplotlib plot embedded in the canvas
		canvas = FigureCanvasTkAgg(f,ImageViewer);
		canvas.show();
		canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,E,W));

		# Toolbar for plot
		toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer);
		toolbar.grid(column=1,row=2,sticky=(N,S,E,W));


def Masking(img,img_d):
	""" This function takes the downsampled class map and img, upsamples the class map and adds 3 channels to the map"""

	# Set the mask to be a logical array based upon the current class and upscale using nearest neighbours and by 5 orders (2225,3015)
	# image is (2448,)  so
	image_mask			= scipy.ndimage.zoom(img_d == CurrentClass.get(), 5, order=0)

	# Preallocate the mask
	mask 				= np.zeros_like(img);# Create a 3D Mask


	# Save Pickle/Image difference dimensions to reading in the first coordinate
	starty 	= (img.height - image_mask.shape[0])/2;
	endy 	= image_mask.shape[0] + starty;

	startx 	= (img.width - image_mask.shape[1])/2;
	endx	= image_mask.shape[1] + startx;

	# Fill in each of the
	for i in range(3):
		mask[starty:endy,startx:endx,i] = image_mask;

	return mask


def ClassForward():
	""" Move forward through the classes."""

	global img_d, img, mask


	if ClassState.get() == 1:

		ClassDrawButton.config(state=NORMAL)

		if CurrentClass.get() < int(ClassNumber.get()):
			# Set the class number to be equal to the first class
			CurrentClass.set(CurrentClass.get() +  1);
		elif CurrentClass.get() == int(ClassNumber.get()):
			CurrentClass.set(1);

		ttk.Label(Class, text="Class: {} / {}".format(CurrentClass.get(),ClassNumber.get())).grid(column=1, row=1, sticky=W)


		## Update canvas environment ##


		# Create the mask
		mask 	= Masking(img,img_d);

		a.imshow(img*mask);
		a.axis('off');

		# Matplotlib plot embedded in the canvas
		canvas = FigureCanvasTkAgg(f,ImageViewer);
		canvas.show();
		canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,E,W));

		# Toolbar for plot
		toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer);
		toolbar.grid(column=1,row=2,sticky=(N,S,E,W));


		## Upate Cover Metric ##
		CoverageVar.set( format((img_d==CurrentClass.get()).sum()/np.float((img_d!=0).sum()),'.4f') );
		ttk.Label(Class, text="{}".format(CoverageVar.get())).grid(column=1, row=12, sticky=W)


	elif ClassState.get() == 0:
		print 'You have to turn on the Class first to access a pickle file of classes.'


def ClassBackwards():
	""" Move backwards through the classes. """

	global img_d, mask

	if ClassState.get() == 1:

		ClassDrawButton.config(state=NORMAL)

		if CurrentClass.get() > 1:
			# Set the class number to be equal to the first class
			CurrentClass.set(CurrentClass.get() -  1);
		elif CurrentClass.get() == 1 or CurrentClass.get() == 0:
			# Set the class number to be back to the top
			CurrentClass.set(int(ClassNumber.get()));

		ttk.Label(Class, text="Class: {} / {}".format(CurrentClass.get(),ClassNumber.get())).grid(column=1, row=1, sticky=W)

		# Create the mask
		mask 	= Masking(img,img_d);

		f 	= plt.figure(figsize=(6,6), dpi=175, tight_layout=True)
		f.subplots_adjust(wspace=0, hspace=0)
		a 	= f.add_subplot(111)

		a.imshow(img*mask);
		a.axis('off');

		# Matplotlib plot embedded in the canvas
		canvas = FigureCanvasTkAgg(f,ImageViewer);
		canvas.show();
		canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,E,W));
		canvas.update()

		# Toolbar for plot
		toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer);
		toolbar.grid(column=1,row=2,sticky=(N,S,E,W));


		## Upate Cover Metric ##
		CoverageVar.set( format((img_d==CurrentClass.get()).sum()/np.float((img_d!=0).sum()),'.4f') );
		ttk.Label(Class, text="{}".format(CoverageVar.get())).grid(column=1, row=12, sticky=W)


	elif ClassState.get() == 0:
		print 'You have to turn on the Class first to access a pickle file of classes.'


def TrainingDataExtract(LabelCoords):
	"""This function takes the coordinates from the drawn line and saves a number of windowsize X windowsave Images
	into either GoodData directroy or BadData directory. """

	global PlotNumber, img, HomeDirectory, GlobalTrainingImageCounter

    # Define LineNumbers coordinates, extract them from the Label class
	LineNumbers  				= LabelCoords[0::3];
	coords      				= LabelCoords[2::3];
	GlobalTrainingImageCounter 	= 0;

    # Write to a directoy folder depending upon the input of the function.
	if LabelState.get() == 'Excellent':
		Directory = 'ExcellentData'
		print 'Labelled data will go into ExcellentData directory.'
	elif LabelState.get() == 'Misclassified':
		Directory = 'MisclassifiedData';
		print 'Labelled data will go into MisclassifiedData directory.'
	else:
		print 'errrrm...'

	try:
		# For each line number
		for LineNumber in xrange(0,len(coords)):
			xcoord = coords[LineNumber][0];
			ycoord = coords[LineNumber][1];

	        # size of window
			windowsize = 21;

	        # Define the x and y coordinates
			xTop    = [x + windowsize/2.0 for x in xcoord];
			xBottom = [x - windowsize/2.0 for x in xcoord];
			yTop    = [x + windowsize/2.0  for x in ycoord];
			yBottom = [x - windowsize/2.0  for x in ycoord];

			#print(xTop,xBottom)
	        # for each coordinate
			for i in xrange(0,len(xcoord)):
				TrainingData = img.crop(box=(xBottom[i],yBottom[i],xTop[i],yTop[i]));
				TrainingData.save('{}/{}_line{}_img{}_lab{}_class{}.png'.format(Directory,PlotNumber,LineNumber,GlobalTrainingImageCounter,LabelState.get(),CurrentClass.get()));
				GlobalTrainingImageCounter += 1;

	except IOError:
		pass


def DrawOnOff():
	""" DrawOnOff changes the drawing state, so that it is on or off."""

	global mask, img, line, line_coords

	if DrawState.get() == 0:
		DrawState.set(1);
		print('DrawState = {}'.format(DrawState.get()))

		#Disable stuff

		f 	= plt.figure(figsize=(6,6), dpi=175, tight_layout=True);
		f.subplots_adjust(wspace=0, hspace=0);
		a 	= f.add_subplot(111);
		a.imshow(img*mask);
		line, = a.plot([155],[118],'red'); # empty line
		a.axis('off');

		# Matplotlib plot embedded in the canvas
		canvas = FigureCanvasTkAgg(f,ImageViewer);
		canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,W,E));

		line_coords, click_coords = ClassLabeller(line);

		# Toolbar for plot
		toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer);
		toolbar.grid(column=1,row=2,sticky=(N,S,E,W));


	elif DrawState.get() == 1:
		DrawState.set(0);

		TrainingDataExtract(line_coords);

		print('DrawState = {}'.format(DrawState.get()))


def on_press(event):
	""" Define events on_press of the mouse when draw is on."""

	global xclick, yclick, NumberOfLines, press, line

	#if event.inaxes!=line.axes: return
	press = 1;
	xclick.append(event.xdata);
	yclick.append(event.ydata);
	NumberOfLines += 1;
	print('The number of lines: {}'.format(NumberOfLines))
	print('Button clicked at:{}'.format((event.xdata,event.ydata)))


def on_motion(event):
	""" Define events on_motion of the mouse when draw is on."""

	global xs, ys, line, NumberOfLines, press

	print NumberOfLines, press,

	#if event.inaxes!=line.axes: return
	if press is not None and NumberOfLines > 1:
		xs.append(event.xdata)
		ys.append(event.ydata)
		line.set_data(xs,ys)
		line.figure.canvas.draw()
		print((event.xdata,event.ydata), event)


def on_release(event):
	""" Define events on_release of the mouse when draw is on."""

	global press,line,NumberOfLines,coords,xs,ys

	print(event)
	press = None;
	line.figure.canvas.draw();

	if NumberOfLines > 1:
		coords.append(('Line Number:',NumberOfLines-1))
		coords.append(('xcoord','ycoord'))
		coords.append((xs,ys))

	xs = [];
	ys = [];


def ClassLabeller(line):
	""" This function takes the line class and extracts the x,y canvas coordinates of the drag motion,
	and also the click coordinates"""

	global xs, ys, coords, xclick, yclick, press, NumberOfLines

	xs     = list(line.get_xdata());
	ys     = list(line.get_ydata());
	coords = [];

	# Mouse click coordinates
	xclick = list(line.get_xdata());
	yclick = list(line.get_ydata());
	press  = None;

	# Store the number of lines
	NumberOfLines = 0;

	linecoords = line.figure.canvas.mpl_connect('button_press_event', on_press);
	cidrelease = line.figure.canvas.mpl_connect('button_release_event', on_release);
	cidmotion  = line.figure.canvas.mpl_connect('motion_notify_event', on_motion);

	return coords, (xclick, yclick)


def _quit():
	root.quit();     	# stops mainloop
	root.destroy();  	# this is necessary on Windows to
	plt.close('all');	# close all plots that seem to be coming up.


def calculate():
	print 1+1

# GUI BEGINS HERE
########################################################################################################################
HomeDirectory = '/Users/jamesbeattie/Documents/Research/2017_McCool_Clustering/GUI'


# OS command calls
########################################################################################################################

os.system('ls {}/*.png > imagefiles.txt'.format(HomeDirectory))

# Create a Bad Label Directory
if os.path.exists('MisclassifiedData') != True: os.system('mkdir MisclassifiedData');
# Create a Good Label Directory
if os.path.exists('ExcellentData') != True: os.system('mkdir ExcellentData')

########################################################################################################################

root = Tk()
root.title("Automated Ecologist Labelling GUI")
root.resizable(width=False, height=False)
root.geometry('1515x1125')

mainframe = ttk.Frame(root)
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
mainframe.columnconfigure(0, weight=1)
mainframe.rowconfigure(0, weight=1)


# Variables
ClassState 		= IntVar(); 			# 0: Class off, 1 Class On
DrawState		= IntVar();
LabelState 		= StringVar();	# this is the state of the label that will be draw on: Excellent or Misclassified
directory 		= StringVar(); 	# directory variable
listitem 		= StringVar(); 	# list item variable
label 			= StringVar(); 	# label variable
image_var		= StringVar();  # image variable
ClassNumber 	= StringVar(); 	# the total number of classes
CurrentClass 	= IntVar();	# the current class
rs 				= 10; 			# row span of the list and the image
CoverageVar 	= DoubleVar();	# the amount of coverage in the class


# Default Variable Values
directory.set(HomeDirectory);
label.set("forb");
image_var.set(HomeDirectory + "/Q1_2009.JPG_ENCODED_008.png");
ClassNumber.set('-');
CurrentClass.set(0);
CoverageVar.set(0.0);
DrawState.set(0);
ClassState.set(0);

# Directory handling ###################################################################################################
ttk.Label(mainframe, text="Directory:",font='Helvetica 14 bold').grid(column=1, row=1, sticky=W)
directory_entry = ttk.Entry(mainframe,width=154, textvariable=directory)
directory_entry.grid(column=2,row=1,columnspan=20)
dir_change = ttk.Button(mainframe, text="Change Directory", command=return_entry)
dir_change.grid(column=37, row=1,columnspan=2, sticky=(N,S,W,E))


#  List box  ###########################################################################################################
lbox 	= Listbox(mainframe, height=10)
files = os.listdir(directory.get())

for name in files:
    lbox.insert('end', name)

# list box
lbox.grid(column=1,row=3,rowspan=rs,columnspan=4,sticky=(N,S,E,W))

# Define scroll bar
scrollbar = Scrollbar(lbox)
scrollbar.pack(side=RIGHT, fill=Y)
lbox.config(yscrollcommand=scrollbar.set);
scrollbar.config(command=lbox.yview)

#scrollbar.config(command=lbox.yview)

# Next Class and Image Buttons #########################################################################################
Class = ttk.Label(mainframe);
Class.grid(column=37,row=3,sticky=W);
ttk.Label(Class, text="Class: {} / {}".format(CurrentClass.get(),ClassNumber.get()),font='Helvetica 14 bold').grid(column=1, row=1, sticky=W);
ClassOnOffButton = ttk.Button(Class, text="Class On / Off", command=ClassOnOff);
ClassOnOffButton.grid(column=1, row=2,columnspan=2, sticky=(N,S,E,W));
NextClassButton = ttk.Button(Class, text="Next Class", command=ClassForward,state=DISABLED);
NextClassButton.grid(column=1, row=3,columnspan=2, sticky=(N,S,E,W));
PreviousClassButton = ttk.Button(Class, text="Previous Class", command=ClassBackwards,state=DISABLED);
PreviousClassButton.grid(column=1, row=4,columnspan=2,sticky=(N,S,E,W));

ImageParms = ttk.Label(mainframe);
ImageParms.grid(column=37, row=4,sticky=(N,S,E,W));
ttk.Label(ImageParms, text="Image 0 / 90:",font='Helvetica 14 bold').grid(column=1, row=1,sticky=(N,S,E,W));
SelectImageButton 	= ttk.Button(ImageParms, text="Select Image", command=image_press);
SelectImageButton.grid(column=1, row=2,sticky=(N,S,E,W));
NextImageButton		= ttk.Button(ImageParms, text="Next Image", command=image_forward,width=20);
NextImageButton.grid(column=1, row=3,sticky=(N,S,E,W));
PreviousImageButton = ttk.Button(ImageParms, text="Previous Image", command=image_backward,width=20);
PreviousImageButton.grid(column=1, row=4,sticky=(N,S,E,W));


# Class label drop down menu ###########################################################################################
ttk.Label(Class, text="Class Label:",font='Helvetica 14 bold').grid(column=1, row=8, sticky=W);
ClassLabels = OptionMenu(Class, label, "forb", "grass", "litter","forb-grass","forb-litter","grass-litter",'ground');
ClassLabels.grid(column = 1, row = 9,sticky=(N,S,E,W),columnspan=2);
#ConfirmClassButton = ttk.Button(Class, text="Confirm", command=calculate,state=DISABLED)
#ConfirmClassButton.grid(column=1, row=10,sticky=(N,S,E,W))
ttk.Label(Class, text="Coverage:",font='Helvetica 14 bold').grid(column=1, row=11, sticky=W);
CoverageLab = ttk.Label(Class, text="{}".format(CoverageVar.get()),width=20);
CoverageLab.grid(column=1, row=12, sticky=W);
ttk.Label(Class, text="Class Drawing:",font='Helvetica 14 bold').grid(column=1, row=13, sticky=W);
ClassDrawButton = ttk.Button(Class, text="Draw On/Off", command=DrawOnOff,state=DISABLED,width=20);
ClassDrawButton.grid(column=1, row=14,sticky=(N,S,E,W));
LearningLabels = OptionMenu(Class, LabelState, 'Excellent','Misclassified');
LearningLabels.grid(column = 1, row = 15,sticky=(N,S,E,W),columnspan=2);

# Image viewer #########################################################################################################
ImageViewer = ttk.Label(mainframe);
ImageViewer.grid(column=5,row=3,rowspan=rs,columnspan=22);

f 	= plt.figure(figsize=(6,6), dpi=175, tight_layout=True);
f.subplots_adjust(wspace=0, hspace=0);
a 	= f.add_subplot(111);
img = Image.open(image_var.get());
a.imshow(img);
a.axis('off');

# Matplotlib plot embedded in the canvas
canvas = FigureCanvasTkAgg(f,ImageViewer);
canvas.show()
canvas.get_tk_widget().grid(column=1,row=1,sticky=(N,S,W,E));
# Toolbar for plot
toolbar = NavigationToolbar2TkAgg(canvas,ImageViewer);
toolbar.grid(column=1,row=2,sticky=(N,S,E,W));


# Export CSV ###########################################################################################################
DatasetParms = ttk.Label(mainframe)
DatasetParms.grid(column=37, row=5,sticky=(N,S,E,W))
ttk.Label(DatasetParms, text="Dataset:",font='Helvetica 14 bold').grid(column=1, row=2, sticky=(N,S,E,W))
ExportCSVButton = ttk.Button(DatasetParms, text="Export CSV",command=calculate,state=DISABLED,width=20)
ExportCSVButton.grid(column=1, row=3, sticky=(N,S,E,W))

Quit_ 	= ttk.Label(mainframe)
Quit_.grid(column=37, row=12,sticky=(N,S,E,W))
button 	= ttk.Button(Quit_, text='Quit', command=_quit,width=20).grid(column=1, row=1, sticky=(N,S,E,W))


########################################################################################################################
# Some kind of padding applied everything
for child in mainframe.winfo_children():
	child.grid_configure(padx=1, pady=1)


lbox.bind('<<ListboxSelect>>', onselect)
#root.bind('<Return>', image_press)


root.mainloop()
