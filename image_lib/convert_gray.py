#Bijan Mapar, 2010
#Needs Python Image Library (PIL) Installed
import Image
import math
import Tkinter
import tkFileDialog
import sys
filename=''
while filename=='':
    root = Tkinter.Tk()
    filename = tkFileDialog.askopenfilename()
    root.destroy()
    
f = open(filename,'rb')
data = []
try:
    header = f.read(256/8)
    byte=1
    while byte!="":
        byte = f.read(1)
        #print byte
        if byte=="":
            break
        pixel = int(ord(byte))
        data.append((pixel,pixel,pixel))
finally:
    f.close()

dimensions = int(round(math.sqrt(len(data))))
im = Image.new("RGB",(dimensions,dimensions))
im.putdata(data[0:dimensions*dimensions])
im.save(filename[0:-4]+"png")
