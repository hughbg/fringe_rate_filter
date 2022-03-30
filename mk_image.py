import os


# Input is x.uvfits 

print("IMPORT ----------\n")
importuvfits(fitsfile="x.uvfits", vis="x.ms")
print("CLEAN-------------\n")
tclean(vis="x.ms", pbcor=True, weighting="briggs", imagename="x", cell="0.14deg", imsize=512, gridder="widefield",  pblimit=-0.01)

print("EXPORT-------------\n")
exportfits(imagename="x.image", fitsimage="x.fits")
print("PNG----------------\n")
os.system("convert x.fits output.png")
os.system("mv x.fits output_im.fits")
os.system("mv x.image output.image")

