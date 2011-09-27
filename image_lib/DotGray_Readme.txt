<FILES>
DotGrayConverter.m
DotGrayFileWriter.m

Attached are MatLab routines to convert .gray & import the
files into Matlab as well at to write a matlab matrix into a .gray
format. In addition the coverter routine can also convert the file   into
a BMP format.

Both programs can be used as stand alones or incorporated into other
programs.

NOTE the directory should be changed  under the "Set Defaults" section
in each code to prevent from having to enter it each time the program   is
run. The directory should be the directory where the images are   stored
on the users computers. If IsOSX/isWin functions do not come with Matlab 
remove the if/else clause.