#author narumeena
#description devlop automaticaly shell script of variaent calling pipeline







# purpose:   read a pipe-delimited text file (i.e., the "fields" are separated by
#            the tab \t
# usage:     perl main.pl


#!usr/bin/perl

$filename = 'sample.txt';

# use the perl open function to open the file
open(FILE, $filename) or die "Could not read from $filename, program halting.";

# loop through each line in the file
# with the typical "perl file while" loop
while(<FILE>)
{
    # get rid of the pesky newline character
    chomp;
    
    # read the fields in the current line into an array
    @fields = split('|', $_);
    
    # print the first field
    print "$fields[0]\n";
}
close FILE;
