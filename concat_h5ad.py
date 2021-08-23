#!/usr/bin/env python3

import sys, getopt
import scanpy as sc, anndata as ad, os
from numpy import loadtxt

def main(argv):
   infile = ''
   outfile = ''
   try:
      opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
   except getopt.GetoptError:
      print('concat_h5ad.py -i <file of h5ad files> -o <outfile h5ad>')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print('concat_h5ad.py -i <file of h5ad files> -o <outfile h5ad>')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         infile = arg
      elif opt in ("-o", "--ofile"):
         outfile = arg

   if infile == None:
      print("Must specify infile")
      sys.exit(2)      

   if outfile == '':
      print("Must specify outfile")
      sys.exit(2)   

   # Run code data
   ################

   # read list of h5a files
   h5adfiles = open(infile, "r").read().split('\n')

   # remove blank lines
   h5adfiles = list(filter(None, h5adfiles))

   print( "Combine", len(h5adfiles), "H5AD files")

   # Create anndata array
   adatas = [ad.read(file) for file in h5adfiles] 

   # Concatenate data
   print(" Concatenating data...")
   adata = ad.concat(adatas[:])

   # Write to disk
   print(" Writing to H5AD...")
   adata.write( outfile )


if __name__ == "__main__":
   main(sys.argv[1:])


