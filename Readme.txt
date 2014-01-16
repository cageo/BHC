;***************************************************************************************************
;***************************************************************************************************
;
; Supplementary IDL code for:
; 
; Beam Hardening Correction for X-ray Computed Tomography of Heterogeneous Natural Materials, submitted
;     to Computers and Geosciences, 2014
;     
;     Richard A. Ketcham and Romy D. Hanna
;     ketcham@jsg.utexas.edu
;
; Introduction:
;   This code represents a very simplified version of the original program written as a user-guided 
;   beam hardening correction (BHC) application.  The original program includes an extensive GUI as the
;   main driver of the application, with the code presented here representing the central core
;   of the BHC optimization algorithm.  It consists of three IDL files:
;         
;   BHC_Example.pro:  A simple shell program that sets up the required input variables and structures 
;                     required by the BHC progam (BHC_FindCorrection)
;                           
;   BHC_FindCorrection.pro:  The main BHC program, with BHC_FindCorrection() as the entry point.  All supporting
;                            functions are also included in this file, listed in the correct compile order. It returns
;                            a floating point array of the BHC coefficients (Y's) and defining the transform function for
;                            the supplied sinogram file and parameters.
;                                 
;   BHC_Amoeba.pro:   This is a modified version of the original IDL Amoeba function (from IDL 7.1 but likely 
;                     unchanged in later IDL releases). IDL(r) is a registered trademark of Exelis Inc. and IDL(r) 
;                     software is copyrighted by Exelis Visual Information Solutions, Inc.
;     
;
;   Before compiling, the IDL file [IDL_Distribution_Path]/lib/obsolete/tiff_dump.pro needs to be edited as follows:
;      1.  comment out lines 249 and 250
;         ; if (typ le 0) or (typ gt 5) then $
;         ;		message,'Illegal type code, ifd = '+string(index)
;      2. change line 347 to "if arg_present(txt) then txt = line ; else print,line" (other print statements
;         in this file can also be commented out if desired)
;			
;   This is a change specific to the BIR TIFF sinogram file, and may not necessarily be needed for other sinogram
;   file types (see below for more information).
;
;   IDL compile order of included files: 1. BHC_FindCorrection.pro, 2. BHC_Amoeba.pro 3. BHC_Example.pro
;     
;   Although it should compile with no issues (in IDL 7.0 and later), the code cannot be run "as is".  At the very least it 
;   requires updated paths to a sinogram file and command-line reconstruction program in BHC_Example.pro. It was 
;   specifically tailored to a circa 2001 BIR reconstruction program and its TIFF file format for sinograms, so it is 
;   likely that certain program elements related to I/O will have to be modified (in BHC_FindCorrection.pro) to work
;   with a different reconstruction program and its associated file types. 
;	
;   This simplified code was not as extensively tested as the full program, so we apologize in advance for any
;   missing or extraneous variables that resulted from the stripping down of the original code.  Please contact us if 
;   you have questions.
;       
; Additional Notes:
;   - The full BHC optimization program has several reporting functions and GUI elements tied into it 
;         which have been removed here.
;   - In addition, for the sake of brevity, almost all error checking (I/O, argument checks) was removed.
;   - 'Raw' file refers to the uncorrected, raw sinogram file (which is a TIFF file out of the BIR scanner). 
;   - This program requires a separate reconstruction program that is callable via command line.  
;
;***************************************************************************************************  
;**************************************************************************************************
