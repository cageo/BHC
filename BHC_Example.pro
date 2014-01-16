;------------------BHC_Example--------------------------------------------------
; PURPOSE:
;    A simple shell program that sets up the required input variables and structures 
;        required by the BHC progam (BHC_FindCorrection).  This is the first place to 
;        start for modification for running with another reconstructor, although it is
;        likely that other changes will be required in the core BHC_FindCorrection.pro code.
;        
; SIDE EFFECTS:
;    - A new directory is created (IDLBHC_AUTO) and populated with reconstructed files. 
;    - detailed_bhc_results.txt lists transformation functions attempted during reconstruction.  
;    - The following are displayed: 
;       - The initial uncorrected image with ROI locations
;       - Window showing reconstructed images as beam hardening corrections are applied - these images
;           will be downsampled as indicated by bhcParams_in.downsample parameter
;       - The final corrected image (full resolution)
;       - Final graph of BHC transform function
; -------------------------------------------------------------------------------

PRO BHC_Example
  
  ; Path to command line reconstructor
  reconutct_in = 'U:\SoftwareLibrary\Recon\TestBed\reconutct.exe'
  
  ; Path to a sinogram file that requires beam hardening correction
  ; BHC_FindCorrection assumes this is in a TIFF file format specific to a BIR scanner.
  raw_f_in = "U:\GroupShare\Romy\FOR_RICH\BHC_Paper\shellTest\apatite0098-Raw"
  divide = strpos(raw_f_in, '-', /REVERSE_SEARCH)    ; assumes filename ends with -Raw
  outfile = strmid(raw_f_in, 0, divide) + '.tiff'
  
  ; Read the sinogram file to get data min and max
  BHC_ReadBIRSinogramSmart, raw_f_in, data, head1, head2
  STACK_MAX = Max(data,MIN=STACK_MIN)
  
  ; Create the Reconutct_Params and initialize
  reconParams_in = Ptr_New({BHC_Reconutct_Params})
  *reconParams_in = BHC_Init_Params(*reconParams_in)
  (*reconParams_in).RawFile = raw_f_in
  (*reconParams_in).ImageFile = outfile
  
  ; Create the BHC_Params and initialize to single ROI mode, 5 points, and reconstruct file
  bhcParams_in = Ptr_New({BHC_Params})
  (*bhcParams_in).roiMode = 0       ; 0: use single, selected ROI for merit function
                                    ; 1: treat all ROIs as separate materials for merit function
                                    ; 2: treat all ROIs as one material (combine) for merit function
                                    ; 3: treat ROIs as two materials for merit function
                                    ;    Implementation of this mode requires examining all regions and 
                                    ;     and splitting them into two materials (for instance, based on their
                                    ;     average CT value).  Not included in this source code example.
  (*bhcParams_in).numPts = 7        ; Number of coefficients to use in function fit (usually 3-15)
  (*bhcParams_in).pos_2nd_der = 1   ; Enforce positive second derivative at all points (0:NO 1:YES)
  (*bhcParams_in).downsample = 2    ; 0: no downsampling of image during reconstructions
                                    ; 2: 2x downsampling of image during reconstructions
                                    ; 4: 4x downsampling of image during reconstructions
                                    ; 8: 8x downsampling of image during reconstructions 
  (*bhcParams_in).reconstruct = 1   ; Reconstruct raw file(s) (0:NO 1:YES) Must be set to 1 for this example.
           
  ; Write ascii file for command line call of reconstructor - this is specific to the BIR reconstruction program
  dir = Path_Sep()
  divide = strpos(raw_f_in, dir, /REVERSE_SEARCH) + 1
  params_fn = strmid(raw_f_in, 0, divide) + 'params_idlbhc.txt'
  if (file_test(params_fn) eq 1) then file_delete, params_fn 

  BHC_Write_Reconutct_File, params_fn, (*reconParams_in)

  ; Reconstruct and load the uncorrected image  
  if (!VERSION.OS_FAMILY eq 'Windows') then begin
     SPAWN, /HIDE, EXIT_STATUS=err, reconutct_in + ' -f ' + params_fn, result, error_result
  endif else begin
      SPAWN, /NOSHELL, EXIT_STATUS=err, [reconutct_in, ' -f ', params_fn], result, error_result
  endelse
  if (err ne 0) then begin
    if (error_result ne '') then err_msg = error_result else err_msg = result
    result = dialog_message(['Reconutct Error:', err_msg], /Error)
  endif
  orig_image = read_tiff((*reconParams_in).ImageFile)
  sz = size(orig_image)
  (*reconParams_in).Matrix = sz[1]  ; record reconstructed field of view size; assumes square image

  ; Create an ROI object so the program has an area to minimize on (in merit function)
  ; This ROI represents an area in the reconstructed sinogram file that should be homogeneous in 
  ;   greyscale but is not due to beam-hardening.
  ; It will have to be modified for different input files (raw_f_in)
  ; Easiest to set and modify using a GUI with a drawing area (i.e. using IDL graphics objects) 
  roi_in = Ptr_New({BHC_ROI_Params})
  (*roi_in).group = Obj_New('IDLgrROIGroup')
  Xcoord_conv = [0.5, 1.0]
  Ycoord_conv = [0.5, 1.0]
  regionColor = [ 54,173,122] ; dark pale green
  (*roi_in).selected= OBJ_NEW('IDLgrROI', COLOR=regionColor, $
                  XCOORD_CONV=Xcoord_conv, YCOORD_CONV=ycoord_conv)
  xmin = 542.5
  xmax = 626.5
  ymin = 695.5
  ymax = 866.5
  (*roi_in).selected->ReplaceData,[[xmin,ymin],[xmin,ymax],[xmax,ymax],[xmax,ymin]]    
  zero = 0     
  (*roi_in).group->Add, (*roi_in).selected, POSITION=zero  
  
  ; Create an object graphics display showing original image and ROI locations
  imageSize = [sz[1], sz[1]]
  oWindow = OBJ_NEW('IDLgrWindow', RETAIN = 2, DIMENSIONS = imageSize, $
        TITLE = 'Uncorrected Image and ROI Locations')
  oView = OBJ_NEW('IDLgrView', VIEWPLANE_RECT = [0., 0., imageSize])
  oModel = OBJ_NEW('IDLgrModel')
  image = bytscl(orig_image)
  oImage = OBJ_NEW('IDLgrImage', image, /GREYSCALE)
  oModel -> Add, oImage
  oView -> Add, oModel
  oWindow -> Draw, oView
  oModel->Add, (*roi_in).group
  oWindow->Draw, oView
    
  ; Finally, launch the optimization routine to find the correction
  coefs = BHC_FindCorrection(raw_f_in, roi_in, reconParams_in, bhcParams_in, stack_min, stack_max, reconutct_in, orig_image)

  ; final result is a transform function representing the beam hardening correction.
  ; It is defined by a set of points (y's = coefs and x's = original CT data value) that define a cubic spline function 
  ;   in BHC_FindCorrection.
  print, "Optimal BHC defined by points: "
  print, "X:", (*(*bhcParams_in).x)
  print, "Y:", coefs
    
  ; clean up
  Obj_Destroy, (*roi_in).group
  Obj_Destroy, (*roi_in).selected
  ; These objects not destroyed to keep original image up after program completion
;  Obj_Destroy, oImage
;  Obj_Destroy, oModel
;  Obj_Destroy, oView
;  Obj_Destroy, oWindow
  Ptr_Free, roi_in
  Ptr_Free, bhcParams_in
  Ptr_Free, reconParams_in
END
