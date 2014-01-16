;------------------BHC_Cubic_BHC-----------------------------
; PURPOSE:
;    Given a 2D sinogram, transforms data by doing a cubic spline interpolation
;      on input [x,y] points which define transformation function.  So, in the end
;      the data will be transformed using a cubic spline function
;
; INPUTS:
;  data:  2D array of raw sinogram data.  Data type shouldn't matter, it is recast to double before
;         calculation.
;  x:  vector of x  (transformation function defined as [x,y])
;  y:  vector of y 
;  
; OUTPUTS:
;  data:  replaces input data array with BH corrected data
; ---------------------------------------------------------
PRO BHC_Cubic_BHC, data, x, y
       
    sz = SIZE(data)
    num_data = sz[4]
    
    x_sz = SIZE(x)
    y_sz = SIZE(y)
    
    data = Double(data)
    
    ; initialize cubic spline func
    y2 = spl_init(x,y, /DOUBLE)
   
    data = spl_interp(x,y,y2,data, /DOUBLE)
  
END


;------------------BHC_Log_And_Scale-----------------------------
;
; PURPOSE:
;    Given a 2D array of data, scales the data (if SHORT) and 
;    calculates the log before performing BHC.  If the inverse argument is set to 1, 
;    it calculates the inverse log (exp) and reverses the scale.
;
; INPUTS:
;   data:  2D array of data.  Should be 16 bit short int or 32 bit float.
;   scale:  value to scale short (uint) data by.  This is usually retrieved from the raw
;          sinogram file (tif tag Raw_Gain_Scale). 
;   bits_per_sample:  either 16 (short) or 32 (float)
;   inverse:  if 0, scale the data and take the log.  If 1, do the 
;            inverse log and then inverse scale
;  
; OUTPUTS:
;   data:  replaces input data array with scaled/log or inverse log/scaled data
; ---------------------------------------------------------
PRO BHC_Log_And_Scale, data, scale, bits_per_sample, inverse

  if ( scale lt 0.0001 ) then scale = 1.0D
  
  if (inverse eq 0) then begin      
        ; scale data and take log
        if (bits_per_sample ne 32) then begin    ; UINT case
           data = DOUBLE(data)
           dtmp = DOUBLE(data / scale)           
        endif else begin                         ; FLOAT case
           dtmp = DOUBLE(data)                  
        endelse
        if (Min(dtmp) GE 1.e-10) then data = DOUBLE(-alog(dtmp)) $
        else begin
          ; First zero out values in <data> less than 1.e-10
          zeroes = where(dtmp lt 1E-10, numZeroes, COMPLEMENT=nonZeroes, NCOMPLEMENT=numNonZeroes)     
          if (numZeroes GT 0) then data[zeroes] = 0       
          ; Zero is a strange value for this case, as that's the alog of 1, which is an ostensibly reasonable data value
          ; Possible explanation is that scaled data might stay in range 0-1 -- but in that case,zeroes would be converted to maximum value!
           if (numNonZeroes GT 0) then data[nonZeroes] = DOUBLE(-alog(dtmp[nonZeroes]))
        endelse
  endif else begin   
        ; take inverse log and scale back to short data type if applicable
        dtmp = DOUBLE(data)
        replace = where(dtmp ne 0, count)
        ; Again: This leaves zeroes intact, whereas they could in theory be converted back to 1, possibly a reasonable data value
        ; This logic also presumes that any number that starts as zero remains a zero after processing; it's unclear whether this is a safe assumption.
        if (count ne 0) then data[replace] = DOUBLE((exp(-dtmp[replace])))
        ; if data not float scale data and convert to short
        if (bits_per_sample ne 32) then begin
            dtmp = DOUBLE(scale*data)
            data = UINT(round(dtmp))
        endif   
  endelse
  
END

;------------------BHC_Construct_Transform-----------------------------
;
; PURPOSE:
;   Compiles the y coordinates of the point-based transform from the
;   output of the iterative function.
;   The first and last coefficients are the y-function values, and 
;   the intermediate coefficients are the fractional increment of 
;   rise in y coordinate value (monotonically increasing).
;
; INPUTS:
;  coefs:  Coefficients of the iterative function.   
;  INVERSE: Optional keyword that will transform real numbers to mono space. 
;           Default behavior is transforming from mono space to real numbers.
;  
; OUTPUTS:
;  y: An array with the y coordinates
; ---------------------------------------------------------

Function BHC_Construct_Transform, coefs, INVERSE=inverse

    if keyword_set(inverse) then begin    ; inverse transforms real numbers into mono space
      y = coefs 
      n = N_Elements(y)
      for i=1,n-2 do y[i] =(coefs[i]- coefs[i-1])/(coefs[n-1]-coefs[i-1])
    endif else begin    ; regular mode transforms from mono to regular space
      y = coefs 
      n = N_Elements(y)
      for i=1,n-2 do y[i] = y[i-1] + coefs[i]*(y[n-1]-y[i-1])
    endelse
    
    return, y
End                          


;----------------------------BHC_ReadBIRSinogramSmart------------------------------------
; PURPOSE:
;   Reads a sinogram produced by BIR ACTIS software.  Sinograms may be 16-bit (from II system) or
;   32-bit (from high-energy system); 32-bit sinogram TIFF headers indicate that they're longword
;   integers, but they're really floating point, so this routine compensates.  Also reads and
;   returns the header.
;   This version is more robust than the original in that it can cope with files without the
;   header all at the beginning of the file.  For this case, it returns 2 headers: an initial
;   header and a final one.
;
; INPUT PARAMETERS:
;   filename: Name of sinogram, including path
;   data: Variable into which the sinogram data will be placed
;   head1, head2: Variables into which the sinogram header information will be placed.
;     head1 is at the beginning of the file and head2 is at the end (making it a footer,
;     I suppose). 
;   sz: Optional. Array returned by "Size" function for a sinogram with the same type and
;     dimensions as the one being read; it makes file reading more efficient.  If an unused
;     variable is passed in, it is filled with the "Size" results for the file being read.
;
; KEYWORD PARAMETERS:
;   TAGS: Set this to an array of tag values to read.  This array will be replaced by an array
;     of double-precision numbers with the corresponding tag values.  Note that if the tag is not
;     a number, this will return an error and stop everything.
;   BAD_TAGS: This variable will be replaced by a byte array that corresponds to tags not found 
;        (1 (true) = bad tag value or tag not defined).  This array matches the TAGS array in
;        size and order.  
;
; ------------------------------BHC_DeallocBIRTags--------------------------------
; Deallocate the tags array returned by BIRTiffDump TAG_VALUE
Pro BHC_DeallocBIRTags, tags
  sz1 = Size(tags)
  if (sz1[1] EQ 0) then return
  for i=0,sz1[1]-1 do Ptr_Free, tags[i]
  tags = 0
End

Pro BHC_ReadBIRSinogramSmart, fileName, data, head1, head2, sz, TAGS=tags, BAD_TAGS=bad_tags
  if (N_Params() LT 4) OR (N_Elements(sz) LE 1) then begin
    tiff = Read_Tiff(fileName)
    sz = Size(tiff)
    tiff = 0
  endif
  makeFloat = (sz[3] EQ 3) OR (sz[3] EQ 13)
  if Keyword_Set(tags) then getTags = [273, tags] else getTags = [273]

  tiff_dump, filename, GET_TAG_VALUE=getTags, TAG_VALUE=vals
  offset = *(vals[0])
  if Keyword_Set(tags) then begin
    tags = DblArr(N_Elements(tags))
    bad_tags = bytarr(N_elements(tags))
    pointer_valid = ptr_valid(vals)
    for i=1,N_Elements(tags) do begin
      if (pointer_valid[i] eq 1) then begin
        tags[i-1] = *(vals[i]) 
      endif else begin 
        tags[i-1] = 0
        bad_tags[i-1] = 1
      endelse
    endfor
  endif
  BHC_DeallocBIRTags, vals

  openr, 1, fileName
  totHeadSize = (fstat(1)).size-(sz[1]*sz[2]*(makeFloat ? 4 : 2))
  head1Size = offset
  head2Size = totHeadSize - head1Size
  head1 = BytArr(head1Size)
  head2 = (head2Size GT 0) ? BytArr(head2Size) : -1
  data = makeFloat ? FltArr(sz[1], sz[2]) : UIntArr(sz[1], sz[2])
  readu, 1, head1
  readu, 1, data
  if (head2Size GT 0) then readu, 1, head2
  close, 1
  free_lun, 1
End

; ------------------------------BHC_WriteBIRSinogramSmart--------------------------------
; PURPOSE:
;   Writes a BIR ACTIS sinogram so the reconstruction software can still understand it.
;
; INPUTS:
;   filename: Name of sinogram file, including path
;   data: sinogram data (unsigned int or float array)
;   head1, head2: Headers at beginning and end of file; these are obtained from
;     ReadBIRSinogramSmart
;   sz: Output of IDL Size() function on original data, to make sure data is
;     written in same format
; ---------------------------------------------------------
Pro BHC_WriteBIRSinogramSmart, fileName, data, head1, head2, sz
  openw, 1, filename
  if (sz[3] EQ 12) then writeu, 1, head1, UINT(data) $
  else writeu, 1, head1, FLOAT(data)
  if (N_Elements(head2) GT 1) then writeu, 1, head2
  close, 1
  free_lun, 1
End

;------------------BHC_Write_TIFF_Tag-----------------------------
;
; PURPOSE:
;   Writes a tiff tag value to the designated raw sinogram file (which is
;   incidentally a tiff file).  Currently only writes UINT and Double tag values.
;
; INPUTS:
;  filename:  Name of existing tiff file (sinogram) to modify
;  tag_id:  decimal id of tag to modify  (most specific to BIR tiff file)
;  value:   value to write
;  type:  Integer value type (implemented as needed).  For simplicity 
;           just following IDLs type codes returned by size() function: 
;           12 : uint
;           5: double
;  
; OUTPUTS:
;  none; a tiff file header is modified
; ---------------------------------------------------------
PRO BHC_Write_TIFF_Tag, filename, tag_id, value, type

    forward_function tiff_ulong, tiff_uint
    
    ; open file, read header
    openu, lun, filename, error = i, /GET_LUN
    if i lt 0 then begin ;OK?
      message, 'Error opening ' + filename
      return
    endif
    hdr = bytarr(8)    ;Read the header
    readu, lun, hdr

    typ = string(hdr[0:1])   ;Either MM or II file type flag
    if (typ ne 'MM') and (typ ne 'II') then begin
      print,'TIFF_READ: File is not a Tiff file: ', filename
      return
    endif

    ; read ifd
    offs = tiff_ulong(hdr, 4)   ;Offset to IFD
    point_lun, lun, offs        ;point to IFD
    a = bytarr(2)               ;Entry count array
    readu, lun, a               
    count = tiff_uint(a,0)      ;count of tif tags

    ifd = bytarr(count * 12 + 4) ;Array for IFD's + Offset of next IFD
    readu,lun, ifd              ; read it
    
    ; search for tif tag in ifd
    for i=0, count-1 do begin 
        tag = tiff_uint(ifd, i*12)
        if (tag eq tag_id) then break
    endfor
    if (i ge count-1) then begin 
      message, 'Error, could not find tif tag ' + string(tag_id)
      return
    endif
    
    ; write value to file at appropriate offset in ifd
    if (type eq 12) then begin
      point_lun, lun, 2 + offs + (i*12+8)       ;repoint to value in IFD
      writeu, lun, uint(value)
    endif else if (type eq 5) then begin
      point_lun, lun, 2 + offs + (i*12+8)       ;repoint to value in IFD
      offset = 0UL
      readu, lun, offset ; since double is > 4 bytes, tag field has offset to value location
      point_Lun, lun, offset
      writeu, lun, double(value)
    endif else begin
       message, 'Warning, unknown type in BHC_Write_TIFF_Tag: ' + string(type)
    endelse
    
    close, lun
    free_lun, lun
    
END

                            
;------------------BHC_Do_BHC-----------------------------
; PURPOSE:
;   Takes the given raw file and does a beam hardening correction
;   with the coefficients indicated.  Before doing BHC, scales
;   and takes the log of the raw sinogram data.  The corrected raw data is 
;   written back to the file as float, and the tiff tag BITS_PER_SAMPLE is modified if neccessary.
;   
; INPUTS:
;  raw_f:  Name of existing raw sinogram file to correct
;  coefs:  vector (of any size) of BHC coefficients  
;  outfile: Name of file that will be written out (fully qualified)
;  
; OUTPUTS:
;  a new, bhc'd raw file is written to a new directory IDLBHC
; ---------------------------------------------------------
PRO BHC_Do_BHC, raw_f_in, coefs, outfile

    COMMON DO_BHC_ARGS, last_raw_f_opened, scaled_data, head1, head2, sz, bits_per_sample
          
    COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
          LUN, show, y0, yN, bhcParams, matMeans, reconutct
  
    recon = bhcParams.reconstruct 
    x = *(bhcParams.x)
    
    skip = 0  
    ; Check to see if this file is already open and converted (log/scale)
    if (n_elements(last_raw_f_opened) ne 0) then $
        if (last_raw_f_opened eq raw_f_in) then begin
          data = scaled_data 
          skip = 1
        endif
    if (skip ne 1) then begin
        ; Read file with TIFF reader and get header size
        tiff = Read_Tiff(raw_f_in)
        sz = Size(tiff)
      
        ; read raw file, these tiff tags are specific to BIR sinogram file
        ; tags:  258 = bits_per_sample, if==32 float, otherwise short 
        ;             (should be 16, but BIR code doesn't tie it to this...)
        ;        32832 = raw_gain_scale
        ;        32830 = image scale
        tags = [258,32832,32830]
        BHC_ReadBIRSinogramSmart, raw_f_in, data, head1, head2, sz, TAGS=tags, BAD_TAGS=bad_tags
        
        ; get bits_per_sample
        bits_per_sample = tags[0]
        
        ; might be possible that for float data scale may be undefined - yes
        if (bad_tags[1] ne 1) then raw_scale = tags[1] else raw_scale = 1
        
        ; get image scale or error out - required
        if (bad_tags[2] ne 1) then orig_scale = tags[2] else begin
            print, 'BHC_Do_BHC ERROR: Could not get image scale from sinogram file: ' + raw_f_in
            print, 'BHC NOT APPLIED'
            return
        endelse
        
        BHC_Log_And_Scale, data, raw_scale, bits_per_sample, 0
       
        ; update data sz[3] to reflect float data (type 4)
        sz[3] = 4
     
        last_raw_f_opened = raw_f_in
        scaled_data = data
        
    endif
  
    y = BHC_Construct_Transform(coefs) 
    BHC_Cubic_BHC, data, x, y
    
    ; if not reconstructing immediately afterwards, apply inverse log/scale
    if not (recon) then begin 
        BHC_Log_And_Scale, data, raw_scale, bits_per_sample, 1
        ; update data sz[3] to reflect original data type
        if (bits_per_sample eq 32) then sz[3] = 4 else sz[3] = 12
    endif 
       
    BHC_WriteBIRSinogramSmart, outfile, data, head1, head2, sz 
    
    ; if sinogram file was not 32 bit float, then modify tiff tag (Bits_per_sample, #258)
    ;   to indicate it is now so.
    if ((bits_per_sample ne 32) AND recon) then BHC_Write_TIFF_Tag, outfile, 258, 32, 12
  
 END
 
;------------------BHC_ROI_Merit_Func-----------------------------
; PURPOSE:
;   Finds the standard deviation of given ROI(s) after correcting beam hardening with given
;   function, defined in terms of function points (x,y).  Uses cubic spline interpolation for data transformation
;   Also, since this is a function to be used by an the IDL optimization function amoeba() everything
;   else except for the BHC coefs has to be accessed via a COMMON block during use with the amoeba() function.
;
; INPUTS:
;   y:  vector that, along with corresponding x's saved in ROI_Merit_Func_Args (fixed, based input data range), define
;       the BHC correction function
;   init_image: this is an optional parameter used only when setting up image statistics in BHC_FindCorrection().  Not
;               used during amoeba() calls.
;
; OUTPUTS:
;    standard deviation of the ROI(s)
; ---------------------------------------------------------
FUNCTION BHC_ROI_Merit_Func, y, INIT_IMAGE=initImage

  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
  
  if Keyword_Set(initImage) then begin
        slice = initImage
  endif else begin
      ; append endpoints to y
      new_y = [y0, y, yN]
      
      BHC_Do_Recon, new_y
      iter = iter + 1
      
      ; read in image
      slice= Read_Tiff(reconParams.ImageFile)  
  endelse
    
  sz = size(slice)
  
  ; the objRefs retrieved here is a collection of IDLgrROI objects (managed by the GUI)
  ; mode = 0 signifies only one ROI used for guiding the correction, otherwise it's a 
  ;    set of ROIs that define one or multiple materials
  if (bhcParams.roiMode EQ 0) then begin      ; selected ROI
    objRefs = [roi.Selected]
    nRegions = 1
  endif else objrefs = roi.group->Get(/ALL, COUNT=nRegions)
  
  ; calculate number of separate regions (materials) for image statistics
  if ((bhcParams.roiMode EQ 2) OR $  ; Joint ROI mode - assumes all ROIs are in the same material
    (bhcParams.roiMode EQ 0)) then nMaterials = 1  $  ; selected ROI only
  else if (bhcParams.roiMode EQ 3) then nMaterials = 2 $ ; 2 material ROI mode
  else nMaterials = nRegions   ; independent mode - each ROI a separate material
  
  ; Build Final Image Mask for image statistics
  final_mask = replicate(0B, sz[1], sz[2])
  for i=0, nRegions-1 do begin
    ; get ROI mask
    objrefs[i]->GetProperty, INTERIOR=interior
    if (interior eq 1) then continue  ; We'll deal with these later
    mask = objrefs[i]->ComputeMask(DIMENSIONS=[sz[1],sz[2]], MASK_RULE=2)
    if (nMaterials eq 1) then $     ; one material, add mask to final_mask
      final_mask[where(mask eq 255)] = 255 $
    else if (nMaterials eq nRegions) then $  ; independent mode, all different materials
      final_mask[where(mask eq 255)] = i + 1  $
    else begin   ; add mask according to material's uvalue (material index from 1-n)
      objrefs[i]->GetProperty, UVALUE=materialID
      final_mask[where(mask eq 255)] = materialID
    endelse
  endfor
  
  ; Now go through again and just delete any interior masks
  for i=0, nRegions-1 do begin
    objrefs[i]->GetProperty, INTERIOR=interior
    if (interior eq 1) then begin   ; region is interior, have to invert mask before adding to current mask
      objrefs[i]->SetProperty, INTERIOR=0     ; interior ROIs don't work well with ComputeMask
      mask = objrefs[i]->ComputeMask(DIMENSIONS=[sz[1],sz[2]], MASK_RULE=2)
      final_mask[where(mask eq 255)] = 0
      objrefs[i]->SetProperty, INTERIOR=1
    endif
  endfor
  
  ; create image statistics arrays for nMaterials
  indMerits = replicate(0.0D, nMaterials)
  
  ; compute image statistics
  labeled = (nMaterials eq 1)? 0:1
  ; FYI, for multiple materials and thus regions (LABELED=1) in mask, rest of mask region (0) counts as first region
  IMAGE_STATISTICS, slice, MASK=final_mask, STDDEV=stddevs, MEAN=means, LABELED=labeled, COUNT=counts
  if Keyword_Set(initImage) then begin
      matMeans = means 
      return, 1.e20
   endif
   
  ; Adjust the image so the means match original
  if (N_Elements(means) EQ 1) then calcSlice = slice*matMeans/means $
  else begin
    line = LINFIT(matMeans,means)
    calcSlice = (slice-line[0])/line[1]
  endelse
  IMAGE_STATISTICS, calcSlice, MASK=final_mask, STDDEV=stddevs, MEAN=means, LABELED=labeled, COUNT=counts
  
  if (nMaterials GT 1) then begin  ; remember, ignoring first region (0) for nMaterials >1
    mean_zero = where(means[1:nMaterials] eq 0, count, COMPLEMENT=mean_not_zero)
    if (count ne 0) then begin
      indMerits[mean_zero] = 1.e20 
      if (count NE nMaterials) then indMerits[mean_not_zero] = stddevs[mean_not_zero+1]/means[mean_not_zero+1]
    endif else indMerits = stddevs[1:nMaterials]/means[1:nMaterials]
  endif else if (means EQ 0) then indMerits = 1.e20 else indMerits = stddevs/means
  
  ; now calculate merit 
  if (nMaterials GT 1) then merit = Total(indMerits) else merit = indMerits
  
  ; Report progress
  BHC_Inverse_Reporting_Iterate, merit, new_y, slice
    
  return, merit
END

;------------------BHC_Write_Reconutct_File-----------------------------
;
; PURPOSE:
;   Given a filename and BHC_Reconutct_Params structure, writes the parameters
;   to a tab-deliminated text file that can be supplied to reconutct with -f.
;   If the file already exists, it appends the parameters as a new line to the end
;   of the file.  If it does not exist, it creates a new file.
;
; INPUTS:
;  filename:  Name of file to write.  Can be a new file or existing file.
;  params:  BHC_Reconutct_Params struct
; ---------------------------------------------------------
PRO BHC_Write_Reconutct_File, filename, params

  ; writes a single param struct to a tab deliminated text file
    if (file_test(filename)) then begin 
      openu, lun, filename, /GET_LUN, /APPEND 
    endif else begin
      openw, lun, filename, /GET_LUN 
      writeu, lun, ';parameter file created by BHC_Write_Reconutct_File.pro' + string(10B)
      writeu, lun, ';'
      writeu, lun,  tag_names(params) + string(9B)
    endelse
    
    writeu, lun, string(10B)
    
    for i=0, N_TAGS(params)-1 do begin
       writeu, lun, strtrim(params.(i),2) + string(9B)
    endfor
    
    close, lun
    free_lun, lun  
END

;------------------BHC_Do_Recon-----------------------------
; PURPOSE:
;   Applies the beam hardening correction and reconstructs the image with call to
;       external reconstruction program (path to program save in reconutct variable
;       in ROI_MERIT_FUNC_ARGS common block)
;
; INPUTS:
;   y:  vector of y's that, along with corresponding x's (fixed, based on
;       input data range), define the BH correction function
;   rescale: Flag to signal that image should be rescaled (only done after optimal BHC is determined)
;           calculates a modified scale and offset so scaling of the image post-BHC is similar to the pre-BHC image
;   
; ---------------------------------------------------------
Pro BHC_Do_Recon, y, RESCALE=rescale 

  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  fn = fn_base + "_BHC" + StrTrim(String(iter,format='(I04)'),2) + "-Raw" 
  outfile = new_dir + fn
  
  ; do BHC on input raw file and write resulting raw sinogram file
  BHC_Do_BHC, raw_f, y, outfile
  
  ; write parameter file for reconstruction
  reconParams.SkipScale = 1
  params_fn = new_dir + 'params_idlbhc.txt'
  
  ; if parameter file already exists, delete
  if (file_test(params_fn) eq 1) then file_delete, params_fn
  reconParams.RawFile = outfile
  ext = strpos(outfile, '-', /REVERSE_SEARCH)
  reconParams.ImageFile = strmid(outfile, 0, ext) + '.tif'
  BHC_Write_Reconutct_File, params_fn, reconParams
  
  ;  call Reconutct to reconstruct bhc file with tab-deliminated file
  if (!VERSION.OS_FAMILY eq 'Windows') then begin
    SPAWN, /HIDE, EXIT_STATUS=err, reconutct + ' -f ' + params_fn, result, error_result
  endif else begin
    SPAWN, /NOSHELL, EXIT_STATUS=err, [reconutct, ' -f ', params_fn], result, error_result
  endelse
  if (err ne 0) then begin
    if (error_result ne '') then err_msg = error_result else err_msg = result
    message, string(13B) + 'RECONTUCT Error:'  + strjoin(err_msg, string(13B), /SINGLE)
  endif
  
  ;  read in new reconstructed image to get new image mean, min, max (only done once final BHC determined)
  if keyword_set(rescale) then begin
    slice= Read_Tiff(reconParams.ImageFile)
    IMAGE_STATISTICS, slice, MEAN=new_image_mean, MAXIMUM=new_image_max, MINIMUM=new_image_min
    ; get original offset and scale
    if (reconParams.ImgScl ne 'NA') then orig_scale = Float(reconParams.ImgScl)
    if (reconParams.ImgOff ne 'NA') then orig_offset = Float(reconParams.ImgOff)
    if ((n_elements(orig_offset) eq 0) OR (n_elements(orig_scale) eq 0)) then begin  ; get from raw file
      tags = [32830, 32831]
      BHC_ReadBIRSinogramSmart, reconParams.RawFile, data, TAGS=tags, BAD_TAGS=bad_tags
      ; get image scale 
      if (n_elements(orig_scale) eq 0) then if (bad_tags[0] ne 1) then orig_scale = tags[0]
      ; get image offset 
      if (n_elements(orig_offset) eq 0) then if (bad_tags[1] ne 1) then orig_offset = tags[1]
    endif
    
    ; get 'real min' i.e. min value not 0 (canvas)
    hist = histogram(slice, MAX=20000)
    ind = where(hist[1:*] GT 0)
    new_image_min = ind[0]+1
    
    ; calculate new scale (that avoids saturation and clipping)
    mod_scale  = (orig_scale * 1.0)* (61000.0 - 2000.0)/(new_image_max - new_image_min)
    mod_scale = (floor(mod_scale/100.0))*100  ; round down to nearest 100
    BHC_Write_TIFF_Tag, reconParams.RawFile, 32830, mod_scale, 5
    
    ; calculate new offset (that avoids saturation and clipping)
    mod_offset = 2000.0 - (mod_scale*1.0)*(new_image_min-orig_offset)/(orig_scale*1.0)
    mod_offset = (ceil(mod_offset/100.0))*100  ; round up to nearest 100
    BHC_Write_TIFF_Tag, reconParams.RawFile, 32831, mod_offset, 5
    
    if (mod_offset LE 0) then mod_offset = 1
    
    ; delete original reconustructed image
    file_delete, reconParams.ImageFile
    ; rewrite params file so scale and offset read from raw file header instead of param file
    file_delete, params_fn
    imgScl_save = reconParams.ImgScl
    imgOff_save = reconParams.ImgOff
    reconParams.ImgScl = 'NA'
    reconParams.ImgOff = 'NA'
    BHC_Write_Reconutct_File, params_fn, reconParams
    
    ; reconstruct image again with modified scale
    if (!VERSION.OS_FAMILY eq 'Windows') then begin
      SPAWN, /HIDE, EXIT_STATUS=err, reconutct + ' -f ' + params_fn, result, error_result
    endif else begin
      SPAWN, /NOSHELL, EXIT_STATUS=err, [reconutct, ' -f ', params_fn], result, error_result
    endelse
    if (err ne 0) then begin
      if (error_result ne '') then err_msg = error_result else err_msg = result
      message, string(13B) + 'RECONTUCT Error:'  + strjoin(err_msg, string(13B), /SINGLE)
    endif
    
    ; save scale and offset
    bhcParams.mod_scale = mod_scale
    bhcParams.mod_offset = mod_offset
    ; restore scale and offset in reconParams
    reconParams.ImgScl = imgScl_save
    reconParams.ImgOff = imgOff_save
    
  endif   ; end of Rescale test case
  
End

;------------------BHC_Eval_Simple_Poly-----------------------------
;
; PURPOSE:
;    ; evaluates simple polynomial
;         y = A + Bx + Cx^2 + Dx^3....
;
; INPUTS:
;  data: vector of x values to evaluate
;  coefs: vector of coefficients of polynomial
;
; OUTPUTS:
;  vector of resulting values
; ---------------------------------------------------------
FUNCTION BHC_Eval_Simple_Poly, data, coefs

  ; data is a vector of data
  sz = SIZE(data)
  
  ; coefs are the coefficients of the polynomial to evaluate for each data point in vector data
  c_sz = SIZE(coefs)
  
  ; create sum vector and data polynomial arrays
  sum = DBLARR(sz[1])
  runX = REPLICATE(1.0D, sz[1], sz[2])
  
  for x = 0,  c_sz[1]-1, 1 do begin
    sum += (coefs[x] * runX)
    runX *= data
  endfor
  
  return, sum
END

;------------------BHC_Define_X-----------------------------
; PURPOSE:
;       Defines x which is the range of values over which BHC_Merit_Func will be minimized. 
;       Saved in common block ROI_MERIT_FUNC_ARGS variable bhcParams.
;
; INPUTS:
;   stack_min: Assuming that the input image (sinogram) is one of a stack of 2D slices that make up a full 
;             3D volume, this is the minimum of the entire stack so a BHC that applies to the entire stack
;             data space can be found.
;             
;   stack_max: Assuming that the input image (sinogram) is one of a stack of 2D slices that make up a full 
;             3D volume, this is the maximum of the entire stack so a BHC that applies to the entire stack
;             data space can be found.
;
; OUTPUTS:
;     none
;
; SIDE EFFECTS:
;     x is calculated and saved in COMMON BLOCK DO_BHC_ARGS.  Also parameters
;     from raw file saved in same common block so don't have to open raw file and
;     scale/log data again.  Scaled data saved in scaled_data.
; ---------------------------------------------------------
PRO BHC_Define_X, stack_min, stack_max
  
  COMMON DO_BHC_ARGS, last_raw_f_opened, scaled_data, head1, head2, sz, bits_per_sample
    
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  num_pts = bhcParams.numPts
  
  ; Read file with TIFF reader and get header size
  tiff = Read_Tiff(raw_f)
  sz = Size(tiff)
  
  ; read raw sinogram file
  ; tags:  258 = bits_per_sample, if==32 float, otherwise short
  ;              (should be 16, but BIR code doesn't tie it to this...)
  ;        32832 = raw_gain_scale
  tags = [258,32832]
  BHC_ReadBIRSinogramSmart, raw_f, data, head1, head2, sz, TAGS=tags, BAD_TAGS=bad_tags
  
  ; log and scale data for bhc
  if (bad_tags[0] ne 1) then bits_per_sample = tags[0] 
   
  ; might be possible that for float data scale may be undefined 
  if (bad_tags[1] ne 1) then raw_scale = tags[1] else raw_scale = 1
  BHC_Log_And_Scale, data, raw_scale, bits_per_sample, 0
  
  ; Use global max and min (this is so spline will encompass all sinogram values in data set)
  dmax = DOUBLE(-alog(STACK_MIN/raw_scale))
  dmin = DOUBLE(-alog(STACK_MAX/raw_scale))
  
  dscale = (dmax-dmin)/(num_pts-1) 
  x = dmin+IndGen(num_pts)*dscale  
  
  ; save x scale to bhcParams
  bhcParams.x = Ptr_new(x)
END

;;------------------BHC_Initial_Y-----------------------------
;
; PURPOSE:
;   Finds the initial y range (initial point P0 for amoeba) based on
;       a BH correction (defined as coefs to a simple polynomial.
;
; INPUTS:
;   coefs:  Coefficients of a simple polynomial BHC that is evaluated on x range to get
;           reasonble initial y's
;
; OUTPUTS:
;    intial y's to be used in the amoeba merit function (BHC_ROI_Merit_Func)
; ---------------------------------------------------------

FUNCTION BHC_Initial_Y, coefs
    
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  x = *(bhcParams.x)
  num_pts = bhcParams.numPts
  
  c_sz = size(coefs)
  sum = replicate(0.0D, num_pts)
  runX = replicate(1.0D, num_pts)
  
  ;evaluate simple polynomial with coefs
  for i = 0,  c_sz[1]-1, 1 do begin
    sum += (coefs[i] * runX)
    runX *= x
  endfor
  
  sum = sum*x[num_pts-1]/sum[num_pts-1]  
  p = sum
  for i=1,num_pts-2 do p[i] = (sum[i]-sum[i-1])/(sum[num_pts-1]-sum[i-1])
  sum = p
  
  ; since ends fixed save y0 and yN and trim endpoints from y
  y0 = x[0]
  yN = x[num_pts-1]
  
  return, sum[1:num_pts-2]
 
END


;------------------BHC_Define_Y_Bounds-----------------------------
;
; PURPOSE:
;   Finds the global bounds scale to bound the BHC_ROI_Merit_Func minimization function
;
;
; INPUTS:
;  num_pts: the number of points that the BHC_ROI_Merit_Func will be
;             fitting.  So, the scale bounds matrix will be size num_pts*2
;  upr_bnd_coefs: this is the simple polynomial that GLOBALLY defines the upper
;             bound of the function.  Note that for <0, this will turn into the lower
;             bound
;  lwr_bnd_coefs: this is the simple polynomial that GLOBALLY defines the lower
;             bound of the function.  Note that for <0, this turns into the upper bound
;
; OUTPUTS:
;    vector of bounds for y, size num_pts*2
; ---------------------------------------------------------
FUNCTION BHC_Define_Y_Bounds, num_pts, upr_bnd_coefs, lwr_bnd_coefs

  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  x = *(bhcParams.x)
  
  upr_bnd = BHC_Eval_Simple_Poly(x, upr_bnd_coefs)
  lwr_bnd = BHC_Eval_Simple_Poly(x, lwr_bnd_coefs)
  
  mono_lower_def = 0.2 
  
  ; place into scale vector
  ; ends fixed so scale shorter by 2
  scale = REPLICATE(0.0D, (num_pts-2)*2)
  ; all points set to alternative monotonic function form
  for i=0, num_pts-3, 1 do begin
    scale[i*2] = mono_lower_def/(num_pts-1-i)    
    scale[i*2+1] = 1.0/(num_pts-1-i)
  endfor
  
  return, scale
END

;;------------------BHC_Define_P0_and_Scale-----------------------------
;
; PURPOSE:
;       Defines the original starting point (initial BHC function to try) and scale for the amoeba function.
;       Note that we have slightly modified the original IDL BHC amoeba function to use scale as the 
;       upper and lower bounds of BHC coefficient space.
;
; INPUTS:
;   stack_min: Assuming that the input image (sinogram) is one of a stack of 2D slices that make up a full 
;             3D volume, this is the minimum of the entire stack so a BHC that applies to the entire stack
;             data space can be found.
;             
;   stack_max: Assuming that the input image (sinogram) is one of a stack of 2D slices that make up a full 
;             3D volume, this is the maximum of the entire stack so a BHC that applies to the entire stack
;             data space can be found.
;
; OUTPUTS:
;    initial guess:  The initial BHC function to try (P0 for amoeba function)
;    scale: bounds of the BHC function (scale for amoeba function)
;     
; SIDE EFFECTS: 
;     various ROI_MERIT_FUNC_ARGS Common block variables initialized for the merit function
;
; ---------------------------------------------------------
PRO BHC_Define_P0_and_Scale, stack_min, stack_max, initial_guess, scale 

  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  ; upper and lower bound BHCs
  upr_bnd_coefs = [0.0,0.98,0.0,0.35]
  lwr_bnd_coefs = [0.0,0.01,0.0,0.0]
  
  ; initial guess BHC
  init_coefs = [0.0, 0.7, 0.15, 0.0]
  
  BHC_Define_X, stack_min, stack_max
  initial_guess = BHC_Initial_Y(init_coefs)
  scale = BHC_Define_Y_Bounds(bhcParams.numPts, upr_bnd_coefs, lwr_bnd_coefs) 
  for i=1,N_Elements(initial_guess)-2 do initial_guess[i] = 0.5*(scale[i*2]+scale[i*2+1]) 
  
END

;------------------BHC_Inverse_Reporting_Setup-----------------------------
; PURPOSE:
;   Create a file to record progress of BHC inversion 
;   
; INPUTS:
;  num_pts:  number of parameters being iterated on
;  
; OUTPUTS:
;    Common block values for LUN, iter, and show set
;    
; SIDE EFFECTS:
;   - File "detailed_bhc_results" is opened with number LUN and header  
;     lines written
; ---------------------------------------------------------
Pro BHC_Inverse_Reporting_Setup, num_pts      
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
        LUN, show, y0, yN, bhcParams, matMeans, reconutct 

  x = *(bhcParams.x)
  
  OPENW, LUN,  new_dir + 'detailed_bhc_results.txt', /GET_LUN
  tab = string(9B)
  header = 'Raw File' + tab + "Merit Val"
  for i=0, num_pts-1 do header = header + tab + 'Coef #' + strtrim(string(i),2)
  for i=0, num_pts-1 do header = header + tab + 'y #' + strtrim(string(i),2)
  printf, LUN, header
  header = 'X coordinates' + tab
  for i=0, num_pts-1 do header = header + tab + strtrim(string(x[i]),2)
  for i=0, num_pts-1 do header = header + tab + strtrim(string(x[i]),2)
  printf, LUN, header
  
  iter = 0
  show = 1
End

;------------------BHC_Inverse_Reporting_Iterate-----------------------------
;
; PURPOSE:
;   Report a single iteration BHC inversion 
;   
; INPUTS:
;  merit:  Merit function value for this iteration
;  coefs:  BHC function coefficients for this iteration (y values)
;  slice:  Slice image for this iteration
;    
; OUTPUTS:
;    none
;    
; SIDE EFFECTS:
;   - One line written to file "detailed_bhc_results" 
;   - Window 0 created if necessary, and slice name shown in window
; ---------------------------------------------------------
Pro BHC_Inverse_Reporting_Iterate, merit, coefs, slice
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
        LUN, show, y0, yN, bhcParams, matMeans, reconutct 
 
  tab = string(9B)
  outln = reconParams.RawFile + tab + StrTrim(String(merit),2)
  for i=0, N_Elements(coefs)-1 do outln = outln + tab + strtrim(string(coefs[i]),2)
  y = BHC_Construct_Transform(coefs)
  for i=0, N_Elements(y)-1 do outln = outln + tab + strtrim(string(y[i]),2)
    
  printf, LUN, outln

  sz_s = Size(slice)
  if (show EQ 1) then begin
    Window, 0, TITLE = "BHC progress", XSIZE=sz_s[1], YSIZE=sz_s[2]
    show=2
  endif
  if (show EQ 2) then TVSCL, slice
  
End

;------------------BHC_Inverse_Reporting_Finalize-----------------------------
; PURPOSE:
;   Shuts down reportiong of BHC inversion 
;   
; INPUTS:
;  none
;    
; OUTPUTS:
;  none
;    
; SIDE EFFECTS:
;   - File LUN is closed
;   - Window 0 is deleted
; ---------------------------------------------------------
Pro BHC_Inverse_Reporting_Finalize, merit, coefs
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
       LUN, show, y0, yN, bhcParams, matMeans, reconutct 
        
  ; print final winning merit value and coefs
  tab = string(9B)
  outln = reconParams.RawFile + tab + StrTrim(String(merit),2)
  for i=0, N_Elements(coefs)-1 do outln = outln + tab + strtrim(string(coefs[i]),2)
  y = BHC_Construct_Transform(coefs) 
  for i=0, N_Elements(y)-1 do outln = outln + tab + strtrim(string(y[i]),2)
    
  printf, LUN, outln
        
  WDelete, 0
  close, LUN
  free_lun, LUN
End

;------------------BHC_Inverse_Reduce_Resolution-----------------------------
;
; PURPOSE:
;   Reduces resolution of image and ROI to speed up finding result
;
; INPUTS:
;  factor: Factor by which to reduce resolution.  For example, 2 reduces resolution by 1/2
;
; OUTPUTS:
;  none
;
; SIDE EFFECTS:
;   - reconParams.Matrix in common block reduced by <factor>
;   - roi.selected scaled down by <factor>
; ---------------------------------------------------------
Pro BHC_Inverse_Reduce_Resolution, factor
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  reconParams.Matrix = FIX(reconParams.Matrix/factor)
  objrefs = roi.group->Get(/ALL, COUNT=nRegions)
  for i=0, nRegions-1 do objrefs[i]->Scale, 1./factor, 1./factor
End


;------------------BHC_Inverse_Restore_Resolution-----------------------------
;
; PURPOSE:
;   Restore resolution of slice and ROI
;
; INPUTS:
;  factor: Factor by which resolution was reduced; use same number passed to BHC_Inverse_Reduce_Resolution
;
; OUTPUTS:
;  none
;
; SIDE EFFECTS:
;   - reconParams.Matrix in common block restored
;   - roi.selected scaled up by <factor>
; ---------------------------------------------------------
Pro BHC_Inverse_Restore_Resolution, factor
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
     LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  reconParams.Matrix = FIX(reconParams.Matrix * factor)
  objrefs = roi.group->Get(/ALL, COUNT=nRegions)
  for i=0, nRegions-1 do objrefs[i]->Scale, factor, factor
  show=1   ; used as a flag to signal resolution restored
End

;------------------BHC_Graph_BHC_transform-----------------------------
;
; PURPOSE:
;   Displays a new window with the final BHC transform function.   Also graphs the 
;       BHC functions used to the globally bound the merit function for comparison
;   
; INPUTS:
;   - bhc_Param structure
;     
; SIDE EFFECTS:
;   - Window 0 is opened
; ---------------------------------------------------------
Pro BHC_Graph_BHC_transform, in_bhcParams

  if (N_Params() EQ 1) then bhcParams = in_bhcParams

  x = *(bhcParams.x)
  
  sz_x = size(x)
  
  upr_bnd_coefs = [0.0,0.98,0.0,0.35]
  lwr_bnd_coefs = [0.0,0.01,0.0,0.0]
  num_pts = bhcParams.numPts
  scale = BHC_define_y_bounds(num_pts, upr_bnd_coefs, lwr_bnd_coefs)
  final_y = *(bhcParams.bhc_coefs)
  
  new_scale = [x[0],x[0],scale,x[num_pts-1],x[num_pts-1]]
   
  ; plot upper bounds of new_scale to get yrange
  y = findgen(num_pts)
  for i=0,num_pts-1 do y[i] = new_scale[i*2+1]
  y = BHC_Construct_Transform(y)
  white = [255,255,255]
  black = [0,0,0]
  iplot, x, y,  XRANGE=[x[0], x[sz_x[1]-1]], YRANGE=[y[0], (y[sz_x[1]-1])+0.2], COLOR=black, $
      WINDOW_TITLE = 'BHC Transformation Function', DIMENSIONS=[900,600], BACKGROUND_COLOR=white, $
      FONT_COLOR=black, XTEXT_COLOR=black, YTEXT_COLOR=black, /ZOOM_ON_RESIZE, IDENTIFIER=ip, /NO_SAVEPROMPT, $
      XTITLE='Original Data Value (-ln(I/Io))', YTITLE='Corrected Data Value (-ln(I/Io))', TITLE='BHC Data Transformation'
      
  ; have to change axes color programatically 
  oSys = _IDLitSys_GetSystem()
  oTool = oSys->GetByIdentifier(ip)
  void = oTool->DoSetProperty('WINDOW', 'AUTO_RESIZE', '1')
  oTool->CommitActions
  
  ; just plot with original x for now; white dotted = original data
  iplot, x, x ,/OVERPLOT, LINESTYLE=1, COLOR=black
  
  y = findgen(num_pts)
  for i=0,num_pts-1 do y[i] = new_scale[i*2]
  y = BHC_Construct_Transform(y) ; RAK 2/8/2010
  iplot, x, y , /OVERPLOT, COLOR=black  ; black - lwr bound
 
  ; Get function 
  y_plot = BHC_Construct_Transform(final_y) ; RAK 2/8/2010
 
  ; Make spline version
  y2 = spl_init(x,y_plot, /DOUBLE)
  xdense = Indgen(101)*0.01*(x[num_pts-1]-x[0])+x[0]
  ydense = spl_interp(x,y_plot,y2,xdense, /DOUBLE)
  iplot, xdense, ydense, /OVERPLOT, COLOR=[255,0,0]

  ; Plot nodes
  iplot, x, y_plot, /OVERPLOT, COLOR=[255,0,0], LINESTYLE=6, SYM_INDEX=4  ; red - final transform
  
  ; deselect all selected items so graph doesn't look odd when user first sees it
  selItems = oTool->GetSelectedItems(COUNT=num)
  for i=0, num-1, 1 do selItems[i]->Select,0
  
End

;;------------------BHC_FindCorrection_Error-----------------------------
;
; PURPOSE:
;       Error function for FindCorrection function.  Allows smooth recovery
;       (restores downsample, if any) on error.
;
; INPUTS:
;     downsample: current downsample setting
;
; OUTPUTS:
;     none
;
; ---------------------------------------------------------
PRO BHC_FindCorrection_Error, downsample

  COMMON DO_BHC_ARGS, last_raw_f_opened, scaled_data, head1, head2, sz, bits_per_sample
    
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
  
  ; present error
  print, 'Error generated at ' + !ERROR_STATE.BLOCK
  msg = (!ERROR_STATE.MSG NE '')? !ERROR_STATE.MSG : 'null'
  sys_msg = (!ERROR_STATE.SYS_MSG NE '')? !ERROR_STATE.SYS_MSG : 'null'
  full_error_msg = 'FIND_BHC detected an error during ''FIND BHC for Slice''.'
  if (msg NE 'null') then full_error_msg= [[full_error_msg], '','Error msg: ', msg]
  if (sys_msg NE 'null') then full_error_msg= [[full_error_msg], '','System error: ', sys_msg]
  full_error_msg= [[full_error_msg], '','BHC could not be found.']
  
  result = dialog_message(full_error_msg, /Error)
  
  ; Restore old ROIs, resolution
  if (n_elements(downsample) GT 0) then if (downsample GT 1 AND SHOW EQ 2) then BHC_Inverse_Restore_Resolution, downsample
  
  ; reset file read flags so reread for each subsequent run
  last_raw_f_opened = ''
  scaled_data = 0

END

;------------------Structure Definitions-----------------------------
PRO BHC_Params__Define
struct = {BHC_Params, $
   bhc_coefs: Ptr_New(), $    ; BHC coefficients that represent the BHC transform function
   merit_value: 0.0, $        ; latest merit function value (BHC correction value - lower number is better)
   merit_valuePrev: 0.0, $    ; previous merit function value
   numPts: 0, $               ; number of coefficients to use in function fit
   downsample: 0, $           ; downsample data for finding BHC
   pos_2nd_der: 0, $          ; enforce positive second derivative at all points
   roiMode: 0, $              ; roi mode for BHC
   x: Ptr_New(), $            ; X (CT Data) range for BHC transform function
   reconstruct: 0, $          ; reconstruct raw file(s)
   mod_scale: 0.0, $          ; modified image scale after BHC, calculated in BHC_Do_Recon
   mod_offset: 0.0 $          ; modified image offset after BHC, calculated in BHC_Do_Recon
  }
END


Pro BHC_ROI_Params__Define
  struct = {BHC_ROI_Params, $
    group : Obj_New(), $      ; Container for region-of-interest objects (IDLgrROI)
    selected : Obj_New() $    ; Currently selected region (null object if none selected)
    }  
End


PRO BHC_Reconutct_Params__Define
  temp = {BHC_Reconutct_Params, $   ; A list of reconstruction parameters required by the 
                                    ;   command line BIR reconstruction program.  
    RawFile:'', $
    ImageFile:'', $
    Quality:'NA', $
    Matrix:'NA', $
    DoRfc:'NA', $
    RfcSlcWd:'NA',  $
    RfcSect:'NA', $
    RfcRngAmp:'NA', $
    RfcRngMlt:'NA', $
    RfcRngAdj:'NA', $
    ImgScl:'NA',  $
    ImgOff:'NA',  $
    ImgRot:'NA',  $
    ImgCrayAdj:'NA',  $
    ImgSOD:'NA',  $
    ImgFOV:'NA',  $
    ImgCentX:'NA',  $
    ImgCentY:'NA',  $
    FanAngle:'NA',  $
    MedFiltIter:'NA', $
    ConvFiltFile:'NA',  $
    FiltOrder:'NA', $
    RmvRaw:'NA',  $
    WedgeFile:'NA', $
    BHCCoefCnt:'NA',  $
    BHCVal1:'NA', $
    BHCVal2:'NA', $
    BHCVal3:'NA', $
    BHCVal4:'NA', $
    BHCVal5:'NA', $
    BHCVal6:'NA', $
    BHCVal7:'NA', $
    BHCVal8:'NA', $
    BHCVal9:'NA', $
    BHCVal10:'NA',  $
    NumViews:'NA',  $
    SkipScale:'NA', $
    NoOverwrite:'NA',  $
    Reserved1:'NA',  $
    Reserved2:'NA'}
END

FUNCTION BHC_Init_Params, params
    params.Quality = 'NA'
    params.Matrix = 'NA'
    params.DoRfc = 'NA'
    params.RfcSlcWd = 'NA'
    params.RfcSect = 'NA'
    params.RfcRngAmp='NA'
    params.RfcRngMlt='NA'
    params.RfcRngAdj='NA'
    params.ImgScl='NA'
    params.ImgOff='NA'
    params.ImgRot='NA'
    params.ImgCrayAdj='NA'
    params.ImgSOD='NA'
    params.ImgFOV='NA'
    params.ImgCentX='NA'
    params.ImgCentY='NA'
    params.FanAngle='NA'
    params.MedFiltIter='NA'
    params.ConvFiltFile='NA'
    params.FiltOrder='NA'
    params.RmvRaw='NA'
    params.WedgeFile='NA'
    params.BHCCoefCnt='NA'
    params.BHCVal1='NA'
    params.BHCVal2='NA'
    params.BHCVal3='NA'
    params.BHCVal4='NA'
    params.BHCVal5='NA'
    params.BHCVal6='NA'
    params.BHCVal7='NA'
    params.BHCVal8='NA'
    params.BHCVal9='NA'
    params.BHCVal10='NA'
    params.NumViews='NA'
    params.SkipScale='NA'
    params.NoOverwrite='NA'
    params.Reserved1='NA'
    params.Reserved2='NA'
    
    params.RawFile = 'NA'
    params.ImageFile = 'NA'
    
    return, params
END


;------------------BHC_FindCorrection-----------------------------
;
; PURPOSE:
;   Finds the optimum beam hardening correction based on user-defined ROI(s) 
;   within an image.  Uses an optimization routine (IDL_AMOEBA) to find the
;   best correction based on minimizing the standard deviation of the given
;   return (BHC_ROI_Merit_Func()).
;
; INPUTS:
;  raw_f_in: (string) path to raw file (sinogram file) to correct and reconstruct
;  
;  roi_in:   (pointer to BHC_ROI_Params) pointer pass from viewer (GUI) which 
;         holds all ROI information.  This structure has been greatly simplied from it's
;         original version since most members are not used in this example code.
;     
;  reconParams_in: (Pointer to Reconutct_Params struct) - A list of reconstruction
;                   parameters required by the command line reconstruction program.  
;  
;  bhcParams_in:  (pointer to BHC_Params struct) pointer pass from viewer (GUI) which holds
;                 all of the beam hardening correction options
;                 
;  stack_min: (float or int) Assuming that the input image (sinogram) is one of a stack of 2D slices that make up a full 
;             3D volume, this is the minimum of the entire stack so a BHC that applies to the entire stack
;             data space can be found.
;             
;  stack_max: (float or int) Assuming that the input image (sinogram) is one of a stack of 2D slices that make up a full 
;             3D volume, this is the maximum of the entire stack so a BHC that applies to the entire stack
;             data space can be found.
;             
;  reconutct_in: (string) path to the command line reconstruction program
;
;  orig_image: (2D numeric array) image reconstructed without a beam hardening correction.
;  
; OUTPUTS:
;   Floating point array of beam hardening coefficients that were found to provide the optimal correction
;
; SIDE EFFECTS:
;   - an IDLBHC directory is created and filled with several raw (sinogram) and reconstructed (image) files
;   
;  COMMON BLOCK ROI_MERIT_FUNC_ARGS: 
;   Variables required for the IDL amoeba() function.  The IDL Amoeba function is 
;     structured as to only allow a single variable to be passed to the merit function.  So here we pass the rest 
;     of the required variables via a common block.  However, the Amoeba function is written in IDL and could be 
;     easily modified to pass other arguments to the function (called with the IDL call_function, which does
;     allow an array of arguments to be passed).  
;
;  COMMON BLOCK DO_BHC_ARGS: 
;   All of these variables are used in an effort make I/O more efficient and
;     speed up the execution of the merit function (called up to 500 times with amoeba()).  
;     BHC_ReadBIRSinogramSmart function in particular uses several of these
;     persistent variables to speed up file I/O.  
; ---------------------------------------------------------

FUNCTION BHC_FindCorrection, raw_f_in, roi_in, reconParams_in, bhcParams_in, stack_min, stack_max, $
            reconutct_in, orig_image
  
  ; This is a common block for variables shared among support functions for BHC_Correction to reduce 
  ;     the number of variables passed in and out of functions.
  COMMON ROI_MERIT_FUNC_ARGS, raw_f, roi, new_dir, fn_base, reconParams, iter, $
    LUN, show, y0, yN, bhcParams, matMeans, reconutct
    
  COMMON DO_BHC_ARGS, last_raw_f_opened, scaled_data, head1, head2, sz, bits_per_sample
 
  ; begin timer
  start_time = systime(/SECONDS)
    
  ; Define error handler (needed to recover downsampling of image)
  Catch, error
  if(error NE 0) then begin
    Catch, /CANCEL
    BHC_FindCorrection_Error, downsample
    return, -1
  endif
  
  ; set COMMON variables (ROI_MERIT_FUNC_ARGS) for merit function
  bhcParams = (*bhcParams_in)
  raw_f = raw_f_in
  roi = (*roi_in)
  reconParams = (*reconParams_in)
  reconutct = reconutct_in
  iter = 0
  
  ; set local variables
  num_pts = bhcParams.numPts
  
  ; create new IDLBHC dir and delete any files that might already be in directory
  dir = Path_Sep()
  divide = strpos(raw_f, dir, /REVERSE_SEARCH) + 1
  new_dir = strmid(raw_f, 0, divide) + 'IDLBHC_AUTO' + dir
  FILE_MKDIR, new_dir, /NOEXPAND_PATH
  
  ; create data folder and delete any existing files in it
  fn = strmid(raw_f, divide)               ; get filename without directory structure
  fn_base = strmid(fn, 0, strlen(fn)-4)    ; strip -Raw from end
  new_dir = new_dir + fn_base + dir
  FILE_MKDIR, new_dir, /NOEXPAND_PATH
  old_files = file_search(new_dir + '*', COUNT=num_old_files)
  if (num_old_files gt 0) then file_delete, old_files, /QUIET
  
  ; get initial_guess and scale
  BHC_Define_P0_and_Scale, stack_min, stack_max, initial_guess, scale
  scale_save = scale
  
  ; Set up reporting of iteration
  BHC_Inverse_Reporting_Setup, num_pts
  
  ; Initialize means
  dummy = BHC_ROI_Merit_Func(initial_guess,INIT_IMAGE=orig_image)
  
  ; Reduce reconstruction size if requested
  downsample = bhcParams.downsample
  if (downsample GT 1) then BHC_Inverse_Reduce_Resolution, downsample
  
  ; call optimization function to begin.
  ; nmax is max possible function calls (reconstructions). Default is 5000
  nmax = 500
  ; ftol is function tolerance used to stop optimization. 
  ;   should never be less than sqrt(fp_prec) (function should return floating point result)
  prec = MACHAR()
  ftol = sqrt(prec.EPS)*2.^(downsample) 
  ; We are calling a custom version of the Amoeba function, which has been minorly modified to 
  ;     accomodate for enforcing a positive second derivative, if that option is selected, and for
  ;     setting upper and lower bounds on the coefficients.
  final_y = BHC_AMOEBA(ftol, FUNCTION_NAME='BHC_ROI_Merit_Func', FUNCTION_VALUE=fvalue, NCALLS=ncalls, $
    NMAX=nmax, P0=initial_guess, SCALE=scale)
  
  if (n_elements(final_y) eq 0) OR (n_elements(final_y) eq 1) then begin
    BHC_Inverse_Reporting_Finalize, 0, 0
    message, 'Could not converge on a suitable BHC after ' + strtrim(string(ncalls),2)+ $
              'iterations.  Please adjust the parameters and/or ROIs and try again.'
  endif
  
  ;  Restore old ROIs, resolution
  if (downsample GT 1) then BHC_Inverse_Restore_Resolution, downsample
  
  ; ends fixed, append y0 and yN to answer returned by BHC_amoeba
  final_y = [y0,final_y,yN]
  
  ; Generate full-resolution slice
  BHC_Do_Recon, final_y, /RESCALE
  
  BHC_ReadBIRSinogramSmart, reconParams.ImageFile, slice, fhead1, fhead2, fsz
  sz = size(slice)
  
  ; update reconstruction parameters
  reconParams.SkipScale = 'NA'
  reconParams.RawFile = 'NA'
  reconParams.ImageFile = 'NA'
  (*reconParams_in) = reconParams
  
  ; update bhc_params
  bhcParams.bhc_coefs = Ptr_new(final_y)
  bhcParams.merit_value = fvalue[0]
  (*bhcParams_in) = bhcParams
  
  ;Shut down iterative reporting and display BHC transform function and corrected image
  BHC_Inverse_Reporting_Finalize, fvalue[0], final_y
  BHC_Graph_BHC_transform, bhcParams
  total_min = ceil((systime(/SECONDS) - start_time)/60)  
  print, "Number of Iterations: " + strtrim(string(iter),2)
  print, "Time to find solution (minutes): " + strtrim(string(total_min),2)
  Window, 0, TITLE = "Final Corrected Image", XSIZE=reconParams.matrix, YSIZE=reconParams.matrix
  TVSCL, slice
  
  
  ; reset file read flags
  last_raw_f_opened = ''
  scaled_data = 0
  
  return, final_y
END
