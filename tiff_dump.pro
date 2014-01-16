; $Id: //depot/Release/IDL_81/idl/idldir/lib/obsolete/tiff_dump.pro#1 $
;
; Copyright (c) 1991-2011. ITT Visual Information Solutions. All
;       rights reserved. Unauthorized reproduction is prohibited.
;
PRO TiffMakeKeyValue, defs, TagIndex, TagName
; Given a TagValue/TagName table in the form of a string array, return
; two pointers to arrays containing the Tag Indices and the Tag Names,
; respectively.
;
names = defs                    ;Separate tag names....
tags = uint(defs)               ;Extract tag indices
for i=0, n_elements(defs)-1 do begin
    blank = strpos(defs[i], ' ')
    names[i] = strmid(defs[i], blank+1, 100)
endfor
order = sort(tags)              ;Sort by tag value
TagIndex = ptr_new(tags[order])
TagName = ptr_new(names[order])
end

function tiff_slong, a, i, len=len ;return singed longword(s) from array a(i)
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs

if n_elements(len) le 0 then len = 1
if len gt 1 then result = long(a,i,len) $
else result = long(a,i)
if order then byteorder, result, /lswap
return, result
end

function tiff_ulong,a,i,len=len	;return longword(s) from array a(i)
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs

if n_elements(len) le 0 then len = 1
if len gt 1 then result = ulong(a,i,len) $
else result = ulong(a,i)
if order then byteorder, result, /lswap
return, result
end


function tiff_urational,a,i, len = len ; return unsigned rational from a(i)
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs

if n_elements(len) le 0 then len = 1
tmp = tiff_ulong(a, i, len = 2 * len) ;1st, cvt to ulongwords
if len gt 1 then begin
    subs = lindgen(len)
    rslt = float(tmp[subs*2]) / tmp[subs*2+1]
endif else rslt = float(tmp[0]) / tmp[1]
return, rslt
end

function tiff_srational, a, i, len = len ; return signed rational from a(i)
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs

if n_elements(len) le 0 then len = 1
tmp = tiff_slong(a, i, len = 2 * len)	;1st, cvt to longwords
if len gt 1 then begin
	subs = lindgen(len)
	rslt = float(tmp[subs*2]) / tmp[subs*2+1]
endif else rslt = float(tmp[0]) / tmp[1]
return, rslt
end

function tiff_sint,a,i, len=len	;return signed short int from TIFF short int
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs
if n_elements(len) le 0 then len = 1
if len gt 1 then begin	;Array?
	result = fix(a,i,len)
endif else begin	;Scalar
	result = fix(a,i)
endelse
if order then byteorder, result, /sswap
return, result
end


function tiff_uint,a,i, len=len	;return unsigned short from TIFF short int
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs
if n_elements(len) le 0 then len = 1
if len gt 1 then begin	;Array?
	result = uint(a,i,len)
endif else begin	;Scalar
	result = uint(a,i)
endelse
if order then byteorder, result, /sswap
return, result
end


function tiff_sbyte, a,i,len=len ;return sbytes from array a(i)
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs

if n_elements(len) le 0 then len = 1
if len gt 1 then begin
    result = fix(a[i:i+len-1])
    tmp = where(result ge 128, count)
    if count gt 0 then result[tmp] = 256-result[tmp]
endif else begin
    result = fix(a[i])
    if result ge 128 then result = 256-result
endelse

return, result
end


function tiff_float, a, i, len=len
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs
if n_elements(len) le 0 then len = 1
if len gt 1 then result = float(a,i,len) $
else result = float(a,i)
if order then byteorder, result, /lswap
return, result
end

function tiff_double, a, i, len=len
common tiff_com, order, ifd, count

on_error,2              ;Return to caller if an error occurs
if n_elements(len) le 0 then len = 1
if len gt 1 then result = double(a,i,len) $
else result = double(a,i)
if order then byteorder, result, /l64swap
return, result
end



PRO define_tiff_tags, Dummy ;Define the TIFF baseline tag values & names.
; Save values in common.
common tiff_dump_com, TagIndex, TagName, TypeName, TypeLen ;Database...
defs = [ $
         '254 NewSubfileType', $
         '255 SubfileType', $
         '256 ImageWidth', $
         '257 ImageLength', $
         '258 BitsPerSample', $
         '259 Compression', $
         '262 PhotometricInterpretation', $
         '263 Thresholding', $
         '264 CellWidth', $
         '265 CellLength', $
         '266 FillOrder', $
         '269 DocumentName', $
         '270 ImageDescription', $
         '271 Make', $
         '272 Model', $
         '273 StripOffsets', $
         '274 Orientation', $
         '277 SamplesPerPixel', $
         '278 RowsPerStrip', $
         '279 StripByteCounts', $
         '280 MinSampleValue', $
         '281 MaxSampleValue', $
         '282 XResolution', $
         '283 YResolution', $
         '284 PlanarConfiguration', $
         '285 PageName', $
         '286 XPosition', $
         '287 YPosition', $
         '288 FreeOffsets', $
         '289 FreeByteCounts', $
         '290 GrayResponseUnit', $
         '291 GrayResponseCurve', $
         '292 T4Options', $
         '293 T6Options', $
         '296 ResolutionUnit', $
         '297 PageNumber', $
         '301 TransferFunction', $
         '305 Software', $
         '306 DateTime', $
         '315 Artist', $
         '316 HostComputer', $
         '317 Predictor', $
         '318 WhitePoint', $
         '319 PrimaryChromaticities', $
         '320 ColorMap', $
         '321 HalftoneHints', $
         '322 TileWidth', $
         '323 TileLength', $
         '324 TileOffsets', $
         '325 TileByteCounts', $
         '330 SubIFDs', $
         '338 ExtraSamples', $
         '339 SampleFormat', $
         '340 SMinSampleValue', $
         '341 SMaxSampleValue ', $
         '342 TransferRange', $
         '343 ClipPath', $
         '344 XClipPathUnits', $
         '345 YClipPathUnits', $
         '346 Indexed', $
         '347 JPEGTables', $
         '529 YCbCrCoefficients', $
         '530 YCbCrSubSampling', $
         '531 YCbCrPositioning', $
         '532 ReferenceBlackWhite']
; Note: we ignore tags 512-521, for the old JPEG format which are
; obsoleted by TIFF Technical Note #2 which specifies a revised
; JPEG-in-TIFF scheme. 

TiffMakeKeyValue, defs, TagIndex, TagName
TypeName = [ 'Undefined', 'Byte', 'Ascii', 'UShort', 'ULong', 'URatnl', $
             'SByte', 'Undefined', 'SShort', 'Slong', 'SRatnl', $
             'Float', 'Double' ]
TypeLen = [0, 1, 1, 2, 4, 8, 1, 1, 2, 4, 8, 4, 8] ;lengths of tiff types
end





function tiff_get_field, index, lun
                                ;Return contents of field index

common tiff_com, order, ifd, count
common tiff_dump_com, TagIndex, TagName, TypeName, TypeLen ;Database...

on_error,2                      ;Return to caller if an error occurs
if n_elements(TagIndex) eq 0 then begin ;Load tables?
    define_tiff_tags
endif

ent = ifd[index * 12: index * 12 + 11] ;Extract the ifd
tag = tiff_uint(ent, 0)
typ = tiff_uint(ent, 2)
if typ gt 12 or typ lt 1 then typ = 7 ;Unknown type if out of range
tname = TypeName[typ]
cnt = tiff_ulong(ent, 4)
;if (typ le 0) or (typ gt 5) then $
;  message,'Illegal type code, ifd = '+string(index)

nbytes = cnt * TypeLen[typ]
IF (nbytes GT 4) THEN BEGIN 	;value size > 4 bytes ?
    offset = tiff_ulong(ent, 8) ;field has offset to value location
    Point_Lun, lun, offset
    val = BytArr(nbytes, /NOZERO) ;buffer will hold value(s)
    Readu, lun, val
endif else begin
    val = ent[8:*]              ;Rest of value field
endelse

CASE typ OF                     ;Ignore bytes, as there is nothing to do
    1: val = fix(tiff_byte(val, 0, len=cnt)) ;For printing
    2: begin                    ;tiff ascii type
        val = String(val) & cnt = 1
    Endcase
    3: val = tiff_uint(val,0, len = cnt)
    4: val = tiff_ulong(val,0, len = cnt)
    5: val = tiff_urational(val,0, len = cnt)
    6: val = tiff_sbyte(val, 0, len=cnt)
    7: begin val = "<Undefined>" & cnt = 1  & endcase
    8: val = tiff_sint(val, 0, len=cnt)
    9: val = tiff_slong(val, 0, len=cnt)
    10: val = tiff_srational(val, 0, len=cnt)
    11: val = tiff_float(val, 0, len=cnt)
    12: val = tiff_double(val, 0, len=cnt)
ENDCASE
return, val
end

pro tiff_dump_field, index, lun, PRIVATE=priv, TEXT=txt
                                ;Return contents of field index

common tiff_com, order, ifd, count
common tiff_dump_com, TagIndex, TagName, TypeName, TypeLen ;Database...

on_error,2                      ;Return to caller if an error occurs
if n_elements(TagIndex) eq 0 then begin ;Load tables?
    define_tiff_tags
endif

ent = ifd[index * 12: index * 12 + 11] ;Extract the ifd
tag = tiff_uint(ent, 0)
typ = tiff_uint(ent, 2)
if typ gt 12 or typ lt 1 then typ = 7 ;Unknown type if out of range
tname = TypeName[typ]
cnt = tiff_ulong(ent, 4)
if tag ge 32768L then begin     ;Private tag?
    name = '<PrivateTag>'
    if keyword_set(priv) then begin ;Search private tags?
        i = where(tag eq *Priv[0], j)
        if j gt 0 then name = (*Priv[1])[i[0]]
    endif
endif else begin
  i = where(tag eq *TagIndex, j) ;Look up name...
  if j gt 0 then name = (*TagName)[i[0]] else name = '<NoName>'
endelse

;if (typ le 0) or (typ gt 5) then $
; message,'Illegal type code, ifd = '+string(index)

line = string(name, FORMAT='(a, t24)') ;Format hdr
line = line + ' (' + strtrim(tag,2) + ')' + tname

nbytes = cnt * TypeLen[typ]
IF (nbytes GT 4) THEN BEGIN 	;value size > 4 bytes ?
    offset = tiff_ulong(ent, 8) ;field has offset to value location
    Point_Lun, lun, offset
    val = BytArr(nbytes, /NOZERO) ;buffer will hold value(s)
    Readu, lun, val
endif else begin
    val = ent[8:*]              ;Rest of value field
endelse

CASE typ OF                     ;Ignore bytes, as there is nothing to do
    1: val = fix(tiff_byte(val, 0, len=cnt)) ;For printing
    2: begin                    ;tiff ascii type
        val = String(val) & cnt = 1
    Endcase
    3: val = tiff_uint(val,0, len = cnt)
    4: val = tiff_ulong(val,0, len = cnt)
    5: val = tiff_urational(val,0, len = cnt)
    6: val = tiff_sbyte(val, 0, len=cnt)
    7: val = fix(tiff_byte(val, 0, len=cnt)) ;For printing
    8: val = tiff_sint(val, 0, len=cnt)
    9: val = tiff_slong(val, 0, len=cnt)
    10: val = tiff_srational(val, 0, len=cnt)
    11: val = tiff_float(val, 0, len=cnt)
    12: val = tiff_double(val, 0, len=cnt)
    else: begin val = "<UnknownType>" & cnt = 1  & endcase
ENDCASE
if cnt gt 1 then begin
    a = ''
    for i=0, cnt-1 < 15 do a = a + ' '  + strtrim(val[i],2)
    line = line +'<Cnt='+strtrim(cnt,2)+'> = '+ a
endif else line = line + ' = '+ strtrim(val, 2)
if arg_present(txt) then txt = line ; else print,line
end


Pro TIFF_READ_DIRECTORY, file, tags, names, values, PRIVATE=priv, $
        XTAGS=xtags, DIRECTORY_OFFSET=doffs
common tiff_com, order, ifd, count
common tiff_dump_com, TagIndex, TagName, TypeName, TypeLen ;Database...

if n_elements(TagIndex) eq 0 then begin ;Load tables?
    define_tiff_tags
endif

openr,lun,file, error = i, /GET_LUN, /BLOCK

if i lt 0 then begin ;OK?
    message, 'Error opening ' + file
endif

hdr = bytarr(8)			;Read the header
readu, lun, hdr

typ = string(hdr[0:1])		;Either MM or II
if (typ ne 'MM') and (typ ne 'II') then begin
    print,'TIFF_READ: File is not a Tiff file: ', file
    return
endif

endian = byte(1,0,2)		;What endian is this?
endian = endian[0] eq 0		;1 for big endian, 0 for little
order = (typ eq 'MM') xor endian ;1 to swap...

if n_elements(doffs) eq 0 then begin
    offs = tiff_ulong(hdr, 4)	;Get Offset to IFD
endif else offs = doffs

point_lun, lun, offs            ;Read it
a = bytarr(2)			;Entry count array
readu, lun, a
count = tiff_uint(a,0)		;count of entries
if count le 0 or count gt 1000 then $
  message, 'Invalid directory count field'

ifd = bytarr(count * 12)	;Array for IFD's
readu,lun, ifd			;read it

old_tag = 0			;Prev tag...
tags = uintarr(count)
names = strarr(count)           ;Get names & values
values = ptrarr(count)

if n_elements(xtags) eq 2 then begin ;Use separate tag defs?
    TTagIndex = xtags[0]
    TTagName = xtags[1]
endif else begin
    TTagIndex = TagIndex
    TTagName = TagName
endelse

for index = 0L,count-1 do begin	;Print each directory entry
    tag = tiff_uint(ifd, index*12)
    tags[index] = tag
    ent = ifd[index * 12: index * 12 + 11] ;Extract the ifd
    tag = tiff_uint(ent, 0)
    typ = tiff_uint(ent, 2)
    if typ gt 12 or typ lt 1 then typ = 7 ;Unknown type if out of range
    cnt = tiff_ulong(ent, 4)    ;Count field
    name = ''                   ;Get tag's name

    if tag ge 32768L and keyword_set(priv) then begin ;Private tag?
        i = where(tag eq *Priv[0], j)
        if j gt 0 then name = (*Priv[1])[i[0]]
    endif else begin
        i = where(tag eq *TTagIndex, j) ;Look up name...
        if j gt 0 then name = (*TTagName)[i[0]]
    endelse                     ;Private tag

    if strlen(name) eq 0 then begin ;Unknown tag...
        name = 'Unk_'+strtrim(tag,2)
    endif

    nbytes = cnt * TypeLen[typ]
    IF (nbytes GT 4) THEN BEGIN ;value size > 4 bytes ?
        offset = tiff_ulong(ent, 8) ;field has offset to value location
        Point_Lun, lun, offset
        val = BytArr(nbytes, /NOZERO) ;buffer will hold value(s)
        Readu, lun, val
    endif else begin
        val = ent[8:*]          ;Rest of value field
    endelse

    CASE typ OF                 ;Ignore bytes, as there is nothing to do
        1: dummy = 0            ;Already in bytes.
        2: begin                ;tiff ascii type
            val = String(val) & cnt = 1
        Endcase
        3: val = tiff_uint(val,0, len = cnt)
        4: val = tiff_ulong(val,0, len = cnt)
        5: val = tiff_urational(val,0, len = cnt)
        6: val = tiff_sbyte(val, 0, len=cnt)
        7: begin val = "<Undefined>" & cnt = 1  & endcase
        8: val = tiff_sint(val, 0, len=cnt)
        9: val = tiff_slong(val, 0, len=cnt)
        10: val = tiff_srational(val, 0, len=cnt)
        11: val = tiff_float(val, 0, len=cnt)
        12: val = tiff_double(val, 0, len=cnt)
    ENDCASE
    names[index] = name
    values[index] = ptr_new(val)
endfor

free_lun, lun
lun = -1
end





pro tiff_dump, file, PRIVATE=priv, TEXT=text, GET_TAG_VALUE=gettag, $
        TAG_VALUE=tagv
;+
; NAME:
;	TIFF_DUMP
;
; PURPOSE:
;	Dump the Image File Directories of a TIFF File to the terminal
;	or into a string array.  This procedure is used mainly for debugging.
;
; CATEGORY:
;	Input/output.
;
; CALLING SEQUENCE:
;	TIFF_DUMP, Filename
;
; INPUTS:
;    Filename:	string containing the name of file to read.
;		The default extension is ".TIF".
;
; OUTPUTS:
;	By default, all output is to the terminal.  Each TIFF Image
;	File Directory entry is printed.
;
; KEYWORDS:
;   PRIVATE: If set, contains a vector of two pointers.  The first
;     pointer points to a vector containing the
;     values of the private keys of interest.  Private keys are those whose
;     directory key value is larger than 32767, and which have a special
;     meaning outside the TIFF standard.  The second pointer points to a
;     string vector, containing the Names of the private keys, described
;     in the first pointer.
;   TEXT: If set, the output of this routine will be placed in a
;     string array, returned in this parameter, rather than on the console.
;   GET_TAG_VALUE: If set, an array containing directory key
;     values for the fields whose value is to be returned in the
;     keyword parameter TAG_VALUE.
;   TAG_VALUE: If GET_TAG_VALUE is set, this keyword parameter will
;     return a vector of pointers, corresponding to the keys in
;     GET_TAG_VALUE.  If a directory entry in the file equal to an
;     entry in GET_TAG_VALUE exists, a pointer to its value is returned in the
;     corresponding element in TAG_VALUE.
; COMMON BLOCKS:
;	TIFF_COM.  Only for internal use.
;
; SIDE EFFECTS:
;	A file is read.
;
; RESTRICTIONS:
;	Not all of the tags have names encoded.
;	In particular, Facsimile, Document Storage and Retrieval,
;	and most no-longer recommended fields are not encoded.
;
; PROCEDURE:
;	The TIFF file Header and the IFD (Image File Directory) are read
;	and listed.
;
; MODIFICATION HISTORY:
;	DMS, Apr, 1991.  Extracted from TIFF_READ.
;	DMS, Dec, 1994.	 Added private tags
;	DMS, Jul, 1999.  Added keywords PRIVATE, TEXT, GET_TAG_VALUE,
;		         and TAG_VALUE.  Also added new unsigned data types.
;
;-

common tiff_com, order, ifd, count

on_error,2                      ;Return to caller if an error occurs

if arg_present(text) then begin
    text = strarr(1000)
    textpntr = 0L
endif

openr,lun,file, error = i, /GET_LUN, /BLOCK

if i lt 0 then begin ;OK?
    message, 'Error opening ' + file
endif

hdr = bytarr(8)			;Read the header
readu, lun, hdr

typ = string(hdr[0:1])		;Either MM or II
if (typ ne 'MM') and (typ ne 'II') then begin
	print,'TIFF_READ: File is not a Tiff file: ', file
	return
	endif
order = typ eq 'MM'  		;1 if Motorola 0 if Intel (LSB first or vax)
endian = byte(1,0,2)		;What endian is this?
endian = endian[0] eq 0		;1 for big endian, 0 for little
order = order xor endian	;1 to swap...

if arg_present(text) then begin
    text[textpntr] = string('Tiff File: byte order=',typ, $
                            ',  Version = ', tiff_uint(hdr,2))
    textpntr = textpntr + 1
endif else print,'Tiff File: byte order=',typ, ',  Version = ', tiff_uint(hdr,2)

offs = tiff_ulong(hdr, 4)	;Offset to IFD
nFid = 0L

while offs ne 0 do begin        ;Process each ifd
    point_lun, lun, offs        ;Read it
    a = bytarr(2)               ;Entry count array
    readu, lun, a
    count = tiff_uint(a,0)      ;count of entries
    print, '*** IFD', nFid, ' starting at offset:', offs
    nFid = nFid + 1
    if arg_present(text) then begin
        text[textpntr] = string(count, ' directory entries')
        textpntr = textpntr+1
    endif else print,count, ' directory entries'
    ifd = bytarr(count * 12 + 4) ;Array for IFD's + Offset of next IFD
    readu,lun, ifd              ;read it
    offs = tiff_ulong(ifd, count*12)

    if n_elements(gettag) gt 0 then begin
        tagv = ptrarr(n_elements(gettag))
    endif
    
    old_tag = 0			;Prev tag...
    for i=0, count-1 do begin	;Print each directory entry
        tag = tiff_uint(ifd, i*12)
        if tag lt old_tag then $
          print,'*** Error in TIFF directory: entries out of order ****'
        old_tag = tag
        if n_elements(gettag) ne 0 then begin
            j = where(tag eq gettag, count)
            if count gt 0 then tagv[j] = ptr_new(tiff_get_field(i,lun))
        endif
        
        if arg_present(text) then begin
            tiff_dump_field, i, lun, PRIVATE=priv, TEXT=str
            text[textpntr] = str
            textpntr = textpntr+1
        endif else tiff_dump_field, i, lun, PRIVATE=priv
    endfor
endwhile

free_lun, lun
if arg_present(text) then text = text[0:textpntr-1]
lun = -1
end
