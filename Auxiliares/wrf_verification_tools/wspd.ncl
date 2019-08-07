;*****************************************************************************************
;
; Wind Speed Plotting (wspd.ncl)
; Created by Jason W. Godwin (jgodwin@rsmas.miami.edu)
; Date modified: 18 Jul 2013
;
; Description: This program takes in wrfout files and plots 10-meter wind speeds in knots
;
;*****************************************************************************************

;-------- Edit this block only! -----------------------------------------

fig_title	= "19 km WRF - Irene (2011)"

;-------------------------------------------------------------------------

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

; Define colormap

cmap = (/(/255, 255, 255/),(/000, 000, 000/),(/000, 230, 255/),(/000, 255, 255/),(/000, 255, 191/),(/000, 255, 135/),(/000, 255, 000/),(/167, 255, 000/),(/207, 255, 000/),(/255, 255, 000/),(/255, 135, 000/),(/255, 095, 000/),(/255, 039, 000/),(/255, 000, 119/),(/255, 000, 159/),(/255, 000, 199/)/) / 255.

; Plot resources

res 		= True
res@MainTitle	= fig_title
;res@Footer	= False

pltres	= True

mpres	= True

;------------- Set map limits -----------------------------------------------
x_start	= 290
x_end	= 450
y_start	= 170
y_end	= 340

;------------- Map resources ------------------------------------------------

mpres@ZoomIn	= True
mpres@Xstart 	= x_start
mpres@Xend 	= x_end
mpres@Ystart 	= y_start
mpres@Yend 	= y_end

mpres@mpGeophysicalLineColor      = "Black"
mpres@mpNationalLineColor         = "Black"
mpres@mpUSStateLineColor          = "Black"
mpres@mpGridLineColor             = "Black"
mpres@mpLimbLineColor             = "Black"
mpres@mpPerimLineColor            = "Black"
mpres@mpGeophysicalLineThicknessF = 2.0
mpres@mpGridLineThicknessF        = 2.0
mpres@mpLimbLineThicknessF        = 2.0
mpres@mpNationalLineThicknessF    = 2.0
mpres@mpUSStateLineThicknessF     = 2.0

; Load in the wrfout files

fList	= systemfunc("ls -1 wrfout*")
nFiles	= dimsizes(fList)

do iFile = 0, nFiles - 1

; Open a workstation

filename	= sprinti("wspd_plot_%03d",iFile)
print(filename)
wks		= gsn_open_wks("png",filename)
gsn_define_colormap(wks,cmap)
print("Workstation opened.")


; Open the wrfout file

a = addfile(fList(iFile),"r")
print("File added.")

; Get WRF variables

print("Getting WRF variables...")
u	= wrf_user_getvar(a,"ua",-1)
v	= wrf_user_getvar(a,"va",-1)
p	= wrf_user_getvar(a,"pressure",-1)
uvmet10	= wrf_user_getvar(a,"uvmet10",-1)
slp	= wrf_user_getvar(a,"slp",-1)
lat	= wrf_user_getvar(a,"lat",-1)
lon	= wrf_user_getvar(a,"lon",-1)

u10	= uvmet10(0,:,:,:)
v10	= uvmet10(0,:,:,:)

; Compute the wind speed

print("Computing wind speed...")
wspd	= (u10 * u10 + v10 * v10)^(0.5)

wspd@description	= "Wind Speed"
wspd@units		= "m/s"

; Convert from meters per second to knots

u_plane		= (u10 * 1.94386)
v_plane 	= (v10 * 1.94386)
u_plane@units	= "kts"
v_plane@units	= "kts"

; Apply wind speed mask

wspd_masked	= mask(wspd,(wspd.gt.35),True)

descr	= "Wind speed at 10 meters"

wspd_masked@description	= descr
wspd_masked@units	= "kts"

; Find location of minimum sea level pressure
print("Finding storm center...")
dsize = dimsizes(slp)
jmax = dsize(1)
imax = dsize(2)
slp_min = min(min(slp(0,y_start:y_end,x_start:x_end)))

do i = x_start, x_end
	do j = y_start, y_end
		if(slp(0,j,i).eq.slp_min) then
			min_lat = lat(0,j,i)
			min_lon = lon(0,j,i)
		end if
	end do
end do

; Plotting options for wind speed
print("Plotting...")

res                      = True
res@gsnDraw              = False
res@gsnFrame             = False
res@cnFillOn             = True
res@cnLinesOn            = False
res@lbLabelAutoStride    = True
times                    = wrf_user_list_times(a)
res@ValidTime            = True
res@TimeLabel            = times

opts					= res
opts@cnFillMode				= "RasterFill"
opts@cnRasterSmoothingOn		= True
opts@ContourParameters			= (/35, 95, 5/)

minpres = "Min SLP" + sprintf("%5.0f",slp_min) + " mb (Lat: " + sprintf("%4.1f",min_lat) + ", Lon: " + sprintf("%4.1f",min_lon) + ")"

pressureres		= True
pressureres@FieldTitle 	= minpres
pressureres@cnLineColor	= "Black"

wind_max = max(max(wspd(0,y_start:y_end,x_start:x_end)))

maxspd	= sprintf("Max winds = %.1f", max(max(wspd(0,170:340,290:450))))
maxspd = maxspd + " kt"

txres			= True
txres@txFontHeightF	= 0.015
gsn_text_ndc(wks,maxspd,0.3,0.15,txres)

contour_slp  				= wrf_contour(a, wks, slp(0,170:340,290:450), pressureres)
contour_spd				= wrf_contour(a,wks,wspd_masked(0,170:340,290:450),opts)
plot					= wrf_map_overlays(a,wks,(/contour_slp,contour_spd/),pltres,mpres)

; Write track data to text file
print("Writing WRF track data...")
alist	= [/times, min_lat, min_lon, slp_min, wind_max/]
write_table("wrfrun_table.txt", "a", alist,"%s%5.2f%5.2f%5.0f%5.0f") 

print("wrfout complete.")
end do

print("Great success!!!")
end
