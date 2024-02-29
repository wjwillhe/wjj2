load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
;;--------------------------------------------------------------------
; wavelet_3.ncl
;
; Concepts illustrated:
;   - Computing wavelets
;   - Overlaying a stipple pattern to show area of interest
;   - Overlaying a hatch pattern to show area of interest
;   - Changing the size of the dot fill pattern in gsn_contour_shade
;   - Increasing the density of fill patterns
;--------------------------------------------------------------------
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;--------------------------------------------------------------------
; CREATE SAME WAVELET FIGURE AS TORRENCE & COMPO using seasonal Nino3
; from 1871-1997
;--------------------------------------------------------------------

begin
 
; data1=asciiread("/public/home/wangji19/data/NCEP_TRACK/mg1958-2020/JRA55/JRA55_smg_fre.txt",-1,"float")

data1=asciiread("/public/home/wangji19/data/ERA5_TRACK/20220830_MAM/year_fre_1950_2022.txt",-1,"integer")

  time = ispan(1950,2022,1)
 c_dtrend  =  dtrend_msg_n (time,data1,False,True,0);;;;;;;;;;;;;year_maxvor;;;去趋势
 ninoseas = dim_standardize(c_dtrend,0) ;;;(1:71);;;;;;;;;;;;;;;标准化
; ninoseas = data1
  ninoseas!0 = "time"
  ninoseas&time = time
  N     = dimsizes(c_dtrend);;;;(1:71)
;************************************
; compute wavelet
;************************************
  mother  = 0
  param   = -1
  dt      = 1    ;timesteps in units of years
  s0      = dt
  dj      = 0.25
  jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.))
  npad    = N
  nadof   = 0
  noise   = 1
  siglvl  = .10
  isigtest= 0

  w = wavelet(ninoseas,mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)
;************************************
; create coodinate arrays for plot
;************************************
  power            = onedtond(w@power,(/jtot,N/))
  power!0          = "period"                        ; Y axis
  power&period     = w@period                      ; convert period to units of years

  power!1          = "time"                          ; X axis
  power&time       = time

  power@long_name  = "Power Spectrum"
  power@units      = "1/unit-freq"

; compute significance ( >= 1 is significant)
  SIG              = power                            ; transfer meta data
  SIG              = power/conform (power,w@signif,0)
  SIG@long_name    = "Significance"
  SIG@units        = " "
  q=w@signif
  printVarSummary(power)
  ; print(w@signif)
;********************************************************************************
; initial resource settings
;****************************************************************maxvor****year_maxvor-33_fre************
  wks = gsn_open_wks("eps","/public/home/wangji19/program2/the year_fre wavelet_MAM 1950-2022")  ; send graphics to PNG file
  
 gsn_define_colormap(wks,"MPL_YlOrRd"); ;"BlAqGrYeOrReVi200" 
  YLValues          = (/1,3,6,10,15,20,30,40,60,70/)
  YLLabels          = (/"1","3","6","10","15","20","30","40","60","70"/)
  res                     = True                  ; plot mods desired
  res@gsnMaximize         = False
  res@gsnDraw             = False                 ; Do not draw plot
  res@gsnFrame            = False                 ; Do not advance frome
  res@gsnLeftString=""
  res@gsnRightString=""
  res@tmXBLabelFontHeightF    = 0.02    ;set the size of x-Axis words
  res@tmYLLabelFontHeightF    = 0.02   ;set the size of y-Axis words

  ; res@vpHeightF           = .4                    ;
  ; res@vpWidthF            = .5
 ; res@tiMainString  = "year_fre" ;;;;;;70_maxvor
  
  ; res@gsnRightString       = "Wavelet Power"
  
  res@trYReverse          = True 
   res@tmYLMode            = "Explicit"
  res@tmYLValues          = YLValues
  res@tmYLLabels          = YLLabels
  res@tmLabelAutoStride   = True
   res@trYMaxF             =max(YLValues)
  ; res@trYMinF             =min(YLValues)
   res@cnLinesOn           = False                 ; turn off contour lines
  res@cnInfoLabelOn       = False
 res@cnLineLabelsOn      = False       ; do not draw contour labels
    res@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res@cnMinLevelValF      = 1.00        ; set min contour level
  res@cnMaxLevelValF      = 7.00        ; set max contour level
  res@cnLevelSpacingF     = 1.00 
res2 = res                            ; res2 probability plots

  ; res@tiXAxisString        = "Year"
  res@tiXAxisOffsetYF        = 0.135
 res@tiYAxisString        = "Years"
 
  res@cnFillOn            = True                  ; turn on color
  
  
 res@cnFillMode          = "RasterFill"          ; turn on raster mode
  res@cnRasterSmoothingOn = True 
  ; res@cnFillOpacityF      = 0.50                  ; make partially transparent
  ; res@cnLineLabelsOn      = False
  ; res@cnInfoLabelOn       = False
  ; res@cnRasterSmoothingOn = True                  ; turn on raster smoothing
; res@gsnLeftString        = "NINO3: GISST"
 
  ; res@tiMainString = "Wavelet of the es path cyclone annual frequency"
 ; res2@gsnDraw             = False       ; Do not draw plot
  ; res2@gsnFrame            = False       ; Do not advance frome
  ; res2@vpHeightF           = .4                    ;
  ; res2@vpWidthF            = .5

  res2@cnLevelSelectionMode = "ManualLevels" ; set manual contour levels
  res2@cnMinLevelValF      = 0.00        ; set min contour level
  res2@cnMaxLevelValF      = 6.00        ; set max contour level
  res2@cnLevelSpacingF     = 1.00       ; set contour spacing
  res2@cnFillScaleF        = 0.5
;---Create the same plot three times so we can compare different fill scales and dot densities.

  ; res@tiMainString = "Wavelet of the es path cyclone annual frequency"
  
;---Create the same plot three times so we can compare different fill scales and dot densities.
 res@lbLabelFontHeightF   = 0.020
 res@pmLabelBarOrthogonalPosF  =0.01
 res@pmLabelBarHeightF    =0.12
 res@pmLabelBarWidthF     =0.6
  plot = gsn_csm_contour(wks,power,res)
 

  iplot = gsn_csm_contour(wks,SIG,res2)

;
; Apply three different types of stippling to the three plots:
;   - Shading with hatch pattern
;   - Shading with increased density of hatch pattern
;   - Shading with increased density of stipple and dot size.
;
  opt                      = True
  opt@gsnShadeHigh         = 17           ; cross hatch pattern
  opt@gsnShadeFillType     = "pattern"
  iplot = gsn_contour_shade(iplot,-999., 1.0, opt)

  overlay(plot,iplot)
  plot=ShadeCOI(wks,plot,w,time,False) 
  
  ;;;;
  XBLabels=(/1.0,2.0,3.0/)
  XBValues=(/"1.0","2.0","3.0"/)
  gws = w@gws
  signif=w@signif
  ; print(max(signif))
  ; print(max(gws))
  
  ; scale = w@scale
  ; Cdelta = w@cdelta
  ; powernorm = power
  ; powernorm = power/conform(power,scale,0)
  ; scaleavg = dj*dt/Cdelta*dim_sum_Wrap(powernorm(time|:,{period|2.:8.}))
  
  resl = True
  resl@vpHeightF           = .4                    ;
  resl@vpWidthF            = .2
   resl@gsnMaximize         = False

  resl@gsnFrame = False
  resl@gsnDraw = False
  resl@tmXTOn = False
  resl@tmXBOn = True
  resl@tmXBMode="Explicit"
  resl@tmXBLabels=XBLabels;;;;想标记在坐标轴上的字符
  resl@tmXBValues=XBValues;;;;对应字符的位置
 resl@trYAxisType = "LogAxis"
  resl@trYReverse          = True    
  resl@tmYLMode = "Explicit"
  resl@tmYLValues          = YLValues
  resl@tmYLLabels          = YLLabels              ; reverse y-axis
  resl@trYMaxF             =max(YLValues)
  resl@trYMinF             =min(YLValues)
  ; plotg = gsn_csm_xy(wks,gws,power&period,resl)
  

  resl2=  resl  
  
  resl2@tmXBMode="Explicit"
  resl2@tmXBLabels=XBLabels;;;;想标记在坐标轴上的字符
  resl2@tmXBValues=XBValues;;;;对应字符的位置
  resl2@xyLineColor="red"
  resl2@xyLineThicknessF=1.5
  resl2@xyDashPattern=0
  ; plotg1 = gsn_csm_xy(wks,signif,power&period,resl2)
  ; plotg = gsn_csm_x2y(wks,gws,signif,power&period,resl,resl2)
  
  ; overlay(plotg,plotg1)
  
  
  ; plotc = gsn_attach_plots(plot,plotg,res,resl)
  
  draw(plot)
  
  frame(wks)


end