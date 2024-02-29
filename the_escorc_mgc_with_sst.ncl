load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;************************************************
begin
;;;;;;;气旋
;;;;;各类路径气旋的位置

 year=ispan(1950,2022,1)
data1=asciiread("/public/home/wangji19/data/ERA5_TRACK/grid1/year_strong_mg_fre_1950_2022.txt",-1,"float")
c_dtrend  =  dtrend_msg_n (year,data1(0:72),False,True,0);;year_maxvor;year_maxvor33_fre;;;;;;year_maxvor;;;;;;;去趋势
x = dim_standardize(c_dtrend,0) ;;;;;;;;;;;;;;;;;;标准化

; filep = "/public/home/wangji19/data/index/lijianping-NAO/NAO1948-present.txt"
; nrs = numAsciiRow(filep)   ;;;; contributed.ncl
; ncs = numAsciiCol(filep)
; filepdat  = asciiread(filep,(/nrs,ncs/),"float")
; filepdat@_FillValue = -99999.0
; ; ; printVarSummary(dat)
; a_MAM=filepdat(0:72,13)

; x=dim_standardize(a_MAM,0);a_MAM
 ; year=ispan(1948,2020,1)
 ; p_dtrend  =  dtrend_msg_n (year,a_MAM,True,True,0);;;;;;;;;;;;;;;;去趋势
; x=dim_standardize(p_dtrend,0)

 ; ft1=addfile("/public/home/wangji19/program2/data/fre_NAsst_index_MAM.nc","r") 
 ; na1    = ft1->index 
  ; ; ; c_dtrend  =  dtrend_msg_n (year,na1,False,True,0);;;year_maxvor;;年际; ;;33_fre;;;;;year_maxvor;;;;;;;;;去趋势
; ; ; x = dim_standardize(c_dtrend,0) ;;;;;data1(1:71);;;;;;;;;;;;;标准化
 ; x = dim_standardize(na1,0) ;;;;;;;;;;;;;;;;;;标准化
 
 ; ft2=addfile("/public/home/wangji19/program2/data/fre/fre_BS_SST_Jan.nc","r") 
 ; na2    = ft2->index 
; x2 = dim_standardize(na2,0) ;;;;;;;;;;;;;;;;;;标准化
; x=(dim_standardize(na2,0)-dim_standardize(na1,0))/2.0
; ;;;;;;
; f=addfile("/cygdrive/E/data/mondata/ERA5/sst.mon.1950-2021.0-180.nc","r") 
; f=addfile("/public/home/wangji19/data/mondata/slp.mon.mean.nc","r") 
f=addfile("/public/home/wangji19/data/mondata/ERA5/sst.1949_202306.nc","r") 
; f=addfile("/public/home/wangji19/data/mondata/HadISST_sst.nc","r") 
; f=addfile("/public/home/wangji19/data/mondata/ERSST5.sst.mnmean.nc","r") 
it_s=194901
 it_e=202212
 tim    = f->time        ; time:units = "hours since 1-1-1 00:00:0.0"                               
                           ;;time:units = "hours since 1-1-1 00:00:0.0"                               
 TIME1    = cd_calendar(tim, -1)  
 rec_s=ind(it_s.eq.TIME1)
 rec_e=ind(it_e.eq.TIME1)  
 
  delete(TIME1)
  delete(tim)
   ; hgt     =  f->slp(rec_s:rec_e,:,:)
   ; hgtt     = new((/3,74,73,144/),float) 
   ; hgta     = new((/3,74,73,144/),float)
   
  hgt     = short2flt( f->sst(rec_s:rec_e,:,:))
  hgt = lonFlip(hgt) 
  ; hgt@_FillValue=-32767doubletofloat()
  hgt=hgt-273.15
  
; copy_VarMeta(hg,hgt)
  delete(f)
  printVarSummary(hgt)
dims=dimsizes(hgt)
    ; exit

   hgtt     =  new((/3,73,dims(1),dims(2)/),float) ; new((/3,72,321,1440/),float) ;
   hgtt(0,:,:,:) = dtrend_msg_n(year,hgt(14:881:12,:,:),False,True,0);hgt(14:859:12,:,:);;842 ;
   hgtt(1,:,:,:) = dtrend_msg_n(year,hgt(15:881:12,:,:),False,True,0);hgt(15:859:12,:,:) ;
   hgtt(2,:,:,:) = dtrend_msg_n(year,hgt(16:881:12,:,:),False,True,0);hgt(13:859:12,:,:) ;
  hgtcl = dim_avg_n (hgtt,0)
  
 ; hgtcl     =  new((/71,dims(1),dims(2)/),float) ; new((/3,72,321,1440/),float) ;
 ; hgtcl(:,:,:) = hgt(12:856:12,:,:);;;
  hgtcl!0="year"
   hgtcl&year=year
copy_VarMeta(hgt(0,:,:),hgtcl(0,:,:))
; hhgtcl     =  new((/71,dims(1),dims(2)/),float)
; hhgtcl=dtrend_msg_n(year,hgtcl(:,:,:),False,True,0)
; copy_VarMeta(hgtcl(:,:,:),hhgtcl(:,:,:))

  printVarSummary(hgtcl)
confi=new((/dims(1),dims(2)/),float)
  n    = dimsizes(x)                ; n=11
  df   = n-2
 t=new((/dims(1),dims(2)/),float)
t@_FillValue=-1e+30
  ccr=new((/dims(1),dims(2)/),float)
  ccr@_FillValue=-1e+30
 
do i=0,dims(1)-1
do j=0,dims(2)-1
ccr(i,j)=escorc(x,hgtcl(:,i,j)) 
t(i,j)=ccr(i,j)*sqrt(df/(1-ccr(i,j)^2)) ;;;;;;;;;;;;;;;;;;;;;escorc
confi(i,j)=student_t(t(i,j),df)

; ccr(i,j)=regCoef(x,hgtcl(:,i,j)) 
 ; confi(i,j)  = student_t(ccr@tval, ccr@nptxy-2);cor_h;;;cor_h;;;;;;;;;;regcof
end do
end do
copy_VarMeta(hgtcl(0,:,:),ccr(:,:))
copy_VarMeta(hgtcl(0,:,:),confi(:,:))
print(max(ccr))


printVarSummary(confi)
print(max(confi))
print(min(confi))
; delete(hgtt)
delete(hgt)
printVarSummary(ccr)
; outfile=addfile("/public/home/wangji19/program2/cesm/dsst_fre_ccr_MAM.nc","c")  
  ; outfile->index=ccr
  ; exit
 ;;;;;;;;;;;;;
 ; hClm = clmMonTLL( hgt )
      ; printVarSummary(hClm)
  ; hAnom = calcMonAnomTLL (hgt,hClm)
     ; ; hgtc     =   new((/3,71,dims(1),dims(2)/),float) ; new((/3,72,321,1440/),float)  ;
   ; ; hgtc(0,:,:,:) = hAnom(14:859:12,:,:);dim_avg_n( dtrend_msg_n(year,hgt(2:863:12,:,:),True,True,0),0)
   ; ; hgtc(1,:,:,:) = hAnom(15:859:12,:,:);dim_avg_n(  dtrend_msg_n(year,hgt(3:863:12,:,:),True,True,0),0)
   ; ; hgtc(2,:,:,:) = hAnom(16:859:12,:,:);dim_avg_n(  dtrend_msg_n(year,hgt(4:863:12,:,:),True,True,0),0)
 ; ; hanoc= dim_avg_n (hgtc,0)
 
 ; hanoc     =   new((/71,dims(1),dims(2)/),float) ; new((/3,72,321,1440/),float)  ;
   ; hanoc(:,:,:) = hAnom(12:859:12,:,:);dim_avg_n( dtrend_msg_n(year,hgt(2:863:12,:,:),True,True,0),0)
  
   ; hanoc!0="year"
   ; hanoc&year=year
; copy_VarMeta(hgt(0,:,:),hanoc(0,:,:))

  ; ; ; ; ;;;;;;计算区域平均指数
rad = 3.1415926/180
lat1 = hgtcl&latitude;;
clat = cos(rad*lat1)
clat!0 = "lat"
clat&lat = lat1 
 

; inda=new(73,float)
 ; inda= wgt_areaave(hgtcl(:,{40:50},{325:350}),clat({40:50}),1.0,0)
  ; ; inda= wgt_areaave(hanoc(:,{40:50},{-45:-30}),clat({40:50}),1.0,0)

; outfile=addfile("/public/home/wangji19/program2/data/fre_NAsst_index_MAM.nc","c")  
  ; outfile->index=-inda
  ; exit
; indb=new(71,float)
 ; indb= wgt_areaave(hanoc(:,{77.5:82.5},{40:60}),clat({77.5:82.5}),1.0,0)

 ; outf=addfile("/public/home/wangji19/program2/data/fre/fre_BS_SST_Jan.nc","c")  
  ; outf->index=indb
  
 ; indc=new(71,float)
 ; indc= wgt_areaave(hanoc(:,{60:65},{-15:-5}),clat({60:65}),1.0,0)
  ; outf2=addfile("/public/home/wangji19/program2/data/fre/fre_BS_SST_Jan.nc","c")  
  ; outf2->index=indc
 
 ; indd=new(71,float)
 ; indd= wgt_areaave(hanoc(:,{60:65},{-15:-5}),clat({60:65}),1.0,0)
  ; outf3=addfile("/public/home/wangji19/program2/data/fre/fre_BS_SST_Jan.nc","c")  
  ; outf3->index=indd
 
 delete(hgtcl)
 ; delete(hhgtcl)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;HadISST;;;;;;;;;;year_maxvor;fre;;;;regcof;year_maxvor33_fre;;year_maxvor-33_fre;;;;maxvor;;;
 wks = gsn_open_wks("eps","/public/home/wangji19/program2/The year_strong_mgc_fre dtrend sst escorc MAM")  
 ;wks = gsn_open_wks("eps","/cygdrive/E/program2/The_maxspeed_sst_regcof_mam1");;;fre_NAT_BS_SST_Jan;;;;;;;;  
 gsn_define_colormap(wks,"cmp_b2r" )

 res=True
res@gsnAddCyclic = True
;res@gsnMaximize  = True
res@gsnDraw=False
res@gsnFrame=False
res@gsnLeftString=""
res@gsnRightString=""
res@tmXBLabelFontHeightF    = 0.02     ;set the size of x-Axis words
res@tmYLLabelFontHeightF    = 0.02   ;set the size of y-Axis words
; res@pmTickMarkDisplayMode="Always"

mpres=res
 mpres@mpOutlineOn=True
mpres@mpFillOn=False
mpres@mpPerimOn = 	True
mpres@mpPerimDrawOrder = "PostDraw"
mpres@mpPerimLineThicknessF = 0.1
mpres@mpGeophysicalLineThicknessF=0.1
mpres@mpLimitMode="LatLon"
mpres@mpCenterLonF=180;;-35;;90 ;;-35;
mpres@mpMinLatF = -35;;10;
mpres@mpMaxLatF = 90;75;;
mpres@mpMinLonF = 0;-90;;-90;;;
mpres@mpMaxLonF =360;30;;270;;20;;
;;;;;;;;;;;xy经纬度间隔
mpres@tmXTOn=False
mpres@tmYROn=False
; mpres@tmXBMode="Explicit"
; mpres@tmXBValues=(/-90,-60,-30,0,30/)
; mpres@tmXBLabels=(/"90~S~o~N~W","60~S~o~N~W","30~S~o~N~W","0~S~o~","30~S~o~N~E"/);
; mpres@tmYLMode="Explicit"
; mpres@tmYLValues=(/60,30/);;;"40~S~o~N~N",,"10~S~o~N~N"
; mpres@tmYLLabels=(/"60~S~o~N~N","30~S~o~N~N"/);
; mpres@gsnMinorLonSpacing=10
; mpres@gsnMajorLonSpacing=90
; mpres@gsnMinorLatSpacing=10
; mpres@gsnMajorLatSpacing=20
; mpres@mpDataSetName              = "/public/home/wangji19/data/MAP/NCL-Chinamap-master/database/Earth..4";;;"/cygdrive/d/MAP/NCL-Chinamap-master/database/Earth..4"
; ;;
; mpres@mpDataBaseVersion          = "Ncarg4_1" ; or "MediumRes"
; mpres@mpAreaMaskingOn            = True
; mpres@mpMaskAreaSpecifiers       = (/"China"/)
; mpres@mpOutlineSpecifiers        = (/"China","China:Provinces"/)

; res@mpGeophysicalLineThicknessF  = 0.2
; res@mpNationalLineThicknessF   = 0.2
; res@mpProvincialLineThicknessF = 0.05
  ; res@mpProvincialLineThicknessF = 0.05    
; res@mpLandFillColor            = "white"
; res@mpInlandWaterFillColor     = "white"
; res@mpOceanFillColor           = "white"
; res@mpOutlineBoundarySets      = "NoBoundaries"
; res@mpNationalLineColor        = "black"
; res@mpProvincialLineColor      = "black"
; res@mpGeophysicalLineColor     = "black"
map=gsn_csm_map(wks,mpres)

; mpres=res
; mpres@gsnPolar="NH"
; mpres@mpFillOn=False
; mpres@mpPerimOn = 	True
; mpres@mpPerimDrawOrder = "PostDraw"
; mpres@mpPerimLineThicknessF = 0.1
; mpres@mpGridAndLimbOn  = 	True
; mpres@mpGridLatSpacingF    = 45.                ; change latitude  line spacing
; mpres@mpGridLonSpacingF    = 90.
; mpres@mpGridLineThicknessF=0.05
; mpres@mpGeophysicalLineThicknessF=0.1

; mpres@mpCenterLonF=90
; mpres@mpMinLatF=10
; mpres@mpMaxLatF=90

; map=gsn_csm_map_polar(wks,mpres)
resc=res
resc@gsnAddCyclic = True
resc@cnFillOn=True
resc@cnLinesOn=False
; resc@cnLineDashPattern=0
; resc@cnLineThicknessF=1
;resc@cnLineLabelsOn=False
; resc@cnInfoLabelOn=False
;;;;;;;;;;;;;;;;;;;;;;色标设置
resc@lbLabelBarOn         =  True
resc@cnInfoLabelOn          = False;;
resc@lbLabelFontHeightF   =0.02
resc@pmLabelBarOrthogonalPosF  =0.2;latlon0.1;;nh;0.01
resc@pmLabelBarHeightF    =0.15
resc@pmLabelBarWidthF     =0.6
  resc@gsnLeftString="(a)" 
  resc@gsnLeftStringParallelPosF=0.02
  resc@gsnLeftStringOrthogonalPosF=-0.11
  resc@gsnLeftStringFontHeightF   =0.04;;;0.05
resc@cnLevelSelectionMode = "ManualLevels"       ; set manual contour levels
resc@cnMinLevelValF       = -1.0              ; set min contour level
resc@cnMaxLevelValF       = 1.0              ; set max contour level
resc@cnLevelSpacingF      = 0.2;;0.08   
; resc@cnLevelSelectionMode="ExplicitLevels"
; resc@cnLevels=(/-2,-1.6,-1.2,-0.8,-0.4,0,0.4,0.8,1.2,1.6,2/)
; resc@cnLevels=(/-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0/);,0,-0.1,0,0.1-0.1,0.1,
;resc@cnLevels=(/-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5/);0,
; resc@cnLevels=(/-0.25,-0.2,-0.15,-0.1,-0.05,0,0.05,0.1,0.15,0.2,0.25/)
; resc@cnFillColors=(/2,7,13,18,24,30,35,45,50,55,60,65/);
plot1=gsn_csm_contour(wks,ccr(:,:),resc)
; plot1=gsn_csm_contour(wks,rc_slp(:,:),resc)

overlay(map,plot1)

;;;;;;;;;;;;;;;显著区域打点 ------------------
res2 = res
res2@cnFillOn                = True
res2@cnLinesOn               = False
res2@cnLevelSelectionMode   = "ExplicitLevels"
res2@cnLevels               = (/0.1/) ; significance level0.01,
; res2@gsnLeftString          ="sig.@0.1"

res2@cnFillPalette          ="GMT_gray"
res2@cnFillColors           =(/2,-1/)
res2@cnFillPatterns         = (/17,-1/);;;;填充形状4,17,,0.5
res2@cnFillScales           = (/0.5,1/)
res2@cnFillDotSizeF         = 0.001;0.003

res2@cnMonoFillPattern      = False
res2@cnMonoFillScale        = False
res2@cnInfoLabelOn          = False
res2@cnLineLabelsOn         = False
res2@lbLabelBarOn           = False
plot_p=gsn_csm_contour(wks,confi(:,:),res2)
overlay(map,plot_p)

;;;;;;添加方框
gonres=True
gonres@gsFillOpacityF =0
gonres@gsLineThicknessF=1
gonres@gsLineColor="black"
txres=True
txres@txFontHeightF=0.015
txres@txFontColor="black"

; x_lat0=67.5
; x_lon0=360
; box_lat0=(/65,70,70,65,65/)
; box_lon0=(/355,355,5,5,355/)
; dum0=gsn_add_polyline(wks,map,box_lon0,box_lat0,gonres) 
; dum00=gsn_add_text(wks,map,"A",x_lon0,x_lat0,txres) 

x_lat1=50
x_lon1=320
; box_lat1=(/40,50,50,40,40/)
; box_lon1=(/325,325,350,350,325/)
; dum1=gsn_add_polyline(wks,map,box_lon1,box_lat1,gonres)
; dum11=gsn_add_text(wks,map,"B",x_lon1,x_lat1,txres)  

; x_lat2=30
; x_lon2=305
; box_lat2=(/27.5,32.5,32.5,27.5,27.5/)
; box_lon2=(/-62.5,-62.5,-52.5,-52.5,-62.5/)
; dum2=gsn_add_polyline(wks,map,box_lon2,box_lat2,gonres)
; dum22=gsn_add_text(wks,map,"C",x_lon2,x_lat2,txres)   

; x_lat3=12.5
; x_lon3=330
; box_lat3=(/10,15,15,10,10/)
; box_lon3=(/325,325,335,335,325/)
; dum3=gsn_add_polyline(wks,map,box_lon3,box_lat3,gonres) 
; dum33=gsn_add_text(wks,map,"D",x_lon3,x_lat3,txres)  

; mkres=True
; mkres@gsMarkerColor="red"
; mkres@gsMarkerSizeF=0.05
; mkres@gsMarkerIndex=1
; mk1=gsn_add_polymarker(wks,map,lon1,lat1,mkres)
; mk2=gsn_add_polymarker(wks,map,lon2,lat2,mkres)

draw(map)
frame(wks)
end

exit
;;;;;;;;;;;;;;;
; it_s=194801
 ; it_e=202012
 ; tim    = f->time        ; time:units = "hours since 1-1-1 00:00:0.0"                               
                           ;;;time:units = "hours since 1-1-1 00:00:0.0"                               
 ; TIME1    = cd_calendar(tim, -1)  
 ; rec_s=ind(it_s.eq.TIME1)
 ; rec_e=ind(it_e.eq.TIME1)  
  
   ; hgt     =  f->slp(rec_s:rec_e,:,:)
   ; hgt     =  f->sst(rec_s:rec_e,:,:)

 ; hgt    = where(abs(hgt).gt.100,hgt@_FillValue,hgt) 

 ; hcl=month_to_season(hgt,"MAM")
 ; hanoc=dim_rmvmean_n_Wrap(hcl,0)
  ; print(hanoc)

