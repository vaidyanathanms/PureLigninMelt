# Plot densities
set key right bottom
p "dens_npt.xvg" u 1:2 title "300",\
  "../T_320/dens_npt.xvg" u 1:2 title "320",\
  "../T_340/dens_npt.xvg" u 1:2 title "340",\
  "../T_360/dens_npt.xvg" u 1:2 title "360",\
  "../T_380/dens_npt.xvg" u 1:2 title "380",\
  "../T_400/dens_npt.xvg" u 1:2 title "400",\
  "../T_420/dens_npt.xvg" u 1:2 title "420",\
  "../T_440/dens_npt.xvg" u 1:2 title "440",\
  "../T_460/dens_npt.xvg" u 1:2 title "460",\
  "../T_480/dens_npt.xvg" u 1:2 title "480",\
  "../T_500/dens_npt.xvg" u 1:2 title "500",\
  "../../run_5/T_300/dens_npt.xvg" u 1:2 title "300_new"
pause -1 
