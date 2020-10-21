# Nature_Comm_CSD
Code used to process imaging and EEG data as described in Smith et al., 2020.
CSD_Image_EEG_Wrapper.m processes imaging and EEG data for a mouse and calls 3 sub-routines: 
  makemask.m
  proc_Imaging.m
  proc_EEG.m
Supporting .txt files include spectra for LEDs used on imaging system or extinction coefficients and are called within proc_Imaging.m. 
Supporting .tif file is called in proc_Imaging.m to subtract ambient light around imaging system.
Supporting .m files are called in proc_Imaging.m and perform various data processing procedures.
