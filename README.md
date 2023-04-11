# LULC-PBS
Regional LULC maps of three target areas located in the main ecoregions of Ecuador at a resolution of 10 m using Google Earth Engine (GEE) cloud-based computing. 

How to use?

Users with access to the repository can add it to the Code Editor using:

https://code.earthengine.google.com/?accept_repo=users/SI_Testing/LULC_EC

Find the next section where you can select the map of the region you want to process. 
To do this you must uncomment the var scriptClassifier and AOI in the option of your preferred.

//********** Visualising LULC-PBS in the COAST Region of Ecuador******//
// If you'd like to visualise LULC-PBS in coast region, you should uncomment next line 
//**********************************************************************************
// Set scripts SDC and PBS from Coast
//var scriptClassifier = lulcClassCoast
// Set Area of Interest (AOI)
//var AOI= AOI_Co.geometry(); 


//********** Visualising LULC-PBS in the ANDEAN Region of Ecuador******//
// If you'd like to visualise LULC-PBS in Andean region, you should uncomment next line 
//************************************************************************************
// Set scripts SDC and PBS from Andean
//var scriptClassifier = lulcClassAndean
// Set Area of Interest (AOI)
//var AOI= AOI_Hi; 


//********** Visualising LULC-PBS in the AMAZON Region of Ecuador******//
// If you'd like to visualise LULC-PBS in Amazon region, you should uncomment next line 
//*************************************************************************************
// Set scripts SDC and PBS from Amazon
//var scriptClassifier = lulcClassAmazon
// Set Area of Interest (AOI)
//var AOI= AOI_Am; 

