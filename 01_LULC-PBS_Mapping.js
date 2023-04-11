/* Users with access to the repository can add it to the Code Editor using:
https://code.earthengine.google.com/?accept_repo=users/SI_Testing/LULC_EC
*/

// ##########################   Calling scripts: #####################################
// ###################################################################################
var lulcClassCoast = require('users/SI_Testing/LULC_EC:Scripts/SDC_PBS_Coast');
var lulcClassAndean = require('users/SI_Testing/LULC_EC:Scripts/SDC_PBS_Andean');
var lulcClassAmazon = require('users/SI_Testing/LULC_EC:Scripts/SDC_PBS_Amazon');

// ##########################   Definitions: #########################################
// ################################################################################### 
// Definitions of the regions 
// Coast area
var ecuCoast = ee.FeatureCollection("users/SI_Testing/ECU_adm2");
var AOI_Co = ee.FeatureCollection(ecuCoast.filter(ee.Filter.stringContains('NAME_2', "Daule")));
// Andean Highland area
var AOI_Hi =  /* color: #52cfd6 */ee.Geometry.Polygon(          
[[[-78.27937421266317, 0.14859741950997463],
          [-78.27937421266317, -0.012420785388960836],
          [-77.94635114137411, -0.012420785388960836],
          [-77.94635114137411, 0.14859741950997463]]]),       
// Amazon area
AOI_Am = /* color: #d63000 */ee.Geometry.Polygon(
[[[-76.96032781586344, -0.19799757372273413],
[-76.96032781586344, -0.5056085930148337],
[-76.63314123139078, -0.5056085930148337],
[-76.63314123139078, -0.19799757372273413]]])
;

// Set RGB visualization parameters for S2 images
var rgbVis = {
  min: [800, 800, 800],
  max: [3500, 3500, 3500],
  bands: ['B4', 'B3', 'B2']
};

// ##########################   User Inputs: #########################################
// ###################################################################################
// Set dates
var start_date='2019-01-01';//
var   end_date='2020-12-30';

// Set the percentage of clouds in the collection
var max_cloud_percent=75;


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
var scriptClassifier = lulcClassAmazon
// Set Area of Interest (AOI)
var AOI= AOI_Am; 
Map.centerObject(AOI,12);

// ####################### Make Image Collections ####################################
// ###################################################################################

var S2_filtered_Q = ee.ImageCollection("COPERNICUS/S2")
                            .filterBounds(AOI)
                            .filterDate(start_date,end_date)
                            .filterMetadata("CLOUDY_PIXEL_PERCENTAGE",'less_than',max_cloud_percent); 
                            print(S2_filtered_Q, 'filtered2');
        
// Bits 10 and 11 are clouds and cirrus, respectively.
var cloudBitMask = ee.Number(2).pow(10).int();
var cirrusBitMask = ee.Number(2).pow(11).int();
          
function maskS2clouds(image) {
  var qa = image.select('QA60');
  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0).and(
             qa.bitwiseAnd(cirrusBitMask).eq(0));
return image.updateMask(mask);}

// Remove clouds
var ST2_nocloud_Q = S2_filtered_Q.filterBounds(AOI).map(maskS2clouds);
Map.addLayer(ST2_nocloud_Q.median().clip(AOI),rgbVis, 'S2_nocloud',false);

//*********************************** VIIRS dataset *************************//
var dataset = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG')
                  .filterDate(start_date,end_date)
var nighttime = dataset.select('avg_rad');
// 1: resample annual VIIRS to match Sentinel-2 projection
// NOTE: there is a PROBLEM with the RESAMPLING METHOD -- SMOOTHS OUTPUT TOO MUCH
// 2:  get projection and mask from Sentinel-2 
var S2_filtered_py = S2_filtered_Q.first(); 
var imViirs = nighttime.first().select('avg_rad').projection().nominalScale();
var S2_proj = S2_filtered_py.select('B4').projection();
var viirs_annual_resamp = nighttime.map(function(im) {
                                            return im
                                          .resample('bilinear')
                                          .reproject(S2_proj.atScale(10))
                                          .set('rad',im.get('avg_rad'));  });        
var imVIIRS= viirs_annual_resamp.select('avg_rad').first();


// ####################### Pre- processing monthly images#############################
// ###################################################################################
// Calculate the number of images per month
var months = ee.List.sequence( 1, 12);

//Create a stack of images per month
var imagesMonthly = months.map( function(month ){
                                   month = ee.Number( month );
                                   var monthImg = ST2_nocloud_Q
                                                         .filter( ee.Filter.calendarRange(  month, month, "month") )
                                                         .median();
                                   // Define a kernel.
                                   var kernel = ee.Kernel.circle({radius: 2});
                                   // Perform an erosion followed by a dilation, display.
                                   monthImg = monthImg
                                                    .focal_min({kernel: kernel, iterations: 1})
                                                    .focal_max({kernel: kernel, iterations: 1});
                                  monthImg =monthImg.set('month', month);
                                  return( ee.Algorithms.If( ST2_nocloud_Q.size().gt(0), monthImg, 0 ));
                                  });

var monthlyDataset = ee.ImageCollection(imagesMonthly);

// ##################### Processing Single Date Classification (SDC)###################
// ###################################################################################  
var bandsS2 = ['B2','B3','B4','B8','B11','B12'];
var colIm = monthlyDataset.map(function(image)
                                {return scriptClassifier.SDC_dataset(image.clip(AOI).divide(10000),bandsS2,imVIIRS)});
var lisImg = ee.List([]);
var cl = colIm.toList(colIm.size());
for (var i = 0; i < cl.length().getInfo(); i++) {
  var im = ee.Image(cl.get(i)).set('month', i);
  im = im.unmask(1);
  lisImg = lisImg.add(im);
  } 
var collSDC = ee.ImageCollection.fromImages(lisImg);

// ##################### Processing Phenology Based Synthesis Classification (PBS)###################
// ################################################################################################## 

var PBS_OUT=scriptClassifier.PBS_classification(collSDC);
          Map.addLayer(PBS_OUT.clip(AOI), scriptClassifier.visPBS, 'PBS_S2',false);
          
// ####################################### Export LULC-PBS map ####################################
// ################################################################################################## 
       
 Export.image.toDrive(
    {image:PBS_OUT.multiply(100).round().toInt16(),   
    description: 'PBS_S2_2019_2020Amazon',
    fileNamePrefix:'PBS_S2_1920Am',
     folder: 'PBStest', 
     region:AOI,
     scale:10, 
     skipEmptyTiles:true,
     maxPixels: 1e13
    }); 

