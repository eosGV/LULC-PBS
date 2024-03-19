/* Users with access to the repository can add it to the Code Editor using:
https://code.earthengine.google.com/?accept_repo=users/SI_Testing/LULC_EC
*/

/**
 *  * SDC and PBS classification have been modified from Simonetti D, 
 *      https://publications.jrc.ec.europa.eu/repository/handle/JRC99695
 *  * Adding classes such as annual and perennial crops, urban areas, and Greenhouses.
 *  * Spectral index have been testing by region
 * 
 /**
*/
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------------------------------------------------------------------------------------------------
//  ----------     Visualization functions                         -----------------------------------------------
//-------------------------------------------------------------------------------------------------------------------------------

function rgb(r,g,b){
          var bin = r << 16 | g << 8 | b;
          return (function(h){
          return new Array(7-h.length).join("0")+h
          })(bin.toString(16).toUpperCase())
}


var CLASS_PBSvis = {
  'bands':'Class',
  'min': 1,
  'max': 11,
  palette: [
            
            rgb(67, 27, 224), //1  Water -----
            rgb(12, 64, 13),    //2  Forest
            rgb(255,153,255),   //3   Annual Crop
            rgb(255,0,127),  //4 Perennial Crop
            rgb(27,229,5),   //5 grassland 
            rgb(255,128,0),  //6 shrub
            rgb(255,213,0),   //7 bare Soil
            rgb(252, 8, 8),   //8 urban  
            rgb(192,192,192),   //9 snow
            rgb(20,20,20),   //10  no data    //9  Agri
            rgb(112,40,143)     //11 Greenhouse
            
    ] //my palette
};
          
exports.visPBS = CLASS_PBSvis;

// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ------------------  SINGLE DATE CLASSIFICATION (SDC--------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

// INPUT 1: Monthly image to be classified 
// INPUT 2: bands conbination (B,G,R,NIR,SWIR1,SWIR2) ['B2','B3','B4','B8','B11','B12']
// OUT :  classified input(same number of layers); Class code is not the same as the PBS  

exports.SDC_dataset =  function (image,BANDS,imageVIIRS){

    var th_NDVI_MAX_WATER= 0;
    var BLU=image.select(BANDS[0]);
    var GREEN=image.select(BANDS[1]);
    var RED=image.select(BANDS[2]);
    var NIR=image.select(BANDS[3]);
    var SWIR1=image.select(BANDS[4]);
    var SWIR2=image.select(BANDS[5]);
             
    var OUT=ee.Image(0);
    
    var th_NDVI_SATURATION=0.0040;
    var th_NDVI_MIN_CLOUD_BARE=0.39;//-0.3
    var th_NDVI_MIN_VEGE=0.4;
    
    var th_SHALLOW_WATER=-0.12;
    var th_RANGELAND=0.53;//0.53
    var th_GRASS=0.58;//0.58
    var th_SHRUB=0.69;//0.69
    var th_TREES=0.80 ;
       
    var min123=BLU.min(GREEN).min(RED);
     
    var min1234=min123.min(NIR);
    var min234=GREEN.min(RED).min(NIR);
    
    var max234=GREEN.max(RED).max(NIR);
    var max1234=max234.max(BLU);
    
    var max57=SWIR1.max(SWIR2);
    var max457=max57.max(NIR);
    
    var max123457= max1234.max(max57);
    
    
    var BLUgtGREEN  = BLU.gt(GREEN);
    var BLUgteGREEN = BLU.gte(GREEN);
    var BLUlteNIR   = BLU.lte(NIR);
    
    var GREENgtRED  = GREEN.gt(RED);
    var GREENlteRED = GREEN.lte(RED);
    var GREENgteRED = GREEN.gte(RED);
    var REDlteNIR= RED.lte(NIR);
    
    var REDsubtractGREEN = (RED.subtract(GREEN)).abs();
    var BLUsubtractNIR   = BLU.subtract(NIR)
    
    var BLUgtGREENgtRED=BLUgtGREEN.and(GREENgtRED)
    
    var growing14=(BLU.lte(GREEN)).and(GREENlteRED).and(REDlteNIR);
    var growing15=growing14.and(NIR.lte(SWIR1));
    
    var decreasing2345=(GREENgteRED).and(RED.gte(NIR)).and(NIR.gte(SWIR1));
    
    
    var SATURATION=(max234.subtract(min234)).divide(max234);

    var WETNESS= image.expression('byte(b("'+BANDS[0]+'")*255)*0.2626 + byte(b("'+BANDS[1]+'")*255)*0.21 + byte(b("'+BANDS[2]+'")*255)*0.0926 + byte(b("'+BANDS[3]+'")*255)*0.0656 - byte(b("'+BANDS[4]+'")*255)*0.7629 - byte(b("'+BANDS[5]+'")*255)*0.5388');
    
    var NDVI= image.normalizedDifference(['B8', 'B4']);//(NIR.subtract(RED)).divide(NIR.add(RED));
    var NDSI=(BLU.subtract(SWIR1)).divide(GREEN.add(SWIR1));
    
    var BRIGTSOIL=((BLU.lt(0.168)).and(growing15)).or((BLU.lt(0.168)).and(growing14).and(  ((NIR.subtract(SWIR1)).gt(-0.079)))); 
    var VIIRSAREA = imageVIIRS.select('avg_rad').gt(10);
    var VIIRSAREA_VOLCAN = imageVIIRS.select('avg_rad').lt(10);
    
    var WATERSHAPE= ((BLU.subtract(GREEN)).gt(-0.2)).and(decreasing2345).and(WETNESS.gt(0)); //add other cond
    var OTHERWATERSHAPE= (BLUgteGREEN).and(GREENgteRED).and(NIR.gte(RED)).and(SWIR1.lt(NIR)).and(SWIR2.lte(SWIR1)).and(NIR.lt((RED).multiply(1.3)).and(NIR.lt(0.12)).and(SWIR1.lt(RED)).and(NIR.lte(GREEN)).and(NIR.gt(0.039)).and(WETNESS.gt(0))  ); //add other cond  07/10 (add replaced with and  :) and(NIR.lte(GREEN))
    
    var SNOWSHAPE=(min1234.gt(0.30)).and(NDSI.gt(0.65));
    
    var CLOUDSHAPE = ((SNOWSHAPE.eq(0)).and(BRIGTSOIL.eq(0))).and(
                  ((max123457.gt(0.47)).and(min1234.gt(0.37))).or(
                    
                    ((min123.gt(0.17)).and((SWIR1).gt(min123))).and(
                          ((SATURATION.gte(0.2)).and(SATURATION.lte(0.4)).and(max234.gte(0.35)) ).or ((NDSI.lt(0.65)).and(max1234.gt(0.30)).and( (NIR.divide(RED)).gte(1.3) ).and((NIR.divide(GREEN)).gte(1.3)).and( (NIR.divide(SWIR1)).gte(0.95)  )) 
                                                                   )
                                                                   
                                                                  ) 
                                                              ) 
    
    min123=0
    
    BRIGTSOIL=0
    SATURATION=0
    decreasing2345=0
    var th_NDVI_river = 0.13 //river
    // main groups based on ndvi
    var ndvi_1 = NDVI.lte(th_NDVI_MAX_WATER);
    var ndvi_2 = NDVI.lt(th_NDVI_MIN_VEGE).and(ndvi_1.eq(0)); 
    var ndvi_3 = NDVI.gte(th_NDVI_MIN_VEGE);
    var ndvi_4 = NDVI.lte(th_NDVI_river)
    
    
    //NDBI = (SWIR2 – NIR)/(SWIR2 + NIR)
    var NDBI = (SWIR2.subtract(NIR)).divide(SWIR2.add(NIR));
    var bareThreshold = -0.13;
    var waterNdbiThreshold = -0.1;
    var ndbi_th = NDBI.gt(bareThreshold)
    var ndbi_Water_th = NDBI.lte(waterNdbiThreshold);
    var waterNdbiThreshold2 = 0;
    var ndbi_Water_th2 = NDBI.gte(waterNdbiThreshold2);
    var waterNdbiThreshold3 = 0.2;
    var ndbi_Water_th3 = NDBI.lte(waterNdbiThreshold3);
    
    //NDTI
    var NDTI = (SWIR1.subtract(SWIR2)).divide(SWIR1.add(SWIR2));
    var ndtiUrbanThreshold = 0.09;
    var ndti_Urban_th =NDTI.gte(ndtiUrbanThreshold);
    var ndtiUrbanThreshold2 = 0.11;
    var ndti_Urban_th2 =NDTI.lte(ndtiUrbanThreshold2);
    
    
    
    //UI = (SWIR2 – NIR)/(SWIR2 + NIR)
    var UI = (SWIR2.subtract(NIR)).divide(SWIR2.add(NIR));
    var urbanThreshold = -0.05;
    var ui_th = UI.gt(urbanThreshold);
 
    
    //DBSI
    var DB = (SWIR1.subtract(GREEN)).divide(SWIR1.add(GREEN));
    
    
    ///NDWI
    var NDWI = image.normalizedDifference(['B3', 'B8']) 
    
    var ndwiThreshold = -0.1 ;
    var ndwi_th = NDWI.gte(ndwiThreshold);
    var ndwiThreshold2 = -0.3;//river
    var ndwi_th2 = NDWI.lte(ndwiThreshold2);
    var ndwiSoilsThreshold = -0.16;
    var ndwi_Soils_th = NDWI.lte(ndwiSoilsThreshold);
    var ndwiUrbanThreshold = -0.2
    var ndwi_Urban_th = NDWI.gt(ndwiUrbanThreshold);
    var ndwiUrbanThreshold1 = -0.08
    var ndwi_Urban_th1 = NDWI.lte(ndwiUrbanThreshold1);
     var ndwiVulcanThreshold = -0.15
    var ndwi_Vulcan_th = NDWI.gt(ndwiVulcanThreshold);
    var ndwiVulcanThreshold1 = 0.11
    var ndwi_Vulcan_th1 = NDWI.lte(ndwiVulcanThreshold1);

    
      ///---------------------------------------BSI Index------------------------------------------------------//
    var BSI = ((RED.add(SWIR1)).subtract((NIR.add(BLU)))).divide((RED.add(SWIR1)).add((NIR.add(BLU))));
    var bsiUrbanThreshold1 = 0.149;
    var bsi_Urban_th1 = BSI.lte(bsiUrbanThreshold1);
    var bsiUrbanThreshold2 = 0;
    var bsi_Urban_th2 = BSI.gt(bsiUrbanThreshold2);
    var bsiSoilThreshold1 = 0.149;//0.12
    var bsi_Soil_th1 = BSI.gt(bsiSoilThreshold1);
    var bsiWaterThreshold = 0;
    var bsi_Water_Th = BSI.lte(bsiWaterThreshold)
    var bsiWaterThreshold2 = 0.5; //river
    var bsi_Water_Th2 = BSI.lte(bsiWaterThreshold2)
    var bsiWaterThreshold3 = 0.15; //river
    var bsi_Water_Th3 = BSI.gte(bsiWaterThreshold3)
    var bsiWay1 = -0.1;
    var bsiWay2 = 0.1;
    var bsiWay_th1 = BSI.gt(bsiWay1);
    var bsiWay_th2 = BSI.lt(bsiWay2);
    
    ///---------------------------------------GNDVI Index------------------------------------------------------//
    //GNDVI = (B8-B3/B8+B3)
    var GNDVI = (NIR.subtract(GREEN)).divide(NIR.add(GREEN));
    var greenThreshold = 0.5;
    var gndvi_th = GNDVI.gte(greenThreshold)
    
      ///---------------------------------------EVI Index------------------------------------------------------//
    //EVI
    var EVI = image.expression(
     '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'RED': image.select('B4'),
      'BLUE': image.select('B2'),
      'NIR': image.select('B8') //calculate NDVI and add an NDVI band
      
    })
    
    ///---------------------------------------SAVI Index------------------------------------------------------//
    //SAVI
    var SAVI = image.expression(
        '((NIR-RED) / (RED+NIR+0.5))*(1.5)', {
      'RED': image.select('B4'),//b4
      'NIR': image.select('B8') //calculate SAVI and add an SAVI band b8
        
    })
    
    var SAVI1 = ((NIR.subtract(RED)).divide(RED.add(NIR).add(0.5))).multiply(1.5)
    
    var saviThreshold1 = -0.1;
    var savi_th1 = SAVI1.lte(saviThreshold1);//MinForest
    var saviThreshold2 = -0.3;
    var savi_th2 = SAVI1.gte(saviThreshold2);//MaxForest
      
    //-------------------------------------------------------------------------------------------------------------
		//----------------------  SECTION 1 : WATER  ---------------------------------------------------------
		//-------------------------------------------------------------------------------------------------------------
    
    OUT=(ndvi_1.and(SNOWSHAPE)).multiply(3);
    
    OUT=OUT.where(( (ndvi_1).and(BLU.gt(0.94)).and(GREEN.gt(0.94)).and(RED.gt(0.94)).and(NIR.gt(0.94)) ),1);  // TEST CLOUDS 
                   
		OUT=OUT.where(( (OUT.eq(0)).and(ndvi_4).and(ndbi_Water_th).and(bsi_Water_Th).and(ndwi_th)),2); //water
		
		OUT=OUT.where(( (OUT.eq(0)).and(ndvi_4).and(ndbi_Water_th2).and(ndbi_Water_th3).and(bsi_Water_Th2).and(bsi_Water_Th3).and(ndwi_th2)),2); //river
    
    
    //-------------------------------------------------------------------------------------------------------------
		//---------------------  SECTION 2 : CLOUDS or SOIL  ---------------------------------------------------------
		//------------------------------------------------------------------------------------------------------------
    
     OUT=OUT.where(( (ndvi_2).and(SNOWSHAPE)),3);
     OUT=OUT.where(( (ndvi_2).and(OTHERWATERSHAPE).and(BLU.gt(0.078)).and(max57.lt(0.058))),2 );//4
     OUT=OUT.where(( (ndvi_2).and(
                      CLOUDSHAPE.or(
                      (BLUgtGREENgtRED.and(NIR.gt(0.254)).and( BLU.gt(0.165)).and(NDVI.lt(0.40))).or(
                      (BLUgtGREEN.and(BLU.gt(0.27)).and(GREEN.gt(0.21)).and( REDsubtractGREEN.lte(0.1)).and(NIR.gt(0.35)))).or(
                      (BLU.gt(0.94)).and(GREEN.gt(0.94)).and(RED.gt(0.94)).and(NIR.gt(0.94))))
                      
                      )),1);

    CLOUDSHAPE=0
    WETNESS=0
   
     OUT=OUT.where( (OUT.eq(0)).and(NDWI.gt(-0.20)).and(NDWI.lte(-0.09)).and(VIIRSAREA) ,14 );//urban***** //.and(EVI.lt(0.09))//ndvi 0.07
     OUT= OUT.where((NDVI.gt(-0.0001)).and(NDVI.lte(0.075)),7);
     OUT= OUT.where( (NDWI.gt(-0.08)).and(NDWI.lte(0.05)),7);
     
     OUT=OUT.where(( (ndvi_2).and(BLU.lt(0.13)).and(BLUgtGREENgtRED).and(RED.lt(0.09)).and( BLUsubtractNIR.lt(-0.04))  ),5   );//DARK_VEG
     
     OUT=OUT.where(( (OUT.eq(0)).and(ndvi_2).and(
                    ((BLU.lt(0.17)).and(BLU.gt(0.133)).and(BLUgtGREENgtRED).and(RED.lt(0.17)).and(NIR.lte(0.27)).and( ((NIR).subtract(BLU)).lte(0.06))).or(
                    ( ((((NIR.subtract(GREEN)).abs().lte(0.06)).add( BLUsubtractNIR.gte(-0.06))).gt(0)).and(BLUgtGREENgtRED).and(NIR.gte(0.22)) )).or(
                    ( (OUT.eq(0)).and(ndvi_2).and(NDVI.lte(0.24)).and(NIR.lte(0.23)).and(GREENlteRED).and(REDlteNIR)) )
                   )),6);//DARK_SOIL
  
     OUT=OUT.where(( (OUT.eq(0)).and(ndvi_2).and(NDVI.lte(0.19)).and(NIR.gt(0.330)).and(growing14)   ),7 );//BRIGHT_SOIL

     OUT=OUT.where(( (OUT.eq(0)).and(ndvi_2).and(NDVI.gte(0.51)).and(BLUgteGREEN).and(REDsubtractGREEN.lt(0.04)).and(NDWI.lt(0.05))    ),8 );//GRASS
     
     OUT=OUT.where(( (OUT.eq(0)).and(ndvi_2).and(NDVI.gte(0.20)).and( REDsubtractGREEN.lt(0.05))   ),11 );//SPARSE
     
     OUT=OUT.where(( (OUT.eq(0)).and(ndvi_2)),12);//SOIL

     
     REDsubtractGREEN=0
     BLUgteGREEN=0
     
    //-------------------------------------------------------------------------------------------------------------
		//----------------------  SECTION 3 : VEGETATION  -------------------------------------------------------------
		//-------------------------------------------------------------------------------------------------------------
    OUT=OUT.where( (NDVI.gte(0.55)).and(OUT.eq(0)).and(NDSI.gt(-0.25))  ,10);//PURE_FOREST
    
    var MyCOND=(ndvi_3).and(NDVI.lt(th_RANGELAND));
    OUT=OUT.where(( (MyCOND).and(NIR.gte(0.15)).and(NDWI.lt(0.05)) ),8); //GRASS 
    OUT=OUT.where(( (MyCOND).and(NIR.lt(0.15))  ),4); //DARK_VEG

    MyCOND=(ndvi_3).and(NDVI.lt(th_GRASS));
    OUT=OUT.where(( (MyCOND).and(BLUlteNIR).and(NIR.lt(0.15))  ),4);//DARK_VEG
    OUT=OUT.where(( (OUT.eq(0)).and(
                                    ((MyCOND).and(BLUlteNIR)).or( (NDVI.lt(th_SHRUB) ).and(NIR.gt(0.22)))).and(NDSI.lt(-0.35))  ),16);//SHRUB
  
                                    
    OUT=OUT.where(( (MyCOND).and(BLU.gt(NIR)) ),4);//DARK_VEG
    OUT=OUT.where(( (OUT.eq(0)).and(MyCOND)).and(NDSI.lt(-0.3)),16); //SHRUB 
    
    OUT=OUT.where(((OUT.eq(0)).and(ndvi_3).and(NDVI.gt(th_TREES)) ),9);//BRIGHT_FOREST
    
    OUT=OUT.where(( (OUT.eq(0)).and(NDVI.lt(th_GRASS))).and(NDWI.lt(0.05)),8);//GRASS
    OUT=OUT.where(( (OUT.eq(0)).and(NDSI.lt(-0.25))),13);//DEGRADED_FOREST
    OUT=OUT.where(((OUT.eq(0)).and(NDWI).gte(0.05)),2);//water
    
    //Agriculture
    OUT = OUT.where((SAVI1.gte(0.10).and (SAVI1.lte(0.25)) ).and((GNDVI.gte(0.30)).and ((GNDVI.lte(0.45))) ),15);
    
    ////****spetial conditions***///PARA LOS VIVEROS
    var COND_VI1 =  (NDVI.gte(0.29)).and(NDVI.lte(0.43));//0.43 
    var COND_VI2 =  (EVI.gte(0.36)).and(EVI.lte(0.48));
    var COND_VI3 =  (SAVI.gte(0.2)).and(SAVI.lte(0.35));
    var COND_VI6 =  (NDWI.gte(-0.50)).and(NDWI.lte(-0.30));
    
    OUT=OUT.where( COND_VI1.and(COND_VI3).and(COND_VI6).and(COND_VI2) ,17);//PURE_FOREST
    
    OUT=OUT.where(( (OUT.eq(7)).and(EVI.lt(0.07))).and(VIIRSAREA),14);//urban  
    OUT= OUT.where((NDVI.gt(-0.0001)).and(NDVI.lte(0.075)).and(VIIRSAREA_VOLCAN),7);//soil
    OUT= OUT.where((NDWI.gt(-0.08)).and(NDWI.lte(0.05)).and(VIIRSAREA_VOLCAN),7);//soil
    OUT=OUT.where( (OUT.eq(6).or(OUT.eq(12)).or(OUT.eq(14))).and(SAVI1.gt(0.25)).and(SAVI1.lte(0.70)),16);//Shrubs
     

    return (OUT.select([0],["Class"]).toByte());
    
   
   
}   // SINGLE DATE CLASSIFICATION



// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ------------------   PHENOLOGY BASED SYNTHESIS  ----------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------

// INPUT: collection of single band image classification (monthly SDC classification)
// OUT :  sinle image with Land Cover labels

exports.PBS_classification =  function (COLL){  

     var tot=COLL.map(function(res){return (res.gt(0)) }).sum().divide(100);
     var relTot=COLL.map(function(res){return (res.gte(2)) }).sum().divide(100);

     print('PBS classification...');
      
    //----------------------------------------------------------------------------------------------------------
    
     var CLOUDS=COLL.map(function(res){return res.eq(1)}).sum().divide(1);//Set clouds index
     tot=0;
     var DARK_VEG=COLL.map(function(res){return res.eq(4)}).sum().divide(relTot);
     var LOW_IL = COLL.map(function(res){return res.eq(5)}).sum().divide(relTot);
     
     var DARK_SOIL= COLL.map(function(res){return res.eq(6)}).sum().divide(relTot);
     var BRIGHT_SOIL=COLL.map(function(res){return res.eq(7)}).sum().divide(relTot);
     var SPARSE=COLL.map(function(res){return res.eq(11)}).sum().divide(relTot);
     var SOIL=COLL.map(function(res){return res.eq(12)}).sum().divide(relTot);
     
      var URBAN=COLL.map(function(res){return res.eq(14)}).sum().divide(relTot);

     var BRIGHT_FOREST=COLL.map(function(res){return res.eq(9)}).sum().divide(relTot);
     
     var DEGRADED_FOREST=COLL.map(function(res){return res.eq(13)}).sum().divide(relTot);
     var PURE_FOREST=COLL.map(function(res){return res.eq(10)}).sum().divide(relTot);
     var FOREST=PURE_FOREST.add(DEGRADED_FOREST)
     var SHRUB=COLL.map(function(res){return (res.eq(16))}).sum().divide(relTot);
     
     var GRASS= COLL.map(function(res){return res.eq(8)}).sum().divide(relTot);
     
     var EVERGREEN_FOREST= (BRIGHT_FOREST.add(FOREST).add(DARK_VEG));
     
     var CROP = COLL.map(function(res){return res.eq(15)}).sum().divide(relTot);
     
     var WATER=COLL.map(function(res){return res.eq(2)}).sum().divide(relTot);
     var SNOW=COLL.map(function(res){return res.eq(3)}).sum().divide(relTot);
     var GREENHOUSE =COLL.map(function(res){return res.eq(17)}).sum().divide(relTot);
     
     var DRY = SPARSE.add(SOIL).add(BRIGHT_SOIL).add(DARK_SOIL).add(LOW_IL);
      
     print('Valiables are now set ... ');
      
     var OUT =relTot.multiply(0).select([0],["Class"]).toByte(); 
     
     
     OUT=OUT.where( (CLOUDS.gt(9)),10); //no data
     
     OUT=OUT.where( (URBAN.gte(80))  ,8);  // urban
     
     OUT=OUT.where( (WATER.gte(80))  ,1);  // water
     
     OUT=OUT.where( (FOREST.gte(80))  ,2);  // FOREST
     
     OUT=OUT.where( (SNOW.gte(80))  ,9);  // snow
     
     var GREEN=(EVERGREEN_FOREST.add(SHRUB).add(GRASS));
     
     var BARE= (BRIGHT_SOIL.add(SOIL).add(DARK_SOIL));
     OUT = OUT.where( (GREEN.eq(0)).and(CROP.eq(0)).and(WATER.eq(0)).and(BARE.gt(5)).and(URBAN.gte(50)), 8 ); //urban
    
     OUT = OUT.where( (BARE.gt(90)).and(URBAN.lt(10)), 7 ); //baresoils
     
     OUT = OUT.where( (BARE.gt(70)).and(URBAN.gt(10)), 7 ); //baresoils
     
     OUT = OUT.where( (BARE.gt(50)).and(URBAN.gt(20)), 7 ); //baresoils
     
     OUT = OUT.where( (BARE.gt(40)).and(URBAN.gt(15)), 7 ); //baresoils
     
     OUT = OUT.where( (BARE.gt(80)).and(WATER.gt(15)), 7 ); //baresoils
     
     OUT = OUT.where( (BARE.gt(20)).and(URBAN.gt(50)), 8 ); //urban
    
     OUT = OUT.where( (BARE.add(URBAN).gt(65)).and ((GREEN.gt(10)).and(GREEN.lt(30))), 7 ); //baresoil
     
     OUT=OUT.where( (GREENHOUSE.gte(80))  ,3);  // FOREST
     
     OUT = OUT.where( (CLOUDS.gt(9)), 10 ); //nodata
     
     ////////**********volcano*******////////////////
     OUT=OUT.where( (SNOW.gt(10))  ,9);  // snow
 
     var COND_VULCANO1 = WATER.gte(10);
     var COND_VULCANO2 = WATER.lte(60);
     var COND_VULCANO3 = BARE.gte(20);
     var COND_VULCANO4 = BARE.lte(80);
     var COND_VULCANO5 = WATER.lte(80);
     var COND_VULCANO6 = BARE.lte(40);
     
     OUT = OUT.where( (OUT.eq(0)).and(COND_VULCANO1).and(COND_VULCANO2).and(COND_VULCANO3).and(COND_VULCANO4).and(GREEN.eq(0)), 7 ); //soil vulcano
     OUT = OUT.where( (OUT.eq(0)).and(COND_VULCANO1).and(COND_VULCANO5).and(COND_VULCANO3).and(COND_VULCANO6).and(GREEN.lt(20)), 1 ); //water soil
     
    
     OUT=OUT.where( (OUT.eq(0)).and(GREEN.gt(65)).and((GRASS.add(SHRUB)).gt(20)).and(DRY.lt(35)) ,5);  // DECIDE IF OPEN DRY or HUMID
    
     OUT=OUT.where( ((GREEN.gt(90)).and(DRY.lt(5))).and(FOREST.lt(50)),2);//forest
     
     OUT=OUT.where( ((GREEN.gt(90)).and(DRY.lt(5))).and(FOREST.gt(70)),2);//forest
    
     OUT=OUT.where( (GREEN.gt(70)).and(SHRUB.gt(20)).and(FOREST.gt(50)).and((GRASS.add(SHRUB)).gt(30)).and(WATER.lt(5)).and(URBAN.eq(0))   ,2);//forest 
     
     
     OUT=OUT.where( (GREEN.gt(85)).and((SHRUB.add(DEGRADED_FOREST)).gt(70)).and(DRY.lt(5)) ,2);
     
     OUT=OUT.where(  (GREEN.gte(95)).and(FOREST.gte(30)).and(SHRUB.gt(40)).and(EVERGREEN_FOREST.gt(40)).and(WATER.lt(5)).and(URBAN.eq(0))   ,2);  // almost evergreen   
     
     OUT=OUT.where( (FOREST.gt(60)).and(EVERGREEN_FOREST.gt(80)).and(WATER.lt(5)).and(URBAN.eq(0)),2);
     OUT=OUT.where( (FOREST.gt(70)).and(GREEN.eq(FOREST)).and(WATER.lt(5)).and(URBAN.eq(0)),2);
     OUT=OUT.where( (FOREST.gt(50)).and(GREEN.gt(60)).and(((DARK_VEG).add(LOW_IL).add(GREEN)).gt(90)).and(WATER.lt(5)).and(URBAN.eq(0)),2);
     
     OUT=OUT.where( (FOREST.gt(50)).and((EVERGREEN_FOREST.add(SHRUB)).gt(80)).and(DRY.eq(0)).and(WATER.lt(5)).and(URBAN.eq(0)),2);
     
     /////paramo conditions
     var condition1 = (GREEN.gte(10).and(GREEN.lt(85))) ;
     var condition2 = (CROP.gte(0)) ;
     var condition3 = (GRASS.gte(10).and(GRASS.lt(85))) ;
     var condition4 = (SPARSE.gte(10).and(SPARSE.lt(50))) ;
     var condition5 = ((DRY.add(URBAN)).gte(0)).and((DRY.add(URBAN)).lt(60))
     var condition6 = ((DEGRADED_FOREST.add(SPARSE)).gte(10)).and((DEGRADED_FOREST.add(SPARSE)).lt(100))
     
     OUT=OUT.where( (OUT.eq(0)).and(condition1).and(condition2).and(condition3).and(WATER.lt(20)) ,5);  // paramo vegetation-grassland
     OUT=OUT.where( (OUT.eq(0)).and(GREEN.eq(0)).and(DRY.eq(0)).and(WATER.eq(0)).and(URBAN.gte(0).and(condition2)) ,5);  // paramo vegetation-grassland
     OUT=OUT.where( (OUT.eq(0)).and(GREEN.eq(0)).and(WATER.lt(20)).and(URBAN.lte(30).and(condition2).and(condition4)) ,5);  // paramo vegetation-grassland
     OUT=OUT.where((GRASS.gte(70)).and(DRY.eq(0)).and(FOREST.eq(0)).and(URBAN.eq(0)).and(CROP.eq(0)),5)// // paramo vegetation-grassland
     OUT=OUT.where( (OUT.eq(0)).and(condition6).and(condition5).and(condition2).and(PURE_FOREST.eq(0)).and(WATER.lte(25)),5)// // paramo vegetation-grassland
     
     
     OUT=OUT.where( (OUT.eq(0)).and(GREEN.gt(0)).and(FOREST.lt(20)).and(condition4).and(condition5).and(condition6) ,6);  // paramo vegetation-shrubs
     OUT=OUT.where( (OUT.eq(0)).and(GREEN.eq(0)).and(DRY.eq(0)).and(FOREST.eq(0)).and(GRASS.eq(0)).and(WATER.lte(40)).and(condition2) ,6);  // paramo vegetation-shrubs 
     OUT=OUT.where( (OUT.eq(0)).and(GREEN.gt(0)).and(FOREST.lt(20)).and(condition2).and(condition5).and(condition6) ,6);  // paramo vegetation-shrubs
     OUT=OUT.where( (OUT.eq(0)).and(DRY.lt(20)).and(FOREST.eq(0)).and(GRASS.eq(0)).and(SHRUB.gte(10)).and(URBAN.lte(50)).and(CROP.gt(10)) ,6);  // paramo vegetation-shrubs
     OUT=OUT.where( (OUT.eq(0)).and(condition5).and(FOREST.eq(0)).and(GRASS.eq(0)).and(SHRUB.eq(0)).and(WATER.eq(0)).and(CROP.gt(10)) ,6);  // paramo vegetation-shrubs
     OUT=OUT.where( (OUT.eq(0)).and(((DRY.add(URBAN)).gte(0)).and((DRY.add(URBAN)).lt(90))).and(FOREST.eq(0)).and(GRASS.eq(0)).and(SHRUB.eq(0)).and(WATER.eq(0)).and(CROP.gt(10)) ,6);  // paramo vegetation-shrubs
      
     OUT=OUT.where((OUT.eq(0)).and(FOREST.eq(0)).and(GREEN.eq(0)).and((DRY.add(URBAN)).gte(50)).and(CROP.gte(20)),7);//paramo soil
     
     
     OUT=OUT.where((OUT.eq(0)).and(FOREST.gte(40)).and(GREEN.gte(40)).and(condition2).and(condition6).and(condition5),2);//paramo forest
     OUT=OUT.where( (PURE_FOREST.gte(40)).and((GRASS.add(SHRUB)).eq(0)).and(DEGRADED_FOREST.eq(0)).and(WATER.lt(5)).and((DRY.add(URBAN)).lte(30)), 2);//paramo forest   ,2);//forest 
     
     
    ///////////////*************************Annual Crop********************************///////////////
     var condition1AC = (CROP.gt(0)).and(CROP.lt(80));
     var condition2AC = ((DRY.add(URBAN)).gt(0)).and(((DRY.add(URBAN)).lte(50)));
     var condition3AC = ((SHRUB.add(GRASS)).gt(0)).and(((SHRUB.add(GRASS)).lte(60)));//*/*****
     var condition4AC = (DEGRADED_FOREST.gt(0)).and(DEGRADED_FOREST.lte(80))
     
     OUT=OUT.where((CLOUDS.lte(5)).and(condition1AC).and(condition2AC).and(condition3AC).and(condition4AC),3);
     OUT=OUT.where((CLOUDS.lte(5)).and(CROP.eq(0)).and(condition2AC).and(condition3AC).and(condition4AC),3);
     OUT=OUT.where((CLOUDS.lte(5)).and(condition3AC).and(condition4AC).and(PURE_FOREST.eq(0)),3);
     OUT=OUT.where((CLOUDS.lte(11)).and(GREENHOUSE.gte(10)).and((SPARSE.add(URBAN)).lt(70)).and(GRASS.lt(70))
                    .and(FOREST.eq(0)).and(SHRUB.eq(0))
                    .and(CROP.eq(0)).and(BARE.eq(0)).and(WATER.eq(0)),11);//3 floricultura 
             
     
     ////////////////*******************Grassland - Pastos*******************////////////////////////////
     
     OUT=OUT.where((CROP.eq(0)).and(DEGRADED_FOREST.gte(70)).and(condition3AC).and((DRY.add(URBAN)).eq(0)),5);
     OUT=OUT.where((CROP.eq(0)).and(DEGRADED_FOREST.gte(50)).and((FOREST.gte(50)).and(FOREST.lte(90))).and(condition3AC).and((DRY.add(URBAN)).eq(0)),5);
     OUT=OUT.where((CROP.eq(0)).and(DEGRADED_FOREST.gte(50)).and((FOREST.gte(50)).and(FOREST.lte(90)))
                      .and(((SHRUB.and(GRASS)).gte(0)).and(((SHRUB.and(GRASS)).lte(50)))).and((DRY.add(URBAN)).eq(0)),5);
     
      
      ////////////////*******************Forest - Shrubs*******************////////////////////////////
     
     var condition1FS = ((SPARSE.add(BARE)).gt(0));
     var condition2FS = (((SPARSE.add(BARE)).lte(70)));
     var condition3FS = ((FOREST.add(CROP)).gt(0));
     var condition4FS = (((FOREST.add(CROP)).lte(70)));
       
     OUT=OUT.where((OUT.eq(0)).and(condition1FS).and(condition2FS)
                            .and(condition3FS).and(condition4FS).and((DRY.gte(20))),6); 
      
    
     OUT=OUT.where((OUT.eq(0)),5);
      
    
     return OUT
}  //  END PBS


