// Define start and end ee.Dates.
var startDate = ee.Date('2000-01-01');
var endDate = ee.Date('2010-01-01');


// Define the regional bounds 
var region = ee.Geometry.Polygon(
        [[[139.81147500927887, -22.584034738321385],
          [139.81147500927887, -39.52984428265981],
          [155.19233438427887, -39.52984428265981],
          [155.19233438427887, -22.584034738321385]]], null, false);

// Subset the to a Country .
var congo = ee.Feature(
  region
);

// Load a FeatureCollection from a table dataset: 'RESOLVE' ecoregions.
var ecoregions = ee.FeatureCollection('RESOLVE/ECOREGIONS/2017');


// Subset o the bounds of the ecoregion feature
// and other criteria. Clip to the intersection with congo.
var protectedAreas = ecoregions
  .filter(ee.Filter.and(
    ee.Filter.bounds(congo.geometry()),
    ee.Filter.eq('BIOME_NUM', 4)             // 4 = TBMF AUstralia, 1=TMBF Amazonia
  ))
  .map(function(feat){
    return congo.intersection(feat,ee.ErrorMargin(1));
  });

// Map.addLayer(protectedAreas, {}, 'Eco region ');

var col =  ee.ImageCollection("MODIS/006/MCD64A1")
             .filterDate('2000-01-01', '2001-01-01')
             .select('BurnDate');

var fire =  ee.ImageCollection("MODIS/006/MCD64A1")
//var fire = ee.ImageCollection('FIRMS')
             .filterDate(startDate, endDate);

//print(fire)
// Clip and add a date band
var clipToRegion = function(img) {
  var clipped = img.clip(protectedAreas);
  return clipped;
};
             
var fire_clipped =fire.map(clipToRegion)

var scale = fire_clipped.first().projection().nominalScale().getInfo();
var crs = fire_clipped.first().projection().crs();

Map.addLayer(fire_clipped.first(), {}, 'Fire clipped');



var rasterFootprint = ee.Image(1).clip(protectedAreas)

Map.addLayer(protectedAreas, {}, 'Eco region ');
Map.addLayer(rasterFootprint, {}, 'Raster Footprint');

Map.centerObject(protectedAreas)
Export.image.toDrive({
    image: rasterFootprint,
//    description: 'FireCountsAustralia_rasterFootprint',
    description: 'BurnedAreaAustralia_rasterFootprint',
    scale: scale,
    region: region,
    crs: crs,
    fileFormat: 'GeoTIFF',
    formatOptions: {
      cloudOptimized: true
    }
  });
