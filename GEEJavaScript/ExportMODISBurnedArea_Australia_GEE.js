// Export Burned Area Australia 
// Export geotif images with 12 layer representing burned pixels by year

// Define start and end ee.Dates.
var startDate = ee.Date('2000-01-01');
var endDate = ee.Date('2021-01-01');

var col =  ee.ImageCollection("MODIS/006/MCD64A1")
             .filterDate(startDate, endDate )
             .select('BurnDate');

// Define a mask to clip the data.
var countries = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017');
  
var country = ee.Feature(
  countries
    .filter(ee.Filter.eq('country_na', 'Australia'))
    .first()
);
var ecoregions = ee.FeatureCollection('RESOLVE/ECOREGIONS/2017');

var protectedAreas = ecoregions
  .filter(ee.Filter.and(
    ee.Filter.bounds(country.geometry()),
    ee.Filter.eq('BIOME_NUM', 4)
  ))
  .map(function(feat){
    return country.intersection(feat);
  });
// Add Uganda outline to the Map as a layer.
Map.centerObject(protectedAreas, 5);
Map.addLayer(protectedAreas);  
//print('mask', mask)

// Define the regional bounds of animation frames.
var region = ee.Geometry.Polygon(
        [[[139.81147500927887, -22.584034738321385],
          [139.81147500927887, -39.52984428265981],
          [155.19233438427887, -39.52984428265981],
          [155.19233438427887, -22.584034738321385]]], null, false);

Map.addLayer(region); 

// create years and export for every year
var years = Array.apply(null, {length: 21}).map(Number.call, Number) // sequence of 16 numbers
            .map(function(number){
              return exportImagePerYear(col, number + 2000)}); // add 2000 for each year


var scale = col.first().projection().nominalScale().getInfo();
print('Scale', scale);

// function to export
function exportImagePerYear(col, startYear) {

  var colFilt = col.filterDate(String(startYear)+'-01-01', String(startYear+1)+'-01-01');

  // Filter fire with more than 50% confidence and add a new band representing areas where confidence of fire > 50%
  var clipToRegion = function(img) {
    var dateString = ee.Date(img.get('system:time_start')).format('yyyy-MM-dd');
    var clipped = img.clip(protectedAreas);
  //  var burned = clipped.gt(0);
    return clipped.rename(dateString);
  };
  var burned_area = colFilt.map(clipToRegion);

  //var check = ee.Image(burned_area.first());
  //Map.addLayer(check, {palette: ['000000', '00FFFF'], max: 366}, 'check');
  //print('Scale', scale )

  // Stack One layer by year 
  //
  var stackCollection = function(collection) {
    // Create an initial image.
    var first = ee.Image(collection.first()).select([]);

    // Write a function that appends a band to an image.
    var appendBands = function(image, previous) {
        var dateString = ee.Date(image.get('system:time_start')).format('yyyy-MM-dd');
        return ee.Image(previous).addBands(image);
    };
    return ee.Image(collection.iterate(appendBands, first));
  };

  var evi_img = stackCollection(burned_area);
  //print("EVI image stack",evi_img);

  Export.image.toDrive({
    image: evi_img,
    description: 'BurnedAreaAustralia'+String(startYear),
    scale: scale,
    region: region,
    fileFormat: 'GeoTIFF',
    formatOptions: {
      cloudOptimized: true
    }
  });
}

