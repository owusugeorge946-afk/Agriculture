// ============================================================
// FIGURE 6: ASI + HOTSPOT PANEL — LONGER CLEAN LAYOUT
// ============================================================

// -------------------- 1. LOAD STUDY AREA --------------------
var ghana = ee.FeatureCollection("projects/ee-owusugeorge946/assets/Ghana");
var region = ghana.geometry().bounds().buffer(50000);

// -------------------- 2. DATASETS --------------------
var chirps = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY")
  .filterBounds(ghana)
  .filterDate("2000-01-01", "2024-12-31")
  .select("precipitation");

var ndviCol = ee.ImageCollection("MODIS/061/MOD13Q1")
  .filterBounds(ghana)
  .filterDate("2000-01-01", "2024-12-31")
  .select("NDVI");

// -------------------- 3. GSL FUNCTION --------------------
function calculateGSL(year) {
  year = ee.Number(year);

  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year, 12, 31);

  var yearData = chirps.filterDate(start, end);

  var withDOY = yearData.map(function(img) {
    var doy = ee.Date(img.get("system:time_start"))
      .getRelative("day", "year")
      .add(1);

    var doyImg = ee.Image.constant(doy).rename("DOY").toFloat();

    return img.toFloat()
      .addBands(doyImg)
      .clipToCollection(ghana)
      .copyProperties(img, ["system:time_start"]);
  });

  var wetDOY = withDOY.map(function(img) {
    var wetMask = img.select("precipitation")
      .gt(5)
      .selfMask()
      .toFloat();

    return img.select("DOY")
      .updateMask(wetMask)
      .rename("wetDOY")
      .toFloat();
  });

  var onset = wetDOY.min().rename("onset").toFloat();
  var cessation = wetDOY.max().rename("cessation").toFloat();

  return cessation.subtract(onset)
    .rename("GSL")
    .toFloat()
    .clipToCollection(ghana);
}

// -------------------- 4. MEAN GSL --------------------
var years = ee.List.sequence(2000, 2024);

var gslCollection = ee.ImageCollection.fromImages(
  years.map(function(y) {
    return calculateGSL(y);
  })
);

var meanGSL = gslCollection.mean()
  .rename("GSL")
  .toFloat()
  .clipToCollection(ghana);

// -------------------- 5. MEAN NDVI --------------------
var meanNDVI = ndviCol.mean()
  .multiply(0.0001)
  .rename("NDVI")
  .toFloat()
  .clipToCollection(ghana);

// -------------------- 6. NORMALIZATION --------------------
function normalize(img, bandName) {
  var stats = img.reduceRegion({
    reducer: ee.Reducer.minMax(),
    geometry: ghana.geometry(),
    scale: 5000,
    maxPixels: 1e13,
    bestEffort: true
  });

  var minVal = ee.Number(stats.get(bandName + "_min"));
  var maxVal = ee.Number(stats.get(bandName + "_max"));

  return img.subtract(minVal)
    .divide(maxVal.subtract(minVal))
    .rename(bandName + "_norm")
    .toFloat()
    .clipToCollection(ghana);
}

var gslNorm = normalize(meanGSL, "GSL");
var ndviNorm = normalize(meanNDVI, "NDVI");

// -------------------- 7. ASI --------------------
var ASI = gslNorm.add(ndviNorm)
  .divide(2)
  .rename("ASI")
  .toFloat()
  .clipToCollection(ghana);

// -------------------- 8. HOTSPOT --------------------
var kernel = ee.Kernel.square({
  radius: 50000,
  units: "meters",
  normalize: false
});

var localMean = ASI.reduceNeighborhood({
  reducer: ee.Reducer.mean(),
  kernel: kernel
}).rename("Local_ASI_Mean");

var hotspot = ASI.subtract(localMean)
  .rename("ASI_Hotspot")
  .toFloat()
  .clipToCollection(ghana);

// -------------------- 9. VISUALIZATION --------------------
var asiVis = {
  min: 0,
  max: 1,
  palette: ["#d7191c", "#fdae61", "#ffffbf", "#a6d96a", "#1a9641"]
};

var hotspotVis = {
  min: -0.15,
  max: 0.15,
  palette: ["#2166ac", "#f7f7f7", "#b2182b"]
};

var boundary = ghana.style({
  color: "000000",
  fillColor: "00000000",
  width: 2
});

var white = ee.Image.constant(1).visualize({
  min: 0,
  max: 1,
  palette: ["ffffff"]
});

function makeComposite(img, vis) {
  return ee.ImageCollection([
    white,
    img.visualize(vis),
    boundary
  ]).mosaic();
}

var asiComposite = makeComposite(ASI, asiVis);
var hotspotComposite = makeComposite(hotspot, hotspotVis);

// -------------------- 10. PANEL FUNCTION --------------------
function makePanel(composite, label) {

  var thumb = ui.Thumbnail({
    image: composite,
    params: {
      region: region,
      dimensions: "480x380",
      format: "png"
    },
    style: {
      width: "480px",
      height: "380px",
      margin: "0px",
      padding: "0px"
    }
  });

  var panelTitle = ui.Label(label, {
    fontWeight: "bold",
    fontSize: "12px",
    color: "black",
    margin: "2px 0 4px 4px"
  });

  return ui.Panel({
    widgets: [panelTitle, thumb],
    layout: ui.Panel.Layout.Flow("vertical"),
    style: {
      width: "505px",
      height: "420px",
      border: "1px solid #777",
      margin: "3px",
      padding: "5px",
      backgroundColor: "white"
    }
  });
}

var panelASI = makePanel(asiComposite, "(a) Agricultural Sustainability Index");
var panelHot = makePanel(hotspotComposite, "(b) ASI Hotspot / Coldspot");

// -------------------- 11. LEGENDS --------------------
function makeLegend(titleText, palette, leftText, midText, rightText) {

  var legendTitle = ui.Label(titleText, {
    fontWeight: "bold",
    fontSize: "11px",
    margin: "2px 0 3px 4px"
  });

  var legendBar = ui.Thumbnail({
    image: ee.Image.pixelLonLat().select(0),
    params: {
      bbox: [0, 0, 1, 0.1],
      dimensions: "480x16",
      format: "png",
      min: 0,
      max: 1,
      palette: palette
    },
    style: {
      width: "480px",
      height: "16px",
      margin: "0 4px"
    }
  });

  var legendLabels = ui.Panel({
    layout: ui.Panel.Layout.Flow("horizontal"),
    style: {
      width: "480px",
      margin: "2px 4px 2px 4px"
    }
  });

  legendLabels.add(ui.Label(leftText, {
    fontSize: "9px",
    textAlign: "left",
    stretch: "horizontal"
  }));

  legendLabels.add(ui.Label(midText, {
    fontSize: "9px",
    textAlign: "center",
    stretch: "horizontal"
  }));

  legendLabels.add(ui.Label(rightText, {
    fontSize: "9px",
    textAlign: "right",
    stretch: "horizontal"
  }));

  return ui.Panel({
    widgets: [legendTitle, legendBar, legendLabels],
    layout: ui.Panel.Layout.Flow("vertical"),
    style: {
      width: "505px",
      backgroundColor: "white",
      margin: "3px"
    }
  });
}

var asiLegend = makeLegend(
  "ASI: Low to High Sustainability",
  asiVis.palette,
  "Low",
  "Moderate",
  "High"
);

var hotLegend = makeLegend(
  "Hotspot Pattern",
  hotspotVis.palette,
  "Coldspot / Vulnerable",
  "Neutral",
  "Hotspot / Resilient"
);

// -------------------- 12. FINAL LAYOUT --------------------
var title = ui.Label("Agricultural Sustainability Patterns (ASI)", {
  fontWeight: "bold",
  fontSize: "15px",
  color: "black",
  margin: "4px 0 4px 4px"
});

var mapRow = ui.Panel({
  widgets: [panelASI, panelHot],
  layout: ui.Panel.Layout.Flow("horizontal"),
  style: {
    backgroundColor: "white",
    margin: "0px",
    padding: "0px"
  }
});

var legendRow = ui.Panel({
  widgets: [asiLegend, hotLegend],
  layout: ui.Panel.Layout.Flow("horizontal"),
  style: {
    backgroundColor: "white",
    margin: "0px",
    padding: "0px"
  }
});

var mainPanel = ui.Panel({
  widgets: [title, mapRow, legendRow],
  layout: ui.Panel.Layout.Flow("vertical"),
  style: {
    width: "1040px",
    backgroundColor: "white",
    padding: "4px"
  }
});

ui.root.clear();
ui.root.add(mainPanel);

// -------------------- 13. PRINT STATISTICS --------------------
print("ASI summary statistics:", ASI.reduceRegion({
  reducer: ee.Reducer.mean()
    .combine({reducer2: ee.Reducer.minMax(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
  geometry: ghana.geometry(),
  scale: 5000,
  maxPixels: 1e13,
  bestEffort: true
}));

print("Hotspot summary statistics:", hotspot.reduceRegion({
  reducer: ee.Reducer.mean()
    .combine({reducer2: ee.Reducer.minMax(), sharedInputs: true})
    .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true}),
  geometry: ghana.geometry(),
  scale: 5000,
  maxPixels: 1e13,
  bestEffort: true
}));
