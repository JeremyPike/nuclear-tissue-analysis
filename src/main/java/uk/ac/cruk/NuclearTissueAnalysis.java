
package uk.ac.cruk;

import java.awt.Rectangle;
import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.Future;

import org.scijava.command.Command;
import org.scijava.command.CommandModule;
import org.scijava.module.ModuleItem;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
import fiji.plugin.trackmate.features.FeatureFilter;
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.Roi;
import ij.measure.ResultsTable;
import ij.plugin.RoiEnlarger;
import ij.process.ImageStatistics;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;

/**
 * This is a ImageJ command for the ROI based analysis of a sub-populations of
 * cells from tissue sections. It was developed for the Narita lab at the
 * CRUK-CI. To use this command nuclei should first be segmented and classified
 * using an Ilastik pixel and object classification project.
 * 
 * The TrackMate Fiji plugin is used for spot detection: Tinevez, Jean-Yves, et
 * al.
 * "TrackMate: An open and extensible platform for single-particle tracking."
 * Methods (2016).
 * 
 * @author Jeremy Pike
 */

@Plugin(type = Command.class, menuPath = "Plugins>Nuclear Tissue Analysis")
public class NuclearTissueAnalysis<T extends RealType<T>> implements Command {

	@Parameter(label = "Select the directory containing the raw data tiles (.tifs files): ", persist = true, style = "directory")
	private File rawDataDir;

	@Parameter(label = "Select the directory containing the labelled object maps from islatik (.tif files): ", persist = true, style = "directory")
	private File objectMapsDir;

	@Parameter(label = "Radius for spot detection (microns): ", persist = false, min = "0.01", max = "20")
	private double spotRadius = 1.125;

	@Parameter(label = "Quality threshold for spot detection: ", persist = false, min = "0.01", max = "20")
	private double minSpotQuality = 40;

	@Parameter(label = "Define width of nucleus edge (microns): ", persist = false, min = "0.01", max = "20")
	private double nucleusEdgeWidth = 1.5;

	@Parameter(label = "Specify object label for cell class of interest: ", persist = false, min = "1")
	private int cellLabel = 1;

	@Parameter(label = "Specify the tile overlap percentage: ", persist = false, min = "0", max = "50")
	private double tileOverlapPercent = 10;

	@Override
	public void run() {

		// create a new Results table to hold the ROI based measurements
		ResultsTable rt = new ResultsTable();

		// list of all files in the directory containing the raw data (.tif
		// files)
		String[] rawDataDirList = rawDataDir.list();

		// configure the Trackmate detector settings
		Settings tmSettings = new Settings();
		tmSettings.detectorFactory = new LogDetectorFactory();
		tmSettings.addSpotFilter(new FeatureFilter("QUALITY", minSpotQuality, true));
		Map<String, Object> map = tmSettings.detectorFactory.getDefaultSettings();
		map.put("RADIUS", spotRadius);
		map.put("DO_MEDIAN_FILTERING", false);
		map.put("DO_SUBPIXEL_LOCALIZATION", true);
		// set to Dapi channel
		map.put("TARGET_CHANNEL", 1);
		map.put("THRESHOLD", 0.0d);
		tmSettings.detectorSettings = map;

		// loop through all files
		for (int i = 0; i < rawDataDirList.length; i++) {

			// if current file is a .tif and isn't a object map
			if (rawDataDirList[i].endsWith(".tif") && !rawDataDirList[i].contains("Object Predictions")) {

				System.out.println("file " + i + " of " + rawDataDirList.length);

				// open the raw image data
				ImagePlus rawSeriesImp = IJ.openImage(rawDataDir.getPath() + "\\" + rawDataDirList[i]);
				// open the corresponding object map which is produced by
				// Ilastik
				ImagePlus objectMapImp = IJ.openImage(objectMapsDir.getPath() + "\\"
						+ rawDataDirList[i].substring(0, rawDataDirList[i].length() - 4) + "_Object Predictions.tif");

				// produce a binary image from the object map where only objects
				// with the specified label are kept.
				IJ.setThreshold(objectMapImp, cellLabel, cellLabel);

				// find all connected components in the objectMap and output as
				// an Overlay
				IJ.run("Set Measurements...", "  redirect=None decimal=3");
				IJ.run(objectMapImp, "Analyze Particles...", "  show=[Overlay Outlines]");
				// get the overlay which contains an ROI for each connected
				// component
				Overlay overlay = objectMapImp.getOverlay();

				// if there is at least one object in the specified class
				if (overlay != null) {

					// set Trackmate to the raw image data
					tmSettings.setFrom(rawSeriesImp);
					// detect spots using Trackmate
					SpotCollection dapiSpots = countSpotsTrackmate(tmSettings, false);

					// find the X and Y bounds for which objects should be
					// included in the analysis. This is based on tile region
					// overlap percentage.
					double boundX = rawSeriesImp.getWidth() - rawSeriesImp.getWidth() * tileOverlapPercent / 100;
					double boundY = rawSeriesImp.getHeight() - rawSeriesImp.getHeight() * tileOverlapPercent / 100;

					// to store if each ROI is in the analysis region
					boolean[] roiCheck = new boolean[overlay.size()];
					// to store the spot counts for each ROI
					int[] roiSpotCountsAll = new int[overlay.size()];
					// to store the spot counts in the centre region of each ROI
					int[] roiSpotCountsCenter = new int[overlay.size()];

					for (int r = 0; r < overlay.size(); r++) {

						Roi roi = overlay.get(r);

						// get the bounding box of the current ROI
						Rectangle boundingBox = roi.getBounds();

						// shrink the ROU by a specified number of microns to
						// get the center region ROI
						Roi roiShrink = RoiEnlarger.enlarge(roi,
								-1 * nucleusEdgeWidth / rawSeriesImp.getCalibration().pixelWidth);

						// if lower bounds of bounding box are not outside of
						// the analysis region
						if (boundingBox.getMinX() != 0 && boundingBox.getMinY() != 0 && boundingBox.getMinX() <= boundX
								&& boundingBox.getMinY() <= boundY) {

							// include this object in analysis
							roiCheck[r] = true;

							// iterate though all spot detections in current
							// image
							Iterator<Spot> spotIterator = dapiSpots.iterator(true);
							while (spotIterator.hasNext()) {
								Spot spot = spotIterator.next();
								// get spot position in pixel coordinates
								int positionXPix = (int) Math.floor(
										spot.getFeature("POSITION_X") / rawSeriesImp.getCalibration().pixelWidth);
								int positionYPix = (int) Math.floor(
										spot.getFeature("POSITION_Y") / rawSeriesImp.getCalibration().pixelWidth);
								// if the spot is in the object ROI count
								if (roi.contains(positionXPix, positionYPix)) {
									roiSpotCountsAll[r]++;
								}
								// if the spot is in the object center ROI count
								if (roiShrink.contains(positionXPix, positionYPix)) {
									roiSpotCountsCenter[r]++;
								}
							}
						}
					}
					// to store area of each ROI
					double[] roiAreas = new double[overlay.size()];
					// to store intensity mean of each ROI in every channel
					double[][] roiMeans = new double[rawSeriesImp.getNChannels()][overlay.size()];
					// to store intensity standard deviation (std) of each ROI
					// in every channel
					double[][] roiStds = new double[rawSeriesImp.getNChannels()][overlay.size()];

					for (int c = 0; c < rawSeriesImp.getNChannels(); c++) {
						// set raw data ImagePlus to current channel
						rawSeriesImp.setC(c + 1);
						for (int r = 0; r < overlay.size(); r++) {

							// if current ROI in the analysis region
							if (roiCheck[r]) {
								// set the raw data to this ROI
								rawSeriesImp.setRoi(overlay.get(r));
								// measure some statistics for this ROI
								ImageStatistics statistics = rawSeriesImp.getStatistics();
								// store required measures
								roiAreas[r] = statistics.area;
								roiMeans[c][r] = statistics.mean;
								roiStds[c][r] = statistics.stdDev;

							}
						}

					}

					for (int r = 0; r < overlay.size(); r++) {
						// if current ROI in the analysis region
						if (roiCheck[r]) {
							// increment the output table
							rt.incrementCounter();
							// add the file name
							rt.addValue("Name", rawDataDirList[i]);
							// add all measures statistics
							rt.addValue("Area", roiAreas[r]);
							for (int c = 0; c < rawSeriesImp.getNChannels(); c++) {
								rt.addValue("Channel" + (c + 1) + " mean intensity", roiMeans[c][r]);
								rt.addValue("Channel" + (c + 1) + " std", roiStds[c][r]);
							}
							rt.addValue("Edge dapi spots", roiSpotCountsAll[r] - roiSpotCountsCenter[r]);
							rt.addValue("Center dapi spots", roiSpotCountsCenter[r]);
						}
					}
				}
			}
		}
		// show ResultsTable with name "Output"
		rt.show("Output");
	}

	/**
	 * This main function serves for development purposes.
	 *
	 * @param args
	 *            whatever, it's ignored
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {

		// create imageJ and show UI
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// get default setting input parameters for command and save in Map
		Map<String, Object> pluginSettings = new HashMap<String, Object>();
		Iterator<ModuleItem<?>> inputs = ij.command().getCommand(NuclearTissueAnalysis.class).inputs().iterator();
		while (inputs.hasNext()) {
			ModuleItem<?> input = inputs.next();
			if (input.getDefaultValue() != null)
				pluginSettings.put(input.getName(), input.getDefaultValue());
		}
		// specify where the raw data and object maps are stored
		pluginSettings.put("rawDataDir", new File(
				"V:\\core\\light_microscopy\\data\\group_folders\\2016_light_microscopy\\Jeremy\\Side_Projects\\Data\\Yoko_Ito\\2017-03-09_tiffs"));
		pluginSettings.put("objectMapsDir", new File(
				"V:\\core\\light_microscopy\\data\\group_folders\\2016_light_microscopy\\Jeremy\\Side_Projects\\Data\\Yoko_Ito\\2017-03-09_tiffs"));

		// run command with default settings
		Future<CommandModule> fc = ij.command().run(NuclearTissueAnalysis.class, true, pluginSettings);

	}

	/**
	 * Uses trackmate to detect spots
	 *
	 * @param settings
	 *            Settings for trackmate spot detection
	 * 
	 * @param display
	 *            if true displays the results
	 * 
	 * @return The collection of spots from trackmate
	 */
	public static SpotCollection countSpotsTrackmate(Settings settings, boolean display) {

		Model model = new fiji.plugin.trackmate.Model();
		TrackMate trackmate = new TrackMate(model, settings);
		// Check input is ok
		boolean ok = trackmate.checkInput();
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// Find spots
		ok = trackmate.execDetection();
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// Compute spot features
		ok = trackmate.computeSpotFeatures(true);
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}

		// Filter spots
		ok = trackmate.execSpotFiltering(true);
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// display spot detections
		if (display) {
			SelectionModel selectionModel = new SelectionModel(model);
			HyperStackDisplayer displayer = new HyperStackDisplayer(model, selectionModel, settings.imp);
			displayer.render();
			displayer.refresh();
		}
		// Return spot collection
		return model.getSpots();
	}
}