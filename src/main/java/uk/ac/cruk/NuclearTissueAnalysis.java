
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

	@Parameter(label = "Define width of neculus edge (microns): ", persist = false, min = "0.01", max = "20")
	private double nucleusEdgeWidth = 1.5;

	@Override
	public void run() {
		// get current ResultsTable
		ResultsTable rt = new ResultsTable();
		String[] rawDataDirList = rawDataDir.list();
		// Use Trackmate to count the number of satellites
		Settings tmSettings = new Settings();
		tmSettings.detectorFactory = new LogDetectorFactory();
		tmSettings.addSpotFilter(new FeatureFilter("QUALITY", minSpotQuality, true));
		Map<String, Object> map = tmSettings.detectorFactory.getDefaultSettings();
		map.put("RADIUS", spotRadius);
		map.put("DO_MEDIAN_FILTERING", false);
		map.put("DO_SUBPIXEL_LOCALIZATION", true);
		map.put("TARGET_CHANNEL", 1);
		map.put("THRESHOLD", 0.0d);
		tmSettings.detectorSettings = map;

		for (int i = 0; i < rawDataDirList.length; i++) {
			if (rawDataDirList[i].endsWith(".tif") && !rawDataDirList[i].contains("Object Predictions")) {

				System.out.println("file " + i + " of " + rawDataDirList.length);
				ImagePlus rawSeriesImp = IJ.openImage(rawDataDir.getPath() + "\\" + rawDataDirList[i]);
				ImagePlus objectMapImp = IJ.openImage(objectMapsDir.getPath() + "\\"
						+ rawDataDirList[i].substring(0, rawDataDirList[i].length() - 4) + "_Object Predictions.tif");

				tmSettings.setFrom(rawSeriesImp);
				SpotCollection dapiSpots = countSpotsTrackmate(tmSettings, false);

				IJ.setThreshold(objectMapImp, 1.0d, 1.0d);
				IJ.run("Set Measurements...", "  redirect=None decimal=3");
				IJ.run(objectMapImp, "Analyze Particles...", "  show=[Overlay Outlines]");
				Overlay overlay = objectMapImp.getOverlay();
				if (overlay != null) {
					double boundX = rawSeriesImp.getWidth() - rawSeriesImp.getWidth() * 0.1;
					double boundY = rawSeriesImp.getHeight() - rawSeriesImp.getHeight() * 0.1;

					boolean[] roiCheck = new boolean[overlay.size()];
					int[] roiSpotCountsAll = new int[overlay.size()];
					int[] roiSpotCountsCenter = new int[overlay.size()];
					// roiCheck

					for (int r = 0; r < overlay.size(); r++) {
						Roi roi = overlay.get(r);
						Rectangle boundingBox = roi.getBounds();
						Roi roiShrink = RoiEnlarger.enlarge(roi,
								-1 * nucleusEdgeWidth / rawSeriesImp.getCalibration().pixelWidth);

						if (boundingBox.getMinX() != 0 && boundingBox.getMinY() != 0 && boundingBox.getMinX() <= boundX
								&& boundingBox.getMinY() <= boundY) {
							roiCheck[r] = true;

							Iterator<Spot> spotIterator = dapiSpots.iterator(true);
							while (spotIterator.hasNext()) {
								Spot spot = spotIterator.next();
								int positionXPix = (int) Math.floor(
										spot.getFeature("POSITION_X") / rawSeriesImp.getCalibration().pixelWidth);
								int positionYPix = (int) Math.floor(
										spot.getFeature("POSITION_Y") / rawSeriesImp.getCalibration().pixelWidth);
								if (roi.contains(positionXPix, positionYPix)) {
									roiSpotCountsAll[r]++;
								}
								if (roiShrink.contains(positionXPix, positionYPix)) {
									roiSpotCountsCenter[r]++;
								}
							}

						} else {
							roiCheck[r] = false;
						}

					}

					double[] roiAreas = new double[overlay.size()];
					double[][] roiMeans = new double[rawSeriesImp.getNChannels()][overlay.size()];
					double[][] roiStds = new double[rawSeriesImp.getNChannels()][overlay.size()];
					for (int c = 0; c < rawSeriesImp.getNChannels(); c++) {
						rawSeriesImp.setC(c + 1);
						for (int r = 0; r < overlay.size(); r++) {

							Roi roi = overlay.get(r);
							Rectangle boundingBox = roi.getBounds();

							if (roiCheck[r]) {
								rawSeriesImp.setRoi(roi);

								ImageStatistics statistics = rawSeriesImp.getStatistics();
								roiAreas[r] = statistics.area;
								roiMeans[c][r] = statistics.mean;
								roiStds[c][r] = statistics.stdDev;

							}

						}

					}

					for (int r = 0; r < overlay.size(); r++) {
						if (roiCheck[r]) {
							rt.incrementCounter();
							rt.addValue("Name", rawDataDirList[i]);
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