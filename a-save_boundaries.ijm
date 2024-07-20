File.makeDirectory("/home/cylinderman/Experiments/tiff/");
run("Image Sequence... ", "format=TIFF use save=/home/cylinderman/Experiments/tiff/orig0000.tif");
run("Split Channels");
imageCalculator("Subtract create stack", "orig (blue)","orig (green)");
selectWindow("Result of orig (blue)");
setAutoThreshold("Default dark");
run("Threshold...");

waitForUser("adjust threshold", "test the threshold parameters in the window that next appears. note the desired Min and Max");
run("Threshold...");
waitForUser("Click OK when you are done");
//setThreshold(98, 255);

setOption("BlackBackground", false);
run("Convert to Mask", "method=Default background=Dark");
run("Analyze Particles...", "size=500-5000 circularity=0.30-1.00 show=Ellipses display exclude clear include stack");
saveAs("Results", "/home/cylinderman/Experiments/Results_blue.csv");
imageCalculator("Subtract create stack", "orig (red)","orig (green)");
selectWindow("Result of orig (red)");
setAutoThreshold("Default dark");
run("Threshold...");

waitForUser("adjust threshold", "test the threshold parameters in the window that next appears. note the desired Min and Max");
run("Threshold...");
waitForUser("Click OK when you are done");
//setThreshold(49, 255);

setOption("BlackBackground", false);
run("Convert to Mask", "method=Default background=Dark");
run("Analyze Particles...", "size=500-5000 circularity=0.30-1.00 show=Ellipses display exclude clear include stack");
saveAs("Results", "/home/cylinderman/Experiments/Results_red.csv");
close();
