File.makeDirectory("/home/cylinderman/Experiments/tiff/");
run("Image Sequence... ", "format=TIFF use save=/home/cylinderman/Experiments/tiff/orig0000.tif");
run("Split Channels");
imageCalculator("Subtract create stack", "orig (blue)","orig (green)");
selectWindow("Result of orig (blue)");
setAutoThreshold("Default dark");
run("Threshold...");

//waitForUser("adjust threshold", "test the threshold parameters in the window that next appears. note the desired Min and Max");
waitForUser("Click OK when you are done");
run("Threshold...");
//setThreshold(98, 255);

setOption("BlackBackground", false);
run("Convert to Mask", "method=Default background=Dark");
run("Analyze Particles...", "size=500-5000 circularity=0.30-1.00 show=Ellipses display exclude clear include stack");
saveAs("Results", "/home/cylinderman/Experiments/Results_blue.csv");
imageCalculator("Subtract create stack", "orig (red)","orig (green)");
selectWindow("Result of orig (red)");
setAutoThreshold("Default dark");
run("Threshold...");

//waitForUser("adjust threshold", "test the threshold parameters in the window that next appears. note the desired Min and Max");
waitForUser("Click OK when you are done");
run("Threshold...");
//setThreshold(49, 255);

setOption("BlackBackground", false);
run("Convert to Mask", "method=Default background=Dark");
run("Analyze Particles...", "size=500-5000 circularity=0.30-1.00 show=Ellipses display exclude clear include stack");
saveAs("Results", "/home/cylinderman/Experiments/Results_red.csv");
selectWindow("orig (green)");
setAutoThreshold("Default dark");
run("Threshold...");

//waitForUser("adjust threshold", "test the threshold parameters in the window that next appears. note the desired Min and Max");
waitForUser("Click OK when you are done");
run("Threshold...");
//setThreshold(76, 255);

setOption("BlackBackground", false);
print("1");
run("Convert to Mask", "method=Default background=Dark");
print("2");
run("Watershed", "stack");
print("3");
run("Erode", "stack");
print("4");
run("Erode", "stack");
run("Analyze Particles...", "size=10000-Inf circularity=0.30-1.00 show=Ellipses display exclude clear include stack");
print("5");
saveAs("Results", "/home/cylinderman/Experiments/Results_green.csv");
close("Results");
close("orig (green)");
close("orig (blue)");
close("orig (red)");
close();
