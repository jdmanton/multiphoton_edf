num_bands = getNumber("Number of pupil zones", 5);
num_slices = getNumber("Number of PSF slices", 301);

run("Stack to Hyperstack...", "order=xyczt(default) channels=1 slices=" + num_slices + " frames=" + num_bands + " display=Color");
run("Reslice [/]...", "output=1.000 start=Top avoid");
run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
run("Z Project...", "projection=[Average Intensity] all");
setSlice(floor((num_slices + 1) / 2));
resetMinAndMax();
