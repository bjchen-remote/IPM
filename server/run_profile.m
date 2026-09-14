% Standard server launcher. The numerical solver still has one entry: ipm.solve.
serverDirectory = fileparts(mfilename('fullpath'));
projectRoot = fileparts(serverDirectory);
restoredefaultpath;
addpath(projectRoot,serverDirectory);
launch_profile(profile_settings(projectRoot),projectRoot);
