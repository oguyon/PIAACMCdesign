This directory contains scripts to run PIAACMC design tool

- install inotifywait (required for script waitforfile)
- copy file waitforfile to executable path (for example /usr/local/bin/)
- copy scripts run, runclean, runopt, runPIAACMC and sim1024 to working directory (or use sym links)
- edit variable "execname" in scripts runPIAACMC and runopt to point to executable (not required if you sym link ./bin/PIAACMCdesign to /usr/local/bin/PIAACMCdesign)
- edit script run (main script) and execute it

