---
title: 01 - UA HPC Tutorial
type: docs
prev: docs/folder/
next: docs/folder/2-prevariant-filtering
---

## Logging into the UA HPC
### Using your computer’s command line prompt
*Find your terminal*
* On a Mac, you can open the Launchpad and search for “Terminal”; alternatively, you can press Command+Space Bar, and then search for “Terminal”
* On a PC, you will need to install a separate program to access a shell prompt. You can find more information on the options on the [UA HPC documentation website](https://hpcdocs.hpc.arizona.edu/quick_start/logging_in/#system-access)
  
In the terminal, type `ssh netid@hpc.arizona.edu` and press enter, then enter your password and confirm your dual authentication when prompted. This lands you on the gatekeeper/bastion host.
Now, type `shell` and press enter. This moves you to one of the login nodes (junonia/wentletrap).

### Using the UA’s Open OnDemand interface
Go to [ood.hpc.arizona.edu](ood.hpc.arizona.edu), then select Clusters -> Shell access from the top menu bar. This will open a terminal window within your browser. When it opens, you’ll see that you’re automatically logged into a login node (junonia/wentletrap).

## Using the UA HPC
Now that you’re logged into the UA HPC and on a login node (either `junonia` or `wentletrap`), we can practice some basic Linux commands.
* `pwd` to print the directory you currently are in (you should land in your home directory when you first login)
* `ls` to see what files are in your current directory
* `hostname` (to find out which node you’re currently on)
* Moving around: `cd` (to move to a different directory), e.g. `cd /xdisk/jrick/consbioinf` to move into our shared class directory
* If we wanted make a new directory, we can use `mkdir new_directory`
* To check permissions, `ls -l `
* To look preview the contents of a file, we can use `less filename`, where `filename` is the file we want to look at
* To copy files, we use `cp`
* To request an interactive session, `interactive -a consbioinf -t 2:00:00`, where the `-t` flag indicates how long you want your interactive session to be good for (e.g., `2:00:00` indicates 2 hours, 0 minutes, 0 seconds)

## Modules!
* `module list`
* `module avail`
* `module spider`

## Moving files to/from the HPC 
There are several options: `ssh` or Globus or file transfer in Open OnDemand

