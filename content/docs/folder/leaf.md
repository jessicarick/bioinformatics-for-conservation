---
title: UA HPC Tutorial
type: docs
prev: docs/folder/
---

# Logging into the UA HPC
## Using your computer’s command line prompt
*Find your terminal*
* On a Mac, you can open the Launchpad and search for “Terminal”; alternatively, you can press Command+Space Bar, and then search for “Terminal”
* On a PC, you will need to install a separate program to access a shell prompt. You can find more information on the options on the [UA HPC documentation website](https://hpcdocs.hpc.arizona.edu/quick_start/logging_in/#system-access)
  
In the terminal, type `ssh netid@hpc.arizona.edu` and press enter, then enter your password and confirm your dual authentication when prompted. This lands you on the gatekeeper/bastion host.
Now, type `shell` and press enter. This moves you to one of the login nodes (junonia/wentletrap).

## Using the UA’s Open OnDemand interface
Go to [ood.hpc.arizona.edu](ood.hpc.arizona.edu), then select Clusters -> Shell access. This will open a terminal window within your browser. When it opens, you’ll see that you’re automatically logged into a login node (junonia/wentletrap).

# Using the UA HPC
Now that you’re logged into the UA HPC and on a login node (either `junonia` or `wentletrap`), we can practice some basic Linux commands.
* `pwd` to print the directory you currently are in (you should land in your home directory when you first login)
* `ls` to see what files are in your current directory 
 

—------
hostname (to find out which node you’re currently on)
Now, some basic Linux commands
pwd (you should land in your home directory)
ls (there likely won’t be much in your home directory for now)
Moving around: cd (to move to a different directory) → cd /xdisk/jrick/consbioinf → how can we verify that we moved to the directory that we wanted to move to? what is there in this directory? → at this point, draw out a map of file structures
If we wanted make a new directory → mkdir new_directory 
To check permissions, ls -l 
Discuss permissions + how to change them (if you are the owner of the file/directory, otherwise will get error)
chmod chgrp
There is a file already in this directory → to see what is inside of it, we can use less
Copying files → cp
Creating files → there are several different text editors available to use in Linux, but all of them look different from what you’re used to for editing text files on your own computer. I personally use vim, so that’s what I’ll be demonstrating to you today. However, if you are already comfortable with something else, you’re welcome to use that. (demonstrate creating a new text file in our shared directory, then save & close, then demonstrate editing an already-existing file)
Also, tab completion
Modules!
module list
module avail
module spider
Moving files to/from the HPC → ssh or Globus or file transfer in Open OnDemand

