################################################################################
## pfoundMultimedia.tcl                                                       ##
################################################################################
## This package is based on vmdmovie1.5. I was getting some problems		  ##
## to make movies using the original package, when calling the procedure      ##
## ::MovieMaker::makemovie {......}											  ##
## The problem is that when running in batch mode, not all of the variables   ##
## are updated, and I get a movie/picture with only one frame. Removed all tk ##
## stuff and commands for other OSes other than linux.						  ##
################################################################################
## File location: $pfound/lib/pfoundMultimedia1.0							  ##
################################################################################
## Author: Nuno Loureiro-Ferreira / nunolf (at) ci uc pt                      ##
##																			  ##
## 11 Dez 2006 file version 1.0												  ##
################################################################################
## 2Do:																		  ##
## clean up this package, it's an authentic mess							  ##
## 																			  ##
## POVRAY beta version 3.7 (presently working with 3.6) notes:				  ##
## 	The most significant change from the end-user point of view between		  ##
##	versions 3.6 and 3.7 is the addition of SMP(symmetric multiprocessing)	  ##
##	support, which in a nutshell allows the renderer to run on as many		  ##
##	CPU's (or cores) as you have installed on your computer. This will be	  ##
##	particularly useful for those users who intend purchasing a dual-core	  ##
##	CPU or who already have a two (or more) processor machine.				  ##
##																			  ##
## from : http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/4501.html    ##
## 	cutting down the file size of an animated gif from ImageMagick's		  ##
##	'convert' is easy.use 'gifsicle' instead:http://www.lcdf.org/gifsicle/	  ##
##	convert just glues the gifs together, whereas you can save a			  ##
##	significant bit of the file size by storing only the changes between	  ##
##	frames, as it is done by gifsicle. 										  ##
##																			  ##
## play with the scale command, or the following commands (default values)    ##
##	display nearclip set 0.500000											  ##
##	display farclip set 10.000000											  ##
##	display eyesep 0.060000													  ##
##	display focallength 2.000000											  ##
##	display height 6.000000													  ##
##	display distance -2.000000												  ##
##																			  ##
## try to implement a temporal bar during the trajectory movies	              ##
##																			  ##
## povray renderer sometimes gives problems (does not write the scene to .ppm)##
## 																			  ##
## think about the value to be given to the following variables: 			  ##
##  	trjstep, movieduration												  ##
##																			  ##
################################################################################

################################################################################
## This is the original (more or less) header of vmdmovie.tcl file:
##
## VMD Movie Generator version 1.4
##
## A script to generate movies in batch mode processing under machines
## running Linux OS.
## Supports a number of built-in movie types, can be easily modified
## or customized by the user to make new kinds of movies.
##
## Authors: John E. Stone, Justin Gullingsrud
##          vmd@ks.uiuc.edu
##
## $Id: vmdmovie.tcl,v 1.98 2006/08/19 20:56:11 johns Exp $
##
## Home Page
## ---------
##   http://www.ks.uiuc.edu/Research/vmd/plugins/vmdmovie
##
## Instructions:
## -------------
##   Make sure you have the NetPBM 9.23 (or later) package
##
##   Make sure you have ImageMagick in your command path 
##   (required for "animated GIF" movies)
##
##   Make sure you have "ffmpeg" in your command path
##   (optional MPEG-1 video encoder)
##
## If you have problems:
## ---------------------
##   Several of the back-end video encoders are very touchy/flaky
##   about the input files they can handle, and many of them have 
##   little or no error handling.  Try changing the size of the 
##   VMD window you're rendering from, or changing the output resolution
##   if the video encoders give you trouble.
##
################################################################################

## Tell Tcl that we're a package and any dependencies we may have
package require pfoundBasic 1.0
package require exectool 1.2
package require Orient 1.01
package provide pfoundMultimedia 1.0

namespace eval ::pfoundMultimedia:: {
  namespace export genMovie

  # global settings for work directories etc
  global tmpdir
  variable workdir     ""            ;# area for temp and output files
  variable basename    "untitled"    ;# basename for animation files 

  # parameters for the movie generator routines
  variable movietype   "rockandroll" ;# what kind of movie we're making
  	# "rockandroll" 	--> Rock and Roll (XY lemniscate)
  	# "rotation"    	--> Rotation about Y axis
  	# "trajectory"  	--> Trajectory
  	# "trajectoryrock"	--> Trajectory Rock
  	# "userdefined"		--> User Defined Procedure
  variable numframes     1           ;# number of real frames in animation
  variable anglestep     1           ;# change in angle per frame step 
  variable trjstart      0           ;# first trajectory frame in anim
  variable trjend        1           ;# last trajectory frame in anim
  variable trjstep       1           ;# step size for trajectory playback
  variable movieduration 10          ;# desired movie duration
  variable rotateangle   180         ;# how much rotation

  # parameters for renderers and video encoders
  variable framerate    24           ;# playback frame rate of animation
	# "24"			--> 24 (Film)
  	# "25"			--> 25 (PAL Video)
  	# "30"			--> 30 (NTSC Video)
  variable renderer    "snapshot"    ;# which rendering method to use
	# "snapshot"		--> Snapshot (Screen Capture)
  	# "libtachyon"		--> Internal Tachyon (Ray Tracer)
  	# "tachyon"		--> Tachyon (Ray Tracer)
  	# "povray"		--> POV-Ray (Ray Tracer)
  variable imgformat   "tga"         ;# which format generated images will use
  variable movieformat "ppmtompeg"   ;# which compress video format to use
	# "imgif"		--> Animated GIF (ImageMagick)
  	# "jpegframes"		--> JPEG frames (ImageMagick)
  	# "targaframes"		--> Targa frames (ImageMagick)
  	# "ppmtompeg"		--> MPEG-1 (ppmtompeg)
  variable bitrate     "8000000"     ;# bitrate for movie encoders

  # post-processing, pre-compression settings and options
  variable presmooth    "0"           ;# smooth images prior to compression
  variable prescale     "0"           ;# post blur enlargement/decimation
  variable scalefactor  "0.5"         ;# post blur enlargement/decimation factor
  ;# default string for  image text labels
  variable prelabel     "0"
  variable labeltext    {Produced by VMD: http://www.ks.uiuc.edu/Research/vmd/}
  variable labelsize    "7"           ;# height of characters above baseline
  variable labelrow     "auto"        ;# baseline pixel row
  variable labelcolor   "gray"        ;# color to use, black, gray, or white
  variable labelbgcolor "transparent" ;# background color or transparent

  # rendering/compression status information
  variable statusmsg   "Ready      " ;# most recent status message
  variable statusstep  "0"           ;# what step we're on in the process
  variable statussteps "0"           ;# how many steps to complete in process
  variable statusfrac  "0"           ;# fraction of total work to complete
  variable statustotal "0"           ;# total work units to complete for phase
  variable usercancel  "0"           ;# user cancellation of running anim

  # post-compression steps
  variable cleanfiles   "1"           ;# delete intermediate image frames etc

  # place to stuff exception info from catch
  variable foo                       ;# for catching exceptions
  variable viewpoints                ;# used by the rock+roll code

  # variables meant to be used by user-defined movie procs
  variable userframe     "0"         ;# traceable variable
}

################################################################################
# Main procedure															   #
# examples:							 										   #
# genMovie simID rotation imgif tachyon										   #
# genMovie simID trajectory ppmtompeg tachyon								   #
################################################################################
proc ::pfoundMultimedia::genMovie { bn mt mf rm } {
  global   tmpdir
  variable workdir 
  variable basename
  variable movietype
  variable movieformat
  variable renderer
  variable movieduration
  variable framerate
  variable rotateangle
  variable trjstep
  variable framerate
  variable presmooth
  variable prescale
  variable prelabel
  variable labeltext
  variable labelsize
  variable labelcolor

  set workdir $tmpdir	;# Temporary working dir
  set basename $bn
  set movietype $mt
  set movieformat $mf
  set renderer $rm
  set framerate 24		;# Animation frame-rate
  set presmooth 1		;# Image smoothing turned on
  set prescale 1		;# Half-size rescaling turned on (scalefactor == 0.5)
  set prelabel 1		;# Text labelling turned on
  
  # Rotation angle for rotate movies
  # Updates the following variables, if $movietype == "rotation"
  # --> numframes & anglestep
  set rotateangle 360
  ::pfoundMultimedia::angleChanged $rotateangle

  # Trajectory step size
  # Updates the following variables, if $movietype == "trajectory"
  # --> numframes & movieduration
  set trjstep 50
  ::pfoundMultimedia::trajstepChanged $trjstep

  # Movie duration. Editable only for rotate and rock 'n roll
  # Updates the following variables:
  # 	if $movietype == "rockandroll" 		--> numframes
  # 	if $movietype == "rotation" 		--> numframes & anglestep
  # 	if $movietype == "trajectory" 		--> numframes & movieduration
  # 	if $movietype == "trajectoryrock" 	--> numframes & movieduration
  set movieduration 10
  ::pfoundMultimedia::durationChanged $movieduration

  #::pfoundMultimedia::ChangeRep

  switch $movietype {
    # want to produce here a
    rotation {
      set basename ${bn}.rot
      set labelsize    "7"
      set labelcolor   "black"
      #if {$renderer == "povray"} {
      #    set presmooth 0
      #    set prescale 0
      #    set prelabel 0
      #}
      buildmovie
      exec mv $workdir/${basename}.gif ${::processedPath}/$bn
    }
    trajectory {
      set labelsize    "7"
      set labelcolor   "white"
      color Display Background black
      display update ui
      if {$movieformat == "imgif"} {
        set basename ${bn}.trj
        #if {$renderer == "povray"} {
        #  set presmooth 0
        #  set prescale 0
        #  set prelabel 0
        #}
        buildmovie
        exec mv $workdir/${basename}.gif ${::processedPath}/$bn
      }
      if {$movieformat == "ppmtompeg"} {
        set basename ${bn}.trj
        #if {$renderer == "povray"} {
        #  set presmooth 0
        #  set prescale 0
        #  set prelabel 0
        #}
        buildmovie
        exec mv $workdir/${basename}.mpg ${::processedPath}/$bn
      }
    }
  }
}

################################################################################
# Movie Rotation angle - only editable if movietype == rotation				   #
################################################################################
proc ::pfoundMultimedia::angleChanged { newval } {
  variable movietype
  variable numframes
  variable anglestep

  if { $movietype != "rotation" } { return 1 }
  if { ![string length $newval] } { return 1 }
  if { ![string is integer $newval] } { return 0 }
  set anglestep [expr double($newval) / $numframes]
  if { $anglestep == 0 } {
    set numframes 1
  } else {
    set numframes [expr round($newval / $anglestep)]
  }
  return 1
}

################################################################################
# Trajectory step - only editable if movietype == trajectory				   #
################################################################################
proc ::pfoundMultimedia::trajstepChanged { newval } {
  variable movietype
  variable numframes
  variable movieduration
  variable framerate

  # XXX trajectory_rock ignore trjstep for no apparent reason, so we only
  # recalculate numframes and movieduration for 'trajectory'
  if { $movietype != "trajectory" } { return 1 }
  if { ![string length $newval] } { return 1 }
  if { ![string is integer $newval] } { return 0 }
  if { $newval < 1 } { return 0 }
  set totframes [molinfo top get numframes]
  set numframes [expr round(double($totframes / $newval))]
  set movieduration [expr round(double($numframes / $framerate))]
  return 1
}

################################################################################
# Movie duration was edited.  Check if that's allowed, and, if so,			   #
# change any other settings that depend on that variable.					   #
# This proc also gets called when we change movie types; we detect this		   #
# situation when args has zero length.										   #
################################################################################
proc ::pfoundMultimedia::durationChanged { args } {
  variable movietype
  variable framerate
  variable numframes
  variable anglestep
  variable rotateangle
  variable movieduration
  variable trjstep

  set oldval $movieduration
  if [llength $args] { 
    set iarg [lindex $args 0] 
    if { [catch { set tmp [expr 1 * $iarg] } ] } {
      set newval $oldval
    } else {
      if { $tmp > 0 } {
        set newval $tmp
      } else {
        set newval 1
      }
    }
  } else {
    set newval $oldval
  }

  switch $movietype {
    rockandroll { 
      if [llength $args] {
        if { ![string length $newval] } { return 1 }
        if { ![string is integer $newval] } { return 0 }
      }
      set numframes [expr $newval * $framerate] 
      return 1
    }
    rotation { 
      if [llength $args] {
        if { ![string length $newval] } { return 1 }
        if { ![string is integer $newval] } { return 0 }
      }
      set numframes [expr $newval * $framerate]
      set anglestep [expr double($rotateangle) / $numframes]
      return 1
    }
    trajectory {
      if [llength $args] {
        if { [expr $oldval] != [expr $newval] } {
          puts "Sorry, can't change duration for movie of type $movietype."
          return 0
        }
      }
      if {[llength [molinfo list]] > 0} {
        set totframes [molinfo top get numframes]
        if { ![string length $trjstep] } {
          set numframes $totframes
        } else {
          set numframes [expr round(double($totframes)/$trjstep)]
        }
      } else {
        set numframes 1
      }
      set movieduration [expr {round($numframes / $framerate)}]
    }
    trajectoryrock {
      if [llength $args] {
        if { [expr $oldval] != [expr $newval] } {
          puts "Sorry, can't change duration for movie of type $movietype."
          return 0
        }
      }
      if {[llength [molinfo list]] > 0} {
        set numframes [expr ([molinfo top get numframes] * 2)]
      } else {
        set numframes 1
      }
      set movieduration [expr {round($numframes / $framerate)}]
    }
    default {
      set numframes [expr $newval * $framerate]
    }
  }

  # return success by default
  return 1
}

################################################################################
##
################################################################################
proc ::pfoundMultimedia::buildmovie {} {
  global tcl_platform
  variable imgformat 
  variable renderer
  variable statusmsg
  variable statusstep
  variable statussteps
  variable statusfrac
  variable statustotal
  variable usercancel
  variable foo
  variable movietype
  variable movieformat
  variable workdir
  variable numframes
  variable presmooth
  variable cleanfiles
  variable prescale
  variable prelabel
  variable scalefactor

  # begin process
  set usercancel "0"
  set statusmsg "Preparing  " 
  set statusstep  "1"
  set statussteps "8"
  update    ;# update the Tk window, callbacks, etc

  # check to make sure the destination filesystem is writable
  # and that output files can be written.
  if {[::pfoundMultimedia::testfilesystem] != 1} {
    puts "Temporary working directory $workdir is not usable."
    puts "Please double check file permissions and available"
    puts "space, and try again, or choose a different directory."
    return;
  }

  # set image format according to platform and renderer
  switch $renderer {
    snapshot {
      set imgformat "ppm"
    }
    libtachyon {
      set imgformat "ppm"
    }
    tachyon {
      set imgformat "ppm"
    }
    povray {
      set imgformat "ppm"
    }
    default {
      set imgformat "tga"
    }
  }

  if {$usercancel == "0"} {
    set statusmsg "Rendering  " 
    set statusstep  "2"
    update    ;# update the Tk window, callbacks, etc

    switch $movietype {
      trajectory {
        genframes_trajectory       ;# Generate frames from VMD 
      }
      trajectoryrock {
        genframes_trajectory_rock  ;# Generate frames from VMD 
      }
      rotation {
        genframes_rotation         ;# Generate frames from VMD 
      }
      rockandroll {
        genframes_rockandroll      ;# Generate frames from VMD 
      }
      userdefined {
        genframes_userdefined      ;# Generate frames from VMD 
      }
      default {
        genframes_rockandroll      ;# Generate frames from VMD 
      }
    }
  }

  if {$usercancel == "0"} {
    set statusmsg "Converting  " 
    set statusstep  "3"
    update    ;# update the Tk window, callbacks, etc
    convertframes ppm 0   ;# Convert frames to PPM format
  }

  if {$usercancel == "0"} {
    set statusmsg "Smoothing   " 
    set statusstep  "4"
    update    ;# update the Tk window, callbacks, etc

    if {$presmooth == "1"} {
      smoothframes           ;# Smooth frames prior to compression
    }
  }

  if {$usercancel == "0"} {
    set statusmsg "Rescaling   " 
    set statusstep  "5"
    update    ;# update the Tk window, callbacks, etc

    if {$prescale == "1"} {
      rescaleframes          ;# Rescale frames prior to compression
    }
  }

  if {$usercancel == "0"} {
    set statusmsg "Text Labels " 
    set statusstep  "6"
    update    ;# update the Tk window, callbacks, etc

    if {$prelabel == "1"} {
      labelframes            ;# Add text labels to frames prior to compression
    }
  }

  if {$usercancel == "0"} {
    set statusmsg "Encoding    " 
    set statusstep  "7"
    update    ;# update the Tk window, callbacks, etc

    switch $movieformat {
      ffmpeg {
        ffmpeg           ;# Convert frame sequence to MPEG-1 video file
      }
      jpegframes {
        jpegframes       ;# Convert frame sequence to JPEG format
      }
      targaframes {
        targaframes      ;# Convert frame sequence to Truevision Targa format
      }
      ppmtompeg {
        ppmtompeg        ;# Convert frame sequence to MPEG-1 video file
      }
      imgif -
      default {
        imgif            ;# Convert frame sequence to animated GIF file
      }
    }
  }

  set statusmsg "Cleaning    " 
  set statusstep  "8"
  if {$cleanfiles == "1"} {
    cleanframes          ;# delete temporary files and single frames
  }
  update    ;# update the Tk window, callbacks, etc

  if {$usercancel == "1"} {
    set statusmsg "Cancelled  "
    puts "Movie generation Cancelled."
  } else {
    set statusmsg "Completed   "
    puts "Movie generation complete."
  }
  update    ;# update the Tk window, callbacks, etc

  ## reset status area
  set statusmsg   "Ready      " ;# most recent status message
  set statusstep  "0"           ;# what step we're on in the process
  set statussteps "0"           ;# how many steps to complete in process
  set statusfrac  "0"           ;# fraction of total work to complete
  set statustotal "0"           ;# total work units to complete for phase
  update    ;# update the Tk window, callbacks, etc
}

################################################################################
# Test for file creation capability for work areas							   #
################################################################################
proc ::pfoundMultimedia::testfilesystem {} {
  variable workdir;
 
  # test access permissions on working directory
  if {[file isdirectory $workdir] != 1} {
    return 0; # failure 
  }    
  if {[file readable  $workdir] != 1} {
    return 0; # failure 
  }    
  if {[file writable  $workdir] != 1} {
    return 0; # failure 
  }    

  return 1; # success  
}

################################################################################
## Generate all of the frames												   #
################################################################################
proc ::pfoundMultimedia::genframes_rotation {} {
  variable workdir
  variable numframes
  variable anglestep
  variable renderer
  variable basename 
  variable basefilename
  variable statusfrac
  variable statustotal
  variable usercancel

  set statusfrac "0"
  set statustotal $numframes  
  update    ;# update the Tk window, callbacks, etc

  # Loop over all frames and generate images.
  set frame 0
  for {set i 0} {$i < $numframes} {incr i 1} {
    if {$usercancel == "1"} {
      set statusmsg "Cancelled  "
      return
    }

    set basefilename [format "%s/$basename.%04d" $workdir $frame]
    display update ui
    renderframe $basefilename

    rotate y by $anglestep

    incr frame
    set statusfrac $frame
  }
}

################################################################################
# Generate all of the frames												   #
################################################################################
proc ::pfoundMultimedia::genframes_trajectory {} {
  variable workdir
  variable numframes
  variable trjstep
  variable renderer
  variable basename 
  variable basefilename
  variable statusfrac
  variable statustotal
  variable usercancel
  variable frame
  variable trjframe

  puts "numframes: $numframes"
  puts "Generating image frames using renderer: $renderer"
  set statusfrac "0"
  set statustotal $numframes  
  update    ;# update the Tk window, callbacks, etc

  # Loop over all frames and generate images.
  set frame 0
  set trjframe 0
  for {set i 0} {$i < $numframes} {incr i 1} {
    if {$usercancel == "1"} {
      set statusmsg "Cancelled  "
      return
    }

    set basefilename [format "%s/$basename.%04d" $workdir $frame]

    animate goto $trjframe

    display update ui
    # updates the display of the SS	
    mol ssrecalc top
    renderframe $basefilename

    incr frame
    incr trjframe $trjstep
    set statusfrac $frame
  }
}

################################################################################
# Render one frame using existing settings									   #
################################################################################
proc ::pfoundMultimedia::renderframe { basefilename } {
  global env
  variable renderer
  variable rendercmd
  variable scenefilename 
  variable imgfilename
  variable imgformat

  # set platform-specific executable suffix
  set archexe ""

  set imgfilename $basefilename.$imgformat

  switch $renderer {
    snapshot {
      render snapshot $imgfilename
    } 
    libtachyon {
      render TachyonInternal $imgfilename
    } 
    tachyon {
      set scenefilename $basefilename.dat
      set tachyonexe [format "tachyon%s" $archexe];
      set tachyoncmd \
        [::ExecTool::find -interactive -description "Tachyon Ray Tracer" \
    -path [file join $env(VMDDIR) "tachyon_[vmdinfo arch]$archexe"] $tachyonexe]
      if {$tachyoncmd == {}} {
        error "Cannot find Tachyon, aborting"
      }
      set rendercmd [ format "\"%s\"" $tachyoncmd]

      switch $imgformat {
        ppm {
          set rendercmd [concat $rendercmd \
         "-mediumshade $scenefilename -format PPM -aasamples 4 -o $imgfilename"]
        }

        rgb {
          set rendercmd [concat $rendercmd \
         "-mediumshade $scenefilename -format RGB -aasamples 4 -o $imgfilename"]
        }

        tga {
          set rendercmd [concat $rendercmd \
       "-mediumshade $scenefilename -format Targa -aasamples 4 -o $imgfilename"]
        }

        default { 
          puts "Image format unsupported, aborting"
        } 
      }
      render Tachyon $scenefilename $rendercmd
    } 
    povray {
      set scenefilename $basefilename.pov
      set povrayexe [format "povray%s" $archexe];
      set povraycmd \
    [::ExecTool::find -interactive -description "POV-Ray Ray Tracer" $povrayexe]
      if {$povraycmd == {}} {
        error "Cannot find POV-Ray, aborting"
      }
      set rendercmd [ format "\"%s\"" $povraycmd]

      # get image resolution so we can pass it into POV command line
      set imagesize [display get size]
      set xsize [lindex $imagesize 0]
      set ysize [lindex $imagesize 1]

      set rendercmd [concat $rendercmd \
             "-I$scenefilename -O$imgfilename +X +A +FP +W$xsize +H$ysize"]
            ;# "-I$scenefilename -O$imgfilename -D +X +A +FP +W$xsize +H$ysize"]

      render POV3 $scenefilename $rendercmd
    }
    default {
      puts "Unsupported renderer"
    }
  }
  # deletes ascii rendering scene file (they are very big !!)
  if {$renderer == "tachyon"} { file delete -force $scenefilename }
}

################################################################################
# Convert to MPEG using "ppmtompeg" 										   #
################################################################################
proc ::pfoundMultimedia::ppmtompeg {} {
  variable numframes
  variable framerate
  variable basename
  variable workdir
  variable mybasefilename
  variable statusfrac
  variable statustotal
  variable foo
  variable parfile
  variable frame
  variable framename
  variable lastframe

  set statusfrac  "0"
  set statustotal "1"
  update    ;# update the Tk window, callbacks, etc

  puts "Converting frames to MPEG-1 video format"
  
  set mybasefilename [format "%s/$basename" $workdir] 

  # generate MPEG-1 encoder parameter file
  set parfile [open "$mybasefilename.par" w]
  puts $parfile "PATTERN    IBBPBBPBBPBBPBB"
  puts $parfile "FORCE_ENCODE_LAST_FRAME"          ;# force anim loopable
  puts $parfile "OUTPUT     $mybasefilename.mpg"
  puts $parfile "INPUT_DIR  $workdir"
  puts $parfile "INPUT"

  set lastframe [format "%04d" [expr $numframes - 1]]
  puts $parfile "$basename.*.ppm \[0000-$lastframe\]"

  puts $parfile "END_INPUT"
  puts $parfile "BASE_FILE_FORMAT PPM"
  puts $parfile "INPUT_CONVERT *"
  puts $parfile "GOP_SIZE 15"
  puts $parfile "SLICES_PER_FRAME 1"
  puts $parfile "PIXEL HALF"
  puts $parfile "RANGE 32"
  puts $parfile "PSEARCH_ALG LOGARITHMIC"
  puts $parfile "BSEARCH_ALG CROSS2"
  puts $parfile "IQSCALE 8"
  puts $parfile "PQSCALE 10"
  puts $parfile "BQSCALE 25"
  puts $parfile "REFERENCE_FRAME DECODED"
  close $parfile

  puts "ppmtompeg $mybasefilename.par"
  file delete -force $mybasefilename.mpg
  set foo [catch { ::ExecTool::exec ppmtompeg $mybasefilename.par >@ stdout }]
}

################################################################################
# Convert to animated GIF using ImageMagick									   #
################################################################################
proc ::pfoundMultimedia::imgif {} {
  variable numframes
  variable framerate
  variable basename
  variable workdir
  variable mybasefilename
  variable statusfrac
  variable statustotal
  variable foo
  variable delay

  set statusfrac  "0"
  set statustotal "1"
  update    ;# update the Tk window, callbacks, etc

  puts "Converting frames to animated GIF format"
  set mybasefilename [format "%s/$basename" $workdir] 
  set delay [format "%5.2f" [expr 100.0 / $framerate]] 
 
  # flag loop = 0, means that the animated gif will always be looping
  puts "convert -delay $delay -loop 0 $mybasefilename.*.ppm $mybasefilename.gif"
  file delete -force $mybasefilename.gif
  # I'm having problems on behemoth (but not on yose)
  # On behemoth, ::ExecTool::exec convert gives an error saying that convert is
  # not found. Check it when i have time.
  ::ExecTool::exec convert -delay $delay -loop 0 \
					$mybasefilename.*.ppm $mybasefilename.gif >@ stdout
  #exec /usr/local/bin/convert -delay $delay -loop 0 \
#		$mybasefilename.*.ppm $mybasefilename.gif >@ stdout
}

################################################################################
# Convert frames to the format required by the selected video encoder	       #
################################################################################
proc ::pfoundMultimedia::convertframes { newformat finalformat } {
  variable imgformat
  variable numframes
  variable oldfilename
  variable newfilename
  variable basename
  variable workdir
  variable statusfrac
  variable statustotal
  variable usercancel
  variable foo

  set statusfrac  "0"
  set statustotal $numframes
  update    ;# update the Tk window, callbacks, etc

  if { $finalformat } {
    set prefix "final."
    if { $imgformat == $newformat } {
      puts "No frame format conversion necessary, copying to final destination."
      return
    }
  } else {
    set prefix ""
    if { $imgformat == $newformat } {
      puts "No frame format conversion necessary, continuing."
      return
    }
  }


  puts "Converting frames from $imgformat format to $newformat..."

  # Loop over all frames and convert images.
  switch $newformat {
    ppm {
      switch $imgformat {
        tga {
          for {set i 0} {$i < $numframes} {incr i 1} {
            if {$usercancel == "1"} {
              set statusmsg "Cancelled  "
              return
            }
            set oldfilename [format "%s/$basename.%04d.$imgformat" $workdir $i]
            set newfilename \
				[format "%s/$prefix$basename.%04d.$newformat" $workdir $i]
            file delete -force $newfilename
            ::ExecTool::exec tgatoppm $oldfilename > $newfilename
            set statusfrac $i
            update    ;# update the Tk window, callbacks, etc
          }
        }
        rgb {
          for {set i 0} {$i < $numframes} {incr i 1} {
            if {$usercancel == "1"} {
              set statusmsg "Cancelled  "
              return
            }
            set oldfilename [format "%s/$basename.%04d.$imgformat" $workdir $i]
            set newfilename \
				[format "%s/$prefix$basename.%04d.$newformat" $workdir $i]
            file delete -force $newfilename
            set foo \
				[catch {::ExecTool::exec sgitopnm $oldfilename > $newfilename}]
            set statusfrac $i
            update    ;# update the Tk window, callbacks, etc
          }
        }
        ppm {
          for {set i 0} {$i < $numframes} {incr i 1} {
            if {$usercancel == "1"} {
              set statusmsg "Cancelled  "
              return
            }
            set oldfilename [format "%s/$basename.%04d.$imgformat" $workdir $i]
            set newfilename \
				[format "%s/$prefix$basename.%04d.$newformat" $workdir $i]
            file delete -force $newfilename
            file copy -force $oldfilename $newfilename
            set statusfrac $i
            update    ;# update the Tk window, callbacks, etc
          }
        }
        default { 
          puts "Conversion from $imgformat to $newformat unsupported, aborting"
        } 
      }
    }


    jpg {
      switch $imgformat {
        ppm -
        tga -
        rgb {
          for {set i 0} {$i < $numframes} {incr i 1} {
            if {$usercancel == "1"} {
              set statusmsg "Cancelled  "
              return
            }
            set oldfilename [format "%s/$basename.%04d.$imgformat" $workdir $i]
            set newfilename \
				[format "%s/$prefix$basename.%04d.$newformat" $workdir $i]
            file delete -force $newfilename
            ::ExecTool::exec \
				convert -quality 100% $oldfilename JPEG:$newfilename >@ stdout
            set statusfrac $i
            update    ;# update the Tk window, callbacks, etc
          }
        }
        jpg {
          for {set i 0} {$i < $numframes} {incr i 1} {
            if {$usercancel == "1"} {
              set statusmsg "Cancelled  "
              return
            }
            set oldfilename [format "%s/$basename.%04d.$imgformat" $workdir $i]
            set newfilename \
				[format "%s/$prefix$basename.%04d.$newformat" $workdir $i]
            file delete -force $newfilename
            file copy -force $oldfilename $newfilename
            set statusfrac $i
            update    ;# update the Tk window, callbacks, etc
          }
        }
        default { 
          puts "Conversion from $imgformat to $newformat unsupported, aborting"
        } 
      }
    }


    tga {
      switch $imgformat {
        ppm -
        jpg -
        rgb {
          for {set i 0} {$i < $numframes} {incr i 1} {
            if {$usercancel == "1"} {
              set statusmsg "Cancelled  "
              return
            }
            set oldfilename [format "%s/$basename.%04d.$imgformat" $workdir $i]
            set newfilename \
				[format "%s/$prefix$basename.%04d.$newformat" $workdir $i]
            file delete -force $newfilename
            ::ExecTool::exec convert $oldfilename $newfilename >@ stdout
            set statusfrac $i
            update    ;# update the Tk window, callbacks, etc
          }
        }
        tga {
          for {set i 0} {$i < $numframes} {incr i 1} {
            if {$usercancel == "1"} {
              set statusmsg "Cancelled  "
              return
            }
            set oldfilename [format "%s/$basename.%04d.$imgformat" $workdir $i]
            set newfilename \
				[format "%s/$prefix$basename.%04d.$newformat" $workdir $i]
            file delete -force $newfilename
            file copy -force $oldfilename $newfilename
            set statusfrac $i
            update    ;# update the Tk window, callbacks, etc
          }
        }
        default { 
          puts "Conversion from $imgformat to $newformat unsupported, aborting"
        } 
      }
    }

    default {
      puts "Conversion unsupported, aborting"
    }
  }
}

################################################################################
# Convert frames to the format required by the selected video encoder	       #
################################################################################
proc ::pfoundMultimedia::smoothframes { } {
  variable numframes
  variable oldfilename
  variable newfilename
  variable basename
  variable workdir
  variable statusfrac
  variable statustotal
  variable usercancel

  set statusfrac  "0"
  set statustotal $numframes
  update    ;# update the Tk window, callbacks, etc

  puts "Smoothing image frames."

  # Loop over all frames and process images.
  for {set i 0} {$i < $numframes} {incr i 1} {
    if {$usercancel == "1"} {
      set statusmsg "Cancelled  "
      return
    }
    set oldfilename [format "%s/$basename.%04d.ppm" $workdir $i]
    set newfilename [format "%s/orig.$basename.%04d.ppm" $workdir $i]
    file delete -force $newfilename
    file rename -force -- $oldfilename $newfilename
    ::ExecTool::exec pnmsmooth $newfilename > $oldfilename
    file delete -force $newfilename
    set statusfrac $i
    update    ;# update the Tk window, callbacks, etc
  }
}

################################################################################
# Convert frames to the format required by the selected video encoder	       #
################################################################################
proc ::pfoundMultimedia::rescaleframes { } {
  variable numframes
  variable oldfilename
  variable newfilename
  variable basename
  variable workdir
  variable statusfrac
  variable statustotal
  variable usercancel
  variable scalefactor

  set statusfrac  "0"
  set statustotal $numframes
  update    ;# update the Tk window, callbacks, etc

  puts "Rescaling image frames."

  # Loop over all frames and process images.
  for {set i 0} {$i < $numframes} {incr i 1} {
    if {$usercancel == "1"} {
      set statusmsg "Cancelled  "
      return
    }
    set oldfilename [format "%s/$basename.%04d.ppm" $workdir $i]
    set newfilename [format "%s/orig.$basename.%04d.ppm" $workdir $i]
    file delete -force $newfilename
    file rename -force -- $oldfilename $newfilename
    ::ExecTool::exec pnmscale $scalefactor  $newfilename > $oldfilename
    file delete -force $newfilename
    set statusfrac $i
    update    ;# update the Tk window, callbacks, etc
  }
}

################################################################################
# Tag frames with copyright or other text, prior to video encoding			   #
################################################################################
proc ::pfoundMultimedia::labelframes { } {
  variable numframes
  variable oldfilename
  variable newfilename
  variable basename
  variable workdir
  variable statusfrac
  variable statustotal
  variable usercancel
  variable labeltext
  variable labelcolor
  variable labelbgcolor
  variable labelsize
  variable labelrow
  variable localrow

  set statusfrac  "0"
  set statustotal $numframes
  # Text label to put on movies
  set labeltext "Pfound - Generated by VMD\n[::pfoundBasic::timeFormat]"
  update    ;# update the Tk window, callbacks, etc

  puts "Tagging image frames with user-specified text."

  if {$labelrow == "auto"} {
    set localrow [expr $labelsize + 3];
  } else {
    set localrow $labelrow
  }

  # Loop over all frames and process images.
  for {set i 0} {$i < $numframes} {incr i 1} {
    if {$usercancel == "1"} {
      set statusmsg "Cancelled  "
      return
    }
    set oldfilename [format "%s/$basename.%04d.ppm" $workdir $i]
    set newfilename [format "%s/orig.$basename.%04d.ppm" $workdir $i]
    file delete -force $newfilename
    file rename -force -- $oldfilename $newfilename

    ##
    ## image compositing method for image-based logos etc...
    ## ::ExecTool::exec  pnmcomp -invert -alpha stamp-30.pgm black-30.ppm \
    ##		$newfilename > $oldfilename
    ##

    ::ExecTool::exec ppmlabel -size $labelsize -y $localrow -color $labelcolor \
          -background $labelbgcolor -text $labeltext $newfilename > $oldfilename
    file delete -force $newfilename
    set statusfrac $i
    update    ;# update the Tk window, callbacks, etc
  }
}

################################################################################
# Composite frames with copyright or other images, prior to video encoding     #
################################################################################
proc ::pfoundMultimedia::compositeframes { } {
  variable numframes
  variable oldfilename
  variable newfilename
  variable basename
  variable workdir
  variable statusfrac
  variable statustotal
  variable usercancel

  set statusfrac  "0"
  set statustotal $numframes
  update    ;# update the Tk window, callbacks, etc

  puts "Compositing image frames with user-specified image."

  # Loop over all frames and process images.
  for {set i 0} {$i < $numframes} {incr i 1} {
    if {$usercancel == "1"} {
      set statusmsg "Cancelled  "
      return
    }
    set oldfilename [format "%s/$basename.%04d.ppm" $workdir $i]
    set newfilename [format "%s/orig.$basename.%04d.ppm" $workdir $i]
    file delete -force $newfilename
    file rename -force -- $oldfilename $newfilename
    
    # image compositing method for image-based logos etc...
    ::ExecTool::exec  pnmcomp -invert -alpha stamp-30.pgm black-30.ppm \
		      $newfilename > $oldfilename
    file delete -force $newfilename
    set statusfrac $i
    update    ;# update the Tk window, callbacks, etc
  }
}

################################################################################
# Clean up all of the temporary frames once the whole						   #
# animation is done.														   #
################################################################################
proc ::pfoundMultimedia::cleanframes {} {
  variable numframes
  variable basename
  variable workdir 
  variable filespec
  variable fileset
  variable statusfrac
  variable statustotal

  set statusfrac  "0"
  set statustotal $numframes
  update    ;# update the Tk window, callbacks, etc

  puts "Cleaning generated data, frames, and encoder parameter files"

  file delete -force $workdir/$basename.par
  for {set i 0} {$i < $numframes} {incr i 1} {
    set filespec [format "$workdir/$basename.%04d" $i]
    set files $filespec 
    file delete -force $files.ppm $files.tga $files.jpg \
               $files.rgb $files.bmp $files.dat $files.pov
    set statusfrac  $i
    update ;# update the Tk window, callbacks, etc
  }
}

################################################################################
# Change default rep for already loaded molecules							   #
# code core taken from VMD-L:												   #
# http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/8071.html			   #
################################################################################

proc ::pfoundMultimedia::ChangeRep { } {
    after idle {
	mol color Structure
        # Got to take care of this
		# I must do a reset before making a new movie/picture
#		scale by 1.5
        foreach mid [molinfo list] {
            mol delrep 0 $mid
            mol addrep $mid
        }
    }
}

################################################################################
# Writes a PDB from the 1st frame of the trajectory							   #
################################################################################
proc ::pfoundMultimedia::writePdb {output} {
  animate goto 0
  display update ui
  set sel [atomselect top all]
  $sel frame 0
  $sel update
  $sel writepdb ./${output}.pdb
}

# em fases de teste; implementação do pacote orient
proc ::pfoundMultimedia::genSnapshot {} {
	#::pfoundMultimedia::ChangeRep
	render Tachyon plot.dat [render options Tachyon]
	set sel [atomselect top all]
	# show/calc the principal axes
	set I [draw principalaxes $sel]
	# rotate axis 2 to match X
	set A [::Orient::orient $sel [lindex $I 2] {1 0 0}] 
	$sel move $A
	# recalc principal axes to check
	set I [draw principalaxes $sel]
	# rotate axis 0 to match Z
	set A [::Orient::orient $sel [lindex $I 0] {0 0 1}] 
	$sel move $A
	set I [draw principalaxes $sel]

# check first if tachyon is found on the PATH
render Tachyon plot1.dat [render options Tachyon]
}
