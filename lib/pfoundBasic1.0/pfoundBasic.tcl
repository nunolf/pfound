################################################################################
## pfoundBasic.tcl                                                            ##
################################################################################
## Survey:                                                                    ##
## This is the main package													  ##
##																			  ##
## This package has all the procedures to calculate							  ##
## several time relevant properties on the trajectories						  ##
## uploaded to pfound data repository.										  ##
##																			  ##
## The core of the following procedures were written by jrui:				  ##
## calcRgyr calcRmsd calcHb CalcHbNative calcSs calcNatContacts distMap	      ##
## calcSasa																	  ##
################################################################################
## 2do:																		  ##
##     Check variables/ start them at namespace eval						  ##
##     Check procs to be exported (gives a clue of what are the				  ##
##     the procs available to be accessed from outside the namespace)	      ##
##     Implement SASA; generate movies; ...									  ##
##     Check robustness of procs that calculate the properties,				  ##
##     namely, are they working if the protein has more than one chain?	      ##
##     Check file extension types											  ##
##     Have to take care of AMBER trajs (with or without box information)     ##
################################################################################
## File location: $pfound/lib/pfoundBasic1.0								  ##
################################################################################
## Author: Nuno Loureiro-Ferreira / nunolf (at) ci uc pt                      ##
##																			  ##
## Acknowledgments: J. Rui Rodrigues / jrui (at) ci uc pt					  ##
##                  John Stone / vmd-l (at) ks uiuc edu						  ##
##																			  ##
## 11 Dez 2006 file version 1.0												  ##
################################################################################


# Tell Tcl that we're a package and any dependencies we may use
package require La
package provide pfoundBasic 1.0

namespace eval ::pfoundBasic:: {
    
    global tmpdir       ;# sets global temp dir
    global dataBasename ;# sets the simid being processed

    variable trajFrames
    variable nbonds
    variable topType
    variable trajType

    # Resnames are defined by a string of up to 4 characters.
    # These resnames differ depending on the f.f. used to generate
    # the trajectory. The codes array, assigns to each residue a unique
    # 3 letter code. This is necessary in order to get the primary
    # sequence of the protein, with the a.a. common names.
    # Nomenclature taken into consideration, came from
    # AMBER, GROMACS/GROMOS, OPLSA and CHARMM forcefields.
    # If for a given resname X, ![info exists $codes(X] T, then $codes(X) == UNK
    variable codes
    array set codes {
        ALA ALA
        ARG ARG  ARGN ARG  ARN ARG
        ASN ASN  ASN1 ASN
        ASP ASP  ASPH ASP  ASH ASP
        CYS CYS  CYS1 CYS  CYS2 CYS  CYSH CYS  CYN CYS  CYX CYS
        GLN GLN
        GLU GLU  GLUH GLU  GLH GLU
        GLY GLY
        HIS HIS  HISA HIS  HISB HIS  HISH HIS  HSD HIS  HSE HIS  HSP HIS  \
        HID HIS  HIE HIS  HIP HIS
        ILE ILE
        LEU LEU
        LYS LYS  LYSH LYS  LYN LYS  LYP LYS
        MET MET
        PHE PHE
        PRO PRO
        SER SER
        THR THR
        TRP TRP
        TYR TYR
        VAL VAL
    }

    # topology/coordinate & trajectory extensions recognized by VMD
    variable valExtensions
    
        # add these
    	# MMTK
	#nc (structure + coordinates)
	#LAMMPS
	#dump
	#CPMD
    array set valExtensions  {
	# Amber
        mdcrd top    parm top    parm7 top    prmtop top
        crd   trj    nc   trj
        # NAMD/CHARMM
        cor top    psf top
        dcd trj
        # GROMOS
        gro top    g96 top    top top
        xtc trj    trr trj
	# General
        pdb top
    }

    # Make the procedures visible to namespace import
    namespace export loadTraj infoMolecule calcRgyr calc calcRmsf calcHb \
                     calcHbNativeTrans calcSs calcNatContacts distMap \
                     createSomeDirs timeFormat checkScripts checkSims \
                     createSimDir getExt extId extUniq endCronTime
}

################################################################################
# Creates some of the working dirs                                             #
################################################################################

proc ::pfoundBasic::createSomeDirs {} {
    set cDir 0
    foreach item $::workdirList {
        if {![file exists $item]} {
			puts "Creating directory $item ..."
            file mkdir $item

			# Group users can write on logsPath dir. This will allow the DW
			# admin to delete/create the simId.* which will instructed VMD what 
			# to do next
			if {$item == $::logsPath || $item == $::toProcessPath} {
				file attributes $item -permissions rwxrwxr--
			} else {
				file attributes $item -permissions rwxr-xr--
			}
    	    set cDir 1
		}
    }
    return $cDir
}

################################################################################
# Time format                                                                  #
################################################################################

proc ::pfoundBasic::timeFormat {} {
    return [clock format [clock seconds] -format "%a %b %d %T %Z %Y" -gmt true]
}

################################################################################
# Checks existence of needed dirs and scripts                                  #
################################################################################

proc ::pfoundBasic::checkScripts {path} {
    foreach item $::scriptsList {
        if {![file exists ${path}/$item]} {
	    	puts $::cronFile [::pfoundBasic::ErrCron $item]
	    	puts $::cronFile [::pfoundBasic::endCronTime]
            exit
		}
    }
    puts $::cronFile "Work scripts checked"
    return
}

################################################################################
# Outputs aborting job                                                         #
################################################################################

proc ::pfoundBasic::ErrCron {value} {
	return "ERROR: $value does not exist! Aborting cronjob."
}

################################################################################
# Returns cron end time                                                        #
################################################################################

proc ::pfoundBasic::endCronTime { } {
	return "\nCronjob ended at [::pfoundBasic::timeFormat]\n###---###\n\n"
}

################################################################################
# Which simId's are present?                                                   #
################################################################################

proc ::pfoundBasic::checkSims { } {
    variable baseName
    cd $::toProcessPath
    if {[catch {set fileList [glob *.*]} res]} {
	puts $::cronFile "No files to be processed"
	puts $::cronFile [::pfoundBasic::endCronTime]
	exit
    } else {
	puts $::cronFile "Available files: $fileList"
	set baseName [list]
	foreach i $fileList {
	    lappend baseName [file rootname $i]
	}
	set baseName [lsort -unique -dictionary $baseName]
	puts $::cronFile "Possible simID's to be processed: $baseName"
    }
    return
}

################################################################################
# Checks files and creates a dir $processedPath/simId                          #
# for each SimId eligible to be processed by VMD, moving the relevant          #
# trajectory files into that dir.                                              #
#																			   #
# 2Do: presently, vmd checks if there is at least one trajectory and one       #
#      topology file, and does that by the file extension. This proc must be   #
#      improved in order to give more flexibility to the system (think in the  #
#      AMBER case for example).												   #
################################################################################

proc ::pfoundBasic::createSimDir { } {
    variable baseName
    #variable workSims
    variable base
    set workSims [list]
    foreach base $baseName {
		puts $::cronFile "\t--> $base <--"
		set workFiles [glob $base.*]
    
		# case1: 1 file  --> not enough to run vmd
    	# case2: 2 files --> topology/coordinates + traj
    	if {[llength $workFiles] == 2} {
			puts $::cronFile "Seems eligible to proceed."
	    	puts $::cronFile "Creating working dir $::processedPath/$base"
	    	lappend workSims $base

	    	# checks if simID directory exists
	    	if {[file exists ${::processedPath}/$base]} {
				file delete -force ${::processedPath}/$base
	    	}
    		file mkdir ${::processedPath}/$base

    		# checks file size
    		foreach elem $workFiles {
				puts $::cronFile "Checking if $elem file is being modified"
        		set answer [::pfoundBasic::CheckFileInfo $elem]

				# 2do: what about if pfounddb is not able to deliver
				#      the full file? Probably should be a good idea
				#      to have here a timeout operation to get out of the loop
				while {$answer == "notok"} {
            		set answer [::pfoundBasic::CheckFileInfo $elem]
        		}
    		}
		} else {
    
    		# not all needed files are present
    		puts $::cronFile "SimID $base will be processed on next cron."
		}
    }
    if {[llength $workSims] == 0} {
		puts $::cronFile "No eligible simId's to be processed!"
		puts $::cronFile [::pfoundBasic::endCronTime]
		exit
    } else {
		return $workSims
    }
}

################################################################################
# Test size file. Is the file being modified while it is copied to behemoth?   #
################################################################################

proc ::pfoundBasic::CheckFileInfo {elem} {
    variable nSeconds $::nSeconds
    variable base
    set sizeA [file size $elem]
    ::pfoundBasic::Sleep $nSeconds
    set sizeB [file size $elem]
    if {$sizeA == $sizeB} {
        puts $::cronFile "Moving $elem to ${::processedPath}/$base"
	exec mv $elem ${::processedPath}/$base
        return ok
    } else {
        puts $::cronFile "File $elem is being modified"
        puts $::cronFile "Sleeping for $nSeconds second(s)"
        return notok
    }
}

################################################################################
# Sleep for N seconds                                                          #
################################################################################

proc ::pfoundBasic::Sleep {N} {
    after [expr {int($N * 1000)}]
}

################################################################################
# Deletes all loaded molecules                                                 #
################################################################################

proc ::pfoundBasic::deleteMolecules {} {
    if {[molinfo num] != 0} {
        set molList [molinfo list]
        foreach item $molList {
            mol delete $item
        }
    }
    return
}

################################################################################
# Loads relevant files to a VMD session                                        #
# 2do: test with different formats                                             #
#      test topology/coordinate                                                #
#      test top2psf.pl                                                         #
#      take care of crd/crdbox like trajectories (AMBER)                       #
################################################################################

proc ::pfoundBasic::loadTraj {topology trajectory} {
    variable trajFrames
    
    # deals with .top/.itp gromacs files --> .psf
    if {([file extension $topology] == ".top") || \
	 ([file extension $topology] == ".itp")} {
	puts $::cronFile "Converting [file extension $topology] to .psf"
	::pfoundBasic::Top2psf $topology  
	set topology [file rootname $topology].psf
    }
    mol new $topology waitfor all
    set topFrames [molinfo top get numframes]
    
    # ensures that frame 0 from traj, is the 1st frame
    if {$topFrames > 0} {
	puts $::cronFile "$topology has coordinates. Going to delete them!"
        mol addfile $trajectory waitfor all
        animate delete beg 0 end 0 top
    } else {
        mol addfile $trajectory waitfor all
    }
    variable trajFrames [molinfo top get numframes]
    return
}

################################################################################
# Call top2psf.pl script                                                       #
# Marc Baaden is acknowledged                                                  #
################################################################################

proc ::pfoundBasic::Top2psf {filename} {
    exec ${::scriptsPath}/top2psf.pl -p $filename \
									 -o [file rootname $filename].psf
}

################################################################################
# Gather some information about the loaded molecule                            #
#2do: get structure vs resid (frame 0)                                         #
#     calculate number residues                                                #
#     calculate number of chains                                               #
################################################################################

proc ::pfoundBasic::infoMolecule {filename} {
    variable codes
    set fileId [open ./${filename}.info w 0640]

    #puts $fileId "Number of frames  : [molinfo top get numframes]"
    #puts $fileId "Number of atoms   : [molinfo top get numatoms]"
    #set sel0 [atomselect top "protein"]
    #set resids [$sel0 get resid] 
    #puts $fileId "Number of residues: [llength $resids]\n"
    # get some molecule info (uses last frame)

    set sel [atomselect top "protein name CA"]
    set ids [$sel get {resid resname chain segname}]
    puts $fileId "Resid\tResname\tChain\tSegment"
    foreach item $ids {
        set resnameVar [::pfoundBasic::CheckCodes [lindex $item 1]]
        puts $fileId \
			"[lindex $item 0]\t$resnameVar\t[lindex $item 2]\t[lindex $item 3]"
    }
    close $fileId
    return
}

################################################################################
# Defines resname in accordance with codes array                               #
################################################################################

proc ::pfoundBasic::CheckCodes {resname} {
    variable codes
    if {![catch {set result $codes($resname)}]} {
		return $codes($resname)
    } else {
		return "UNK"
    }
}

################################################################################
# Calculates the gyration radius for $selection                                #
################################################################################

proc ::pfoundBasic::calcRgyr {selection filename} {
    variable trajFrames
    set sel [atomselect top $selection frame 0]
    if {[$sel num] <= 0} {
        error "calcRgyr: must have at least one atom in selection"
    }
    set numfrms $trajFrames
    set nframe 0
    set fileId [open ./$filename w 0640]
    puts $fileId [::pfoundBasic::XyHeader $::outExt]
    while {$nframe < $numfrms} {
        set sel [atomselect top $selection frame $nframe]
        set rgyr [measure rgyr $sel]
        puts $fileId "$nframe\t$rgyr"
        incr nframe
    }
    close $fileId
}

################################################################################
# Calculates the RMSD for $selection                                           #
################################################################################

proc ::pfoundBasic::calcRmsd {selection filename} {
    variable trajFrames

    set sel0 [atomselect top $selection frame 0]
    if {[$sel0 num] <= 0} {
        error "calcRmsd: must have at least one atom in selection"
    }
    set numfrms $trajFrames
    set nframe 0
    set fileId [open ./$filename w 0640]
    puts $fileId [::pfoundBasic::XyHeader $::outExt]
    while {$nframe < $numfrms} {
        set sel1 [atomselect top $selection frame $nframe]
        set tm [measure fit $sel1 $sel0]
        set move_sel [atomselect top "all" frame $nframe]
        $move_sel move $tm
        set rmsd [measure rmsd $sel0 $sel1]
        puts $fileId "$nframe\t$rmsd"
        incr nframe
    }
    close $fileId
}

################################################################################
# Calculates the cumulative RMSF for $selection                                #
# 2do: implement the RMSF in time intervals --> RMSF vs time                   #
#      check if vmd >= 1.8.4                                                   #
################################################################################

proc ::pfoundBasic::calcRmsf {selection filename start {step 1}} {
    variable trajFrames
    set sel [atomselect top $selection frame 0]
    set selnum [$sel num] 
    set numframes $trajFrames
    if {[$sel num] <= 0} {
        error "calcRmsf: must have at least one atom in selection"
    }
    set fileId [open ./$filename w 0640]
    puts $fileId [::pfoundBasic::XyHeader $::outExt]
    set rmsf \
    	[measure rmsf $sel first $start last [expr $numframes - 1] step $step]
    set fr 1
    foreach i $rmsf {
        puts $fileId "$fr\t$i"
        incr fr
    }
    close $fileId
}

################################################################################
# Calculates the total number of H-bonds on $selection                         #
################################################################################

proc ::pfoundBasic::calcHb {output {selection "all"} {dist 3.0} {angle 20.0}} {
    variable trajFrames
    set nframes $trajFrames
    set fileId [open ./$output w 0644]
    puts $fileId [::pfoundBasic::XyHeader $::outExt]
    set sel [atomselect top $selection]
    for {set i 0} {$i < $nframes} {incr i} {
        $sel frame $i
        $sel update
        lassign [measure hbonds $dist $angle $sel] ldonors lacceptors lhindex
        set nbonds [llength $ldonors]
        puts $fileId "$i\t$nbonds"
    }
    close $fileId
}

################################################################################
# Calculates the number of H-bonds on $selection (frame 0),                    #
# so called native H-bonds (NHBs), and counts their persistence                #
# For those  NHBs, a matrix is built (H-bond ID vs frame),                     #
# to track the persistence of each NHB.                                        #
# 1 = NHB present    0 = NHB not present                                       #
# Distance = distancia entre os atomos dador e receptor                        #
# Angle = desvio ao alinhamento D-H:::R                                        #
# 2do: check if its possible to simplify this proc                             #
################################################################################

proc ::pfoundBasic::CalcHbNative {output {selection "all"} {dist 3.0} \
	{angle 20.0}} {
    variable nbonds
    variable trajFrames
    
    # do I need this animate comand?
    #animate goto start
    
    set nframes $trajFrames
    set refsel [atomselect top "$selection" frame 0]
    # do I need this 2 next lines? Check it!
    #$refsel frame 0
    #$refsel update
    
    lassign [measure hbonds $dist $angle $refsel] refdad refrec refhyd
    set refpairs [::pfoundBasic::GenPairs $refdad $refrec]
    set fileId [open ./${output}.nhbi w 0644]
    variable nbonds [llength $refdad]
    set ilist [list]
    for {set i 0} {$i < $nbonds} {incr i} {
        set atmdad [lindex $refdad $i]
        set atmrec [lindex $refrec $i]
        set seldad [atomselect top "index $atmdad"]
        set selrec [atomselect top "index $atmrec"]
        set iddad \
			[$seldad get "chain"]:[$seldad get "resname"]:[$seldad get "resid"]:[$seldad get "name"]
        set idrec \
			[$selrec get "chain"]:[$selrec get "resname"]:[$selrec get "resid"]:[$selrec get "name"]
        set j [expr $i+1]
        puts $fileId "$j $iddad $idrec"
        lappend ilist $j
    }
    close $fileId
    set fileId [open ./${output}.nhbm_tmp w 0644]
    set fileId1 [open ./${output}.nhb w 0644]
    puts $fileId1 [::pfoundBasic::XyHeader nhb]
    set sel [atomselect top "$selection"]
    for {set i 0} {$i < $nframes} {incr i} {
        $sel frame $i
        $sel update
        lassign [measure hbonds $dist $angle $sel] dad rec hyd
        set pairs [::pfoundBasic::GenPairs $dad $rec]
        set hblist [list]
        set counter 0
        foreach value $refpairs {
            
	    # hb=0 => bond does not exists; hb=1 => bond exists	
            if {[lsearch $pairs $value] == -1} {
                set hb 0
            } else {
                set hb 1
                set counter [expr $counter + 1]
            }
            lappend hblist $hb
        }
        puts $fileId "$hblist"
        puts $fileId1 "$i\t$counter"
    }
    close $fileId
    close $fileId1
}

################################################################################
# Proc used by CalcHbNative                                                    #
################################################################################

proc ::pfoundBasic::GenPairs {{list1} {list2}} {
    set pairs [list]
    foreach ai $list1 bi $list2 {lappend pairs "$ai $bi"}
    return "$pairs"
}

################################################################################
# Transposes a matrix using the La package                                     #
# Did not used transtranspose because the                                      #
# La:show output is easier to work with.                                       #
################################################################################

proc ::pfoundBasic::Transpose {infile fileExt nrows ncols} {
    set inFile [open $infile r]
    set outFile [open [file rootname ${infile}].$fileExt w]
    lappend lineList 2 $nrows $ncols
    
    # The read loop command is faster than the gets loop command
    # www.beedub.com.com/book/2nd/unix.d???????
    foreach line [split [read $inFile] \n] {
        lappend lineList $line
    }
    set matrix [join $lineList]
    set matrix [La::transpose $matrix]
    puts $outFile [La::show $matrix]
    close $outFile
}

################################################################################
# Calls CalcHbNative and transposes the outputed matrix                        #
################################################################################

proc ::pfoundBasic::calcHbNativeTrans {output {selection "all"} {dist 3.0} \
	{angle 20.0}} {
    variable trajFrames
    variable nbonds
    ::pfoundBasic::CalcHbNative $output $selection $dist $angle
    ::pfoundBasic::Transpose ${output}.nhbm_tmp nhbm $trajFrames $nbonds
    file delete ${output}.nhbm_tmp
}

################################################################################
# Calculates the secondary structure for $selection using stride               #
# 2do: do a parser to change output, or rewrite this proc                      #
#      Vitaliy asks to change the output letters: check e-mail 27abr2006       #
################################################################################

proc ::pfoundBasic::CalcSs {filename selection} {
    variable trajFrames
    set fileId [open ./${filename}.ss_tmp w 0640]
    set sel [atomselect top "$selection and name CA"]
    for {set frame 0} {$frame < $trajFrames} {incr frame} {
        animate goto $frame
		display update ui
        $sel frame $frame
		$sel update
        mol ssrecalc top
        set structList [$sel get structure]
		puts $fileId "$structList"
        unset structList
    }
    close $fileId
}

################################################################################
# Calls CalcSs and transposes the outputed matrix                              #
################################################################################

proc ::pfoundBasic::calcSsTrans {filename selection} {
    variable trajFrames
    set sel [atomselect top "$selection and name CA"]
    set selNum [$sel num]
    ::pfoundBasic::CalcSs $filename $selection
    ::pfoundBasic::Transpose ${filename}.ss_tmp ss $trajFrames $selNum 
    file delete ${filename}.ss_tmp
}

################################################################################
# Calculates native contacts                                                   #
################################################################################

proc ::pfoundBasic::calcNatContacts {selection output {cutoff 4.2} \
	{neighbours 4}} {
    set numfrms [molinfo top get numframes]
    set sel [atomselect top "$selection"]
    set refcontacts [::pfoundBasic::GetResContacts $sel 0 $cutoff $neighbours]
    set nrefcontacts [llength $refcontacts]
    set fileID [open ${output}.nc w 0644]
    puts $fileID [::pfoundBasic::XyHeader nc]
#    puts $fileID "#Selection: $selection     Cutoff: $cutoff   Exclude \
#	neighbours: $neighbours"
    set nframe 0
    while {$nframe < $numfrms} {
        set count 0
        for {set i 0} {$i < [llength $refcontacts]} {incr i} {
            set selA [atomselect top "resid [lindex [lindex $refcontacts $i] 0]\
	        	and $selection" frame $nframe]
            set selB [atomselect top "resid [lindex [lindex $refcontacts $i] 1]\
	        	and $selection" frame $nframe]
            lassign [measure contacts $cutoff $selA $selB] la lb
            $selA delete
            $selB delete
            if {[llength $la] > 0} {incr count}
        }
        # outputs the relative percentage of native contacts to frame 0 vs time
	#puts $fileID "$nframe\t[expr 1.0 * $count / $nrefcontacts]"
	# outputs the number of native contacts vs time
	puts $fileID "$nframe\t$count"
        incr nframe
    }
    close $fileID
}

################################################################################
# Proc used by calcNatContacts                                                 #
# 2do: check if its possible to simplify this proc                             #
################################################################################

proc ::pfoundBasic::GetResContacts {sel frame cutoff neighbours} {
    $sel frame $frame
    $sel update
    set index [$sel get index]
    set resid [$sel get resid]
    set seq [lsort -integer -unique $resid]
    set nres [llength $seq]

    #Retrieve lists with atoms in contact
    lassign [measure contacts $cutoff $sel] la lb
    set nla [llength $la]
    set fullreslist [list]

    #Change atom-atom list -> residue-residue list
    for {set i 0} {$i < $nla} {incr i} {
        set resA [lindex $resid [lsearch $index [lindex $la $i]]]
        set resB [lindex $resid [lsearch $index [lindex $lb $i]]]
        lappend fullreslist [lsort -integer "$resA $resB"]
    }
    set reslist [lsort -unique $fullreslist]

    #Eliminate neighbours from residue-residue list
    set rescontacts [list]
    for {set i 0} {$i < $nres} {incr i} {
        set resA [lindex $seq $i]
        for {set j $i} {$j<$nres} {incr j} {
            set resB [lindex $seq $j]
            set val [lsearch $reslist "$resA $resB"]
            if {$val != -1 && [expr abs($resA-$resB)] > $neighbours} {
                lappend rescontacts "$resA $resB"
            }
        }
    }
    return $rescontacts
}

################################################################################
# Calculates a CA-CA distance map on $selection                                #
################################################################################

proc ::pfoundBasic::distMap {selection out {fr 0}} {
    set sel [atomselect top "$selection" frame $fr]
    set reslist [lsort -integer -unique [$sel get resid]]
    set resnum [llength $reslist]
    #puts "Computing contacts map for $resnum residues..."
    set fileId [open ./${out}.dm w 0644]
    #puts $fileId "resA resB dist"
    for {set i 0} {$i < $resnum} {incr i} {
        set resA [lindex $reslist $i]
        set selA [atomselect top "$selection and resid $resA" frame $fr]
        set coordsA [$selA get "x y z"]
        set distlist1 [list]
        for {set j 0} {$j < $resnum} {incr j} {
            set resB [lindex $reslist $j]
            set selB [atomselect top "$selection and resid $resB" frame $fr]
            set coordsB [$selB get "x y z"]
            
	  	  	# calculates all distances between atoms of 2 residues
            set distlist [list]
            foreach coordA $coordsA {
                foreach coordB $coordsB {
                    lappend distlist [vecdist $coordA $coordB]
                }
            }
	    
	    	# Selects the minimum distance from $distlist
            set mindist [format %6.2f [lindex [lsort -real $distlist] 0]]
            lappend distlist1 $mindist
        }
        set distlist1 [join $distlist1]
        puts $fileId $distlist1
    }
    close $fileId
}

################################################################################
# Produces a XY graphic using xmgrace                                          #
# Image produced: PNG image data, 792 x 612, 8-bit colormap, non-interlaced    #
################################################################################

proc ::pfoundBasic::xyToGrace {infile {hdevice PNG}} {
    ::pfoundBasic::exec gracebat -hardcopy $infile -hdevice $hdevice \
    	                         -printfile ${infile}.[string tolower $hdevice]
}


################################################################################
# Drop-in replacement for exec; before running the program, it first	       #
# searches for the executable and prompts the user for its location if it's    #
# not found.																   #
# Based on proc ::ExecTool::exec { args } [vmd plugin exectools 1.2]	       #
################################################################################
proc ::pfoundBasic::exec { args } {
  #puts "\n$args"
  if {[llength $args] < 1} {
    error "Insufficient arguments."
  }
  set exec_name [lindex $args 0]
  set exec_path [auto_execok $exec_name]
  #puts "exec path = $exec_path"
  if {$exec_path == {}} {
    error "Couldn't find program `$exec_name'."
  }

  eval ::exec [list $exec_path] [lrange $args 1 end]
}

proc ::pfoundBasic::PfoundTag {} {
    return "@ with string\
    	  \n@   string on\
    	  \n@   string 0.28, 0.02\
    	  \n@   string def \"\\-Generated by \\f{3}Pfound\\f{}: \
		[::pfoundBasic::timeFormat] - Simulation ID: $::sim\""
}

################################################################################
# Writes a header to XY property file (xmgrace format)                         #
################################################################################

proc ::pfoundBasic::XyHeader {extension} {
    switch -exact -- $extension {
	rmsd {
	      return "@ title \"Atom positional RMSD from 1\\Sst\\N\
			frame structure\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"RMS deviation (\\#{C5})\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	     }

	rmsf {
	      return "@ title \"Displacement of C\\s\\f{12}a\\N\\f{} atoms\"\
	            \n@ xaxis label \"Residue\"\
		    \n@ yaxis label \"RMS fluctuation (\\#{C5})\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}
	
	rgyr {
	      return "@ title \"Radius of gyration\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"Rgyr (\\#{C5})\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}
	
	hb   {
	      return "@ title \"Hydrogen-bonds\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"Number of H-bonds\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}
	
	nhb  {
	      return "@ title \"Native Hydrogen-bonds\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"Number of H-bonds\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}

	nc   {
	      return "@ title \"Native Contacts\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"Number of native contacts\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}
	
	sasa_aa   {
	      return "@ title \"All-atoms SASA\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"SASA (\\#{C5}\\S2\\N)\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}
	
	sasa_anp   {
	      return "@ title \"Non-Polar SASA\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"SASA (\\#{C5}\\S2\\N)\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}
	
	sasa_ap   {
	      return "@ title \"Polar SASA\"\
	            \n@ xaxis label \"Time (ps)\"\
		    \n@ yaxis label \"SASA (\\#{C5}\\S2\\N)\"\
		    \n[::pfoundBasic::PfoundTag]\
		    \n@ TYPE xy\
		    \n@ s0 line linewidth 2\
		    \n@ s0 line color 2"
	}

	default {
	      return
	}
    }
}

#	ssc  {
#	      return "# Produced by P-found: [::pfoundBasic::timeFormat]\
#		    \n@ title \"Counting secondary structure motifs\"\
#	            \n@ xaxis label \"Time (ps)\"\
#		    \n@ yaxis label \"Number of residues\"\
#		    \n@ TYPE xy"
#	}

################################################################################
# Outputs file extensions of nameDir files                                     #
################################################################################

proc ::pfoundBasic::CheckFileExtensions {nameDir} {
    set listExtensions [list]
    set listDir [glob -type f *{${nameDir}.}*]
    foreach i $listDir {
        lappend listExtensions [::pfoundBasic::GetExtension $i]
    }
    return $listExtensions
}

################################################################################
# Gets infile extension, without the dot                                       #
################################################################################

proc ::pfoundBasic::GetExtension {infile} {
    return [lindex [split [file extension $infile] .] 1]
}

################################################################################
# Parser between file extension and top || trj                                 #
################################################################################

proc ::pfoundBasic::AvailableExtensions {fileTag} {
    variable valExtensions
    if {[info exists valExtensions($fileTag)]} {
        return [list $fileTag $valExtensions($fileTag)]
    } else {
    	return
    }
}

################################################################################
# Gets file extension of files on $Dir                                         #
################################################################################

proc ::pfoundBasic::getExt {Dir} {
    cd ${::processedPath}/$Dir
    #puts $::cronFile "\t--> $Dir <--"
    return [::pfoundBasic::CheckFileExtensions $Dir]
}

################################################################################
# Parsing extension to file type: top || trj                                   #
################################################################################

proc ::pfoundBasic::extId {sim item} {
    variable topType  
    variable trajType 
    set extType [::pfoundBasic::AvailableExtensions $item]
    if {[lindex $extType 1] == "top"} {
		variable topType [lindex $extType 0]
		puts $::cronFile "${sim}.$item: identified as a topology file"
    } elseif {[lindex $extType 1] == "trj"} {
        variable trajType [lindex $extType 0]
		puts $::cronFile "${sim}.$item: identified as a trajectory file"
    } else {
    	puts $::cronFile "File(s) have unknown extension!"
    }
    return
}

################################################################################
# Outputs the topology and trajectory file, in this order                      #
################################################################################

proc ::pfoundBasic::extUniq { } {
    variable topType
    variable trajType

    # Check if there's a topology || coordinates & trajectory
    puts -nonewline $::cronFile "Checking if top\/traj are uniquely defined ..."
    foreach item {topType trajType} {
	if {![info exists $item]} {
	    puts $::cronFile " notok\nSkipping this data set."
	    return "notok"
	}
    }
    puts $::cronFile " ok"
    return [list $topType $trajType]
}

################################################################################
# Not in use                                                                   #
# tests if $nameSpace is available                                             #
################################################################################

proc ::pfoundBasic::searchNamespace {nameSpace} {
    if {[lsearch ::$nameSpace:: [namespace children]] == -1} {
		return
    } else {
		continue
    }
}

################################################################################
# Have to compare performance with measure sasa                                #
# Calculates residue solvent accessibility (naccess)                           #
# 2do: think about how to deal with the produced data						   #
#      - put abs/rel sasa per molecule in just one file?					   #
#      - how about the sasa per residue?									   #
################################################################################

proc ::pfoundBasic::calcSasa {output {selection protein} {zslice 0.1}} {
    variable trajFrames
    global tmpdir
    global dataBasename
    # SASA per residue
    #set file1Id [open ./${output}.sasa_abs w 0640]
    #set file2Id [open ./${output}.sasa_rel w 0640]
    #set file4Id [open ./${output}.sasa_abs_bb w 0640]
    #set file5Id [open ./${output}.sasa_abs_sc w 0640]

    # SASA per molecule (1 value/frame) 
    #set file3Id [open ./${output}.sasa_tot w 0640]
    # absolute all-atoms sasa
    set file6Id [open ./${output}.sasa_aa w 0640]
    puts $file6Id [::pfoundBasic::XyHeader sasa_aa]
    # absolute polar-atoms sasa 
    set file7Id [open ./${output}.sasa_anp w 0640]
    puts $file7Id [::pfoundBasic::XyHeader sasa_anp]
    # absolute non-polar atoms sasa
    set file8Id [open ./${output}.sasa_ap w 0640]
    puts $file8Id [::pfoundBasic::XyHeader sasa_ap]

    set sel [atomselect top "$selection and not hydrogen"]
    set count 0

    # Going to work on global temp dir
    set tmpDir $tmpdir
    set currDir [exec pwd]
    cd $tmpDir
    while {$count < $trajFrames} {
	#puts "Working on frame $count..."
	$sel frame $count
	$sel update
	set outfile [format naccess_${dataBasename}_%04d [expr $count]]
	$sel writepdb ${outfile}.pdb

	# Edit pdb file to match naccess conventions for residue and atom names
	# 2do: work on a parser
	exec mv ${outfile}.pdb naccess_${dataBasename}.pdb
	exec sed -e "s/HSD/HIS/g"   -e "s/HSE/HIS/g"   -e "s/CD  ILE/CD1 ILE/g"\
                 -e "s/OT1 /O   /g" -e "s/OT2 /O   /g" \
                 naccess_${dataBasename}.pdb > ${outfile}.pdb
	
        # starts sasa calc on ${outfile}.pdb
	::pfoundBasic::exec naccess ${outfile}.pdb -z $zslice

        #set abs [pfoundBasic::GetSasa ${outfile}.rsa 4]
	#set rel [pfoundBasic::GetSasa ${outfile}.rsa 5]
	#set bb  [pfoundBasic::GetSasa ${outfile}.rsa 8]
	#set sc  [pfoundBasic::GetSasa ${outfile}.rsa 9]
	#set accesstot  [exec cat ${outfile}.rsa | awk \
    # "/^TOTAL/ {print \"\t\" \$2 \"\t\" \$3 \"\t\" \$4 \"\t\" \$5 \"\t\" \$6}"]
	set access_aa  [exec cat ${outfile}.rsa | awk "/^TOTAL/ {print \$2}"]
    	set access_anp [exec cat ${outfile}.rsa | awk "/^TOTAL/ {print \$5}"]
	set access_ap  [exec cat ${outfile}.rsa | awk "/^TOTAL/ {print \$6}"]

	#puts $file1Id "$count $abs"
	#puts $file2Id "$count $rel"
	#puts $file3Id "$count $accesstot"
	#puts $file4Id "$count $bb"
	#puts $file5Id "$count $sc"
	puts $file6Id "$count\t$access_aa"
	puts $file7Id "$count\t$access_anp"
	puts $file8Id "$count\t$access_ap"
	
	# Será mm necessário apagar estes files? Será que um overwrite não é mais rápido
	file delete naccess_${dataBasename}.pdb ${outfile}.pdb ${outfile}.asa \
                    ${outfile}.log ${outfile}.rsa
	incr count
    }
    #close $file1Id
    #close $file2Id
    #close $file3Id
    #close $file4Id
    #close $file5Id
    close $file6Id
    close $file7Id
    close $file8Id
    cd $currDir
}

################################################################################
# Proc used by calcSasa                                                        #
################################################################################

proc ::pfoundBasic::GetSasa {rsa_file column} {
    set fileID [open $rsa_file r]
    set col [list]
    while {[gets $fileID line] >= 0} {
        if ([string compare "RES" [lindex $line 0]]==0) {
		append col "\x9 [lindex $line $column]"
	}
    }
    close $fileID
    return "$col"
}

################################################################################
# Calculates SS content for each frame (helix/sheet)						   #
# code core taken from VMD-L:												   #
# http://www.ks.uiuc.edu/Research/vmd/mailing_list/vmd-l/5107.html			   #
################################################################################

proc ::pfoundBasic::calcSsCounts {output} {
    set numframes [molinfo top get numframes]
    set sel [atomselect top "name CA"]
    set fileId [open ./${output}.ssc w 0640]
    for {set frame 0} {$frame < $numframes} {incr frame} {
	animate goto $frame
	vmd_calculate_structure top
	$sel frame $frame
	$sel update
	set helixlist [$sel get alpha_helix]
	set sheetlist [$sel get sheet]
	set helixcount 0
	foreach i $helixlist { incr helixcount $i }
	set sheetcount 0
	foreach i $sheetlist { incr sheetcount $i }
	puts $fileId "$frame\t$helixcount\t$sheetcount"
        display update ui
    }
    $sel delete
    close $fileId
}
