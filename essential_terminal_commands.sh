zip -r varmodel2.zip . -x *.git* varmodel2.xcodeproj/\*

# SLURM get job info
sacct -j 21777695 --format=jobid,jobname,partition,account,alloccpus,state,cputime,maxrss,maxvmsize --state=RUNNING --starttime 2018-04-30T18:15:00

# copy all sqlite files across directories from scratch to local computer
for i in {7..12}
do
rsync -rav -e ssh --include '*/' --include='mtn_'$i'*.sqlite' --exclude='*' \.  \  shai@192.170.193.140:/home/shai/Documents/sqlite
done

rsync -rav -e ssh --include '*/' --include='*.RData' --exclude='*' \.  \  pilosofs@midway.rcc.uchicago.edu:/scratch/midway/pilosofs

# Comapre folders
rsync -vn --include='*.csv' /scratch/midway/pilosofs/mtn_4 shai@ee-jazz.uchicago.edu:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_4


# Delete files across folders
find . -name "*.RData" -type f -delete

# Count files across folders
find . -name "*.zip" -type f | wc -l

# Copy all RData files from multiple folders to one location
find . -name *.RData -type f -exec cp '{}' /home/pilosofs/malaria_temporal_networks/ \;

# Remove only folders but not files
rm -r */

# zip selected files
zip -r sqlite mtn_4 -i '*.sqlite'



# Delete files in folders
for i in {4..15}
do
cd /home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i
rm *.pdf
rm *.csv
done

# Copy files from a given date
for i in {4..15}
do
	echo 'mtn_'$i
	localDir=/scratch/midway/pilosofs/mtn_$i
	cd $localDir
	touch --date "2016-12-20" start
	touch --date "2016-12-21" end
	scp `find -type f -newer start -not -newer end -name '*.pdf'` shai@128.135.150.227:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i
	scp `find -type f -newer start -not -newer end -name '*.csv'` shai@128.135.150.227:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i
done


# loop for zip
for i in {25..28}
do
echo 'mtn_'$i
localDir=/scratch/midway/pilosofs/mtn_$i
cd $localDir
for r in {1..6}
do
echo '-- run_'$r
cd localDir/$r
zip 'mtn_"$i'_run_'$r'.zip -i '*.sqlite'
mv 'mtn_"$i'_run_'$r'.zip  /scratch/midway/pilosofs/sqlite
done
done



for i in {2..10}
do
echo 'run_'$i
zip 'mtn_12_seasonality_empirical_irs_2_run_'$i'.zip' *.* -i '*run_'$i'_*.*' -m
done



# Create a copy of sqlites while renaming the experiment
for i in {4..15}
do
	echo 'mtn_'$i
	localDir=/scratch/midway/pilosofs/mtn_$i
	let x=i+14
	remoteDir=/scratch/midway/pilosofs/mtn_$x
	mkdir $remoteDir
	cd $localDir
	for r in {1..20}
	do
		echo '-- run_'$r
		mkdir $remoteDir/run_$r
		cd $localDir/run_$r
		cp *.sqlite $remoteDir/run_$r
		cd $remoteDir/run_$r
		mv 'mtn_'$i'_run_'$r'.sqlite' 'mtn_'$x'_run_'$r'.sqlite' 
	done
done


# Copy files with chaning names
for i in {1..24};  do
echo 'mtn_'$i
localDir=/scratch/midway/pilosofs/mtn_$i
let x=i+24
remoteDir=/scratch/midway/pilosofs/mtn_$x
#mkdir $remoteDir
cd $localDir
cp 'mtn_'$i'.sqlite' $remoteDir
cd $remoteDir
mv 'mtn_'$i'.sqlite' 'mtn_'$x'.sqlite'
cd ..

done

##
for i in {1..20}
do
cd mtn_$i
rm -r */
cd ..
done


##
localDir=/scratch/midway/pilosofs/mtn_7
remoteDir=/scratch/midway/pilosofs/mtn_17
for i in {1..20}
do
echo $i
cd $remoteDir
rm -rf run_$i # remove the run folder if it exists
mkdir run_$i
cd $localDir/run_$i
cp *.sqlite $remoteDir/run_$i
cd $remoteDir/run_$i
mv mtn_7_run_$i.sqlite mtn_17_run_$i.sqlite
done

#
for i in {4..15}
do
cd /home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i
mkdir 300_layers
mv *.csv 300_layers
mv *.pdf 300_layers
done



for i in {1..25}
do
echo $i
for r in {1..10}
do
../../Infomap 'layer_'$i'_realization_'$r'.txt' output/ --undirected -i pajek -s 123 -N 5 --two-level --map --tree --silent
done
done
	


for i in {4..5}
do
echo 'mtn_'$i
cd '/scratch/midway/pilosofs/mtn_'$i
for r in {1..20}
do
echo '-- run_'$r
unzip 'mtn_'$i'_run_'$r'_0.95.zip'
done
rm *.zip
done



type='mda'
for i in {27..28}
do
	echo 'mtn_'$i$type
	cd '/scratch/midway/pilosofs/mtn_'$i$type
	scp *run_1* shai@192.170.193.140:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i$type
	#scp *.csv shai@192.170.193.140:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i$type
	#scp *.pdf shai@128.135.150.227:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i$type
	#scp *.zip shai@128.135.150.227:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/zip
	#rm *.zip
done
Agam&Ofek215


type=''
for i in {25..28}
do
cd '/scratch/midway/pilosofs/mtn_'$i$type 
zip -r 'mtn_'$i$type'.zip' . -i 'run_*'
done


for i in 4 5 6 7
do
	echo 'mtn_'$i
	cd '/scratch/midway/pilosofs/mtn_'$i
	scp -r permutations/*.pdf shai@128.135.150.227:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i
	scp -r permutations/*.csv shai@128.135.150.227:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i
	scp -r permutations/*.tree shai@128.135.150.227:/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_$i
done


cd /scratch/midway/pilosofs; zip -r accumCurves.zip mtn_25 mtn_27 mtn_28 mtn_30 -i '*accumCurve*.*'


# BULK RENAME
for i in 26irs 27irs 28irs; 
	do cd mtn_$i;
	for r in {1..5};
		do cd run_$r;
		rename '0.95' '0.95_300' *.*;
		cd ..; 
	done; 
	cd ..;
done;


#This will replace the string SEARCH with REPLACE in every file (that is, *). The /g means global, so if you had a SEARCH SEARCH.jpg, it would be renamed REPLACE REPLACE.jpg. If you didn't have /g, it would have only done substitution once, and thus now named REPLACE SEARCH.jpg. If you want case-insensitive, add /i (that would be, /gi or /ig at the end).
rename "s/SEARCH/REPLACE/g"  *


type=''
for i in {25..28}
do
cd 'mtn_'$i$type 
rm *0.95_t*.*
rm *0.95_g*.*
rm *0.95_i*.*
rm *0.95.pdf
cd ..
done


## Copy from Midway without the run folders
cd /home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments
for i in {1..12}
do
echo $i
mkdir 'mtn_'$i'Nt'
done

for i in {1..12}
do
echo '----------  '$i'  -------------'
cd '/scratch/midway/pilosofs/mtn_'$i'Nt'
scp *.* shai@192.170.193.140:'/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_'$i'Nt'
done

Agam&Ofek215





for i in  {1..12}; 
do cd 'mtn_'$i;
for r in {6..10};
do cd 'run_'$r;
cp 'mtn_'$i'_run_'$r'_0.95_300_inter_edges.csv' ..; 
cp 'mtn_'$i'_run_'$r'_0.95_300_intra_edges.csv' ..; 
cp 'mtn_'$i'_run_'$r'_0.95_300_layerInfo.csv' ..; 
cp 'mtn_'$i'_run_'$r'_0.95_300_strain_composition.csv' ..;
cp 'mtn_'$i'_run_'$r'_0.95_300_followedHosts.csv' ..;
cp 'mtn_'$i'_run_'$r'_0.95_300_strain_compositionFH.csv' ..;
cp 'mtn_'$i'_run_'$r'_0.95_300_similarityMatrix.csv' ..;
cd ..
done; 
cd ..;
done;

# Loop through folders and do something.
for f in *;
do
    if [[ -d $f ]];
    then
		if [ $f != "mtn_12_seasonality_empirical_irs_2" ];
		then
			echo $f
			cd $f
			rm -r run_1
			cd ..
		fi
    fi
done



## Loop and extract particular files
for i in 2 2Nt
do
echo 'mtn_'$i
cd '/home/shai/Dropbox/RESEARCH/PROJECTS/Malaria_interventions/current_experiments/mtn_'$i
for r in {1..10}
do
echo '-- run_'$r
unzip -j 'mtn_'$i'_run_'$r'_0.25_300.zip' '*/*.tree'
done
done



## Move err and out files to their respective folders
for i in 2 7 12
do
for r in {11..50}
do mv 'mtn_'$i'_'*$r.* '/scratch/midway2/pilosofs/mtn_'$i'/run_'$r;
done;
done;

for i in 2 7 12
do
for r in {1..10}
do mv 'mtn_'$i'_Sn_'*$r.* '/scratch/midway2/pilosofs/mtn_'$i'_Sn/run_'$r;
done;
done;

for i in 2 7 12
do
for r in {1..10}
do mv 'mtn_'$i'Nt_'*$r.* '/scratch/midway2/pilosofs/mtn_'$i'Nt/run_'$r;
done;
done;



for i in 2_Sn 7_Sn
do
    cd '/home/shai/Documents/mtn_data/current_experiments/mtn_'$i
    zip 'mtn_'$i'_additional_result_files.zip' *.* -i '*.csv' -x '*diagnostics*.csv' '*EIR*.csv' '*MOI*.csv' '*PTS*.csv' '*layerInfo*.csv' '*moduleSummary*.csv' '*temporalData*.csv' '*durationOfInfection*.csv' -m -u
done




for i in 2_Sn 2Nt 2Nt_Sn
do
cd '/home/shai/Documents/mtn_data/current_experiments/mtn_'$i
unzip 'mtn_'$i'_additional_result_files.zip' '*temporalData*.csv'
zip -d 'mtn_'$i'_additional_result_files.zip' '*temporalData*.csv*'
done



for i in 2 7 12
do
for r in {1..50}
do cp 'mtn_'$i'_run_'$r'.json' '/scratch/midway2/pilosofs/mtn_'$i'/run_'$r;
done;
done;


for i in 2_Sn 7_Sn 12_Sn
do
for r in {1..50}
do
cd '/scratch/midway2/pilosofs/mtn_'$i'/run_'$r
cp *.tree /scratch/midway2/pilosofs/calculate_PTS
cp *strain_composition.csv /scratch/midway2/pilosofs/calculate_PTS
done
done

for f in *;
do
    if [[ -d $f ]];
    then
	cd $f
	cp *0.95_Infomap*.txt ../Infomap_input
	cd ..
	fi
done


for d in {1..100}
do
echo $d
cd '/project/pascualmm/shai/ghana_season/'$d
rm *_Infomap*.txt
rm *.R
rm NAbuggy_instances_infomap.txt
rm *repertoire_persistence.csv
rm *.pdf
rm Infomap*
done
cd '/project/pascualmm/shai/ghana_season/'
rm slurm_output/*.*
rm Results/*.*
rm Infomap_input/*.*
sbatch networks_for_data_job_all_scenarios.sbatch



cd '/project/pascualmm/shai/ghana_season'
zip repertoire_persistence.zip . -i '*repertoire_persistence*' -r
# Then open:
unzip -j repertoire_persistence.zip


# Extract all json files
for i in {1..50}
do
unzip 'mtn_12_run_'$i'_0.9_300.zip' '*.json'
done

for i in {1..50}
do
cd 'run_'$i
cp *.* ../
cd ..
done


for i in {1..50}
do
unzip -l 'mtn_12Nt_run_'$i'.zip'
done


for i in {1..50}
do
unzip -j 'mtn_12_Sn_'$i'_0.9_300.zip' '*.json'
done


for i in {1..50}
do
cd 'run_'$i
cp *.json ../../sqlite
cd ..
done

for i in {1..50}
do
zip 'mtn_7Gn_run_'$i'.zip' *.* -i 'mtn_7Gn_run_'$i'.*'
done

for r in {11..20}
do
echo $r
cd '/scratch/midway2/pilosofs/mtn_12Gn/run_'$r
cp *.tree /scratch/midway2/pilosofs/epi
cp *strain_composition.csv /scratch/midway2/pilosofs/epi
cp *EIR.csv /scratch/midway2/pilosofs/epi
done


for r in {11..25}
do
unzip -j 'mtn_7Gn_run_'$r'_0.8_300.zip' '*.tree' '*strain_composition.csv'
done


for i in {1..50}
do
unzip -j 'mtn_2_run_'$i'.zip'
done


zip mtn_7Gn_additional_result_files.zip *.* -i '*.csv' -x '*diagnostics*.csv' '*EIR*.csv' '*MOI*.csv' '*PTS*.csv' '*layerInfo*.csv' '*moduleSummary*.csv' '*temporalData*.csv' '*durationOfInfection*.csv' -m -u
