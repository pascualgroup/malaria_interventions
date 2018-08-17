cd /scratch/midway2/pilosofs/malaria_interventions/sqlite

# Move all E000 files to a folder
rm e000files.txt
find *E000* ! -name '*CP*' >> e000files.txt
xargs mv -t E000 < e000files.txt

# Move all CP files to a folder
mv *_CP*.sqlite CP
