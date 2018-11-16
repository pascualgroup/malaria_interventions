module load git
cd /home/pilosofs/varmodel2/
git pull
git log -1 > 'version.txt'
rm *.zip
zip -r varmodel2.zip . -x *.git* varmodel2.xcodeproj/\* tests/\* parameters-example.py sqlite3o_err.txt run_tests.py 
cp varmodel2.zip /scratch/midway2/pilosofs/malaria_interventions/
cp varmodel2.zip /scratch/midway2/pilosofs/PLOS_Biol/
cd /scratch/midway2/pilosofs/PLOS_Biol/
