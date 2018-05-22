module load git
cd /home/pilosofs/varmodel2/
git pull
git log -1
rm *.zip
zip -r varmodel2.zip . -x *.git* varmodel2.xcodeproj/\*
zip -d varmodel2.zip "parameters-example.py"
cp varmodel2.zip ../malaria_interventions/
cd ../malaria_interventions/
