
# create the dirs where a file README.sh is on /work/project/fragencode
dir=/work/project/fragencode
find $dir -name README.sh | while read f 
do 
file=${f#$dir"/"}
echo $file | awk '{s=""; n=split($1,a,"/"); for(k=1; k<=(n-1); k++){s=(s)(a[k])("/"); print s} }' 
done | while read d
do
mkdir -p ~/save/READMEs/$d
done

# copy the readmes in the proper places
find $dir -name README.sh | while read f
do 
file=${f#$dir"/"}
cd ~/save/READMEs/${file%/README.sh}
cp $dir/$file .
done

# put the correct rights back
chmod 755 -R ~/save/READMEs
