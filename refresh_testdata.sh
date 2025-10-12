DIR_DATA=/fast/peijinz/day-peel/data
DIR_DATA_SRC=/fast/peijinz/day-peel/data.tar # tar file of the testdata
DIR_OUTPUT=/fast/peijinz/day-peel/

# remove the directory
rm -rf $DIR_DATA
# untar the testdata
tar -xvf $DIR_DATA_SRC -C $DIR_OUTPUT