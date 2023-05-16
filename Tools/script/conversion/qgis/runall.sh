

#Directory script
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"


echo 'Run meshToShapefile.py'
python $DIR/meshToShapefile.py

echo ''
echo 'Run resultToShapefileAll.py'
python $DIR/resultToShapefileAll.py

echo ''
echo 'Run stationObservationToShapefile.py'
python $DIR/stationObservationToShapefile.py
